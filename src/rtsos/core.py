# ===============================================================
#  RTSOS – Radiative Transfer model based on Successive Orders of Scattering
#  Copyright (C) 2025  Pengwang Zhai

#  Licensed under the Creative Commons Attribution–NonCommercial 4.0
#  International License (CC BY-NC 4.0).
#  You may use, modify, and share this code for research and
#  educational purposes with proper attribution.
#  Commercial use requires written permission from the author.

#  Full license: https://creativecommons.org/licenses/by-nc/4.0/
#  Contact: Pengwang Zhai  |  [pwzhai@gmail.com]
# ===============================================================

import subprocess
from argparse import Namespace
from collections.abc import Iterable
from pathlib import Path
from platform import system
from shutil import copy
from tempfile import TemporaryDirectory
from uuid import uuid4

import dask.config
import numpy as np
import xarray as xr
from dask.tokenize import tokenize

from .execution import (
    Identity,
    ExecutionFailure,
    archive_output,
    atomic_write,
    output_lock,
    publish_copy,
    read_json,
    records_path,
    report_progress,
    validate_radiometry,
    write_json,
)

expect_exe = system() == "Windows"
# TODO Slurm, threading, and https://github.com/pydata/xarray/issues/7549
dask.config.set(scheduler="synchronous")


def groupby(dataset: xr.Dataset, groups: xr.DataArray) -> Iterable[tuple]:
    """Return subsets selected from `dataset` at the coordinates of each unique
    value present in `groups`.
    """
    if groups.dims:
        # yield from the iterable created by xr.DataArray.groupby
        for key, value in groups.groupby(groups):
            subset = dataset.sel(value.unstack().indexes)
            yield key, subset
    else:
        # compensate for xr.DataArray.groupby's inability to handle no dims
        yield groups.item(), dataset


def split_list_arg(arg: str) -> dict:
    """Split the comma separated elements of `arg`, which are key, value
    pairs for the returned dictionary, parsing the value for later use
    in xr.Dataset indexing.
    """
    split = {}
    for item in arg.split(","):
        key, *value = item.split(":")
        if value:
            start, *stop = value[0].split("-")
            start = int(start)
            if stop:
                stop = int(stop[0])
            else:
                stop = start
            stop = stop + 1
        else:
            start = stop = None
        split[key] = slice(start, stop)
    return split


def reduce_by_coords(
    coords: dict,
    inputs: xr.Dataset,
    outdirs: Path,
) -> tuple[xr.Dataset, xr.DataArray]:
    """Select inputs and build paths along split dimensions using original indices."""
    selected = inputs.isel(coords)
    paths = xr.DataArray(outdirs)
    for dim, selection in coords.items():
        positions = np.arange(inputs.sizes[dim])[selection]
        subdirs = xr.DataArray(
            [Path(dim) / str(i) for i in positions],
            dims=dim,
            coords={dim: selected[dim]} if dim in selected.coords else None,
        )
        paths = paths / subdirs
    return selected, paths


class Program:
    def __init__(
        self,
        defaults: xr.Dataset,
        params: tuple = (),
    ) -> None:
        self.program = self.program + (".exe" if expect_exe else "")
        self.defaults = defaults
        # initialize the params @property
        self._params = tuple(i.__name__ for i in params)

    def __call__(self, args: Namespace) -> None:
        # with the `--pre` argument, write inputs and return
        if args.pre or not args.inputs.exists():
            # build coordinates dataset and write it to netCDF
            path = args.inputs
            path.parent.mkdir(parents=True, exist_ok=True)
            self.defaults.to_netcdf(path=path)
            if args.pre:
                return
        # read existing inputs
        inputs = xr.load_dataset(args.inputs)
        # begin construction of output path(s) relative to args.outputs.parent
        outdirs = xr.DataArray(Path())
        prefix = Path(args.outputs.stem)
        # with the `--split` argument, prepare to process a subset of inputs
        # in subdirectories, defined by the dimension(s) given as the split
        if args.split:
            coordinates = split_list_arg(args.split)
            inputs, outdirs = reduce_by_coords(coordinates, inputs, prefix)
            prefix = Path()
        identity = Identity()
        selections = list(groupby(dataset=inputs, groups=outdirs))
        identity.add_output(args.outputs, tree=args.outputs.parent / args.outputs.stem)
        for key, _ in selections:
            identity.add_output(args.outputs.parent / key / args.outputs.name)
        dispatch = dict(
            records=[
                dict(
                    name=str(key),
                    output=str(args.outputs.parent / key / args.outputs.name),
                    status="unattempted",
                    execution=dict(
                        records=[
                            dict(name=str(i), status="unattempted")
                            for i, _ in enumerate(
                                [None]
                                if args.combine and not args.dry
                                else self.execution_inputs(value)
                            )
                        ]
                    ),
                )
                for key, value in selections
            ]
        )
        dispatch_file = records_path(args.outputs) / f"dispatch-{uuid4().hex}.json"
        if not args.dry:
            write_json(dispatch_file, dispatch)
        reports = []
        for record, (key, value) in zip(dispatch["records"], selections):
            path = args.outputs.parent / key / args.outputs.name
            try:
                if args.combine and not args.dry:
                    self._combine(value, key, prefix, args, identity)
                    result = dict(
                        complete=True,
                        records=[dict(name=str(path), status="successful")],
                    )
                else:
                    result = self.run(
                        value,
                        path,
                        args.tmp,
                        args.dry,
                        resume=getattr(args, "resume", False),
                        retry_failed=getattr(args, "retry_failed", False),
                        identity=identity,
                    )
            except KeyboardInterrupt as error:
                record.update(
                    status="cancelled",
                    reason="User cancellation",
                    execution=getattr(error, "report", record["execution"]),
                )
                write_json(dispatch_file, dispatch)
                raise
            except (OSError, ValueError, NotImplementedError) as error:
                if args.dry:
                    raise
                result = dict(
                    complete=False,
                    records=[
                        dict(
                            name=str(path),
                            status=(
                                "blocked"
                                if isinstance(error, (OSError, ExecutionFailure))
                                else "scientifically_failed"
                            ),
                            reason=f"{type(error).__name__}: {error}",
                        )
                    ],
                )
            if not args.dry:
                reports.append(result)
                record.update(
                    status=(
                        "successful"
                        if result["complete"]
                        else next(
                            (
                                r["status"]
                                for r in result["records"]
                                if r["status"] != "successful"
                            ),
                            "retryable",
                        )
                    ),
                    execution=result,
                )
                if result["complete"]:
                    report_progress(dispatch_file, dispatch)
                else:
                    write_json(dispatch_file, dispatch)
        if not args.dry:
            records = [r for report in reports for r in report["records"]]
            if any(not r["complete"] for r in reports):
                raise ExecutionFailure(
                    dict(
                        records=records,
                        reason="; ".join(r.get("reason", "") for r in reports),
                    )
                )

    def _combine(self, value, key, prefix, args, identity):
        """Publish one complete selected destination; siblings remain independent."""
        coordinates = split_list_arg(args.combine)
        selected, combine_dirs = reduce_by_coords(coordinates, value, key)
        destination = args.outputs.parent / key / args.outputs.name
        with output_lock(records_path(destination)):
            archive_output(destination)
            paths = []
            for child, subset in groupby(selected, combine_dirs):
                path = args.outputs.parent / prefix / child / args.outputs.name
                binding, _ = identity.inputs(self.prepare_inputs(subset), self.program)
                accepted = records_path(path) / binding / "combined.nc"
                if (
                    not path.is_file()
                    or not accepted.is_file()
                    or identity.path(path) != identity.path(accepted)
                ):
                    raise ExecutionFailure(
                        dict(
                            records=[
                                dict(
                                    name=str(path),
                                    status="unattempted",
                                    reason="Selected split product is missing, incomplete or belongs to another request/runtime/data identity",
                                )
                            ]
                        )
                    )
                paths.append(path)
            with xr.open_mfdataset(
                paths,
                data_vars="different",
                compat="no_conflicts",
                join="outer",
            ) as dataset:
                try:
                    dataset = self.post(dataset)
                except NotImplementedError:
                    pass
                validate_radiometry(dataset)
                binding, details = identity.inputs(
                    self.prepare_inputs(value), self.program
                )
                directory = records_path(destination) / binding
                directory.mkdir(exist_ok=True)
                write_json(
                    directory / "combine.json",
                    dict(identity=details, sources=[str(p) for p in paths]),
                )
                atomic_write(directory / "combined.nc", lambda p: dataset.to_netcdf(p))
                publish_copy(directory / "combined.nc", destination)

    def run(
        self,
        inputs: xr.Dataset,
        outputs: Path,
        tmp: Path,
        dry: bool,
        *,
        resume=False,
        retry_failed=False,
        identity=None,
        request=None,
        finalize=None,
    ) -> dict | None:
        """Publish each native execution group before starting the next one.

        Complete group files are the recovery authority. JSON summaries can be
        rebuilt from them after reporting failures or interruption. Native work
        uses durable attempt directories; ``tmp`` remains the dry-run workspace.
        """
        inputs = self.prepare_inputs(inputs)
        outputs = Path(outputs).expanduser().resolve()
        outputs.parent.mkdir(parents=True, exist_ok=True)
        if dry:
            with TemporaryDirectory(dir=tmp) as work:
                for one_input in self.execution_inputs(inputs):
                    infile, _ = self.infile(Path(work), one_input.squeeze())
                    copy(Path(work) / infile, outputs.parent)
            return None
        identity = identity or Identity()
        identity.add_output(outputs)
        with output_lock(records_path(outputs)):
            key, binding = identity.inputs(inputs, self.program, request)
            directory = records_path(outputs) / key
            directory.mkdir(exist_ok=True)
            groups = list(self.execution_inputs(inputs))
            report = dict(
                request_id=key,
                output=str(outputs),
                records=[
                    dict(name=str(i), status="unattempted") for i in range(len(groups))
                ],
            )
            write_json(directory / "plan.json", dict(identity=binding, **report))
            archive_output(outputs)
            products = []
            blocked = None
            for index, one_input in enumerate(groups):
                record = report["records"][index]
                record_dir = directory / str(index)
                record_dir.mkdir(exist_ok=True)
                state_file = record_dir / "state.json"
                product = record_dir / "product.nc"
                previous = read_json(state_file, {})
                if blocked is not None:
                    record.update(status="blocked", reason=blocked)
                    write_json(state_file, record)
                    continue
                # A completed atomic product survives missing/stale optional status.
                if resume and product.exists():
                    try:
                        receipt = read_json(record_dir / "product.json", {})
                        if receipt.get("sha256") != identity.path(product):
                            raise ValueError(
                                "Completed product hash differs from publication receipt"
                            )
                        with xr.open_dataset(product) as dataset:
                            self.validate_output(dataset, one_input)
                    except (OSError, ValueError) as error:
                        previous = dict(
                            status="retryable", reason=f"Stored product: {error}"
                        )
                    else:
                        record.update(
                            status="successful", product=str(product), reused=True
                        )
                        products.append(product)
                        continue
                if (
                    resume
                    and previous.get("status") == "scientifically_failed"
                    and not retry_failed
                ):
                    record.update(previous)
                    continue
                # Never overwrite a completed or interrupted attempt's evidence.
                attempt = record_dir / uuid4().hex
                attempt.mkdir()
                if product.exists():
                    product.replace(attempt / "previous-product.nc")
                record.update(
                    status="retryable",
                    reason="Interrupted before publication",
                    attempt=str(attempt),
                )
                write_json(state_file, record)
                stage = "input"
                try:
                    infile, outfile = self.infile(attempt, one_input.squeeze())
                    copy(attempt / infile, outputs.parent)
                    with (attempt / "native.log").open("wb") as log:
                        stage = "native"
                        subprocess.run(
                            args=[self.program, str(attempt / infile)],
                            check=True,
                            stdout=log,
                            stderr=subprocess.STDOUT,
                        )
                    stage = "read"
                    mie_only = (
                        "MIE_TABLE_CAL" in one_input
                        and one_input.MIE_TABLE_CAL.item() == 1
                    )
                    if not mie_only:
                        with self.read_output(attempt / outfile, one_input) as native:
                            one_output = self._merge_output(one_input, native.load())
                        self.validate_output(one_output, one_input)
                        one_output = self.result_dataset(one_output)
                        stage = "publish"

                        def write(path):
                            one_output.to_netcdf(path)
                            if finalize is not None:
                                finalize(path)
                            write_json(
                                record_dir / "product.json",
                                dict(sha256=identity.path(path)),
                            )

                        atomic_write(product, write)
                        products.append(product)
                except KeyboardInterrupt as error:
                    record.update(status="cancelled", reason="User cancellation")
                    write_json(attempt / "state.json", record)
                    write_json(state_file, record)
                    report_progress(records_path(outputs) / "status.json", report)
                    error.report = report
                    raise
                except Exception as error:
                    if isinstance(error, subprocess.CalledProcessError):
                        status = (
                            "retryable"
                            if error.returncode < 0
                            else "scientifically_failed"
                        )
                    elif stage == "native" and isinstance(
                        error, (FileNotFoundError, PermissionError)
                    ):
                        status = "blocked"
                    elif stage == "publish" or isinstance(
                        error, (FileNotFoundError, PermissionError)
                    ):
                        status = "retryable"
                    elif stage == "input" and isinstance(error, OSError):
                        status = "retryable"
                    else:
                        status = "scientifically_failed"
                    record.update(
                        status=status, reason=f"{type(error).__name__}: {error}"
                    )
                    write_json(attempt / "state.json", record)
                    write_json(state_file, record)
                    if status == "blocked":
                        blocked = record["reason"]
                else:
                    record.update(
                        status="successful",
                        reason=(
                            "Native Mie-table-only execution; no radiometric product"
                            if mie_only
                            else ""
                        ),
                        product=str(product) if not mie_only else None,
                    )
                    report_progress(state_file, record)
                report_progress(records_path(outputs) / "status.json", report)
            report["complete"] = all(
                r["status"] == "successful" for r in report["records"]
            )
            if report["complete"] and products:
                try:
                    if len(products) == 1:
                        publish_copy(products[0], directory / "combined.nc")
                        publish_copy(directory / "combined.nc", outputs)
                    else:
                        datasets = [xr.load_dataset(p) for p in products]
                        ds = xr.combine_by_coords(
                            datasets,
                            data_vars="different",
                            coords="different",
                            compat="equals",
                        )
                        atomic_write(
                            directory / "combined.nc", lambda p: ds.to_netcdf(p)
                        )
                        publish_copy(directory / "combined.nc", outputs)
                except (OSError, ValueError) as error:
                    report.update(
                        complete=False, reason=f"Combined publication: {error}"
                    )
            report_progress(records_path(outputs) / "status.json", report)
            return report

    def _merge_output(self, one_input, one_output):
        if "phony_dim_0" in one_output.dims:
            return one_output
        one_input = self.result_inputs(one_input, one_output)
        for item in one_input.coords:
            if item not in one_output:
                continue
            if (one_input[item] == one_output[item]).all():
                one_output = one_output.drop_vars(item)
            elif (
                np.issubdtype(one_input[item].dtype, np.number)
                and np.isnan(one_input[item]).all()
            ):
                one_input = one_input.drop_vars(item)
            else:
                raise ValueError(f"Inputs/outputs differ for {item}")
        return xr.merge((one_input, one_output))

    def validate_output(self, dataset, inputs):
        validate_radiometry(dataset)

    def infile(self, path: Path, dataset: xr.Dataset) -> tuple[Path, Path]:
        """Write parameters to a text file, and return its path."""
        dataset = self.prepare_inputs(dataset)
        self.params = dataset
        lines = []
        for item in self.params:
            if item not in dataset:
                continue
            param = dataset[item]
            name = param.attrs.get("long_name", item)
            value = param.item()
            lines.append(f"{value:<24} # {name}")
        outfile = tokenize(dataset, ensure_deterministic=True)
        outfile = Path(outfile).with_suffix(".outfile")
        lines += [str(path / outfile), ""]
        infile = outfile.with_suffix(".infile.txt")
        with (path / infile).open("w") as stream:
            stream.write("\n".join(lines))
        return infile, outfile.with_suffix(".outfile.h5")

    def post(self, dataset: xr.Dataset) -> xr.Dataset:
        raise NotImplementedError

    def execution_inputs(self, inputs: xr.Dataset) -> Iterable[xr.Dataset]:
        """Native execution groups; by default one scalar parameter combination."""
        shape = tuple(inputs.sizes.values())
        each_input = xr.DataArray(
            # See https://github.com/pydata/xarray/issues/7423 for the float dtype.
            data=np.arange(np.prod(shape), dtype=float).reshape(shape),
            coords=inputs.coords,
        )
        for _, one_input in groupby(dataset=inputs, groups=each_input):
            yield one_input

    def read_output(self, path: Path, inputs: xr.Dataset) -> xr.Dataset:
        return xr.open_dataset(path)

    def result_dataset(self, dataset: xr.Dataset) -> xr.Dataset:
        """Form the scientific product before publishing the combined output."""
        return dataset

    def prepare_inputs(self, dataset: xr.Dataset) -> xr.Dataset:
        """Prepare coordinates for native serialization and output metadata."""
        return dataset

    def result_inputs(self, inputs: xr.Dataset, outputs: xr.Dataset) -> xr.Dataset:
        """Return input coordinates to check and merge with a native result."""
        return inputs

    @property
    def params(self) -> tuple:
        return self._params

    @params.setter
    def params(self, dataset: xr.Dataset) -> None:
        pass
