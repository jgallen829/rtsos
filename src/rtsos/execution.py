"""Ordinary durable records shared by case and legacy native execution."""

from contextlib import contextmanager
from hashlib import sha256
import json
import os
from pathlib import Path
from shutil import copyfile, which
import sys
from uuid import uuid4

from dask.tokenize import tokenize
import numpy as np
import xarray as xr


class ExecutionFailure(ValueError):
    """A requested population has records without accepted products."""

    def __init__(self, report):
        self.report = report
        reasons = "; ".join(
            f"{r['name']}: {r['status']}: {r.get('reason', '')}"
            for r in report["records"]
            if r["status"] != "successful"
        )
        super().__init__(reasons or report.get("reason", "Incomplete execution"))


def records_path(output):
    return output.with_name(output.name + ".records")


def read_json(path, default=None):
    if not path.exists():
        return default
    return json.loads(path.read_text())


def atomic_write(path, write):
    """Only the rename publishes a finished file; failed writes stay as evidence."""
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(f".{path.stem}.{uuid4().hex}.partial{path.suffix}")
    write(partial)
    with partial.open("rb") as stream:
        os.fsync(stream.fileno())
    partial.replace(path)
    if os.name != "nt":
        descriptor = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(descriptor)
        finally:
            os.close(descriptor)


def write_json(path, value):
    atomic_write(path, lambda p: p.write_text(json.dumps(value, indent=2) + "\n"))


def report_progress(path, value):
    """The summary is convenient; immutable products and record states are primary."""
    try:
        write_json(path, value)
    except OSError as error:
        # Reporting after publication must not change execution eligibility, even
        # when a caller configures warnings as errors.
        try:
            print(f"RTSOS: could not update {path}: {error}", file=sys.stderr)
        except OSError:
            pass  # An unavailable reporting stream cannot invalidate a product.


def archive_output(output):
    """An earlier public product cannot masquerade as a failed replacement."""
    if output.exists():
        history = records_path(output) / "history" / uuid4().hex
        history.mkdir(parents=True)
        output.replace(history / output.name)


@contextmanager
def output_lock(directory):
    """OS locks release on interruption; one output never has competing writers."""
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / "writer.lock").open("a+b") as stream:
        if os.name == "nt":
            import msvcrt

            stream.write(b"0")
            stream.flush()
            stream.seek(0)
            msvcrt.locking(stream.fileno(), msvcrt.LK_NBLCK, 1)
        else:
            import fcntl

            fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        try:
            yield
        finally:
            if os.name == "nt":
                stream.seek(0)
                msvcrt.locking(stream.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                fcntl.flock(stream, fcntl.LOCK_UN)


class Identity:
    """Authenticate immutable data/runtime once per invocation, shared by cases.

    Directory identities include names and file contents, following data links.
    Declared output files/sidecars are excluded when sharing an input directory.
    Git revisions, version labels and documentation do not control reuse.
    """

    def __init__(self, *, ignored=()):
        self.paths = {}
        self.outputs = set()
        self.output_trees = {}
        self.ignored = {Path(p).resolve() for p in ignored}

    def add_output(self, output, *, tree=None):
        output = Path(output).expanduser().resolve()
        self.outputs.add(output)
        self.ignored.update((output, records_path(output)))
        if tree is not None:
            self.output_trees[Path(tree).resolve()] = output.name

    def _generated(self, path):
        if path.resolve() in self.ignored or path.name.endswith(".infile.txt"):
            return True
        for root, name in self.output_trees.items():
            if path.resolve().is_relative_to(root) and (
                path.name in (name, name + ".records")
                or (
                    path.name.startswith(f".{Path(name).stem}.")
                    and path.name.endswith(f".partial{Path(name).suffix}")
                )
            ):
                # Standard legacy split products may be nested or staged. Only
                # their declared filenames/sidecars are artifacts, not the tree.
                return True
        return any(
            path.parent == p.parent
            and path.name.startswith(f".{p.stem}.")
            and path.name.endswith(f".partial{p.suffix}")
            for p in self.outputs
        )

    def path(self, path):
        path = Path(path).expanduser().resolve()
        if path not in self.paths:
            if path.is_dir():
                entries = []
                for p in sorted(path.iterdir()):
                    if self._generated(p):
                        continue
                    fingerprint = self.path(p)
                    # Empty/generated-only parents carry no input file bytes.
                    # Creating split-output directories must not change data identity.
                    if p.is_dir() and fingerprint == tokenize([]):
                        continue
                    entries.append((p.name, fingerprint))
                value = tokenize(entries)
            elif path.is_file():
                with path.open("rb") as stream:
                    digest = sha256()
                    while chunk := stream.read(1024 * 1024):
                        digest.update(chunk)
                value = digest.hexdigest()
            else:
                value = "missing"
            self.paths[path] = value
        return self.paths[path]

    def runtime(self, program):
        package = Path(__file__).parent
        modules = (
            "core",
            "execution",
            "simulation",
            "parameters",
            "physical",
            "PACE_Simulator_DoubleK",
            "monochromatic",
            "bands",
            "responses",
        )
        specific = program.removeprefix("rtsos_").removesuffix(".exe")
        sources = {
            n: self.path(package / f"{n}.py")
            for n in (*modules, specific)
            if (package / f"{n}.py").is_file()
        }
        exe = which(program)
        return dict(
            program=program,
            executable=str(Path(exe).resolve()) if exe else None,
            executable_sha256=self.path(exe) if exe else "missing",
            sources=sources,
            python=sys.version,
            numpy=np.__version__,
            xarray=xr.__version__,
        )

    def inputs(self, inputs, program, extra=None):
        paths = {}
        for name in (
            "Aux_Dir",
            "Atmos_Dir",
            "Mie_Database_Dir",
            "Aerosol_Phasematrix_File_Hi",
            "Aerosol_Phasematrix_File_Low",
            "land_ref_spectra_filename",
            "CFILE_INSTRUMENT",
        ):
            if (
                name == "Mie_Database_Dir"
                and "MIE_TABLE_CAL" in inputs
                and (inputs.MIE_TABLE_CAL != 3).all()
            ):
                continue  # mode 1 writes this database; mode 2 does not consume it
            if name in inputs:
                for value in np.asarray(inputs[name]).flat:
                    if str(value):
                        path = Path(str(value)).expanduser().resolve()
                        paths[str(path)] = self.path(path)
        # Coordinate documentation does not alter a native request. Dataset attrs
        # include actual solar/phase/geometry policies and are part of its identity.
        request = [
            (n, v.dims, v.dtype.str, v.values) for n, v in inputs.variables.items()
        ]
        identity = dict(
            runtime=self.runtime(program),
            data=paths,
            request=tokenize(request, inputs.attrs, extra),
        )
        return tokenize(identity), identity


def validate_radiometry(dataset):
    """Narrow acceptance: finite measured fields, including valid signed zeros."""
    if not dataset.data_vars:
        raise ValueError("Native output contains no scientific variables")
    for name, variable in dataset.data_vars.items():
        field = name.removeprefix("Node_")
        if field.startswith(
            ("Stokes_Vector", "Radiance", "Irrad", "Solar_Irradiance", "Q_", "U_", "V_")
        ):
            if not variable.size or not np.all(np.isfinite(variable)):
                raise ValueError(f"nonfinite native {name} or empty radiometric field")


def publish_copy(source, destination):
    atomic_write(destination, lambda p: copyfile(source, p))
