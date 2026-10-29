"""Readable physical cases and the shared Python/CLI simulation entry point."""

from copy import deepcopy
from dataclasses import asdict, dataclass
from os import PathLike
from pathlib import Path
from tempfile import gettempdir
import tomllib
from time import time_ns

from netCDF4 import Dataset
import xarray as xr

from .execution import (
    ExecutionFailure,
    Identity,
    archive_output,
    output_lock,
    read_json,
    records_path,
    report_progress,
    write_json,
)


@dataclass
class Case:
    """Plain case tables; direct paths are absolute after construction.

    Use ``from_mapping`` for Python cases and ``load_case`` for TOML. Relative
    paths use the supplied base directory, or the case file's directory. Loading
    represents a request; the adapter separately checks execution capabilities.
    """

    sensor: dict
    scene: dict
    geometry: dict
    physics: dict
    numerics: dict
    illumination: dict
    spectral: dict
    data: dict

    @classmethod
    def from_mapping(cls, tables: dict, *, base_dir: str | Path) -> "Case":
        tables = deepcopy(tables)
        expected = cls.__dataclass_fields__.keys()
        if tables.keys() != expected or any(
            not isinstance(value, dict) for value in tables.values()
        ):
            raise ValueError(f"A case needs these tables: {', '.join(expected)}")
        for name, table in tables["scene"].items():
            if not isinstance(table, dict):
                raise ValueError(f"scene.{name} must be a table")
        base = Path(base_dir).expanduser().resolve()
        path_fields = [
            (tables["data"], name) for name in ("root", "auxiliary", "gas_absorption")
        ]
        path_fields += [(tables["sensor"], "response_file")]
        path_fields += [(tables["scene"].get("surface", {}), "spectrum")]
        path_fields += [
            (tables["scene"].get("aerosol", {}), name)
            for name in ("phase_matrix_hi", "phase_matrix_low")
        ]
        for table, key in path_fields:
            if key in table:
                if not isinstance(table[key], (str, PathLike)):
                    raise ValueError(f"{key} must be a filesystem path")
                path = Path(table[key]).expanduser()
                table[key] = str((base / path).resolve())
        return cls(**tables)

    def load_responses(self):
        """Load/select complete response channels independently of execution."""
        from .responses import load_oci_rsr, load_responses

        path = self.sensor["response_file"]
        fmt = self.sensor.get("response_format", "toml")
        if fmt == "toml":
            response = load_responses(path)
        elif fmt == "oci":
            response = load_oci_rsr(path)
        else:
            raise ValueError(f"Unknown response_format {fmt!r}; use toml or oci")
        return response.select(
            self.sensor.get("channels"),
            detector=self.sensor.get("detector"),
            mirror=self.sensor.get("mirror"),
        )


def load_case(path: str | Path) -> Case:
    """Read TOML, resolving direct paths relative to this file, never to cwd."""
    path = Path(path).expanduser().resolve()
    with path.open("rb") as stream:
        return Case.from_mapping(tomllib.load(stream), base_dir=path.parent)


def simulate(
    case: Case | str | Path,
    outputs: str | Path,
    *,
    tmp: str | Path | None = None,
    dry_run: bool = False,
    resume: bool = False,
    retry_failed: bool = False,
    _identity: Identity | None = None,
) -> xr.Dataset | None:
    """Translate one case and execute through the existing Program mechanics.

    ``outputs`` and ``tmp`` are invocation paths relative to cwd. A dry run writes
    the same native input records beside ``outputs`` and returns None.
    """
    outputs = Path(outputs).expanduser().resolve()
    try:
        program, inputs, request, finalize = _prepare_simulation(case, dry_run)
    except (ValueError, NotImplementedError, OSError) as error:
        if not dry_run:
            with output_lock(records_path(outputs)):
                archive_output(outputs)
                report = dict(
                    records=[
                        dict(
                            name="case",
                            status=(
                                "blocked"
                                if isinstance(error, OSError)
                                else "scientifically_failed"
                            ),
                            reason=f"{type(error).__name__}: {error}",
                        )
                    ]
                )
                write_json(
                    records_path(outputs) / f"preparation-{time_ns()}.json", report
                )
                report_progress(records_path(outputs) / "status.json", report)
        raise
    report = program.run(
        inputs,
        outputs,
        Path(tmp or gettempdir()).expanduser().resolve(),
        dry_run,
        resume=resume,
        retry_failed=retry_failed,
        identity=_identity,
        request=request,
        finalize=finalize,
    )
    if dry_run:
        return None
    if not report["complete"]:
        raise ExecutionFailure(report)
    return xr.load_dataset(outputs)


def _prepare_simulation(case, dry_run):
    from .PACE_Simulator_DoubleK import case_inputs, init, solar_source

    if not isinstance(case, Case):
        case = load_case(case)
    method = case.spectral.get("method")
    responses = None
    if method in ("monochromatic", "response_integration"):
        from .monochromatic import Monochromatic, case_inputs, solar_source, check_data

        if method == "response_integration":
            from . import bands

            responses = case.load_responses()
            inputs = bands.case_inputs(case, responses)
            program = bands.ResponseIntegration(responses)
        else:
            inputs = case_inputs(case)
            program = Monochromatic()
    else:
        program = init()
        inputs = case_inputs(case)
    illumination = None if dry_run else solar_source(case)
    if method in ("monochromatic", "response_integration") and not dry_run:
        check_data(case, inputs, illumination)
    attributes = dict(
        spectral_method=method,
        illumination_policy=case.illumination["source"],
        illumination_source=str(illumination),
        sensor=case.sensor.get("name", "monochromatic"),
    )
    if method == "pace_native":
        attributes["channel_selection"] = case.sensor["channels"]

    def finalize(path):
        with Dataset(path, "a") as result:
            result.setncatts(attributes)
            if method in ("monochromatic", "response_integration"):
                _monochromatic_metadata(result)
            if responses is not None:
                bands.write_calibrations(result, responses)

    return program, inputs, asdict(case), finalize


def simulate_cases(
    cases: str | Path,
    outputs: str | Path,
    *,
    resume: bool = False,
    retry_failed=(),
    tmp: str | Path | None = None,
    dry_run: bool = False,
) -> dict:
    """Run an ordered TOML list of named existing case files, independently.

    ``[[cases]]`` entries contain only ``name`` and ``file``. File paths resolve
    relative to the list; each case retains its own path/physics conventions.
    Products are ``outputs/name.nc``. Return accounting for the full population.
    On resume, scientific failures require an explicit name in ``retry_failed``;
    operational blocks, interruptions and incomplete products are reconsidered.
    """
    from dask.tokenize import tokenize
    from .core import expect_exe

    source = Path(cases).expanduser().resolve()
    with source.open("rb") as stream:
        tables = tomllib.load(stream)
    entries = tables.get("cases", [])
    if set(tables) != {"cases"} or not isinstance(entries, list) or not entries:
        raise ValueError("A case list needs a nonempty [[cases]] array")
    names = []
    for entry in entries:
        if set(entry) != {"name", "file"}:
            raise ValueError("Each case-list entry needs only name and file")
        name = entry["name"]
        if (
            not isinstance(name, str)
            or not name
            or name in (".", "..")
            or any(c in name for c in "/\\")
        ):
            raise ValueError("Case names must be nonempty filename components")
        if name in names:
            raise ValueError(f"Duplicate case name {name!r}")
        names.append(name)
        if not isinstance(entry["file"], str):
            raise ValueError("Case file must be a path")
    unknown = set(retry_failed) - set(names)
    if unknown:
        raise ValueError(f"Unknown retry case names: {sorted(unknown)}")
    output_dir = Path(outputs).expanduser().resolve()
    records_dir = output_dir / "case-runs"
    identity = Identity(ignored=(records_dir,))
    identity.add_output(output_dir / "cases.json")
    for name in names:
        identity.add_output(output_dir / f"{name}.nc")
    with output_lock(records_dir):
        histories = sorted(records_dir.glob("*.json"))
        previous = read_json(histories[-1], {}) if resume and histories else {}
        previous = {r["name"]: r for r in previous.get("records", [])}
        report = dict(
            source=str(source),
            records=[
                dict(
                    name=e["name"],
                    file=str((source.parent / e["file"]).resolve()),
                    output=str(output_dir / f"{e['name']}.nc"),
                    status="unattempted",
                )
                for e in entries
            ],
        )
        run_file = records_dir / f"{time_ns()}.json"
        write_json(run_file, report)
        for record in report["records"]:
            output = Path(record["output"])
            try:
                case = load_case(record["file"])
                program = (
                    "rtsos_PACE_Simulator_DoubleK"
                    if case.spectral.get("method") == "pace_native"
                    else "rtsos_Monochromatic"
                )
                program += ".exe" if expect_exe else ""
                dependencies = dict(case.data)
                dependencies.update(
                    {k: v for k, v in case.sensor.items() if k == "response_file"}
                )
                for section in case.scene.values():
                    dependencies.update(
                        {
                            k: v
                            for k, v in section.items()
                            if k in ("spectrum", "phase_matrix_hi", "phase_matrix_low")
                        }
                    )
                binding = tokenize(
                    asdict(case),
                    identity.runtime(program),
                    {str(p): identity.path(p) for p in dependencies.values()},
                )
                record["request_id"] = binding
                old = previous.get(record["name"], {})
                if (
                    resume
                    and old.get("request_id") == binding
                    and old.get("status") == "scientifically_failed"
                    and "execution" not in old
                    and record["name"] not in retry_failed
                ):
                    record.update(status=old["status"], reason=old["reason"])
                else:
                    simulate(
                        case,
                        output,
                        tmp=tmp,
                        dry_run=dry_run,
                        resume=resume,
                        retry_failed=record["name"] in retry_failed,
                        _identity=identity,
                    )
                    record.update(
                        status="unattempted" if dry_run else "successful", reason=""
                    )
            except KeyboardInterrupt:
                record.update(status="cancelled", reason="User cancellation")
                write_json(run_file, report)
                report_progress(output_dir / "cases.json", report)
                raise
            except ExecutionFailure as error:
                states = {r["status"] for r in error.report["records"]}
                status = (
                    "blocked"
                    if "blocked" in states
                    else (
                        "scientifically_failed"
                        if "scientifically_failed" in states
                        else "retryable"
                    )
                )
                record.update(status=status, reason=str(error), execution=error.report)
            except OSError as error:
                record.update(
                    status="blocked", reason=f"{type(error).__name__}: {error}"
                )
            except (ValueError, NotImplementedError, KeyError) as error:
                record.update(
                    status="scientifically_failed",
                    reason=f"{type(error).__name__}: {error}",
                )
            if record["status"] in ("blocked", "scientifically_failed", "retryable"):
                try:
                    with output_lock(records_path(output)):
                        archive_output(output)
                except OSError as error:
                    record.update(
                        status="blocked",
                        reason=f"{record['reason']}; cannot retire public output: {error}",
                    )
                write_json(run_file, report)
            else:
                # A missing success summary is repaired from each atomic case product.
                report_progress(run_file, report)
            report_progress(output_dir / "cases.json", report)
        return report


def _monochromatic_metadata(result):
    result.setncatts(
        dict(
            stokes_convention="Native RTSOS I,Q,U,V; meridian-plane reference; no sign changes",
            angle_convention="ThetaV=acos(mu); 0 upward, 180 downward; PhiV native solar-plane azimuth",
            detector_convention="TOA top of atmosphere; AirSensor at sensor_height_km; BOA above interface; TOO below ocean interface",
            gas_method="Native cross-section searches at requested wavelength, then trapezoidal height integration; species cutoffs retained",
        )
    )
    for name, variable in result.variables.items():
        name = name.removeprefix("Node_")
        if name.startswith("Stokes_Vector"):
            variable.units = "W m-2 sr-1 um-1"
        elif name.startswith("Irradiance_") or name == "Solar_Irradiance":
            variable.units = "W m-2 um-1"
        elif name == "Wavelength_Evaluated":
            variable.units = "nm"
        elif name in (
            "ThetaV",
            "PhiV",
            "viewing_zenith_requested",
            "relative_azimuth_requested",
        ):
            variable.units = "degree"
        elif name.startswith("Tau_"):
            variable.units = "1"
    result.variables["N_Stokes_Parameter"].stokes_labels = "I Q U V"
    if "observation" in result.dimensions:
        result.setncatts(
            dict(
                geometry_mode="explicit_directions",
                angle_convention="Degrees; propagation zenith from upward normal (0 upward, 180 downward); solar propagation azimuth=0; PhiV evaluated modulo 360, no 0..180 folding",
                geometry_scope="One solar_zenith_deg (sun position) and sensor_height_km per physical case; ordered pairs at all detector levels",
            )
        )
        result.variables["ThetaV"].long_name = (
            "Native evaluated propagation zenith acos(mu)"
        )
        result.variables["PhiV"].long_name = (
            "Native evaluated propagation azimuth relative to solar propagation plane"
        )
