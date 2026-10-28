"""Readable physical cases and the shared Python/CLI simulation entry point."""

from copy import deepcopy
from dataclasses import dataclass
from os import PathLike
from pathlib import Path
from tempfile import gettempdir
import tomllib

from netCDF4 import Dataset
import xarray as xr


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
) -> xr.Dataset | None:
    """Translate one case and execute through the existing PACE Program.

    ``outputs`` and ``tmp`` are invocation paths relative to cwd. A dry run writes
    the same native input records beside ``outputs`` and returns None.
    """
    from .PACE_Simulator_DoubleK import case_inputs, init, solar_source

    if not isinstance(case, Case):
        case = load_case(case)
    inputs = case_inputs(case)
    illumination = None if dry_run else solar_source(case)
    outputs = Path(outputs).expanduser().resolve()
    init().run(
        inputs, outputs, Path(tmp or gettempdir()).expanduser().resolve(), dry_run
    )
    if dry_run:
        return None
    attributes = dict(
        spectral_method="pace_native",
        illumination_policy=case.illumination["source"],
        illumination_source=str(illumination),
        sensor="pace",
        channel_selection=case.sensor["channels"],
    )
    with Dataset(outputs, "a") as result:
        result.setncatts(attributes)
    return xr.load_dataset(outputs)
