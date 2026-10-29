"""Verify a non-editable installation using tests/assets from an unpacked sdist.

Run with that installation's Python, outside the checkout, with PYTHONPATH unset:
    python -I /path/to/sdist/validation/check_installed.py /path/to/sdist /path/to/work
"""

from argparse import ArgumentParser
from importlib import import_module, metadata
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import sysconfig


def check_identity(source):
    distribution = metadata.distribution("rtsos")
    files = {
        str(path): distribution.locate_file(path).resolve()
        for path in distribution.files
    }
    modules = {}
    for path in (source / "src/rtsos").glob("*.py"):
        name = "rtsos" if path.stem == "__init__" else f"rtsos.{path.stem}"
        installed = Path(import_module(name).__file__).resolve()
        assert installed == files[f"rtsos/{path.name}"], (name, installed)
        assert not installed.is_relative_to(source), installed
        modules[name] = str(installed)
    from rtsos._version import __version__

    assert __version__ == distribution.version
    assert not any(
        name == "zhai_rt" or name.startswith("zhai_rt.") for name in sys.modules
    )
    scripts = Path(sysconfig.get_path("scripts"))
    commands = {}
    for name in (
        "rtsos",
        "rtsos_Monochromatic",
        "rtsos_PACE_Simulator_DoubleK",
        "rtsos_GSFC_AC_LUT",
        "rtsos_Aerosol_Phmx_Cal",
        "rtsos_rao_dg",
    ):
        installed = (scripts / name).resolve()
        resolved = shutil.which(name)
        assert resolved is not None and Path(resolved).resolve() == installed, (
            name,
            resolved,
        )
        assert installed in files.values(), installed
        commands[name] = str(installed)
    for name in ("PACE_Simulator_DoubleK", "GSFC_AC_LUT", "Aerosol_Phmx_Cal"):
        program = import_module(f"rtsos.{name}").init()
        assert program.program == f"rtsos_{name}"
    from rtsos.monochromatic import Monochromatic

    assert Monochromatic().program == "rtsos_Monochromatic"
    print(
        json.dumps(
            dict(
                python=sys.executable,
                version=distribution.version,
                modules=modules,
                commands=commands,
            ),
            indent=2,
        ),
        flush=True,
    )
    return commands


def smoke(commands, source, work):
    def run(*args):
        print("+", *map(str, args), flush=True)
        result = subprocess.run(
            args,
            cwd=work,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        if result.returncode:
            print(result.stdout, flush=True)
        result.check_returncode()
        return result.stdout

    for name, command in commands.items():
        output = run(command, "--help") if name == "rtsos" else run(command)
        assert ("--help" if name == "rtsos" else "TOO FEW INPUTS") in output, output
        print(output, flush=True)

    import xarray as xr

    for program in ("PACE_Simulator_DoubleK", "GSFC_AC_LUT", "Aerosol_Phmx_Cal"):
        assert "--help" in run(commands["rtsos"], program, "--help")
        folder = work / program
        folder.mkdir()
        inputs = folder / "inputs.nc"
        run(commands["rtsos"], program, "--pre", inputs)
        dataset = xr.load_dataset(inputs)
        dataset.isel({dim: slice(0, 1) for dim in dataset.dims}).to_netcdf(inputs)
        run(
            commands["rtsos"],
            program,
            inputs,
            folder / "outputs.nc",
            "--dry-run",
            "--tmp",
            work,
        )
        assert len(list(folder.glob("*.infile.txt"))) == 1
        assert not (folder / "outputs.nc").exists()

    for example in ("single_direction", "response_directions", "pace_ocean"):
        folder = work / example
        folder.mkdir()
        run(
            commands["rtsos"],
            "simulate",
            source / f"examples/{example}.toml",
            folder / "outputs.nc",
            "--dry-run",
            "--tmp",
            work,
        )
        assert list(folder.glob("*.infile.txt"))
        assert not (folder / "outputs.nc").exists()
    run(commands["rtsos"], "simulate-cases", "--help")
    run(
        commands["rtsos"],
        "simulate-cases",
        source / "examples/cases.toml",
        work / "cases",
        "--dry-run",
        "--tmp",
        work,
    )
    assert list((work / "cases").rglob("*.infile.txt"))


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="unpacked sdist root")
    parser.add_argument("work", type=Path, help="new directory for smoke/test outputs")
    args = parser.parse_args()
    source, work = args.source.resolve(), args.work.resolve()
    assert "PYTHONPATH" not in os.environ, "Unset PYTHONPATH for installed checks"
    assert not (source / ".git").exists(), "Use the sdist's tests and assets"
    work.mkdir(parents=True)
    commands = check_identity(source)
    smoke(commands, source, work)
    environment = dict(os.environ, PYTEST_DISABLE_PLUGIN_AUTOLOAD="1")
    subprocess.run(
        [
            sys.executable,
            "-I",
            "-m",
            "pytest",
            "--import-mode=importlib",
            "-c",
            str(source / "pyproject.toml"),
            str(source / "src/rtsos/tests"),
            "-q",
            "-o",
            f"cache_dir={work / 'pytest-cache'}",
            f"--basetemp={work / 'pytest-tmp'}",
        ],
        cwd=work,
        env=environment,
        check=True,
    )


if __name__ == "__main__":
    main()
