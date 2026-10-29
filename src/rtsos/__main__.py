from argparse import ArgumentParser
from importlib import import_module
from pathlib import Path
from tempfile import gettempdir

parser = ArgumentParser(
    description=(
        "Run a program that uses RTSOS to solve the multiple scattering radiative"
        " transfer equation from UV, visible, to infrared."
    ),
)
common = ArgumentParser(add_help=False)
common.add_argument(
    "--resume", action="store_true", help="reuse completed matching native groups"
)
common.add_argument(
    "--retry-failed",
    action="store_true",
    help="with --resume, explicitly retry scientific failures as well",
)
common.add_argument(
    "--dry-run",
    default=False,
    dest="dry",
    action="store_true",
    help="only generate parameter files",
)
common.add_argument(
    "--tmp",
    type=Path,
    default=gettempdir(),
    help="path where temporary directories will be created",
)
common.add_argument(
    "--pre",
    action="store_true",
    help="write default parameters to the location given as `inputs`",
)
common.add_argument(
    "--split",
    type=str,
    help=(
        "comma-separated list of dimensions by which inputs are split into "
        "nested subdirectories, with optional positions or inclusive ranges as "
        "`dim:index` or `dim:start-stop` respectively"
    ),
)
common.add_argument(
    "--combine",
    type=str,
    help=(
        "dimensions by which outputs from a previous call using "
        "`--split` are combined, formatted as for `--split`"
    ),
)
common.add_argument(
    "inputs",
    type=Path,
    help=(
        "path for input file, a NetCDF whose coordinates are parameters "
        "over which the model will generate outputs"
    ),
)
common.add_argument(
    "outputs",
    nargs="?",
    type=Path,
    help=("path for output file(s) (ignored, if given, with `--pre`)"),
)
subparsers = parser.add_subparsers(
    dest="program",
    title="program (required)",
    required=True,
    help="use `rtsos <program> --help` for more information",
)
subparsers.add_parser(
    name="GSFC_AC_LUT",
    parents=[common],
    description=(
        "Generate look-up-tables (LUTs) used in atmospheric correction for"
        " top-of-atmosphere (i.e. satellite based) measurements of the Earth's surface."
    ),
)
pace_parser = subparsers.add_parser(
    name="PACE_Simulator_DoubleK",
    parents=[common],
    description=(
        "The PACE simulator runs the monochromatic RTSOS, with a  a list of built-in"
        " aerosol and ocean inherent optical properties."
    ),
)
pace_parser.add_argument(
    "--data-root",
    type=Path,
    help="directory containing Data and Gas_Absorption_Coefficients (overrides input coordinates)",
)
subparsers.add_parser(
    name="Aerosol_Phmx_Cal",
    parents=[common],
    description=(
        "Generate a file containing the aerosol scattering matrix to be read by RTSOS."
    ),
)

simulate_parser = subparsers.add_parser(
    "simulate",
    description="Simulate one physical case with PACE, monochromatic or full-response integration.",
)
simulate_parser.add_argument(
    "case", type=Path, help="TOML case file; data paths are relative to its directory"
)
simulate_parser.add_argument("outputs", type=Path, help="output NetCDF path")
simulate_parser.add_argument(
    "--tmp",
    type=Path,
    default=Path(gettempdir()),
    help="parent directory for native temporary files",
)
simulate_parser.add_argument(
    "--dry-run",
    action="store_true",
    help="write native input records without execution",
)
simulate_parser.add_argument("--resume", action="store_true")
simulate_parser.add_argument("--retry-failed", action="store_true")

cases_parser = subparsers.add_parser(
    "simulate-cases",
    description="Run an ordered list of named case files independently.",
)
cases_parser.add_argument(
    "cases", type=Path, help="TOML file containing [[cases]] name/file entries"
)
cases_parser.add_argument(
    "outputs", type=Path, help="directory for name.nc products and accounting"
)
cases_parser.add_argument("--resume", action="store_true")
cases_parser.add_argument("--retry-failed", action="append", default=[], metavar="NAME")
cases_parser.add_argument("--tmp", type=Path, default=Path(gettempdir()))
cases_parser.add_argument("--dry-run", action="store_true")


def main(argv: str | None = None) -> None:
    args = parser.parse_args(argv)
    from .execution import ExecutionFailure

    if args.program == "simulate-cases":
        from .simulation import simulate_cases

        try:
            report = simulate_cases(
                args.cases,
                args.outputs,
                resume=args.resume,
                retry_failed=args.retry_failed,
                tmp=args.tmp,
                dry_run=args.dry_run,
            )
        except (ValueError, OSError) as error:
            cases_parser.error(str(error))
        expected_status = "unattempted" if args.dry_run else "successful"
        if any(r["status"] != expected_status for r in report["records"]):
            cases_parser.exit(1, str(ExecutionFailure(report)) + "\n")
        return
    if args.program == "simulate":
        from .simulation import simulate

        try:
            simulate(
                args.case,
                args.outputs,
                tmp=args.tmp,
                dry_run=args.dry_run,
                resume=args.resume,
                retry_failed=args.retry_failed,
            )
        except (ValueError, NotImplementedError, OSError) as error:
            simulate_parser.error(str(error))
        return
    module = import_module(f"{__package__}.{args.program}")
    program = (
        module.init(data_root=args.data_root)
        if args.program == "PACE_Simulator_DoubleK"
        else module.init()
    )
    try:
        program(args)
    except ExecutionFailure as error:
        parser.exit(1, str(error) + "\n")


if __name__ == "__main__":
    main()
