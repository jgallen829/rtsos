"""Recovery at real native file/process and product publication boundaries.

The executable below is a small synthetic I/O producer, not an RT model. Existing
analytical/native tests own physics; these tests exercise ordering and durability.
"""

import json
import errno
import os
from pathlib import Path
import subprocess
import sys

from netCDF4 import Dataset
import numpy as np
import pytest
import xarray as xr

from rtsos import load_case, simulate, simulate_cases
from rtsos.__main__ import main
from rtsos.execution import ExecutionFailure, Identity, read_json, records_path
from rtsos.PACE_Simulator_DoubleK import init

EXAMPLES = Path(__file__).resolve().parents[3] / "examples"

SYNTHETIC = """
import json, os, sys
from pathlib import Path
import numpy as np
import xarray as xr
rows = Path(sys.argv[1]).read_text().splitlines()
mono = 'thuillier' in rows
values = [r.split('#')[0].strip() for r in rows]
key = int(float(values[1 if mono else 15]))
control = Path(os.environ['RTSOS_TEST_CONTROL'])
mode = json.loads(control.read_text()).get(str(key), '')
with control.with_suffix('.calls').open('a') as log:
    log.write(str(key) + '\\n')
if mono:
    p = rows.index('thuillier')
    wave = np.array([float(w) for w in rows[p+2:p+2+int(rows[p+1])]])
    dims = ['evaluation']
    coords = dict(evaluation=np.arange(len(wave)), N_Stokes_Parameter=[1,2,3,4])
    fields = dict(Wavelength_Evaluated=('evaluation', wave), Solar_Irradiance=('evaluation', wave+10))
    shape = [len(wave)]
    if 'explicit_directions' in rows:
        j = rows.index('explicit_directions')
        angles = np.loadtxt(rows[j+2:], ndmin=2)
        coords['observation'] = np.arange(len(angles))
        fields['ThetaV'] = ('observation', angles[:,0])
        fields['PhiV'] = ('observation', angles[:,1] % 360)
        dims.append('observation')
        shape.append(len(angles))
    dims.append('N_Stokes_Parameter')
    shape.append(4)
    signal = np.zeros(shape)
    signal[...,0] = key
    signal[...,1] = -1
    signal[...,3] = -0.0
    fields['Stokes_Vector_TOA'] = (dims, signal)
    result = xr.Dataset(fields, coords=coords)
    for level in ['TOA','AirSensor','BOA'] + (['TOO'] if float(values[12]) >= 0 else []):
        if 'Stokes_Vector_' + level not in result:
            result['Stokes_Vector_' + level] = result.Stokes_Vector_TOA.copy()
        for direction in ['Downwelling','Upwelling']:
            result['Irradiance_' + level + '_' + direction] = result.Solar_Irradiance.copy()
    if mode == 'incomplete':
        result = result.drop_vars('Stokes_Vector_BOA')
    outfile = Path(rows[p-1] + '.h5')
else:
    result = xr.Dataset(dict(Radiance_TOA=('Wavelength', [float(key), 0.0]),
                             Q_TOA=('Wavelength', [-1.0, -0.0])),
                        coords=dict(Wavelength=[400.,500.]))
    outfile = Path(rows[-1] + '.h5')
if mode == 'nan':
    field = 'Stokes_Vector_TOA' if mono else 'Radiance_TOA'
    result[field].values.flat[0] = np.nan
if mode != 'missing':
    result.to_netcdf(outfile)
print('synthetic record', key, mode, flush=True)
if mode == 'exit':
    sys.exit(7)
"""


@pytest.fixture
def native_io(tmp_path, monkeypatch):
    scripts = tmp_path / "bin"
    scripts.mkdir()
    for name in ("rtsos_PACE_Simulator_DoubleK", "rtsos_Monochromatic"):
        path = scripts / name
        path.write_text(f"#!{sys.executable}\n" + SYNTHETIC)
        path.chmod(0o755)
    control = tmp_path / "control.json"
    control.write_text("{}")
    monkeypatch.setenv("RTSOS_TEST_CONTROL", str(control))
    monkeypatch.setenv("PATH", str(scripts) + os.pathsep + os.environ["PATH"])
    data = tmp_path / "data"
    data.mkdir()
    solar = data / "f0.txt"
    solar.write_text("# solar\n" * 13 + "300 100\n2500 100\n")
    gas = tmp_path / "gas"
    gas.mkdir()
    return control, data, gas


def calls(control):
    path = control.with_suffix(".calls")
    return [int(v) for v in path.read_text().splitlines()] if path.exists() else []


def sweep():
    return (
        init()
        .defaults.isel(Solar_Zenith_Angle=0, tau_ref_hi=0, chla=0)
        .assign_coords(chla=("chla", np.array([10, 20, 30], dtype=np.float32)))
    )


def case_list(
    tmp_path, native_io, *, methods=("monochromatic",) * 3, data_override=None
):
    _, data, gas = native_io
    entries = []
    for i, (sun, method) in enumerate(zip((10, 20, 30), methods)):
        filename = (
            "response_directions.toml"
            if method == "response_integration"
            else "single_direction.toml"
        )
        text = (EXAMPLES / filename).read_text()
        text = text.replace(
            'root = "../RTSOS_data"',
            f'auxiliary = "{data_override.get(i, data) if data_override else data}"\ngas_absorption = "{gas}"',
        )
        text = text.replace("solar_zenith_deg = 20.0", f"solar_zenith_deg = {sun}.0")
        for name in ("response.toml", "land_spectrum.txt"):
            text = text.replace(f'"{name}"', f'"{EXAMPLES / name}"')
        path = tmp_path / f"case{i}.toml"
        path.write_text(text)
        entries.append(f'[[cases]]\nname = "case{i}"\nfile = "{path.name}"\n')
    listing = tmp_path / "cases.toml"
    listing.write_text("\n".join(entries))
    return listing


def states(report):
    return [r["status"] for r in report["records"]]


def test_legacy_valid_failing_valid_resume_and_explicit_retry(tmp_path, native_io):
    control, _, _ = native_io
    control.write_text('{"20": "exit"}')
    output = tmp_path / "sweep.nc"
    program = init()
    report = program.run(sweep(), output, tmp_path, False)
    assert states(report) == ["successful", "scientifically_failed", "successful"]
    assert calls(control) == [10, 20, 30]
    assert not output.exists()  # no partial population posed as a complete sweep
    directory = records_path(output) / report["request_id"]
    assert (directory / "0/product.nc").is_file()
    assert (directory / "2/product.nc").is_file()
    failed = report["records"][1]
    assert "7" in failed["reason"]
    assert list(Path(failed["attempt"]).glob("*.outfile.h5"))
    assert "20 exit" in (Path(failed["attempt"]) / "native.log").read_text()
    control.write_text("{}")
    again = program.run(sweep(), output, tmp_path, False, resume=True)
    assert states(again) == states(report)
    assert calls(control) == [10, 20, 30]
    repaired = program.run(
        sweep(), output, tmp_path, False, resume=True, retry_failed=True
    )
    assert repaired["complete"] and calls(control) == [10, 20, 30, 20]
    ds = xr.load_dataset(output)
    np.testing.assert_array_equal(ds.chla, [10, 20, 30])
    assert (ds.Radiance_TOA.sel(Wavelength=500) == 0).all()
    assert np.signbit(ds.Q_TOA.sel(Wavelength=500)).all()


@pytest.mark.parametrize(
    "mode,status", [("nan", "scientifically_failed"), ("missing", "retryable")]
)
def test_invalid_or_missing_native_output_continues(tmp_path, native_io, mode, status):
    control, _, _ = native_io
    control.write_text(json.dumps({"20": mode}))
    output = tmp_path / "sweep.nc"
    report = init().run(sweep(), output, tmp_path, False)
    assert states(report) == ["successful", status, "successful"]
    assert calls(control) == [10, 20, 30]
    assert report["records"][1]["reason"]
    control.write_text("{}")
    again = init().run(sweep(), output, tmp_path, False, resume=True)
    assert calls(control) == ([10, 20, 30] if mode == "nan" else [10, 20, 30, 20])
    assert again["complete"] == (mode == "missing")


def test_case_list_products_accounting_cli_and_selected_retry(tmp_path, native_io):
    control, _, _ = native_io
    listing = case_list(
        tmp_path,
        native_io,
        methods=("monochromatic", "monochromatic", "response_integration"),
    )
    control.write_text('{"20": "nan"}')
    output = tmp_path / "results"
    report = simulate_cases(listing, output)
    assert states(report) == ["successful", "scientifically_failed", "successful"]
    assert calls(control) == [10, 20, 30]
    for name in ("case0", "case2"):
        ds = xr.load_dataset(output / f"{name}.nc")
        assert ds.attrs["illumination_policy"] == "thuillier"
        assert ds.observation.dtype == np.dtype("int64")
        assert np.isfinite(ds.Stokes_Vector_TOA).all()
    with Dataset(output / "case2.nc") as result:
        assert list(result["calibration"].groups) == ["0", "1", "2"]
        case = load_case(tmp_path / "case2.toml")
        for i, channel in enumerate(case.load_responses().channels):
            for name in ("wavelength", "response"):
                stored = result[f"calibration/{i}/{name}"][:]
                np.testing.assert_array_equal(stored, getattr(channel, name))
                assert stored.dtype == getattr(channel, name).dtype
    with pytest.raises(SystemExit) as error:
        main(["simulate-cases", str(listing), str(output), "--resume"])
    assert error.value.code == 1
    assert calls(control) == [10, 20, 30]
    control.write_text("{}")
    main(
        [
            "simulate-cases",
            str(listing),
            str(output),
            "--resume",
            "--retry-failed",
            "case1",
        ]
    )
    assert calls(control) == [10, 20, 30, 20]
    assert states(read_json(output / "cases.json")) == ["successful"] * 3


def test_shared_data_block_only_dependents_and_repaired_resume(tmp_path, native_io):
    control, _, _ = native_io
    missing = tmp_path / "missing-data"
    missing.mkdir()
    listing = case_list(tmp_path, native_io, data_override={0: missing, 1: missing})
    output = tmp_path / "results"
    report = simulate_cases(listing, output)
    assert states(report) == ["blocked", "blocked", "successful"]
    assert calls(control) == [30]
    assert all("f0.txt" in r["reason"] for r in report["records"][:2])
    (missing / "f0.txt").write_text("# solar\n" * 13 + "300 100\n2500 100\n")
    repaired = simulate_cases(listing, output, resume=True)
    assert states(repaired) == ["successful"] * 3
    assert calls(control) == [30, 10, 20]


def test_missing_executable_blocks_sweep_and_becomes_schedulable(tmp_path, native_io):
    control, _, _ = native_io
    executable = tmp_path / "bin/rtsos_PACE_Simulator_DoubleK"
    saved = executable.with_suffix(".saved")
    executable.rename(saved)
    # Exclude any real installed binary for this synthetic failure.
    original_path = os.environ["PATH"]
    os.environ["PATH"] = str(executable.parent)
    try:
        output = tmp_path / "sweep.nc"
        report = init().run(sweep(), output, tmp_path, False)
        assert states(report) == ["blocked"] * 3
        assert not calls(control)
        assert len(list(records_path(output).rglob("native.log"))) == 1
        saved.rename(executable)
        repaired = init().run(sweep(), output, tmp_path, False, resume=True)
        assert repaired["complete"] and calls(control) == [10, 20, 30]
    finally:
        os.environ["PATH"] = original_path


def test_cancellation_accounts_remaining_and_resumes(tmp_path, native_io, monkeypatch):
    control, _, _ = native_io
    listing = case_list(tmp_path, native_io)
    output = tmp_path / "results"
    native = subprocess.run

    def cancel(args, **kwargs):
        rows = Path(args[1]).read_text().splitlines()
        if float(rows[1].split("#")[0]) == 20:
            raise KeyboardInterrupt
        return native(args, **kwargs)

    monkeypatch.setattr(subprocess, "run", cancel)
    with pytest.raises(KeyboardInterrupt):
        simulate_cases(listing, output)
    assert states(read_json(output / "cases.json")) == [
        "successful",
        "cancelled",
        "unattempted",
    ]
    assert calls(control) == [10]
    assert (output / "case0.nc").exists()
    monkeypatch.setattr(subprocess, "run", native)
    assert states(simulate_cases(listing, output, resume=True)) == ["successful"] * 3
    assert calls(control) == [10, 20, 30]


def test_optional_accounting_failure_does_not_strand_work(
    tmp_path, native_io, monkeypatch, capsys
):
    control, _, _ = native_io
    listing = case_list(tmp_path, native_io)
    output = tmp_path / "results"

    def unavailable(*args):
        raise OSError("optional reporting unavailable")

    monkeypatch.setattr("rtsos.execution.write_json", unavailable)
    assert states(simulate_cases(listing, output)) == ["successful"] * 3
    assert "optional reporting unavailable" in capsys.readouterr().err
    assert calls(control) == [10, 20, 30]
    assert states(simulate_cases(listing, output, resume=True)) == ["successful"] * 3
    assert calls(control) == [10, 20, 30]


def test_changed_request_data_identity_and_preserved_history(tmp_path, native_io):
    control, data, gas = native_io
    case = load_case(EXAMPLES / "single_direction.toml")
    case.data = dict(auxiliary=str(data), gas_absorption=str(gas))
    output = tmp_path / "single.nc"
    original = simulate(case, output)
    assert calls(control) == [20]
    # Unrelated documentation is outside the runtime/data identity.
    (tmp_path / "README.md").write_text("new documentation")
    simulate(case, output, resume=True)
    assert calls(control) == [20]
    case.geometry["solar_zenith_deg"] = 30
    changed = simulate(case, output, resume=True)
    assert calls(control) == [20, 30]
    assert not original.Stokes_Vector_TOA.equals(changed.Stokes_Vector_TOA)
    assert len(list(records_path(output).glob("*/0/product.nc"))) == 2
    assert list(records_path(output).glob("history/*/single.nc"))
    (data / "f0.txt").write_text("# solar\n" * 13 + "300 101\n2500 101\n")
    simulate(case, output, resume=True)
    assert calls(control) == [20, 30, 30]


def test_partial_metadata_write_retains_raw_and_never_reuses_old_product(
    tmp_path, native_io, monkeypatch
):
    control, data, gas = native_io
    case = load_case(EXAMPLES / "response_directions.toml")
    case.data = dict(auxiliary=str(data), gas_absorption=str(gas))
    output = tmp_path / "bands.nc"
    simulate(case, output)
    writer = __import__(
        "rtsos.bands", fromlist=["write_calibrations"]
    ).write_calibrations

    def interrupted(result, responses):
        result.createGroup("calibration")
        raise RuntimeError("interrupted calibration write")

    monkeypatch.setattr("rtsos.bands.write_calibrations", interrupted)
    with pytest.raises(
        ExecutionFailure, match="interrupted calibration write"
    ) as error:
        simulate(case, output)  # deliberate replacement of the same request
    assert not output.exists()
    assert states(error.value.report) == ["retryable"]
    directory = records_path(output) / error.value.report["request_id"]
    assert not (directory / "0/product.nc").exists()
    assert list(directory.glob("0/.*.partial.nc"))
    assert list(directory.glob("0/*/*.outfile.h5"))
    assert list(directory.glob("0/*/previous-product.nc"))
    monkeypatch.setattr("rtsos.bands.write_calibrations", writer)
    simulate(case, output, resume=True)
    assert calls(control) == [20, 20, 20]
    with Dataset(output) as result:
        assert len(result["calibration"].groups) == 3


def test_public_copy_failure_reuses_durable_success_and_later_cases_run(
    tmp_path, native_io, monkeypatch
):
    control, _, _ = native_io
    listing = case_list(tmp_path, native_io)
    output = tmp_path / "results"
    from rtsos.core import publish_copy

    def fail_first(source, destination):
        if destination == output / "case0.nc":
            raise OSError("public destination unavailable")
        return publish_copy(source, destination)

    monkeypatch.setattr("rtsos.core.publish_copy", fail_first)
    report = simulate_cases(listing, output)
    assert states(report) == ["retryable", "successful", "successful"]
    assert calls(control) == [10, 20, 30]
    monkeypatch.setattr("rtsos.core.publish_copy", publish_copy)
    assert states(simulate_cases(listing, output, resume=True)) == ["successful"] * 3
    assert calls(control) == [10, 20, 30]


def test_selected_combine_rejects_missing_or_changed_population(tmp_path, native_io):
    control, _, _ = native_io
    inputs = tmp_path / "inputs.nc"
    sweep().to_netcdf(inputs)
    output = tmp_path / "sweep.nc"
    args = ["PACE_Simulator_DoubleK", str(inputs), str(output)]
    main([*args, "--split=chla"])
    main([*args, "--combine=chla"])
    child = tmp_path / "sweep/chla/1/sweep.nc"
    child.rename(child.with_suffix(".saved"))
    with pytest.raises(SystemExit) as error:
        main([*args, "--combine=chla"])
    assert error.value.code == 1 and not output.exists()
    child.with_suffix(".saved").rename(child)
    changed = sweep().assign_coords(chla=("chla", np.array([10, 21, 30], np.float32)))
    changed.to_netcdf(inputs)
    with pytest.raises(SystemExit):
        main([*args, "--combine=chla"])
    assert not output.exists() and calls(control) == [10, 20, 30]


def test_shared_identity_hashes_file_bytes_once(tmp_path, monkeypatch):
    root = tmp_path / "data"
    root.mkdir()
    table = root / "table.dat"
    table.write_bytes(b"immutable data")
    counts = []
    original = Path.open

    def tracked(path, *args, **kwargs):
        if path == table and args and args[0] == "rb":
            counts.append(path)
        return original(path, *args, **kwargs)

    monkeypatch.setattr(Path, "open", tracked)
    identity = Identity()
    identity.path(root)
    identity.path(root)
    identity.path(table)
    assert counts == [table]


def test_changed_stored_product_is_recomputed_independently(tmp_path, native_io):
    control, _, _ = native_io
    output = tmp_path / "sweep.nc"
    report = init().run(sweep(), output, tmp_path, False)
    product = records_path(output) / report["request_id"] / "1/product.nc"
    with Dataset(product, "a") as dataset:
        dataset["Radiance_TOA"][0] = 999.0  # finite corruption still fails the receipt
    repaired = init().run(sweep(), output, tmp_path, False, resume=True)
    assert repaired["complete"]
    assert calls(control) == [10, 20, 30, 20]
    assert (
        xr.load_dataset(output).Radiance_TOA.sel(chla=20, Wavelength=400).item() == 20
    )


def test_locked_split_does_not_block_unrelated_selected_work(tmp_path, native_io):
    from rtsos.execution import output_lock

    control, _, _ = native_io
    inputs = tmp_path / "inputs.nc"
    sweep().to_netcdf(inputs)
    output = tmp_path / "sweep.nc"
    child = tmp_path / "sweep/chla/1/sweep.nc"
    with output_lock(records_path(child)):
        with pytest.raises(SystemExit) as error:
            main(["PACE_Simulator_DoubleK", str(inputs), str(output), "--split=chla"])
    assert error.value.code == 1
    assert calls(control) == [10, 30]
    main(
        ["PACE_Simulator_DoubleK", str(inputs), str(output), "--split=chla", "--resume"]
    )
    assert calls(control) == [10, 30, 20]


def test_locked_case_does_not_strand_later_cases(tmp_path, native_io):
    from rtsos.execution import output_lock

    control, _, _ = native_io
    listing = case_list(tmp_path, native_io)
    output = tmp_path / "results"
    with output_lock(records_path(output / "case1.nc")):
        report = simulate_cases(listing, output)
    assert states(report) == ["successful", "blocked", "successful"]
    assert calls(control) == [10, 30]
    assert states(simulate_cases(listing, output, resume=True)) == ["successful"] * 3
    assert calls(control) == [10, 30, 20]


def test_direct_repair_supersedes_old_case_failure_summary(tmp_path, native_io):
    control, _, _ = native_io
    listing = case_list(tmp_path, native_io)
    output = tmp_path / "results"
    control.write_text('{"20": "exit"}')
    assert states(simulate_cases(listing, output)) == [
        "successful",
        "scientifically_failed",
        "successful",
    ]
    control.write_text("{}")
    simulate(
        tmp_path / "case1.toml", output / "case1.nc", resume=True, retry_failed=True
    )
    assert states(simulate_cases(listing, output, resume=True)) == ["successful"] * 3
    assert calls(control) == [10, 20, 30, 20]


def test_partial_combine_finishes_independent_destinations(tmp_path, native_io):
    control, _, _ = native_io
    inputs = tmp_path / "inputs.nc"
    sweep().assign_coords(
        tau_ref_hi=("tau_ref_hi", np.array([0.1, 0.2], np.float32))
    ).to_netcdf(inputs)
    output = tmp_path / "sweep.nc"
    args = ["PACE_Simulator_DoubleK", str(inputs), str(output)]
    main([*args, "--split=chla,tau_ref_hi"])
    missing = tmp_path / "sweep/chla/0/tau_ref_hi/1/sweep.nc"
    missing.rename(missing.with_suffix(".saved"))
    with pytest.raises(SystemExit):
        main([*args, "--split=chla", "--combine=tau_ref_hi"])
    assert not (tmp_path / "sweep/chla/0/sweep.nc").exists()
    assert (tmp_path / "sweep/chla/1/sweep.nc").exists()
    assert (tmp_path / "sweep/chla/2/sweep.nc").exists()
    assert calls(control) == [10, 10, 20, 20, 30, 30]


def test_mie_only_preserves_dispatch_without_claiming_radiometric_reuse(
    tmp_path, monkeypatch
):
    from rtsos.GSFC_AC_LUT import init as ac_init

    program = ac_init()
    inputs = program.defaults.isel(
        {dim: 0 for dim in program.defaults.dims}
    ).assign_coords(MIE_TABLE_CAL=np.int32(1))
    calls = []
    monkeypatch.setattr(subprocess, "run", lambda **kwargs: calls.append(kwargs))
    output = tmp_path / "mie.nc"
    first = program.run(inputs, output, tmp_path, False)
    second = program.run(inputs, output, tmp_path, False, resume=True)
    assert first["complete"] and second["complete"] and len(calls) == 2
    assert not output.exists()
    assert all(r["product"] is None for r in second["records"])


def test_split_cancellation_accounts_all_native_groups_and_resumes(
    tmp_path, native_io, monkeypatch
):
    control, _, _ = native_io
    inputs = tmp_path / "inputs.nc"
    sweep().assign_coords(
        tau_ref_hi=("tau_ref_hi", np.array([0.1, 0.2], np.float32))
    ).to_netcdf(inputs)
    output = tmp_path / "sweep.nc"
    args = ["PACE_Simulator_DoubleK", str(inputs), str(output), "--split=chla"]
    native = subprocess.run

    def interrupt(args, **kwargs):
        if calls(control) == [10]:
            raise KeyboardInterrupt
        return native(args, **kwargs)

    monkeypatch.setattr(subprocess, "run", interrupt)
    with pytest.raises(KeyboardInterrupt):
        main(args)
    manifest = read_json(next(records_path(output).glob("dispatch-*.json")))
    assert states(manifest["records"][0]["execution"]) == ["successful", "cancelled"]
    assert all(
        states(record["execution"]) == ["unattempted"] * 2
        for record in manifest["records"][1:]
    )
    monkeypatch.setattr(subprocess, "run", native)
    main([*args, "--resume"])
    assert calls(control) == [10, 10, 20, 20, 30, 30]


def test_resume_with_outputs_inside_input_data_directory(tmp_path, native_io):
    control, data, gas = native_io
    case = load_case(EXAMPLES / "single_direction.toml")
    case.data = dict(auxiliary=str(data), gas_absorption=str(gas))
    output = data / "single.nc"
    simulate(case, output)
    simulate(case, output, resume=True)
    assert calls(control) == [20]


def test_incomplete_native_schema_retains_reason_and_later_case_runs(
    tmp_path, native_io
):
    control, _, _ = native_io
    listing = case_list(tmp_path, native_io)
    control.write_text('{"20": "incomplete"}')
    output = tmp_path / "results"
    report = simulate_cases(listing, output)
    assert states(report) == ["successful", "scientifically_failed", "successful"]
    assert "Stokes_Vector_BOA" in report["records"][1]["reason"]
    assert calls(control) == [10, 20, 30] and not (output / "case1.nc").exists()


def test_record_log_permission_failure_keeps_later_work_runnable(
    tmp_path, native_io, monkeypatch
):
    control, _, _ = native_io
    original = Path.open

    def inaccessible(path, *args, **kwargs):
        if path.name == "native.log" and path.parent.parent.name == "1":
            raise PermissionError(errno.EACCES, "record log unavailable", str(path))
        return original(path, *args, **kwargs)

    monkeypatch.setattr(Path, "open", inaccessible)
    output = tmp_path / "sweep.nc"
    report = init().run(sweep(), output, tmp_path, False)
    assert states(report) == ["successful", "retryable", "successful"]
    assert calls(control) == [10, 30]
    assert "record log unavailable" in report["records"][1]["reason"]
    monkeypatch.setattr(Path, "open", original)
    assert init().run(sweep(), output, tmp_path, False, resume=True)["complete"]
    assert calls(control) == [10, 30, 20]


def test_output_read_permission_failure_is_retryable(tmp_path, native_io, monkeypatch):
    control, _, _ = native_io
    program = init()
    original = program.read_output

    def inaccessible(path, inputs):
        if inputs.chla.item() == 20:
            raise PermissionError(
                errno.EACCES, "temporary output access error", str(path)
            )
        return original(path, inputs)

    monkeypatch.setattr(program, "read_output", inaccessible)
    output = tmp_path / "sweep.nc"
    report = program.run(sweep(), output, tmp_path, False)
    assert states(report) == ["successful", "retryable", "successful"]
    monkeypatch.setattr(program, "read_output", original)
    assert program.run(sweep(), output, tmp_path, False, resume=True)["complete"]
    assert calls(control) == [10, 20, 30, 20]


@pytest.mark.parametrize("combine", ["resume", "direct", "staged"])
def test_split_products_share_input_directory_without_changing_identity(
    tmp_path, native_io, combine
):
    control, data, _ = native_io
    inputs = tmp_path / "inputs.nc"
    sweep().assign_coords(
        Aux_Dir=str(data),
        tau_ref_hi=("tau_ref_hi", np.array([0.1, 0.2], np.float32)),
    ).to_netcdf(inputs)
    output = data / "sweep.nc"
    args = ["PACE_Simulator_DoubleK", str(inputs), str(output)]
    selection = "chla,tau_ref_hi"
    main([*args, "--split=" + selection])
    if combine == "staged":
        main([*args, "--split=chla", "--combine=tau_ref_hi"])
        main([*args, "--combine=chla"])
    elif combine == "direct":
        main([*args, "--combine=" + selection])
    else:
        main([*args, "--split=" + selection, "--resume"])
    if combine != "resume":
        np.testing.assert_array_equal(xr.load_dataset(output).chla, [10, 20, 30])
    assert calls(control) == [10, 10, 20, 20, 30, 30]
    # Genuine input changes still invalidate reuse, even under this layout.
    (data / "f0.txt").write_text("# solar\n" * 13 + "300 101\n2500 101\n")
    main([*args, "--split=" + selection, "--resume"])
    assert calls(control) == [10, 10, 20, 20, 30, 30] * 2
