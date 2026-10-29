"""Public-dispatch checks that need no native executable or auxiliary data."""

from itertools import product
from pathlib import Path
import subprocess

import numpy as np
import pytest
import xarray as xr

from rtsos import Aerosol_Phmx_Cal, GSFC_AC_LUT, PACE_Simulator_DoubleK
from rtsos.__main__ import main
from rtsos.core import groupby, reduce_by_coords, split_list_arg

PROGRAMS = ("GSFC_AC_LUT", "PACE_Simulator_DoubleK", "Aerosol_Phmx_Cal")


@pytest.mark.parametrize("program", [None, *PROGRAMS, "simulate", "simulate-cases"])
def test_help(program, capsys):
    with pytest.raises(SystemExit) as exc:
        main(([program] if program else []) + ["--help"])
    assert exc.value.code == 0
    assert "--help" in capsys.readouterr().out


@pytest.mark.parametrize("program", PROGRAMS)
def test_template_and_dry_run(program, tmp_path, monkeypatch):
    def unexpected_native(*args, **kwargs):
        pytest.fail("a dry run attempted native execution")

    monkeypatch.setattr(subprocess, "run", unexpected_native)
    inputs = tmp_path / "inputs.nc"
    main([program, "--pre", str(inputs)])
    ds = xr.load_dataset(inputs)
    assert not ds.data_vars
    # Retain a single state for this smoke check, including indexed dimensions.
    ds.isel({dim: slice(0, 1) for dim in ds.dims}).to_netcdf(inputs)
    main([program, "--dry-run", str(inputs), str(tmp_path / "outputs.nc")])
    assert len(list(tmp_path.glob("*.infile.txt"))) == 1


def test_default_pace_sweep(tmp_path):
    inputs = tmp_path / "inputs.nc"
    main(["PACE_Simulator_DoubleK", "--pre", str(inputs)])
    ds = xr.load_dataset(inputs)
    assert ds.sizes == {"Solar_Zenith_Angle": 2, "tau_ref_hi": 2, "chla": 5}
    assert ds.wv_pace_ref.item() == np.float32(532.0)
    assert ds.tau_ref_hi.dtype == np.float32
    assert ds.OCEAN_RAMAN_FLAG.dtype == np.int32
    main(
        [
            "PACE_Simulator_DoubleK",
            "--dry-run",
            str(inputs),
            str(tmp_path / "outputs.nc"),
        ]
    )
    records = [p.read_text().splitlines() for p in tmp_path.glob("*.infile.txt")]
    assert len(records) == 20
    # These are the original native record positions, independent of self.params.
    actual = {
        (
            float(r[1].split("#")[0]),
            float(r[3].split("#")[0]),
            float(r[15].split("#")[0]),
        )
        for r in records
    }
    assert actual == set(
        product(
            ds.Solar_Zenith_Angle.values.tolist(),
            ds.tau_ref_hi.values.tolist(),
            ds.chla.values.tolist(),
        )
    )
    assert all(len(r) == 58 for r in records)


def test_invalid_legacy_dry_run_is_not_silently_successful(tmp_path):
    inputs = (
        PACE_Simulator_DoubleK.init()
        .defaults.isel(Solar_Zenith_Angle=0, tau_ref_hi=0, chla=0)
        .assign_coords(OCEAN_CASE_SELECT=np.int32(99))
    )
    inputs.to_netcdf(tmp_path / "inputs.nc")
    with pytest.raises(ValueError, match="Unsupported OCEAN_CASE_SELECT"):
        main(
            [
                "PACE_Simulator_DoubleK",
                str(tmp_path / "inputs.nc"),
                str(tmp_path / "outputs.nc"),
                "--dry-run",
            ]
        )


def test_case_list_dry_run_reports_invalid_preparation(tmp_path):
    (tmp_path / "bad.toml").write_text('[sensor]\nname = "bad"\n')
    listing = tmp_path / "cases.toml"
    listing.write_text('[[cases]]\nname = "bad"\nfile = "bad.toml"\n')
    with pytest.raises(SystemExit) as error:
        main(["simulate-cases", str(listing), str(tmp_path / "results"), "--dry-run"])
    assert error.value.code == 1


@pytest.mark.parametrize(
    "spec",
    [
        "chla:0",
        "chla:2",
        "chla:1-3",
        "chla",
        "chla:0-4",
        "chla:1-3,tau_ref_hi:1,Solar_Zenith_Angle:0-1",
        "Solar_Zenith_Angle:0,tau_ref_hi:0,chla:0",
    ],
)
def test_pace_selection(spec):
    inputs = PACE_Simulator_DoubleK.init().defaults
    selection = split_list_arg(spec)
    selected, dirs = reduce_by_coords(selection, inputs, Path("outputs"))
    xr.testing.assert_identical(selected, inputs.isel(selection))
    assert dirs.sizes == {dim: selected.sizes[dim] for dim in selection}
    for path, subset in groupby(selected, dirs):
        expected = Path("outputs")
        for dim, indexer in selection.items():
            positions = np.arange(inputs.sizes[dim])[indexer]
            value = subset[dim].item()
            original = next(int(i) for i in positions if inputs[dim][i].item() == value)
            expected /= Path(dim) / str(original)
            assert subset.sizes[dim] == 1
        assert path == expected


def test_indexed_selection():
    inputs = GSFC_AC_LUT.init().defaults
    selected, dirs = reduce_by_coords({"sza-dt": slice(2, 5)}, inputs, Path("outputs"))
    xr.testing.assert_identical(selected, inputs.isel({"sza-dt": slice(2, 5)}))
    assert dirs.sizes == {"sza-dt": 3}
    assert selected["sza-dt"].values.tolist() == [2, 3, 4]
    for path, subset in groupby(selected, dirs):
        assert path == Path("outputs/sza-dt") / str(subset["sza-dt"].item())
        assert subset.sizes["sza-dt"] == 1
        assert subset.Tau_NIR.min().item() == 0.0
        assert subset.Wave_Mean_Square_Slope.min().item() == 0.0


@pytest.mark.parametrize(
    "model,nmode,names",
    [
        (11, 2, "Aerosol_Model IRH WAVEBAND_SEG_FLAG Aux_Dir"),
        (
            -96,
            2,
            "Aerosol_Model Reff_Cloud Veff_Cloud Mr_Cloud Mi_Cloud WAVEBAND_SEG_FLAG Aux_Dir",
        ),
        (-98, 2, "Aerosol_Model Reff_Cloud Veff_Cloud WAVEBAND_SEG_FLAG Aux_Dir"),
        (
            -99,
            2,
            "Aerosol_Model nmode fmfrac cmsfracs rf1 rf2 rf3 Aerosol_Mixing_Flag dust_ivar dust_ireff Relative_Humidity r0f r0c sigf sigc r0_dl r0_ws r0_bc r0_s r0_ss WAVEBAND_SEG_FLAG Aux_Dir",
        ),
        (
            -99,
            3,
            "Aerosol_Model nmode fmfrac dustfrac dsfrac rf2 rf3 Aerosol_Mixing_Flag dust_ivar dust_ireff Relative_Humidity r0f r0c sigf sigc r0_dl r0_ws r0_bc r0_s r0_ss WAVEBAND_SEG_FLAG Aux_Dir",
        ),
    ],
)
def test_repeated_aerosol_records(model, nmode, names, tmp_path):
    program = Aerosol_Phmx_Cal.init()
    inputs = program.defaults.assign_coords(
        Aerosol_Model=np.int32(model),
        nmode=np.int32(nmode),
        WAVEBAND_SEG_FLAG=("WAVEBAND_SEG_FLAG", np.array([0, 1], np.int32)),
    )
    inputs.to_netcdf(tmp_path / "inputs.nc")
    main(
        [
            "Aerosol_Phmx_Cal",
            "--dry-run",
            str(tmp_path / "inputs.nc"),
            str(tmp_path / "outputs.nc"),
        ]
    )
    records = list(tmp_path.glob("*.infile.txt"))
    assert len(records) == 2
    expected_names = names.split()
    segments = set()
    for path in records:
        rows = path.read_text().splitlines()
        values = [row.split("#")[0].strip() for row in rows[:-1]]
        assert len(values) == len(expected_names)
        segment = int(values[expected_names.index("WAVEBAND_SEG_FLAG")])
        segments.add(segment)
        state = inputs.isel(WAVEBAND_SEG_FLAG=segment)
        assert values == [str(state[name].item()) for name in expected_names]
        # The public parameter list is reusable after serialization as well.
        program.params = state
        assert program.params == tuple(expected_names)
        assert tuple(program.params) == tuple(program.params)
    assert segments == {0, 1}


def test_aerosol_model_sweep(tmp_path):
    inputs = Aerosol_Phmx_Cal.init().defaults.assign_coords(
        Aerosol_Model=("Aerosol_Model", np.array([-96, -98, 11], np.int32))
    )
    inputs.to_netcdf(tmp_path / "inputs.nc")
    main(
        [
            "Aerosol_Phmx_Cal",
            "--dry-run",
            str(tmp_path / "inputs.nc"),
            str(tmp_path / "outputs.nc"),
        ]
    )
    rows = [p.read_text().splitlines() for p in tmp_path.glob("*.infile.txt")]
    assert {int(r[0].split("#")[0]): len(r) for r in rows} == {-96: 8, -98: 6, 11: 5}


@pytest.fixture
def synthetic_native(monkeypatch):
    """Replace only the native subprocess, retaining real input/output file I/O."""
    calls = []

    def run(args, check, **kwargs):
        assert check is True
        executable, infile = args
        rows = Path(infile).read_text().splitlines()
        values = [row.split("#")[0].strip() for row in rows[:-1]]
        if executable.startswith("rtsos_PACE"):
            sz, tau, chla = (float(values[i]) for i in (1, 3, 15))
            output = xr.Dataset(
                {
                    "Radiance_TOA": (
                        "Wavelength",
                        np.array([chla + 10 * tau + sz / 100, 0.0]),
                    )
                },
                coords={"Wavelength": [400.0, 500.0]},
            )
        else:
            assert executable.startswith("rtsos_GSFC_AC_LUT")
            model, rh, sz, tau, flag = (float(values[i]) for i in (7, 9, 11, 12, 18))
            signal = model * 10000 + rh * 1000 + sz * 10 + tau * 100 + flag
            output = xr.Dataset(
                coords={"ThetaV": [0.0, 45.0, 100.0], "Wavelength": [400.0, 500.0]}
            )
            for name in (
                "Radiance_TOA",
                "Radiance_TOA_Glint",
                "Q_TOA",
                "Q_TOA_Glint",
                "U_TOA",
                "U_TOA_Glint",
                "Irrad_Down_TOA",
            ):
                output[name] = (("ThetaV", "Wavelength"), np.full((3, 2), signal))
            output["Radiance_BOA"] = (("ThetaV", "Wavelength"), np.full((3, 2), 2.0))
            output["Tau_Aerosol_Extinction"] = ("NTLYERA", [tau, 0.0])
            output["Tau_Rayleigh_Extinction"] = ("NTLYERA", [0.1, 0.2])
            output["Rayleigh_Depolarization_Ratio"] = ("NTLYERA", [0.03, 0.03])
            output["height"] = ("Altitude", [0.0, 1.0])
        output.to_netcdf(Path(rows[-1]).with_suffix(".outfile.h5"))
        calls.append(tuple(values))

    monkeypatch.setattr(subprocess, "run", run)
    return calls


@pytest.mark.parametrize("partial", [False, True])
def test_pace_split_combine(tmp_path, synthetic_native, partial):
    inputs = PACE_Simulator_DoubleK.init().defaults
    inputs.to_netcdf(tmp_path / "inputs.nc")
    args = [str(tmp_path / "inputs.nc"), str(tmp_path / "outputs.nc")]
    selection = "chla:1-3,tau_ref_hi:1,Solar_Zenith_Angle"
    main(["PACE_Simulator_DoubleK", "--split=" + selection, *args])
    assert len(synthetic_native) == 6
    products = list(
        (tmp_path / "outputs").glob(
            "chla/*/tau_ref_hi/*/Solar_Zenith_Angle/*/outputs.nc"
        )
    )
    assert len(products) == 6
    for path in products:
        ds = xr.load_dataset(path)
        assert (
            ds.sizes["chla"]
            == ds.sizes["tau_ref_hi"]
            == ds.sizes["Solar_Zenith_Angle"]
            == 1
        )
        i = int(path.relative_to(tmp_path).parts[2])
        assert ds.chla.item() == inputs.chla[i].item()
    if partial:
        main(
            [
                "PACE_Simulator_DoubleK",
                "--split=chla:1-3",
                "--combine=tau_ref_hi:1,Solar_Zenith_Angle",
                *args,
            ]
        )
        assert len(list((tmp_path / "outputs").glob("chla/*/outputs.nc"))) == 3
        main(["PACE_Simulator_DoubleK", "--combine=chla:1-3", *args])
    else:
        main(["PACE_Simulator_DoubleK", "--combine=" + selection, *args])
    result = xr.load_dataset(tmp_path / "outputs.nc")
    expected = PACE_Simulator_DoubleK.init().prepare_inputs(
        inputs.isel(chla=slice(1, 4), tau_ref_hi=slice(1, 2))
    )
    xr.testing.assert_identical(
        result.drop_vars("Radiance_TOA").drop_vars("Wavelength"), expected
    )
    signal = (
        expected.chla + 10 * expected.tau_ref_hi + expected.Solar_Zenith_Angle / 100
    )
    actual, expected_signal = xr.broadcast(
        result.Radiance_TOA.sel(Wavelength=400, drop=True), signal
    )
    xr.testing.assert_allclose(actual, expected_signal.transpose(*actual.dims))
    assert (result.Radiance_TOA.sel(Wavelength=500) == 0.0).all()
    assert len(synthetic_native) == 6  # recombination performs no native work


@pytest.mark.parametrize("partial", [False, True])
def test_ac_lut_split_combine(tmp_path, synthetic_native, partial):
    inputs = GSFC_AC_LUT.init().defaults.isel(
        Wave_Mean_Square_Slope=slice(0, 1), Tau_NIR=slice(0, 2)
    )
    inputs.to_netcdf(tmp_path / "inputs.nc")
    args = [str(tmp_path / "inputs.nc"), str(tmp_path / "outputs.nc")]
    selection = "Aerosol_Model:1-2,Relative_Humidity:1,sza-dt:1-3"
    main(["GSFC_AC_LUT", "--split=" + selection, *args])
    assert len(synthetic_native) == 12
    products = list(
        (tmp_path / "outputs").glob(
            "Aerosol_Model/*/Relative_Humidity/1/sza-dt/*/outputs.nc"
        )
    )
    assert len(products) == 6
    for path in products:
        ds = xr.load_dataset(path)
        index = int(path.parent.name)
        model_index = int(path.relative_to(tmp_path / "outputs").parts[1])
        assert ds.Aerosol_Model.item() == inputs.Aerosol_Model[model_index].item()
        assert ds["sza-dt"].item() == index
        assert ds.Solar_Zenith_Angle.item() == inputs.Solar_Zenith_Angle[index].item()
        assert ds.sizes["sza-dt"] == 1
        assert ds.Tau_NIR.values.tolist() == inputs.Tau_NIR.values.tolist()
    if partial:
        main(
            [
                "GSFC_AC_LUT",
                "--split=Aerosol_Model:1-2,Relative_Humidity:1",
                "--combine=sza-dt:1-3",
                *args,
            ]
        )
        paths = sorted(
            (tmp_path / "outputs").glob(
                "Aerosol_Model/*/Relative_Humidity/1/outputs.nc"
            )
        )
        assert len(paths) == 2
    else:
        main(["GSFC_AC_LUT", "--combine=" + selection, *args])
        paths = [tmp_path / "outputs.nc"]
    for path in paths:
        ds = xr.load_dataset(path)
        model_indices = (
            [int(path.relative_to(tmp_path / "outputs").parts[1])]
            if partial
            else [1, 2]
        )
        assert (
            ds.Aerosol_Model.values.tolist()
            == inputs.Aerosol_Model[model_indices].values.tolist()
        )
        assert ds.Solar_Zenith_Angle.values.tolist() == [2.0, 4.0]
        assert ds.ThetaV.values.tolist() == [0.0, 45.0]
        assert ds.Relative_Humidity.values.tolist() == [
            inputs.Relative_Humidity[1].item()
        ]
        assert ds.Tau_NIR.values.tolist() == inputs.Tau_NIR.values.tolist()
        expected = (
            ds.Aerosol_Model * 10000
            + ds.Relative_Humidity * 1000
            + ds.Solar_Zenith_Angle * 10
            + ds.Tau_NIR * 100
        )
        expected, actual = xr.broadcast(expected, ds.Radiance_TOA)
        xr.testing.assert_allclose(actual, expected.transpose(*actual.dims))
        assert (ds.aerosol_optical_thickness.sel(Tau_NIR=0) == 0).all()
    if partial:
        expected = xr.combine_by_coords(
            [xr.load_dataset(path) for path in paths],
            data_vars="different",
            compat="equals",
        )
        main(
            [
                "GSFC_AC_LUT",
                "--combine=Aerosol_Model:1-2,Relative_Humidity:1",
                *args,
            ]
        )
        xr.testing.assert_identical(xr.load_dataset(tmp_path / "outputs.nc"), expected)
    assert len(synthetic_native) == 12
