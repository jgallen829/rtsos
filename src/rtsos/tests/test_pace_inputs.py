"""PACE records checked against main_PACE_Simulator_DoubleK.f90's READ order."""

from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from rtsos.PACE_Simulator_DoubleK import init
from rtsos.__main__ import main


def records(path):
    return [line.split("#", 1)[0].strip() for line in path.read_text().splitlines()]


@pytest.mark.parametrize("surface", [0, 1, 2, 3, -200, -201, -202, -203])
def test_native_record_order(surface, tmp_path):
    program = init()
    state = program.defaults.isel(Solar_Zenith_Angle=0, tau_ref_hi=0, chla=0)
    state = state.assign_coords(
        OCEAN_CASE_SELECT=np.int32(surface),
        land_ref_spectra_filename="land spectrum.txt",
        CHLA_HOMOGENEITY=np.int32(0),
        OCEAN_RAMAN_FLAG=np.int32(1),
        OCEAN_FCHLA_FLAG=np.int32(0),
        OCEAN_FCDOM_FLAG=np.int32(1),
        ap_select=np.int32(2),
        MONOCHROMATIC_FLAG=np.int32(1),
        hyspectral_flag=np.int32(0),
        NonPhotochemicalQuenching_FLAG=np.int32(1),
        OCEAN_PHMX_ONE=np.int32(0),
        ATMOS_ZERO=np.int32(1),
        SUN_GLINT_FLAG=np.int32(1),
        gas_abs_flag=np.int32(0),
    )
    infile, outfile = program.infile(tmp_path, state)
    actual = records(tmp_path / infile)
    # Explicit numeric records, following the native READs, not program.params.
    expected = [5, 20, 532, 0.1, 15, 2, 0, -1, -1, -1, 0.3, 0.3, surface]
    if surface >= 0:
        expected += [
            0.3,
            200,
            0.03,
            1.02,
            3,
            1.2,
            4,
            0,
            1,
            0.05,
            0.01,
            0.015,
            0.3,
            0.01,
        ]
    expected += [40, 40, 60, 40, 36, 19]
    if surface >= 0:
        expected += [0, 1, 0, 1, 2]  # homogeneity, Raman, Fchla, Fcdom, ap_select
    expected += [1, 0]  # monochromatic, hyperspectral
    if surface >= 0:
        expected += [1, 0]  # NPQ, ocean phase matrix
    expected += [1]  # atmosphere zero
    if surface >= 0:
        expected += [1]  # sun glint
    expected += [0, 345.66, 1.4387, 0.2075, 1013, 0, 2.2, 0]
    np.testing.assert_array_equal(
        [float(v) for v in actual[: len(expected)]],
        np.array(expected, dtype=np.float32).astype(float),
    )
    expected_paths = [
        str(Path("Auxiliary_Files").resolve()) + "/",
        str(Path("Gas_Absorption_Coefficients").resolve()),
        "afglus.dat",
        "",
        "",
    ]
    if surface < 0:
        expected_paths += [str(Path("land spectrum.txt").resolve())]
    assert actual[len(expected) : -1] == expected_paths
    assert actual[-1] == str(tmp_path / outfile.with_suffix(""))
    assert len(actual) == (58 if surface >= 0 else 37)


def test_repeated_surface_selection(tmp_path):
    program = init()
    state = program.defaults.isel(Solar_Zenith_Angle=0, tau_ref_hi=0, chla=0)
    for surface in [1, -200, 0, -203, 2, 3, -201, 1]:
        state = state.assign_coords(
            OCEAN_CASE_SELECT=np.int32(surface),
            land_ref_spectra_filename="surface.txt",
        )
        infile, _ = program.infile(tmp_path, state)
        rows = records(tmp_path / infile)
        assert int(rows[12]) == surface
        if surface >= 0:
            assert float(rows[13]) == np.float32(0.3)
            assert rows[27:33] == ["40", "40", "60", "40", "36", "19"]
            assert len(rows) == 58
        else:
            assert rows[13:19] == ["40", "40", "60", "40", "36", "19"]
            assert rows[-2] == str(Path("surface.txt").resolve())
            assert len(rows) == 37
        assert tuple(program.params) == tuple(program.params)


@pytest.mark.parametrize("missing", [True, False])
def test_land_requires_spectrum(tmp_path, missing):
    state = init().defaults.isel(Solar_Zenith_Angle=0, tau_ref_hi=0, chla=0)
    state = state.assign_coords(OCEAN_CASE_SELECT=np.int32(-200))
    if missing:
        state = state.drop_vars("land_ref_spectra_filename")
    state.to_netcdf(tmp_path / "input.nc")
    with pytest.raises(ValueError, match="land_ref_spectra_filename is required"):
        main(
            [
                "PACE_Simulator_DoubleK",
                "--dry-run",
                str(tmp_path / "input.nc"),
                str(tmp_path / "out.nc"),
            ]
        )
    assert not list(tmp_path.glob("*.infile.txt"))


def test_data_root_roundtrip_and_cli(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    root = Path("longer data directory")
    program = init(data_root=root)
    original = init().defaults
    for name in ["Aux_Dir", "Atmos_Dir"]:
        assert program.defaults[name].attrs == original[name].attrs
    program.defaults.to_netcdf("template.nc")
    loaded = xr.load_dataset("template.nc")
    xr.testing.assert_identical(loaded, program.defaults)
    assert loaded.Aux_Dir.item() == str(root.resolve() / "Data") + "/"
    assert loaded.Atmos_Dir.item() == str(
        root.resolve() / "Gas_Absorption_Coefficients"
    )
    # Also exercise overriding a pre-existing template with legacy short strings.
    state = original.isel(Solar_Zenith_Angle=0, tau_ref_hi=0, chla=0)
    state.to_netcdf("input.nc")
    main(
        [
            "PACE_Simulator_DoubleK",
            "--data-root",
            str(root),
            "--dry-run",
            "input.nc",
            "output.nc",
        ]
    )
    rows = records(next(tmp_path.glob("*.infile.txt")))
    assert rows[52:55] == [loaded.Aux_Dir.item(), loaded.Atmos_Dir.item(), "afglus.dat"]
    main(
        ["PACE_Simulator_DoubleK", "--data-root", str(root), "--pre", "cli-template.nc"]
    )
    xr.testing.assert_identical(xr.load_dataset("cli-template.nc"), loaded)
    assert original.Aux_Dir.item() == "Auxiliary_Files"
    assert original.Atmos_Dir.item() == "Gas_Absorption_Coefficients"


def test_direct_path_replacement(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    state = init().defaults.isel(Solar_Zenith_Angle=0, tau_ref_hi=0, chla=0)
    attrs = state.land_ref_spectra_filename.attrs
    state = state.assign_coords(
        OCEAN_CASE_SELECT=np.int32(-200),
        land_ref_spectra_filename=(
            (),
            "spectra/a longer land spectrum filename.txt",
            attrs,
        ),
        Aerosol_Phasematrix_File_Hi=(
            (),
            "aerosols/top layer.h5",
            state.Aerosol_Phasematrix_File_Hi.attrs,
        ),
    )
    state.to_netcdf("land.nc")
    loaded = xr.load_dataset("land.nc")
    assert loaded.land_ref_spectra_filename.attrs == attrs
    program = init()
    prepared = program.prepare_inputs(loaded)
    infile, _ = program.infile(tmp_path, loaded)
    rows = records(tmp_path / infile)
    assert rows[-2] == str(tmp_path / "spectra/a longer land spectrum filename.txt")
    assert rows[-4] == str(tmp_path / "aerosols/top layer.h5")
    assert prepared.land_ref_spectra_filename.attrs == attrs
    xr.testing.assert_identical(program.prepare_inputs(prepared), prepared)


@pytest.mark.parametrize("surface", [1, -200])
def test_native_path_comment_marker(tmp_path, surface):
    program = init(data_root=tmp_path)
    state = program.defaults.isel(Solar_Zenith_Angle=0, tau_ref_hi=0, chla=0)
    state = state.assign_coords(
        OCEAN_CASE_SELECT=np.int32(surface),
        land_ref_spectra_filename="surface#1.txt",
    )
    if surface < 0:
        with pytest.raises(ValueError, match="cannot contain #"):
            program.infile(tmp_path, state)
    else:
        program.infile(
            tmp_path, state
        )  # An unused land path does not block ocean work.


@pytest.mark.parametrize("surface", [-1, -204, 4])
def test_unsupported_surface(tmp_path, surface):
    program = init()
    state = program.defaults.isel(Solar_Zenith_Angle=0, tau_ref_hi=0, chla=0)
    state = state.assign_coords(
        OCEAN_CASE_SELECT=np.int32(surface), land_ref_spectra_filename="surface.txt"
    )
    with pytest.raises(ValueError, match="Unsupported OCEAN_CASE_SELECT"):
        program.infile(tmp_path, state)
    assert not list(tmp_path.glob("*.infile.txt"))


@pytest.mark.parametrize("mixed", [False, True])
@pytest.mark.parametrize("surface", [0, -200])
def test_output_keeps_native_inactive_flags(tmp_path, monkeypatch, mixed, surface):
    """Land omits ocean flags; case zero reads them but disables their processes."""
    import subprocess

    def native(args, check, **kwargs):
        assert check
        rows = records(Path(args[1]))
        selected = int(rows[12])
        land = selected < 0
        assert len(rows) == (37 if land else 58)
        output = xr.Dataset(
            {
                "OCEAN_RAMAN_FLAG": np.int32(selected > 0),
                "OCEAN_FCHLA_FLAG": np.int32(selected > 0),
                "OCEAN_FCDOM_FLAG": np.int32(selected > 0),
                "CHLA_HOMOGENEITY": np.int32(not land),
                "SUN_GLINT_FLAG": np.int32(land),
                "Stokes_Vector_TOA": ("Stokes", [0.1, -0.01, 0.0, 0.0]),
            }
        )
        output.to_netcdf(Path(rows[-1]).with_suffix(".outfile.h5"))

    monkeypatch.setattr(subprocess, "run", native)
    program = init()
    state = program.defaults.isel(Solar_Zenith_Angle=0, tau_ref_hi=0, chla=0)
    state = state.assign_coords(
        OCEAN_CASE_SELECT=(
            ("OCEAN_CASE_SELECT", np.array([1, surface], np.int32))
            if mixed
            else np.int32(surface)
        ),
        land_ref_spectra_filename="surface.txt",
    )
    assert state.OCEAN_RAMAN_FLAG.item() == 1
    program.run(state, tmp_path / "output.nc", tmp_path, False)
    output = xr.load_dataset(tmp_path / "output.nc")
    if mixed:
        assert output.sel(OCEAN_CASE_SELECT=1).OCEAN_RAMAN_FLAG.item() == 1
        output = output.sel(OCEAN_CASE_SELECT=surface)
    assert output.OCEAN_RAMAN_FLAG.item() == 0
    assert output.SUN_GLINT_FLAG.item() == int(surface < 0)
    assert output.CHLA_HOMOGENEITY.item() == int(surface >= 0)
    assert output.Stokes_Vector_TOA.values.tolist() == [0.1, -0.01, 0.0, 0.0]
    assert output.land_ref_spectra_filename.item() == str(Path("surface.txt").resolve())


@pytest.mark.parametrize("gas_abs,atmos_zero", [(0, 0), (1, 1)])
def test_native_spectral_flag_reset(tmp_path, monkeypatch, gas_abs, atmos_zero):
    import subprocess

    def native(args, check, **kwargs):
        assert check
        rows = records(Path(args[1]))
        assert rows[38] == "0"  # Requested monochromatic flag remains unchanged.
        xr.Dataset(
            {
                "MONOCHROMATIC_FLAG": np.int32(1),
                "SUN_GLINT_FLAG": np.int32(atmos_zero),
                "Radiance_TOA": ("Wavelength", [0.0, 1.0]),
            }
        ).to_netcdf(Path(rows[-1]).with_suffix(".outfile.h5"))

    monkeypatch.setattr(subprocess, "run", native)
    program = init()
    state = program.defaults.isel(Solar_Zenith_Angle=0, tau_ref_hi=0, chla=0)
    state = state.assign_coords(
        gas_abs_flag=np.int32(gas_abs), ATMOS_ZERO=np.int32(atmos_zero)
    )
    program.run(state, tmp_path / "output.nc", tmp_path, False)
    output = xr.load_dataset(tmp_path / "output.nc")
    assert output.MONOCHROMATIC_FLAG.item() == 1
    assert output.SUN_GLINT_FLAG.item() == atmos_zero


def test_native_flag_reset_remains_checked(tmp_path, monkeypatch):
    import subprocess

    def wrong_native(args, check, **kwargs):
        assert check
        rows = records(Path(args[1]))
        # Case zero must report Raman disabled, even when the request was one.
        xr.Dataset({"OCEAN_RAMAN_FLAG": np.int32(1)}).to_netcdf(
            Path(rows[-1]).with_suffix(".outfile.h5")
        )

    monkeypatch.setattr(subprocess, "run", wrong_native)
    program = init()
    state = program.defaults.isel(Solar_Zenith_Angle=0, tau_ref_hi=0, chla=0)
    state = state.assign_coords(OCEAN_CASE_SELECT=np.int32(0))
    report = program.run(state, tmp_path / "output.nc", tmp_path, False)
    assert report["records"][0]["status"] == "scientifically_failed"
    assert "OCEAN_RAMAN_FLAG" in report["records"][0]["reason"]
    assert not (tmp_path / "output.nc").exists()
