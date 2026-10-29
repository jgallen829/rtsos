"""Requests through simulate/Program, plus independently solvable native surfaces."""

import os
from pathlib import Path
import subprocess

import numpy as np
import pytest
import xarray as xr

from rtsos import load_case, simulate
from rtsos.__main__ import main
from rtsos.monochromatic import Monochromatic, case_inputs, check_data, solar_source

EXAMPLES = Path(__file__).resolve().parents[3] / "examples"


def ocean():
    return load_case(EXAMPLES / "monochromatic_ocean.toml")


def test_request_records_units_precision_and_order(tmp_path):
    case = ocean()
    case.spectral.update(
        wavelength_unit="um", wavelengths=[0.55125123456789, 0.443125, 0.55125123456789]
    )
    case.scene["water"]["chla"] = 0.123456789123456
    inputs = case_inputs(case)
    assert inputs.chla.item() == case.scene["water"]["chla"]
    assert inputs.wavelength_nm.dtype == np.dtype("float64")
    program = Monochromatic()
    for i in range(3):
        infile, _ = program.infile(tmp_path, inputs.isel(request=i))
        rows = (tmp_path / infile).read_text().splitlines()
        assert float(rows[-1]) == case.spectral["wavelengths"][i] * 1000
        assert rows[-3:-1] == ["thuillier", "1"]
        assert float(rows[15].split("#")[0]) == case.scene["water"]["chla"]


def test_python_cli_real_program_route(tmp_path, monkeypatch, complete_native):
    case = ocean()
    case.spectral["wavelengths"] = [551.25, 443.125, 551.25]
    case.scene["aerosol"]["tau_ref_hi"] = 0.0
    case.data = {"auxiliary": str(tmp_path), "gas_absorption": str(tmp_path)}
    (tmp_path / "f0.txt").touch()
    monkeypatch.setattr("rtsos.monochromatic.check_data", lambda *a: None)
    records = []

    def native(args, check, **kwargs):
        assert args[0] == "rtsos_Monochromatic" and check
        rows = Path(args[1]).read_text().splitlines()
        policy = rows.index("thuillier")
        waves = np.array([float(w) for w in rows[policy + 2 :]])
        assert int(rows[policy + 1]) == 2  # Exact A/B/A reuse in one process.
        records.append(rows[: policy - 1] + rows[policy:])
        output = xr.Dataset(
            {
                "Wavelength_Evaluated": ("evaluation", waves),
                "Stokes_Vector_TOA": (
                    ("evaluation", "N_Stokes_Parameter", "PhiV", "ThetaV"),
                    np.broadcast_to(
                        np.array([[w, -0.02, 0.0, -0.0] for w in waves], np.float64)[
                            :, :, None, None
                        ],
                        (len(waves), 4, 2, 2),
                    ),
                ),
            },
            coords={
                "evaluation": np.arange(len(waves)),
                "N_Stokes_Parameter": [1, 2, 3, 4],
            },
        )
        complete_native(output, rows[policy - 1] + ".h5", ocean=True)

    monkeypatch.setattr(subprocess, "run", native)
    result = simulate(case, tmp_path / "python.nc", tmp=tmp_path)
    # CLI uses the same TOML request and adapter, including repeat identities.
    text = (EXAMPLES / "monochromatic_ocean.toml").read_text()
    text = text.replace(
        'root = "../RTSOS_data"',
        f'auxiliary = "{tmp_path}"\ngas_absorption = "{tmp_path}"',
    )
    text = text.replace("[551.25, 443.125]", "[551.25, 443.125, 551.25]").replace(
        "tau_ref_hi = 0.125", "tau_ref_hi = 0.0"
    )
    path = tmp_path / "case.toml"
    path.write_text(text)
    main(["simulate", str(path), str(tmp_path / "cli.nc"), "--tmp", str(tmp_path)])
    xr.testing.assert_identical(result, xr.load_dataset(tmp_path / "cli.nc"))
    np.testing.assert_array_equal(
        result.Wavelength_Evaluated, [551.25, 443.125, 551.25]
    )
    assert result.sizes["request"] == 3
    assert result.Stokes_Vector_TOA.dtype == np.dtype("float64")
    assert np.signbit(result.Stokes_Vector_TOA.sel(N_Stokes_Parameter=4)).all()
    assert result.tau_ref_hi.item() == 0
    assert len(records) == 2 and records[0] == records[1]


@pytest.mark.parametrize(
    "process,operation",
    [
        ("raman", "Raman"),
        ("chlorophyll_fluorescence", "chlorophyll fluorescence"),
        ("cdom_fluorescence", "CDOM fluorescence"),
        ("nonphotochemical_quenching", "NPQ"),
    ],
)
def test_unsupported_process_before_execution(tmp_path, process, operation):
    case = ocean()
    case.physics[process] = True
    with pytest.raises(NotImplementedError, match=operation):
        simulate(case, tmp_path / "out.nc", dry_run=True)
    assert not list(tmp_path.iterdir())


def test_missing_external_matrix_operation():
    case = ocean()
    case.scene["aerosol"]["Aerosol_Model"] = -99
    with pytest.raises(
        NotImplementedError, match="requested-wavelength matrix loading"
    ):
        case_inputs(case)


@pytest.mark.parametrize("waves", [[], [0], [-1], [float("nan")], [[440]], ["440"]])
def test_wavelength_boundary(waves):
    case = ocean()
    case.spectral["wavelengths"] = waves
    with pytest.raises(ValueError, match="wavelengths"):
        case_inputs(case)


def test_settings_have_no_pace_clamp_or_workload_cap():
    case = ocean()
    case.numerics.update(nquadainput=256, MAXMORDINPUT=80)
    case.spectral["wavelengths"] = list(np.linspace(400, 700, 1001))
    inputs = case_inputs(case)
    assert inputs.sizes["request"] == 1001
    assert inputs.attrs["rt_fourier_order_limit"] == 80
    assert inputs.attrs["phase_expansion_order"] == 256
    assert "pace_rt_fourier_order_limit" not in inputs.attrs


def test_explicit_solar_choice_independent_of_precedence(tmp_path):
    case = ocean()
    case.data = {"auxiliary": str(tmp_path), "gas_absorption": str(tmp_path)}
    (tmp_path / "f0.txt").touch()
    (tmp_path / "tsis1_hsrs_binned_fs.csv").touch()
    assert solar_source(case).name == "f0.txt"
    case.illumination["source"] = "tsis"
    assert solar_source(case).name == "tsis1_hsrs_binned_fs.csv"


def test_active_data_coverage_without_response_files(tmp_path):
    case = load_case(EXAMPLES / "monochromatic_land.toml")
    case.data = {"auxiliary": str(tmp_path), "gas_absorption": str(tmp_path)}
    path = tmp_path / "f0.txt"
    path.write_text("# header\n" * 13 + "400 100\n500 100\n600 100\n")
    inputs = case_inputs(case)
    check_data(case, inputs, path)
    # No instrument calibration files have been created or consulted.
    case.spectral.update(wavelength_unit="nm", wavelengths=[650])
    with pytest.raises(ValueError, match="solar file.*coverage"):
        check_data(case, case_inputs(case), path)
    case.spectral["wavelengths"] = [250]
    with pytest.raises(ValueError, match="aerosol refractive-index.*coverage"):
        check_data(case, case_inputs(case), path)


def test_no_atmosphere_does_not_silently_disable_requested_gases():
    case = ocean()
    case.physics.update(atmosphere=False, gas_absorption=True)
    with pytest.raises(ValueError, match="requires an atmosphere"):
        case_inputs(case)


@pytest.mark.native
@pytest.mark.parametrize("surface", ["land", "ocean"])
def test_native_analytical_surface(surface, tmp_path):
    """Affine solar spectrum; Lambertian land or black surface-only ocean.

    Expected values use radiance=albedo*E0*cos(SZA)/pi and a black-ocean
    zero. The inherited disabled-atmosphere floor is 102*1e-12 in tau;
    atol accounts for that residual scattering, not fitting a solver result.
    """
    auxiliary = Path(os.environ["RTSOS_AUXILIARY"])
    for path in auxiliary.iterdir():
        if path.name != "f0.txt":
            (tmp_path / path.name).symlink_to(path)
    (tmp_path / "f0.txt").write_text(
        "# solar\n" * 13
        + "".join(f"{w} {100 + .01*w:.12f}\n" for w in range(300, 2502))
    )
    case = load_case(EXAMPLES / f"monochromatic_{surface}.toml")
    case.data = {
        "auxiliary": str(tmp_path),
        "gas_absorption": os.environ["RTSOS_GAS_ABSORPTION"],
    }
    case.spectral = {
        "method": "monochromatic",
        "wavelength_unit": "nm",
        "wavelengths": [551.25],
    }
    case.physics = {"atmosphere": False, "gas_absorption": False}
    case.scene["aerosol"]["tau_ref_hi"] = 0
    if surface == "ocean":
        case.scene["surface"]["wind_speed_m_s"] = 0
        case.scene["water"] = {"model": "none"}
        case.physics["sun_glint"] = False
    result = simulate(case, tmp_path / "surface.nc", tmp=tmp_path)
    e0 = 1000 + 0.1 * 551.25  # f0 mW cm-2 um-1 -> W m-2 um-1
    np.testing.assert_allclose(result.Solar_Irradiance, e0, rtol=2e-14)
    expected_i = (0.2 if surface == "land" else 0) * e0 * np.cos(np.deg2rad(20)) / np.pi
    actual = result.Stokes_Vector_TOA.sel(ThetaV=0, N_Stokes_Parameter=1)
    np.testing.assert_allclose(actual, expected_i, rtol=0, atol=e0 * 2e-9)
    np.testing.assert_array_equal(result.Tau_Gas_Absorption, 0)
    np.testing.assert_array_equal(result.Tau_Aerosol_Extinction, 0)


@pytest.mark.native
@pytest.mark.parametrize("mode", ["pace_grid", "explicit_directions"])
def test_native_absorbing_atmosphere_restores_attenuation(tmp_path, mode):
    """A black surface cannot receive more than E0*mu0*exp(-vertical gas tau).

    Every downward photon traverses at least the vertical absorber column. This
    bound catches loss of the attenuation restoration after the native driver
    reduces absorption depth above ten; it needs no reference from that driver.
    """
    case = load_case(EXAMPLES / "monochromatic_land.toml")
    case.data = {
        "auxiliary": os.environ["RTSOS_AUXILIARY"],
        "gas_absorption": os.environ["RTSOS_GAS_ABSORPTION"],
    }
    spectrum = tmp_path / "black.txt"
    spectrum.write_text("black\n2\n300 0 1 0 0 1.5 0\n2500 0 1 0 0 1.5 0\n")
    case.scene["surface"]["spectrum"] = str(spectrum)
    case.spectral.update(wavelength_unit="nm", wavelengths=[1400.0])
    if mode == "explicit_directions":
        case.geometry = dict(
            mode=mode,
            solar_zenith_deg=20.0,
            sensor_height_km=2.0,
            viewing_zenith_deg=[37.25, 0, 145, 37.25],
            relative_azimuth_deg=[315, 0, 120, -45],
        )
    result = simulate(case, tmp_path / "absorbing.nc", tmp=tmp_path)
    gas_tau = result.Tau_Gas_Absorption.sum("Atmosphere_Layer")
    assert (gas_tau > 10).all()  # Exercise the native absorption-reduction branch.
    upper_bound = (
        result.Solar_Irradiance
        * np.cos(np.deg2rad(case.geometry["solar_zenith_deg"]))
        * np.exp(-gas_tau)
    )
    assert (result.Irradiance_BOA_Downwelling >= 0).all()
    assert (result.Irradiance_BOA_Downwelling <= upper_bound).all()
