"""Small opt-in model connections; no sensor averaging or accuracy gate.

Reference provenance and tolerances are documented in fixtures/model_reference.json
and examples/qualification.md. Auxiliary tables are supplied by the caller.
"""

from copy import deepcopy
import json
import os
from pathlib import Path

import numpy as np
import pytest

from rtsos import load_case, simulate
from rtsos.monochromatic import case_inputs

EXAMPLES = Path(__file__).resolve().parents[3] / "examples"
DATA = Path(__file__).parent / "fixtures"
MODELS = {
    "none": 0,
    "chlorophyll": 1,
    "chlorophyll_sediment": 2,
    "seven_parameter": 3,
    "lambertian_fresnel": -200,
    "mrpv": -201,
    "snow": -202,
    "ross_li": -203,
}


def model_case(model):
    """One explicit physical state with nonzero active model parameters."""
    ocean = MODELS[model] >= 0
    case = load_case(EXAMPLES / f"monochromatic_{'ocean' if ocean else 'land'}.toml")
    case.spectral.update(wavelength_unit="nm", wavelengths=[551.25])
    case.physics["gas_absorption"] = False
    case.geometry = dict(
        mode="explicit_directions",
        solar_zenith_deg=60.0 if model == "snow" else 20.0,
        sensor_height_km=2.0,
        viewing_zenith_deg=[30.0, 30.0, 145.0],
        relative_azimuth_deg=[63.0, -63.0, 120.0],
    )
    if ocean:
        case.scene["water"]["model"] = model
        if model == "none":
            case.scene["water"] = {"model": model}
            case.physics = {
                key: value
                for key, value in case.physics.items()
                if key
                in {
                    "scattering",
                    "spherical_shell",
                    "atmosphere",
                    "gas_absorption",
                    "sun_glint",
                }
            }
        elif model == "chlorophyll_sediment":
            case.scene["water"].update(
                phytoplankton_index_refraction=1.02,
                phytoplankton_spectral_slope=3.5,
                sediment_index_refraction=1.2,
                sediment_spectral_slope=4.0,
                sediment_concentration=0.2,  # g m-3; nonzero mineral contribution
            )
        elif model == "seven_parameter":
            case.scene["water"].pop("absorption_model")
            case.scene["water"].update(
                adg440=0.03,
                bbp660_BackscatterCoeff=0.003,
                Bp660_BackscatterFraction=0.015,
                Sdg=0.018,
                Sbp=1.0,
                S_Bp=0.1,
            )
    else:
        spectrum = {
            "lambertian_fresnel": DATA / "mixture_spectrum.txt",
            "mrpv": DATA / "mrpv_spectrum.txt",
            "snow": DATA / "snow_spectrum.txt",
            "ross_li": EXAMPLES / "ross_li_spectrum.txt",
        }[model]
        case.scene["surface"] = dict(model=model, spectrum=str(spectrum))
    return case


def native_data(case):
    case.data = dict(
        auxiliary=os.environ["RTSOS_AUXILIARY"],
        gas_absorption=os.environ["RTSOS_GAS_ABSORPTION"],
    )
    if case.scene["surface"]["model"] == "snow":
        folder = Path(case.data["auxiliary"]) / "SnowBRDF_Data/brdfPaperMaterial"
        required = ["albedo.txt"] + [f"eof_{letter}.txt" for letter in "ABCDEF"]
        missing = [
            str(folder / name) for name in required if not (folder / name).is_file()
        ]
        if missing:
            pytest.skip("Snow BRDF tables unavailable: " + ", ".join(missing))
    return case


@pytest.mark.parametrize("model,code", MODELS.items())
def test_model_case_translation(model, code):
    inputs = case_inputs(model_case(model))
    assert inputs.OCEAN_CASE_SELECT.item() == code
    assert inputs.wavelength_nm.item() == 551.25
    assert inputs.ncolinput.item() == 8
    if model == "chlorophyll_sediment":
        assert inputs.sediment_concentration.item() == 0.2
        assert inputs.phytoplankton_spectral_slope.item() == 3.5
    if model == "seven_parameter":
        assert inputs.adg440.item() == 0.03
        assert inputs.bbp660_BackscatterCoeff.item() == 0.003
        assert inputs.Bp660_BackscatterFraction.item() == 0.015


@pytest.mark.native
@pytest.mark.parametrize("model,code", MODELS.items())
def test_native_model_matrix(model, code, tmp_path):
    case = native_data(model_case(model))
    result = simulate(case, tmp_path / "model.nc", tmp=tmp_path)
    assert result.OCEAN_CASE_SELECT.item() == code
    np.testing.assert_array_equal(result.Wavelength_Evaluated, [551.25])
    np.testing.assert_array_equal(result.N_Stokes_Parameter, [1, 2, 3, 4])
    levels = ["TOA", "AirSensor", "BOA"] + (["TOO"] if code >= 0 else [])
    assert ("Stokes_Vector_TOO" in result) == (code >= 0)
    for level in levels:
        stokes = result[f"Stokes_Vector_{level}"]
        assert stokes.dims == ("request", "observation", "N_Stokes_Parameter")
        assert stokes.shape == (1, 3, 4)
        assert stokes.dtype == np.dtype("float64")
        assert stokes.units == "W m-2 sr-1 um-1"
        assert np.isfinite(stokes).all()
        assert (stokes.sel(N_Stokes_Parameter=1) >= 0).all()
        # Snow's original off-node quadratic azimuth interpolation is not exactly
        # mirror symmetric in I; compare both directions with the public reference
        # below. The other models have I,Q even and U,V odd.
        even = slice(1, 2) if model == "snow" else slice(0, 2)
        np.testing.assert_allclose(
            stokes[0, 0, even], stokes[0, 1, even], rtol=2e-12, atol=1e-13
        )
        np.testing.assert_allclose(
            stokes[0, 0, 2:], -stokes[0, 1, 2:], rtol=2e-12, atol=1e-13
        )
        for direction in ("Downwelling", "Upwelling"):
            field = result[f"Irradiance_{level}_{direction}"]
            assert field.dims == ("request",)
            assert field.dtype == np.dtype("float64")
            assert field.units == "W m-2 um-1"
            assert np.isfinite(field).all() and (field >= 0).all()
    np.testing.assert_array_equal(result.Tau_Gas_Absorption, 0)
    # No incident diffuse light at TOA in a downward direction away from the beam.
    np.testing.assert_array_equal(result.Stokes_Vector_TOA[0, 2], 0)
    for name in (
        "Tau_Rayleigh_Extinction",
        "Tau_Aerosol_Extinction",
        "Tau_Aerosol_Scattering",
    ):
        assert result[name].dims == ("request", "Atmosphere_Layer")
        assert np.isfinite(result[name]).all() and (result[name] >= 0).all()
    rows = next(tmp_path.glob("*.infile.txt")).read_text().splitlines()
    assert int(rows[12].split("#")[0]) == code
    if code > 0:
        assert float(rows[15].split("#")[0]) == 0.5
    elif code < 0:
        assert case.scene["surface"]["spectrum"] in [
            row.split("#")[0].strip() for row in rows
        ]
    references = json.loads((DATA / "model_reference.json").read_text())
    reference = references["models"][model]
    # Frozen pre-extraction public preparation, one exact angle and wavelength.
    # Shared solver/equations: this is connection parity, not physical accuracy.
    for name, expected in reference.items():
        actual = result[name].isel(request=0)
        if name.startswith("Stokes_Vector_"):
            actual = actual.isel(observation=0)
        elif name.startswith("Tau_"):
            actual = actual.sum("Atmosphere_Layer")
        np.testing.assert_allclose(
            actual, expected, rtol=2e-12, atol=1e-13, err_msg=name
        )
    if model == "snow":
        for name, expected in references["snow_mirror"].items():
            np.testing.assert_allclose(
                result[name].isel(request=0, observation=1),
                expected,
                rtol=2e-12,
                atol=1e-13,
                err_msg=name,
            )


@pytest.mark.native
@pytest.mark.parametrize("model", ["mrpv", "ross_li"])
def test_native_land_intensity_equation(model, tmp_path):
    """Independent direct-surface equations with polarization term disabled.

    RossThick/LiSparse-R: MODIS ATBD and Tilstra et al. (2024), Appendix B.
    MRPV: public native exponential phase law plus Rahman et al. (1993) hotspot.
    These are different angular laws with nonzero anisotropy coefficients.
    The disabled-atmosphere tau floor uses the existing E0*2e-9 tolerance.
    """
    case = native_data(model_case(model))
    case.physics.update(atmosphere=False, gas_absorption=False)
    case.scene["aerosol"]["tau_ref_hi"] = 0.0
    table = np.loadtxt(case.scene["surface"]["spectrum"], skiprows=2)
    table[:, 4] = 0.0  # Bpol / pBRDFe: isolate the scalar surface equation.
    spectrum = tmp_path / "surface.txt"
    np.savetxt(spectrum, table, header="surface\n2", comments="")
    case.scene["surface"]["spectrum"] = str(spectrum)
    case.geometry["viewing_zenith_deg"] = [30.0]
    case.geometry["relative_azimuth_deg"] = [63.0]
    result = simulate(case, tmp_path / "equation.nc", tmp=tmp_path)
    sun, view, azimuth = np.deg2rad([20.0, 30.0, 180.0 - 63.0])
    mu0, mu = np.cos(sun), np.cos(view)
    cos_phase = mu0 * mu + np.sin(sun) * np.sin(view) * np.cos(azimuth)
    separation = np.hypot(
        np.tan(sun) - np.tan(view) * np.cos(azimuth), np.tan(view) * np.sin(azimuth)
    )
    if model == "mrpv":
        reflectance = (
            0.08
            * np.exp(0.2 * cos_phase)
            * (mu0 * mu * (mu0 + mu)) ** (-0.2)
            * (1 + 0.92 / (1 + separation))
        )
    else:
        phase = np.arccos(cos_phase)
        volume = ((np.pi / 2 - phase) * cos_phase + np.sin(phase)) / (
            mu0 + mu
        ) - np.pi / 4
        sec_sum = 1 / mu0 + 1 / mu
        cos_t = (
            2
            * np.hypot(separation, np.tan(sun) * np.tan(view) * np.sin(azimuth))
            / sec_sum
        )
        t = np.arccos(
            min(1.0, cos_t)
        )  # LiSparse overlap's analytical piecewise domain.
        overlap = (t - np.sin(t) * np.cos(t)) * sec_sum / np.pi
        geometric = overlap - sec_sum + (1 + cos_phase) / (2 * mu0 * mu)
        reflectance = 0.0399 + 0.0245 * volume + 0.0072 * geometric
    expected = result.Solar_Irradiance.item() * mu0 * reflectance / np.pi
    np.testing.assert_allclose(
        result.Stokes_Vector_TOA[0, 0, 0],
        expected,
        rtol=0,
        atol=result.Solar_Irradiance.item() * 2e-9,
    )
    np.testing.assert_allclose(
        result.Stokes_Vector_TOA[0, 0, 1:],
        0,
        rtol=0,
        atol=result.Solar_Irradiance.item() * 2e-9,
    )


def test_uv_land_example_sensor_does_not_change_scene():
    case = load_case(EXAMPLES / "uv_land.toml")
    expected = case_inputs(case)
    renamed = deepcopy(case)
    renamed.sensor["name"] = "another sensor"
    assert expected.identical(case_inputs(renamed))
    np.testing.assert_array_equal(expected.wavelength_nm, [350.0, 380.0])
    assert expected.OCEAN_CASE_SELECT.item() == -203
