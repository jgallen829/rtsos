"""Independent spectral expectations, public request routing and native grouping."""

import os
from pathlib import Path
import subprocess

from netCDF4 import Dataset
import numpy as np
import pytest
import xarray as xr

from rtsos import load_case, simulate
from rtsos.__main__ import main
from rtsos.bands import case_inputs, integrate, integration_grid
from rtsos.monochromatic import Monochromatic, case_inputs as mono_inputs
from rtsos.responses import ChannelResponse, Responses

EXAMPLES = Path(__file__).resolve().parents[3] / "examples"


def test_full_response_requests_scale_and_reuse_without_thinning():
    case = load_case(EXAMPLES / "response_land.toml")
    case.numerics.update(nquadainput=256, MAXMORDINPUT=80)
    waves = np.linspace(400, 700, 2001)
    responses = Responses(
        tuple(
            ChannelResponse("many", waves, np.ones_like(waves), "nm", mirror=m)
            for m in (0, 1)
        ),
        "supplied",
    )
    inputs = case_inputs(case, responses)
    np.testing.assert_array_equal(inputs.wavelength_nm, waves)
    assert inputs.sizes["request"] == 2001
    assert inputs.attrs["rt_fourier_order_limit"] == 80
    assert inputs.attrs["phase_expansion_order"] == 256


def test_constant_affine_signed_weights_units_order_and_zeros():
    # Hand trapezoids: signed area=212.5 nm, first moment=141250 nm^2.
    # Quadrature-weighted wavelength=141250/212.5 nm, not a band center.
    a = ChannelResponse(
        "a",
        np.array([400, 500, 700], np.float32),
        np.array([1, -0.25, 2], np.float32),
        "nm",
        mirror=0,
    )
    b = ChannelResponse("b", [0.45, 0.6], [-1.0, -2.0], "um", mirror=1)
    c = ChannelResponse("a", [0.4, 0.5, 0.7], [1.0, -0.25, 2.0], "um", mirror=1)
    responses = Responses((b, a, c), "hand samples")
    waves = np.array([400.0, 450.0, 500.0, 600.0, 700.0])
    nodes = xr.Dataset(
        {
            "Solar_Irradiance": ("request", 3 + 2 * waves / 1000),
            "Irradiance_TOA_Downwelling": ("request", np.full(5, 8.0)),
            "Stokes_Vector_TOA": (
                ("request", "N_Stokes_Parameter"),
                np.array([[2 + w / 1000, -3, 0, 0] for w in waves]),
            ),
            "Tau_Gas_Absorption": ("request", np.arange(5.0)),
        },
        coords={
            "request": np.arange(5),
            "wavelength_nm": ("request", waves),
            "N_Stokes_Parameter": [1, 2, 3, 4],
        },
    )
    result = integrate(nodes, responses)
    mean_um = np.array([0.55, 141250 / 212.5 / 1000, 141250 / 212.5 / 1000])
    np.testing.assert_allclose(result.Solar_Irradiance, 3 + 2 * mean_um, rtol=2e-15)
    np.testing.assert_allclose(result.Stokes_Vector_TOA[:, 0], 2 + mean_um, rtol=2e-15)
    np.testing.assert_allclose(result.Stokes_Vector_TOA[:, 1], -3, rtol=2e-15)
    np.testing.assert_array_equal(result.Stokes_Vector_TOA[:, 2:], 0)
    np.testing.assert_allclose(result.Irradiance_TOA_Downwelling, 8, rtol=2e-15)
    np.testing.assert_allclose(
        result.response_area, [-0.225, 0.2125, 0.2125], rtol=2e-15
    )
    assert list(zip(result.channel_id.values, result.mirror.values)) == [
        ("b", "1"),
        ("a", "0"),
        ("a", "1"),
    ]
    assert result.Tau_Gas_Absorption.dims == ("spectral_node",)
    np.testing.assert_array_equal(a.response, [1, -0.25, 2])
    assert a.response.dtype == np.dtype("float32")


@pytest.mark.parametrize(
    "waves,weights,error",
    [
        ([500], [1], "at least two"),
        ([500, 400], [1, 1], "strictly increasing"),
        ([400, 400], [1, 1], "strictly increasing"),
        ([0, 500], [1, 1], "positive"),
        ([400, 500], [0, 0], "degenerate"),
        ([400, 500], [-1, 1], "degenerate"),
        ([400, 500], [-1, 1 + np.finfo(float).eps], "degenerate"),
        ([400, 500], [np.nan, 1], "finite and unmasked"),
        (np.ma.array([400, 500], mask=[0, 1]), [1, 1], "finite and unmasked"),
        ([400, 500], np.ma.array([1, 1], mask=[1, 0]), "finite and unmasked"),
    ],
)
def test_integration_rejects_incomplete_or_degenerate_calibrations(
    waves, weights, error
):
    channel = ChannelResponse("x", waves, weights, "nm")
    with pytest.raises(ValueError, match=error):
        integration_grid(channel)


def test_python_cli_complete_variants_and_exact_calibration(
    tmp_path, monkeypatch, complete_native
):
    case = load_case(EXAMPLES / "response_land.toml")
    case.data = {"auxiliary": str(tmp_path), "gas_absorption": str(tmp_path)}
    (tmp_path / "f0.txt").touch()
    monkeypatch.setattr("rtsos.monochromatic.check_data", lambda *a: None)
    calls = []

    def native(args, check, **kwargs):
        rows = Path(args[1]).read_text().splitlines()
        policy = rows.index("thuillier")
        waves = np.array([float(w) for w in rows[policy + 2 :]])
        calls.append(waves)
        assert len(waves) == int(rows[policy + 1])
        output = xr.Dataset(
            {
                "Wavelength_Evaluated": ("evaluation", waves),
                "Solar_Irradiance": ("evaluation", waves + 10),
                "Stokes_Vector_TOA": (
                    ("evaluation", "N_Stokes_Parameter", "PhiV", "ThetaV"),
                    np.broadcast_to(
                        np.array([[w, -1, 0, 0] for w in waves], dtype=np.float64)[
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
        complete_native(output, rows[policy - 1] + ".h5")

    monkeypatch.setattr(subprocess, "run", native)
    actual = simulate(case, tmp_path / "python.nc")
    text = (
        (EXAMPLES / "response_land.toml")
        .read_text()
        .replace(
            'root = "../RTSOS_data"',
            f'auxiliary = "{tmp_path}"\ngas_absorption = "{tmp_path}"',
        )
        .replace(
            'response_file = "response.toml"',
            f'response_file = "{EXAMPLES / "response.toml"}"',
        )
        .replace(
            'spectrum = "land_spectrum.txt"',
            f'spectrum = "{EXAMPLES / "land_spectrum.txt"}"',
        )
    )
    path = tmp_path / "case.toml"
    path.write_text(text)
    main(["simulate", str(path), str(tmp_path / "cli.nc")])
    xr.testing.assert_identical(actual, xr.load_dataset(tmp_path / "cli.nc"))
    assert len(calls) == 2
    assert actual.channel_id.values.tolist() == ["red", "blue", "blue"]
    assert actual.mirror.values.tolist() == ["0", "0", "1"]
    assert actual.sizes["channel"] == 3
    assert actual.Solar_Irradiance.units == "W m-2 um-1"
    assert actual.Stokes_Vector_TOA.units == "W m-2 sr-1 um-1"
    with Dataset(tmp_path / "python.nc") as output:
        for i, channel in enumerate(case.load_responses().channels):
            group = output[f"calibration/{i}"]
            for name in ["wavelength", "response"]:
                original = getattr(channel, name)
                stored = group[name][:]
                np.testing.assert_array_equal(stored, original)
                assert stored.dtype == original.dtype
    # Dry-run uses exactly one native array record and no process.
    simulate(case, tmp_path / "dry/out.nc", dry_run=True)
    assert len(calls) == 2
    assert len(list((tmp_path / "dry").glob("*.infile.txt"))) == 1

    def nonfinite_native(args, check, **kwargs):
        native(args, check)
        rows = Path(args[1]).read_text().splitlines()
        with Dataset(rows[rows.index("thuillier") - 1] + ".h5", "a") as output:
            output["Stokes_Vector_TOA"][0, 0] = np.nan

    monkeypatch.setattr(subprocess, "run", nonfinite_native)
    with pytest.raises(ValueError, match="nonfinite native Stokes"):
        simulate(case, tmp_path / "failed.nc")
    assert not (tmp_path / "failed.nc").exists()


@pytest.mark.native
@pytest.mark.parametrize("surface", ["land", "ocean"])
def test_native_grouped_matches_single_and_repeated_order(surface, tmp_path):
    case = load_case(EXAMPLES / f"monochromatic_{surface}.toml")
    case.data = {
        "auxiliary": os.environ["RTSOS_AUXILIARY"],
        "gas_absorption": os.environ["RTSOS_GAS_ABSORPTION"],
    }
    # A/absorbing/B/A exercises atmospheric restoration and changing optics.
    waves = [551.25, 1400.0, 443.125, 551.25]
    case.spectral.update(wavelength_unit="nm", wavelengths=waves)
    grouped = simulate(case, tmp_path / "grouped.nc")
    for i, wave in enumerate(waves[:3]):
        case.spectral["wavelengths"] = [wave]
        single = simulate(case, tmp_path / f"single{i}.nc")
        for name in single.data_vars:
            np.testing.assert_allclose(
                (
                    grouped[name].isel(request=i)
                    if "request" in grouped[name].dims
                    else grouped[name]
                ),
                single[name].squeeze(),
                rtol=2e-13,
                atol=0,
                err_msg=name,
            )
    for name in grouped.data_vars:
        if "request" in grouped[name].dims:
            np.testing.assert_array_equal(
                grouped[name].isel(request=0), grouped[name].isel(request=3)
            )
    # The native array protocol itself also evaluates A again, independently of
    # Python's exact-value reuse, so a cached A cannot hide native contamination.
    case.spectral["wavelengths"] = waves
    inputs = mono_inputs(case)
    program = Monochromatic()
    infile, outfile = program.infile(tmp_path, inputs)
    rows = (tmp_path / infile).read_text().splitlines()
    policy = rows.index("thuillier")
    (tmp_path / infile).write_text(
        "\n".join(rows[: policy + 1] + [str(len(waves))] + [str(w) for w in waves])
        + "\n"
    )
    subprocess.run([program.program, str(tmp_path / infile)], check=True)
    direct = xr.load_dataset(tmp_path / outfile)
    for name in direct.data_vars:
        if "evaluation" in direct[name].dims:
            np.testing.assert_array_equal(
                direct[name].isel(evaluation=0),
                direct[name].isel(evaluation=3),
                err_msg=name,
            )


@pytest.mark.native
def test_native_sediment_reference_reuse(tmp_path):
    case = load_case(EXAMPLES / "monochromatic_ocean.toml")
    case.data = {
        "auxiliary": os.environ["RTSOS_AUXILIARY"],
        "gas_absorption": os.environ["RTSOS_GAS_ABSORPTION"],
    }
    case.scene["water"].update(
        model="chlorophyll_sediment", sediment_concentration=0.03
    )
    case.physics["gas_absorption"] = False
    case.spectral["wavelengths"] = [551.25, 443.125]
    grouped = simulate(case, tmp_path / "sediment.nc")
    for i, wave in enumerate([551.25, 443.125]):
        case.spectral["wavelengths"] = [wave]
        single = simulate(case, tmp_path / f"single{i}.nc")
        for name in single.data_vars:
            actual = (
                grouped[name].isel(request=i)
                if "request" in grouped[name].dims
                else grouped[name]
            )
            np.testing.assert_allclose(
                actual, single[name].squeeze(), rtol=2e-13, atol=0, err_msg=name
            )


@pytest.mark.native
@pytest.mark.parametrize("surface", ["land", "ocean"])
def test_native_full_response_analytical_surface(surface, tmp_path):
    auxiliary = Path(os.environ["RTSOS_AUXILIARY"])
    for path in auxiliary.iterdir():
        if path.name != "f0.txt":
            (tmp_path / path.name).symlink_to(path)
    (tmp_path / "f0.txt").write_text(
        "# solar\n" * 13
        + "".join(f"{w} {100 + .01*w:.12f}\n" for w in range(300, 2502))
    )
    case = load_case(EXAMPLES / f"response_{surface}.toml")
    case.data = {
        "auxiliary": str(tmp_path),
        "gas_absorption": os.environ["RTSOS_GAS_ABSORPTION"],
    }
    case.physics = {"atmosphere": False, "gas_absorption": False}
    case.scene["aerosol"]["tau_ref_hi"] = 0
    if surface == "ocean":
        case.scene["surface"]["wind_speed_m_s"] = 0
        case.scene["water"] = {"model": "none"}
        case.physics["sun_glint"] = False
    result = simulate(case, tmp_path / "band.nc")
    # Independent primitive interval sum, not production integration.
    expected_solar = []
    for c in case.load_responses().channels:
        w = [
            float(x) * (1000 if c.wavelength_unit == "um" else 1) for x in c.wavelength
        ]
        numerator = denominator = 0.0
        for j in range(len(w) - 1):
            half_step = (w[j + 1] - w[j]) / 2
            denominator += half_step * (c.response[j] + c.response[j + 1])
            numerator += half_step * (
                c.response[j] * (1000 + 0.1 * w[j])
                + c.response[j + 1] * (1000 + 0.1 * w[j + 1])
            )
        expected_solar.append(numerator / denominator)
    expected_solar = np.array(expected_solar)
    np.testing.assert_allclose(result.Solar_Irradiance, expected_solar, rtol=2e-14)
    expected_i = (
        (0.2 if surface == "land" else 0)
        * expected_solar
        * np.cos(np.deg2rad(20))
        / np.pi
    )
    np.testing.assert_allclose(
        result.Stokes_Vector_TOA.sel(ThetaV=0, N_Stokes_Parameter=1),
        np.broadcast_to(expected_i[:, None], (3, 2)),
        rtol=0,
        atol=max(expected_solar) * 2e-9,
    )
    np.testing.assert_array_equal(result.Tau_Gas_Absorption, 0)
