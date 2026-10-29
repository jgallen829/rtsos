"""Ordered directions across the public interfaces and native angular boundary."""

from copy import deepcopy
import os
from pathlib import Path
import subprocess

import numpy as np
import pytest
import xarray as xr

from rtsos import load_case, simulate
from rtsos.__main__ import main
from rtsos.monochromatic import Monochromatic, case_inputs
from rtsos.PACE_Simulator_DoubleK import case_inputs as pace_inputs

EXAMPLES = Path(__file__).resolve().parents[3] / "examples"


def directions(theta=(37.25, 0, 145, 37.25, 30), phi=(315, 0, 120, -45, 63)):
    return dict(
        mode="explicit_directions",
        solar_zenith_deg=20.0,
        sensor_height_km=2.0,
        viewing_zenith_deg=theta if np.ma.isMaskedArray(theta) else list(theta),
        relative_azimuth_deg=list(phi),
    )


def test_paired_records_precision_zeros_and_no_count_limit(tmp_path):
    case = load_case(EXAMPLES / "monochromatic_directions.toml")
    case.geometry = directions([0, 37.25123456789123, 180], [0, -45, 405])
    inputs = case_inputs(case)
    assert inputs.viewing_zenith_requested.dtype == np.dtype("float64")
    infile, _ = Monochromatic().infile(tmp_path, inputs)
    rows = (tmp_path / infile).read_text().splitlines()
    i = rows.index("explicit_directions")
    assert rows[i + 1] == "3"
    np.testing.assert_array_equal(
        np.loadtxt(rows[i + 2 :]), [[0, 0], [37.25123456789123, -45], [180, 405]]
    )
    # This is an operating size, not a maximum. No Cartesian product in public inputs.
    case.geometry = directions(np.linspace(0, 89, 2001), np.arange(2001))
    assert case_inputs(case).sizes["observation"] == 2001
    case.geometry = directions([89.999999999, 90.000000001], [0, 360])
    assert case_inputs(case).sizes["observation"] == 2


@pytest.mark.parametrize(
    "theta,phi,error,exception",
    [
        ([], [], "nonempty", ValueError),
        ([0, 30], [0], "equal lengths", ValueError),
        ([[0]], [0], "list", ValueError),
        ([float("nan")], [0], "finite", ValueError),
        ([0], [float("inf")], "finite", ValueError),
        ([True], [0], "numeric", ValueError),
        (np.ma.array([0], mask=[1]), [0], "finite", ValueError),
        ([-1], [0], "\\[0, 180\\]", ValueError),
        ([181], [0], "\\[0, 180\\]", ValueError),
        ([90], [0], "mu=0 transport/direct-beam limit", NotImplementedError),
    ],
)
def test_direction_boundary(theta, phi, error, exception):
    case = load_case(EXAMPLES / "monochromatic_directions.toml")
    case.geometry = directions(theta, phi)
    with pytest.raises(exception, match=error):
        case_inputs(case)


def test_legacy_capability_error_and_sensor_independence():
    case = load_case(EXAMPLES / "pace_ocean.toml")
    case.geometry = directions()
    with pytest.raises(NotImplementedError, match="legacy PACE double-k"):
        pace_inputs(case)
    mono = load_case(EXAMPLES / "monochromatic_directions.toml")
    response = load_case(EXAMPLES / "response_directions.toml")
    assert mono.geometry == response.geometry
    assert mono.scene == response.scene
    assert mono.physics == response.physics
    assert mono.illumination == response.illumination
    inputs = case_inputs(mono)
    mono.sensor["name"] = "another sensor"
    xr.testing.assert_identical(inputs, case_inputs(mono))
    mono.geometry["zenith_count"] = 2
    with pytest.raises(ValueError, match="Unsupported geometry"):
        case_inputs(mono)


@pytest.mark.parametrize("method", ["monochromatic", "response"])
@pytest.mark.parametrize("single", [False, True])
def test_python_cli_observation_identity_and_dimensions(
    tmp_path, monkeypatch, method, single, complete_native
):
    filename = f"{method}_directions.toml"
    text = (
        (EXAMPLES / filename)
        .read_text()
        .replace(
            'root = "../RTSOS_data"',
            f'auxiliary = "{tmp_path}"\ngas_absorption = "{tmp_path}"',
        )
    )
    for name in ("response.toml", "land_spectrum.txt"):
        text = text.replace(f'"{name}"', f'"{EXAMPLES / name}"')
    if single:
        text = text.replace("[37.25, 0.0, 145.0, 37.25, 30.0]", "[37.25]").replace(
            "[315.0, 0.0, 120.0, -45.0, 63.0]", "[315.0]"
        )
    path = tmp_path / "case.toml"
    path.write_text(text)
    (tmp_path / "f0.txt").touch()
    monkeypatch.setattr("rtsos.monochromatic.check_data", lambda *a: None)
    calls = []

    def native(args, check, **kwargs):
        assert args[0] == "rtsos_Monochromatic" and check
        rows = Path(args[1]).read_text().splitlines()
        i = rows.index("thuillier")
        j = rows.index("explicit_directions")
        waves = np.array([float(w) for w in rows[i + 2 : j]])
        assert len(waves) == int(rows[i + 1])
        angles = np.loadtxt(rows[j + 2 :], ndmin=2)
        assert len(angles) == int(rows[j + 1])
        theta, phi = angles.T
        # A synthetic native boundary checks transport/order, not RT physics.
        values = np.zeros((len(waves), len(theta), 4), np.float64)
        values[:, :, 0] = waves[:, None] + theta + phi
        values[:, :, 1] = -1
        values[:, :, 3] = -0.0
        output = xr.Dataset(
            {
                "Wavelength_Evaluated": ("evaluation", waves),
                "Solar_Irradiance": ("evaluation", waves + 10),
                "Irradiance_TOA_Downwelling": ("evaluation", waves + 5),
                "Tau_Gas_Absorption": (
                    ("evaluation", "Atmosphere_Layer"),
                    np.zeros((len(waves), 2)),
                ),
                "Stokes_Vector_TOA": (
                    ("evaluation", "observation", "N_Stokes_Parameter"),
                    values,
                ),
                "ThetaV": ("observation", theta),
                "PhiV": ("observation", phi % 360),
            },
            coords={
                "evaluation": np.arange(len(waves)),
                "observation": np.arange(len(theta)),
                "N_Stokes_Parameter": [1, 2, 3, 4],
                "Atmosphere_Layer": [1, 2],
            },
        )
        complete_native(output, rows[i - 1] + ".h5")
        calls.append(angles)

    monkeypatch.setattr(subprocess, "run", native)
    result = simulate(path, tmp_path / "python.nc", tmp=tmp_path)
    main(["simulate", str(path), str(tmp_path / "cli.nc"), "--tmp", str(tmp_path)])
    xr.testing.assert_identical(result, xr.load_dataset(tmp_path / "cli.nc"))
    assert len(calls) == 2
    expected = load_case(path).geometry
    np.testing.assert_array_equal(result.observation, np.arange(1 if single else 5))
    np.testing.assert_array_equal(
        result.viewing_zenith_requested, expected["viewing_zenith_deg"]
    )
    np.testing.assert_array_equal(
        result.relative_azimuth_requested, expected["relative_azimuth_deg"]
    )
    assert result.ThetaV.dims == result.PhiV.dims == ("observation",)
    spectral_axis = "request" if method == "monochromatic" else "channel"
    assert result.Stokes_Vector_TOA.dims == (
        spectral_axis,
        "observation",
        "N_Stokes_Parameter",
    )
    assert result.Stokes_Vector_TOA.dtype == np.dtype("float64")
    assert result.Stokes_Vector_TOA.units == "W m-2 sr-1 um-1"
    assert result.Irradiance_TOA_Downwelling.dims == (spectral_axis,)
    assert "observation" not in result.Tau_Gas_Absorption.dims
    np.testing.assert_array_equal(result.Stokes_Vector_TOA.sel(N_Stokes_Parameter=3), 0)
    if method == "response":
        assert result.channel_id.values.tolist() == ["red", "blue", "blue"]
        assert result.mirror.values.tolist() == ["0", "0", "1"]
        assert result.Node_Stokes_Vector_TOA.dims == (
            "spectral_node",
            "observation",
            "N_Stokes_Parameter",
        )
    else:
        assert np.signbit(result.Stokes_Vector_TOA.sel(N_Stokes_Parameter=4)).all()


def native_case(surface="land"):
    case = load_case(EXAMPLES / f"monochromatic_{surface}.toml")
    case.data = dict(
        auxiliary=os.environ["RTSOS_AUXILIARY"],
        gas_absorption=os.environ["RTSOS_GAS_ABSORPTION"],
    )
    return case


@pytest.mark.native
@pytest.mark.parametrize("surface", ["land", "ocean"])
def test_native_pairs_match_unequal_legacy_axes(surface, tmp_path):
    case = native_case(surface)
    # Two scattering orders use the native unmapped output branch. This isolates
    # grid/coordinate pairing from the mapped scratch-fit correction (tested below).
    case.numerics["ncolinput"] = 2
    case.geometry.update(zenith_count=6, azimuth_count=4)
    case.spectral.update(wavelength_unit="nm", wavelengths=[551.25, 443.125, 551.25])
    grid = simulate(case, tmp_path / "grid.nc")
    assert grid.sizes["ThetaV"] == 6 and grid.sizes["PhiV"] == 4
    # The expected directions come from the existing grid constructor, not the new adapter.
    ti = [4, 0, 2, 1, 4, 1]
    pi = [2, 0, 3, 1, 2, 3]
    case.geometry = directions([150, 0, 60, 30, 150, 30], [120, 0, 180, 60, 120, 180])
    paired = simulate(case, tmp_path / "paired.nc")
    np.testing.assert_allclose(
        paired.ThetaV, grid.ThetaV.values[ti], rtol=0, atol=3e-14
    )
    np.testing.assert_allclose(paired.PhiV, grid.PhiV.values[pi], rtol=0, atol=3e-14)
    for name in grid.data_vars:
        if name.startswith("Stokes_Vector_"):
            expected = (
                grid[name]
                .isel(
                    ThetaV=xr.DataArray(ti, dims="observation"),
                    PhiV=xr.DataArray(pi, dims="observation"),
                )
                .transpose("request", "observation", "N_Stokes_Parameter")
            )
            np.testing.assert_allclose(
                paired[name], expected, rtol=2e-12, atol=1e-13, err_msg=name
            )
            np.testing.assert_array_equal(
                paired[name].isel(observation=0), paired[name].isel(observation=4)
            )
        elif name.startswith(("Irradiance_", "Tau_")) or name == "Solar_Irradiance":
            np.testing.assert_array_equal(paired[name], grid[name], err_msg=name)


@pytest.mark.native
@pytest.mark.parametrize("surface", ["land", "ocean"])
def test_native_off_grid_close_directions_single_and_wavelength_reuse(
    surface, tmp_path
):
    case = native_case(surface)
    case.physics["gas_absorption"] = False
    case.geometry = directions(
        [37.25, 37.250001, 0, 180, 37.25, 37.25], [63, 63, 315, 120, 63, -63]
    )
    waves = [551.25, 443.125, 551.25]
    case.spectral.update(wavelength_unit="nm", wavelengths=waves)
    grouped = simulate(case, tmp_path / "grouped.nc")
    for i in range(4):
        single = deepcopy(case)
        single.geometry = directions(
            [case.geometry["viewing_zenith_deg"][i]],
            [case.geometry["relative_azimuth_deg"][i]],
        )
        single.spectral["wavelengths"] = [waves[i % 2]]
        result = simulate(single, tmp_path / f"single{i}.nc")
        assert result.sizes["observation"] == result.sizes["request"] == 1
        for name in result.data_vars:
            if name.startswith("Stokes_Vector_"):
                np.testing.assert_allclose(
                    grouped[name].isel(request=i % 2, observation=i),
                    result[name].isel(request=0, observation=0),
                    rtol=2e-12,
                    atol=1e-13,
                    err_msg=name,
                )
    assert (
        grouped.ThetaV[1] > grouped.ThetaV[0]
    )  # Below the legacy 1e-6 mu merge tolerance.
    for name in grouped.data_vars:
        if name.startswith("Stokes_Vector_"):
            s = grouped[name].values
            np.testing.assert_array_equal(s[:, 0], s[:, 4])
            np.testing.assert_allclose(s[:, 0, :2], s[:, 5, :2], rtol=2e-12, atol=1e-13)
            np.testing.assert_allclose(
                s[:, 0, 2:], -s[:, 5, 2:], rtol=2e-12, atol=1e-13
            )
        if "request" in grouped[name].dims:
            np.testing.assert_array_equal(
                grouped[name].isel(request=0), grouped[name].isel(request=2)
            )
    assert (
        abs(
            grouped.Stokes_Vector_TOA.sel(N_Stokes_Parameter=4).isel(
                request=0, observation=0
            )
        )
        > 1e-9
    )
    # Force native A/B/A, bypassing Python wavelength reuse while preserving the pair records.
    inputs = case_inputs(case)
    program = Monochromatic()
    infile, outfile = program.infile(tmp_path, inputs)
    path = tmp_path / infile
    rows = path.read_text().splitlines()
    i = rows.index("thuillier")
    j = rows.index("explicit_directions")
    path.write_text(
        "\n".join(rows[: i + 1] + ["3"] + [str(w) for w in waves] + rows[j:]) + "\n"
    )
    subprocess.run([program.program, str(path)], check=True)
    native = xr.load_dataset(tmp_path / outfile)
    for name in native.data_vars:
        if "evaluation" in native[name].dims:
            np.testing.assert_array_equal(
                native[name].isel(evaluation=0),
                native[name].isel(evaluation=2),
                err_msg=name,
            )


@pytest.mark.native
def test_native_rayleigh_mirror_and_meridian_u_sign(tmp_path):
    case = native_case()
    black = tmp_path / "black.txt"
    black.write_text("black\n2\n300 0 1 0 0 1.5 0\n2500 0 1 0 0 1.5 0\n")
    case.scene["surface"]["spectrum"] = str(black)
    case.scene["aerosol"]["tau_ref_hi"] = 0
    case.physics["gas_absorption"] = False
    case.numerics["ncolinput"] = 1  # Pure atmospheric single scattering, black bottom.
    case.spectral.update(wavelength_unit="nm", wavelengths=[443.125])
    case.geometry = directions([30, 30, 30, 30], [45, -45, 315, 405])
    result = simulate(case, tmp_path / "rayleigh.nc")
    s = result.Stokes_Vector_TOA.isel(request=0).values
    np.testing.assert_allclose(s[0, :2], s[1, :2], rtol=3e-14)
    np.testing.assert_allclose(s[0, 2:], -s[1, 2:], rtol=3e-14, atol=0)
    np.testing.assert_array_equal(s[1], s[2])
    np.testing.assert_array_equal(s[0], s[3])
    np.testing.assert_array_equal(s[:, 3], 0)  # Rayleigh V is a meaningful finite zero.
    # Independent Cartesian polarization geometry. Rayleigh's polarized light is
    # perpendicular to the scattering plane. Native +Q is along the meridian;
    # +U=2*E_theta*E_phi. Depolarization/column strength cancel from U/Q.
    solar, theta, phi = np.deg2rad([20, 30, 45])
    incoming = [np.sin(solar), 0, -np.cos(solar)]
    outgoing = [np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)]
    e_theta = [np.cos(theta) * np.cos(phi), np.cos(theta) * np.sin(phi), -np.sin(theta)]
    e_phi = [-np.sin(phi), np.cos(phi), 0]
    perpendicular = np.cross(incoming, outgoing)
    a, b = np.dot(perpendicular, e_theta), np.dot(perpendicular, e_phi)
    expected_ratio = 2 * a * b / (a * a - b * b)
    assert s[0, 2] > 0 and s[1, 2] < 0
    np.testing.assert_allclose(s[0, 2] / s[0, 1], expected_ratio, rtol=3e-13)


@pytest.mark.native
@pytest.mark.parametrize("surface", ["land", "ocean"])
@pytest.mark.parametrize("method", ["monochromatic", "response_integration"])
def test_native_directional_analytical_surface(surface, method, tmp_path):
    # Reuse the spectral-quadrature/calibration tests; here only constant illumination
    # and directional Lambertian/black radiance test the new angular boundary.
    case = native_case(surface)
    auxiliary = Path(os.environ["RTSOS_AUXILIARY"])
    for path in auxiliary.iterdir():
        if path.name != "f0.txt":
            (tmp_path / path.name).symlink_to(path)
    (tmp_path / "f0.txt").write_text(
        "# solar\n" * 13 + "".join(f"{w} 100\n" for w in range(300, 2502))
    )
    case.data["auxiliary"] = str(tmp_path)
    case.physics = {"atmosphere": False, "gas_absorption": False}
    case.scene["aerosol"]["tau_ref_hi"] = 0
    if surface == "ocean":
        case.scene["water"] = {"model": "none"}
        case.scene["surface"]["wind_speed_m_s"] = 0
        case.physics["sun_glint"] = False
    case.geometry = directions([37.25, 0, 65, 37.25], [315, 0, 63, -45])
    if method == "monochromatic":
        case.spectral.update(wavelength_unit="nm", wavelengths=[551.25])
    else:
        case.spectral = {"method": method}
        case.sensor = dict(
            response_file=str(EXAMPLES / "response.toml"), channels=["red"]
        )
    result = simulate(case, tmp_path / "surface.nc")
    expected = (0.2 if surface == "land" else 0) * 1000 * np.cos(np.deg2rad(20)) / np.pi
    np.testing.assert_allclose(result.Solar_Irradiance, 1000, rtol=2e-14)
    np.testing.assert_allclose(
        result.Stokes_Vector_TOA.sel(N_Stokes_Parameter=1),
        expected,
        rtol=0,
        atol=1000 * 2e-9,
    )
    np.testing.assert_allclose(
        result.Stokes_Vector_TOA.sel(N_Stokes_Parameter=[2, 3, 4]),
        0,
        rtol=0,
        atol=1000 * 2e-9,
    )
    np.testing.assert_array_equal(result.Tau_Gas_Absorption, 0)
