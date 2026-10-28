"""Physical cases through the real dispatcher, serializer, and result reader."""

from copy import deepcopy
from dataclasses import asdict
from pathlib import Path
import subprocess
import tomllib

import numpy as np
import pytest
import xarray as xr

from rtsos import Case, load_case, simulate
from rtsos.__main__ import main
from rtsos.PACE_Simulator_DoubleK import case_inputs, init, solar_source

EXAMPLES = Path(__file__).resolve().parents[3] / "examples"


def ocean():
    return load_case(EXAMPLES / "pace_ocean.toml")


@pytest.mark.parametrize(
    "surface,water,code",
    [
        ("ross_li", None, -203),
        ("snow", None, -202),
        ("mrpv", None, -201),
        ("lambertian_fresnel", None, -200),
        ("ocean", "none", 0),
        ("ocean", "chlorophyll", 1),
        ("ocean", "chlorophyll_sediment", 2),
        ("ocean", "seven_parameter", 3),
    ],
)
def test_named_models_native_records(surface, water, code, tmp_path):
    case = ocean()
    case.scene["surface"] = {"model": surface}
    if water is None:
        case.scene.pop("water")
        case.scene["surface"]["spectrum"] = str(EXAMPLES / "land_spectrum.txt")
    else:
        case.scene["water"] = {"model": water}
        if code:
            case.scene["water"].update(chla=0.7, water_depth_max=180)
        if code == 2:
            case.scene["water"].update(
                phytoplankton_index_refraction=1.03,
                phytoplankton_spectral_slope=3.2,
                sediment_index_refraction=1.15,
                sediment_spectral_slope=4.2,
                sediment_concentration=0.0,
            )
        if code == 3:
            case.scene["water"].update(
                adg440=0.0,
                bbp660_BackscatterCoeff=0.02,
                Bp660_BackscatterFraction=0.015,
                Sdg=0.017,
                Sbp=0.2,
                S_Bp=-0.1,
            )
    case.physics = {"scattering": "vector"}
    inputs = case_inputs(case)
    assert not inputs.dims  # Exactly one physical state, including land.
    program = init()
    infile, _ = program.infile(tmp_path, inputs)
    rows = [
        r.split("#", 1)[0].strip() for r in (tmp_path / infile).read_text().splitlines()
    ]
    assert int(rows[12]) == code
    assert len(rows) == (37 if code < 0 else 58)
    start = 13 if code < 0 else 27
    assert [int(v) for v in rows[start : start + 6]] == [8, 12, 24, 3, 2, 2]
    if code < 0:
        assert rows[-2] == str(EXAMPLES / "land_spectrum.txt")
    if code > 0:
        assert float(rows[15]) == float(np.float32(0.7))
    if code == 2:
        np.testing.assert_array_equal(
            [float(v) for v in rows[16:21]],
            np.array([1.03, 3.2, 1.15, 4.2, 0], dtype=np.float32).astype(float),
        )
    if code == 3:
        np.testing.assert_array_equal(
            [float(v) for v in rows[21:27]],
            np.array([0, 0.02, 0.015, 0.017, 0.2, -0.1], dtype=np.float32).astype(
                float
            ),
        )
        assert inputs.Sbp.attrs["units"] == inputs.S_Bp.attrs["units"] == "1"


def test_case_path_roundtrip(tmp_path, monkeypatch):
    work = tmp_path / "cases with spaces"
    work.mkdir()
    source = (
        (EXAMPLES / "pace_land.toml")
        .read_text()
        .replace("land_spectrum.txt", "spectra/a longer spectrum filename.txt")
    )
    path = work / "land.toml"
    path.write_text(source)
    case = load_case(path)
    mapping = asdict(case)
    monkeypatch.chdir(tmp_path)
    from_mapping = Case.from_mapping(mapping, base_dir=tmp_path)
    assert from_mapping == case
    inputs = case_inputs(case)
    inputs.to_netcdf(tmp_path / "inputs.nc")
    loaded = xr.load_dataset(tmp_path / "inputs.nc")
    xr.testing.assert_identical(inputs, loaded)
    infile, _ = init().infile(tmp_path, loaded)
    rows = (tmp_path / infile).read_text().splitlines()
    assert str(work / "spectra/a longer spectrum filename.txt") in rows[-2]
    assert str(tmp_path / "RTSOS_data/Data") in rows[-7]
    assert case.scene["atmosphere"]["profile"] == "afglus.dat"


def test_python_and_cli_same_case(tmp_path, monkeypatch):
    case_file = tmp_path / "ocean.toml"
    case_file.write_text((EXAMPLES / "pace_ocean.toml").read_text())
    case = load_case(case_file)
    auxiliary = tmp_path / "aux"
    auxiliary.mkdir()
    (auxiliary / "f0.txt").touch()
    case.data = {"auxiliary": str(auxiliary), "gas_absorption": str(tmp_path / "gas")}
    # Put the same absolute data paths in the file used by the CLI.
    case_file.write_text(
        case_file.read_text().replace(
            'root = "../RTSOS_data"',
            f'auxiliary = "{auxiliary}"\ngas_absorption = "{tmp_path / "gas"}"',
        )
    )
    records = []

    def native(args, check):
        rows = Path(args[1]).read_text().splitlines()
        records.append(rows[:-1])
        assert check
        xr.Dataset(
            {
                "Stokes_Vector_TOA": ("Stokes", [0.1, -0.02, 0.0, 0.0]),
                "MONOCHROMATIC_FLAG": np.int32(0),
            }
        ).to_netcdf(rows[-1] + ".h5")

    monkeypatch.setattr(subprocess, "run", native)
    result = simulate(case, tmp_path / "python.nc", tmp=tmp_path)
    main(["simulate", str(case_file), str(tmp_path / "cli.nc"), "--tmp", str(tmp_path)])
    xr.testing.assert_identical(result, xr.load_dataset(tmp_path / "cli.nc"))
    assert records[0] == records[1]
    assert result.Stokes_Vector_TOA.values.tolist() == [0.1, -0.02, 0.0, 0.0]
    assert result.attrs["illumination_source"] == str(auxiliary / "f0.txt")
    assert result.attrs["pace_rt_fourier_order_limit"] == 3
    main(["simulate", str(case_file), str(tmp_path / "dry.nc"), "--dry-run"])
    assert len(records) == 2
    assert not (tmp_path / "dry.nc").exists()


@pytest.mark.parametrize("method", ["monochromatic", "response_integration"])
def test_unconnected_spectral_operation(method):
    case = ocean()
    case.spectral = {"method": method, "wavelength_unit": "nm", "wavelengths": [440.0]}
    with pytest.raises(NotImplementedError, match="requested.*wavelengths"):
        case_inputs(case)


def test_custom_selection_is_separate_from_execution():
    case = load_case(EXAMPLES / "custom_response.toml")
    responses = case.load_responses()
    assert [c.id for c in responses.channels] == ["red", "blue"]
    assert responses.channels[1].response.tolist() == [0, -0.001, 1, 0.3, 0]
    with pytest.raises(NotImplementedError, match="custom response integration"):
        case_inputs(case)
    case.spectral = {"method": "pace_native"}
    with pytest.raises(NotImplementedError, match="Custom sensors/responses"):
        case_inputs(case)


def test_cli_reports_missing_response_case_data(tmp_path, capsys):
    path = tmp_path / "case.toml"
    path.write_text(
        (EXAMPLES / "custom_response.toml")
        .read_text()
        .replace('root = "../RTSOS_data"', 'root = "missing"')
        .replace(
            'response_file = "response.toml"',
            f'response_file = "{EXAMPLES / "response.toml"}"',
        )
    )
    with pytest.raises(SystemExit) as error:
        main(
            [
                "simulate",
                str(path),
                str(tmp_path / "out.nc"),
            ]
        )
    assert error.value.code == 2
    assert "f0.txt" in capsys.readouterr().err
    assert not (tmp_path / "out.nc").exists()


@pytest.mark.parametrize(
    "change,error",
    [
        ({"zenith_count": 3}, "middle angle"),
        ({"sensor_height_km": 100}, "reset"),
        ({"zenith_count": 2.5}, "integer"),
    ],
)
def test_native_geometry_boundary(change, error):
    case = ocean()
    case.geometry.update(change)
    with pytest.raises(ValueError, match=error):
        case_inputs(case)


def test_operating_settings_scale_and_fractional_counts_rejected():
    case = ocean()
    case.numerics.update(ncolinput=128, nquadainput=256, MAXMORDINPUT=80)
    case.geometry.update(zenith_count=400, azimuth_count=500)
    inputs = case_inputs(case)
    assert inputs.nquadainput.item() == 256
    assert inputs.MAXMORDINPUT.item() == 80
    assert inputs.attrs["pace_rt_fourier_order_limit"] == 30
    assert inputs.NPHIV.item() == 500
    case.numerics["ncolinput"] = 8.5
    with pytest.raises(ValueError, match="integer"):
        case_inputs(case)


@pytest.mark.parametrize("quadrature,order,bound", [(12, 80, 12), (256, 5, 5)])
def test_native_fourier_bound_is_reported(quadrature, order, bound):
    case = ocean()
    case.numerics.update(nquadainput=quadrature, MAXMORDINPUT=order)
    inputs = case_inputs(case)
    assert inputs.MAXMORDINPUT.item() == order
    assert inputs.attrs["pace_rt_fourier_order_limit"] == bound


@pytest.mark.parametrize("section", ["surface", "water", "aerosol", "atmosphere"])
def test_malformed_scene_table(section):
    mapping = asdict(ocean())
    mapping["scene"][section] = "invalid table"
    with pytest.raises(ValueError, match=rf"scene\.{section} must be a table"):
        Case.from_mapping(mapping, base_dir=EXAMPLES)


def test_malformed_case_path():
    mapping = asdict(ocean())
    mapping["data"]["root"] = 123
    with pytest.raises(ValueError, match="root must be a filesystem path"):
        Case.from_mapping(mapping, base_dir=EXAMPLES)


def test_reject_irrelevant_water_inputs_and_missing_sweep_choice():
    case = ocean()
    case.scene["water"]["sediment_concentration"] = 1.0
    with pytest.raises(ValueError, match="sediment_concentration"):
        case_inputs(case)
    case = ocean()
    del case.scene["water"]["chla"]
    with pytest.raises(ValueError, match="requires chla"):
        case_inputs(case)


def test_illumination_precedence(tmp_path):
    case = ocean()
    case.data = {"auxiliary": str(tmp_path), "gas_absorption": str(tmp_path)}
    (tmp_path / "f0.txt").touch()
    assert solar_source(case).name == "f0.txt"
    case.illumination["source"] = "tsis"
    with pytest.raises(ValueError, match="native solar reader selects"):
        solar_source(case)
    (tmp_path / "tsis1_hsrs_binned_fs.csv").touch()
    assert solar_source(case).name == "tsis1_hsrs_binned_fs.csv"
    case.illumination["source"] = "thuillier"
    with pytest.raises(ValueError, match="native solar reader selects"):
        solar_source(case)


def test_mapping_does_not_change_input():
    with (EXAMPLES / "pace_ocean.toml").open("rb") as stream:
        mapping = tomllib.load(stream)
    before = deepcopy(mapping)
    Case.from_mapping(mapping, base_dir=EXAMPLES)
    assert mapping == before
