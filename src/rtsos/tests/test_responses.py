from pathlib import Path

from netCDF4 import Dataset
import numpy as np
import pytest

from rtsos.responses import ChannelResponse, load_oci_rsr, load_responses


def test_exact_toml_samples_and_variant_selection():
    path = Path(__file__).resolve().parents[3] / "examples/response.toml"
    responses = load_responses(path)
    selected = responses.select(["red", "blue"])
    assert [(c.id, c.mirror) for c in selected.channels] == [
        ("red", 0),
        ("blue", 0),
        ("blue", 1),
    ]
    for channel, original in zip(
        selected.channels, [responses.channels[i] for i in [2, 0, 1]]
    ):
        np.testing.assert_array_equal(channel.wavelength, original.wavelength)
        np.testing.assert_array_equal(channel.response, original.response)
    blue = responses.select(["blue"], mirror=1).channels[0]
    assert blue.wavelength.tolist() == [389, 431, 440.2, 451, 511]
    assert blue.response.tolist() == [0, -0.002, 0.9, 0.4, 0]
    assert blue.nominal_label == "blue 440"
    assert blue.calibration_center_wavelength is None
    assert selected.channels[0].wavelength_unit == "um"
    with pytest.raises(ValueError, match="No supplied calibration"):
        responses.select(["red"], mirror=1)


def test_exact_oci_source_arrays(tmp_path):
    path = tmp_path / "rsr.nc"
    wave = np.array([320, 400, 441.25, 900, 2500], dtype=np.float32)
    weights = np.array([0, -1e-6, 0.8, 1e-7, 0], dtype=np.float32)
    with Dataset(path, "w") as source:
        for detector in ["blue", "red", "SWIR"]:
            group = source.createGroup(detector)
            for dim, count in [("band", 2), ("sample", 5), ("mirror", 2)]:
                group.createDimension(dim, count)
            w = group.createVariable("sampled_wavelength", "f4", ("sample",))
            w.units = "nm"
            w[:] = wave
            group.createVariable("center_wavelength", "f4", ("band", "mirror"))[:] = [
                [441.5, 441.75],
                [442.5, 442.75],
            ]
            r = group.createVariable("RSR", "f4", ("band", "sample", "mirror"))
            for index in range(2):
                for mirror in range(2):
                    r[index, :, mirror] = weights * (index + 1 + mirror)
    responses = load_oci_rsr(path)
    assert len(responses.channels) == 12
    chosen = responses.select(["oci_swir_002", "oci_blue_001"], mirror=1)
    assert [c.detector for c in chosen.channels] == ["swir", "blue"]
    for channel, factor in zip(chosen.channels, [3, 2]):
        np.testing.assert_array_equal(channel.wavelength, wave)
        np.testing.assert_array_equal(channel.response, weights * factor)
        assert channel.response.dtype == weights.dtype
    assert chosen.channels[1].calibration_center_wavelength == 441.75
    assert chosen.channels[1].nominal_label == ""
    assert chosen.channels[1].wavelength[2] == 441.25
    with pytest.raises(ValueError, match="No supplied calibration"):
        responses.select(["oci_blue_001"], mirror=2)


def test_explicit_units():
    with pytest.raises(ValueError, match="wavelength_unit"):
        ChannelResponse("x", [400, 500], [0, 1], "")


def test_oci_masks_and_missing_center_preserved(tmp_path):
    path = tmp_path / "masked_rsr.nc"
    wave = np.ma.array(
        np.array([400, 450, 500], dtype=np.float32), mask=[False, True, False]
    )
    weights = np.ma.array(
        np.array([-0.0, -0.01, 1.0], dtype=np.float32), mask=[False, True, False]
    )
    with Dataset(path, "w") as source:
        for detector in ["blue", "red", "SWIR"]:
            group = source.createGroup(detector)
            for dim, count in [("band", 1), ("sample", 3), ("mirror", 1)]:
                group.createDimension(dim, count)
            w = group.createVariable("sampled_wavelength", "f4", ("sample",))
            w.units = "nm"
            w[:] = wave
            group.createVariable("center_wavelength", "f4", ("band", "mirror"))
            r = group.createVariable("RSR", "f4", ("band", "sample", "mirror"))
            r[0, :, 0] = weights
    channel = load_oci_rsr(path).select(["oci_blue_001"]).channels[0]
    with Dataset(path) as source:
        for actual, variable in [
            (channel.wavelength, source["blue/sampled_wavelength"][:]),
            (channel.response, source["blue/RSR"][0, :, 0]),
        ]:
            assert actual.dtype == variable.dtype
            np.testing.assert_array_equal(actual.data, variable.data)
            np.testing.assert_array_equal(actual.mask, variable.mask)
    assert channel.calibration_center_wavelength is None
    assert np.signbit(channel.response.data[0])


def test_programmatic_mask_and_signed_zero_preserved():
    response = np.ma.array([-0.0, -0.01, 1.0], mask=[False, True, False])
    channel = ChannelResponse("x", [400, 450, 500], response, "nm")
    np.testing.assert_array_equal(channel.response.data, response.data)
    np.testing.assert_array_equal(channel.response.mask, response.mask)
    assert np.signbit(channel.response.data[0])
