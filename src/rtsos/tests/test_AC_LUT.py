import pytest
from pytest import fixture
import xarray as xr

from rtsos.__main__ import main
from rtsos.parameters import Parameters as P


@fixture(scope="module")
def cache_dir(pytestconfig):
    return pytestconfig.cache.mkdir("ac_lut")


@fixture(scope="module")
def tmp_default(cache_dir):
    path = cache_dir / "default.nc"
    main(["GSFC_AC_LUT", "--pre", str(path)])
    yield path
    path.unlink()


@fixture(scope="module")
def tmp_input(tmp_default):
    path = tmp_default.parent / "input.nc"
    ds = xr.load_dataset(tmp_default)
    ds[P.NWV.__name__][...] = 6
    ds = ds.isel(
        {
            "Wave_Mean_Square_Slope": slice(0, 1),
            "sza-dt": slice(0, 3),  # three sza-dt folders
            "Tau_NIR": slice(0, 2),  # two parameter files in each sza-dt folder
        },
    )
    ds.to_netcdf(path)
    yield path
    path.unlink()


@fixture(scope="module")
def cache_output(tmp_input):
    path = tmp_input.parent / "output.nc"
    parent = path.parent / "output"
    if not any(parent.glob("**/output.nc")):
        # run one AM, one RH, and every sza-dt
        main(
            [
                "GSFC_AC_LUT",
                "--split=Aerosol_Model:2,Relative_Humidity:4,sza-dt",
                str(tmp_input),
                str(path),
            ],
        )
        for item in parent.glob("**/*.infile.txt"):
            item.unlink()
    return path


@fixture
def tmp_split_combine(tmp_input, cache_output):
    # aggregate sza-dt only, not AM or RH
    main(
        [
            "GSFC_AC_LUT",
            "--split=Aerosol_Model:2,Relative_Humidity:4",
            "--combine=sza-dt",
            str(tmp_input),
            str(cache_output),
        ],
    )
    yield cache_output
    path = cache_output.parent / "output"
    for item in path.glob("Aerosol_Model/*/Relative_Humidity/*/output.nc"):
        item.unlink()


@fixture
def tmp_combine(tmp_input, cache_output):
    main(
        [
            "GSFC_AC_LUT",
            "--combine=Aerosol_Model:2,Relative_Humidity:4,sza-dt",
            str(tmp_input),
            str(cache_output),
        ],
    )
    yield cache_output
    cache_output.unlink()


def test_default(tmp_default):
    ds = xr.open_dataset(tmp_default)
    assert len(ds.data_vars) == 0


def test_infile(tmp_input):
    path = tmp_input.parent / "output.nc"
    # generate infiles for one AM, one RH, and every sza-dt
    main(
        [
            "GSFC_AC_LUT",
            "--dry-run",
            "--split=Aerosol_Model:2,Relative_Humidity:4,sza-dt",
            str(tmp_input),
            str(path),
        ],
    )
    infile = tuple(path.parent.glob("**/*.infile.txt"))
    for item in infile:
        item.unlink()
    assert len(infile) == 6


@pytest.mark.native
def test_output(cache_output):
    path = cache_output.parent / "output"
    output = tuple(path.glob("Aerosol_Model/*/Relative_Humidity/*/sza-dt/*/output.nc"))
    assert len(output) == 3
    for item in output:
        dataset = xr.open_dataset(item)
        assert dataset.sizes["Tau_NIR"] == 2


@pytest.mark.native
def test_split_combine(tmp_split_combine):
    path = tmp_split_combine.parent / "output"
    output = tuple(path.glob("Aerosol_Model/*/Relative_Humidity/*/output.nc"))
    assert len(output) == 1
    for item in output:
        dataset = xr.open_dataset(item)
        assert dataset.sizes["Tau_NIR"] == 2
        assert dataset.sizes["Solar_Zenith_Angle"] == 2


@pytest.mark.native
def test_combine(tmp_combine):
    dataset = xr.open_dataset(tmp_combine)
    assert dataset.sizes["Tau_NIR"] == 2
    assert dataset.sizes["Solar_Zenith_Angle"] == 2
