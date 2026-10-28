import pytest
from pytest import fixture

from rtsos.__main__ import main


@fixture(scope="session")
def phmx_cache_dir(pytestconfig):
    return pytestconfig.cache.mkdir("phmx")


@fixture(scope="session")
def phmx_tmp_default(phmx_cache_dir):
    path = phmx_cache_dir / "default.nc"
    main(["Aerosol_Phmx_Cal", "--pre", str(path)])
    yield path
    path.unlink()


@fixture(scope="session")
def phmx_cache_output(phmx_tmp_default):
    path = phmx_tmp_default.parent / "output.nc"
    if not path.exists():
        main(
            [
                "Aerosol_Phmx_Cal",
                str(phmx_tmp_default),
                str(path),
            ],
        )
        for item in path.parent.glob("**/*.infile.txt"):
            item.unlink()
    return path


def pytest_addoption(parser):
    parser.addoption(
        "--native",
        action="store_true",
        help="run native integration tests (requires RTSOS executables and auxiliary data)",
    )


def pytest_collection_modifyitems(config, items):
    if not config.getoption("--native"):
        skip = pytest.mark.skip(
            reason="native integration requires --native and auxiliary data"
        )
        for item in items:
            if "native" in item.keywords:
                item.add_marker(skip)
