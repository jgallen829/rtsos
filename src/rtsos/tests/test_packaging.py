"""Source-build version extraction for the project's existing release tags."""

from pathlib import Path
import re
import tomllib

import pytest


@pytest.mark.parametrize(
    "tag,version",
    [
        ("v1.0-public-2025-10-23", "1.0"),
        ("v1.1.0-rc1", "1.1.0-rc1"),
        ("v1.1.0-rc2", "1.1.0-rc2"),
    ],
)
def test_release_tag_version(tag, version):
    config = Path(__file__).resolve().parents[3] / "pyproject.toml"
    with config.open("rb") as stream:
        pattern = tomllib.load(stream)["tool"]["setuptools_scm"]["tag_regex"]
    match = re.fullmatch(pattern, tag)
    assert match is not None
    assert match.group("version") == version
