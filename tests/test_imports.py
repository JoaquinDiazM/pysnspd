"""Basic package import smoke tests."""
from __future__ import annotations

from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:  # Python 3.10, the minimum supported version.
    import tomli as tomllib


def test_import_package():
    import pysnspd

    project_file = Path(__file__).resolve().parents[1] / "pyproject.toml"
    with project_file.open("rb") as stream:
        project = tomllib.load(stream)["project"]
    assert pysnspd.__version__ == project["version"]
