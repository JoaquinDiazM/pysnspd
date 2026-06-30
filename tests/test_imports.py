"""Basic package import smoke tests."""
from __future__ import annotations


def test_import_package():
    import pysnspd

    assert pysnspd.__version__ == "0.0.0"
