"""Basic import smoke tests for the template repository."""


def test_import_package():
    """Import the top-level package."""
    import pysnspd

    assert pysnspd.__version__ == "0.0.0"
