"""Configuration interface for pySNSPD.

This module will eventually validate the full project configuration:
paths, units, material parameters, mesh parameters, catalog resolution,
stationary-run parameters, photon-run parameters, and circuit parameters.
"""


def load_config(config_path):
    """Load a project configuration file.

    Parameters
    ----------
    config_path : str or pathlib.Path
        Path to the project configuration file.

    Returns
    -------
    int
        Placeholder return value. Future implementation should return a
        structured configuration object.
    """
    return 0


def validate_config(config):
    """Validate required configuration fields and physical units.

    This should enforce that ``big_data_root`` is explicitly configured and
    that no heavy raw-data output is written inside the Python package.
    """
    return 0


def summarize_config(config):
    """Create a compact human-readable summary of the configuration."""
    return 0
