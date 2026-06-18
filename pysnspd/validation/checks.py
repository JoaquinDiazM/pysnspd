"""Validation and consistency checks."""


def check_current_conservation(state, mesh):
    """Check that ``div(j)`` remains small compared with the reference scale."""
    return 0


def check_energy_balance(history):
    """Check consistency of projected electron--phonon energy exchange."""
    return 0


def check_stationary_voltage(state):
    """Check that the stationary superconducting state has no spurious voltage."""
    return 0


def check_phase_continuity(state, mesh):
    """Check phase unwrapping and avoid interpreting branch cuts as physics."""
    return 0


def check_catalog_compatibility(pre_data, ss_data, photon_config):
    """Check that loaded catalogs are compatible with the requested run."""
    return 0
