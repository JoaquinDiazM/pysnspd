"""Effective two-temperature thermal projection."""


def electronic_energy_density(Te, delta, q, dos_catalog):
    """Evaluate the superconducting electronic energy density used on the LHS."""
    return 0


def phonon_energy_density(Tph, phonon_dos):
    """Evaluate the projected phonon energy density."""
    return 0


def electronic_heat_capacity(Te, delta, q, dos_catalog):
    """Evaluate or interpolate the effective electronic heat capacity."""
    return 0


def phonon_heat_capacity(Tph, phonon_dos):
    """Evaluate or interpolate the effective phonon heat capacity."""
    return 0


def advance_thermal_state(thermal_state, powers, dt):
    """Advance ``Te`` and ``Tph`` by one thermal substep."""
    return 0
