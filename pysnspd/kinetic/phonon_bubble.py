"""Photon-energy helpers for the future phonon-bubble layer.

Only model-independent photon-energy conversion is provided here. The spatial
phonon-bubble profile and thermal injection rule are future OE targets and
should be implemented only after the coupled no-photon stationary state is
validated.
"""

from __future__ import annotations

from scipy import constants


def photon_energy_from_wavelength(wavelength_m: float) -> float:
    """Return photon energy in joules from wavelength in meters.

    Parameters
    ----------
    wavelength_m:
        Photon wavelength in SI meters.

    Returns
    -------
    float
        Photon energy ``h c / wavelength_m`` in joules.
    """
    wavelength = float(wavelength_m)
    if wavelength <= 0.0:
        raise ValueError("wavelength_m must be positive.")
    return float(constants.h * constants.c / wavelength)


__all__ = ["photon_energy_from_wavelength"]