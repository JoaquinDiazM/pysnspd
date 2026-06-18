"""Projected microscopic powers for the thermal model."""


def compute_scattering_power(Te, Tph, delta, q, phase_space_catalog, alpha2F):
    """Compute ``P_ep_S`` from the projected Simon scattering channel."""
    return 0


def compute_recombination_power(Te, Tph, delta, q, phase_space_catalog, alpha2F):
    """Compute ``P_ep_R`` from the projected recombination / pair-breaking channel."""
    return 0


def compute_escape_power(Tph, Tbath, phonon_dos, tau_esc):
    """Compute the phonon escape power ``P_esc`` to the substrate."""
    return 0


def compute_total_electron_phonon_power(Te, Tph, delta, q, catalogs):
    """Compute the net projected electron--phonon power used by the thermal solver."""
    return 0
