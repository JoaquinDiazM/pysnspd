"""Phase-space integrals for projected electron--phonon powers."""


def compute_J_scattering(omega, Te, delta, q, dos_catalog):
    """Compute the scattering phase-space functional ``mathcal{J}_S``."""
    return 0


def compute_J_recombination(omega, Te, delta, q, dos_catalog):
    """Compute the recombination / pair-breaking functional ``mathcal{J}_R``."""
    return 0


def build_phase_space_catalog(usadel_catalog, kinetic_config, parallel_config):
    """Build fine catalogs for ``mathcal{J}_S`` and ``mathcal{J}_R``.

    This is part of the PRE-run and should be parallelizable independently
    of the gTDGL solver.
    """
    return 0


def save_phase_space_catalog(catalog, path):
    """Save phase-space catalogs to disk."""
    return 0


def load_phase_space_catalog(path):
    """Load phase-space catalogs from disk."""
    return 0
