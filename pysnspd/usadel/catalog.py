"""Catalog construction for DOS, gap, current, and spectral quantities."""


def build_dos_catalog(material_config, bias_config, catalog_config, parallel_config):
    """Build a fine catalog of ``rho(E; |Delta|, q)`` values.

    This is part of the PRE-run and should be parallelizable.
    """
    return 0


def build_current_catalog(material_config, bias_config, catalog_config, parallel_config):
    """Build a catalog relating bias current, superfluid momentum, and gap."""
    return 0


def save_usadel_catalog(catalog, path):
    """Save a Usadel catalog to disk."""
    return 0


def load_usadel_catalog(path):
    """Load a Usadel catalog from disk."""
    return 0
