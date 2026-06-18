"""Generalized TDGL solver placeholders."""


def advance_gtdgl_state(state, mesh, material_fields, dt):
    """Advance the gTDGL fields by one time step."""
    return 0


def run_stationary_gtdgl(initial_state, mesh, catalogs, ss_config):
    """Run gTDGL relaxation toward a stationary pre-photon state."""
    return 0


def run_transient_gtdgl(initial_state, mesh, catalogs, photon_config):
    """Run gTDGL dynamics during the photon-triggered transient."""
    return 0
