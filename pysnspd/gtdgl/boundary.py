"""Boundary conditions for gTDGL and circuit coupling."""


def build_static_boundary_conditions(mesh, bias_config):
    """Build static longitudinal and transverse boundary conditions."""
    return 0


def build_dynamic_boundary_conditions(mesh, circuit_state, bias_config):
    """Build boundary conditions dynamically coupled to the external circuit."""
    return 0


def compute_terminal_voltage(state, mesh):
    """Compute the terminal voltage from the scalar potential and phase dynamics."""
    return 0
