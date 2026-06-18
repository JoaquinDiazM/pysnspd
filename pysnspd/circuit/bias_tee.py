"""Minimal external circuit placeholders."""


def initialize_circuit_state(circuit_config, bias_config):
    """Initialize the bias tee and readout variables before a run."""
    return 0


def advance_circuit_state(circuit_state, V_tdgl, dt, circuit_config):
    """Advance the external circuit state from the internal detector voltage."""
    return 0


def compute_output_voltage(circuit_state, circuit_config):
    """Compute the observable output voltage across the load."""
    return 0
