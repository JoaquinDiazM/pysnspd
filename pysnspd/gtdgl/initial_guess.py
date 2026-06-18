"""Analytic initial guesses for stationary gTDGL states."""


def build_uniform_bias_guess(mesh, material_params, bias_config, usadel_catalog):
    """Build a physically motivated initial guess for the SS-run.

    The future implementation should set ``Te = Tph = Tbath``, choose an
    amplitude compatible with the bias state, impose a longitudinal phase
    ramp, and initialize the scalar potential close to zero.
    """
    return 0


def estimate_phase_ramp_from_bias(bias_config, usadel_catalog, geometry_config):
    """Estimate the phase ramp or superfluid momentum associated with ``I_bias``."""
    return 0
