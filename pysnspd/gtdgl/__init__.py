"""Mesoscopic gTDGL package for pySNSPD.

The active backend is the pyTDGL-like stationary/transient core promoted to the
package root.  The old legacy OE7 solver modules were removed from this package
slice.  Shared SI material, state, seed and finite-volume diagnostic helpers
remain at package level because the promoted backend still uses them.
"""
from __future__ import annotations

from .material import GTDGLMaterial, build_gtdgl_material
from .operators import FVOperators, build_fv_operators
from .options import SolverOptions, SolverOptionsError, SparseSolver
from .solver import SolverResult, TDGLSolver, validate_terminal_currents
from .adapter import solve_stationary_pytdgl_like
from .allmaras import AllmarasForcingFields, compute_allmaras_forcing_dimensionless
from .state import CurrentFields, GTDGLStationaryState, RelaxationResult
from .ss_targets import (
    apply_terminal_proximity_seed,
    contact_recovery_diagnostics,
    continuity_diagnostics,
    stationarity_diagnostics,
)

__all__ = [
    "GTDGLMaterial",
    "build_gtdgl_material",
    "FVOperators",
    "build_fv_operators",
    "CurrentFields",
    "GTDGLStationaryState",
    "RelaxationResult",
    "SolverOptions",
    "SolverOptionsError",
    "SolverResult",
    "SparseSolver",
    "TDGLSolver",
    "validate_terminal_currents",
    "solve_stationary_pytdgl_like",
    "AllmarasForcingFields",
    "compute_allmaras_forcing_dimensionless",
    "apply_terminal_proximity_seed",
    "contact_recovery_diagnostics",
    "continuity_diagnostics",
    "stationarity_diagnostics",
]
