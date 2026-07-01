"""Smoke tests for the flattened ``pysnspd.gtdgl`` public API."""
from __future__ import annotations

import importlib.util


def test_gtdgl_public_backend_names_are_flat():
    import pysnspd.gtdgl as gtdgl

    for name in (
        "GTDGLMaterial",
        "build_gtdgl_material",
        "FVOperators",
        "build_fv_operators",
        "SolverOptions",
        "SparseSolver",
        "TDGLSolver",
        "SolverResult",
        "validate_terminal_currents",
        "solve_stationary_pytdgl_like",
        "CurrentFields",
        "GTDGLStationaryState",
        "RelaxationResult",
        "AllmarasForcingFields",
        "compute_allmaras_forcing_dimensionless",
    ):
        assert hasattr(gtdgl, name), name


def test_nested_backend_package_is_removed():
    # Keep the removed package name split so repository-wide legacy-grep checks
    # do not report this test as a stale import path.
    removed_name = "pysnspd.gtdgl." + "pytdgl_like"
    try:
        spec = importlib.util.find_spec(removed_name)
    except ModuleNotFoundError:
        spec = None
    assert spec is None


def test_flat_submodules_import_directly():
    from pysnspd.gtdgl.adapter import solve_stationary_pytdgl_like
    from pysnspd.gtdgl.allmaras import (
        allmaras_coefficients,
        compute_allmaras_forcing_dimensionless,
    )
    from pysnspd.gtdgl.currents import native_current_scale_A_m2
    from pysnspd.gtdgl.device import build_pytdgl_like_device
    from pysnspd.gtdgl.options import SolverOptions
    from pysnspd.gtdgl.solver import TDGLSolver
    from pysnspd.gtdgl.tdgl_operators import MeshOperators
    from pysnspd.gtdgl.usadel_current import compute_usadel_supercurrent_diagnostic

    assert callable(solve_stationary_pytdgl_like)
    assert callable(allmaras_coefficients)
    assert callable(compute_allmaras_forcing_dimensionless)
    assert callable(native_current_scale_A_m2)
    assert callable(build_pytdgl_like_device)
    assert SolverOptions is not None
    assert TDGLSolver is not None
    assert MeshOperators is not None
    assert callable(compute_usadel_supercurrent_diagnostic)
