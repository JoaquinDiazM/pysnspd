"""Solver option validation for the flattened gTDGL backend."""
from __future__ import annotations

import pytest

from pysnspd.solver.options import SolverOptions, SolverOptionsError, SparseSolver
from pysnspd.solver.core import TDGLSolver, validate_terminal_currents


def test_solver_options_default_superlu_validates():
    opts = SolverOptions(solve_time=1.0e-3, dt_init=1.0e-4, dt_max=1.0e-4)
    opts.validate()
    assert opts.sparse_solver is SparseSolver.SUPERLU


def test_solver_options_reject_bad_dt():
    opts = SolverOptions(solve_time=1.0e-3, dt_init=2.0e-4, dt_max=1.0e-4)
    with pytest.raises(SolverOptionsError):
        opts.validate()


def test_terminal_current_validation_enforces_zero_sum():
    opts = SolverOptions(solve_time=1.0e-3)
    terminals = [type("T", (), {"name": "left"})(), type("T", (), {"name": "right"})()]
    validate_terminal_currents({"left": -1.0e-6, "right": 1.0e-6}, terminals, opts)
    with pytest.raises(ValueError):
        validate_terminal_currents({"left": 1.0e-6, "right": 1.0e-6}, terminals, opts)


def test_solver_class_keeps_pytdgl_like_core_methods():
    for method in (
        "solve_for_psi_squared",
        "adaptive_euler_step",
        "solve_for_observables",
        "update",
        "solve",
    ):
        assert hasattr(TDGLSolver, method), method
