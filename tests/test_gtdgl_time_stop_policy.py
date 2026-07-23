from __future__ import annotations

import inspect

from pysnspd.solver.stationary import solve_stationary_pytdgl_like
from pysnspd.solver.core import TDGLSolver


def test_solver_exposes_explicit_stop_on_convergence_policy():
    sig = inspect.signature(TDGLSolver)
    assert "stop_on_convergence" in sig.parameters
    assert sig.parameters["stop_on_convergence"].default is False


def test_adapter_runs_to_time_by_default_not_eta_stop():
    sig = inspect.signature(solve_stationary_pytdgl_like)
    assert "stop_on_convergence" in sig.parameters
    assert sig.parameters["stop_on_convergence"].default is False
