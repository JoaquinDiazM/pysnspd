from __future__ import annotations

import inspect

from pysnspd.solver.stationary import solve_stationary_pytdgl_like
from pysnspd.solver.core import TDGLSolver
from pysnspd.solver.steady_gate import classify_photon_readiness


def test_solver_removed_eta_residual_stop_policy():
    sig = inspect.signature(TDGLSolver)
    assert "stop_on_convergence" not in sig.parameters
    assert "stop_eta" not in sig.parameters
    assert "stop_min_steps" not in sig.parameters


def test_adapter_exposes_persistent_photon_ready_policy():
    sig = inspect.signature(solve_stationary_pytdgl_like)
    assert "stop_on_convergence" not in sig.parameters
    assert "stationarity_eta" not in sig.parameters
    assert sig.parameters["photon_ready_consecutive_evaluations"].default == 5
    assert sig.parameters["dynamic_stationarity_minimum_tail_ps"].default == 5.0


def test_photon_ready_accepts_dynamic_attractor_as_alternative_to_fixed_point():
    result = classify_photon_readiness(
        strict_stationarity_passes=False,
        dynamic_stationarity_passes=True,
        contact_recovery_passes=True,
        continuity_passes=True,
        thermal_stationarity_passes=True,
        circuit_stationarity_passes=True,
        phase_drive_converged=True,
    )
    assert result == {
        "passes": True,
        "mesoscopic_passes": True,
        "mesoscopic_mode": "weak_dynamic_attractor",
    }


def test_photon_ready_keeps_hard_validity_gates_mandatory():
    result = classify_photon_readiness(
        strict_stationarity_passes=True,
        dynamic_stationarity_passes=False,
        contact_recovery_passes=True,
        continuity_passes=False,
        thermal_stationarity_passes=True,
        circuit_stationarity_passes=True,
        phase_drive_converged=True,
    )
    assert result["passes"] is False
    assert result["mesoscopic_mode"] == "strict_fixed_point"
