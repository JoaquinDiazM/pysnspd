from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from pysnspd.gtdgl.solver import TDGLSolver
from pysnspd.plotting.ss_run import plot_ss_adaptive_timestep_history


def test_adaptive_euler_records_retry_shrink_diagnostics():
    solver = TDGLSolver.__new__(TDGLSolver)
    solver.options = SimpleNamespace(
        adaptive=True,
        max_solve_retries=3,
        adaptive_time_step_multiplier=0.25,
    )
    solver.allmaras_forcing_callback = None
    solver.operators = SimpleNamespace(psi_laplacian=object())
    solver.device = SimpleNamespace(layer=SimpleNamespace(gamma=0.0, u=1.0))
    solver.terminal_psi_value = None
    solver.normal_boundary_index = np.array([], dtype=np.int64)
    solver.apply_terminal_psi = lambda psi: psi

    attempted_dt = []

    def fake_solve_for_psi_squared(**kwargs):
        attempted_dt.append(float(kwargs["dt"]))
        if len(attempted_dt) == 1:
            return None
        psi = np.asarray(kwargs["psi"], dtype=np.complex128)
        return psi, np.abs(psi) ** 2

    solver.solve_for_psi_squared = fake_solve_for_psi_squared

    psi0 = np.array([1.0 + 0.0j, 0.5 + 0.1j], dtype=np.complex128)
    abs_sq0 = np.abs(psi0) ** 2
    mu0 = np.zeros(psi0.size, dtype=float)
    epsilon = np.ones(psi0.size, dtype=float)

    psi1, abs_sq1, dt = TDGLSolver.adaptive_euler_step(
        solver,
        step=7,
        psi=psi0,
        abs_sq_psi=abs_sq0,
        mu=mu0,
        epsilon=epsilon,
        dt=1.0,
    )

    assert np.allclose(psi1, psi0)
    assert np.allclose(abs_sq1, abs_sq0)
    assert dt == 0.25
    assert attempted_dt == [1.0, 0.25]
    assert solver.last_adaptive_dt_attempt == 1.0
    assert solver.last_adaptive_dt_accepted == 0.25
    assert solver.last_adaptive_retries == 1
    assert solver.last_adaptive_rejected_attempts == 1


def test_adaptive_timestep_plot_smoke(tmp_path):
    t_s = np.linspace(0.0, 2.0e-12, 5)
    history = {
        "t_s": t_s,
        "dt_s": np.array([0.2, 0.2, 0.4, 0.8, 1.0]) * 1.0e-15,
        "dt_attempt_s": np.array([0.2, 0.2, 0.4, 0.8, 1.0]) * 1.0e-15,
        "dt_next_s": np.array([0.2, 0.4, 0.8, 1.0, 1.0]) * 1.0e-15,
        "adaptive_target_dt_s": np.array([np.nan, np.nan, 0.9, 1.0, 1.1]) * 1.0e-15,
        "adaptive_retries": np.array([0, 1, 0, 0, 0]),
        "adaptive_rejected_attempts": np.array([0, 1, 0, 0, 0]),
        "adaptive_window_mean_d_abs_sq": np.array([np.nan, 2.0e-2, 1.0e-2, 8.0e-3, 7.0e-3]),
    }
    out = plot_ss_adaptive_timestep_history(history, tmp_path / "adaptive.png", dpi=80)
    assert out.exists()
    assert out.stat().st_size > 0
