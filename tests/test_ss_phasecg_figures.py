from __future__ import annotations

import numpy as np

from pysnspd.plotting.ss_phasecg_figures import make_phasecg_ss_figures


def test_make_phasecg_ss_figures(tmp_path):
    nodes_x_nm = np.array([0.0, 10.0, 0.0, 10.0])
    nodes_y_nm = np.array([0.0, 0.0, 10.0, 10.0])
    triangles = np.array([[0, 1, 2], [1, 3, 2]], dtype=np.int64)
    snapshot_t_ps = np.array([0.0, 1.0])
    history_t_ps = np.linspace(0.0, 1.0, 5)
    snapshot_field = np.array(
        [[0.8, 0.9, 1.0, 0.85], [0.7, 0.8, 0.9, 0.75]],
        dtype=float,
    )

    dataset = {
        "run_name": "phasecg_unit",
        "snapshot_t_ps": snapshot_t_ps,
        "t_ps": history_t_ps,
        "nodes_x_nm": nodes_x_nm,
        "nodes_y_nm": nodes_y_nm,
        "triangles": triangles,
        "delta_snapshot_over_delta0": snapshot_field,
        "phi_snapshot_mV": snapshot_field - 0.8,
        "qxi_snapshot": snapshot_field,
        "js_snapshot_over_javg": snapshot_field,
        "jn_snapshot_over_javg": 1.0 - snapshot_field,
        "jtot_snapshot_over_javg": np.ones_like(snapshot_field),
        "div_j_snapshot_normalized": 1.0e-10 * (snapshot_field - 0.8),
        "target_current_uA": 43.0,
        "current_total_snapshot_uA": np.array([43.0, 43.0]),
        "current_super_snapshot_uA": np.array([35.0, 34.0]),
        "current_normal_snapshot_uA": np.array([8.0, 9.0]),
        "terminal_voltage_mV": np.linspace(0.0, 2.0, 5),
        "voltage_center_snapshot_mV": np.array([0.0, 0.5]),
        "voltage_terminal_snapshot_mV": np.array([0.0, 2.0]),
        "delta_center_min": np.array([0.8, 0.7]),
        "delta_center_mean": np.array([0.9, 0.8]),
        "delta_center_max": np.array([1.0, 0.9]),
        "normal_current_fraction_snapshot": np.array([8.0 / 43.0, 9.0 / 43.0]),
        "dt_attempt_fs": np.linspace(0.1, 0.2, 5),
        "dt_accepted_fs": np.linspace(0.08, 0.15, 5),
        "dt_next_fs": np.linspace(0.1, 0.2, 5),
        "solve_attempts_per_step": np.array([1.0, 2.0, 1.0, 3.0, 1.0]),
        "adaptive_retries": np.array([0.0, 1.0, 0.0, 2.0, 0.0]),
        "cumulative_rejected_attempts": np.array([0.0, 1.0, 1.0, 3.0, 3.0]),
        "estimated_wall_step_s": np.ones(5) * 0.2,
        "estimated_wall_cumulative_s": np.linspace(0.2, 1.0, 5),
        "measured_wall_time_s": 1.0,
        "eta_R": np.logspace(-2, -3, 5),
        "allmaras_update_forcing_max_abs": np.logspace(1, 2, 5),
        "poisson_residual_rel": np.logspace(-12, -11, 5),
        "poisson_tolerance": 1.0e-9,
        "div_j_normalized_max_snapshot": np.array([1.0e-10, 2.0e-10]),
        "div_j_normalized_rms_snapshot": np.array([5.0e-11, 1.0e-10]),
        "allmaras_phase_convergence_residual_rel": np.logspace(-4, -3, 5),
        "phase_convergence_tolerance": 1.0e-3,
        "allmaras_phase_convergence_iterations": np.array([2, 3, 3, 4, 3]),
        "allmaras_phase_direct_node_count": np.ones(5) * 3,
        "allmaras_phase_continued_node_count": np.ones(5),
        "allmaras_phase_zero_amplitude_node_count": np.zeros(5),
        "allmaras_phase_drive_rms_snapshot": np.array([0.1, 0.2]),
        "allmaras_phase_drive_max_snapshot": np.array([0.2, 0.3]),
        "usadel_vs_gl_relative_l2_snapshot": np.array([0.5, 0.4]),
        "allmaras_phase_convergence_converged": np.ones(5, dtype=bool),
        "continuity_passes": True,
        "stationarity_passes": False,
    }

    saved = make_phasecg_ss_figures(dataset=dataset, output_dir=tmp_path, dpi=40)
    assert set(saved) == {"snapshot_fields", "physical_evolution", "numerical_diagnostics"}
    for path in saved.values():
        assert path.exists()
        assert path.stat().st_size > 0
