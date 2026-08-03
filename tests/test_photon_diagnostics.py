"""Smoke tests for single-run E3 photon diagnostics."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from pysnspd.plotting.photon_diagnostics import (
    make_photon_run_diagnostic_figures,
    nearest_unique_snapshot_indices,
    plot_photon_censored_recovery_diagnostics,
)


@dataclass
class DummyMesh:
    nodes: np.ndarray
    triangles: np.ndarray


def test_make_photon_run_diagnostic_figures(tmp_path: Path):
    mesh = DummyMesh(
        nodes=1.0e-9
        * np.array(
            [
                [-50.0, -20.0],
                [50.0, -20.0],
                [-50.0, 20.0],
                [50.0, 20.0],
            ]
        ),
        triangles=np.array([[0, 1, 2], [1, 3, 2]], dtype=np.int64),
    )
    snapshot_t_ps = np.array([0.0, 50.0, 60.0])
    phase = np.array(
        [
            [0.0, 0.1, 0.0, 0.1],
            [0.0, 0.3, 0.0, 0.3],
            [0.0, 0.2, 0.0, 0.2],
        ]
    )
    amplitude_J = 1.2 * 1.602176634e-22 * np.ones_like(phase)
    snapshots = {
        "snapshot_t_ps": snapshot_t_ps,
        "delta_snapshot_meV": amplitude_J / 1.602176634e-22,
        "psi_real_snapshot_J": amplitude_J * np.cos(phase),
        "psi_imag_snapshot_J": amplitude_J * np.sin(phase),
        "phi_snapshot_V": 1.0e-4 * np.array(
            [[-1.0, 1.0, -1.0, 1.0]] * snapshot_t_ps.size
        ),
        "Te_snapshot_K": np.array(
            [[0.9, 0.9, 0.9, 0.9], [1.5, 2.0, 1.5, 2.0], [1.1, 1.2, 1.1, 1.2]]
        ),
        "Tph_snapshot_K": np.array(
            [[0.9, 0.9, 0.9, 0.9], [1.0, 1.2, 1.0, 1.2], [0.95, 1.0, 0.95, 1.0]]
        ),
    }
    time = np.linspace(0.0, 80.0, 81)
    pulse = np.exp(-0.5 * ((time - 58.0) / 8.0) ** 2)
    history = {
        "t_ps": time,
        "photon_applied": time >= 50.0,
        "I_s_A": 30.0e-6 - 2.0e-6 * pulse,
        "I_b_A": np.full(time.size, 30.0e-6),
        "I_rf_A": 2.0e-6 * pulse,
        "V_tdgl_center_V": 0.2e-3 * pulse,
        "V_out_V": 0.1e-3 * pulse,
        "v_c_V": 0.05e-3 * pulse,
        "max_Te_K": 0.9 + 2.0 * pulse,
        "max_Tph_K": 0.9 + 0.6 * pulse,
        "mean_delta_over_delta0": 1.0 - 0.2 * pulse,
        "min_delta_over_delta0": 1.0 - 0.5 * pulse,
        "max_pairbreaking_ratio": 0.1 + 0.3 * pulse,
    }
    selected = np.array([1, 2], dtype=np.int64)
    diagnostic_shape = (selected.size, mesh.nodes.shape[0])
    snapshot_diagnostics = {
        "snapshot_t_ps": snapshot_t_ps[selected],
        "selected_indices": selected,
        "nodes_x_nm": 1.0e9 * mesh.nodes[:, 0],
        "nodes_y_nm": 1.0e9 * mesh.nodes[:, 1],
        "triangles": mesh.triangles,
        "q_abs_snapshot_m_inv": np.ones(diagnostic_shape) * 1.0e7,
        "joule_snapshot_W_m3": np.ones(diagnostic_shape) * 2.0,
        "P_total_snapshot_W_m3": np.ones(diagnostic_shape) * 3.0,
        "P_diff_snapshot_W_m3": np.zeros(diagnostic_shape),
        "P_esc_snapshot_W_m3": np.ones(diagnostic_shape) * 0.1,
        "u_e_snapshot_J_m3": np.ones(diagnostic_shape) * 4.0,
        "u_ph_snapshot_J_m3": np.ones(diagnostic_shape) * 5.0,
        "C_e_snapshot_J_m3_K": np.ones(diagnostic_shape) * 6.0,
        "C_ph_snapshot_J_m3_K": np.ones(diagnostic_shape) * 7.0,
        "snapshot_global_limits": {
            "q_abs_snapshot_m_inv": np.array([0.0, 2.0e7]),
            "joule_snapshot_W_m3": np.array([0.0, 20.0]),
            "P_total_snapshot_W_m3": np.array([-3.0, 30.0]),
            "P_diff_snapshot_W_m3": np.array([-2.0, 2.0]),
            "P_esc_snapshot_W_m3": np.array([0.0, 1.0]),
            "u_e_snapshot_J_m3": np.array([-4.0, 40.0]),
            "u_ph_snapshot_J_m3": np.array([0.0, 50.0]),
            "C_e_snapshot_J_m3_K": np.array([0.0, 60.0]),
            "C_ph_snapshot_J_m3_K": np.array([0.0, 70.0]),
        },
    }
    timing = {
        "latency": {"crossing_time_ps": 53.0, "t_lat_ps": 3.0},
        "recovery": {
            "selected": {
                "mode": "electrical",
                "entry_time_ps": 75.0,
                "t_rec_ps": 25.0,
            }
        },
    }

    saved = make_photon_run_diagnostic_figures(
        mesh=mesh,
        history=history,
        snapshots=snapshots,
        summary={"photon": {"time_s": 50.0e-12, "x_m": 0.0, "y_m": 0.0}},
        delta0_meV=1.2,
        xi_m=5.0e-9,
        requested_times_ps=[49.0, 51.0, 60.0],
        output_dir=tmp_path,
        dpi=40,
        timing=timing,
        snapshot_diagnostics=snapshot_diagnostics,
    )

    assert set(saved) == {
        "scalar_evolution",
        "field_evolution",
        "power_density_snapshots",
        "energy_heat_capacity_snapshots",
    }
    for path in saved.values():
        assert path.exists()
        assert path.stat().st_size > 0

    censored_timing = {
        "baseline": {
            "values": {
                "I_b_A": 30.0e-6,
                "I_s_A": 30.0e-6,
                "I_rf_A": 0.0,
                "V_out_V": 0.0,
                "v_c_V": 0.0,
                "V_tdgl_center_V": 0.0,
            }
        },
        "latency": {"detected": True, "t_lat_ps": 3.0},
        "recovery_criteria": {"hold_s": 10.0e-12},
        "recovery": {
            "selected": {
                "mode": "electrical",
                "recovered": False,
                "t_rec_ps": None,
                "lower_bound_ps": 30.0,
                "absolute_tolerances": {
                    "I_b_A": 0.3e-6,
                    "I_s_A": 0.3e-6,
                    "I_rf_A": 0.05e-6,
                    "V_out_V": 10.0e-6,
                    "v_c_V": 10.0e-6,
                    "V_tdgl_center_V": 10.0e-6,
                },
            }
        },
    }
    recovery_path = plot_photon_censored_recovery_diagnostics(
        history=history,
        summary={
            "photon": {"time_s": 50.0e-12},
            "circuit": {
                "params": {
                    "R_load_ohm": 50.0,
                    "R_bias_ohm": 1.0e4,
                    "L_bias_H": 1.0e-6,
                    "L_k_H": 10.0e-9,
                    "C_couple_F": 100.0e-12,
                }
            },
        },
        timing=censored_timing,
        output_path=tmp_path / "censored.pdf",
        dpi=40,
    )
    assert recovery_path.exists()
    assert recovery_path.stat().st_size > 0


def test_nearest_unique_photon_snapshot_indices():
    indices = nearest_unique_snapshot_indices(
        [0.0, 50.0, 60.0],
        [49.0, 51.0, 60.0],
    )

    np.testing.assert_array_equal(indices, np.array([1, 2]))
