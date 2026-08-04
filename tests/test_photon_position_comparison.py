"""Smoke coverage for the center/edge photon comparison atlas."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from pysnspd.plotting.photon_comparison import make_photon_position_figures


@dataclass
class DummyMesh:
    nodes: np.ndarray
    triangles: np.ndarray


def _history(*, scale: float) -> dict[str, np.ndarray]:
    time = np.linspace(0.0, 60.0, 61)
    pulse = scale * np.exp(-0.5 * ((time - 53.0) / 3.0) ** 2)
    return {
        "t_ps": time,
        "photon_applied": time >= 50.0,
        "I_b_A": np.full(time.size, 30.0e-6),
        "I_s_A": 30.0e-6 - 2.0e-6 * pulse,
        "I_rf_A": 2.0e-6 * pulse,
        "V_out_V": 0.1e-3 * pulse,
        "v_c_V": 0.05e-3 * pulse,
        "V_tdgl_center_V": 0.2e-3 * pulse,
        "max_Te_K": 0.9 + 2.0 * pulse,
        "max_Tph_K": 0.9 + 0.4 * pulse,
        "mean_delta_over_delta0": 1.0 - 0.2 * pulse,
    }


def _snapshots(*, scale: float) -> dict[str, np.ndarray]:
    times = np.array([50.0, 51.0, 52.0, 53.0])
    x = np.linspace(0.0, 1.0, 4)
    phase = scale * times[:, None] * 0.01 * x[None, :]
    amplitude = 1.2e-22 * (1.0 - 0.1 * scale * x[None, :])
    amplitude = np.broadcast_to(amplitude, phase.shape)
    return {
        "snapshot_t_ps": times,
        "delta_snapshot_meV": amplitude / 1.602176634e-22,
        "psi_real_snapshot_J": amplitude * np.cos(phase),
        "psi_imag_snapshot_J": amplitude * np.sin(phase),
        "phi_snapshot_V": scale * 1.0e-4 * np.broadcast_to(x, phase.shape),
        "Te_snapshot_K": 0.9 + scale * np.broadcast_to(x, phase.shape),
        "Tph_snapshot_K": 0.9 + 0.2 * scale * np.broadcast_to(x, phase.shape),
    }


def _diagnostics(*, scale: float) -> dict[str, object]:
    shape = (4, 4)
    positive_keys = (
        "q_abs_snapshot_m_inv",
        "joule_snapshot_W_m3",
        "P_esc_snapshot_W_m3",
        "u_ph_snapshot_J_m3",
        "C_e_snapshot_J_m3_K",
        "C_ph_snapshot_J_m3_K",
    )
    signed_keys = (
        "P_total_snapshot_W_m3",
        "P_diff_snapshot_W_m3",
        "u_e_snapshot_J_m3",
    )
    result: dict[str, object] = {
        "snapshot_t_ps": np.array([50.0, 51.0, 52.0, 53.0]),
        "selected_indices": np.arange(4, dtype=np.int64),
    }
    limits: dict[str, np.ndarray] = {}
    for index, key in enumerate(positive_keys, start=1):
        values = scale * index * np.arange(1, 17, dtype=float).reshape(shape)
        result[key] = values
        limits[key] = np.array([0.0, float(np.max(values)) * 2.0])
    for index, key in enumerate(signed_keys, start=1):
        values = scale * index * np.linspace(-1.0, 1.0, 16).reshape(shape)
        result[key] = values
        limit = float(np.max(np.abs(values))) * 2.0
        limits[key] = np.array([-limit, limit])
    result["snapshot_global_limits"] = limits
    return result


def _timing(*, recovered: bool) -> dict[str, object]:
    baseline = {
        "I_b_A": 30.0e-6,
        "I_s_A": 30.0e-6,
        "I_rf_A": 0.0,
        "V_out_V": 0.0,
        "v_c_V": 0.0,
        "V_tdgl_center_V": 0.0,
    }
    tolerances = {
        "I_b_A": 0.3e-6,
        "I_s_A": 0.3e-6,
        "I_rf_A": 0.05e-6,
        "V_out_V": 10.0e-6,
        "v_c_V": 10.0e-6,
        "V_tdgl_center_V": 10.0e-6,
    }
    return {
        "baseline": {"values": baseline},
        "latency": {"detected": True, "t_lat_ps": 2.0},
        "recovery_criteria": {"hold_s": 10.0e-12},
        "recovery": {
            "selected": {
                "mode": "electrical",
                "recovered": recovered,
                "t_rec_ps": 8.0 if recovered else None,
                "absolute_tolerances": tolerances,
            }
        },
    }


def test_make_photon_position_figures_includes_all_comparison_atlases(
    tmp_path: Path,
) -> None:
    mesh = DummyMesh(
        nodes=1.0e-9
        * np.array(
            [[-50.0, -20.0], [50.0, -20.0], [-50.0, 20.0], [50.0, 20.0]]
        ),
        triangles=np.array([[0, 1, 2], [1, 3, 2]], dtype=np.int64),
    )
    circuit = {
        "params": {
            "R_load_ohm": 50.0,
            "R_bias_ohm": 1.0e4,
            "L_bias_H": 1.0e-6,
            "L_k_H": 10.0e-9,
            "C_couple_F": 100.0e-12,
        }
    }
    saved = make_photon_position_figures(
        mesh=mesh,
        center_history=_history(scale=1.0),
        center_snapshots=_snapshots(scale=1.0),
        center_summary={"photon": {"x_m": 0.0, "y_m": 0.0}, "circuit": circuit},
        edge_history=_history(scale=0.7),
        edge_snapshots=_snapshots(scale=0.7),
        edge_summary={"photon": {"x_m": 20.0e-9, "y_m": 10.0e-9}, "circuit": circuit},
        delta0_meV=1.2,
        xi_m=5.0e-9,
        requested_times_ps=[50.0, 51.0, 52.0, 53.0],
        output_dir=tmp_path,
        dpi=40,
        center_timing=_timing(recovered=False),
        edge_timing=_timing(recovered=True),
        center_snapshot_diagnostics=_diagnostics(scale=1.0),
        edge_snapshot_diagnostics=_diagnostics(scale=0.7),
    )

    assert set(saved) == {
        "field_comparison",
        "circuit_comparison",
        "power_density_comparison",
        "energy_heat_capacity_comparison",
        "censored_recovery_diagnostics",
    }
    for path in saved.values():
        assert path.exists()
        assert path.stat().st_size > 0
