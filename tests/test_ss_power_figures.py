from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from pysnspd.gtdgl.snapshot_diagnostics import compute_snapshot_joule_power_density
from pysnspd.plotting.ss_power_figures import (
    plot_ss_snapshot_power_balance_maps,
    plot_ss_snapshot_power_energy_maps,
    plot_ss_snapshot_profile_comparison,
    plot_ss_snapshot_runtime_metrics,
    plot_ss_snapshot_state_atlas,
)


@dataclass
class DummyMesh:
    nodes: np.ndarray
    triangles: np.ndarray


def _dummy_mesh() -> DummyMesh:
    x = np.array([0.0, 100e-9, 200e-9, 300e-9])
    y0 = np.array([-50e-9, -50e-9, -50e-9, -50e-9])
    y1 = np.array([50e-9, 50e-9, 50e-9, 50e-9])
    nodes = np.column_stack([np.r_[x, x], np.r_[y0, y1]])
    triangles = np.array(
        [
            [0, 1, 4],
            [1, 5, 4],
            [1, 2, 5],
            [2, 6, 5],
            [2, 3, 6],
            [3, 7, 6],
        ],
        dtype=np.int64,
    )
    return DummyMesh(nodes=nodes, triangles=triangles)


def _power_payload(n_snapshots: int, n_nodes: int) -> dict[str, np.ndarray]:
    t_s = np.linspace(0.0, 8.0e-12, n_snapshots)
    x = np.linspace(-1.0, 1.0, n_nodes)
    p_ep = np.vstack([1.0e12 * (0.2 + 0.5 * np.sin(np.pi * x + 0.2 * k)) for k in range(n_snapshots)])
    joule = np.vstack([1.0e15 * (1.2 + x * x + 0.1 * k) for k in range(n_snapshots)])
    kappa = np.vstack([1.0e-6 + 1.0e-3 * (1.0 + x) ** 2 + 2.0e-4 * k for k in range(n_snapshots)])
    p_esc = 0.8 * p_ep
    q_abs = np.vstack([3.0e7 * (1.0 + 0.2 * (x + 1.0) + 0.03 * k) for k in range(n_snapshots)])
    u_e = np.vstack([20.0 + 0.4 * k + 2.0 * x for k in range(n_snapshots)])
    u_ph = np.vstack([7.0 + 0.2 * k + 1.0 * x for k in range(n_snapshots)])
    c_e = np.vstack([5.0e1 + 2.0 * k + x for k in range(n_snapshots)])
    c_ph = np.vstack([2.0e2 + 3.0 * k + 2.0 * x for k in range(n_snapshots)])
    idx = np.tile(np.arange(n_nodes, dtype=float), (n_snapshots, 1))
    return {
        "snapshot_t_s": t_s,
        "P_total_snapshot_W_m3": p_ep,
        "joule_snapshot_W_m3": joule,
        "kappa_s_snapshot_W_m_K": kappa,
        "P_esc_snapshot_W_m3": p_esc,
        "q_abs_snapshot_m_inv": q_abs,
        "u_e_snapshot_J_m3": u_e,
        "u_ph_snapshot_J_m3": u_ph,
        "C_e_snapshot_J_m3_K": c_e,
        "C_ph_snapshot_J_m3_K": c_ph,
        "power_table_iTe": 0.0 * idx,
        "power_table_iTph": 1.0 + 0.0 * idx,
        "power_table_iDelta": 2.0 + 0.0 * idx,
        "power_table_iQ": 3.0 + 0.0 * idx,
    }


def _dataset(mesh: DummyMesh) -> dict[str, np.ndarray | str]:
    n_hist = 21
    t_ps = np.linspace(0.0, 8.0, n_hist)
    return {
        "run_name": "unit_power",
        "x_nm": mesh.nodes[:, 0] * 1.0e9,
        "y_nm": mesh.nodes[:, 1] * 1.0e9,
        "triangles": mesh.triangles,
        "sigma_n_S_m": 4.2e5,
        "t_ps": t_ps,
        "thermal_mean_Te_K_history": np.linspace(0.9, 1.4, n_hist),
        "thermal_max_Te_K_history": np.linspace(0.9, 2.0, n_hist),
        "thermal_mean_Tph_K_history": np.linspace(0.9, 1.1, n_hist),
        "thermal_max_Tph_K_history": np.linspace(0.9, 1.3, n_hist),
        "thermal_max_rate_K_per_ps_history": np.linspace(0.10, 0.01, n_hist),
        "thermal_substeps_history": np.linspace(1.0, 3.0, n_hist),
        "thermal_max_abs_dTe_K_history": np.linspace(0.05, 0.01, n_hist),
        "thermal_max_abs_dTph_K_history": np.linspace(0.02, 0.005, n_hist),
        "thermal_max_P_J_W_m3_history": np.linspace(1.0e15, 2.0e15, n_hist),
        "thermal_max_P_ep_W_m3_history": np.linspace(1.0e12, 3.0e12, n_hist),
        "thermal_max_P_esc_W_m3_history": np.linspace(2.0e11, 8.0e11, n_hist),
        "thermal_max_P_diff_W_m3_history": np.linspace(-5.0e12, 5.0e12, n_hist),
    }


def test_snapshot_joule_power_density_is_positive_definite_jn_squared():
    history = {
        "edge_jn_snapshot_A_m2": np.array([[-2.0, 4.0], [3.0, -5.0]]),
        "edge_i": np.array([0, 1], dtype=np.int64),
        "edge_j": np.array([1, 2], dtype=np.int64),
        "edge_length_m": np.ones(2),
        "dual_face_length_m": np.ones(2),
    }

    joule = compute_snapshot_joule_power_density(
        history,
        sigma_n_S_m=2.0,
        n_snap=2,
        n_nodes=3,
    )

    assert joule is not None
    expected = np.array(
        [
            [2.0, 0.5 * (2.0 + 8.0), 8.0],
            [4.5, 0.5 * (4.5 + 12.5), 12.5],
        ]
    )
    np.testing.assert_allclose(joule, expected)
    assert np.all(joule >= 0.0)


def _snapshot_payload(n_snapshots: int, n_nodes: int) -> dict[str, np.ndarray]:
    t_s = np.linspace(0.0, 8.0e-12, n_snapshots)
    x = np.linspace(0.0, 1.0, n_nodes)
    delta = np.vstack([1.0 + 0.03 * k + 0.05 * x for k in range(n_snapshots)])
    phi = np.vstack([2.0e-3 * (k + 1) * (x - 0.5) for k in range(n_snapshots)])
    phi[0] = 0.2 * np.sin(20.0 * np.pi * x)
    jx = np.vstack([2.0e10 * (1.0 + 0.02 * k + x) for k in range(n_snapshots)])
    jy = np.zeros_like(jx)
    te = np.vstack([0.9 + 0.2 * k + 0.1 * x for k in range(n_snapshots)])
    tph = np.vstack([0.9 + 0.05 * k + 0.02 * x for k in range(n_snapshots)])
    return {
        "snapshot_t_s": t_s,
        "delta_snapshot_meV": delta,
        "phi_snapshot_V": phi,
        "jtot_snapshot_x_A_m2": jx,
        "jtot_snapshot_y_A_m2": jy,
        "Te_snapshot_K": te,
        "Tph_snapshot_K": tph,
        "delta0_meV": np.array([1.0]),
        "javg_A_m2": np.array([2.0e10]),
    }


def test_plot_ss_snapshot_power_energy_maps(tmp_path: Path):
    mesh = _dummy_mesh()
    power = _power_payload(n_snapshots=5, n_nodes=mesh.nodes.shape[0])
    output = plot_ss_snapshot_power_energy_maps(
        mesh,
        power,
        _dataset(mesh),
        tmp_path / "ss_snapshot_power_energy_maps.png",
        dpi=80,
    )
    assert output.exists()
    assert output.stat().st_size > 0


def test_plot_ss_snapshot_power_balance_maps(tmp_path: Path):
    mesh = _dummy_mesh()
    power = _power_payload(n_snapshots=5, n_nodes=mesh.nodes.shape[0])
    output = plot_ss_snapshot_power_balance_maps(
        mesh,
        power,
        _dataset(mesh),
        tmp_path / "ss_snapshot_power_balance_maps.png",
        dpi=80,
    )
    assert output.exists()
    assert output.stat().st_size > 0


def test_plot_ss_snapshot_state_atlas(tmp_path: Path):
    mesh = _dummy_mesh()
    snapshots = _snapshot_payload(n_snapshots=5, n_nodes=mesh.nodes.shape[0])
    output = plot_ss_snapshot_state_atlas(
        mesh,
        snapshots,
        _dataset(mesh),
        tmp_path / "ss_snapshot_state_atlas.png",
        dpi=80,
    )
    assert output.exists()
    assert output.stat().st_size > 0


def test_plot_ss_snapshot_profile_comparison(tmp_path: Path):
    mesh = _dummy_mesh()
    snapshots = _snapshot_payload(n_snapshots=5, n_nodes=mesh.nodes.shape[0])
    power = _power_payload(n_snapshots=5, n_nodes=mesh.nodes.shape[0])
    output = plot_ss_snapshot_profile_comparison(
        mesh,
        snapshots,
        power,
        _dataset(mesh),
        tmp_path / "ss_snapshot_profile_comparison.png",
        dpi=80,
    )
    assert output.exists()
    assert output.stat().st_size > 0


def test_plot_ss_snapshot_runtime_metrics(tmp_path: Path):
    mesh = _dummy_mesh()
    snapshots = _snapshot_payload(n_snapshots=5, n_nodes=mesh.nodes.shape[0])
    power = _power_payload(n_snapshots=5, n_nodes=mesh.nodes.shape[0])
    output = plot_ss_snapshot_runtime_metrics(
        power,
        snapshots,
        _dataset(mesh),
        tmp_path / "ss_snapshot_runtime_metrics.png",
        dpi=80,
    )
    assert output.exists()
    assert output.stat().st_size > 0
