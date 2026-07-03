from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from pysnspd.plotting.ss_power_figures import make_ss_snapshot_power_figures


@dataclass
class DummyMesh:
    nodes: np.ndarray
    triangles: np.ndarray
    length_m: float = 2.0e-7
    width_m: float = 1.0e-7
    target_spacing_m: float = 5.0e-9


def test_make_ss_snapshot_power_figures(tmp_path: Path):
    nodes = np.array(
        [
            [0.0, 0.0],
            [1.0e-7, 0.0],
            [2.0e-7, 0.0],
            [0.0, 1.0e-7],
            [1.0e-7, 1.0e-7],
            [2.0e-7, 1.0e-7],
        ],
        dtype=float,
    )
    triangles = np.array([[0, 1, 3], [1, 4, 3], [1, 2, 4], [2, 5, 4]], dtype=np.int64)
    mesh = DummyMesh(nodes=nodes, triangles=triangles)
    raw_ss = tmp_path / "raw" / "run" / "ss"
    raw_ss.mkdir(parents=True)
    n_snap, n_nodes = 3, nodes.shape[0]
    t = np.array([0.0, 0.5e-12, 1.0e-12])
    base = np.linspace(0.8, 1.0, n_nodes)

    np.savez(
        raw_ss / "stationary_snapshots.npz",
        snapshot_t_s=t,
        delta0_meV=np.array([1.3]),
        javg_A_m2=np.array([2.0e10]),
        delta_snapshot_meV=np.vstack([base, 0.95 * base, 0.9 * base]),
        phi_snapshot_V=np.vstack([base, 2 * base, 3 * base]) * 1.0e-3,
        jtot_snapshot_x_A_m2=np.ones((n_snap, n_nodes)) * 2.0e10,
        jtot_snapshot_y_A_m2=np.zeros((n_snap, n_nodes)),
    )
    np.savez(
        raw_ss / "snapshot_power_energy_diagnostics.npz",
        snapshot_t_s=t,
        P_total_snapshot_W_m3=np.ones((n_snap, n_nodes)) * 1.0e14,
        joule_snapshot_W_m3=np.ones((n_snap, n_nodes)) * 2.0e14,
        P_esc_snapshot_W_m3=np.ones((n_snap, n_nodes)) * 1.0e13,
        kappa_s_snapshot_W_m_K=np.ones((n_snap, n_nodes)) * 0.1,
        q_abs_snapshot_m_inv=np.ones((n_snap, n_nodes)) * 3.0e7,
        u_e_snapshot_J_m3=np.ones((n_snap, n_nodes)) * 1.0e3,
        u_ph_snapshot_J_m3=np.ones((n_snap, n_nodes)) * 1.0e2,
        C_e_snapshot_J_m3_K=np.ones((n_snap, n_nodes)) * 1.0e2,
        C_ph_snapshot_J_m3_K=np.ones((n_snap, n_nodes)) * 1.0e3,
        power_table_iTe=np.tile(np.arange(n_nodes), (n_snap, 1)),
        power_table_iTph=np.tile(np.arange(n_nodes), (n_snap, 1)),
        power_table_iDelta=np.tile(np.arange(n_nodes), (n_snap, 1)),
        power_table_iQ=np.tile(np.arange(n_nodes), (n_snap, 1)),
    )

    saved = make_ss_snapshot_power_figures(
        mesh=mesh,
        dataset={
            "run_name": "dummy",
            "x_nm": nodes[:, 0] * 1.0e9,
            "y_nm": nodes[:, 1] * 1.0e9,
            "triangles": triangles,
            "delta0_meV": 1.3,
            "javg_A_m2": 2.0e10,
        },
        raw_ss=raw_ss,
        output_dir=tmp_path / "figures",
        dpi=80,
    )

    assert "snapshot_state_atlas" in saved
    assert "snapshot_power_energy_maps" in saved
    assert "snapshot_runtime_metrics" in saved
    for path in saved.values():
        assert path.exists()
        assert path.stat().st_size > 0
