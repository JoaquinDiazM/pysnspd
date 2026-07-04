from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from pysnspd.plotting.ss_power_figures import (
    plot_ss_snapshot_power_balance_maps,
    plot_ss_snapshot_power_energy_maps,
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
    t_s = np.linspace(0.0, 60.0e-12, n_snapshots)
    x = np.linspace(-1.0, 1.0, n_nodes)
    p_ep = np.vstack([1.0e12 * (1.0 + 0.5 * np.sin(np.pi * x + 0.2 * k)) for k in range(n_snapshots)])
    joule = np.vstack([1.0e15 * (1.2 + x * x + 0.1 * k) for k in range(n_snapshots)])
    kappa = np.vstack([1.0e-6 + 1.0e-3 * (1.0 + x) ** 2 + 2.0e-4 * k for k in range(n_snapshots)])
    p_esc = 0.8 * p_ep
    return {
        "snapshot_t_s": t_s,
        "P_total_snapshot_W_m3": p_ep,
        "joule_snapshot_W_m3": joule,
        "kappa_s_snapshot_W_m_K": kappa,
        "P_esc_snapshot_W_m3": p_esc,
    }


def _dataset(mesh: DummyMesh) -> dict[str, np.ndarray | str]:
    return {
        "run_name": "unit_power",
        "x_nm": mesh.nodes[:, 0] * 1.0e9,
        "y_nm": mesh.nodes[:, 1] * 1.0e9,
        "triangles": mesh.triangles,
    }


def test_plot_ss_snapshot_power_energy_maps(tmp_path: Path):
    mesh = _dummy_mesh()
    power = _power_payload(n_snapshots=9, n_nodes=mesh.nodes.shape[0])
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
    power = _power_payload(n_snapshots=9, n_nodes=mesh.nodes.shape[0])
    output = plot_ss_snapshot_power_balance_maps(
        mesh,
        power,
        _dataset(mesh),
        tmp_path / "ss_snapshot_power_balance_maps.png",
        dpi=80,
    )
    assert output.exists()
    assert output.stat().st_size > 0


def _snapshot_payload(n_snapshots: int, n_nodes: int) -> dict[str, np.ndarray]:
    t_s = np.linspace(0.0, 60.0e-12, n_snapshots)
    x = np.linspace(0.0, 1.0, n_nodes)
    delta = np.vstack([1.0 + 0.03 * k + 0.05 * x for k in range(n_snapshots)])
    phi = np.vstack([2.0e-3 * (k + 1) * (x - 0.5) for k in range(n_snapshots)])
    phi[0] = 0.2 * np.sin(20.0 * np.pi * x)
    jx = np.vstack([2.0e10 * (1.0 + 0.02 * k + x) for k in range(n_snapshots)])
    jy = np.zeros_like(jx)
    return {
        "snapshot_t_s": t_s,
        "delta_snapshot_meV": delta,
        "phi_snapshot_V": phi,
        "jtot_snapshot_x_A_m2": jx,
        "jtot_snapshot_y_A_m2": jy,
        "delta0_meV": np.array([1.0]),
        "javg_A_m2": np.array([2.0e10]),
    }


def test_plot_ss_snapshot_state_atlas(tmp_path: Path):
    mesh = _dummy_mesh()
    snapshots = _snapshot_payload(n_snapshots=9, n_nodes=mesh.nodes.shape[0])
    output = plot_ss_snapshot_state_atlas(
        mesh,
        snapshots,
        _dataset(mesh),
        tmp_path / "ss_snapshot_state_atlas.png",
        dpi=80,
    )
    assert output.exists()
    assert output.stat().st_size > 0
