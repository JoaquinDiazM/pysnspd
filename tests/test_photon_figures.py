"""Smoke tests for the standard pipeline-03 plotting module."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from pysnspd.plotting.photon_figures import (
    plot_photon_center_scalar_snapshot_rows,
)


@dataclass
class DummyMesh:
    nodes: np.ndarray
    triangles: np.ndarray


def _strip_mesh() -> DummyMesh:
    x = np.array([-50.0, -15.0, 15.0, 50.0]) * 1.0e-9
    lower = np.column_stack([x, np.full(x.size, -50.0e-9)])
    upper = np.column_stack([x, np.full(x.size, 50.0e-9)])
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
    return DummyMesh(nodes=np.vstack([lower, upper]), triangles=triangles)


def test_plot_photon_center_scalar_snapshot_rows(tmp_path: Path):
    mesh = _strip_mesh()
    n_nodes = mesh.nodes.shape[0]
    x = np.linspace(-1.0, 1.0, n_nodes)
    snapshots = {
        "snapshot_t_s": np.array([0.0, 1.0e-12]),
        "delta_snapshot_meV": np.vstack([1.5 + 0.05 * x, 1.2 + 0.10 * x]),
        "delta0_J": np.array([1.5e-3 * 1.602176634e-19]),
        "phi_snapshot_V": np.vstack([1.0e-6 * x, 2.0e-6 * x]),
        "Te_snapshot_K": np.vstack([1.0 + 0.1 * x, 1.5 + 0.2 * x]),
        "Tph_snapshot_K": np.vstack([0.9 + 0.02 * x, 1.0 + 0.04 * x]),
    }

    output = plot_photon_center_scalar_snapshot_rows(
        mesh=mesh,
        snapshots=snapshots,
        summary={},
        requested_times_ps=[0.0, 1.0],
        output_path=tmp_path / "photon_center_scalar_snapshots.png",
        center_width_nm=100.0,
        dpi=60,
    )

    assert output.exists()
    assert output.stat().st_size > 0
