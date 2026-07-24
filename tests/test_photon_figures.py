"""Smoke tests for the standard pipeline-03 plotting module."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from pysnspd.plotting.photon_figures import (
    plot_photon_circuit_response,
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


def test_plot_photon_circuit_response_with_timing_annotations(tmp_path: Path):
    t_ps = np.linspace(0.0, 80.0, 161)
    pulse = np.exp(-0.5 * ((t_ps - 35.0) / 8.0) ** 2)
    output = plot_photon_circuit_response(
        history={
            "t_ps": t_ps,
            "photon_applied": t_ps >= 20.0,
            "I_s_A": 30.0e-6 - 5.0e-6 * pulse,
            "I_rf_A": 5.0e-6 * pulse,
            "V_tdgl_center_V": 0.2e-3 * pulse,
            "V_out_V": 0.25e-3 * pulse,
        },
        summary={},
        timing={
            "latency": {"crossing_time_ps": 24.0, "t_lat_ps": 4.0},
            "recovery": {
                "selected": {
                    "mode": "electrical",
                    "entry_time_ps": 60.0,
                    "t_rec_ps": 40.0,
                }
            },
        },
        output_path=tmp_path / "photon-circuit-timing.png",
        dpi=60,
    )

    assert output.exists()
    assert output.stat().st_size > 0
