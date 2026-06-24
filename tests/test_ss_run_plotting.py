from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from pysnspd.plotting.ss_run import (
    plot_ss_available_snapshots,
    plot_ss_phi_snapshots,
    plot_ss_relaxation_history,
)


def _toy_mesh():
    nodes = np.array(
        [
            [0.0, -0.5],
            [1.0, -0.5],
            [2.0, -0.5],
            [0.0, 0.0],
            [1.0, 0.0],
            [2.0, 0.0],
            [0.0, 0.5],
            [1.0, 0.5],
            [2.0, 0.5],
        ],
        dtype=float,
    )
    triangles = np.array(
        [
            [0, 1, 4],
            [0, 4, 3],
            [1, 2, 5],
            [1, 5, 4],
            [3, 4, 7],
            [3, 7, 6],
            [4, 5, 8],
            [4, 8, 7],
        ],
        dtype=np.int64,
    )
    return SimpleNamespace(
        nodes=nodes,
        triangles=triangles,
        length_m=2.0,
        width_m=1.0,
    )


def _history_with_phi_snapshots():
    n_snap = 3
    n_nodes = 9

    return {
        "t_s": np.array([0.0, 1.0e-15, 2.0e-15]),
        "eta_R": np.array([1.0e-2, 1.0e-3, 1.0e-4]),
        "current_residual": np.array([1.0e-1, 1.0e-2, 1.0e-3]),
        "pairbreaking_max": np.array([0.2, 0.15, 0.1]),
        "terminal_voltage_V": np.array([0.0, 1.0e-6, 2.0e-6]),
        "delta_min_over_delta0": np.array([0.95, 0.97, 0.99]),
        "normal_current_fraction_rms": np.array([0.2, 0.1, 0.05]),
        "phi_snapshot_t_s": np.array([0.0, 1.0e-15, 2.0e-15]),
        "phi_snapshot_V": np.linspace(-1.0e-6, 1.0e-6, n_snap * n_nodes).reshape(n_snap, n_nodes),
    }


def test_phi_snapshot_plot_is_written(tmp_path):
    mesh = _toy_mesh()
    history = _history_with_phi_snapshots()

    out = plot_ss_phi_snapshots(
        mesh,
        history,
        tmp_path / "phi_snapshots.png",
        dpi=80,
        ncols=2,
    )

    assert out.exists()
    assert out.stat().st_size > 0


def test_available_snapshots_skips_missing_fields(tmp_path):
    mesh = _toy_mesh()
    history = _history_with_phi_snapshots()

    out = plot_ss_available_snapshots(
        mesh,
        history,
        tmp_path,
        dpi=80,
        ncols=2,
    )

    assert set(out) == {"phi"}
    assert out["phi"].exists()
    assert out["phi"].stat().st_size > 0


def test_relaxation_history_plot_is_written(tmp_path):
    history = _history_with_phi_snapshots()

    out = plot_ss_relaxation_history(
        history,
        tmp_path / "history.png",
        dpi=80,
    )

    assert out.exists()
    assert out.stat().st_size > 0