"""Smoke tests for PRE diagnostic plotting helpers."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from pysnspd.mesh.delaunay import MeshData
from pysnspd.mesh.edges import build_edge_data
from pysnspd.plotting.pre_diagnostics import write_pre_diagnostic_plots


def test_write_pre_diagnostic_plots_smoke(tmp_path):
    nodes = np.array(
        [
            [0.0, -0.5],
            [1.0, -0.5],
            [1.0, 0.5],
            [0.0, 0.5],
            [0.5, 0.0],
        ],
        dtype=float,
    )
    triangles = np.array(
        [
            [0, 1, 4],
            [1, 2, 4],
            [2, 3, 4],
            [3, 0, 4],
        ],
        dtype=np.int64,
    )

    mesh = MeshData(
        nodes=nodes,
        triangles=triangles,
        length_m=1.0,
        width_m=1.0,
        target_spacing_m=0.5,
        seed=12345,
        triangulation_method="test",
        boundary_guard_layers=0,
    )
    edge_data = build_edge_data(nodes, triangles, length_m=1.0, width_m=1.0)
    usadel_catalog = SimpleNamespace(
        calibration_q_values_m_inv=np.linspace(0.0, 1.0e8, 8),
        calibration_current_values_A=np.array(
            [0.0, 5.0e-6, 1.0e-5, 2.0e-5, 3.0e-5, 3.8e-5, 3.6e-5, 3.2e-5],
            dtype=float,
        ),
        metadata={"Ic_target_A": 3.88e-5},
    )

    outputs = write_pre_diagnostic_plots(
        mesh=mesh,
        edge_data=edge_data,
        usadel_catalog=usadel_catalog,
        output_dir=tmp_path,
        dpi=80,
    )

    assert set(outputs) == {
        "mesh_boundary_tags_png",
        "mesh_triangle_area_hist_png",
        "mesh_edge_length_hist_png",
        "usadel_supercurrent_curve_png",
    }
    for value in outputs.values():
        path = tmp_path / value.split("/")[-1]
        assert path.exists()
        assert path.stat().st_size > 0
