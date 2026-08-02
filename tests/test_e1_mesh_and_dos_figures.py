from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from pysnspd.mesh.delaunay import MeshData
from pysnspd.plotting.mesh import (
    _triangle_quality_metrics,
    plot_mesh_edge_length_histograms,
    plot_mesh_pytdgl_style,
    plot_mesh_triangle_quality,
)
from pysnspd.plotting.usadel_dos_curves import _eta_legend_line


def _strip_mesh() -> MeshData:
    x = np.array([0.0, 120.0, 220.0, 360.0]) * 1.0e-9
    y = np.array([-60.0, 0.0, 60.0]) * 1.0e-9
    xx, yy = np.meshgrid(x, y, indexing="xy")
    nodes = np.column_stack((xx.ravel(), yy.ravel()))
    triangles: list[tuple[int, int, int]] = []
    nx = x.size
    for row in range(y.size - 1):
        for col in range(x.size - 1):
            ll = row * nx + col
            lr = ll + 1
            ul = ll + nx
            ur = ul + 1
            triangles.extend(((ll, lr, ur), (ll, ur, ul)))
    return MeshData(
        nodes=nodes,
        triangles=np.asarray(triangles, dtype=np.int64),
        length_m=360.0e-9,
        width_m=120.0e-9,
        target_spacing_m=4.0e-9,
        seed=1,
        triangulation_method="test",
        boundary_guard_layers=0,
    )


def test_eta_legend_reports_catalog_dynes_broadening() -> None:
    catalog = SimpleNamespace(eta_J=1.315e-3 * 1.602176634e-22, metadata={})

    assert _eta_legend_line(catalog) == r"$\eta=1.315\,\mu\mathrm{eV}$"


def test_e1_mesh_figures_smoke(tmp_path) -> None:
    mesh = _strip_mesh()
    outputs = [
        plot_mesh_pytdgl_style(mesh, tmp_path / "mesh_pytdgl_style.pdf", dpi=60),
        plot_mesh_edge_length_histograms(
            mesh,
            tmp_path / "mesh_edge_length_histograms.pdf",
            dpi=60,
        ),
        plot_mesh_triangle_quality(mesh, tmp_path / "mesh_triangle_quality.pdf", dpi=60),
    ]

    for output in outputs:
        assert output.exists()
        assert output.stat().st_size > 0


def test_triangle_quality_metrics_are_bounded() -> None:
    mesh = _strip_mesh()
    minimum_angle_deg, shape_factor = _triangle_quality_metrics(
        mesh.nodes,
        mesh.triangles,
    )

    assert np.all((minimum_angle_deg > 0.0) & (minimum_angle_deg <= 60.0))
    assert np.all((shape_factor > 0.0) & (shape_factor <= 1.0))
