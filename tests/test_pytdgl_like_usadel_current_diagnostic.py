from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from pysnspd.mesh.delaunay import MeshData
from pysnspd.mesh.edges import build_edge_data
from pysnspd.gtdgl.material import GTDGLMaterial, K_B_J_K
from pysnspd.gtdgl.operators import build_fv_operators
from pysnspd.gtdgl.pytdgl_like.usadel_current import compute_usadel_supercurrent_diagnostic


def _small_mesh():
    length = 2.0e-8
    width = 1.0e-8
    nodes = np.array(
        [
            [0.0, -0.5 * width],
            [0.5 * length, -0.5 * width],
            [length, -0.5 * width],
            [0.0, 0.0],
            [0.5 * length, 0.0],
            [length, 0.0],
            [0.0, 0.5 * width],
            [0.5 * length, 0.5 * width],
            [length, 0.5 * width],
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
    mesh = MeshData(
        nodes=nodes,
        triangles=triangles,
        length_m=length,
        width_m=width,
        target_spacing_m=5.0e-9,
        seed=1,
    )
    return mesh, build_edge_data(nodes, triangles, length_m=length, width_m=width)


def _material():
    Tc = 8.65
    return GTDGLMaterial(
        Tc_K=Tc,
        D_m2_s=1.58e-4,
        sigma_n_S_m=4.2e5,
        delta0_J=1.764 * K_B_J_K * Tc,
        thickness_m=7.0e-9,
        width_m=1.2e-7,
        tau_ee_Tc_s=5.0e-12,
        tau_ep_Tc_s=24.7e-12,
        tau_scale=1.0,
    )


def test_usadel_supercurrent_diagnostic_interpolates_q_table():
    mesh, edge_data = _small_mesh()
    ops = build_fv_operators(mesh, edge_data)
    mat = _material()
    q0 = 2.0e7
    psi = 0.9 * np.exp(1j * q0 * mesh.nodes[:, 0])
    catalog = SimpleNamespace(
        q_axis_m_inv=np.array([0.0, 1.0e7, 3.0e7]),
        js_A_m2=np.array([0.0, 1.0e10, 3.0e10]),
    )

    diag = compute_usadel_supercurrent_diagnostic(
        usadel_catalog=catalog,
        psi_dimensionless=psi,
        material=mat,
        Te_K=np.full(mesh.n_nodes, 0.9),
        ops=ops,
    )

    assert diag.available
    assert diag.edge_js_usadel_A_m2.shape == (ops.n_edges,)
    assert diag.node_js_usadel_x_A_m2.shape == (mesh.n_nodes,)
    assert np.all(np.isfinite(diag.edge_js_usadel_A_m2))
    assert np.nanmax(np.abs(diag.edge_js_usadel_A_m2)) > 0.0


def test_usadel_supercurrent_diagnostic_marks_missing_catalog_unavailable():
    mesh, edge_data = _small_mesh()
    ops = build_fv_operators(mesh, edge_data)
    mat = _material()
    diag = compute_usadel_supercurrent_diagnostic(
        usadel_catalog=None,
        psi_dimensionless=np.ones(mesh.n_nodes, dtype=np.complex128),
        material=mat,
        Te_K=np.full(mesh.n_nodes, 0.9),
        ops=ops,
    )
    assert not diag.available
    assert np.isnan(diag.edge_js_usadel_A_m2).all()
