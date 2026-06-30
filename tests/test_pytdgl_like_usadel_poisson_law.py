from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from pysnspd.mesh.delaunay import MeshData
from pysnspd.mesh.edges import build_edge_data
from pysnspd.gtdgl.material import GTDGLMaterial, K_B_J_K
from pysnspd.gtdgl.operators import build_fv_operators
from pysnspd.gtdgl.pytdgl_like.adapter import solve_stationary_pytdgl_like


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


def _seed(mesh, mat):
    q0 = 2.0e7
    psi = 0.9 * mat.delta0_J * np.exp(1j * q0 * mesh.nodes[:, 0])
    return SimpleNamespace(
        node_psi_real_J=np.real(psi),
        node_psi_imag_J=np.imag(psi),
        node_phi_electric_V=np.zeros(mesh.n_nodes),
        node_Te_K=np.full(mesh.n_nodes, 0.9),
        node_Tph_K=np.full(mesh.n_nodes, 0.9),
    )


def test_solve_stationary_pytdgl_like_can_use_usadel_poisson_supercurrent():
    mesh, edge_data = _small_mesh()
    mat = _material()
    ops = build_fv_operators(mesh, edge_data)
    seed = _seed(mesh, mat)
    catalog = SimpleNamespace(
        q_axis_m_inv=np.array([-1.0e8, 0.0, 1.0e8]),
        delta_axis_J=np.array([0.0, mat.delta0_J]),
        js_A_m2=np.array(
            [
                [-0.0, 0.0, 0.0],
                [-5.0e10, 0.0, 5.0e10],
            ],
            dtype=float,
        ),
    )

    result = solve_stationary_pytdgl_like(
        mesh=mesh,
        edge_data=edge_data,
        seed=seed,
        material=mat,
        ops=ops,
        steps=2,
        dt_s=1.0e-18,
        target_current_A=0.0,
        terminal_psi=None,
        adaptive=False,
        n_snapshots=2,
        usadel_catalog=catalog,
        supercurrent_law="usadel-poisson",
    )

    assert result.summary["supercurrent_law"] == "usadel_poisson"
    assert result.summary["usadel_current_available"]
    assert "edge_js_actual_snapshot_A_m2" in result.history
    assert "edge_js_gl_snapshot_A_m2" in result.history
    assert result.history["edge_js_actual_snapshot_A_m2"].shape == (2, ops.n_edges)
