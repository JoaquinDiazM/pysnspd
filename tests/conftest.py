"""Shared fixtures for the flattened gTDGL backend tests."""
from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from pysnspd.mesh.delaunay import MeshData
from pysnspd.mesh.edges import build_edge_data
from pysnspd.gtdgl.material import GTDGLMaterial, K_B_J_K
from pysnspd.gtdgl.operators import build_fv_operators


@pytest.fixture()
def small_strip_mesh_bundle():
    """Return a tiny but non-degenerate 2D strip mesh, edge data and FV operators."""
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
    edge_data = build_edge_data(nodes, triangles, length_m=length, width_m=width)
    ops = build_fv_operators(mesh, edge_data)
    return mesh, edge_data, ops


@pytest.fixture()
def gtdgl_material():
    Tc = 8.65
    return GTDGLMaterial(
        Tc_K=Tc,
        D_m2_s=1.58e-4,
        sigma_n_S_m=4.2e5,
        delta0_J=1.764 * K_B_J_K * Tc,
        thickness_m=7.0e-9,
        width_m=1.2e-7,
        tau_ee_Tc_s=0.5e-12,
        tau_ep_Tc_s=2.47e-12,
    )


@pytest.fixture()
def stationary_seed_factory():
    def _make(mesh, material, *, q0_m_inv: float = 0.0, amplitude: float = 0.9):
        phase = q0_m_inv * np.asarray(mesh.nodes, dtype=float)[:, 0]
        psi = amplitude * material.delta0_J * np.exp(1j * phase)
        return SimpleNamespace(
            node_psi_real_J=np.real(psi),
            node_psi_imag_J=np.imag(psi),
            node_phi_electric_V=np.zeros(mesh.n_nodes),
            node_Te_K=np.full(mesh.n_nodes, 0.9),
            node_Tph_K=np.full(mesh.n_nodes, 0.9),
        )

    return _make
