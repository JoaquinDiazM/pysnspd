"""OE7 boundary-policy smoke tests."""
from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from pysnspd.gtdgl.material import GTDGLMaterial, K_B_J_K
from pysnspd.gtdgl.operators import build_fv_operators
from pysnspd.gtdgl.relax import relax_stationary_gtdgl
from pysnspd.mesh.delaunay import MeshData
from pysnspd.mesh.edges import build_edge_data


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
    edge_data = build_edge_data(nodes, triangles, length_m=length, width_m=width)
    return mesh, edge_data


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
        tau_scale=0.10,
    )


def _seed(mesh, mat):
    psi = np.full(mesh.n_nodes, 0.7 * mat.delta0_J, dtype=np.complex128)
    return SimpleNamespace(
        node_psi_real_J=np.real(psi),
        node_psi_imag_J=np.imag(psi),
        node_phi_electric_V=np.zeros(mesh.n_nodes),
        node_Te_K=np.full(mesh.n_nodes, 0.9),
        node_Tph_K=np.full(mesh.n_nodes, 0.9),
    )


def test_boundary_policy_switches_smoke_zero_current():
    mesh, edge_data = _small_mesh()
    ops = build_fv_operators(mesh, edge_data)
    mat = _material()
    seed = _seed(mesh, mat)

    for delta_policy, poisson_policy in [
        ("current_inversion", "target_flux"),
        ("current_inversion", "zero_flux"),
        ("vacuum_only", "target_flux"),
        ("normal_terminal", "target_flux"),
    ]:
        result = relax_stationary_gtdgl(
            mesh=mesh,
            edge_data=edge_data,
            seed=seed,
            material=mat,
            ops=ops,
            steps=2,
            min_steps=1,
            dt_s=1.0e-16,
            tolerance_eta=1.0,
            tolerance_current_residual=1.0,
            target_current_A=0.0,
            delta_boundary_policy=delta_policy,
            poisson_terminal_policy=poisson_policy,
        )
        assert result.summary["accepted_steps"] >= 1
        assert result.summary["delta_boundary_policy"] == delta_policy
        assert result.summary["poisson_terminal_policy"] == poisson_policy
        assert np.isfinite(result.summary["terminal_voltage_V"])
