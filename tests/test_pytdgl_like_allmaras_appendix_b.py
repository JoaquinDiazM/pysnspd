from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from pysnspd.mesh.delaunay import MeshData
from pysnspd.mesh.edges import build_edge_data
from pysnspd.gtdgl.material import GTDGLMaterial, K_B_J_K
from pysnspd.gtdgl.operators import build_fv_operators
from pysnspd.gtdgl.pytdgl_like.allmaras import (
    allmaras_coefficients,
    compute_allmaras_appendix_b_diagnostic,
)
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


def test_allmaras_coefficients_are_appendix_b_shapes_and_positive():
    mesh, edge_data = _small_mesh()
    mat = _material()
    psi = _seed(mesh, mat).node_psi_real_J / mat.delta0_J + 0j
    Te = np.full(mesh.n_nodes, 0.9)

    coeff = allmaras_coefficients(psi_dimensionless=psi, material=mat, Te_K=Te)

    assert coeff.gamma_kwt_dimensionless.shape == (mesh.n_nodes,)
    assert coeff.rho_kwt.shape == (mesh.n_nodes,)
    assert np.all(coeff.gamma_kwt_dimensionless > 0.0)
    assert np.all(coeff.rho_kwt >= 1.0)
    assert np.all(coeff.delta_mod_over_delta0 >= 0.0)
    assert np.all(coeff.solver_epsilon >= 0.0)
    assert np.all(coeff.xi_mod2_m2 > 0.0)


def test_allmaras_mismatch_diagnostic_has_bulk_mask_and_finite_drive():
    mesh, edge_data = _small_mesh()
    mat = _material()
    ops = build_fv_operators(mesh, edge_data)
    seed = _seed(mesh, mat)
    psi = (seed.node_psi_real_J + 1j * seed.node_psi_imag_J) / mat.delta0_J
    Te = np.full(mesh.n_nodes, 0.9)
    terminal = np.zeros(mesh.n_nodes, dtype=bool)
    x = mesh.nodes[:, 0]
    terminal[(x == x.min()) | (x == x.max())] = True

    diag = compute_allmaras_appendix_b_diagnostic(
        psi_dimensionless=psi,
        material=mat,
        Te_K=Te,
        ops=ops,
        terminal_node_mask=terminal,
        bulk_guard_layers=1,
    )

    assert diag.edge_js_us_allmaras_A_m2.shape == (ops.n_edges,)
    assert diag.node_mismatch_divergence_A_m3.shape == (ops.n_nodes,)
    assert diag.node_phase_drive_abs_over_delta0.shape == (ops.n_nodes,)
    assert diag.bulk_node_mask.shape == (ops.n_nodes,)
    assert np.any(diag.bulk_node_mask)
    assert np.all(np.isfinite(diag.node_phase_drive_abs_over_delta0))


def test_solver_history_contains_appendix_b_diagnostics_without_wz_rewrite():
    mesh, edge_data = _small_mesh()
    mat = _material()
    ops = build_fv_operators(mesh, edge_data)
    result = solve_stationary_pytdgl_like(
        mesh=mesh,
        edge_data=edge_data,
        seed=_seed(mesh, mat),
        material=mat,
        ops=ops,
        steps=3,
        dt_s=1.0e-18,
        target_current_A=0.0,
        terminal_psi=0.0,
        adaptive=False,
        n_snapshots=3,
    )

    assert result.summary["allmaras_coefficients_backend"] == "appendix_b_coefficients_v1_without_wz_rewrite"
    assert result.summary["allmaras_solver_u"] == 1.0
    assert "allmaras_mismatch_divergence_snapshot_A_m3" in result.history
    assert "allmaras_phase_drive_abs_over_delta0_snapshot" in result.history
    assert "allmaras_bulk_node_mask" in result.history
    assert result.history["allmaras_mismatch_divergence_snapshot_A_m3"].shape[1] == mesh.n_nodes



def test_bulk_mask_falls_back_to_terminal_only_on_tiny_mesh():
    from pysnspd.gtdgl.pytdgl_like.allmaras import bulk_node_mask_from_terminal_mask

    mesh, edge_data = _small_mesh()
    ops = build_fv_operators(mesh, edge_data)
    terminal = np.zeros(mesh.n_nodes, dtype=bool)
    x = mesh.nodes[:, 0]
    terminal[(x == x.min()) | (x == x.max())] = True

    bulk = bulk_node_mask_from_terminal_mask(ops, terminal, guard_layers=1)

    assert np.any(bulk)
    assert np.array_equal(bulk, ~terminal)
