"""Tests for the Appendix-B Allmaras local forcing injected into the solver."""
from __future__ import annotations

import numpy as np

from pysnspd.gtdgl.adapter import solve_stationary_pytdgl_like
from pysnspd.gtdgl.allmaras import compute_allmaras_forcing_dimensionless


def test_allmaras_forcing_dimensionless_is_finite(small_strip_mesh_bundle, gtdgl_material, stationary_seed_factory):
    mesh, _, ops = small_strip_mesh_bundle
    seed = stationary_seed_factory(mesh, gtdgl_material, q0_m_inv=1.0e7)
    psi = (seed.node_psi_real_J + 1j * seed.node_psi_imag_J) / gtdgl_material.delta0_J
    Te = np.full(mesh.n_nodes, 0.9)
    lap = np.zeros(mesh.n_nodes, dtype=np.complex128)
    edge_js = np.zeros(ops.n_edges, dtype=float)

    forcing = compute_allmaras_forcing_dimensionless(
        psi_dimensionless=psi,
        psi_laplacian_dimensionless=lap,
        material=gtdgl_material,
        Te_K=Te,
        ops=ops,
        length_scale_m=5.0e-9,
        edge_js_usadel_A_m2=edge_js,
    )

    assert forcing.forcing_dimensionless.shape == (mesh.n_nodes,)
    assert forcing.edge_js_us_A_m2.shape == (ops.n_edges,)
    assert forcing.edge_js_gl_A_m2.shape == (ops.n_edges,)
    assert np.all(np.isfinite(forcing.forcing_dimensionless))
    assert np.all(np.isfinite(forcing.node_mismatch_divergence_A_m3))


def test_solver_uses_appendix_b_wz_update_backend(small_strip_mesh_bundle, gtdgl_material, stationary_seed_factory):
    mesh, edge_data, ops = small_strip_mesh_bundle
    result = solve_stationary_pytdgl_like(
        mesh=mesh,
        edge_data=edge_data,
        seed=stationary_seed_factory(mesh, gtdgl_material, q0_m_inv=1.0e7),
        material=gtdgl_material,
        ops=ops,
        steps=3,
        dt_s=1.0e-18,
        target_current_A=0.0,
        terminal_psi=0.0,
        adaptive=False,
        n_snapshots=3,
    )

    assert result.summary["allmaras_coefficients_backend"] == "appendix_b_allmaras_wz_update_v1"
    assert result.summary["allmaras_update_backend"] == "appendix_b_explicit_forcing_rho_kwt_wz_v1"
    assert "allmaras_update_forcing_max_abs" in result.history
    assert np.all(np.isfinite(result.history["allmaras_update_forcing_max_abs"]))
