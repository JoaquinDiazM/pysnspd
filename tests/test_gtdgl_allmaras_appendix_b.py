"""Appendix-B Allmaras diagnostic tests with flat imports."""

from __future__ import annotations

import numpy as np

from pysnspd.solver.stationary import solve_stationary_pytdgl_like
from pysnspd.gtdgl.allmaras import (
    allmaras_coefficients,
    compute_allmaras_appendix_b_diagnostic,
    nonterminal_node_mask,
)

ALLMARAS_UPDATE_BACKEND = "appendix_b_normalized_phase_drive_harmonic_continuation_v2"


def test_allmaras_coefficients_are_appendix_b_shapes_and_positive(
    small_strip_mesh_bundle,
    gtdgl_material,
    stationary_seed_factory,
):
    mesh, _, _ = small_strip_mesh_bundle
    seed = stationary_seed_factory(mesh, gtdgl_material)

    psi = (seed.node_psi_real_J + 1j * seed.node_psi_imag_J) / gtdgl_material.delta0_J
    Te = np.full(mesh.n_nodes, 0.9)

    coeff = allmaras_coefficients(
        psi_dimensionless=psi,
        material=gtdgl_material,
        Te_K=Te,
    )

    assert coeff.gamma_kwt_dimensionless.shape == (mesh.n_nodes,)
    assert coeff.rho_kwt.shape == (mesh.n_nodes,)
    assert np.all(coeff.gamma_kwt_dimensionless > 0.0)
    assert np.all(coeff.rho_kwt >= 1.0)
    assert np.all(coeff.delta_mod_over_delta0 >= 0.0)
    assert np.all(coeff.solver_epsilon >= 0.0)
    assert np.all(coeff.xi_mod2_m2 > 0.0)


def test_allmaras_mismatch_diagnostic_has_bulk_mask_and_finite_drive(
    small_strip_mesh_bundle,
    gtdgl_material,
    stationary_seed_factory,
):
    mesh, _, ops = small_strip_mesh_bundle
    seed = stationary_seed_factory(mesh, gtdgl_material, q0_m_inv=2.0e7)

    psi = (seed.node_psi_real_J + 1j * seed.node_psi_imag_J) / gtdgl_material.delta0_J
    Te = np.full(mesh.n_nodes, 0.9)

    terminal = np.zeros(mesh.n_nodes, dtype=bool)
    x = mesh.nodes[:, 0]
    terminal[
        np.isclose(x, x.min(), rtol=0.0, atol=1.0e-15)
        | np.isclose(x, x.max(), rtol=0.0, atol=1.0e-15)
    ] = True

    diag = compute_allmaras_appendix_b_diagnostic(
        psi_dimensionless=psi,
        material=gtdgl_material,
        Te_K=Te,
        ops=ops,
        terminal_node_mask=terminal,
    )

    assert diag.edge_js_us_allmaras_A_m2.shape == (ops.n_edges,)
    assert diag.node_mismatch_divergence_A_m3.shape == (ops.n_nodes,)
    assert diag.node_phase_drive_abs_over_delta0.shape == (ops.n_nodes,)
    assert diag.bulk_node_mask.shape == (ops.n_nodes,)
    assert np.any(diag.bulk_node_mask)
    assert np.all(np.isfinite(diag.node_phase_drive_abs_over_delta0))


def test_solver_history_contains_appendix_b_diagnostics(
    small_strip_mesh_bundle,
    gtdgl_material,
    stationary_seed_factory,
):
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
    assert result.summary["allmaras_update_backend"] == ALLMARAS_UPDATE_BACKEND
    assert result.summary["allmaras_solver_u"] == 1.0
    assert "allmaras_mismatch_divergence_snapshot_A_m3" in result.history
    assert "allmaras_phase_drive_abs_over_delta0_snapshot" in result.history
    assert "allmaras_bulk_node_mask" in result.history
    assert result.history["allmaras_mismatch_divergence_snapshot_A_m3"].shape[1] == mesh.n_nodes


def test_nonterminal_mask_excludes_only_exact_terminals(small_strip_mesh_bundle):
    mesh, _, _ = small_strip_mesh_bundle

    terminal = np.zeros(mesh.n_nodes, dtype=bool)
    x = mesh.nodes[:, 0]
    terminal[
        np.isclose(x, x.min(), rtol=0.0, atol=1.0e-15)
        | np.isclose(x, x.max(), rtol=0.0, atol=1.0e-15)
    ] = True

    bulk = nonterminal_node_mask(terminal, n_nodes=mesh.n_nodes)

    assert bulk.shape == (mesh.n_nodes,)
    assert bulk.dtype == np.bool_
    assert np.any(bulk)

    assert np.array_equal(bulk, ~terminal)


def test_nonterminal_mask_all_terminal_degenerate_fixture_is_nonempty(small_strip_mesh_bundle):
    mesh, _, _ = small_strip_mesh_bundle

    terminal = np.ones(mesh.n_nodes, dtype=bool)

    bulk = nonterminal_node_mask(terminal, n_nodes=mesh.n_nodes)

    assert bulk.shape == (mesh.n_nodes,)
    assert bulk.dtype == np.bool_
    assert np.any(bulk)
    assert np.all(bulk)
