"""Tests for the Appendix-B Allmaras local forcing injected into the solver."""

from __future__ import annotations

import numpy as np
import pytest

from pysnspd.gtdgl.adapter import solve_stationary_pytdgl_like
from pysnspd.gtdgl.allmaras import (
    PhaseDriveContinuationSolver,
    compute_allmaras_forcing_dimensionless,
)

ALLMARAS_UPDATE_BACKEND = "appendix_b_normalized_phase_drive_harmonic_continuation_v2"


def test_normalized_phase_drive_matches_exact_direct_formula(small_strip_mesh_bundle, gtdgl_material):
    mesh, _, ops = small_strip_mesh_bundle
    psi = np.exp(1j * np.linspace(0.0, 0.2, mesh.n_nodes))
    mismatch = np.linspace(1.0, float(mesh.n_nodes), mesh.n_nodes)
    delta0 = float(gtdgl_material.delta0_J)
    continuation = PhaseDriveContinuationSolver.from_operators(ops)

    drive, info = continuation.solve(
        psi_dimensionless=psi,
        mismatch_divergence_A_m3=mismatch,
        correction_C_J2_m3_A=np.full(mesh.n_nodes, delta0 * delta0),
        delta0_J=delta0,
    )

    expected = 1j * mismatch * psi / (np.abs(psi) ** 2)
    assert np.allclose(drive, expected, rtol=1.0e-13, atol=1.0e-13)
    assert info.converged
    assert info.continued_node_count == 0


def test_low_amplitude_phase_drive_is_continued_to_convergence(small_strip_mesh_bundle, gtdgl_material):
    mesh, _, ops = small_strip_mesh_bundle
    psi = np.ones(mesh.n_nodes, dtype=np.complex128)
    low_node = mesh.n_nodes // 2
    psi[low_node] = 1.0e-4
    delta0 = float(gtdgl_material.delta0_J)
    continuation = PhaseDriveContinuationSolver.from_operators(
        ops,
        direct_amplitude_fraction=1.0e-2,
        tolerance=1.0e-8,
        max_iterations=24,
    )

    drive, info = continuation.solve(
        psi_dimensionless=psi,
        mismatch_divergence_A_m3=np.ones(mesh.n_nodes),
        correction_C_J2_m3_A=np.full(mesh.n_nodes, delta0 * delta0),
        delta0_J=delta0,
    )

    assert info.converged
    assert info.continued_node_count == 1
    assert np.isclose(drive[low_node], 1j, rtol=1.0e-6, atol=1.0e-6)


def test_zero_amplitude_node_is_inside_continuation_domain(small_strip_mesh_bundle, gtdgl_material):
    mesh, _, ops = small_strip_mesh_bundle
    psi = np.ones(mesh.n_nodes, dtype=np.complex128)
    zero_node = mesh.n_nodes // 2
    psi[zero_node] = 0.0
    delta0 = float(gtdgl_material.delta0_J)
    continuation = PhaseDriveContinuationSolver.from_operators(
        ops,
        tolerance=1.0e-8,
        max_iterations=24,
    )

    drive, info = continuation.solve(
        psi_dimensionless=psi,
        mismatch_divergence_A_m3=np.ones(mesh.n_nodes),
        correction_C_J2_m3_A=np.full(mesh.n_nodes, delta0 * delta0),
        delta0_J=delta0,
    )

    assert info.converged
    assert info.continued_node_count == 1
    assert info.zero_amplitude_node_count == 1
    assert np.isclose(drive[zero_node], 1j, rtol=1.0e-6, atol=1.0e-6)


def test_phase_drive_refuses_unconverged_harmonic_extension(small_strip_mesh_bundle, gtdgl_material):
    mesh, _, ops = small_strip_mesh_bundle
    psi = np.ones(mesh.n_nodes, dtype=np.complex128)
    middle = np.isclose(
        mesh.nodes[:, 0],
        np.median(np.unique(mesh.nodes[:, 0])),
        rtol=0.0,
        atol=1.0e-15,
    )
    psi[middle] = 0.0
    delta0 = float(gtdgl_material.delta0_J)
    continuation = PhaseDriveContinuationSolver.from_operators(
        ops,
        tolerance=1.0e-14,
        max_iterations=1,
    )

    with pytest.raises(RuntimeError, match="did not converge"):
        continuation.solve(
            psi_dimensionless=psi,
            mismatch_divergence_A_m3=np.linspace(1.0, 3.0, mesh.n_nodes),
            correction_C_J2_m3_A=np.full(mesh.n_nodes, delta0 * delta0),
            delta0_J=delta0,
        )


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
    assert result.summary["allmaras_update_backend"] == ALLMARAS_UPDATE_BACKEND
    assert "allmaras_update_forcing_max_abs" in result.history
    assert np.all(np.isfinite(result.history["allmaras_update_forcing_max_abs"]))
