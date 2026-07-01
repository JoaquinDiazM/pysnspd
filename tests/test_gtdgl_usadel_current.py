"""Strict 3D Usadel supercurrent catalogue diagnostics with flat imports."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from pysnspd.gtdgl.usadel_current import compute_usadel_supercurrent_diagnostic


def _strict_3d_catalog(material, *, q_scale: float = 2.0e3) -> SimpleNamespace:
    q_axis = np.array([0.0, 1.0e7, 2.0e7, 3.0e7], dtype=float)
    delta_axis = np.array([0.0, 0.5, 0.9, 1.0], dtype=float) * material.delta0_J
    Te_axis = np.array([0.5, 0.9, 1.2], dtype=float)

    js = np.empty((Te_axis.size, delta_axis.size, q_axis.size), dtype=float)
    for iT, T in enumerate(Te_axis):
        # A deliberately simple monotone mock law.  The exact values are not
        # physical; the test is about strict [Te, |Delta|, q] interpolation.
        temp_factor = 1.0 - 0.05 * (T - 0.9)
        for iD, delta in enumerate(delta_axis):
            js[iT, iD, :] = temp_factor * (delta / material.delta0_J) * q_scale * q_axis

    return SimpleNamespace(
        js_A_m2=js,
        Te_axis_K=Te_axis,
        delta_axis_J=delta_axis,
        q_axis_m_inv=q_axis,
    )


def test_usadel_current_unavailable_without_catalog(small_strip_mesh_bundle, gtdgl_material, stationary_seed_factory):
    mesh, _, ops = small_strip_mesh_bundle
    seed = stationary_seed_factory(mesh, gtdgl_material, q0_m_inv=1.0e7)
    psi = (seed.node_psi_real_J + 1j * seed.node_psi_imag_J) / gtdgl_material.delta0_J
    diag = compute_usadel_supercurrent_diagnostic(
        usadel_catalog=None,
        psi_dimensionless=psi,
        material=gtdgl_material,
        Te_K=seed.node_Te_K,
        ops=ops,
    )
    assert not diag.available
    assert diag.edge_q_m_inv.shape == (ops.n_edges,)
    assert np.all(np.isnan(diag.edge_js_usadel_A_m2))


def test_usadel_current_strict_3d_table_interpolation(small_strip_mesh_bundle, gtdgl_material, stationary_seed_factory):
    mesh, _, ops = small_strip_mesh_bundle
    catalog = _strict_3d_catalog(gtdgl_material)
    seed = stationary_seed_factory(mesh, gtdgl_material, q0_m_inv=1.5e7, amplitude=0.9)
    psi = (seed.node_psi_real_J + 1j * seed.node_psi_imag_J) / gtdgl_material.delta0_J

    diag = compute_usadel_supercurrent_diagnostic(
        usadel_catalog=catalog,
        psi_dimensionless=psi,
        material=gtdgl_material,
        Te_K=seed.node_Te_K,
        ops=ops,
    )

    assert diag.available, diag.reason
    assert "Te" in diag.backend and "delta" in diag.backend and "q" in diag.backend
    assert diag.edge_js_usadel_A_m2.shape == (ops.n_edges,)
    assert np.all(np.isfinite(diag.edge_js_usadel_A_m2))
    assert diag.node_div_js_usadel_A_m3.shape == (ops.n_nodes,)


def test_usadel_current_legacy_delta_q_table_is_rejected(small_strip_mesh_bundle, gtdgl_material, stationary_seed_factory):
    mesh, _, ops = small_strip_mesh_bundle
    q_axis = np.array([0.0, 1.0e7, 2.0e7, 3.0e7])
    delta_axis = np.array([0.5, 0.9, 1.0]) * gtdgl_material.delta0_J
    js = np.outer(delta_axis / gtdgl_material.delta0_J, 2.0e3 * q_axis)
    catalog = SimpleNamespace(js_A_m2=js, q_axis_m_inv=q_axis, delta_axis_J=delta_axis)
    seed = stationary_seed_factory(mesh, gtdgl_material, q0_m_inv=1.5e7, amplitude=0.9)
    psi = (seed.node_psi_real_J + 1j * seed.node_psi_imag_J) / gtdgl_material.delta0_J

    diag = compute_usadel_supercurrent_diagnostic(
        usadel_catalog=catalog,
        psi_dimensionless=psi,
        material=gtdgl_material,
        Te_K=seed.node_Te_K,
        ops=ops,
    )

    assert not diag.available
    assert "3D" in diag.reason or "Te_axis" in diag.reason


def test_usadel_current_callable_backend_is_rejected_by_strict_policy(
    small_strip_mesh_bundle,
    gtdgl_material,
    stationary_seed_factory,
):
    mesh, _, ops = small_strip_mesh_bundle

    def current_callable(q_m_inv, delta_J, Te_K):
        return 1.0e-4 * q_m_inv * delta_J / gtdgl_material.delta0_J * np.ones_like(Te_K)

    catalog = SimpleNamespace(interpolate_supercurrent_density_A_m2=current_callable)
    seed = stationary_seed_factory(mesh, gtdgl_material, q0_m_inv=1.0e7)
    psi = (seed.node_psi_real_J + 1j * seed.node_psi_imag_J) / gtdgl_material.delta0_J
    diag = compute_usadel_supercurrent_diagnostic(
        usadel_catalog=catalog,
        psi_dimensionless=psi,
        material=gtdgl_material,
        Te_K=seed.node_Te_K,
        ops=ops,
    )

    assert not diag.available
    assert "3D" in diag.reason or "table" in diag.reason
