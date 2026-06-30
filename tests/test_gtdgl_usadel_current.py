"""Usadel supercurrent catalogue diagnostics with flat imports."""
from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from pysnspd.gtdgl.usadel_current import compute_usadel_supercurrent_diagnostic


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


def test_usadel_current_table_delta_q_interpolation(small_strip_mesh_bundle, gtdgl_material, stationary_seed_factory):
    mesh, _, ops = small_strip_mesh_bundle
    q_axis = np.array([0.0, 1.0e7, 2.0e7, 3.0e7])
    delta_axis = np.array([0.5, 0.9, 1.0]) * gtdgl_material.delta0_J
    # Canonical PRE layout: js_A_m2[delta, q].  Use an odd-in-q law through q>=0 table.
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
    assert diag.available
    assert diag.backend == "table:delta,q"
    assert diag.edge_js_usadel_A_m2.shape == (ops.n_edges,)
    assert np.all(np.isfinite(diag.edge_js_usadel_A_m2))
    assert diag.node_div_js_usadel_A_m3.shape == (ops.n_nodes,)


def test_usadel_current_callable_backend(small_strip_mesh_bundle, gtdgl_material, stationary_seed_factory):
    mesh, _, ops = small_strip_mesh_bundle

    def current_callable(q_m_inv, delta_J, Te_K):
        del Te_K
        return 1.0e-4 * q_m_inv * delta_J / gtdgl_material.delta0_J

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
    assert diag.available
    assert diag.backend.startswith("callable:")
    assert np.all(np.isfinite(diag.edge_js_usadel_A_m2))
