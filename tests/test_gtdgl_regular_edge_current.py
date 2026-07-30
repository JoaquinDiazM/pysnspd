"""Gauge-invariant regular-edge tests for the production Usadel closure."""
from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from pysnspd.gtdgl.usadel_current import compute_usadel_supercurrent_diagnostic
from pysnspd.mesh.operators import divergence_from_edge_scalar


def _constant_stiffness_catalog(material, stiffness: float = 1.0e48) -> SimpleNamespace:
    delta_axis = np.array([0.0, material.delta0_J])
    return SimpleNamespace(
        js_stiffness_A_per_m_J2=np.full((1, 2, 2), stiffness, dtype=float),
        Te_axis_K=np.array([0.9]),
        delta2_axis_J2=delta_axis**2,
        q_axis_m_inv=np.array([0.0, 1.0e9]),
    )


def _diagnostic(psi, *, ops, material, link=None):
    return compute_usadel_supercurrent_diagnostic(
        usadel_catalog=_constant_stiffness_catalog(material),
        psi_dimensionless=np.asarray(psi, dtype=np.complex128),
        material=material,
        Te_K=np.full(ops.n_nodes, 0.9),
        ops=ops,
        edge_link_variable=link,
    )


def test_regular_edge_current_is_exactly_zero_at_zero_endpoint(
    small_strip_mesh_bundle,
    gtdgl_material,
):
    _mesh, _edge_data, ops = small_strip_mesh_bundle
    psi = np.exp(1j * np.linspace(0.0, 0.3, ops.n_nodes))
    zero_node = ops.n_nodes // 2
    psi[zero_node] = 0.0

    diag = _diagnostic(psi, ops=ops, material=gtdgl_material)
    incident = (ops.edge_i == zero_node) | (ops.edge_j == zero_node)

    assert np.any(incident)
    assert np.array_equal(diag.edge_pair_flow_J2_m_inv[incident], np.zeros(np.count_nonzero(incident)))
    assert np.array_equal(diag.edge_js_usadel_A_m2[incident], np.zeros(np.count_nonzero(incident)))
    assert np.array_equal(diag.edge_js_gl_A_m2[incident], np.zeros(np.count_nonzero(incident)))


def test_regular_edge_current_has_quadratic_amplitude_power(
    small_strip_mesh_bundle,
    gtdgl_material,
):
    mesh, _edge_data, ops = small_strip_mesh_bundle
    phase = 2.0e7 * np.asarray(mesh.nodes)[:, 0]
    full = _diagnostic(0.8 * np.exp(1j * phase), ops=ops, material=gtdgl_material)
    small = _diagnostic(0.08 * np.exp(1j * phase), ops=ops, material=gtdgl_material)
    active = np.abs(full.edge_js_usadel_A_m2) > 1.0e-12 * np.max(
        np.abs(full.edge_js_usadel_A_m2)
    )

    assert np.any(active)
    assert np.allclose(
        small.edge_js_usadel_A_m2[active] / full.edge_js_usadel_A_m2[active],
        1.0e-2,
        rtol=2.0e-13,
        atol=0.0,
    )


def test_regular_edge_current_matches_discrete_plane_wave(
    small_strip_mesh_bundle,
    gtdgl_material,
):
    mesh, _edge_data, ops = small_strip_mesh_bundle
    qx = 1.5e7
    amplitude = 0.7
    stiffness = 1.0e48
    psi = amplitude * np.exp(1j * qx * np.asarray(mesh.nodes)[:, 0])
    catalog = _constant_stiffness_catalog(gtdgl_material, stiffness=stiffness)

    diag = compute_usadel_supercurrent_diagnostic(
        usadel_catalog=catalog,
        psi_dimensionless=psi,
        material=gtdgl_material,
        Te_K=np.full(ops.n_nodes, 0.9),
        ops=ops,
    )
    dtheta = qx * np.asarray(ops.edge_vec_m)[:, 0]
    expected_pair = (
        (amplitude * gtdgl_material.delta0_J) ** 2
        * np.sin(dtheta)
        / np.asarray(ops.edge_length_m)
    )

    assert np.allclose(diag.edge_pair_flow_J2_m_inv, expected_pair, rtol=3.0e-13, atol=1.0e-50)
    assert np.allclose(
        diag.edge_js_usadel_A_m2,
        stiffness * expected_pair,
        rtol=3.0e-13,
        atol=1.0e-13 * np.max(np.abs(stiffness * expected_pair)),
    )


def test_regular_edge_current_is_gauge_invariant(
    small_strip_mesh_bundle,
    gtdgl_material,
):
    _mesh, _edge_data, ops = small_strip_mesh_bundle
    amplitude = np.linspace(0.2, 1.0, ops.n_nodes)
    phase = np.linspace(-0.4, 0.5, ops.n_nodes)
    psi = amplitude * np.exp(1j * phase)
    gauge = np.linspace(0.7, -0.2, ops.n_nodes)
    transformed = psi * np.exp(1j * gauge)
    link = np.exp(1j * (gauge[ops.edge_i] - gauge[ops.edge_j]))

    reference = _diagnostic(psi, ops=ops, material=gtdgl_material)
    changed = _diagnostic(transformed, ops=ops, material=gtdgl_material, link=link)

    assert np.allclose(changed.edge_q_m_inv, reference.edge_q_m_inv, rtol=2.0e-13, atol=1.0e-8)
    assert np.allclose(changed.edge_pair_flow_J2_m_inv, reference.edge_pair_flow_J2_m_inv)
    assert np.allclose(changed.edge_js_usadel_A_m2, reference.edge_js_usadel_A_m2)


def test_mismatch_is_formed_on_edges_before_one_divergence(
    small_strip_mesh_bundle,
    gtdgl_material,
):
    mesh, _edge_data, ops = small_strip_mesh_bundle
    psi = np.exp(1j * 2.0e7 * np.asarray(mesh.nodes)[:, 0])
    diag = _diagnostic(psi, ops=ops, material=gtdgl_material)
    expected = divergence_from_edge_scalar(
        diag.edge_js_usadel_A_m2 - diag.edge_js_gl_A_m2,
        ops,
    )

    assert np.array_equal(diag.edge_js_mismatch_A_m2, diag.edge_js_usadel_A_m2 - diag.edge_js_gl_A_m2)
    assert np.allclose(diag.node_mismatch_divergence_A_m3, expected, rtol=0.0, atol=0.0)
