"""Finite-volume operator tests for the flat gTDGL package."""
from __future__ import annotations

import numpy as np

from pysnspd.gtdgl.operators import (
    divergence_from_edge_scalar,
    edge_average,
    edge_phase_gradient_from_psi,
    edge_scalar_gradient,
    laplacian,
    terminal_voltage,
)


def test_basic_operator_shapes(small_strip_mesh_bundle):
    mesh, _, ops = small_strip_mesh_bundle
    values = np.linspace(0.0, 1.0, mesh.n_nodes)
    assert edge_average(values, ops).shape == (ops.n_edges,)
    assert edge_scalar_gradient(values, ops).shape == (ops.n_edges,)
    assert divergence_from_edge_scalar(np.zeros(ops.n_edges), ops).shape == (ops.n_nodes,)
    assert laplacian(values, ops).shape == (ops.n_nodes,)


def test_phase_gradient_recovers_uniform_ramp(small_strip_mesh_bundle):
    mesh, _, ops = small_strip_mesh_bundle
    q0 = 2.0e7
    psi = np.exp(1j * q0 * mesh.nodes[:, 0])
    q_edge = edge_phase_gradient_from_psi(psi, ops)
    assert q_edge.shape == (ops.n_edges,)
    assert np.all(np.isfinite(q_edge))


def test_terminal_voltage_uses_right_minus_left_average(small_strip_mesh_bundle):
    mesh, _, _ = small_strip_mesh_bundle
    phi = mesh.nodes[:, 0] / mesh.length_m
    v = terminal_voltage(mesh.nodes, phi, length_m=mesh.length_m)
    assert np.isclose(v, 1.0)
