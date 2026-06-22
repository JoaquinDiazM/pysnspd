"""OE7 smoke tests for stationary gTDGL/Poisson relaxation."""
from __future__ import annotations

from types import SimpleNamespace
import numpy as np

from pysnspd.mesh.delaunay import MeshData
from pysnspd.mesh.edges import build_edge_data
from pysnspd.gtdgl.material import GTDGLMaterial, K_B_J_K
from pysnspd.gtdgl.operators import build_fv_operators
from pysnspd.gtdgl.relax import (
    compute_current_fields,
    current_residual,
    kwt_local_update,
    relax_stationary_gtdgl,
    solve_poisson_potential,
)


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


def test_tau_scale_changes_effective_relaxation_times():
    mat = _material()
    assert np.isclose(mat.tau_ee_s(8.65), 0.5e-12)
    assert np.isclose(mat.tau_ep_s(8.65), 2.47e-12)
    assert mat.tau_sc_s(8.65) < mat.tau_ep_s(8.65)


def test_kwt_local_update_preserves_state_for_zero_forcing():
    mesh, edge_data = _small_mesh()
    ops = build_fv_operators(mesh, edge_data)
    mat = _material()

    psi = np.full(mesh.n_nodes, 0.7 * mat.delta0_J, dtype=np.complex128)
    phi = np.zeros(mesh.n_nodes)
    Te = np.full(mesh.n_nodes, 0.9)
    forcing = np.zeros(mesh.n_nodes, dtype=np.complex128)

    updated, ok, _ = kwt_local_update(
        psi_J=psi,
        phi_V=phi,
        Te_K=Te,
        forcing_J=forcing,
        dt_s=1.0e-16,
        material=mat,
    )

    assert ok
    assert np.allclose(updated, psi, rtol=1.0e-12, atol=1.0e-30)
    assert ops.n_edges > 0


def test_poisson_zero_supercurrent_returns_zero_potential():
    mesh, edge_data = _small_mesh()
    ops = build_fv_operators(mesh, edge_data)
    mat = _material()

    phi = solve_poisson_potential(
        edge_js_us_A_m2=np.zeros(ops.n_edges),
        material=mat,
        ops=ops,
    )

    assert np.all(np.isfinite(phi))
    assert abs(float(np.mean(phi))) < 1.0e-15
    assert np.max(np.abs(phi)) < 1.0e-15


def test_stationary_relaxation_smoke_zero_current_seed():
    mesh, edge_data = _small_mesh()
    ops = build_fv_operators(mesh, edge_data)
    mat = _material()

    psi = np.full(mesh.n_nodes, 0.7 * mat.delta0_J, dtype=np.complex128)
    seed = SimpleNamespace(
        node_psi_real_J=np.real(psi),
        node_psi_imag_J=np.imag(psi),
        node_phi_electric_V=np.zeros(mesh.n_nodes),
        node_Te_K=np.full(mesh.n_nodes, 0.9),
        node_Tph_K=np.full(mesh.n_nodes, 0.9),
    )

    result = relax_stationary_gtdgl(
        mesh=mesh,
        edge_data=edge_data,
        seed=seed,
        material=mat,
        ops=ops,
        steps=5,
        min_steps=1,
        dt_s=1.0e-16,
        tolerance_eta=1.0,
        tolerance_current_residual=1.0,
        target_current_A=0.0,
    )

    assert result.summary["accepted_steps"] >= 1
    assert np.isfinite(result.summary["terminal_voltage_V"])
    assert np.isfinite(result.summary["current_residual"])

    currents = compute_current_fields(
        psi_J=result.state.psi_J,
        phi_V=result.state.phi_V,
        Te_K=result.state.Te_K,
        material=mat,
        ops=ops,
    )
    assert np.isfinite(current_residual(currents, mesh))