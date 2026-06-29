"""Smoke tests for the pyTDGL-like OE7 comparison backend."""
from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from pysnspd.mesh.delaunay import MeshData
from pysnspd.mesh.edges import build_edge_data
from pysnspd.gtdgl.material import GTDGLMaterial, K_B_J_K
from pysnspd.gtdgl.operators import build_fv_operators
from pysnspd.gtdgl.pytdgl_like import SolverOptions, SparseSolver, TDGLSolver
from pysnspd.gtdgl.pytdgl_like.device import build_pytdgl_like_device
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
        tau_scale=1.0,
    )


def _seed(mesh, mat):
    psi = np.full(mesh.n_nodes, 0.9 * mat.delta0_J, dtype=np.complex128)
    return SimpleNamespace(
        node_psi_real_J=np.real(psi),
        node_psi_imag_J=np.imag(psi),
        node_phi_electric_V=np.zeros(mesh.n_nodes),
        node_Te_K=np.full(mesh.n_nodes, 0.9),
        node_Tph_K=np.full(mesh.n_nodes, 0.9),
    )


def test_pytdgl_like_public_names_and_options():
    opts = SolverOptions(solve_time=1.0e-3, dt_init=1.0e-4, dt_max=1.0e-4)
    opts.validate()
    assert opts.sparse_solver is SparseSolver.SUPERLU
    assert hasattr(TDGLSolver, "solve_for_psi_squared")
    assert hasattr(TDGLSolver, "adaptive_euler_step")
    assert hasattr(TDGLSolver, "solve_for_observables")
    assert hasattr(TDGLSolver, "update")
    assert hasattr(TDGLSolver, "solve")


def test_pytdgl_like_device_and_solver_zero_current_smoke():
    mesh, edge_data = _small_mesh()
    mat = _material()
    ops = build_fv_operators(mesh, edge_data)
    seed = _seed(mesh, mat)
    device = build_pytdgl_like_device(
        mesh=mesh,
        edge_data=edge_data,
        material=mat,
        ops=ops,
        Te_K=seed.node_Te_K,
        target_current_A=0.0,
    )
    options = SolverOptions(
        solve_time=2.0e-4,
        dt_init=1.0e-4,
        dt_max=1.0e-4,
        adaptive=False,
        terminal_psi=None,
    )
    solver = TDGLSolver(device, options, terminal_currents={"left": 0.0, "right": 0.0})
    sol = solver.solve()
    assert sol is not None
    assert np.all(np.isfinite(sol.tdgl_data.psi))
    assert np.all(np.isfinite(sol.tdgl_data.mu))


def test_solve_stationary_pytdgl_like_returns_relaxation_result():
    mesh, edge_data = _small_mesh()
    mat = _material()
    ops = build_fv_operators(mesh, edge_data)
    seed = _seed(mesh, mat)
    result = solve_stationary_pytdgl_like(
        mesh=mesh,
        edge_data=edge_data,
        seed=seed,
        material=mat,
        ops=ops,
        steps=2,
        dt_s=1.0e-18,
        target_current_A=0.0,
        terminal_psi=None,
        adaptive=False,
        n_snapshots=2,
    )
    assert result.summary["backend"] == "pytdgl_like_minimal_no_screening"
    assert result.summary["accepted_steps"] >= 1
    assert np.all(np.isfinite(result.state.psi_J))
    assert np.all(np.isfinite(result.state.phi_V))
    assert "snapshot_t_s" in result.history
    assert "phi_snapshot_t_s" in result.history
    assert "supercurrent_density_snapshot_A_m2" in result.history
    assert "normal_current_density_snapshot_A_m2" in result.history
    assert "edge_i" in result.history
    assert "edge_j" in result.history
    assert "pytdgl_like_poisson_residual_rel" in result.history
    assert "pytdgl_like_poisson_residual_snapshot" in result.history
    assert "pytdgl_like_native_supercurrent_snapshot" in result.history
    assert result.history["psi_snapshot_real_J"].shape == (2, mesh.n_nodes)
    assert result.history["edge_js_us_snapshot_A_m2"].shape == (2, ops.n_edges)
    assert result.history["pytdgl_like_poisson_residual_snapshot"].shape == (2, mesh.n_nodes)
    assert result.history["pytdgl_like_native_supercurrent_snapshot"].shape == (2, ops.n_edges)
    assert np.isfinite(result.summary["native_poisson_residual_rel_final"])


def test_pytdgl_like_adapter_keeps_terminal_currents_in_amperes():
    mesh, edge_data = _small_mesh()
    mat = _material()
    ops = build_fv_operators(mesh, edge_data)
    seed = _seed(mesh, mat)
    target_current_A = 3.5e-6

    result = solve_stationary_pytdgl_like(
        mesh=mesh,
        edge_data=edge_data,
        seed=seed,
        material=mat,
        ops=ops,
        steps=2,
        dt_s=1.0e-18,
        target_current_A=target_current_A,
        terminal_psi=None,
        adaptive=False,
        n_snapshots=2,
    )

    bc = result.summary["boundary_currents_A"]
    assert np.isclose(bc["left_A"], -target_current_A)
    assert np.isclose(bc["right_A"], target_current_A)
    assert np.isclose(bc["net_A"], 0.0)
    assert result.summary["terminal_neumann_current_unit_A"] > 0
