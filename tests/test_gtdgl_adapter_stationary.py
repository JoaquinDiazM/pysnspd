"""Stationary adapter smoke tests for the flattened backend."""
from __future__ import annotations

import numpy as np

from pysnspd.circuit.readout import CircuitParams, CircuitRuntimeConfig
from pysnspd.solver.stationary import solve_stationary_pytdgl_like


def test_solve_stationary_returns_relaxation_result(small_strip_mesh_bundle, gtdgl_material, stationary_seed_factory):
    mesh, edge_data, ops = small_strip_mesh_bundle
    seed = stationary_seed_factory(mesh, gtdgl_material)
    result = solve_stationary_pytdgl_like(
        mesh=mesh,
        edge_data=edge_data,
        seed=seed,
        material=gtdgl_material,
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
    for key in (
        "snapshot_t_s",
        "phi_snapshot_t_s",
        "supercurrent_density_snapshot_A_m2",
        "normal_current_density_snapshot_A_m2",
        "edge_i",
        "edge_j",
        "pytdgl_like_poisson_residual_rel",
        "pytdgl_like_poisson_residual_snapshot",
        "pytdgl_like_native_supercurrent_snapshot",
        "pytdgl_like_native_si_current_scale_A_m2",
        "pytdgl_like_native_si_residual_plus_boundary_rms_A_m3",
    ):
        assert key in result.history
    assert result.history["psi_snapshot_real_J"].shape == (2, mesh.n_nodes)
    assert result.history["edge_js_us_snapshot_A_m2"].shape == (2, ops.n_edges)
    assert result.history["pytdgl_like_poisson_residual_snapshot"].shape == (2, mesh.n_nodes)
    assert result.history["pytdgl_like_native_supercurrent_snapshot"].shape == (2, ops.n_edges)
    assert np.isfinite(result.summary["native_poisson_residual_rel_final"])
    assert result.summary["native_si_current_scale_A_m2"] > 0.0


def test_adapter_keeps_terminal_currents_in_amperes(small_strip_mesh_bundle, gtdgl_material, stationary_seed_factory):
    mesh, edge_data, ops = small_strip_mesh_bundle
    seed = stationary_seed_factory(mesh, gtdgl_material)
    target_current_A = 3.5e-6
    result = solve_stationary_pytdgl_like(
        mesh=mesh,
        edge_data=edge_data,
        seed=seed,
        material=gtdgl_material,
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
    assert result.summary["terminal_neumann_current_unit_A"] > 0.0
    assert result.summary["native_si_current_scale_A_m2"] > 0.0
    assert "native_si_boundary_currents_from_total_A" in result.summary


def test_stationary_adapter_serializes_circuit_after_final_step(
    small_strip_mesh_bundle,
    gtdgl_material,
    stationary_seed_factory,
):
    mesh, edge_data, ops = small_strip_mesh_bundle
    seed = stationary_seed_factory(mesh, gtdgl_material)
    result = solve_stationary_pytdgl_like(
        mesh=mesh,
        edge_data=edge_data,
        seed=seed,
        material=gtdgl_material,
        ops=ops,
        steps=2,
        dt_s=1.0e-18,
        target_current_A=0.0,
        terminal_psi=None,
        adaptive=False,
        n_snapshots=2,
        circuit_enabled=True,
        circuit_params=CircuitParams(),
        circuit_runtime_config=CircuitRuntimeConfig(start_time_s=0.0),
    )

    circuit = result.summary["circuit_runtime"]
    assert circuit["enabled"] is True
    assert circuit["params"]["R_load_ohm"] == 50.0
    assert circuit["params"]["V_bias_V"] == 0.0
    assert circuit["final_state"]["I_s_A"] == 0.0
    state_circuit = result.state.metadata["circuit_runtime"]
    assert state_circuit["params"] == circuit["params"]
    assert state_circuit["final_state"] == circuit["final_state"]


def test_adaptive_solver_stops_exactly_at_requested_chunk_boundary(
    small_strip_mesh_bundle,
    gtdgl_material,
    stationary_seed_factory,
):
    mesh, edge_data, ops = small_strip_mesh_bundle
    seed = stationary_seed_factory(mesh, gtdgl_material)
    requested_time_s = 3.25e-18
    result = solve_stationary_pytdgl_like(
        mesh=mesh,
        edge_data=edge_data,
        seed=seed,
        material=gtdgl_material,
        ops=ops,
        total_time_s=requested_time_s,
        dt_s=1.0e-18,
        target_current_A=0.0,
        terminal_psi=None,
        adaptive=True,
        n_snapshots=2,
    )

    assert np.isclose(
        result.summary["final_time_ps"],
        requested_time_s / 1.0e-12,
        rtol=1.0e-12,
        atol=1.0e-18,
    )


def test_stationary_adapter_builds_total_stationarity_gate(
    small_strip_mesh_bundle,
    gtdgl_material,
    stationary_seed_factory,
):
    mesh, edge_data, ops = small_strip_mesh_bundle
    seed = stationary_seed_factory(mesh, gtdgl_material)
    result = solve_stationary_pytdgl_like(
        mesh=mesh,
        edge_data=edge_data,
        seed=seed,
        material=gtdgl_material,
        ops=ops,
        steps=2,
        dt_s=1.0e-18,
        target_current_A=0.0,
        terminal_psi=None,
        adaptive=False,
        n_snapshots=2,
        stop_on_total_stationarity=True,
        stationarity_phase_gradient_rel=0.123,
        stationarity_phi_gradient_rel=0.234,
    )

    assert result.summary["accepted_steps"] >= 1
    assert "photon_ready" in result.summary


def test_adapter_rejects_usadel_poisson_without_catalog(small_strip_mesh_bundle, gtdgl_material, stationary_seed_factory):
    mesh, edge_data, ops = small_strip_mesh_bundle
    seed = stationary_seed_factory(mesh, gtdgl_material, q0_m_inv=1.0e7)
    try:
        solve_stationary_pytdgl_like(
            mesh=mesh,
            edge_data=edge_data,
            seed=seed,
            material=gtdgl_material,
            ops=ops,
            steps=1,
            dt_s=1.0e-18,
            target_current_A=0.0,
            terminal_psi=None,
            adaptive=False,
            n_snapshots=2,
            supercurrent_law="usadel_poisson",
        )
    except RuntimeError as exc:
        assert "Usadel" in str(exc) or "catalogue" in str(exc)
    else:
        raise AssertionError("usadel_poisson without a catalogue must fail loudly")
