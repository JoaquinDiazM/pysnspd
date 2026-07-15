"""Tests for stationary-state target diagnostics."""
from __future__ import annotations

import numpy as np

from pysnspd.gtdgl.material import GTDGLMaterial
from pysnspd.gtdgl.ss_targets import (
    apply_terminal_proximity_seed,
    contact_recovery_diagnostics,
    dynamic_stationarity_diagnostics,
    stationarity_diagnostics,
)


def _material() -> GTDGLMaterial:
    return GTDGLMaterial(
        Tc_K=8.65,
        D_m2_s=1.58e-4,
        sigma_n_S_m=4.2e5,
        delta0_J=2.10667708314e-22,
        thickness_m=7e-9,
        width_m=120e-9,
        tau_ee_Tc_s=0.5e-12,
        tau_ep_Tc_s=2.47e-12,
    )


def test_terminal_proximity_seed_reaches_bulk_after_few_xi():
    material = _material()
    x = np.linspace(0.0, 240e-9, 401)
    nodes = np.column_stack([x, np.zeros_like(x)])
    psi = np.ones(x.size, dtype=np.complex128)
    Te = np.full(x.size, 0.9)

    out, seed_diag = apply_terminal_proximity_seed(
        psi,
        nodes_m=nodes,
        material=material,
        Te_K=Te,
        healing_target_xi=2.5,
        target_bulk_fraction=0.95,
    )
    rec = contact_recovery_diagnostics(
        psi_dimensionless=out,
        nodes_m=nodes,
        material=material,
        Te_K=Te,
        threshold_fraction=0.95,
        min_allowed_xi=1.5,
        max_allowed_xi=4.0,
        bin_width_m=0.5e-9,
    )

    assert seed_diag.enabled
    assert np.isclose(np.abs(out[0]), 0.0)
    assert np.max(np.abs(out)) <= 1.0 + 1e-12
    assert rec.passes
    assert 1.5 <= rec.left_recovery_length_xi <= 4.0
    assert 1.5 <= rec.right_recovery_length_xi <= 4.0


def test_stationarity_diagnostics_uses_gradients_not_global_fields():
    material = _material()
    delta0 = material.delta0_J
    edge_i = np.array([0, 1, 2], dtype=np.int64)
    edge_j = np.array([1, 2, 3], dtype=np.int64)
    edge_length_m = np.full(3, 1.0e-9)

    # Same physical edge Q and same grad(phi), but different global phase and
    # constant potential offset. This must pass because those offsets are gauge
    # choices in the present A=0 gauge-fixed backend.
    psi0 = np.ones(4, dtype=np.complex128) * delta0
    psi1 = psi0 * np.exp(1j * 0.7)
    phi0 = np.array([0.0, 1.0e-6, 2.0e-6, 3.0e-6])
    phi1 = phi0 + 123.0
    history = {
        "psi_snapshot_real_J": np.vstack([np.real(psi0), np.real(psi1)]),
        "psi_snapshot_imag_J": np.vstack([np.imag(psi0), np.imag(psi1)]),
        "phi_snapshot_V": np.vstack([phi0, phi1]),
        "edge_Q_snapshot_m_inv": np.vstack(
            [
                np.array([1.0e7, 1.1e7, 1.2e7]),
                np.array([1.0e7, 1.1e7, 1.2e7]),
            ]
        ),
        "edge_i": edge_i,
        "edge_j": edge_j,
        "edge_length_m": edge_length_m,
        "eta_R": np.array([1.0e-3, 1.0e-2]),  # info-only now
    }

    diag = stationarity_diagnostics(
        history=history,
        material=material,
        phase_gradient_rel_tol=1.0e-6,
        phi_gradient_rel_tol=1.0e-6,
        phase_gradient_abs_tol_m_inv=1.0,
        phi_gradient_abs_tol_V_m=1.0,
    )

    assert diag.passes
    assert diag.phase_gradient_rel_change == 0.0
    assert diag.phi_gradient_rel_change < 1.0e-8
    assert diag.eta_R_window_max > 1.0e-3


def test_stationarity_diagnostics_detects_changing_phase_gradient():
    material = _material()
    delta0 = material.delta0_J
    history = {
        "psi_snapshot_real_J": np.ones((2, 4)) * delta0,
        "psi_snapshot_imag_J": np.zeros((2, 4)),
        "phi_snapshot_V": np.zeros((2, 4)),
        "edge_Q_snapshot_m_inv": np.vstack(
            [
                np.array([1.0e7, 1.0e7, 1.0e7]),
                np.array([1.2e7, 1.0e7, 0.8e7]),
            ]
        ),
        "edge_i": np.array([0, 1, 2], dtype=np.int64),
        "edge_j": np.array([1, 2, 3], dtype=np.int64),
        "edge_length_m": np.full(3, 1.0e-9),
        "eta_R": np.array([1.0e-8, 1.0e-8]),
    }

    diag = stationarity_diagnostics(
        history=history,
        material=material,
        phase_gradient_rel_tol=1.0e-4,
        phi_gradient_rel_tol=1.0e-4,
        phase_gradient_abs_tol_m_inv=1.0e3,
        phi_gradient_abs_tol_V_m=1.0,
    )

    assert not diag.passes
    assert diag.phase_gradient_abs_change_m_inv > 1.0e6


def test_stationarity_diagnostics_excludes_zero_delta_terminal_edges():
    material = _material()
    delta0 = material.delta0_J
    psi0 = np.array([0.0, 0.0, 1.0, 1.0], dtype=np.complex128) * delta0
    psi1 = psi0.copy()
    history = {
        "psi_snapshot_real_J": np.vstack([np.real(psi0), np.real(psi1)]),
        "psi_snapshot_imag_J": np.vstack([np.imag(psi0), np.imag(psi1)]),
        "phi_snapshot_V": np.zeros((2, 4)),
        "edge_Q_snapshot_m_inv": np.vstack(
            [
                np.array([9.0e99, 1.0e7, 1.0e7]),
                np.array([-9.0e99, 1.0e7, 1.0e7]),
            ]
        ),
        "edge_i": np.array([0, 1, 2], dtype=np.int64),
        "edge_j": np.array([1, 2, 3], dtype=np.int64),
        "edge_length_m": np.full(3, 1.0e-9),
        "normal_terminal_edge_mask": np.array([True, True, False]),
        "eta_R": np.array([1.0e-8, 1.0e-8]),
    }

    diag = stationarity_diagnostics(
        history=history,
        material=material,
        phase_gradient_rel_tol=1.0e-6,
        phi_gradient_rel_tol=1.0e-6,
    )

    assert diag.passes
    assert diag.active_edge_count == 1


def test_stationarity_diagnostics_excludes_contact_conversion_region():
    material = _material()
    delta0 = material.delta0_J

    # Edge 0 is inside the metallic-contact conversion region and changes a lot.
    # Edges 1 and 2 are bulk edges and remain stationary.
    history = {
        "psi_snapshot_real_J": np.ones((2, 4)) * delta0,
        "psi_snapshot_imag_J": np.zeros((2, 4)),
        "phi_snapshot_V": np.zeros((2, 4)),
        "edge_Q_snapshot_m_inv": np.vstack(
            [
                np.array([1.0e7, 1.0e7, 1.0e7]),
                np.array([5.0e8, 1.0e7, 1.0e7]),
            ]
        ),
        "edge_phi_gradient_snapshot_V_m": np.vstack(
            [
                np.array([0.0, 10.0, 12.0]),
                np.array([2.0e5, 10.0, 12.0]),
            ]
        ),
        "edge_i": np.array([0, 1, 2], dtype=np.int64),
        "edge_j": np.array([1, 2, 3], dtype=np.int64),
        "edge_length_m": np.full(3, 1.0e-9),
        "edge_distance_from_contact_m": np.array([5.0e-9, 50.0e-9, 60.0e-9]),
        "stationarity_xi_m": np.array([10.0e-9]),
        "eta_R": np.array([1.0e-8, 1.0e-8]),
    }

    bulk_diag = stationarity_diagnostics(
        history=history,
        material=material,
        phase_gradient_rel_tol=1.0e-6,
        phi_gradient_rel_tol=1.0e-6,
        phase_gradient_abs_tol_m_inv=1.0e3,
        phi_gradient_abs_tol_V_m=1.0,
        bulk_exclusion_xi=4.0,
    )

    assert bulk_diag.passes
    assert bulk_diag.active_edge_count == 2
    assert np.isclose(bulk_diag.bulk_exclusion_length_m, 40.0e-9)

    full_domain_diag = stationarity_diagnostics(
        history=history,
        material=material,
        phase_gradient_rel_tol=1.0e-6,
        phi_gradient_rel_tol=1.0e-6,
        phase_gradient_abs_tol_m_inv=1.0e3,
        phi_gradient_abs_tol_V_m=1.0,
        bulk_exclusion_xi=0.0,
    )

    assert not full_domain_diag.passes
    assert full_domain_diag.active_edge_count == 3


def _dynamic_history(*, add_final_band: bool) -> tuple[dict[str, np.ndarray], np.ndarray]:
    material = _material()
    x_axis = np.linspace(0.0, 200.0e-9, 81)
    y_axis = np.linspace(-60.0e-9, 60.0e-9, 9)
    xx, yy = np.meshgrid(x_axis, y_axis, indexing="xy")
    nodes = np.column_stack([xx.reshape(-1), yy.reshape(-1)])
    snapshot_t_s = np.linspace(0.0, 4.0e-12, 5)
    frames = []
    for index, _ in enumerate(snapshot_t_s):
        profile = np.ones_like(x_axis)
        profile -= 0.65 * np.exp(-0.5 * ((x_axis - 75.0e-9) / 5.0e-9) ** 2)
        profile -= 0.60 * np.exp(-0.5 * ((x_axis - 130.0e-9) / 5.0e-9) ** 2)
        profile *= 1.0 + 1.0e-3 * np.sin(index)
        if add_final_band and index == len(snapshot_t_s) - 1:
            profile -= 0.65 * np.exp(-0.5 * ((x_axis - 102.0e-9) / 4.0e-9) ** 2)
        frames.append(np.tile(profile, y_axis.size))
    amplitude = np.asarray(frames) * material.delta0_J
    history_t_s = np.linspace(0.0, 4.0e-12, 101)
    voltage = 10.0e-3 * (1.0 + 2.0e-3 * np.sin(2.0 * np.pi * history_t_s / 0.4e-12))
    history = {
        "snapshot_t_s": snapshot_t_s,
        "psi_snapshot_real_J": amplitude,
        "psi_snapshot_imag_J": np.zeros_like(amplitude),
        "stationarity_xi_m": np.array([10.0e-9]),
        "t_s": history_t_s,
        "terminal_voltage_V": voltage,
    }
    return history, nodes


def test_dynamic_stationarity_accepts_stable_psl_morphology_with_small_oscillation():
    material = _material()
    history, nodes = _dynamic_history(add_final_band=False)
    diag = dynamic_stationarity_diagnostics(
        history=history,
        nodes_m=nodes,
        delta0_J=material.delta0_J,
        tail_snapshots=4,
        minimum_tail_duration_ps=2.0,
    )

    assert diag.passes
    assert diag.topology_count_stable
    assert not diag.new_suppressed_band_in_tail
    assert diag.psl_count_final == 2
    assert diag.voltage_relative_span < diag.tolerance_voltage_relative


def test_dynamic_stationarity_rejects_new_psl_in_tail():
    material = _material()
    history, nodes = _dynamic_history(add_final_band=True)
    diag = dynamic_stationarity_diagnostics(
        history=history,
        nodes_m=nodes,
        delta0_J=material.delta0_J,
        tail_snapshots=4,
        minimum_tail_duration_ps=2.0,
    )

    assert not diag.passes
    assert diag.new_suppressed_band_in_tail
    assert not diag.topology_count_stable
