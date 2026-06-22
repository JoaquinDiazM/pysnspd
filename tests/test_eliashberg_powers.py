from types import SimpleNamespace

import numpy as np

from pysnspd.kinetic.eliashberg import build_debye_reference_spectrum
from pysnspd.kinetic.phase_space import PhaseSpaceCatalog
from pysnspd.kinetic.powers import (
    TAU0_OVER_TAU_EP_TC,
    bose_difference,
    build_usadel_self_consistent_trajectory,
    compute_power_curve_at_usadel_state,
    compute_projected_powers,
    compute_usadel_q_power_scan,
    compute_vodolazov_debye_power_density,
    electronic_density_of_states_from_sigma_D,
    select_usadel_state_by_current_fraction,
    tau0_from_tau_ep_Tc,
    tau_ep_Tc_from_tau0,
)


def _toy_phase_catalog():
    omega = np.linspace(0.0, 20.0, 80) * 1.602176634e-22
    Te_axis = np.array([1.0, 20.0])
    delta_axis = np.array([0.0, 1.0e-22])
    q_axis = np.array([0.0, 1.0])

    shape = (Te_axis.size, delta_axis.size, q_axis.size, omega.size)
    JS = np.ones(shape) * 1.0e-22
    JR = np.ones(shape) * 2.0e-22
    JR[:, 0, :, :] = 0.0

    return PhaseSpaceCatalog(
        Te_values_K=Te_axis,
        omega_values_J=omega,
        delta_values_J=delta_axis,
        gamma_values_J=np.array([0.0, 1.0e-22]),
        q_values_m_inv=q_axis,
        J_S_TdqO_J=JS,
        J_R_TdqO_J=JR,
        delta_indices=np.array([0, 1]),
        q_indices=np.array([0, 1]),
        metadata={"backend": "test"},
    )


def _toy_usadel_catalog():
    gamma = np.linspace(0.0, 1.0e-22, 6)
    q = np.linspace(0.0, 1.0, 6)
    delta = np.array([1.0, 0.95, 0.85, 0.70, 0.40, 0.10]) * 1.0e-22
    current = np.array([0.0, 1.0, 2.0, 2.5, 2.2, 1.0])
    current_density = current * 10.0

    return SimpleNamespace(
        calibration_gamma_values_J=gamma,
        calibration_q_values_m_inv=q,
        calibration_delta_eq_values_J=delta,
        calibration_current_values_A=current,
        calibration_current_density_values_A_m2=current_density,
        metadata={"T_bias_K": 0.9},
    )


def test_density_of_states_from_sigma_D_is_positive():
    N0 = electronic_density_of_states_from_sigma_D(
        sigma_n_S_m=4.2e5,
        D_m2_s=1.58e-4,
    )
    assert np.isfinite(N0)
    assert N0 > 0.0


def test_tau0_is_not_tau_ep_Tc():
    tau_ep = 24.7e-12
    tau0 = tau0_from_tau_ep_Tc(tau_ep)

    assert np.isclose(tau0 / tau_ep, TAU0_OVER_TAU_EP_TC)
    assert tau0 > 70.0 * tau_ep
    assert np.isclose(tau_ep_Tc_from_tau0(tau0), tau_ep)


def test_bose_difference_equilibrium_is_zero():
    omega = np.linspace(0.0, 10.0, 100) * 1.602176634e-22
    diff = bose_difference(omega, 5.0, 5.0)
    assert np.all(np.isfinite(diff))
    assert np.allclose(diff, 0.0)


def test_vodolazov_debye_power_has_T5_scaling():
    N0 = 5.0e47
    tau0_s = 25.0e-12
    Tc_K = 10.0
    Tph_K = 1.0

    p1 = compute_vodolazov_debye_power_density(
        5.0,
        Tph_K,
        N0_J_m3=N0,
        tau0_s=tau0_s,
        Tc_K=Tc_K,
    )
    p2 = compute_vodolazov_debye_power_density(
        10.0,
        Tph_K,
        N0_J_m3=N0,
        tau0_s=tau0_s,
        Tc_K=Tc_K,
    )

    expected = (10.0**5 - Tph_K**5) / (5.0**5 - Tph_K**5)
    assert np.isclose(p2 / p1, expected)


def test_projected_power_vanishes_at_thermal_equilibrium():
    phase = _toy_phase_catalog()
    spectrum = build_debye_reference_spectrum(
        phase.omega_values_J,
        lambda_ep=1.0,
        omega_D_J=float(phase.omega_values_J[-1]),
    )

    result = compute_projected_powers(
        5.0,
        5.0,
        1.0e-22,
        0.0,
        phase,
        spectrum,
        N0_J_m3=5.0e47,
    )

    assert abs(result.P_S_W_m3) < 1.0e-30
    assert abs(result.P_R_W_m3) < 1.0e-30
    assert abs(result.P_total_W_m3) < 1.0e-30


def test_projected_power_sign_follows_temperature_bias():
    phase = _toy_phase_catalog()
    spectrum = build_debye_reference_spectrum(
        phase.omega_values_J,
        lambda_ep=1.0,
        omega_D_J=float(phase.omega_values_J[-1]),
    )

    hot_e = compute_projected_powers(
        10.0,
        2.0,
        0.0,
        0.0,
        phase,
        spectrum,
        N0_J_m3=5.0e47,
    )
    hot_ph = compute_projected_powers(
        2.0,
        10.0,
        0.0,
        0.0,
        phase,
        spectrum,
        N0_J_m3=5.0e47,
    )

    assert hot_e.P_total_W_m3 > 0.0
    assert hot_ph.P_total_W_m3 < 0.0


def test_usadel_self_consistent_trajectory_uses_stable_branch():
    usadel = _toy_usadel_catalog()

    traj = build_usadel_self_consistent_trajectory(usadel, n_q=10)

    assert traj.q_values_m_inv.size == 10
    assert np.max(traj.current_values_A) <= 2.5
    assert np.isclose(traj.metadata["Ic_A"], 2.5)
    assert np.isclose(traj.q_values_m_inv[-1], 0.6)
    assert np.all(np.diff(traj.q_values_m_inv) >= 0.0)


def test_select_usadel_state_by_current_fraction():
    usadel = _toy_usadel_catalog()
    traj = build_usadel_self_consistent_trajectory(usadel, n_q=20)

    state = select_usadel_state_by_current_fraction(traj, 0.8)

    assert 0.0 <= state["current_fraction"] <= 1.0
    assert state["delta_J"] > 0.0
    assert state["q_m_inv"] >= 0.0


def test_usadel_q_power_scan_shapes():
    phase = _toy_phase_catalog()
    usadel = _toy_usadel_catalog()
    traj = build_usadel_self_consistent_trajectory(usadel, n_q=7)
    spectrum = build_debye_reference_spectrum(
        phase.omega_values_J,
        lambda_ep=1.0,
        omega_D_J=float(phase.omega_values_J[-1]),
    )

    out = compute_usadel_q_power_scan(
        np.array([2.0, 5.0]),
        Tph_K=1.0,
        trajectory=traj,
        phase_space_catalog=phase,
        spectrum=spectrum,
        N0_J_m3=5.0e47,
    )

    assert out["P_S_W_m3"].shape == (2, 7)
    assert out["P_R_W_m3"].shape == (2, 7)
    assert out["P_total_W_m3"].shape == (2, 7)
    assert np.all(np.isfinite(out["P_total_W_m3"]))


def test_power_curve_at_usadel_state_shapes():
    phase = _toy_phase_catalog()
    usadel = _toy_usadel_catalog()
    traj = build_usadel_self_consistent_trajectory(usadel, n_q=7)
    state = select_usadel_state_by_current_fraction(traj, 0.5)

    spectrum = build_debye_reference_spectrum(
        phase.omega_values_J,
        lambda_ep=1.0,
        omega_D_J=float(phase.omega_values_J[-1]),
    )

    out = compute_power_curve_at_usadel_state(
        np.array([2.0, 5.0]),
        Tph_K=1.0,
        state=state,
        phase_space_catalog=phase,
        spectrum=spectrum,
        N0_J_m3=5.0e47,
        tau0_s=1.0e-9,
        Tc_K=10.0,
    )

    assert out["P_total_W_m3"].shape == (2,)
    assert np.all(out["delta_values_J"] == state["delta_J"])
    assert np.all(out["q_values_m_inv"] == state["q_m_inv"])