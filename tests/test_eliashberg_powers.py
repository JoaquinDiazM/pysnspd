from types import SimpleNamespace

import numpy as np

from pysnspd.kinetic.eliashberg import build_debye_reference_spectrum
from pysnspd.kinetic.phase_space import PhaseSpaceCatalog
from pysnspd.kinetic.powers import (
    TAU0_OVER_TAU_EP_TC,
    bose_difference,
    build_thermal_usadel_grid,
    compute_power_curve_thermal_usadel_state,
    compute_power_scan_thermal_usadel,
    compute_projected_powers,
    compute_vodolazov_debye_power_density,
    electronic_density_of_states_from_sigma_D,
    select_thermal_usadel_q_state,
    tau0_from_tau_ep_Tc,
    tau_ep_Tc_from_tau0,
    thermal_usadel_delta_at_state,
)


def _toy_phase_catalog():
    omega = np.linspace(0.0, 20.0, 80) * 1.602176634e-22
    Te_axis = np.array([0.9, 5.0, 10.0, 20.0])
    delta_axis = np.array([0.0, 1.0e-22, 3.0e-22])
    q_axis = np.array([0.0, 0.5, 1.0])

    shape = (Te_axis.size, delta_axis.size, q_axis.size, omega.size)
    JS = np.ones(shape) * 1.0e-22
    JR = np.ones(shape) * 2.0e-22
    JR[:, 0, :, :] = 0.0

    return PhaseSpaceCatalog(
        Te_values_K=Te_axis,
        omega_values_J=omega,
        delta_values_J=delta_axis,
        gamma_values_J=np.array([0.0, 0.5e-22, 1.0e-22]),
        q_values_m_inv=q_axis,
        J_S_TdqO_J=JS,
        J_R_TdqO_J=JR,
        delta_indices=np.array([0, 1, 2]),
        q_indices=np.array([0, 1, 2]),
        metadata={"backend": "test"},
    )


def _toy_usadel_catalog():
    q = np.linspace(0.0, 1.0, 6)
    gamma = 0.5 * 1.054571817e-34 * 1.0e-4 * q * q
    delta = np.array([1.0, 0.95, 0.85, 0.70, 0.40, 0.10]) * 2.0e-22
    current = np.array([0.0, 1.0, 2.0, 2.5, 2.2, 1.0])
    current_density = current * 10.0

    return SimpleNamespace(
        calibration_gamma_values_J=gamma,
        calibration_q_values_m_inv=q,
        calibration_delta_eq_values_J=delta,
        calibration_current_values_A=current,
        calibration_current_density_values_A_m2=current_density,
        metadata={
            "T_bias_K": 0.9,
            "Tc_K": 10.0,
            "D_m2_s": 1.0e-4,
            "sigma_n_S_m": 4.2e5,
            "width_m": 120e-9,
            "thickness_m": 7e-9,
        },
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


def test_thermal_usadel_grid_collapses_gap_above_Tc():
    usadel = _toy_usadel_catalog()
    grid = build_thermal_usadel_grid(
        usadel,
        np.array([1.0, 5.0, 12.0]),
        n_q=5,
        n_matsubara=80,
    )

    assert grid.delta_eq_Tq_J.shape == (3, 5)
    assert np.max(grid.delta_eq_Tq_J[0, :]) > 0.0
    assert np.max(grid.delta_eq_Tq_J[-1, :]) == 0.0


def test_select_thermal_usadel_q_state():
    usadel = _toy_usadel_catalog()
    grid = build_thermal_usadel_grid(
        usadel,
        np.array([1.0, 5.0, 12.0]),
        n_q=8,
        n_matsubara=80,
    )

    state = select_thermal_usadel_q_state(grid, 0.8)

    assert 0.0 <= state["reference_current_fraction"] <= 1.0
    assert state["q_m_inv"] >= 0.0
    assert state["gamma_J"] >= 0.0


def test_thermal_usadel_delta_interpolation():
    usadel = _toy_usadel_catalog()
    grid = build_thermal_usadel_grid(
        usadel,
        np.array([1.0, 5.0, 12.0]),
        n_q=8,
        n_matsubara=80,
    )

    d_low = thermal_usadel_delta_at_state(
        grid,
        Te_K=1.0,
        q_m_inv=float(grid.q_values_m_inv[2]),
    )
    d_high = thermal_usadel_delta_at_state(
        grid,
        Te_K=12.0,
        q_m_inv=float(grid.q_values_m_inv[2]),
    )

    assert d_low >= 0.0
    assert d_high == 0.0


def test_power_scan_thermal_usadel_shapes():
    phase = _toy_phase_catalog()
    usadel = _toy_usadel_catalog()
    grid = build_thermal_usadel_grid(
        usadel,
        np.array([1.0, 5.0, 12.0]),
        n_q=7,
        n_matsubara=80,
    )
    spectrum = build_debye_reference_spectrum(
        phase.omega_values_J,
        lambda_ep=1.0,
        omega_D_J=float(phase.omega_values_J[-1]),
    )

    out = compute_power_scan_thermal_usadel(
        np.array([2.0, 5.0]),
        Tph_K=1.0,
        thermal_grid=grid,
        phase_space_catalog=phase,
        spectrum=spectrum,
        N0_J_m3=5.0e47,
    )

    assert out["P_S_W_m3"].shape == (2, 7)
    assert out["P_R_W_m3"].shape == (2, 7)
    assert out["P_total_W_m3"].shape == (2, 7)
    assert np.all(np.isfinite(out["P_total_W_m3"]))


def test_power_curve_thermal_usadel_state_shapes():
    phase = _toy_phase_catalog()
    usadel = _toy_usadel_catalog()
    grid = build_thermal_usadel_grid(
        usadel,
        np.array([1.0, 5.0, 12.0]),
        n_q=7,
        n_matsubara=80,
    )
    state = select_thermal_usadel_q_state(grid, 0.5)

    spectrum = build_debye_reference_spectrum(
        phase.omega_values_J,
        lambda_ep=1.0,
        omega_D_J=float(phase.omega_values_J[-1]),
    )

    out = compute_power_curve_thermal_usadel_state(
        np.array([2.0, 5.0]),
        Tph_K=1.0,
        state=state,
        thermal_grid=grid,
        phase_space_catalog=phase,
        spectrum=spectrum,
        N0_J_m3=5.0e47,
        tau0_s=1.0e-9,
        Tc_K=10.0,
    )

    assert out["P_total_W_m3"].shape == (2,)
    assert out["P_normal_Eliashberg_W_m3"].shape == (2,)
    assert np.all(out["q_values_m_inv"] == state["q_m_inv"])