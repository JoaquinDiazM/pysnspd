import numpy as np

from pysnspd.kinetic.eliashberg import build_debye_reference_spectrum
from pysnspd.kinetic.phase_space import PhaseSpaceCatalog
from pysnspd.kinetic.powers import (
    TAU0_OVER_TAU_EP_TC,
    bose_difference,
    compute_power_curve,
    compute_projected_powers,
    compute_vodolazov_debye_power_density,
    diagnostic_bcs_gap_factor,
    electronic_density_of_states_from_sigma_D,
    tau0_from_tau_ep_Tc,
    tau_ep_Tc_from_tau0,
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
    omega = np.linspace(0.0, 20.0, 80) * 1.602176634e-22
    Te_axis = np.array([1.0, 10.0])
    delta_axis = np.array([0.0, 1.0e-22])
    q_axis = np.array([0.0])

    shape = (Te_axis.size, delta_axis.size, q_axis.size, omega.size)
    JS = np.ones(shape) * 1.0e-22
    JR = np.ones(shape) * 1.0e-22
    JR[:, 0, :, :] = 0.0

    phase = PhaseSpaceCatalog(
        Te_values_K=Te_axis,
        omega_values_J=omega,
        delta_values_J=delta_axis,
        gamma_values_J=np.array([0.0]),
        q_values_m_inv=q_axis,
        J_S_TdqO_J=JS,
        J_R_TdqO_J=JR,
        delta_indices=np.array([0, 1]),
        q_indices=np.array([0]),
        metadata={"backend": "test"},
    )

    spectrum = build_debye_reference_spectrum(
        omega,
        lambda_ep=1.0,
        omega_D_J=float(omega[-1]),
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
    omega = np.linspace(0.0, 20.0, 80) * 1.602176634e-22
    Te_axis = np.array([1.0, 20.0])
    delta_axis = np.array([0.0])
    q_axis = np.array([0.0])

    shape = (Te_axis.size, delta_axis.size, q_axis.size, omega.size)
    JS = np.ones(shape) * 1.0e-22
    JR = np.zeros(shape)

    phase = PhaseSpaceCatalog(
        Te_values_K=Te_axis,
        omega_values_J=omega,
        delta_values_J=delta_axis,
        gamma_values_J=np.array([0.0]),
        q_values_m_inv=q_axis,
        J_S_TdqO_J=JS,
        J_R_TdqO_J=JR,
        delta_indices=np.array([0]),
        q_indices=np.array([0]),
        metadata={"backend": "test"},
    )

    spectrum = build_debye_reference_spectrum(
        omega,
        lambda_ep=1.0,
        omega_D_J=float(omega[-1]),
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


def test_bcs_gap_factor_collapses_above_Tc():
    Te = np.array([1.0, 5.0, 10.0, 20.0])
    factor = diagnostic_bcs_gap_factor(Te, Tc_K=10.0)

    assert factor[0] > factor[1]
    assert factor[2] == 0.0
    assert factor[3] == 0.0


def test_bcs_gap_policy_suppresses_recombination_above_Tc():
    omega = np.linspace(0.0, 20.0, 80) * 1.602176634e-22
    Te_axis = np.array([1.0, 20.0])
    delta_axis = np.array([0.0, 1.0e-22])
    q_axis = np.array([0.0])

    shape = (Te_axis.size, delta_axis.size, q_axis.size, omega.size)
    JS = np.ones(shape) * 1.0e-22
    JR = np.ones(shape) * 2.0e-22
    JR[:, 0, :, :] = 0.0

    phase = PhaseSpaceCatalog(
        Te_values_K=Te_axis,
        omega_values_J=omega,
        delta_values_J=delta_axis,
        gamma_values_J=np.array([0.0]),
        q_values_m_inv=q_axis,
        J_S_TdqO_J=JS,
        J_R_TdqO_J=JR,
        delta_indices=np.array([0, 1]),
        q_indices=np.array([0]),
        metadata={"backend": "test"},
    )

    spectrum = build_debye_reference_spectrum(
        omega,
        lambda_ep=1.0,
        omega_D_J=float(omega[-1]),
    )

    Te_values = np.array([5.0, 12.0])
    delta_values = 1.0e-22 * diagnostic_bcs_gap_factor(Te_values, Tc_K=10.0)

    curve = compute_power_curve(
        Te_values,
        Tph_K=1.0,
        delta_J=1.0e-22,
        q_m_inv=0.0,
        phase_space_catalog=phase,
        spectrum=spectrum,
        N0_J_m3=5.0e47,
        tau0_s=1.0e-9,
        Tc_K=10.0,
        delta_values_J=delta_values,
    )

    assert curve["delta_values_J"][0] > 0.0
    assert curve["delta_values_J"][1] == 0.0
    assert curve["P_R_W_m3"][1] == 0.0