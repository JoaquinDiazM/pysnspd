from pathlib import Path

import numpy as np

from pysnspd.kinetic.eliashberg import (
    EliashbergSpectrum,
    build_debye_reference_spectrum,
    load_simon_eliashberg_dat,
    mev_to_j,
)
from pysnspd.kinetic.phase_space import PhaseSpaceCatalog
from pysnspd.kinetic.powers import (
    bose_difference,
    compute_projected_powers,
    compute_vodolazov_debye_power_density,
    electronic_density_of_states_from_sigma_D,
)


def test_density_of_states_from_sigma_D_is_positive():
    N0 = electronic_density_of_states_from_sigma_D(
        sigma_n_S_m=4.2e5,
        D_m2_s=1.58e-4,
    )
    assert np.isfinite(N0)
    assert N0 > 0.0


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