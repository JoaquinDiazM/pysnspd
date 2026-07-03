from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from pysnspd.kinetic.eliashberg import build_debye_reference_spectrum
from pysnspd.kinetic.phase_space import PhaseSpaceCatalog
from pysnspd.kinetic.power_table import (
    build_power_table_catalog,
    load_power_table_catalog_npz,
    power_table_summary,
    save_power_table_catalog_npz,
)


def _tiny_config() -> dict:
    return {
        "project": {"name": "test", "big_data_root": "/tmp", "default_run_name": "test"},
        "parallel": {"enabled": False, "workers": 1, "backend": "serial"},
        "material": {
            "name": "NbN",
            "Tc_K": 8.65,
            "sigma_n_S_m": 4.2e5,
            "lambda_L_m": 5.4e-7,
            "thickness_m": 7.0e-9,
            "width_m": 120.0e-9,
            "tau_ee_Tc_ps": 0.5,
            "tau_ep_Tc_ps": 2.47,
            "tau_esc_ps": 20.0,
            "ion_density_nm3": 48.0,
        },
        "calibration": {"Ic_target_A": 38.8e-6},
        "bias": {"T_bias_K": 0.9, "I_bias_A": 35.0e-6},
        "mesh": {"type": "delaunay", "target_spacing_m": 4.1e-9, "seed": 1},
        "catalogs": {
            "dos": {"n_delta": 2, "n_q": 2, "n_energy": 8, "n_matsubara": 20},
            "phase_space": {"n_Te": 2, "n_Tph": 2, "n_delta": 2, "n_q": 2, "n_omega": 16},
        },
        "ss_run": {"max_steps": 1, "dt_s": 1e-15, "convergence_tol": 1e-7},
        "photon_run": {"photon_wavelength_m": 1064e-9, "max_steps": 1, "dt_s": 1e-15, "bubble_radius_m": 10e-9},
        "circuit": {"R_load_ohm": 50.0, "L_bias_H": 1e-6, "C_rf_F": 1e-12},
    }


def test_power_table_catalog_smoke(tmp_path):
    Te = np.array([0.9, 3.0], dtype=float)
    Tph = np.array([0.9, 3.0], dtype=float)
    delta = np.array([0.0, 2.0e-22], dtype=float)
    gamma = np.array([0.0, 1.0e-23], dtype=float)
    q = np.array([0.0, 4.0e7], dtype=float)
    omega = np.linspace(0.0, 8.0e-22, 16)
    energy = np.linspace(0.0, 9.0e-22, 24)

    JS = np.ones((Te.size, delta.size, q.size, omega.size), dtype=float) * 1.0e-22
    JR = np.ones_like(JS) * 2.0e-23
    JR[:, 0, :, :] = 0.0

    phase = PhaseSpaceCatalog(
        Te_values_K=Te,
        omega_values_J=omega,
        delta_values_J=delta,
        gamma_values_J=gamma,
        q_values_m_inv=q,
        J_S_TdqO_J=JS,
        J_R_TdqO_J=JR,
        delta_indices=np.array([0, 1], dtype=np.int64),
        q_indices=np.array([0, 1], dtype=np.int64),
        metadata={"backend": "test"},
    )

    rho = np.ones((delta.size, q.size, energy.size), dtype=float)
    usadel = SimpleNamespace(
        energy_values_J=energy,
        delta_values_J=delta,
        q_values_m_inv=q,
        rho_delta_gamma_E=rho,
        calibration_delta_eq_values_J=delta,
        metadata={
            "D_m2_s": 1.58e-4,
            "sigma_n_S_m": 4.2e5,
            "delta0_J": 2.0e-22,
        },
    )

    spectrum = build_debye_reference_spectrum(
        omega,
        lambda_ep=0.1,
        omega_D_J=float(omega[-1]),
    )

    catalog = build_power_table_catalog(
        phase_space_catalog=phase,
        usadel_catalog=usadel,
        spectrum=spectrum,
        config=_tiny_config(),
        n_Tph=Tph.size,
        Tph_min_K=float(Tph[0]),
        Tph_max_K=float(Tph[-1]),
        workers=1,
        parallel_backend="serial",
        progress=False,
    )

    assert catalog.P_total_W_m3.shape == (2, 2, 2, 2)
    assert catalog.u_e_J_m3.shape == (2, 2, 2)
    assert catalog.kappa_s_W_m_K.shape == (2, 2)
    assert catalog.u_ph_J_m3.shape == (2,)
    assert catalog.C_ph_J_m3_K.shape == (2,)
    assert catalog.P_esc_W_m3.shape == (2,)
    assert np.all(np.isfinite(catalog.P_total_W_m3))
    assert np.all(np.isfinite(catalog.u_e_J_m3))
    assert np.all(np.isfinite(catalog.kappa_s_W_m_K))
    assert np.all(np.isfinite(catalog.u_ph_J_m3))
    assert np.all(np.isfinite(catalog.C_ph_J_m3_K))
    assert np.all(np.isfinite(catalog.P_esc_W_m3))
    assert catalog.kappa_s_W_m_K[0, 0] > catalog.kappa_s_W_m_K[0, 1]
    assert catalog.u_ph_J_m3[1] > catalog.u_ph_J_m3[0]
    assert abs(catalog.P_esc_W_m3[0]) < 1.0e-20
    assert catalog.P_esc_W_m3[1] > 0.0
    assert np.max(np.abs(catalog.P_total_W_m3[0, 0])) < 1.0e-20
    assert np.max(np.abs(catalog.P_total_W_m3[1, 1])) < 1.0e-20
    assert np.max(np.abs(catalog.P_total_W_m3[1, 0])) > 0.0

    summary = power_table_summary(catalog)
    assert summary["n_Te"] == 2
    assert summary["n_Tph"] == 2
    assert summary["P_total_is_finite"] is True
    assert summary["kappa_s_is_finite"] is True
    assert summary["P_esc_is_finite"] is True

    path = save_power_table_catalog_npz(catalog, tmp_path / "power_table_catalog.npz")
    loaded = load_power_table_catalog_npz(path)
    assert loaded.P_total_W_m3.shape == catalog.P_total_W_m3.shape
    assert np.allclose(loaded.P_total_W_m3, catalog.P_total_W_m3)
    assert np.allclose(loaded.kappa_s_W_m_K, catalog.kappa_s_W_m_K)
    assert np.allclose(loaded.u_ph_J_m3, catalog.u_ph_J_m3)
    assert np.allclose(loaded.P_esc_W_m3, catalog.P_esc_W_m3)
