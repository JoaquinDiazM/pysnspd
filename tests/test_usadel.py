from pathlib import Path

import numpy as np

from pysnspd.config import validate_config
from pysnspd.usadel.catalog import (
    build_usadel_catalog_from_config,
    catalog_summary,
    save_usadel_catalog_npz,
    load_usadel_catalog_npz,
)
from pysnspd.usadel.parameters import (
    bcs_gap_J,
    bcs_gap_zero_J,
    depairing_energy_grid_J,
    q_axis_from_depairing_energy_m_inv,
)
from pysnspd.usadel.solver import dynes_bcs_dos


def _minimal_config(tmp_path: Path) -> dict:
    return {
        "project": {
            "name": "test_project",
            "big_data_root": str(tmp_path / "big_data"),
            "default_run_name": "usadel_test",
        },
        "parallel": {
            "enabled": True,
            "workers": 2,
            "backend": "process",
        },
        "material": {
            "name": "NbN",
            "Tc_K": 8.65,
            "D_m2_s": 1.58e-4,
            "sigma_n_S_m": "4.2e5",
            "lambda_L_m": 5.4e-7,
            "thickness_m": 7.0e-9,
            "width_m": 120.0e-9,
        },
        "bias": {
            "T_bias_K": 0.9,
            "I_bias_A": 35.0e-6,
        },
        "mesh": {
            "type": "delaunay",
            "length_m": 240.0e-9,
            "target_spacing_m": 20.0e-9,
            "seed": 12345,
        },
        "catalogs": {
            "dos": {
                "n_delta": 6,
                "n_q": 5,
                "n_energy": 80,
                "n_matsubara": 20,
            },
            "phase_space": {
                "n_Te": 5,
                "n_Tph": 5,
                "n_delta": 6,
                "n_q": 5,
                "n_omega": 50,
            },
        },
        "ss_run": {
            "max_steps": 10,
            "dt_s": 1.0e-15,
            "convergence_tol": 1.0e-7,
        },
        "photon_run": {
            "photon_wavelength_m": 1064.0e-9,
            "max_steps": 10,
            "dt_s": 1.0e-15,
            "bubble_radius_m": 10.0e-9,
        },
        "circuit": {
            "R_load_ohm": 50.0,
            "L_bias_H": 1.0e-6,
            "C_rf_F": 1.0e-12,
        },
    }


def test_bcs_gap_scales_are_reasonable():
    Tc = 8.65
    delta0 = bcs_gap_zero_J(Tc)

    assert delta0 > 0.0
    assert bcs_gap_J(0.0, Tc) == delta0
    assert bcs_gap_J(Tc, Tc) == 0.0
    assert 0.0 < bcs_gap_J(0.9, Tc) <= delta0


def test_depairing_and_q_axes():
    delta0 = bcs_gap_zero_J(8.65)
    gammas = depairing_energy_grid_J(delta_ref_J=delta0, n_q=5)
    q_values = q_axis_from_depairing_energy_m_inv(gammas, D_m2_s=1.58e-4)

    assert gammas.shape == (5,)
    assert q_values.shape == (5,)
    assert gammas[0] == 0.0
    assert q_values[0] == 0.0
    assert np.all(np.diff(gammas) >= 0.0)
    assert np.all(np.diff(q_values) >= 0.0)


def test_dynes_dos_normal_state_and_finite_values():
    E = np.linspace(0.0, 10.0, 100)
    rho_normal = dynes_bcs_dos(E, delta_J=0.0)

    assert np.allclose(rho_normal, 1.0)

    rho_sc = dynes_bcs_dos(E, delta_J=1.0, eta_J=1.0e-3)
    assert rho_sc.shape == E.shape
    assert np.all(np.isfinite(rho_sc))
    assert np.min(rho_sc) >= 0.0


def test_build_save_load_usadel_catalog(tmp_path):
    cfg = validate_config(
        _minimal_config(tmp_path),
        require_big_data_root_exists=False,
    )

    catalog = build_usadel_catalog_from_config(cfg)

    assert catalog.rho_delta_gamma_E.shape == (6, 5, 80)
    assert catalog.anomalous_delta_gamma_E.shape == (6, 5, 80)
    assert catalog.energy_values_J.shape == (80,)
    assert catalog.delta_values_J.shape == (6,)
    assert catalog.gamma_values_J.shape == (5,)
    assert catalog.q_values_m_inv.shape == (5,)
    assert np.all(np.isfinite(catalog.rho_delta_gamma_E))

    summary = catalog_summary(catalog)
    assert summary["shape"] == [6, 5, 80]
    assert summary["rho_is_finite"] is True

    path = save_usadel_catalog_npz(catalog, tmp_path / "usadel_dos_catalog.npz")
    loaded = load_usadel_catalog_npz(path)

    assert np.allclose(catalog.energy_values_J, loaded.energy_values_J)
    assert np.allclose(catalog.delta_values_J, loaded.delta_values_J)
    assert np.allclose(catalog.gamma_values_J, loaded.gamma_values_J)
    assert np.allclose(catalog.rho_delta_gamma_E, loaded.rho_delta_gamma_E)
    assert loaded.metadata["backend"] == "dynes_bcs_proxy_for_usadel_oe3"