from __future__ import annotations

from pysnspd.config import validate_config


def _base_config() -> dict:
    return {
        "project": {
            "name": "test",
            "default_run_name": "run",
            "big_data_root": "/tmp",
        },
        "parallel": {"enabled": False, "workers": 1, "backend": "serial"},
        "material": {
            "name": "NbN",
            "Tc_K": 8.65,
            "sigma_n_S_m": 4.2e5,
            "lambda_L_m": 5.4e-7,
            "thickness_m": 7.0e-9,
            "width_m": 120.0e-9,
            "D_m2_s": 1.58e-4,
        },
        "calibration": {"Ic_target_A": 38.8e-6},
        "bias": {"T_bias_K": 0.9, "I_bias_A": 35.0e-6},
        "mesh": {"type": "delaunay", "target_spacing_m": 4.0e-9, "seed": 1},
        "catalogs": {
            "dos": {"n_delta": 2, "n_q": 2, "n_energy": 8, "n_matsubara": 8},
            "phase_space": {"n_Te": 2, "n_Tph": 2, "n_delta": 2, "n_q": 2, "n_omega": 8},
        },
        "ss_run": {"total_time_ps": 20.0, "dt_s": 1.0e-15, "convergence_tol": 1.0e-7},
        "photon_run": {
            "photon_wavelength_m": 1064.0e-9,
            "max_steps": 10,
            "dt_s": 1.0e-15,
            "bubble_radius_m": 10.0e-9,
        },
        "circuit": {"R_load_ohm": 50.0, "L_bias_H": 1.0e-6, "C_rf_F": 1.0e-12},
    }


def test_ss_run_total_time_does_not_require_max_steps() -> None:
    cfg = validate_config(_base_config(), require_big_data_root_exists=False)
    assert cfg["ss_run"]["total_time_ps"] == 20.0
    assert cfg["ss_run"]["snapshots_per_ps"] == 10.0
    assert "max_steps" not in cfg["ss_run"]


def test_legacy_ss_run_max_steps_is_converted_to_physical_time() -> None:
    old = _base_config()
    old["ss_run"] = {"max_steps": 1000, "dt_s": 2.0e-15, "convergence_tol": 1.0e-7}
    cfg = validate_config(old, require_big_data_root_exists=False)
    assert cfg["ss_run"]["total_time_ps"] == 2.0
