from types import SimpleNamespace

import pytest

from pysnspd.config import ConfigError, validate_config
from pysnspd.gtdgl.material import build_gtdgl_material


def _minimal_config(material_extra=None):
    material = {
        "name": "NbN",
        "Tc_K": 8.65,
        "sigma_n_S_m": 4.2e5,
        "lambda_L_m": 5.4e-7,
        "thickness_m": 7.0e-9,
        "width_m": 120.0e-9,
    }
    if material_extra:
        material.update(material_extra)
    return {
        "project": {
            "name": "test",
            "big_data_root": "/tmp",
            "default_run_name": "run",
        },
        "parallel": {"enabled": True, "workers": 1, "backend": "process"},
        "material": material,
        "calibration": {"Ic_target_A": 38.8e-6},
        "bias": {"T_bias_K": 0.9, "I_bias_A": 20.0e-6},
        "mesh": {"type": "delaunay", "target_spacing_m": 4.0e-9, "seed": 1},
        "catalogs": {
            "dos": {"n_delta": 3, "n_q": 4, "n_energy": 8, "n_matsubara": 5},
            "phase_space": {"n_Te": 3, "n_Tph": 3, "n_delta": 3, "n_q": 4, "n_omega": 8},
        },
        "ss_run": {"max_steps": 10, "dt_s": 1.0e-15, "convergence_tol": 1.0e-7},
        "photon_run": {"photon_wavelength_m": 1064.0e-9, "max_steps": 10, "dt_s": 1.0e-15, "bubble_radius_m": 10.0e-9},
        "circuit": {"R_load_ohm": 50.0, "L_bias_H": 1.0e-6, "C_rf_F": 1.0e-12},
    }


def test_material_relaxation_times_are_validated_from_yaml_ps():
    cfg = validate_config(
        _minimal_config({"tau_ee_Tc_ps": 6.0, "tau_ep_Tc_ps": 31.0}),
        require_big_data_root_exists=False,
    )
    material = build_gtdgl_material(
        cfg,
        SimpleNamespace(metadata={"D_m2_s": 1.58e-4}),
    )

    assert material.tau_ee_Tc_s == pytest.approx(6.0e-12)
    assert material.tau_ep_Tc_s == pytest.approx(31.0e-12)


def test_material_relaxation_time_duplicate_aliases_are_rejected():
    with pytest.raises(ConfigError, match="multiple aliases"):
        validate_config(
            _minimal_config({"tau_ee_Tc_ps": 5.0, "tau_ee_Tc_s": 5.0e-12}),
            require_big_data_root_exists=False,
        )


def test_material_relaxation_times_are_used_without_hidden_multiplier():
    cfg = validate_config(
        _minimal_config({"tau_ee_Tc_ps": 0.5, "tau_ep_Tc_ps": 2.47}),
        require_big_data_root_exists=False,
    )
    material = build_gtdgl_material(cfg, SimpleNamespace(metadata={"D_m2_s": 1.58e-4}))

    assert "tau_" + "scale" not in getattr(material, "__dataclass_fields__", {})
    assert material.tau_ee_s(8.65) == pytest.approx(0.5e-12)
    assert material.tau_ep_s(8.65) == pytest.approx(2.47e-12)
