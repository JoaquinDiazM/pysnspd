from types import SimpleNamespace

import numpy as np
import pytest

from pysnspd.gtdgl.material import build_gtdgl_material


def _config():
    return {
        "material": {
            "Tc_K": 8.65,
            "D_m2_s": 1.58e-4,
            "sigma_n_S_m": 4.2e5,
            "thickness_m": 7.0e-9,
            "width_m": 120.0e-9,
            "tau_ee_Tc_ps": 0.5,
            "tau_ep_Tc_ps": 2.47,
        }
    }


def _catalog():
    return SimpleNamespace(
        metadata={
            "Tc_K": 8.65,
            "D_m2_s": 1.58e-4,
            "sigma_n_S_m": 4.2e5,
            "delta0_J": 2.10e-22,
            "thickness_m": 7.0e-9,
            "width_m": 120.0e-9,
        }
    )


def test_build_gtdgl_material_uses_effective_diffusion_factor():
    material = build_gtdgl_material(_config(), _catalog(), diffusion_factor=0.68)

    assert material.D_base_m2_s == pytest.approx(1.58e-4)
    assert material.D_effective_factor == pytest.approx(0.68)
    assert material.D_m2_s == pytest.approx(1.58e-4 * 0.68)
    assert np.all(material.xi_mod_squared_m2(0.9) > 0.0)


def test_build_gtdgl_material_rejects_invalid_diffusion_factor():
    with pytest.raises(ValueError, match="diffusion_factor"):
        build_gtdgl_material(_config(), _catalog(), diffusion_factor=0.0)
