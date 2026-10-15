"""Calibration provenance, exact native scaling and callback recalibration."""
import functools
import json
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.constants import Boltzmann as kB, elementary_charge as e, hbar

from pysnspd.experimental.depairing_calibration import calibration_result


@functools.lru_cache(None)
def nominal(rs=425., width=80e-9):
    return calibration_result(sheet_resistance_ohm=rs, width_m=width,
                              requested_current_A=21.5e-6)


def test_nominal_425_and_608_reproduce_target_with_finite_quadrature():
    first, alternative = nominal(), nominal(608.)
    assert first["diffusion_m2_s"] == pytest.approx(1.4104e-4, rel=8e-5)
    assert alternative["diffusion_m2_s"] == pytest.approx(6.89137e-5, rel=8e-5)
    for value in (first, alternative):
        assert value["target_depairing_current_A"] == pytest.approx(38.8e-6*80/120)
        assert value["achieved_depairing_current_A"] == pytest.approx(
            value["target_depairing_current_A"], rel=2e-8)
        assert abs(value["relative_target_residual"]) < 2e-8
        assert value["quadrature"]["electronic_states"] == 54
        assert value["requested_current_below_depairing"]


def test_material_derived_from_fixed_sheet_conductance_and_new_diffusion():
    first, alternative = nominal(), nominal(608.)
    for result in (first, alternative):
        diffusion = result["diffusion_m2_s"]
        assert result["N0_per_J_m3"] == pytest.approx(1/(2*e**2*diffusion*
            result["sheet_resistance_ohm"]*result["thickness_m"]), rel=2e-15)
        assert result["ell0_m"] == pytest.approx(np.sqrt(hbar*diffusion/(2*kB*8.65)))
        assert result["N0_per_J_m3"]/result["reference"]["N0_per_J_m3"] == pytest.approx(
            result["reference"]["diffusion_m2_s"]/diffusion)
    assert first["diffusion_m2_s"]/alternative["diffusion_m2_s"] == pytest.approx((608/425)**2, rel=2e-8)


def test_width_scaled_target_keeps_same_calibrated_material():
    wide = nominal(425., 120e-9)
    assert wide["target_depairing_current_A"] == pytest.approx(38.8e-6)
    assert wide["diffusion_m2_s"] == pytest.approx(nominal()["diffusion_m2_s"], rel=2e-8)


def test_provenance_distinguishes_inferences_and_requires_new_closure_check():
    result = nominal()
    json.dumps(result, allow_nan=False)
    assert result["calibration"]["target_is_width_scaled_inference"]
    assert result["calibration"]["diffusion_is_fitted_effective_parameter"]
    assert result["calibration"]["requires_verification_with_new_closure"]
    assert not result["calibration"]["physical_admission"]
    assert "conditional" in result["calibration"]["sheet_resistance_status"]
    assert result["geometry"]["rebuild_physical_geometry_with_calibrated_ell0"]
    assert not result["geometry"]["reuse_old_normalized_mesh_or_ell0"]
    assert "PDF13" in result["sources"]["K20"]["depairing_locator"]
    assert "independent" in result["assumptions"]["sigma_nm"]


def test_actual_custom_closure_is_rebuilt_and_root_recalibrated():
    calls = []
    target = 20e-6
    def factory(**kwargs):
        calls.append(kwargs)
        # Deliberately break the native square-root scaling, so the utility must
        # solve this supplied closure rather than attach an unchecked D estimate.
        current = 30e-6*(5e-5/kwargs["diffusion_m2_s"])**.75
        return SimpleNamespace(bulk=lambda **ignored:
                               SimpleNamespace(depairing_current_A=current))
    result = calibration_result(target_depairing_current_A=target, bulk_factory=factory)
    assert result["diffusion_m2_s"] == pytest.approx(5e-5*(30/20)**(1/.75), rel=1e-8)
    assert result["achieved_depairing_current_A"] == pytest.approx(target, rel=2e-8)
    assert "root" in result["calibration"]["method"]
    assert not result["calibration"]["requires_verification_with_new_closure"]
    assert not result["calibration"]["target_is_width_scaled_inference"]
    assert len({id(item["graph"]) for item in calls}) == len(calls)
    assert all(item["graph"].n_nodes == 2 for item in calls)
    for item in calls:
        assert item["ell0_m"] == pytest.approx(np.sqrt(hbar*item["diffusion_m2_s"]/(2*kB*8.65)))


@pytest.mark.parametrize("name,value", [
    ("sheet_resistance_ohm", 0.), ("width_m", -1.), ("thickness_m", float("nan")),
    ("critical_temperature_K", float("inf")), ("bath_temperature_K", 8.65),
    ("reference_diffusion_m2_s", True), ("target_depairing_current_A", 1j),
    ("requested_current_A", -1.), ("requested_current_A", 38.8e-6*80/120),
    ("count_order", 1), ("count_cutoff_bar", -1.), ("bulk_factory", 42)])
def test_invalid_inputs_fail_before_calibration(name, value):
    with pytest.raises(ValueError):
        calibration_result(**{name: value})