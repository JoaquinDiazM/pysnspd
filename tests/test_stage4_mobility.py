"""Stage 4 KWT scenario changes, checked against the dimensional equation."""
from dataclasses import FrozenInstanceError, replace
import json

import numpy as np
import pytest

from pysnspd.experimental.cell_closures import CellScales, KWTMobility
from pysnspd.experimental.spatial_dynamics import kwt_spatial_response


HBAR = 1.054571817e-34
KB = 1.380649e-23
SCALES = CellScales(1.764 * KB * 8.65, 1.1e47, 8.65, .9, 1.)


def dimensional_response(z, force, temperature_bar, scales, tau_ee_ps, tau_ep_ps):
    """Invert the dimensional temporal matrix instead of its eigenformula."""
    delta = scales.delta0_J * np.asarray(z)
    temperature = max(temperature_bar * scales.delta0_J / KB, scales.Tb_K)
    ratio = temperature / scales.Tc_K
    tau_s = 1e-12 / (ratio / tau_ee_ps + ratio**3 / tau_ep_ps)
    tau_gl_s = np.pi * HBAR / (8 * KB * scales.Tc_K)
    prefactor = scales.N0_per_J_m3 * np.sqrt((1 + ratio) / 2)
    temporal = prefactor * tau_gl_s * (
        np.eye(2) + 4 * tau_s**2 / HBAR**2 * np.outer(delta, delta)
    ) / np.sqrt(1 + 4 * tau_s**2 / HBAR**2 * np.dot(delta, delta))
    force_SI = scales.N0_per_J_m3 * scales.delta0_J * np.asarray(force)
    velocity_SI = np.linalg.solve(temporal, -force_SI / 2)
    heat_SI = 2 * np.dot(velocity_SI, temporal @ velocity_SI)
    return (
        velocity_SI * scales.t_ref_s / scales.delta0_J,
        heat_SI * scales.t_ref_s / scales.energy_density_scale_J_m3,
    )


@pytest.mark.parametrize("times", [(.50, 2.47), (5., 24.7), (6., 24.7)])
@pytest.mark.parametrize("amplitude", [0., .05, .8])
@pytest.mark.parametrize("temperature_bar", [0., .4])
def test_actual_relaxation_times_match_dimensional_equation(times, amplitude, temperature_bar):
    model = KWTMobility(SCALES, tau_ee_Tc_ps=times[0], tau_ep_Tc_ps=times[1])
    z = amplitude * np.array([.6, .8])
    force = np.array([.7, -.45])
    response = model.tensor_response(z, force, temperature_bar)
    velocity, heat = dimensional_response(z, force, temperature_bar, SCALES, *times)
    # The independent matrix solve loses a few digits at large KWT anisotropy.
    np.testing.assert_allclose(response.velocity, velocity, rtol=2e-10, atol=1e-13)
    assert response.heat == pytest.approx(heat, rel=2e-10, abs=1e-13)
    assert response.heat == pytest.approx(-force @ response.velocity, rel=1e-12)
    assert response.heat > 0
    assert np.linalg.eigvalsh(response.mobility).min() > 0


def test_changed_times_are_not_a_global_clock_multiplier():
    inherited = KWTMobility(SCALES)
    changed = KWTMobility(SCALES, tau_ee_Tc_ps=6., tau_ep_Tc_ps=24.7)
    theta = SCALES.temperature_to_bar(.9)
    old = inherited.tensor_response([.8, 0.], [1., 1.], theta)
    new = changed.tensor_response([.8, 0.], [1., 1.], theta)
    assert new.mobility[0, 0] < old.mobility[0, 0]
    assert new.mobility[1, 1] > old.mobility[1, 1]
    radial = changed.amplitude_response(.8, 1., theta)
    assert radial.velocity == pytest.approx(new.velocity[0], rel=1e-13)
    scaled_clock = replace(changed, scales=replace(SCALES, t_ref_ps=3.))
    scaled = scaled_clock.tensor_response([.8, 0.], [1., 1.], theta)
    np.testing.assert_allclose(scaled.velocity / 3, new.velocity, rtol=1e-13)
    assert scaled.taupsi_ps == new.taupsi_ps


def test_changed_mobility_reaches_spatial_assembly_with_power_balance():
    model = KWTMobility(SCALES, tau_ee_Tc_ps=6., tau_ep_Tc_ps=24.7)
    z = np.array([.3 + .4j])
    gradient = np.array([[.7, -.2]])
    load = np.array([[.01, .02]])
    result = kwt_spatial_response(
        z, gradient, [2.], model, [.2], [1e-5], boundary_load_bar=load,
    )
    direct = model.tensor_response([z[0].real, z[0].imag], (gradient[0] - load[0]) / 2, .2)
    np.testing.assert_allclose(
        [result.material_velocity_bar[0].real, result.material_velocity_bar[0].imag],
        direct.velocity, rtol=1e-12, atol=1e-13,
    )
    assert result.condensate_heat_rate_bar == pytest.approx(2 * direct.heat, rel=1e-12)
    assert result.identity_residual_bar == pytest.approx(0., abs=1e-12)


def test_parameters_are_immutable_and_metadata_records_actual_values():
    default = KWTMobility(SCALES)
    assert (default.tau_ee_Tc_ps, default.tau_ep_Tc_ps) == (.50, 2.47)
    model = replace(default, tau_ee_Tc_ps=6., tau_ep_Tc_ps=24.7)
    with pytest.raises(FrozenInstanceError):
        model.tau_ep_Tc_ps = 1.
    metadata = json.loads(json.dumps(model.metadata()))
    assert metadata["tau_ee_Tc_ps"] == 6.
    assert metadata["tau_ep_Tc_ps"] == 24.7
    assert metadata["scales"]["t_ref_ps"] == SCALES.t_ref_ps
    metadata["tau_ee_Tc_ps"] = 1.
    metadata["scales"]["t_ref_ps"] = 99.
    assert model.metadata()["tau_ee_Tc_ps"] == 6.
    assert model.scales.t_ref_ps == SCALES.t_ref_ps


@pytest.mark.parametrize("name", ["tau_ee_Tc_ps", "tau_ep_Tc_ps"])
@pytest.mark.parametrize("value", [0., -1., np.nan, np.inf, -np.inf])
def test_invalid_relaxation_times_fail_at_construction(name, value):
    with pytest.raises(ValueError, match=name):
        KWTMobility(SCALES, **{name: value})
