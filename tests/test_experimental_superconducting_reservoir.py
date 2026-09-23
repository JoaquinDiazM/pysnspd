"""Cheap thermal and branch tests; synthetic spectra are algebraic oracles only."""
import numpy as np
import pytest
from scipy.special import expit, xlogy

from pysnspd.experimental.energy_catalog import build_vacuum_catalog, HBAR_J_S, E_CHARGE_C
from pysnspd.experimental.superconducting_reservoir import (
    SuperconductingReservoir, ReservoirBranchError)


class AlgebraCatalog:
    """Analytic BCS vacuum plus a deliberately simple, nonphysical test spectrum."""
    def __init__(self):
        self.vacuum = build_vacuum_catalog(np.array([.08, .5, 1., 1.5]),
            np.array([0., .2, .6, 1.2]), Tc_K=8., N0_per_J_m3=1e47,
            D_m2_s=5e-5, analytic=True)
        self.count_nodes = np.array([.03, .2, .7, 1.3, 2.5, 6.])
        self.count_weights = np.array([.1, .2, .4, .6, 1., 1.7])

    def energy_kernel(self, a, gamma):
        self.vacuum.evaluate(a, gamma)
        e0 = np.sqrt(self.count_nodes**2+a*a)
        return e0+.03*gamma, a/e0, np.full_like(e0, .03)


@pytest.fixture
def catalog():
    return AlgebraCatalog()


def reservoir(catalog, theta=.12, **kwargs):
    return SuperconductingReservoir(catalog, theta, cross_section_m2=120e-9*7e-9,
                                    resolved_length_m=360e-9, **kwargs)


def independent_free_energy(catalog, model, a, q):
    # Closed vacuum expression, separate from the reservoir implementation.
    m = a*a/(a*a+.01)
    gamma = (m*q)**2/model.gap_ratio
    assert gamma < a
    vacuum = a*a*(np.log(a)-.5)+np.pi*a*gamma/2-gamma*gamma/3
    energies = np.sqrt(catalog.count_nodes**2+a*a)+.03*gamma
    f = vacuum-4*model.bath_theta*np.dot(catalog.count_weights,
                                        np.logaddexp(0., -energies/model.bath_theta))
    return f+np.pi/4*a*a*q*q*(1-m*m)


def fd5(function, step=2e-5):
    return (function(-2*step)-8*function(-step)+8*function(step)-function(2*step))/(12*step)


def test_thermal_free_energy_and_gradients_match_independent_grand_potential(catalog):
    model = reservoir(catalog, theta=.35)
    a, q = .87, .09
    point = model.evaluate(a, q)
    expected = independent_free_energy(catalog, model, a, q)
    assert point.free_energy_bar == pytest.approx(expected, abs=3e-15)
    da = fd5(lambda h: independent_free_energy(catalog, model, a+h, q))
    dq = fd5(lambda h: independent_free_energy(catalog, model, a, q+h))
    assert point.amplitude_force_bar == pytest.approx(da, abs=2e-10)
    assert point.current_bar == pytest.approx(dq, abs=2e-10)
    energies, ea, eg = catalog.energy_kernel(a, point.gamma_bar)
    p = expit(-energies/model.bath_theta)
    np.testing.assert_array_equal(point.occupation, p)
    u0, ua0, ug0 = catalog.vacuum.evaluate(a, point.gamma_bar)
    expected_moments = (u0+4*np.dot(catalog.count_weights*energies, p),
                        ua0+4*np.dot(catalog.count_weights*ea, p),
                        ug0+4*np.dot(catalog.count_weights*eg, p))
    np.testing.assert_allclose(point.electronic_moments_bar, expected_moments, atol=2e-15)
    entropy = -4*np.dot(catalog.count_weights, xlogy(p, p)+xlogy(1-p, 1-p))
    m = a*a/(a*a+.01)
    from_entropy = expected_moments[0]-model.bath_theta*entropy+model.kappa*a*a*q*q*(1-m*m)
    assert point.free_energy_bar == pytest.approx(from_entropy, abs=2e-15)


def test_zero_temperature_branch_is_BCS_and_skips_excitation_kernel(catalog, monkeypatch):
    model = reservoir(catalog, theta=0.)
    def forbidden(*args):
        raise AssertionError("the vacuum branch requires no excitation solve")
    monkeypatch.setattr(catalog, "energy_kernel", forbidden)
    branch = model.solve_at_q(0.)
    assert branch.amplitude_bar == pytest.approx(1., abs=2e-11)
    assert branch.free_energy_bar == pytest.approx(-.5, abs=2e-14)
    assert abs(branch.amplitude_force_bar) < 5e-11
    assert branch.current_bar == 0.
    assert branch.amplitude_curvature == pytest.approx(2., abs=2e-8)
    m = 1/(1+.01)
    expected_slope = np.pi*m*m/model.gap_ratio+2*model.kappa*(1-m*m)
    assert branch.differential_current_bar == pytest.approx(expected_slope, abs=2e-8)
    expected_inductance = (HBAR_J_S/(2*E_CHARGE_C))*model.length_bar/(model.current_scale_A*expected_slope)
    assert branch.L_res_diff_H == pytest.approx(expected_inductance, rel=2e-8)
    assert branch.checked_qs == (0.,)
    assert branch.stable


def test_reduced_current_slope_includes_the_amplitude_response(catalog):
    model = reservoir(catalog, theta=0.)
    q = .1
    branch = model.solve_at_q(q, continuation_qs=(.03, .07))
    h = 3e-4
    plus = model.solve_at_q(q+h)
    minus = model.solve_at_q(q-h)
    measured = (plus.current_bar-minus.current_bar)/(2*h)
    assert branch.differential_current_bar == pytest.approx(measured, rel=2e-7)
    assert branch.differential_current_bar < branch.thermal_hessian[1, 1]
    assert branch.amplitude_curvature > branch.amplitude_curvature_uncertainty
    assert branch.differential_current_bar > branch.differential_current_uncertainty
    assert branch.maxwell_disagreement < 2e-7
    assert branch.checked_qs == (0., .03, .07, .1)


def test_explicit_current_inverse_and_fixed_exterior_partition(catalog):
    model = reservoir(catalog, theta=0.)
    expected = model.solve_at_q(.075)
    inverse = model.solve_at_current(expected.current_A, q_bracket=(0., .12))
    assert inverse.q_bare_bar == pytest.approx(.075, abs=3e-10)
    assert inverse.amplitude_bar == pytest.approx(expected.amplitude_bar, abs=3e-10)
    partition = model.inductance_partition(inverse)
    assert partition.total_reference_H == 10e-9
    assert partition.exterior_fixed_H > 0
    assert partition.exterior_fixed_H+partition.resolved_differential_H == pytest.approx(10e-9)
    assert partition.reference_current_A == inverse.current_A
    saved_exterior = partition.exterior_fixed_H
    model.solve_at_q(.1)  # A later field does not mutate the reference partition.
    assert partition.exterior_fixed_H == saved_exterior
    with pytest.raises(ReservoirBranchError, match="exterior"):
        model.inductance_partition(inverse, total_reference_H=inverse.L_res_diff_H/2)
    with pytest.raises(ReservoirBranchError, match="outside"):
        model.solve_at_current(10*expected.current_A, q_bracket=(0., .12))


@pytest.mark.parametrize("kind", ["amplitude", "current"])
def test_unstable_roots_are_rejected_instead_of_selected(catalog, kind):
    source = catalog.vacuum
    class VacuumOracle:
        def __getattr__(self, key):
            return getattr(source, key)
        def evaluate(self, a, gamma):
            source.evaluate(a, gamma)  # Real support guards remain operative.
            if kind == "amplitude":
                return -(a-1)**2+gamma, -2*(a-1), 1.
            return (a-1)**2-gamma, 2*(a-1), -1.
    catalog.vacuum = VacuumOracle()
    model = reservoir(catalog, theta=0.)
    with pytest.raises(ReservoirBranchError, match=("curvature" if kind == "amplitude" else "differential")):
        model.solve_at_q(0.)


def test_brackets_paths_temperature_and_support_are_explicit(catalog):
    model = reservoir(catalog, theta=0.)
    for theta in (-.1, np.nan, np.inf):
        with pytest.raises(ValueError):
            reservoir(catalog, theta=theta)
    for bracket in ((.5, .7),):
        with pytest.raises(ReservoirBranchError, match="bracket"):
            model.solve_at_q(0., amplitude_bracket=bracket)
    for bracket in ((1., .8), (.01, 1.02), (0., 1.02)):
        with pytest.raises(ValueError, match="bracket"):
            model.solve_at_q(0., amplitude_bracket=bracket)
    for path in ((.07, .03), (.2,), (-.02,), (.1,)):
        with pytest.raises(ValueError, match="continuation"):
            model.solve_at_q(.1, continuation_qs=path)
    with pytest.raises(ValueError, match="Gamma"):
        model.evaluate(.9, 4.)
    opposite = model.evaluate(.9, -.08)
    forward = model.evaluate(.9, .08)
    assert opposite.free_energy_bar == forward.free_energy_bar
    assert opposite.amplitude_force_bar == forward.amplitude_force_bar
    assert opposite.current_bar == -forward.current_bar


def test_warm_root_is_stationary_for_the_actual_thermal_free_energy(catalog):
    model = reservoir(catalog, theta=.12)
    branch = model.solve_at_q(.06)
    measured = fd5(lambda h: independent_free_energy(catalog, model, branch.amplitude_bar+h, .06))
    assert abs(measured) < 3e-9
    assert abs(branch.amplitude_force_bar) < 3e-9
    assert branch.amplitude_bar < 1.
    assert branch.stable
    assert np.all((branch.occupation > 0) & (branch.occupation < .5))
