"""Physical native-energy inverse, wide temperature jumps and domain limits."""
import numpy as np
import pytest
from scipy.optimize import brentq
from scipy.special import expit
from test_thermal_stiffness_coupled_candidate import small_system


def fd(functional, gap, temperature):
    en = functional.energies_bar(gap)
    return expit(-en/(np.asarray(temperature)/functional.Tc_K)[:, None])


@pytest.mark.parametrize("temperature", [
    np.array([.03, .9, 7.]),
    np.array([40., 100., 1000.]),
    np.array([.1, 8.65, 300.]),
])
def test_recovers_known_FD_temperature_without_bath_floor(temperature):
    system, state = small_system()
    f = system.functional
    # Zero gap is a normal-state limiting case, not a singular temperature.
    state.delta_bar[1] = 0.
    p = fd(f, state.delta_bar, temperature)
    actual = f.equivalent_temperature_K(state.delta_bar, p)
    np.testing.assert_allclose(actual, temperature, rtol=4e-13, atol=2e-14)
    en = f.energies_bar(state.delta_bar)
    np.testing.assert_allclose(
        np.sum(f.weights_bar*en*fd(f, state.delta_bar, actual), axis=1),
        np.sum(f.weights_bar*en*p, axis=1), rtol=4e-13, atol=0.)


def test_nonthermal_population_matches_independent_scalar_root():
    system, state = small_system()
    f = system.functional
    p = fd(f, state.delta_bar, [1.5, 6., 15.])
    p *= np.linspace(.4, 1.6, p.shape[1])[None, :]
    en = f.energies_bar(state.delta_bar)
    target = np.sum(f.weights_bar*en*p, axis=1)
    def root(row, target_value):
        return brentq(
            lambda T: np.sum(f.weights_bar*row*expit(-row*f.Tc_K/T))-target_value,
            .001, 10000., xtol=1e-14, rtol=4*np.finfo(float).eps)
    expected = np.array([root(row, value) for row, value in zip(en, target)])
    actual = f.equivalent_temperature_K(state.delta_bar, p)
    np.testing.assert_allclose(actual, expected, rtol=2e-13, atol=2e-14)


def test_cached_guess_does_not_remember_old_physical_temperature():
    system, state = small_system()
    f = system.functional
    for temperatures in ([.03, .9, 10.], [1000., 300., 800.], [.04, .1, .08],
                         [7., 5., 20.]):
        p = fd(f, state.delta_bar, temperatures)
        actual = f.equivalent_temperature_K(state.delta_bar, p)
        np.testing.assert_allclose(actual, temperatures, rtol=4e-13, atol=2e-14)


@pytest.mark.parametrize("occupation", [0., .5, .6, 1.])
def test_invalid_positive_FD_energy_domain_is_not_clipped(occupation):
    system, state = small_system()
    p = np.full_like(state.p, occupation)
    with pytest.raises(ValueError, match="temperature"):
        system.functional.equivalent_temperature_K(state.delta_bar, p)
