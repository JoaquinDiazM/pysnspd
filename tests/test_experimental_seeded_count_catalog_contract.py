"""The rejected polynomial is a Newton seed, never a returned energy/force.

These tiny causal inversions exercise the wrapper contract independently of
the saved benchmark. They are not another material or field-domain admission.
"""
from types import SimpleNamespace

import numpy as np
import pytest

from pysnspd.experimental.refined_cells import ComplementaryCountCatalog
from pysnspd.experimental.seeded_count_catalog import SeededCountCatalog, invert_count_seeded
from pysnspd.experimental.energy_catalog import energy_at_count_batch, retarded_spectrum


def source_and_predictor():
    # This analytic vacuum is only an algebraic oracle for moment derivatives.
    def vacuum(a, g):
        if not .9 <= a <= 1.1 or not 0 <= g <= .1:
            raise ValueError('test vacuum support')
        return .5*a*a+.7*g+.1*g*g, a, .7+.2*g
    base = SimpleNamespace(eta=1e-8, vacuum=SimpleNamespace(evaluate=vacuum))
    source = ComplementaryCountCatalog(base, np.array([1e-7, 1e-4, .03, .2, 1., 3.]),
                                       np.array([2e-5, .001, .03, .15, .5, 1.]))
    coeff = np.zeros((3, 3, len(source.count_nodes)))
    coeff[0, 0] = [6., 5., 4., 3., 2., 1.]  # Invalid as an energy: reversed.

    def forbidden(*args):
        raise AssertionError('predictor energy or force must not enter production values')
    predictor = SimpleNamespace(count_nodes=source.count_nodes,
        count_weights=source.count_weights, eta=source.eta, coefficients=coeff,
        amplitude_bounds=(.985, 1.005), gamma_bounds=(.003, .009),
        energy_kernel=forbidden, evaluate=forbidden,
        _da=np.full_like(coeff, 1e20), _dg=np.full_like(coeff, -1e20))
    return source, predictor


def test_wrapper_corrects_invalid_predictor_and_uses_original_implicit_derivatives():
    source, predictor = source_and_predictor()
    model = SeededCountCatalog(source, predictor)
    a, g = .997, .0055
    values, diagnostics = model.energy_kernel_with_diagnostics(a, g)
    expected = source.energy_kernel(a, g)
    np.testing.assert_allclose(values, expected, rtol=2e-9, atol=2e-12)
    assert diagnostics['seed_bracket_replacements'] > 0
    assert diagnostics['spectrum_points'] >= len(source.count_nodes)
    assert np.all(values[0] > 0) and np.all(np.diff(values[0]) > 0)
    # Re-evaluate the original causal equations at the actual returned E.
    c, s = retarded_spectrum(values[0], delta=a, gamma=g, eta=source.eta)
    w = 1j*a/s
    count = np.real(w-1j*g*a*a/(2*w*w))
    np.testing.assert_allclose(count, source.count_nodes, rtol=1e-8, atol=2e-12)
    np.testing.assert_array_equal(values[1], s.imag/c.real)
    np.testing.assert_array_equal(values[2], -s.real*(s.imag/c.real))
    assert model.count_nodes is source.count_nodes
    assert model.count_weights is source.count_weights
    assert model.vacuum is source.vacuum


def test_fixed_population_moment_forces_differentiate_corrected_energy():
    source, predictor = source_and_predictor()
    model = SeededCountCatalog(source, predictor)
    a, g = .997, .0055
    p = .13*np.exp(-source.count_nodes/.7)
    original = p.copy()
    moments = model.evaluate(a, g, p)
    h = 1e-5
    def derivative(da, dg):
        def energy(t):
            return model.evaluate(a+t*da, g+t*dg, p)[0]
        return (energy(-2*h)-8*energy(-h)+8*energy(h)-energy(2*h))/(12*h)
    assert moments[1] == pytest.approx(derivative(1., 0.), rel=2e-7, abs=2e-9)
    assert moments[2] == pytest.approx(derivative(0., 1.), rel=2e-7, abs=2e-9)
    np.testing.assert_array_equal(p, original)


@pytest.mark.parametrize('field', [( .984, .006), (1.006, .006), (.997, .0029),
                                  (.997, .0091), (.997, 0.), (np.nan, .006)])
def test_wrapper_never_extrapolates_predictor_box(field):
    source, predictor = source_and_predictor()
    with pytest.raises(ValueError, match='predictor box'):
        SeededCountCatalog(source, predictor).energy_kernel(*field)


@pytest.mark.parametrize('field', ['count_nodes', 'count_weights', 'eta'])
def test_wrapper_rejects_incompatible_count_measure_or_regulator(field):
    source, predictor = source_and_predictor()
    setattr(predictor, field, getattr(predictor, field)*1.01)
    with pytest.raises(ValueError):
        SeededCountCatalog(source, predictor)


@pytest.mark.parametrize('invalid', [-.01, 1.01, np.nan, np.inf, 1j])
def test_wrapper_rejects_nonphysical_populations_without_solving(invalid, monkeypatch):
    source, predictor = source_and_predictor()
    model = SeededCountCatalog(source, predictor)
    def forbidden(*args):
        raise AssertionError('invalid populations must be rejected before the inverse')
    monkeypatch.setattr(SeededCountCatalog, 'energy_kernel', forbidden)
    p = np.full(len(source.count_nodes), invalid)
    with pytest.raises(ValueError):
        model.evaluate(.997, .0055, p)


def test_nonfinite_trial_is_replaced_by_original_upper_bracket_not_returned():
    x = np.array([1e-7, .03, .4, 2.])
    guess = np.array([np.nan, np.inf, -np.inf, -1.])
    actual, diagnostics = invert_count_seeded(x, .997, .0055, 1e-8, guess)
    expected = energy_at_count_batch(x, delta=.997, gamma=.0055, eta=1e-8)
    np.testing.assert_array_equal(actual, expected)
    assert diagnostics['seed_bracket_replacements'] == len(x)
