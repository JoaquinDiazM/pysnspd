"""Cell contracts use the admitted R2 artefact, not an analytic mock DOS."""
from pathlib import Path

import numpy as np
import pytest

from pysnspd.experimental.energy_catalog import OccupationEnergyCatalog
from pysnspd.experimental.cell_validation import ElectronicCell, EnergyFacePair, linear_count_moments, rk4_trajectory


@pytest.fixture(scope="module")
def catalog():
    root = Path(__file__).resolve().parents[1]
    return OccupationEnergyCatalog.load(
        root / "docs/implementation/stage1_r2/catalogs/occupation_catalog.npz")


@pytest.mark.parametrize("amplitude,gamma", [(1.0, 0.0), (.72, .2), (.35, .65)])
def test_bgk_conserves_catalogue_energy_and_increases_entropy(catalog, amplitude, gamma):
    cell = ElectronicCell(catalog, amplitude, gamma)
    x = catalog.count_nodes
    p = .12*np.exp(-x/.25) + .064*np.exp(-((x-.95)/.3)**2)
    rhs, temperature = cell.bgk(p, .7)
    assert abs(4*np.dot(cell.weights*cell.energies, rhs)) < 2e-12
    target = cell.fermi_dirac(temperature)
    relaxed = target + (p-target)*np.exp(-2/.7)
    assert abs(cell.excitation_energy(relaxed)-cell.excitation_energy(p)) < 2e-12
    assert cell.entropy(relaxed) > cell.entropy(p)
    assert np.max(abs(cell.bgk(target, .7)[0])) < 2e-13


def test_rk4_bgk_matches_independent_exact_solution(catalog):
    cell = ElectronicCell(catalog, .72, .2)
    p = .2*np.exp(-catalog.count_nodes/.25)
    target = cell.fermi_dirac(cell.equivalent_temperature(p))
    exact = target + (p-target)*np.exp(-3)
    errors = []
    for steps in (30, 60, 120):
        _, trajectory = rk4_trajectory(lambda time, state: cell.bgk(state, 1)[0], p, 3, steps)
        assert np.min(trajectory) >= 0 and np.max(trajectory) <= 1
        errors.append(np.max(abs(trajectory[-1]-exact)))
    assert errors[1] < errors[0]/12
    assert errors[2] < errors[1]/12
    assert errors[-1] < 1e-8


def test_source_defined_at_vacuum_and_has_exact_energy_moment(catalog):
    cell = ElectronicCell(catalog, .72, .2)
    p = np.zeros_like(cell.energies)
    source = cell.heating(p, .05, .08)
    assert np.min(source) >= 0
    assert abs(4*np.dot(cell.weights*cell.energies, source)-.05) < 1e-14
    p[0] = 1
    assert cell.heating(p, .05, .08)[0] == 0
    with pytest.raises(ValueError, match="saturated"):
        cell.heating(np.ones_like(p), .05, .08)


def test_source_has_documented_thermal_direction(catalog):
    cell = ElectronicCell(catalog, .72, .2)
    p = cell.fermi_dirac(.25)
    expected = cell.energies*p*(1-p)
    expected *= .05/(4*np.dot(cell.weights*cell.energies, expected))
    assert np.allclose(cell.heating(p, .05, .08), expected, rtol=1e-11, atol=1e-15)


def test_positive_temperature_support_is_not_silently_extended(catalog):
    cell = ElectronicCell(catalog, .72, .2)
    assert cell.equivalent_temperature(np.zeros_like(cell.energies)) == 0
    with pytest.raises(ValueError, match="positive-temperature"):
        cell.equivalent_temperature(np.full_like(cell.energies, .8))


@pytest.mark.parametrize("value", [-1e-15, 1+1e-15, np.nan, np.inf])
def test_no_pauli_or_nonfinite_repair(catalog, value):
    cell = ElectronicCell(catalog, .72, .2)
    p = np.zeros_like(cell.energies)
    p[0] = value
    with pytest.raises(ValueError):
        cell.bgk(p, 1)


def test_cell_uses_catalogue_support_guard(catalog):
    with pytest.raises(ValueError):
        ElectronicCell(catalog, .01, .2)
    with pytest.raises(ValueError):
        ElectronicCell(catalog, .72, 1.3)


def test_linear_count_moments_have_correct_number_and_energy():
    # x(E)=2E on [1,4], integrated on three unequal hat intervals.
    mass, energetic = linear_count_moments(np.array([2., 8.]), np.array([1., 4.]),
                                           np.array([1., 1.7, 2.3, 4.]))
    assert abs(np.sum(mass)-6.) < 2e-14
    assert abs(np.sum(energetic)-15.) < 3e-14
    assert np.min(mass) > 0 and np.min(energetic) > 0


@pytest.fixture
def pair(catalog):
    left, right = ElectronicCell(catalog, .43, .02), ElectronicCell(catalog, 1.2, .2)
    boundary = min(left.energies[-1], right.energies[-1])
    grid = np.r_[np.linspace(1e-6, boundary, 65), max(left.energies[-1], right.energies[-1])]
    return EnergyFacePair(left, right, grid)


def test_unequal_spectra_exchange_conserves_energy_and_pauli(pair):
    rng = np.random.default_rng(2131)
    p = rng.uniform(0, 1, pair.capacities.shape)
    p[0, 0] = 0; p[1, -1] = 1
    derivative = pair.rhs(p)
    assert abs(4*np.sum(pair.energy_moments*derivative)) < 2e-14
    updated = pair.exact(p, 12.)
    assert np.min(updated) >= 0 and np.max(updated) <= 1
    assert abs(np.sum(pair.reconstructed_energy(updated)-pair.reconstructed_energy(p))) < 5e-13
    assert np.array_equal(updated[:, ~pair.active], p[:, ~pair.active])


def test_same_energy_distribution_is_stationary_for_different_gaps(pair):
    f = 1/(1+np.exp(pair.shared_energy/.3))
    state = np.tile(f, (2, 1))
    assert np.array_equal(pair.rhs(state), np.zeros_like(state))
    native = pair.native_populations(state)
    assert np.max(abs(native[0]-native[1])) > .01


def test_transport_refuses_native_population_extrapolation(catalog):
    cell = ElectronicCell(catalog, .72, .2)
    with pytest.raises(ValueError, match="extrapolation"):
        EnergyFacePair(cell, cell, np.linspace(.5, 2, 20))


def test_exact_exchange_does_not_cancel_a_cold_population_at_zero_time(pair):
    p = np.tile(np.linspace(.01, .2, len(pair.shared_energy)), (2, 1))
    p[1] = 1e-44
    assert np.array_equal(pair.exact(p, 0.), p)
    assert np.min(pair.exact(p, 1e-14)) >= 0


def test_private_tail_cannot_exchange_nonthermal_population(pair):
    p = np.zeros_like(pair.capacities)
    outside = pair.shared_energy > pair.common_upper_energy
    p[0, outside] = 1
    assert np.array_equal(pair.conductance[outside], np.zeros(np.count_nonzero(outside)))
    assert np.array_equal(pair.rhs(p), np.zeros_like(p))
    assert np.array_equal(pair.exact(p, 10), p)
    boundary = pair.shared_energy == pair.common_upper_energy
    assert np.all(pair.conductance[boundary] > 0)


@pytest.mark.parametrize("mobility", [0, -1, np.nan, np.inf])
def test_condensate_mobility_must_be_positive(catalog, mobility):
    cell = ElectronicCell(catalog, .72, 0)
    with pytest.raises(ValueError, match="mobility"):
        cell.condensate_relaxation(np.zeros_like(cell.energies), mobility)
