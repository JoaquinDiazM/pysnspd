"""Small independent invariants/accuracy checks; no production-sized network."""
from pathlib import Path
import sys
from types import SimpleNamespace
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT), str(Path(__file__).resolve().parent)]
from limited_ssp import step, integrate, _ratio, _electron_stoichiometry
from pysnspd.experimental.cell_validation import ElectronicCell
from pysnspd.experimental.kinetic_events import PhononGrid
from pysnspd.experimental.coupled_cells import CoupledCellSystem


class NormalCatalogue:
    count_nodes = np.array([.2, .5, 1., 2.])
    count_weights = np.array([.15, .3, .5, .6])
    eta = 1e-8

    def energy_kernel(self, delta, gamma):
        assert delta == gamma == 0
        return self.count_nodes.copy(), np.zeros(4), np.zeros(4)

    def evaluate(self, delta, gamma, population):
        return float(4*np.dot(self.count_weights*self.count_nodes, population)), 0., 0.


class FrozenMobility:
    scales = SimpleNamespace(bath_temperature_bar=.2)

    def amplitude_response(self, amplitude, force, temperature):
        assert force == 0
        return SimpleNamespace(velocity=0., heat=0.)


def make_system(cells=2, reaction_rate=5., face=.4, power=0., tau=np.inf):
    omega = np.array([.1, .3, .7, 1.5, 3., 4.])
    grid = PhononGrid(omega, np.full(6, .4),
        lambda e: np.zeros_like(e) if reaction_rate == 0 else np.ones_like(e),
        (.1, 4.), "four-state invariant oracle")
    return CoupledCellSystem(NormalCatalogue(), grid, FrozenMobility(),
        gammas=(0.,)*cells, tau_kin=tau, tau_escape=np.inf,
        diffusion_over_length_squared=face, external_powers=(power,)*cells,
        rate_prefactor=1. if reaction_rate == 0 else reaction_rate, reaction_quadrature_order=2,
        reaction_method="projected", reaction_layout="resolved_panels")


def energy(system, state):
    _, p, n = system.unpack(state)
    return float(np.sum(p*(4*system.catalog.count_weights*system.catalog.count_nodes))
                 +np.sum(n*(system.phonons.capacities*system.phonons.energies)))


def test_coupled_competing_reaction_face_budget_is_physical_and_conservative():
    system = make_system(reaction_rate=200., face=200.)
    initial = system.pack([0., 0.], [[.1, .4, .2, 1e-44], [0., 0., 0., 0.]], np.zeros((2, 6)))
    stats = {}
    obtained = step(system, 0., initial, 1., stats)
    _, p, n = system.unpack(obtained)
    assert stats["reaction_limited_evaluations"] > 0
    assert stats["transport_limited_evaluations"] > 0
    assert np.all((p >= 0) & (p <= 1)) and np.all(n >= 0)
    assert energy(system, obtained) == pytest.approx(energy(system, initial), abs=3e-14)
    assert stats["integrated_energy_weighted_flux_defect"] > 0


def test_empty_full_and_phonon_vacuum_faces_survive_ssp_composition():
    system = make_system(reaction_rate=50., face=30.)
    # Keep total energy in the positive-temperature FD range for the regular
    # closure, while including exact individual Pauli boundaries.
    initial = system.pack([0., 0.], [[1., 0., .01, 0.], [0., .2, 0., 0.]], np.zeros((2, 6)))
    _, states, stats = integrate(system, initial, .2, 4)
    for state in states:
        system.unpack(state)
        assert energy(system, state) == pytest.approx(energy(system, initial), abs=5e-14)
    assert stats["status"] == "COMPLETED_REQUIRES_CONVERGENCE_ASSESSMENT"


def test_source_energy_and_ledger_use_same_unlimited_heat():
    system = make_system(cells=1, reaction_rate=0., face=0., power=.03)
    cell = ElectronicCell(system.catalog, 0., 0.)
    initial = system.pack([0.], [cell.fermi_dirac(.2)], np.zeros((1, 6)))
    _, states, stats = integrate(system, initial, .1, 4)
    expected = .03*.1
    assert energy(system, states[-1])-energy(system, initial) == pytest.approx(expected, abs=2e-14)
    assert system.ledgers(states[-1])[0][0, 1] == pytest.approx(expected, abs=2e-15)
    assert stats["integrated_external_input"] == pytest.approx(expected, abs=2e-15)


def test_bgk_has_third_order_when_limiter_is_inactive():
    system = make_system(cells=1, reaction_rate=0., face=0., tau=.7)
    p0 = np.array([.3, .05, .015, .001])
    cell = ElectronicCell(system.catalog, 0., 0.)
    equilibrium = cell.fermi_dirac(cell.equivalent_temperature(p0))
    reference = equilibrium+(p0-equilibrium)*np.exp(-.2/.7)
    errors = []
    for steps in (5, 10, 20):
        initial = system.pack([0.], [p0], np.zeros((1, 6)))
        _, states, stats = integrate(system, initial, .2, steps)
        actual = system.unpack(states[-1])[1][0]
        errors.append(float(np.max(abs(actual-reference))))
        assert stats["limited_event_evaluations"] == 0
    assert 7 < errors[0]/errors[1] < 9
    assert 7 < errors[1]/errors[2] < 9


def test_baseline_infeasibility_is_rejected_without_repair():
    system = make_system(cells=1, reaction_rate=0., face=0., tau=.01)
    initial = system.pack([0.], [[.2, .03, .005, .001]], np.zeros((1, 6)))
    unchanged = initial.copy()
    with pytest.raises(ValueError, match="convex interval"):
        step(system, 0., initial, .02, {})
    assert np.array_equal(initial, unchanged)


def test_ratios_do_not_limit_available_interior_or_create_boundary_population():
    actual = _ratio(np.array([0., .2, 2.]), np.array([1., 1., 1.]), .5)
    assert actual == pytest.approx([0., .36, 1.], abs=1e-16)


def test_safety_ramp_is_continuous_at_raw_one_and_at_inactive_transition():
    eps = 1e-9
    raw = np.array([1-eps, 1., 1+eps, 1/.9-eps, 1/.9, 1/.9+eps])
    actual = _ratio(raw, np.ones(len(raw)), 1.)
    assert actual[:3] == pytest.approx(.9*raw[:3], abs=2e-16)
    assert abs(actual[2]-actual[0]) <= 2*eps
    assert actual[3] == pytest.approx(1-.9*eps, abs=2e-16)
    assert actual[4] == pytest.approx(1., abs=2e-16)
    assert actual[5] == 1.


def test_unlimited_map_equals_complete_rhs_with_moving_condensate():
    class AffineCatalogue(NormalCatalogue):
        # Analytic common-energy oracle with a nonzero amplitude derivative.
        # Native-pair reactions avoid claiming this synthetic dispersion is BCS.
        def energy_kernel(self, delta, gamma):
            return self.count_nodes+.1*delta, np.full(4, .1), np.zeros(4)

        def evaluate(self, delta, gamma, population):
            e = self.count_nodes+.1*delta
            return (.5*delta**2+float(4*np.dot(self.count_weights*e, population)),
                    delta+float(.4*np.dot(self.count_weights, population)), 0.)

    class MovingMobility(FrozenMobility):
        def amplitude_response(self, amplitude, force, temperature):
            return SimpleNamespace(velocity=-.03*force, heat=.03*force**2)

    system = make_system(cells=2, reaction_rate=2., face=.4, power=.03, tau=.7)
    system.catalog = AffineCatalogue()
    system.mobility = MovingMobility()
    system.tau_escape = 1.5
    system.reaction_method = "native_pairs"
    initial = system.pack([.03, .04], [[.2, .1, .02, .002], [.25, .08, .01, .001]],
        np.broadcast_to(system.bath_phonons+.01, (2, 6)))
    rhs = system.rhs(0., initial)
    stats = {}
    h = 1e-6
    obtained = step(system, 0., initial, h, stats)
    assert stats["limited_event_evaluations"] == 0
    assert np.max(abs(obtained-initial-h*rhs)) < 5e-16
    assert obtained[0] != initial[0]
    assert rhs[-1] != 0  # face transport is active
    ledger_rhs = rhs[system.population_size:-1].reshape(2, 4)
    assert np.all(ledger_rhs[:, 0] > 0)  # escape
    assert np.all(ledger_rhs[:, 1] > 0)  # external heat
    assert np.all(ledger_rhs[:, 2] != 0)  # reactions
    assert np.all(ledger_rhs[:, 3] > 0)  # condensate dissipation


def test_repeated_native_stoichiometric_slots_are_merged():
    event = SimpleNamespace(index_i=np.array([0, 0]), index_j=np.array([0, 1]),
        recombination=np.array([True, False]), omega=np.array([.4, .3]))
    _, coefficients = _electron_stoichiometry(event)
    assert np.array_equal(coefficients[0], [-2., 1.])
    assert np.array_equal(coefficients[1], [0., -1.])


def test_callback_is_small_metadata_at_first_tenth_last_steps():
    system = make_system(cells=1, reaction_rate=0., face=0.)
    initial = system.pack([0.], [[.2, .03, .005, .001]], np.zeros((1, 6)))
    messages = []
    integrate(system, initial, .1, 11, callback=messages.append)
    assert [m["macrostep"] for m in messages] == [1, 10, 11]
    assert all(not isinstance(value, np.ndarray) for message in messages for value in message.values())
