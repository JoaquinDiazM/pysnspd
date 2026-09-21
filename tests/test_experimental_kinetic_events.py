"""Physical event contracts; no production transient is activated."""
from pathlib import Path
import numpy as np
import pytest
from scipy.special import expit

from pysnspd.experimental.cell_validation import ElectronicCell
from pysnspd.experimental.energy_catalog import OccupationEnergyCatalog
from pysnspd.experimental.kinetic_events import ElectronPhononEvents, ProjectedElectronPhononEvents, HybridElectronPhononEvents, PhononGrid


def phonons(low=.01, high=4., nodes=65):
    energy = np.linspace(low, high, nodes)
    return PhononGrid(energy, .2*np.ones(nodes), lambda omega: .03*(omega/4)**2,
                      (low, high), "synthetic positive capacities and acoustic coupling")


@pytest.fixture(scope="module")
def catalog():
    return OccupationEnergyCatalog.load(Path(__file__).resolve().parents[1] /
        "docs/implementation/stage1_r2/catalogs/occupation_catalog.npz")


@pytest.fixture(params=["native_pairs", "projected", "hybrid", "native_intervals", "energy_panels", "resolved_panels"])
def events(catalog, request):
    cell = ElectronicCell(catalog, .72, .2)
    if request.param != "native_pairs":
        constructor = HybridElectronPhononEvents if request.param == "hybrid" else ProjectedElectronPhononEvents
        if request.param in ("native_intervals", "energy_panels", "resolved_panels"):
            return constructor.from_cell(cell, phonons(), rate_prefactor=1., quadrature_order=4,
                                         integration_layout=request.param)
        return constructor.from_cell(cell, phonons(), rate_prefactor=1., quadrature_order=24)
    return ElectronPhononEvents.from_cell(cell, phonons(), rate_prefactor=1.)


def test_recombination_diagonal_half_weight_and_double_stoichiometry():
    grid = PhononGrid(np.array([.5, 1., 1.5]), np.ones(3), lambda q: np.ones_like(q),
                     (.5, 1.5), "diagonal oracle")
    network = ElectronPhononEvents.from_arrays([.5], [.4], [.3], grid, rate_prefactor=2., label="one electronic bin")
    expected_rate = 2*.4**2*(1+.3**2)/2*.2**2
    dp, dn = network.rhs(np.array([.2]), np.zeros(3))
    assert len(network.omega) == 1
    assert network.rates(np.array([.2]), np.zeros(3))[0] == pytest.approx(expected_rate)
    assert dp[0] == pytest.approx(-2*expected_rate/(4*.4))
    assert dn[1] == pytest.approx(expected_rate)
    assert network.energy_rate(dp, dn) == pytest.approx(0, abs=1e-17)


def test_scattering_weight_and_particle_count():
    grid = PhononGrid(np.array([.39, .4, .41]), np.ones(3), lambda q: np.ones_like(q),
                     (.39, .41), "single scattering oracle")
    network = ElectronPhononEvents.from_arrays([.5, .9], [.2, .7], [.3, .4], grid,
                                               rate_prefactor=3., label="two electronic bins")
    dp, dn = network.rhs(np.array([0., 1.]), np.zeros(3))
    expected_rate = 3*.2*.7*(1-.3*.4)
    assert network.rates(np.array([0., 1.]), np.zeros(3))[0] == pytest.approx(expected_rate)
    assert np.dot(network.electron_capacities, dp) == pytest.approx(0, abs=1e-15)
    assert network.energy_rate(dp, dn) == pytest.approx(0, abs=1e-15)


@pytest.mark.parametrize("temperature", [.1, .3, 1.])
def test_same_temperature_fd_be_is_detailed_balance(events, temperature):
    p = expit(-events.electron_energies/temperature)
    n = 1/np.expm1(events.phonons.energies/temperature)
    forward, reverse = events.log_activities(p, n)
    scale = events.coefficients*np.exp(np.maximum(forward, reverse))
    rates = events.rates(p, n)
    assert np.max(abs(rates)/scale) < 1e-10
    dp, dn = events.rhs(p, n)
    assert abs(events.energy_rate(dp, dn)) < 1e-12


def test_barycentric_event_energy_and_count(events):
    recovered = events.beta_lower*events.phonons.energies[events.phonon_lower]
    recovered += events.beta_upper*events.phonons.energies[events.phonon_upper]
    assert np.max(abs(recovered-events.omega)) < 1e-14
    assert np.array_equal(events.beta_lower+events.beta_upper, np.ones_like(events.omega))


def test_rhs_energy_and_recombination_count(events):
    rng = np.random.default_rng(5321)
    p = rng.uniform(.01, .9, len(events.electron_energies))
    n = rng.uniform(.01, 2., len(events.phonons.energies))
    rates = events.rates(p, n)
    dp, dn = events.rhs_from_rates(rates)
    assert abs(events.energy_rate(dp, dn)) < 1e-11
    assert np.dot(events.electron_capacities, dp) == pytest.approx(-2*np.sum(rates[events.recombination]), abs=1e-12)


def test_all_physical_faces_point_inward_without_clipping(events):
    rng = np.random.default_rng(218)
    p = rng.choice([0., .25, 1.], size=len(events.electron_energies))
    n = rng.choice([0., .4], size=len(events.phonons.energies))
    dp, dn = events.rhs(p, n)
    assert np.all(dp[p == 0] >= 0)
    assert np.all(dp[p == 1] <= 0)
    assert np.all(dn[n == 0] >= 0)
    assert np.array_equal(p, np.round(p*4)/4)


def test_all_vacuum_is_stationary(events):
    dp, dn = events.rhs(np.zeros_like(events.electron_energies), np.zeros_like(events.phonons.energies))
    assert np.array_equal(dp, np.zeros_like(dp))
    assert np.array_equal(dn, np.zeros_like(dn))


@pytest.mark.parametrize("bad", [-1e-15, 1+1e-15, np.nan, np.inf])
def test_bad_pauli_data_rejected(events, bad):
    p = np.zeros_like(events.electron_energies); p[0] = bad
    with pytest.raises(ValueError, match="electronic occupations"):
        events.rhs(p, np.zeros_like(events.phonons.energies))


@pytest.mark.parametrize("bad", [-1e-15, np.nan, np.inf])
def test_bad_phonon_data_rejected(events, bad):
    n = np.zeros_like(events.phonons.energies); n[0] = bad
    with pytest.raises(ValueError, match="phonon occupations"):
        events.rhs(np.zeros_like(events.electron_energies), n)


def test_missing_coupled_phonon_bracket_is_not_discarded():
    grid = PhononGrid(np.array([.2, 1.]), np.ones(2), lambda q: np.ones_like(q),
                     (0., 4.), "deliberately insufficient mode support")
    with pytest.raises(ValueError, match="lack a phonon bracket"):
        ElectronPhononEvents.from_arrays([.5, .9], [.2, .7], [0., 0.], grid,
                                          rate_prefactor=1., label="support oracle")


def test_compact_zero_coupling_is_recorded():
    grid = phonons()
    network = ElectronPhononEvents.from_arrays([.5, 3., 6.], [.2, .7, .4], [0., 0., 0.], grid,
                                              rate_prefactor=1., label="explicit cutoff")
    assert network.metadata["outside_compact_coupling_events"] > 0
    assert np.max(network.omega) <= 4.


@pytest.mark.parametrize("ratio", [1+1e-14, np.nan])
def test_incompatible_spectral_ratio_rejected(ratio):
    with pytest.raises(ValueError, match="R2/N1"):
        ElectronPhononEvents.from_arrays([.5], [.2], [ratio], phonons(), rate_prefactor=1., label="invalid spectrum")


def test_projected_electron_energy_matches_its_exact_native_stoichiometry(catalog):
    cell = ElectronicCell(catalog, .72, 0.)
    network = ProjectedElectronPhononEvents.from_cell(cell, phonons(), rate_prefactor=1., quadrature_order=32)
    for suffix in ("i", "j"):
        represented = getattr(network, "electron_beta_lower_"+suffix)*cell.energies[getattr(network, "electron_lower_"+suffix)]
        represented += getattr(network, "electron_beta_upper_"+suffix)*cell.energies[getattr(network, "electron_upper_"+suffix)]
        assert np.max(abs(represented-getattr(network, "target_energy_"+suffix))) < 2e-15
    p = expit(-cell.energies/.2)
    n = 1/np.expm1(network.phonons.energies/.2)
    assert abs(network.energy_rate(*network.rhs(p, n))) < 1e-14


def test_projected_mixed_empty_full_endpoints_remain_defined(catalog):
    network = ProjectedElectronPhononEvents.from_cell(ElectronicCell(catalog, .72, .2), phonons(),
                                                      rate_prefactor=1., quadrature_order=24)
    p = np.arange(len(network.electron_energies)) % 2
    dp, dn = network.rhs(p.astype(float), np.zeros_like(network.phonons.energies))
    assert np.all(np.isfinite(dp)) and np.all(np.isfinite(dn))
    assert np.all(dp[p == 0] >= 0) and np.all(dp[p == 1] <= 0) and np.all(dn >= 0)


@pytest.mark.parametrize("amplitude", [.55, .72, .95, 1.47])
def test_exact_bcs_virtual_support_moves_with_amplitude(catalog, amplitude):
    cell = ElectronicCell(catalog, amplitude, 0.)
    network = ProjectedElectronPhononEvents.from_cell(cell, phonons(), rate_prefactor=1.,
        quadrature_order=2, integration_layout="resolved_panels")
    assert "exact causal" in network.metadata["electron_reconstruction"]
    for suffix in ("i", "j"):
        target = getattr(network, "target_energy_"+suffix)
        assert np.all(target >= cell.energies[0])
        assert np.all(target <= cell.energies[-1])
        beta = getattr(network, "electron_beta_upper_"+suffix)
        assert np.all((beta >= 0) & (beta <= 1))
    p = expit(-cell.energies/.23)
    n = 1/np.expm1(network.phonons.energies/.23)
    assert abs(network.energy_rate(*network.rhs(p, n))) < 1e-13


def test_resolved_partitions_at_nearly_coincident_energy_knots():
    # An interval narrower than one quadrature-node ULP is a legitimate
    # vanishing-measure interval, not permission to extrapolate a reaction.
    energy = np.array([.125, np.nextafter(.125, 1.), .25, .375, .5, .75, 1.])
    network = ProjectedElectronPhononEvents.from_arrays(energy, np.ones(len(energy)),
        energy, np.zeros(len(energy)), phonons(.125, 1., 33), rate_prefactor=1.,
        quadrature_order=2, label="near-coincident partition regression",
        integration_layout="resolved_panels", spectral_parameters=(0., 0., 1e-8))
    assert np.all(np.isfinite(network.coefficients))
    assert np.all(network.coefficients > 0)
    assert np.all(network.target_energy_i >= energy[0])
    assert np.all(network.target_energy_j <= energy[-1])
    assert np.all(network.omega >= .125)
    assert np.all(network.omega <= 1.)
