"""Instantaneous algebraic oracles; no admission of a spatial transient."""
import numpy as np
import pytest

from pysnspd.experimental.cell_closures import CellScales, KWTMobility
from pysnspd.experimental.cell_validation import ElectronicCell
from pysnspd.experimental.electrical_ports import solve_potential, ThesisCircuitParameters
from pysnspd.experimental.energy_catalog import build_vacuum_catalog, E_CHARGE_C, HBAR_J_S, K_B_J_K
from pysnspd.experimental.spatial_open import OpenSpatialFunctional
from pysnspd.experimental.mixed_spatial import MixedSpatialFunctional
from pysnspd.experimental.spatial_dynamics import (
    equivalent_temperatures, kwt_spatial_response, gauge_power_check,
    normal_heat_distribution, deposit_spatial_heat, spatial_energy_rate,
    cm9_power_balance,
)


class AlgebraCatalog:
    """Small analytic excitation oracle, deliberately not a material model."""
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

    def evaluate(self, a, gamma, p):
        vacuum = self.vacuum.evaluate(a, gamma)
        return tuple(v+4*np.dot(self.count_weights*k, p)
                     for v, k in zip(vacuum, self.energy_kernel(a, gamma)))


def scales(catalog, time_ps=1.):
    return CellScales(catalog.vacuum.delta0_J, catalog.vacuum.N0_per_J_m3,
                     8., .12*catalog.vacuum.delta0_J/K_B_J_K, time_ps)


def test_kwt_identification_sums_positive_temporal_matrices_before_inversion():
    cat = AlgebraCatalog()
    mobility = KWTMobility(scales(cat))
    z = np.array([0., .8+.1j, .9-.05j])
    gradient = np.array([[.1, -.2], [-.3, .07], [.2, -.1]])
    load = np.array([[0., .03], [.02, 0.], [-.01, .005]])
    mass = np.array([1., .3, .7, 2.])
    mapping = np.array([0, 1, 1, 2])
    temperature = np.array([.12, .2, .4, .15])
    phi = np.array([1e-5, 3e-6, 0.])
    response = kwt_spatial_response(z, gradient, mass, mobility, temperature, phi,
        quadrature_to_node=mapping, boundary_load_bar=load)
    # Independent real 2x2 tensor construction of D.12, not the radial solve.
    matrices = []
    total = np.zeros((3, 2, 2))
    for j, node in enumerate(mapping):
        c = mobility.coefficients(abs(z[node]), temperature[j])
        zr = np.array([z[node].real, z[node].imag])
        matrix = 2*c['Abar']*c['tau0bar']/c['R']*(np.eye(2)+c['c']*np.outer(zr, zr))
        matrices.append(matrix)
        total[node] += mass[j]*matrix
    vr = np.array([np.linalg.solve(matrix, -(g-f))
                   for matrix, g, f in zip(total, gradient, load)])
    expected = vr[:, 0]+1j*vr[:, 1]
    expected_heat = np.array([vr[node]@matrix@vr[node]
                             for node, matrix in zip(mapping, matrices)])
    np.testing.assert_allclose(response.material_velocity_bar, expected, rtol=3e-14, atol=1e-15)
    np.testing.assert_allclose(response.heat_density_bar, expected_heat, rtol=3e-14)
    np.testing.assert_allclose(response.field_velocity_bar,
        expected-1j*(2*E_CHARGE_C*mobility.scales.t_ref_s/HBAR_J_S)*phi*z, rtol=3e-14)
    assert np.all(response.heat_density_bar > 0)
    assert abs(response.identity_residual_bar) < 1e-14
    np.testing.assert_array_equal(response.node_mass_bar, [1., 1., 2.])


def test_reference_time_changes_normalized_rates_but_not_physical_rates():
    cat = AlgebraCatalog()
    z, g = np.array([.9+.1j]), np.array([[.03, .008]])
    first = kwt_spatial_response(z, g, [2.], KWTMobility(scales(cat, 1.)), [.2], [1e-5])
    second = kwt_spatial_response(z, g, [2.], KWTMobility(scales(cat, 2.)), [.2], [1e-5])
    np.testing.assert_allclose(first.field_velocity_bar/first.time_scale_s,
                               second.field_velocity_bar/second.time_scale_s, rtol=3e-15)
    np.testing.assert_allclose(first.heat_density_bar/first.time_scale_s,
                               second.heat_density_bar/second.time_scale_s, rtol=3e-15)


def test_temperature_inversion_uses_the_instantaneous_spectrum():
    cat = AlgebraCatalog()
    before, after = ElectronicCell(cat, .9, .04), ElectronicCell(cat, .7, .1)
    p = before.fermi_dirac(.24)
    temperatures = equivalent_temperatures([before, after], [p, p])
    assert temperatures[0] == pytest.approx(.24, abs=2e-14)
    assert abs(temperatures[1]-.24) > .01
    assert equivalent_temperatures([before], [np.zeros_like(p)])[0] == 0
    with pytest.raises(ValueError, match='finite positive-temperature'):
        equivalent_temperatures([before], [np.ones_like(p)])


def test_common_energy_kwt_potential_heat_and_thesis_circuit_balance():
    cat = AlgebraCatalog()
    sc = scales(cat)
    model = OpenSpatialFunctional(cat, 360e-9, 120e-9*7e-9, 1, degree=2)
    z = np.array([.9, .94*np.exp(.06j), .91*np.exp(.13j)])
    p = np.array([(.07+.003*i)*np.exp(-cat.count_nodes/.6) for i in range(3)])
    value = model.evaluate(z, p, require_stability=False)
    cells = [ElectronicCell(cat, a, g) for a, g in
             zip(value.amplitude_nodes_bar, value.gamma_nodes_bar)]
    edges = np.array([[0, 1], [1, 2]])
    sigma = 2*E_CHARGE_C**2*cat.vacuum.N0_per_J_m3*cat.vacuum.D_m2_s
    conductance = sigma*model.cross_section_m2/(np.diff(model.x_bar)*model.ell0_m)
    current = 5e-6
    potential = solve_potential(edges, conductance, value.current_A,
        np.array([current, 0., -current]), left_node=0, right_node=2)
    current_scale = 2*E_CHARGE_C/HBAR_J_S*model.energy_scale_J
    # Explicit terminal phase traction; this is work supplied by reservoirs.
    phase_load = np.array([-current/current_scale, 0., current/current_scale])
    complex_load = 1j*z/abs(z)**2*phase_load
    load = np.column_stack((complex_load.real, complex_load.imag))
    response = kwt_spatial_response(z, value.gradient_cartesian_bar, model.mass_bar,
        KWTMobility(sc), equivalent_temperatures(cells, p), potential.phi_V,
        boundary_load_bar=load)
    power = gauge_power_check(z, value.gradient_cartesian_bar, edges,
        value.link_current_bar, potential.phi_V,
        energy_scale_J=model.energy_scale_J, time_scale_s=sc.t_ref_s)
    np.testing.assert_allclose(power.noether_residual_bar, 0., atol=2e-16)
    assert power.field_phase_power_W == pytest.approx(potential.super_power_W, rel=2e-13, abs=1e-23)
    assert abs(power.power_residual_W) < 1e-22
    # Manufactured conservative allocation; no spatial-order claim is made.
    heat = normal_heat_distribution(conductance, potential.delta_phi_V,
        np.array([[.5, 0.], [.5, .5], [0., .5]]))
    deposition = deposit_spatial_heat(cells, p, response, heat.quadrature_power_W,
                                      energy_scale_J=model.energy_scale_J, scales=sc)
    assert abs(deposition.moment_residual_bar) < 2e-14
    rate = spatial_energy_rate(value.gradient_cartesian_bar, response.field_velocity_bar,
                               cells, deposition.electron_rhs, model.mass_bar)
    physical_rate = rate.total_rate_bar*model.energy_scale_J/sc.t_ref_s
    boundary_work = response.boundary_work_rate_bar*model.energy_scale_J/sc.t_ref_s
    assert physical_rate == pytest.approx(potential.port_power_W+boundary_work,
                                         rel=2e-12, abs=2e-22)
    circuit = ThesisCircuitParameters(Lk_ext_H=9.8e-9)
    circuit_balance = circuit.power_balance(np.array([5.1e-6, current, 0.]), potential.Vdev_V)
    balance = cm9_power_balance(rate.total_rate_bar, energy_scale_J=model.energy_scale_J,
        time_scale_s=sc.t_ref_s, circuit_balance=circuit_balance,
        reservoir_power_into_W=boundary_work)
    assert abs(balance.device_port_residual_W) < 2e-22
    assert abs(balance.residual_W) < 2e-21

    # Independent directional derivative of the WHOLE actual common energy.
    # Gamma is recomputed from the perturbed field, p stays on fixed count.
    dp = np.asarray(deposition.electron_rhs)
    h = 2e-5
    def energy(offset):
        return model.evaluate(z+offset*response.field_velocity_bar,
                              p+offset*dp, require_stability=False).energy_bar
    measured = (energy(-2*h)-8*energy(-h)+8*energy(h)-energy(2*h))/(12*h)
    assert measured == pytest.approx(rate.total_rate_bar, rel=2e-7, abs=2e-9)
    moved = model.evaluate(z+h*response.field_velocity_bar, p, require_stability=False)
    assert np.max(abs(moved.gamma_nodes_bar-value.gamma_nodes_bar)) > 1e-10


def test_mixed_trace_shares_field_velocity_but_retains_separate_spectral_heating():
    cat = AlgebraCatalog()
    sc = scales(cat)
    model = MixedSpatialFunctional(cat, 120e-9, 40e-9, 7e-9, 60e-9, 60e-9,
        elements_x=1, elements_y=1, left_elements=1, right_elements=1, degree=2)
    x, y = model.dof_coordinates_bar.T
    z = (.92+.008*np.cos(x/3))*np.exp(1j*(.015*x+.003*y))
    p = np.array([(.065+.0007*i)*np.exp(-cat.count_nodes/.6)
                  for i in range(model.quadrature_size)])
    value = model.evaluate(z, p, require_stability=False)
    fields = value.sampled_fields
    cells = [ElectronicCell(cat, a, g) for a, g in
             zip(fields.amplitude_quadrature_bar, fields.gamma_quadrature_bar)]
    temperature = equivalent_temperatures(cells, p)
    shared = model.quadrature_to_dof == model.interface_dofs['left']
    assert np.ptp(temperature[shared]) > 1e-3
    sigma = 2*E_CHARGE_C**2*cat.vacuum.N0_per_J_m3*cat.vacuum.D_m2_s
    conductance = model.conductance(sigma)
    injections = np.zeros(model.cells)
    injections[0], injections[-1] = 2e-6, -2e-6
    potential = solve_potential(model.graph_edges, conductance, value.current_A,
        injections, left_node=0, right_node=model.cells-1)
    response = kwt_spatial_response(z, value.gradient_cartesian_bar,
        model.quadrature_mass_bar, KWTMobility(sc), temperature, potential.phi_V,
        quadrature_to_node=model.quadrature_to_dof)
    np.testing.assert_allclose(response.node_mass_bar, model.mass_bar, rtol=3e-15)
    assert np.ptp(response.heat_density_bar[shared]) > 1e-8
    allocation = np.zeros((model.quadrature_size, len(model.graph_edges)))
    # Explicit test allocation: half per endpoint, split among its quadratures
    # in proportion to their volume. It does not claim a spatial convergence order.
    for face, (tail, head) in enumerate(model.graph_edges):
        for node in (tail, head):
            rows = model.quadrature_to_dof == node
            allocation[rows, face] += .5*model.quadrature_mass_bar[rows]/model.mass_bar[node]
    heat = normal_heat_distribution(conductance, potential.delta_phi_V, allocation)
    deposition = deposit_spatial_heat(cells, p, response, heat.quadrature_power_W,
                                      energy_scale_J=model.energy_scale_J, scales=sc)
    rate = spatial_energy_rate(value.gradient_cartesian_bar, response.field_velocity_bar,
                               cells, deposition.electron_rhs, model.quadrature_mass_bar)
    assert rate.total_rate_bar*model.energy_scale_J/sc.t_ref_s == pytest.approx(
        potential.port_power_W, rel=2e-11, abs=1e-22)
    power = gauge_power_check(z, value.gradient_cartesian_bar, model.graph_edges,
        value.link_current_bar, potential.phi_V, energy_scale_J=model.energy_scale_J,
        time_scale_s=sc.t_ref_s)
    assert np.max(abs(power.noether_residual_bar)) < 2e-15
    assert abs(power.power_residual_W) < 1e-22
    circuit_balance = ThesisCircuitParameters(Lk_ext_H=9.8e-9).power_balance(
        [2.1e-6, 2e-6, 0.], potential.Vdev_V)
    balance = cm9_power_balance(rate.total_rate_bar, energy_scale_J=model.energy_scale_J,
        time_scale_s=sc.t_ref_s, circuit_balance=circuit_balance)
    assert abs(balance.residual_W) < 2e-21
    h = 2e-5
    dp = np.asarray(deposition.electron_rhs)
    def energy(offset):
        return model.evaluate(z+offset*response.field_velocity_bar,
                              p+offset*dp, require_stability=False).energy_bar
    measured = (energy(-2*h)-8*energy(-h)+8*energy(h)-energy(2*h))/(12*h)
    assert measured == pytest.approx(rate.total_rate_bar, rel=2e-7, abs=2e-9)


@pytest.mark.parametrize('mass,mapping', [([0.], None), ([1., 2.], None),
    ([1.], [0.]), ([1.], [1]), ([1., 2.], np.array([[1., 0.], [0., 1.]]))])
def test_rejects_invalid_mass_and_identification_without_repair(mass, mapping):
    with pytest.raises(ValueError):
        kwt_spatial_response([.9], [[.1, .2]], mass, KWTMobility(scales(AlgebraCatalog())),
                             np.full(len(mass), .2), [0.], quadrature_to_node=mapping)


@pytest.mark.parametrize('allocation', [[[.2], [.7]], [[-.1], [1.1]], [[np.nan], [1.]]])
def test_heat_allocation_requires_a_positive_partition_of_each_face(allocation):
    with pytest.raises(ValueError):
        normal_heat_distribution([1.], [.2], allocation)


def test_explicit_gauge_link_work_is_separate_and_requires_both_arguments():
    cat = AlgebraCatalog()
    cell = ElectronicCell(cat, .9, .1)
    zero = np.zeros_like(cat.count_nodes)
    result = spatial_energy_rate([[0., 0.]], [0j], [cell], [zero], [2.],
        phonon_density_rates_bar=[-.03], link_current_bar=[.2, -.1],
        link_phase_velocity=[.3, -.4])
    assert result.total_rate_bar == pytest.approx(.04)
    assert result.explicit_link_work_bar == pytest.approx(.1)
    with pytest.raises(ValueError, match='both explicit'):
        spatial_energy_rate([[0., 0.]], [0j], [cell], [zero], [2.], link_current_bar=[.1])


def test_heat_rejects_mismatched_spectrum_and_saturated_occupations():
    cat = AlgebraCatalog()
    sc = scales(cat)
    response = kwt_spatial_response([.9], [[.1, .02]], [1.], KWTMobility(sc), [.2], [0.])
    with pytest.raises(ValueError, match='amplitude'):
        deposit_spatial_heat([ElectronicCell(cat, .8, .1)], [np.zeros(6)], response, [0.],
                             energy_scale_J=1e-20, scales=sc)
    with pytest.raises(ValueError, match='saturated'):
        deposit_spatial_heat([ElectronicCell(cat, .9, .1)], [np.ones(6)], response, [0.],
                             energy_scale_J=1e-20, scales=sc)
