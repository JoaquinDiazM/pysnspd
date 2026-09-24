"""Physical algebra controls for a static two-mode Keldysh response."""
import unittest
import numpy as np
from pysnspd.experimental import thermal_spatial_usadel as thermal
from pysnspd.experimental import retarded_spatial_usadel as retarded
from pysnspd.experimental import frozen_kinetic_usadel as kinetic


class FrozenKineticTests(unittest.TestCase):
    def setUp(self):
        self.graph = thermal.rectangular_graph(5, 5, 2., 2.)
        self.xy = self.graph.coordinates_bar
        self.rng = np.random.default_rng(94026)

    def spectral(self):
        d = 1.2*np.exp(.4j*self.xy[:, 0])
        spectrum = retarded.solve(self.graph, d, .05-.8j, tolerance=1e-11)
        operator = kinetic.assemble(self.graph, d, spectrum.g, spectrum.f, spectrum.f_tilde)
        return d, spectrum, operator

    def test_normal_limit_exact_discrete_diffusion(self):
        size = self.graph.n_nodes
        operator = kinetic.assemble(self.graph, np.zeros(size), np.ones(size), np.zeros(size), np.zeros(size))
        h = self.rng.normal(size=(size, 2))
        observed = operator.evaluate(h)
        tail, head = self.graph.edges.T
        expected = self.graph.conductance[:, None]*(h[head]-h[tail])
        np.testing.assert_allclose(observed.edge_flux, expected, atol=1e-14)
        np.testing.assert_array_equal(observed.pair_conversion, np.zeros_like(h))
        np.testing.assert_allclose(operator.matrix@h.ravel(), observed.residual.ravel(), atol=1e-14)

    def test_uniform_charge_conversion_sign_and_diffusion(self):
        d = np.full(self.graph.n_nodes, 1.2)
        a, b = retarded.uniform_fields(d, .05-.8j)
        f, ft, g = retarded.fields(a, b)
        operator = kinetic.assemble(self.graph, d, g, f, ft)
        h = np.column_stack((np.zeros(len(d)), self.rng.normal(size=len(d))))
        value = operator.evaluate(h)
        np.testing.assert_allclose(value.pair_conversion[:, 1], -self.graph.area_weights*d*f.real*h[:, 1], atol=1e-14)
        DT = (1+abs(g[0])**2+abs(f[0])**2)/2
        tail, head = self.graph.edges.T
        np.testing.assert_allclose(value.edge_flux[:, 1], self.graph.conductance*DT*(h[head, 1]-h[tail, 1]), atol=1e-14)

    def test_ward_identity_without_kinetic_or_spectral_stationarity(self):
        d = (1+.2*self.xy[:, 0])*np.exp(.4j*self.xy[:, 1])
        a, b = retarded.uniform_fields(d, .06-.9j)
        f, ft, g = retarded.fields(a, b)
        operator = kinetic.assemble(self.graph, d, g, f, ft)
        h = self.rng.normal(size=(len(d), 2))
        value = operator.evaluate(h)
        self.assertLess(np.max(abs(value.ward_residual)), 2e-14)
        np.testing.assert_allclose(operator.matrix@h.ravel(), value.residual.ravel(), atol=2e-14)

    def test_uniform_energy_mode_at_finite_eta_is_equilibrium(self):
        d, spectrum, operator = self.spectral()
        value = operator.evaluate(np.column_stack((np.ones(len(d)), np.zeros(len(d)))), eta=.05)
        self.assertLess(np.max(abs(value.residual)), 1e-10)
        self.assertGreater(np.max(abs(value.artificial_eta_leakage)), .001)

    def test_charge_response_closes_only_charge_not_energy_equation(self):
        d, spectrum, operator = self.spectral()
        r = np.linalg.norm(self.xy-[1., 0.], axis=1)
        hL = np.maximum(0., 1-r*r)**3
        directions = operator.charge_response(hL, self.graph.boundary_nodes)
        value = operator.evaluate(directions)
        free = np.ones(len(d), bool); free[self.graph.boundary_nodes] = False
        self.assertGreater(np.max(abs(directions[:, 1])), 1e-4)
        self.assertLess(np.max(abs(value.residual[free, 1])), 1e-13)
        self.assertGreater(np.max(abs(value.residual[free, 0])), .01)
        self.assertLess(np.max(abs(value.ward_residual)), 2e-14)
        np.testing.assert_array_equal(directions[self.graph.boundary_nodes], np.zeros((len(self.graph.boundary_nodes), 2)))

    def test_gauge_covariance_of_spectral_response(self):
        d, spectrum, operator = self.spectral()
        chi = self.rng.normal(size=len(d))*.3
        alpha = chi[self.graph.edges[:, 1]]-chi[self.graph.edges[:, 0]]
        gauged = kinetic.assemble(self.graph, d*np.exp(1j*chi), spectrum.g,
            spectrum.f*np.exp(1j*chi), spectrum.f_tilde*np.exp(-1j*chi), alpha)
        h = self.rng.normal(size=(len(d), 2))
        first, second = operator.evaluate(h), gauged.evaluate(h)
        np.testing.assert_allclose(first.residual, second.residual, atol=2e-14)
        np.testing.assert_allclose(first.edge_flux, second.edge_flux, atol=2e-14)
        np.testing.assert_allclose(first.gap_force_increment*np.exp(1j*chi), second.gap_force_increment, atol=2e-14)

    def test_zero_distribution_increment_does_not_double_count_thermal_fields(self):
        d, spectrum, operator = self.spectral()
        value = operator.evaluate(np.zeros((len(d), 2)))
        np.testing.assert_array_equal(value.gap_force_increment, np.zeros(len(d)))
        np.testing.assert_array_equal(value.charge_current_increment, np.zeros(len(self.graph.edges)))

    def test_contacts_and_real_directions_are_explicit(self):
        d, spectrum, operator = self.spectral()
        with self.assertRaises(ValueError):
            operator.charge_response(np.zeros(len(d)), np.array([], int))
        with self.assertRaises(ValueError):
            operator.charge_response(np.ones(len(d)), self.graph.boundary_nodes)
        with self.assertRaises(ValueError):
            operator.evaluate(np.zeros((len(d), 2), complex))


if __name__ == '__main__':
    unittest.main()
