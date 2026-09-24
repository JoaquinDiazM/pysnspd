"""Independent analytic controls for the new retarded spectral coordinates."""
import unittest
import numpy as np

from pysnspd.experimental import thermal_spatial_usadel as thermal
from pysnspd.experimental import retarded_spatial_usadel as retarded


class RetardedSpatialTests(unittest.TestCase):
    def setUp(self):
        self.graph = thermal.rectangular_graph(4, 3, 2., 1.5)
        self.rng = np.random.default_rng(30924)
        self.d = .7*np.exp(.2j*self.graph.coordinates_bar[:, 0])

    def test_uniform_exact_causal_bcs_across_gap(self):
        d = np.full(self.graph.n_nodes, 1.2*np.exp(.3j))
        for energy in (0., .6, 1.2, 2.4):
            z = .03-1j*energy
            result = retarded.solve(self.graph, d, z)
            root = np.sqrt(z*z+abs(d)**2)
            np.testing.assert_allclose(result.g, z/root, atol=2e-13)
            np.testing.assert_allclose(result.f, d/root, atol=2e-13)
            self.assertGreaterEqual(result.minimum_DOS, 0.)
            self.assertLess(result.normalization_residual, 2e-13)

    def test_normal_exact_zero_current(self):
        result = retarded.solve(self.graph, np.zeros(len(self.d)), .01-4j)
        np.testing.assert_array_equal(result.g, np.ones(len(self.d)))
        np.testing.assert_array_equal(result.f, np.zeros(len(self.d)))

    def test_analytic_jacobian_complex_direction(self):
        a, b = retarded.uniform_fields(self.d, .4-.8j)
        a += .03*(self.rng.normal(size=len(a))+1j*self.rng.normal(size=len(a)))
        b += .03*(self.rng.normal(size=len(b))+1j*self.rng.normal(size=len(b)))
        direction = self.rng.normal(size=(len(a), 2))+1j*self.rng.normal(size=(len(a), 2))
        alpha = self.rng.normal(size=len(self.graph.edges))*.1
        _, jac = retarded.residual_jacobian(self.graph, self.d, .4-.8j, a, b, alpha)
        h = 1e-6
        high = retarded.evaluate(self.graph, self.d, .4-.8j, a+h*direction[:, 0], b+h*direction[:, 1], alpha)
        low = retarded.evaluate(self.graph, self.d, .4-.8j, a-h*direction[:, 0], b-h*direction[:, 1], alpha)
        np.testing.assert_allclose((high.residual-low.residual).ravel()/(2*h), jac@direction.ravel(), rtol=2e-8, atol=2e-8)

    def test_action_root_gradient_identity(self):
        a, b = retarded.uniform_fields(self.d, .4-.8j)
        direction = self.rng.normal(size=(len(a), 2))+1j*self.rng.normal(size=(len(a), 2))
        value = retarded.evaluate(self.graph, self.d, .4-.8j, a, b)
        da = 2*value.residual[:, 1]/(1+a*b)**2
        db = 2*value.residual[:, 0]/(1+a*b)**2
        h = 1e-6
        high = retarded.evaluate(self.graph, self.d, .4-.8j, a+h*direction[:, 0], b+h*direction[:, 1]).action
        low = retarded.evaluate(self.graph, self.d, .4-.8j, a-h*direction[:, 0], b-h*direction[:, 1]).action
        self.assertAlmostEqual(abs((high-low)/(2*h)-np.sum(da*direction[:, 0]+db*direction[:, 1])), 0., places=7)

    def test_same_matsubara_action_fields_and_current(self):
        eps = .5
        u = thermal.solve_frequency(self.graph, self.d, eps, tol=1e-11)
        a, b = u.f/(1+u.g), np.conj(u.f)/(1+u.g)
        value = retarded.evaluate(self.graph, self.d, eps, a, b)
        original = thermal.spectral_energy_gradient(self.graph, self.d, eps, u.u)
        self.assertAlmostEqual(value.action.real, original.energy, places=11)
        self.assertAlmostEqual(value.action.imag, 0., places=11)
        np.testing.assert_allclose(value.link_derivative_alpha, original.link_derivative_alpha, atol=1e-12)
        result = retarded.solve(self.graph, self.d, eps, initial_a=a, initial_b=b, tolerance=1e-10)
        np.testing.assert_allclose(result.g, u.g, atol=2e-11)
        np.testing.assert_allclose(result.f, u.f, atol=2e-11)

    def test_gauge_covariance_including_contacts(self):
        z = .3-.6j
        result = retarded.solve(self.graph, self.d, z, tolerance=1e-11)
        chi = self.rng.normal(size=len(self.d))*.3
        alpha = chi[self.graph.edges[:, 1]]-chi[self.graph.edges[:, 0]]
        contacts = np.array([0, len(self.d)-1])
        transformed = retarded.solve(self.graph, self.d*np.exp(1j*chi), z, alpha,
            initial_a=result.a*np.exp(1j*chi), initial_b=result.b*np.exp(-1j*chi), tolerance=1e-10,
            fixed_nodes=contacts, fixed_a=result.a[contacts]*np.exp(1j*chi[contacts]),
            fixed_b=result.b[contacts]*np.exp(-1j*chi[contacts]))
        np.testing.assert_allclose(transformed.g, result.g, atol=1e-10)
        np.testing.assert_allclose(transformed.f, result.f*np.exp(1j*chi), atol=1e-10)

    def test_retarded_partner_is_not_complex_conjugate(self):
        a, b = retarded.uniform_fields(self.d, .02-.8j)
        self.assertGreater(np.max(abs(b-np.conj(a))), .2)
        _, _, good = retarded.fields(a, b)
        _, _, bad = retarded.fields(a, np.conj(a))
        self.assertGreater(np.max(abs(good-bad)), .2)

    def test_malformed_input_and_noncausal_axis_rejected(self):
        with self.assertRaises(ValueError):
            retarded.solve(self.graph, self.d, -.1-.2j)
        with self.assertRaises(ValueError):
            retarded.solve(self.graph, self.d, .1-.2j, initial_a=np.zeros(len(self.d)))
        with self.assertRaises(ValueError):
            retarded.solve(self.graph, self.d, .1-.2j, fixed_nodes=np.array([0]))


if __name__ == '__main__':
    unittest.main()
