"""Independent variational and analytic limits of the thermal graph oracle."""
import unittest

import numpy as np

from pysnspd.experimental.thermal_spatial_usadel import (
    ThermalGraph, rectangular_graph, spectral_energy_gradient,
    spectral_residual_jacobian, solve_frequency, evaluate_thermal)


class TestThermalSpatialUsadel(unittest.TestCase):
    def setUp(self):
        self.graph = rectangular_graph(5, 4, 4., 3.)
        x, y = self.graph.coordinates_bar.T
        self.d = (.8+.12*np.cos(x)*np.cos(y))*np.exp(1j*(.12*x+.06*y))
        self.u = (.3+.03*x)*np.exp(1j*(.1*x+.08*y))
        self.alpha = .03*np.sin(np.arange(len(self.graph.edges)))

    def test_cartesian_geometry_measures_and_inactive_edges(self):
        self.assertAlmostEqual(self.graph.area_weights.sum(), 12.)
        self.assertEqual(len(self.graph.boundary_nodes), 14)
        graph = ThermalGraph(self.graph.area_weights, self.graph.edges,
            np.zeros(len(self.graph.edges)))
        result = solve_frequency(graph, self.d, .7)
        np.testing.assert_allclose(result.u, self.d/.7, atol=0.)

    def test_full_real_gradient_jacobian_and_link_derivative(self):
        graph, h = self.graph, 2e-6
        dz = .15*np.cos(np.arange(graph.n_nodes))+.1j*np.sin(np.arange(graph.n_nodes))
        da = .12*np.sin(np.arange(len(graph.edges))+.3)
        value, jac = spectral_residual_jacobian(graph, self.d, .73, self.u, self.alpha)
        plus = spectral_energy_gradient(graph, self.d, .73, self.u+h*dz, self.alpha+h*da)
        minus = spectral_energy_gradient(graph, self.d, .73, self.u-h*dz, self.alpha-h*da)
        expected = np.real(np.vdot(value.gradient_u, dz))+np.dot(value.link_derivative_alpha, da)
        self.assertAlmostEqual((plus.energy-minus.energy)/(2*h), expected, delta=1e-8)
        plus = spectral_energy_gradient(graph, self.d, .73, self.u+h*dz, self.alpha)
        minus = spectral_energy_gradient(graph, self.d, .73, self.u-h*dz, self.alpha)
        predicted = jac@np.column_stack((dz.real, dz.imag)).ravel()
        measured = (plus.residual-minus.residual)/(2*h)
        np.testing.assert_allclose(predicted, np.column_stack((measured.real, measured.imag)).ravel(), atol=1e-9)

    def test_uniform_and_normal_analytic_solutions(self):
        epsilon, gap = .37, .82+.26j
        for d in (np.full(self.graph.n_nodes, gap), np.zeros(self.graph.n_nodes, complex)):
            result = solve_frequency(self.graph, d, epsilon)
            denominator = np.sqrt(epsilon**2+abs(d)**2)
            np.testing.assert_allclose(result.f, d/denominator, atol=2e-15)
            np.testing.assert_allclose(result.g, epsilon/denominator, atol=2e-15)
            self.assertEqual(result.iterations, 0)
            exact = np.sum(self.graph.area_weights*2*(epsilon-denominator))
            self.assertAlmostEqual(result.energy, exact, delta=1e-13)

    def test_damped_newton_stationarity_and_fixed_boundary_reactions(self):
        fixed = self.graph.boundary_nodes
        result = solve_frequency(self.graph, self.d, .6, self.alpha,
            fixed_nodes=fixed, fixed_u=self.d[fixed]/.6, tol=1e-9)
        self.assertLessEqual(result.residual, 1e-9)
        self.assertLess(result.iterations, 15)
        np.testing.assert_array_equal(result.u[fixed], self.d[fixed]/.6)
        free = np.ones(self.graph.n_nodes, bool)
        free[fixed] = False
        gradient = spectral_energy_gradient(self.graph, self.d, .6, result.u, self.alpha).gradient_u
        self.assertLess(np.max(abs(gradient[free])), 3e-9)
        self.assertGreater(np.max(abs(gradient[fixed])), 1e-3)

    def test_gauge_covariance_and_discrete_noether(self):
        t, count = .35, 3
        solutions = [solve_frequency(self.graph, self.d, 2*np.pi*t*(n+.5), self.alpha, tol=1e-9)
                     for n in range(count)]
        first = evaluate_thermal(self.graph, self.d, t, solutions, self.alpha)
        chi = .7*np.sin(np.arange(self.graph.n_nodes))
        tail, head = self.graph.edges.T
        alpha = self.alpha-chi[tail]+chi[head]
        d = self.d*np.exp(1j*chi)
        other = [solve_frequency(self.graph, d, solution.epsilon, alpha,
                    initial_u=solution.u*np.exp(1j*chi), tol=1e-9) for solution in solutions]
        second = evaluate_thermal(self.graph, d, t, other, alpha)
        self.assertAlmostEqual(first.energy, second.energy, delta=2e-13)
        np.testing.assert_allclose(second.gap_gradient, first.gap_gradient*np.exp(1j*chi), atol=2e-13)
        np.testing.assert_allclose(second.current_bar, first.current_bar, atol=2e-13)
        np.testing.assert_allclose(first.noether_residual, 0., atol=2e-13)

    def test_onshell_gap_force_and_current_are_energy_derivatives(self):
        t, h = .6, 5e-5
        dz = .07*np.cos(np.arange(self.graph.n_nodes))+.04j*np.sin(np.arange(self.graph.n_nodes))
        da = .03*np.sin(np.arange(len(self.graph.edges)))
        def evaluated(amount):
            d, alpha = self.d+amount*dz, self.alpha+amount*da
            solutions = [solve_frequency(self.graph, d, 2*np.pi*t*(n+.5), alpha, tol=1e-9) for n in range(3)]
            return evaluate_thermal(self.graph, d, t, solutions, alpha)
        center = evaluated(0.)
        measured = (evaluated(h).energy-evaluated(-h).energy)/(2*h)
        expected = np.real(np.vdot(center.gap_gradient, dz))+np.dot(center.link_derivative_alpha, da)
        self.assertAlmostEqual(measured, expected, delta=2e-7)
        np.testing.assert_array_equal(center.current_bar, -center.link_derivative_alpha)

    def test_wrong_inputs_and_unconverged_paths_are_not_silently_repaired(self):
        with self.assertRaises(ValueError):
            ThermalGraph(np.ones(2), np.array([[0, 1]]), np.array([-.1]))
        with self.assertRaises(ValueError):
            solve_frequency(self.graph, self.d, .7, fixed_nodes=[0])
        with self.assertRaisesRegex(RuntimeError, 'iteration limit'):
            solve_frequency(self.graph, self.d, .7, self.alpha, max_iterations=1, tol=1e-13)
        solution = solve_frequency(self.graph, self.d, np.pi*.6, self.alpha)
        with self.assertRaisesRegex(ValueError, 'these fields'):
            evaluate_thermal(self.graph, self.d+.001, .6, [solution], self.alpha)


if __name__ == '__main__':
    unittest.main()
