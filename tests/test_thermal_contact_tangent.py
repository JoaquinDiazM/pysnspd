import unittest
import numpy as np

from pysnspd.experimental import thermal_spatial_usadel as thermal
from pysnspd.experimental.thermal_stable_newton import solve_frequency
from pysnspd.experimental.thermal_contact_tangent import ThermalMovingContactTangent, ThermalContactTangent
from sandbox.stage4_core.harmonic_uniform_reference import end_contacts


class ThermalMovingContactTangentTest(unittest.TestCase):
    def setUp(self):
        self.graph = thermal.rectangular_graph(6, 4, 4., 2.)
        self.fixed = np.unique(np.r_[end_contacts(self.graph)])
        x, y = self.graph.coordinates_bar.T
        self.d = (1.2+.03*np.cos(.5*x))*np.exp(1j*(.07*x+.02*np.sin(y)))
        self.epsilon = .7
        self.solution = solve_frequency(self.graph, self.d, self.epsilon,
            fixed_nodes=self.fixed, fixed_u=self.d[self.fixed]/self.epsilon, tol=1e-12)
        self.tangent = ThermalMovingContactTangent(self.graph, self.d, self.solution)

    def test_global_gauge_rotation_moves_contacts_and_preserves_current(self):
        response = self.tangent.apply(1j*self.d)
        np.testing.assert_allclose(response.du, 1j*self.solution.u, atol=2e-12)
        np.testing.assert_allclose(response.df, 1j*self.solution.f, atol=2e-12)
        np.testing.assert_allclose(response.dg, 0., atol=2e-13)
        np.testing.assert_allclose(response.current_derivative, 0., atol=2e-13)
        self.assertLess(response.free_equation_residual_max, 1e-13)
        self.assertEqual(response.fixed_law_residual_max, 0.)
        # A reservoir may sustain a reaction even though its free neighbours
        # solve their spectral equations. Its gauge derivative remains visible.
        np.testing.assert_allclose(response.boundary_residual_derivative,
            1j*self.tangent.value.residual[self.fixed], atol=2e-12)
        self.assertGreater(np.max(abs(response.boundary_residual_derivative)), 1e-3)

    def test_exact_tangent_matches_new_stationary_solutions_with_moving_contacts(self):
        x, y = self.graph.coordinates_bar.T
        direction = .1*np.cos(.4*x)+.07j*np.sin(.5*y+.2*x)
        response = self.tangent.apply(direction)
        resident = ThermalContactTangent(self.graph, self.d, self.epsilon, self.solution.u,
                                        fixed_nodes=self.fixed).apply(direction)
        np.testing.assert_allclose(resident.df, response.df, atol=0.)
        np.testing.assert_allclose(resident.boundary_reaction, response.boundary_reaction, atol=0.)
        self.assertEqual(resident.equation_residual, response.free_equation_residual_max)
        step = 2e-5
        snapshots = []
        for sign in (-1., 1.):
            gap = self.d+sign*step*direction
            solution = solve_frequency(self.graph, gap, self.epsilon,
                initial_u=self.solution.u+sign*step*response.du,
                fixed_nodes=self.fixed, fixed_u=gap[self.fixed]/self.epsilon, tol=1e-13)
            snapshots.append(thermal.spectral_energy_gradient(
                self.graph, gap, self.epsilon, solution.u))
        before, after = snapshots
        np.testing.assert_allclose(response.df, (after.f-before.f)/(2*step), rtol=2e-7, atol=2e-10)
        np.testing.assert_allclose(response.dg, (after.g-before.g)/(2*step), rtol=2e-7, atol=2e-10)
        np.testing.assert_allclose(response.current_derivative,
            -(after.link_derivative_alpha-before.link_derivative_alpha)/(2*step), rtol=2e-7, atol=2e-10)
        np.testing.assert_allclose(response.residual_derivative,
            (after.residual-before.residual)/(2*step), rtol=2e-6, atol=2e-9)
        np.testing.assert_allclose(response.spectral_gradient_derivative,
            (after.gradient_u-before.gradient_u)/(2*step), rtol=2e-6, atol=2e-9)
        self.assertLess(response.free_equation_residual_max, 1e-13)
        self.assertGreater(np.max(abs(response.boundary_gradient_derivative)), 1e-3)


if __name__ == '__main__':
    unittest.main()
