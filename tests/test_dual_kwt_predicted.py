"""A spectral predictor never replaces the nonlinear stationarity test."""
import unittest
from unittest.mock import patch
import numpy as np
from sandbox.stage4_core.dual_kwt_predicted import UniformPredictorMode, base


class DualKWTPredictorTests(unittest.TestCase):
    def test_predictor_is_verified_and_falls_back_when_nonlinear_residual_is_large(self):
        graph = base.thermal.rectangular_graph(6, 5, 5., 4.)
        d0, _, epsilon = base.uniform_reference(graph, .9/8.65, 32)
        plan = dict(spectral_tolerance=1e-7, maximum_newton_iterations=60)
        mode = UniformPredictorMode(graph, d0, epsilon[0], plan)
        x, y = graph.coordinates_bar.T
        shape = np.sin(np.pi*x/5)*np.cos(np.pi*y/4)
        shape[graph.boundary_nodes] = 0
        original = base.solve_frequency
        with patch.object(base, 'solve_frequency', wraps=original) as solve:
            small = d0*(1+1e-5*shape)*np.exp(1e-5j*shape)
            first = mode.evaluate('small', small)
            self.assertFalse(first['newton_solved'])
            self.assertEqual(solve.call_count, 0)
            self.assertLessEqual(first['residual'], plan['spectral_tolerance'])
            exact = original(graph, small, epsilon[0], fixed_nodes=graph.boundary_nodes,
                fixed_u=d0[graph.boundary_nodes]/epsilon[0], tol=1e-11)
            self.assertLess(np.max(abs(first['u']-exact.u)), 2e-7)
            large = d0*(1+.1*shape)*np.exp(.1j*shape)
            second = mode.evaluate('large', large)
            self.assertTrue(second['newton_solved'])
            self.assertEqual(solve.call_count, 1)
            self.assertGreater(second['predicted_residual'], plan['spectral_tolerance'])
            self.assertLessEqual(second['residual'], plan['spectral_tolerance'])
            self.assertTrue(np.array_equal(second['u'][graph.boundary_nodes],
                                           d0[graph.boundary_nodes]/epsilon[0]))


if __name__ == '__main__':
    unittest.main()
