"""Check the invariant-sector shortcut against the existing full graph laws."""
import unittest
import numpy as np
from pysnspd.experimental import thermal_spatial_usadel as thermal
from pysnspd.experimental.thermal_weak_response import SpectralTangent
from sandbox.stage4_core.dual_longitudinal_mode import make_control, evolve


class DualLongitudinalModeTests(unittest.TestCase):
    def test_exact_uniform_hessian_and_nonthermal_exchange_on_existing_operators(self):
        graph = thermal.rectangular_graph(7, 5, 8., 4.)
        control = make_control(graph, matsubara_count=12)
        d = np.full(graph.n_nodes, control['gap'], complex)
        v = control['mode']
        actual = 2*graph.area_weights*np.log(control['temperature_ratio'])*v
        for epsilon in control['eps']:
            solution = thermal.solve_frequency(graph, d, epsilon,
                fixed_nodes=graph.boundary_nodes,
                fixed_u=d[graph.boundary_nodes]/epsilon, tol=1e-10)
            direction = SpectralTangent(graph, d, solution).apply(v)
            actual += 4*np.pi*control['temperature_ratio']*graph.area_weights*(
                v/epsilon-direction.df.real)
        np.testing.assert_allclose(actual, control['hessian']*graph.area_weights*v,
            rtol=1e-10, atol=1e-12)
        self.assertLess(control['embedding_error'], 1e-10)
        self.assertLess(abs(control['reciprocal_power_residual']), 1e-14)
        # Suppressing reciprocal gap/population coupling changes the RHS:
        # the control genuinely exercises exchange, not two uncoupled curves.
        full = control['matrix']@control['initial']
        uncoupled_gap = -control['mobility']*control['hessian']*control['initial'][0]/control['tD_ps']
        self.assertGreater(abs(full[0]-uncoupled_gap), 1e-9)

    def test_euler_refinement_integrated_availability_and_occupation_support(self):
        control = make_control(thermal.rectangular_graph(7, 5, 8., 4.), matsubara_count=12)
        a = evolve(control, steps=512, duration_ps=.1)
        b = evolve(control, steps=1024, duration_ps=.1)
        ratio = a['metric_error_over_initial']/b['metric_error_over_initial']
        self.assertTrue(1.8 < ratio < 2.2)
        self.assertLess(a['integrated_balance_over_initial'], .01)
        self.assertGreaterEqual(a['minimum_occupation'], 0)
        self.assertLessEqual(a['maximum_occupation'], 1)
        self.assertTrue(np.all(np.diff(a['availability']) <= 1e-15))
        self.assertTrue(np.all(a['integrated_loss'][-1] > 0))


if __name__ == '__main__':
    unittest.main()
