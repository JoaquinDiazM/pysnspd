import unittest
import numpy as np
from pysnspd.experimental import thermal_spatial_usadel as thermal
from pysnspd.experimental.thermal_snapshot import spectral_difference


class ThermalSnapshotTests(unittest.TestCase):
    def test_exact_difference_matches_original_action_and_derivatives(self):
        graph=thermal.rectangular_graph(5,4,3.,2.);rng=np.random.default_rng(782)
        d0=.8+.1j*rng.normal(size=graph.n_nodes);u0=.4+.2j*rng.normal(size=graph.n_nodes)
        d=d0+.03*(rng.normal(size=graph.n_nodes)+1j*rng.normal(size=graph.n_nodes))
        u=u0+.02*(rng.normal(size=graph.n_nodes)+1j*rng.normal(size=graph.n_nodes));eps=.7
        a=thermal.spectral_energy_gradient(graph,d0,eps,u0);b=thermal.spectral_energy_gradient(graph,d,eps,u)
        result=spectral_difference(graph,d0,u0,d,u,eps)
        direct=b.energy-a.energy+np.dot(graph.area_weights,abs(d)**2-abs(d0)**2)/eps
        self.assertAlmostEqual(result['renormalized_energy_difference'],direct,delta=2e-14)
        np.testing.assert_allclose(result['gap_gradient_difference'],2*graph.area_weights*((d-d0)/eps-(b.f-a.f)),atol=1e-15)
        np.testing.assert_allclose(result['current_difference'],a.link_derivative_alpha-b.link_derivative_alpha,atol=1e-15)

    def test_equal_states_are_exactly_zero_and_tiny_differences_are_retained(self):
        graph=thermal.rectangular_graph(4,4,2.,2.);d0=np.ones(graph.n_nodes,complex);u0=np.full(graph.n_nodes,.7+.2j)
        zero=spectral_difference(graph,d0,u0,d0,u0,.3)
        self.assertEqual(zero['renormalized_energy_difference'],0.)
        np.testing.assert_array_equal(zero['gap_gradient_difference'],0.)
        result=spectral_difference(graph,d0,u0,d0+1e-10,u0+1e-10,.3)
        self.assertTrue(np.all(np.isfinite(result['f_difference'])))
        self.assertGreater(np.linalg.norm(result['f_difference']),1e-11)


if __name__=='__main__':unittest.main()
