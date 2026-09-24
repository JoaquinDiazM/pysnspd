"""A resident Jacobian must preserve the exact nonlinear spectral solve."""
import unittest

import numpy as np

from pysnspd.experimental import thermal_spatial_usadel as thermal
from pysnspd.experimental.bulk_current_reference import homogeneous_spectrum
from sandbox.stage5_prephoton.dc_hold import _mode_contribution


class PrephotonSpectralCacheTests(unittest.TestCase):
    def setUp(self):
        original=thermal.rectangular_graph(5,4,4.,3.)
        x=original.coordinates_bar[:,0]
        fixed=np.flatnonzero((x==x.min()) | (x==x.max()))
        self.graph=thermal.ThermalGraph(original.area_weights,original.edges,
            original.conductance,original.coordinates_bar,fixed)
        self.fixed=fixed;self.x=x;self.epsilon=.37
        phase=np.exp(.15j*x)
        self.gap=1.4*phase
        scalar=homogeneous_spectrum(1.4,np.array([self.epsilon]),2*(1-np.cos(.15)))[0][0]
        self.u=scalar*phase
        self.plan=dict(spectral_tolerance=1e-10,maximum_newton_iterations=60)

    def compare(self,gap,u,fixed_u,cache):
        cached,observed=_mode_contribution(self.graph,gap,self.epsilon,u,fixed_u,self.plan,cache)
        direct,exact=_mode_contribution(self.graph,gap,self.epsilon,u,fixed_u,self.plan)
        np.testing.assert_allclose(cached,direct,rtol=2e-10,atol=2e-10)
        np.testing.assert_allclose(observed['force'],exact['force'],rtol=0,atol=2e-10)
        np.testing.assert_allclose(observed['current'],exact['current'],rtol=0,atol=2e-10)
        self.assertAlmostEqual(observed['energy'],exact['energy'],delta=2e-10)
        evaluated=thermal.spectral_energy_gradient(self.graph,gap,self.epsilon,cached)
        free=np.ones(self.graph.n_nodes,bool);free[self.fixed]=False
        residual=np.max(abs(evaluated.residual[free])/
            (self.graph.area_weights[free]*np.maximum(1.,abs(gap[free]))))
        self.assertLessEqual(residual,self.plan['spectral_tolerance'])
        np.testing.assert_array_equal(cached[self.fixed],fixed_u)
        return cached,observed

    def test_one_factor_serves_successive_perturbations_at_unchanged_tolerance(self):
        cache={};old=self.u
        for index in (1,2,3):
            gap=self.gap*(1+index*1e-6*np.sin(np.pi*self.x/4.))
            old,part=self.compare(gap,old,self.u[self.fixed],cache)
            self.assertEqual(part['newton_solves'],0)
            self.assertEqual(part['chord_corrections'],1)
            self.assertEqual(part['spectral_factorizations'],int(index==1))
            if index==1:
                initial_factor=cache['factor']
            else:
                self.assertIs(cache['factor'],initial_factor)

    def test_moving_reservoir_phase_is_included_in_the_nonlinear_residual(self):
        cache={};old=self.u
        for angle in (1e-6,2e-6):
            phase=np.exp(1j*angle)
            gap=self.gap*phase
            old,part=self.compare(gap,old,self.u[self.fixed]*phase,cache)
            self.assertEqual(part['newton_solves'],0)
            self.assertEqual(part['chord_corrections'],1)
            np.testing.assert_allclose(old,self.u*phase,rtol=2e-10,atol=2e-10)

    def test_a_bad_cached_predictor_falls_back_to_newton_and_refreshes(self):
        class WrongPredictor:
            def solve(self,rhs):
                return np.zeros_like(rhs)
        bad=WrongPredictor();cache={'factor':bad}
        gap=self.gap*(1+1e-4*np.sin(np.pi*self.x/4.))
        _,part=self.compare(gap,self.u,self.u[self.fixed],cache)
        self.assertEqual(part['chord_corrections'],0)
        self.assertEqual(part['newton_solves'],1)
        self.assertEqual(part['spectral_factorizations'],1)
        self.assertIsNot(cache['factor'],bad)


if __name__=='__main__':
    unittest.main()
