"""Independent finite differences and power identities of the full-node tangent."""
import unittest
import numpy as np
from pysnspd.experimental import thermal_spatial_usadel as thermal
from pysnspd.experimental.thermal_weak_response import (
    SpectralTangent,thermal_hessian_action,ThermalKWTNormal,divergence)


class TestThermalWeakResponse(unittest.TestCase):
    def setUp(self):
        self.graph=thermal.rectangular_graph(6,5,4.,3.)
        x,y=self.graph.coordinates_bar.T
        self.d=(.8+.08*np.sin(x))*np.exp(1j*(.12*x+.08*y))
        self.t=.35;self.count=4
        self.direction=(.1*np.cos(x)+.07j*np.sin(y))*self.d
        self.direction[self.graph.boundary_nodes]=0
        self.fixed=[self.d[self.graph.boundary_nodes]/(2*np.pi*self.t*(n+.5)) for n in range(self.count)]

    def evaluated(self,d):
        spectra=[thermal.solve_frequency(self.graph,d,2*np.pi*self.t*(n+.5),
            fixed_nodes=self.graph.boundary_nodes,fixed_u=self.fixed[n],tol=1e-11)
            for n in range(self.count)]
        return spectra,thermal.evaluate_thermal(self.graph,d,self.t,spectra)

    def test_full_node_hessian_and_current_match_independent_nonlinear_differences(self):
        spectra,base=self.evaluated(self.d)
        tangents=[SpectralTangent(self.graph,self.d,s) for s in spectra]
        actions=[[t.apply(self.direction)] for t in tangents]
        H,I=thermal_hessian_action(self.graph,self.t,[self.direction],actions,[s.epsilon for s in spectra])
        h=1e-4;plus=self.evaluated(self.d+h*self.direction)[1];minus=self.evaluated(self.d-h*self.direction)[1]
        np.testing.assert_allclose(H[0],(plus.gap_gradient-minus.gap_gradient)/(2*h),atol=2e-8,rtol=2e-6)
        np.testing.assert_allclose(I[0],(plus.current_bar-minus.current_bar)/(2*h),atol=2e-9,rtol=2e-6)
        ward=np.imag(np.conj(self.direction)*base.gap_gradient+np.conj(self.d)*H[0])+divergence(self.graph,I[0])
        free=np.ones(self.graph.n_nodes,bool);free[self.graph.boundary_nodes]=False
        self.assertLess(np.max(abs(ward[free])),2e-9)
        # Reuse the very same factors for an unrelated full-node direction.
        rng=np.random.default_rng(37);other=.03*(rng.normal(size=len(self.d))+1j*rng.normal(size=len(self.d)))
        other[~free]=0
        Ho,_=thermal_hessian_action(self.graph,self.t,[other],[[t.apply(other)] for t in tangents],[s.epsilon for s in spectra])
        self.assertAlmostEqual(np.real(np.vdot(other,H[0])),np.real(np.vdot(self.direction,Ho[0])),delta=1e-10)

    def test_kwt_normal_power_and_exact_nonstationary_rhs_tangent(self):
        spectra,base=self.evaluated(self.d)
        actions=[[SpectralTangent(self.graph,self.d,s).apply(self.direction)] for s in spectra]
        H,I=thermal_hessian_action(self.graph,self.t,[self.direction],actions,[s.epsilon for s in spectra])
        make=lambda d:ThermalKWTNormal(self.graph,d,Tc_K=8.65,T_K=.35*8.65)
        model=make(self.d);observed=model.response(base.gap_gradient,base.current_bar)
        self.assertGreater(observed['kwt_loss'],0)
        self.assertGreater(observed['normal_loss'],0)
        self.assertLess(observed['continuity_max'],2e-14)
        self.assertLess(abs(observed['dissipation_residual']),2e-10)
        self.assertGreater(np.linalg.norm(observed['velocity']),1e-2) # Reference deliberately not stationary.
        tangent=model.rhs_tangent(self.direction,H[0],base.gap_gradient,base.current_bar,I[0])
        h=1e-4;vel=[]
        for sign in (1,-1):
            d=self.d+sign*h*self.direction;state=self.evaluated(d)[1]
            vel.append(make(d).response(state.gap_gradient,state.current_bar)['velocity'])
        np.testing.assert_allclose(tangent['velocity_direction'],(vel[0]-vel[1])/(2*h),atol=3e-8,rtol=3e-6)
        self.assertGreater(np.linalg.norm(tangent['baseline_mobility_gauge_correction']),1e-5)
        self.assertLess(tangent['continuity_direction_max'],1e-13)

    def test_local_mobility_derivative_and_zero_gap_are_regular(self):
        d=self.d.copy();d[12]=0
        make=lambda z:ThermalKWTNormal(self.graph,z,Tc_K=8.65,T_K=.9)
        model=make(d);force=.4+.2j+np.zeros(len(d))
        h=1e-5
        expected=(make(d+h*self.direction).local_inverse(force)-make(d-h*self.direction).local_inverse(force))/(2*h)
        np.testing.assert_allclose(model.local_inverse_derivative(self.direction,force),expected,rtol=1e-7,atol=3e-9)
        self.assertTrue(np.all(np.isfinite(model.local_inverse(force))))

    def test_fixed_boundary_motion_is_rejected(self):
        spectra,_=self.evaluated(self.d)
        tangent=SpectralTangent(self.graph,self.d,spectra[0])
        with self.assertRaisesRegex(ValueError,'zero gap directions'):
            tangent.apply(np.ones(len(self.d)))
        with self.assertRaisesRegex(ValueError,'cannot move'):
            ThermalKWTNormal(self.graph,self.d,Tc_K=8.65,T_K=.9).local_inverse_derivative(np.ones(len(self.d)),self.d)

    def test_mobility_matches_inherited_delta0_convention(self):
        from pysnspd.experimental.energy_catalog import K_B_J_K as k
        from pysnspd.experimental.cell_closures import CellScales,KWTMobility
        Tc,Tb,r=8.65,.9,1.7638769425910419
        model=ThermalKWTNormal(self.graph,self.d,Tc_K=Tc,T_K=Tb)
        scales=CellScales(r*k*Tc,1.,Tc,Tb,model.tD_ps)
        old=KWTMobility(scales)
        force=.4+.2j+np.zeros(len(self.d))
        measured=-model.local_inverse(force)
        for i in np.flatnonzero(model.free):
            z=self.d[i]/r
            g=force[i]/(r*self.graph.area_weights[i])
            expected=old.tensor_response([z.real,z.imag],[g.real,g.imag],scales.bath_temperature_bar).velocity
            np.testing.assert_allclose([measured[i].real,measured[i].imag],r*expected,atol=2e-13,rtol=2e-13)


if __name__=='__main__':unittest.main()
