"""Physical invariants of the effective neutral BCS energy closure."""
import unittest
import numpy as np
from scipy.constants import Boltzmann as kB, hbar, elementary_charge as e
from pysnspd.experimental.thermal_spatial_usadel import rectangular_graph
from pysnspd.experimental.effective_bcs_functional import (
    EffectiveBCSFunctional, bcs_count_quadrature)


class TestEffectiveBCSFunctional(unittest.TestCase):
    def setUp(self):
        self.Tc, self.D, self.N0 = 8.65, 5e-5, 9.153297821e46
        self.Eu = kB*self.Tc
        self.ell = np.sqrt(hbar*self.D/(2*self.Eu))
        self.graph = rectangular_graph(4,3,4.,2.)
        self.xi,self.w = bcs_count_quadrature(32.,6)
        self.model = EffectiveBCSFunctional(self.graph,
            energy_scale_J=self.Eu,N0_per_J_m3=self.N0,ell0_m=self.ell,
            thickness_m=7e-9,count_nodes=self.xi,count_weights=self.w,
            gradient_law="bcs_radial",bath_temperature_K=.9,diffusion_m2_s=self.D)
        x,y = self.graph.coordinates_bar.T
        self.d=(1.4+.2*np.cos(x))*np.exp(1j*(.22*x+.04*y))
        self.p=self.model.fermi_dirac(self.d,5.)+.01*np.exp(-self.xi)[None,:]
        self.alpha=.03*np.sin(np.arange(len(self.graph.edges)))
        self.dd=.15*np.cos(np.arange(len(x)))+.09j*np.sin(np.arange(len(x)))
        self.dp=.02*self.p*(1-self.p)*np.cos(np.arange(len(x))[:,None]+self.xi)
        self.da=.02*np.cos(np.arange(len(self.graph.edges)))

    def test_energy_force_current_and_population_derivative(self):
        m,h=self.model,2e-6
        value=m.evaluate(self.d,self.p,self.alpha)
        direct=(m.evaluate(self.d+h*self.dd,self.p+h*self.dp,self.alpha+h*self.da).energy_J
                -m.evaluate(self.d-h*self.dd,self.p-h*self.dp,self.alpha-h*self.da).energy_J)/(2*h)
        direction=m.direction(self.d,self.p,self.dd,self.dp,self.alpha,self.da)
        self.assertAlmostEqual(direct/direction.energy_direction_J,1.,delta=2e-8)
        np.testing.assert_allclose(value.edge_current_A,
            -2*e/hbar*value.link_derivative_alpha_J,rtol=1e-15,atol=0.)

    def test_analytic_full_hessian_and_current_direction(self):
        m,h=self.model,2e-6
        plus=m.evaluate(self.d+h*self.dd,self.p+h*self.dp,self.alpha+h*self.da)
        minus=m.evaluate(self.d-h*self.dd,self.p-h*self.dp,self.alpha-h*self.da)
        direct=m.direction(self.d,self.p,self.dd,self.dp,self.alpha,self.da)
        np.testing.assert_allclose((plus.gap_force_J_per_gapbar-minus.gap_force_J_per_gapbar)/(2*h),
            direct.gap_force_direction_J_per_gapbar,rtol=3e-8,atol=2e-32)
        np.testing.assert_allclose((plus.edge_current_A-minus.edge_current_A)/(2*h),
            direct.edge_current_direction_A,rtol=2e-7,atol=2e-15)

    def test_gauge_covariance_and_ward_identity(self):
        m=self.model
        value=m.evaluate(self.d,self.p,self.alpha)
        chi=.5*np.sin(np.arange(self.graph.n_nodes))
        t,h=self.graph.edges.T
        rotated=m.evaluate(self.d*np.exp(1j*chi),self.p,self.alpha-chi[t]+chi[h])
        self.assertAlmostEqual(rotated.energy_J/value.energy_J,1.,delta=4e-15)
        np.testing.assert_allclose(rotated.gap_force_J_per_gapbar,
            value.gap_force_J_per_gapbar*np.exp(1j*chi),rtol=4e-14,atol=1e-33)
        np.testing.assert_allclose(rotated.edge_current_A,value.edge_current_A,
            rtol=4e-14,atol=1e-21)
        np.testing.assert_allclose(value.current_divergence_A,
            4*e/hbar*np.imag(np.conj(value.gap_force_J_per_gapbar)*self.d),
            atol=4e-21,rtol=1e-13)

    def test_normal_point_is_smooth_without_core_cutoff(self):
        m=self.model
        d=np.zeros(self.graph.n_nodes,complex)
        p=m.fermi_dirac(d,8.65)
        value=m.evaluate(d,p)
        np.testing.assert_array_equal(value.gap_force_J_per_gapbar,0)
        np.testing.assert_array_equal(value.edge_current_A,0)
        kap,kp,kpp=m.rigidity(np.array([0.,1e-22,1e-10,.01,1.]))
        self.assertTrue(np.all(np.isfinite(kap)))
        self.assertTrue(np.all(np.isfinite(kp)))
        self.assertTrue(np.all(np.isfinite(kpp)))
        self.assertTrue(np.all(kap>0))
        self.assertLess(abs(kap[1]-kap[0])/kap[0],1e-14)
        at_Tc=m.bulk_coefficient(0.,8.65)
        self.assertLess(abs(at_Tc),1e-15)
        self.assertLess(m.bulk_coefficient(0.,.9),0)
        self.assertGreater(m.bulk_coefficient(0.,10.),0)

    def test_constant_gradient_limit_and_no_double_pair_breaking(self):
        k=self.N0*np.pi*hbar*self.D/(8*self.Eu)
        m=EffectiveBCSFunctional(self.graph,energy_scale_J=self.Eu,
            N0_per_J_m3=self.N0,ell0_m=self.ell,thickness_m=7e-9,
            count_nodes=self.xi,count_weights=self.w,gradient_stiffness_SI=k)
        d=np.full(self.graph.n_nodes,1.5+0j)
        p=m.fermi_dirac(d,4.)
        value=m.evaluate(d,p)
        self.assertEqual(value.gradient_energy_J,0.)
        np.testing.assert_allclose(m.rigidity(np.array([0.,1.,4.]))[0],np.pi/4,
            rtol=4e-15)
        # Local spectrum has no q/Gamma input: only link energy changes under drive.
        driven=m.evaluate(d,p,np.full(len(self.graph.edges),.1))
        np.testing.assert_array_equal(value.electron_energies_J,driven.electron_energies_J)
        self.assertGreater(driven.gradient_energy_J,0.)

    def test_radial_rigidity_matches_dirty_bcs_current_and_origin(self):
        m=self.model
        amplitudes=np.array([0.,.01,.3,1.6])
        rho=amplitudes**2
        kap,kp,kpp=m.rigidity(rho)
        K=kap*self.N0*self.ell**2
        expected=np.empty_like(amplitudes)
        expected[0]=np.pi*self.N0*hbar*self.D/(8*kB*.9)
        expected[1:]=np.pi*self.N0*hbar*self.D/(4*self.Eu*amplitudes[1:])*np.tanh(
            self.Eu*amplitudes[1:]/(2*kB*.9))
        np.testing.assert_allclose(K,expected,rtol=3e-15)
        for r in (.0001,.01,.7,3.):
            step=r*2e-5
            kplus=m.rigidity(np.array(r+step))
            kminus=m.rigidity(np.array(r-step))
            center=m.rigidity(np.array(r))
            self.assertAlmostEqual(float((kplus[0]-kminus[0])/(2*step)/center[1]),1.,
                                   delta=2e-7)
            self.assertAlmostEqual(float((kplus[1]-kminus[1])/(2*step)/center[2]),1.,
                                   delta=2e-6)

    def test_fresh_bulk_21uA_and_tangent(self):
        m=self.model
        reference=m.bulk(temperature_K=.9,current_A=21.297206e-6,width_m=80e-9)
        self.assertTrue(reference.longitudinal_stable)
        self.assertLess(reference.stationarity_residual_bar,2e-12)
        self.assertGreater(reference.depairing_current_A,reference.current_A)
        self.assertAlmostEqual(reference.current_A/21.297206e-6,1.,delta=1e-10)
        # Independent implicit q continuation checks the declared current tangent.
        dq=reference.q_bar*1e-5
        ip,ap=m.bulk_current(.9,reference.q_bar+dq,80e-9)
        im,am=m.bulk_current(.9,reference.q_bar-dq,80e-9)
        measured_dq_dI=2*dq/(ip-im)/m.ell0_m
        measured_dd_dI=(ap-am)/(ip-im)*m.Eu_J
        self.assertAlmostEqual(measured_dq_dI/reference.q_derivative_m_inv_per_A,1.,delta=2e-7)
        self.assertAlmostEqual(measured_dd_dI/reference.gap_derivative_J_per_A,1.,delta=2e-7)


    def test_sparse_jacobians_match_variation_and_hessian_symmetry(self):
        m=self.model
        hessian,current_jac=m.hessian_and_current_jacobian(self.d,self.p,self.alpha)
        direct=m.direction(self.d,self.p,self.dd,alpha=self.alpha)
        vector=np.column_stack((self.dd.real,self.dd.imag)).ravel()
        measured=hessian@vector
        expected=np.column_stack((direct.gap_force_direction_J_per_gapbar.real,
                                  direct.gap_force_direction_J_per_gapbar.imag)).ravel()
        np.testing.assert_allclose(measured,expected,rtol=3e-14,atol=1e-33)
        np.testing.assert_allclose(current_jac@vector,direct.edge_current_direction_A,
                                   rtol=4e-14,atol=1e-20)
        self.assertLess(np.max(abs((hessian-hessian.T).data),initial=0.),
                        np.max(abs(hessian.data))*3e-15)


    def test_discrete_dc_matches_link_force_and_current(self):
        for spacing in (5e-9,3.5e-9):
            nx,ny=7,17
            graph=rectangular_graph(nx,ny,(nx-1)*spacing/self.ell,80e-9/self.ell)
            m=EffectiveBCSFunctional(graph,energy_scale_J=self.Eu,
                N0_per_J_m3=self.N0,ell0_m=self.ell,thickness_m=7e-9,
                count_nodes=self.xi,count_weights=self.w,gradient_law="bcs_radial",
                bath_temperature_K=.9,diffusion_m2_s=self.D)
            ref=m.bulk(temperature_K=.9,current_A=21.297206e-6,width_m=80e-9,
                       phase_step_m=spacing)
            x=graph.coordinates_bar[:,0]
            d=ref.amplitude_bar*np.exp(1j*ref.q_bar*x)
            value=m.evaluate(d,m.fermi_dirac(d,.9))
            interior=(x>x.min()+1e-12)&(x<x.max()-1e-12)
            force=value.gap_force_J_per_gapbar/m.node_energy_scale_J
            self.assertLess(np.max(abs(force[interior])),5e-13)
            left=x==x.min()
            self.assertAlmostEqual(np.sum(value.current_divergence_A[left])/ref.current_A,
                                   1.,delta=5e-14)
            dq=ref.q_bar*1e-5
            ip,ap=m.bulk_current(.9,ref.q_bar+dq,80e-9,spacing)
            im,am=m.bulk_current(.9,ref.q_bar-dq,80e-9,spacing)
            self.assertAlmostEqual((2*dq/(ip-im)/self.ell)/ref.q_derivative_m_inv_per_A,
                                   1.,delta=3e-7)

    def test_invalid_inputs_fail_without_clipping(self):
        with self.assertRaises(ValueError):
            self.model.evaluate(self.d,-self.p)
        with self.assertRaises(ValueError):
            self.model.evaluate(self.d,self.p+1)
        with self.assertRaises(ValueError):
            self.model.bulk(temperature_K=.9,current_A=.001,width_m=80e-9)
        with self.assertRaises(ValueError):
            self.model.rigidity(np.array(-1.))


if __name__=="__main__":
    unittest.main()
