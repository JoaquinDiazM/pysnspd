"""Independent normal-metal, energy, static/time-gauge moment controls."""
import unittest
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.constants import elementary_charge
from pysnspd.experimental import adiabatic_matrix_usadel as ad
from pysnspd.experimental.nambu_moments import NambuMoments


class NambuMomentsTests(unittest.TestCase):
    def setUp(self):
        points,weights=leggauss(128)
        self.Eu=2e-22;self.E=2.5*(points+1);self.w=2.5*weights
        self.volume=np.array([2e-24,3e-24]);self.N0=8e46;self.I0=4e-6
        self.obs=NambuMoments(self.E*self.Eu,self.w*self.Eu,self.volume,self.N0,self.Eu,self.I0)
        self.shape=(len(points),2,2,2);self.zeros=np.zeros(self.shape,complex)
        self.flux=np.zeros((len(points),1,2,2),complex)
        self.T=.2;self.h0=np.tanh(self.E/(2*self.T))
        self.reference=np.broadcast_to(2*self.h0[:,None,None,None]*ad.TAU3,(len(points),2,2,2)).copy()

    def evaluate(self,K0,K1=None,phi=None,thermal=None,flux0=None,flux1=None):
        return self.obs.evaluate(K_zero=K0,K_first=self.zeros if K1 is None else K1,
            K_thermal_zero=self.reference,phi_V=np.zeros(2) if phi is None else phi,
            thermal_gradient_J_per_gapbar=np.zeros(2) if thermal is None else thermal,
            edge_flux_zero=self.flux if flux0 is None else flux0,
            edge_flux_first=self.flux if flux1 is None else flux1,
            reference_contract='normal thermal reference at same instantaneous spectrum and gauge')

    def test_finite_normal_reservoir_neutrality_not_only_small_voltage(self):
        # Integral of hT(E) from zero to infinity is -v for every finite v,
        # by oddness of tanh and a shift of the two integration intervals.
        v=np.array([.1,-.4])
        minus=np.tanh((self.E[:,None]-v)/(2*self.T))
        plus=np.tanh((self.E[:,None]+v)/(2*self.T))
        hL=(minus+plus)/2;hT=(minus-plus)/2
        K=2*(hL[:,:,None,None]*ad.TAU3+hT[:,:,None,None]*ad.I2)
        r=self.evaluate(K,phi=v*self.Eu/elementary_charge)
        logcosh=lambda z:np.logaddexp(z,-z)-np.log(2.)
        finite_integral=self.T*(logcosh((5-v)/(2*self.T))-logcosh((5+v)/(2*self.T)))
        np.testing.assert_allclose(r.charge_moment_zero_J/self.Eu,finite_integral,rtol=2e-13,atol=1e-15)
        # The omitted thermal tail is known analytically and is reported, not repaired.
        expected=(v+finite_integral)*self.Eu/elementary_charge
        np.testing.assert_allclose(r.neutrality_residual_V,expected,atol=2e-14*self.Eu/elementary_charge)
        self.assertGreater(np.linalg.norm(r.charge_moment_zero_J),.1*self.Eu)

    def test_normal_hot_electrons_energy_has_spin_and_branch_count(self):
        excess=np.array([.02,.03])
        K=self.reference-4*excess[None,:,None,None]*ad.TAU3
        r=self.evaluate(K)
        # Four N0 states per positive energy: hL changes by minus two times f.
        expected=4*self.N0*sum(self.volume*excess)*self.Eu**2*(5**2/2)
        self.assertAlmostEqual(r.energy_increment_zero_J/expected,1.,places=14)
        self.assertEqual(r.energy_increment_first_J,0.)
        np.testing.assert_allclose(r.neutrality_residual_V,0.,atol=0.)

    def test_anomalous_first_moment_is_retained_once(self):
        K1=self.zeros.copy();K1[:,:,0,1]=(.1+.2j)
        r=self.evaluate(self.reference,K1)
        expected=1j*self.N0*self.Eu*self.volume*(.1+.2j)*5*self.Eu
        np.testing.assert_allclose(r.gap_force_J_per_gapbar,expected,rtol=2e-15,atol=0.)
        np.testing.assert_allclose(r.gap_force_increment_zero_J_per_gapbar,0.,atol=0.)
        self.assertGreater(np.linalg.norm(r.gap_force_increment_first_J_per_gapbar),0.)

    def test_static_gauge_rotates_pair_force_preserves_charge_energy(self):
        K=self.reference.copy();K[:,:,0,1]=.07+.09j;K[:,:,1,0]=-.11j
        r=self.evaluate(K)
        chi=np.array([.2,-.5]);U=ad.links(-chi).value
        gauged=U[None]@K@U.conj().transpose(0,2,1)[None]
        q=self.evaluate(gauged)
        np.testing.assert_allclose(q.gap_force_J_per_gapbar,np.exp(1j*chi)*r.gap_force_J_per_gapbar,rtol=3e-15)
        np.testing.assert_allclose(q.neutrality_residual_V,r.neutrality_residual_V,atol=1e-20)
        self.assertLess(abs(q.energy_increment_zero_J-r.energy_increment_zero_J),
            3e-15*self.N0*self.Eu**2*sum(self.volume))

    def test_normal_temporal_gauge_charge_compensates_voltage(self):
        kappa=.07;chi_t=np.array([.13,-.04]);chi=np.array([.2,-.5])
        hE=(1-self.h0**2)/(2*self.T)
        K=ad.Jet.of(self.reference,energy=np.broadcast_to(2*hE[:,None,None,None]*ad.TAU3,self.shape))
        U=ad.links(-chi,-chi_t)
        moved=ad.transport(U,K,kappa)
        r=self.evaluate(moved.value,moved.first,phi=-kappa*chi_t*self.Eu/elementary_charge)
        np.testing.assert_allclose(r.neutrality_residual_V,0.,atol=1e-12*self.Eu/elementary_charge)
        self.assertGreater(np.linalg.norm(r.charge_moment_first_J),0.)
        without=self.evaluate(moved.value,phi=-kappa*chi_t*self.Eu/elementary_charge)
        self.assertGreater(np.linalg.norm(without.neutrality_residual_V),1e-7)

    def test_full_charge_flux_and_first_correction_in_amperes(self):
        flux0=self.flux.copy();flux1=self.flux.copy()
        flux0[:,:,0,0]=.2;flux0[:,:,1,1]=-.2
        flux1[:,:,0,0]=.07;flux1[:,:,1,1]=-.07
        r=self.evaluate(self.reference,flux0=flux0,flux1=flux1)
        np.testing.assert_allclose(r.edge_current_A,[self.I0*.27*5],rtol=2e-15)
        self.assertGreater(r.edge_current_first_A[0],0.)

    def test_reference_contract_and_finite_measures_are_required(self):
        with self.assertRaisesRegex(ValueError,'positive'):
            NambuMoments(self.E*self.Eu,-self.w*self.Eu,self.volume,self.N0,self.Eu,self.I0)
        with self.assertRaisesRegex(ValueError,'real'):
            NambuMoments(self.E*self.Eu,self.w*self.Eu,self.volume.astype(complex),self.N0,self.Eu,self.I0)
        with self.assertRaisesRegex(ValueError,'reference'):
            self.obs.evaluate(K_zero=self.reference,K_first=self.zeros,K_thermal_zero=self.reference,
                phi_V=np.zeros(2),thermal_gradient_J_per_gapbar=np.zeros(2),
                edge_flux_zero=self.flux,edge_flux_first=self.flux,reference_contract='')


    def test_stable_current_preserves_nonzero_thermal_baseline_and_reports_defect(self):
        baseline=np.array([21e-6]);thermalflux=self.flux.copy();thermalflux[:,:,0,0]=.2
        actual=thermalflux.copy();actual[:,:,0,0]+=.03
        first=self.flux.copy();first[:,:,0,0]=.01
        completed=self.obs.complete_current(thermal_edge_current_A=baseline,edge_flux_zero=actual,
            thermal_edge_flux_zero=thermalflux,edge_flux_first=first,
            reference_contract='same instantaneous spectrum; independent Matsubara thermal current with its tail')
        np.testing.assert_allclose(completed.edge_current_A,baseline+self.I0*.04/2*5,rtol=2e-15)
        # A baseline restored by subtraction must not hide its real-axis discrepancy.
        np.testing.assert_allclose(completed.thermal_quadrature_discrepancy_A,
            [self.I0*.2/2*5-21e-6],rtol=2e-15)
        dc=self.obs.complete_current(thermal_edge_current_A=baseline,edge_flux_zero=thermalflux,
            thermal_edge_flux_zero=thermalflux,edge_flux_first=self.flux,reference_contract='same-gap DC reference')
        np.testing.assert_allclose(dc.edge_current_A,baseline,atol=0.)
        self.assertGreater(abs(dc.thermal_quadrature_discrepancy_A[0]),1e-6)
