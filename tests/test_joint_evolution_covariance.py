"""Independent Wilson-time gauge and passive-current controls for stage 5.1.

These tests derive the current from Nambu Jet products. No additive Ohmic
constitutive term, detector trajectory or stage 5.2/5.3 admission is involved.
"""
import unittest
from dataclasses import replace
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.constants import elementary_charge as e, hbar, Boltzmann
from pysnspd.experimental import adiabatic_matrix_usadel as ad
from pysnspd.experimental import causal_kinetic_regulator as regulator
from pysnspd.experimental import thermal_spatial_usadel as thermal
from pysnspd.experimental import bulk_current_reference as bulk
from pysnspd.experimental.nambu_moments import NambuMoments
from pysnspd.experimental.electrical_ports import ThesisCircuitParameters


class JointTemporalGaugeTests(unittest.TestCase):
    def setUp(self):
        self.graph = thermal.rectangular_graph(4, 3, 2., 1.5)
        self.n = self.graph.n_nodes
        self.Eu = 1.2e-22
        self.kappa = hbar/(2*self.Eu)
        self.rng = np.random.default_rng(510091)

    def test_local_temporal_gauge_transports_full_R_H_B_and_flux(self):
        g = self.graph
        n = self.n
        d = .8+.03*self.rng.normal(size=n)
        dt = .03e12*self.rng.normal(size=n)
        E, eta = .6, .09
        z = eta-1j*E
        w = np.sqrt(z*z+d*d)
        tau1 = np.array([[0., 1.], [1., 0.]])
        R0 = (z*ad.TAU3+d[:, None, None]*tau1)/w[:, None, None]
        RE = -1j*ad.TAU3/w[:, None, None]+1j*z*R0/(w*w)[:, None, None]
        Rt = dt[:, None, None]*(tau1/w[:, None, None]-
                                      d[:, None, None]*R0/(w*w)[:, None, None])
        R = ad.Jet.of(R0, RE, Rt)
        H = ad.distribution(.4+.02*self.rng.normal(size=n),
            .01*self.rng.normal(size=n),energy_L=.2+.01*self.rng.normal(size=n),
            energy_T=.01*self.rng.normal(size=n),
            time_L=.01e12*self.rng.normal(size=n),
            time_T=.01e12*self.rng.normal(size=n))
        B = ad.generator(d, E, d_time=dt)
        alpha = .03*self.rng.normal(size=len(g.edges))
        alpha_time = .02e12*self.rng.normal(size=len(g.edges))
        chi = .2*self.rng.normal(size=n)
        chi_time = .1e12*self.rng.normal(size=n)
        U = ad.links(-chi, -chi_time)
        RG, HG, BG = [ad.transport(U, field, self.kappa) for field in (R,H,B)]
        tail, head = g.edges.T
        aG = alpha+chi[head]-chi[tail]
        atG = alpha_time+chi_time[head]-chi_time[tail]
        original, flux, K, correction = regulator.kinetic_residual(g,R,H,B,
            numerical_eta=eta,kappa_s=self.kappa,alpha=alpha,alpha_time_per_s=alpha_time)
        transformed, transformed_flux, KG, transformed_correction = regulator.kinetic_residual(
            g,RG,HG,BG,numerical_eta=eta,kappa_s=self.kappa,
            alpha=aG,alpha_time_per_s=atG)
        for actual, expected in (
            (transformed,ad.transport(U,original,self.kappa)),
            (KG,ad.transport(U,K,self.kappa)),
            (transformed_correction,ad.transport(U,correction,self.kappa)),
            (transformed_flux,ad.transport(U.take(tail),flux,self.kappa))):
            np.testing.assert_allclose(actual.value,expected.value,atol=2e-14)
            np.testing.assert_allclose(actual.first,expected.first,atol=2e-14)
        # Hfirst is transformed, even though its ordinary diagonal value is
        # unchanged by U. Erasing this term in both frames breaks covariance.
        h_expected = self.kappa*chi_time[:,None,None]*(ad.TAU3@H.energy)
        np.testing.assert_allclose(HG.first,h_expected,atol=2e-15)
        broken_H = replace(HG,first=np.zeros_like(HG.first))
        broken, _, _, _ = regulator.kinetic_residual(g,RG,broken_H,BG,
            numerical_eta=eta,kappa_s=self.kappa,alpha=aG,alpha_time_per_s=atG)
        self.assertGreater(np.max(abs(broken.first-transformed.first)),1e-5)

    def _normal_fields(self,phi):
        x, w = leggauss(96)
        max_energy_bar, temperature_bar = 12., .23
        eb = .5*max_energy_bar*(x+1)
        weights = .5*max_energy_bar*w*self.Eu
        hL = np.tanh(eb/(2*temperature_bar))
        hE = (1-hL*hL)/(2*temperature_bar)
        tail,head = self.graph.edges.T
        at = 2*e/hbar*(phi[head]-phi[tail])
        R = ad.Jet.of(np.broadcast_to(ad.TAU3,(self.n,2,2)))
        currents, K0, flux0, flux1, Hs = [],[],[],[],[]
        for energy, hl, he in zip(eb,hL,hE):
            H = ad.distribution(np.full(self.n,hl),np.zeros(self.n),
                               energy_L=he)
            B = ad.generator(np.zeros(self.n),energy)
            _,F,K = ad.kinetic_residual(self.graph,R,H,B,
                alpha=np.zeros(len(tail)),alpha_time=at,kappa=self.kappa)
            K0.append(K.value);flux0.append(F.value);flux1.append(F.first)
            Hs.append(H)
        moments = NambuMoments(eb*self.Eu,weights,np.full(self.n,3e-25),
            9e46,self.Eu,self.Eu/(2*e*608.))
        return moments,R,Hs,np.asarray(K0),np.asarray(flux0),np.asarray(flux1),at

    def test_normal_electrical_conductivity_emerges_from_alpha_time(self):
        phi = 8e-5-7e-5*self.graph.coordinates_bar[:,0]
        obs,R,Hs,K0,F0,F1,at = self._normal_fields(phi)
        value = obs.evaluate(K_zero=K0,K_first=np.zeros_like(K0),
            K_thermal_zero=K0,phi_V=np.zeros(self.n),
            thermal_gradient_J_per_gapbar=np.zeros(self.n,complex),
            edge_flux_zero=F0,edge_flux_first=F1,
            reference_contract='normal Wilson-energy frame; no added Ohmic current')
        tail,head = self.graph.edges.T
        drop = phi[tail]-phi[head]
        expected = self.graph.conductance*drop/608.
        np.testing.assert_allclose(value.edge_current_A,expected,rtol=3e-13,atol=1e-21)
        np.testing.assert_array_equal(value.edge_current_zero_A,np.zeros(len(tail)))
        self.assertGreater(np.dot(drop,value.edge_current_A),0.)
        _,_,_,_,_,reversed_flux,_ = self._normal_fields(-phi+1.3e-3)
        # Reverse the physical field and add an arbitrary common voltage.
        np.testing.assert_allclose(reversed_flux,-F1,rtol=2e-14,atol=2e-15)

    def test_local_energy_coordinate_recovers_lab_neutrality_and_same_current(self):
        phi = 8e-5-7e-5*self.graph.coordinates_bar[:,0]
        obs,R,Hs,K0,F0,F1,at = self._normal_fields(phi)
        # lab -> local temporal: chi_dot=+2e phi/hbar. Consequently the
        # inverse frame carries H1_lab=-(e phi/Eu)*tau3*H_E.
        inverse = ad.links(np.zeros(self.n),2*e*phi/hbar)
        lab_K1,lab_F1 = [],[]
        tail,head = self.graph.edges.T
        for energy,H in zip(obs.energy/self.Eu,Hs):
            lab_H = ad.transport(inverse,H,self.kappa)
            expected = -(e*phi/self.Eu)[:,None,None]*(ad.TAU3@H.energy)
            np.testing.assert_allclose(lab_H.first,expected,atol=2e-16)
            # The scalar generator correction is the inverse gauge
            # transform of temporal B. It commutes ordinarily with R0.
            B = ad.generator(np.zeros(self.n),energy)
            lab_B = ad.transport(inverse,B,self.kappa)
            lab_R = ad.transport(inverse,R,self.kappa)
            _,F,K = ad.kinetic_residual(self.graph,lab_R,lab_H,lab_B,
                alpha=np.zeros(len(tail)),alpha_time=np.zeros(len(tail)),kappa=self.kappa)
            lab_K1.append(K.first);lab_F1.append(F.first)
        lab = obs.evaluate(K_zero=K0,K_first=np.asarray(lab_K1),
            K_thermal_zero=K0,phi_V=phi,
            thermal_gradient_J_per_gapbar=np.zeros(self.n,complex),
            edge_flux_zero=F0,edge_flux_first=np.asarray(lab_F1),
            reference_contract='lab energy frame reconstructed from same Wilson fields')
        np.testing.assert_allclose(lab.neutrality_residual_V,0.,atol=3e-17)
        np.testing.assert_allclose(np.asarray(lab_F1),F1,rtol=2e-13,atol=2e-15)
        wrong = obs.evaluate(K_zero=K0,K_first=np.zeros_like(K0),
            K_thermal_zero=K0,phi_V=phi,
            thermal_gradient_J_per_gapbar=np.zeros(self.n,complex),
            edge_flux_zero=F0,edge_flux_first=F1,
            reference_contract='deliberately incorrect mixed gauge, for control')
        np.testing.assert_allclose(wrong.neutrality_residual_V,phi,rtol=1e-14)

    def test_bank_qdot_is_the_same_passive_inductive_field_and_CM_sign(self):
        params = dict(width_m=80e-9,Tc_K=10.,diffusion_m2_s=5e-5,
                      sheet_resistance_ohm=608.)
        reference = bulk.solve_bulk_reference(.12,256,.27)
        tangent = bulk.bulk_current_tangent(reference,**params)
        scale = bulk.physical_scales(Tc_K=params['Tc_K'],
            diffusion_m2_s=params['diffusion_m2_s'],
            sheet_resistance_ohm=params['sheet_resistance_ohm'])
        I = reference.total_current_A(**params)
        circuit = ThesisCircuitParameters(Lk_ext_H=1.3e-7,V_bias_V=1e4*I)
        voltage = 8e-5
        rate = circuit.rhs([I,I,0.],voltage)
        self.assertLess(rate[1],0.)
        qdot = tangent.q_bar_derivative_per_A*rate[1]
        electric_V_per_m = hbar*qdot/(2*e*scale['ell0_m'])
        self.assertLess(electric_V_per_m,0.)
        L_out = 3e-6
        L_wire = L_out*bulk.differential_inductance_per_length_H_m(reference,**params)
        np.testing.assert_allclose(electric_V_per_m*L_out,L_wire*rate[1],rtol=2e-15)
        # With covariant phase rate zero, alpha_t=-qdot*dX in the bank.
        dX = .2
        phi_head_minus_tail = hbar/(2*e)*(-qdot*dX)
        alpha_t = 2*e/hbar*phi_head_minus_tail
        np.testing.assert_allclose(alpha_t,-qdot*dX,rtol=2e-15)
        balance = circuit.power_balance([I,I,0.],voltage,derivative=rate)
        self.assertLess(abs(balance.residual_W),1e-22)


if __name__ == '__main__':
    unittest.main()
