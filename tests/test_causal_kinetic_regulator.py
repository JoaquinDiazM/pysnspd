import unittest
import numpy as np
from pysnspd.experimental import adiabatic_matrix_usadel as ad
from pysnspd.experimental import causal_kinetic_regulator as reg
from pysnspd.experimental import thermal_spatial_usadel as thermal
from pysnspd.experimental import matrix_kinetic_rhs as old
from tests.test_matrix_kinetic_rhs import uniform_jet

class CausalRegulatorTest(unittest.TestCase):
    def setUp(self):
        self.g=thermal.rectangular_graph(4,3,2.,1.5);self.n=self.g.n_nodes
        self.rng=np.random.default_rng(81026)
    def test_full_rak_algebra_with_independent_noncommuting_jets(self):
        def jet():return ad.Jet(*(self.rng.normal(size=(self.n,2,2))+1j*self.rng.normal(size=(self.n,2,2)) for _ in range(4)))
        R,H=jet(),jet();A=ad.advanced(R);K=ad.keldysh(R,H,.31)
        SR,SA,SK=reg.self_energy_triad(H,.027,kappa_s=.31)
        direct=ad.star(SR,K,.31)-ad.star(K,SA,.31)+ad.star(SK,A,.31)-ad.star(R,SK,.31)
        reduced=reg.regulator_term(R,H,.027,kappa_s=.31)
        for field in ('value','energy','time','first'):
            np.testing.assert_allclose(getattr(direct,field),getattr(reduced,field),atol=2e-15)
    def test_static_thermal_exact_eta_identity(self):
        R=uniform_jet(self.n,1.2,.7,.005);H=ad.distribution(np.full(self.n,.4),np.zeros(self.n))
        B=ad.generator(np.full(self.n,1.2),.7)
        oldr,_,_=ad.kinetic_residual(self.g,R,H,B,kappa=.4)
        reduced,_,_,correction=reg.kinetic_residual(self.g,R,H,B,numerical_eta=.005,kappa_s=.4)
        self.assertGreater(np.max(abs(oldr.value)),.001)
        np.testing.assert_allclose(reduced.value,0,atol=4e-16)
        np.testing.assert_allclose(ad.projections(correction.value),0,atol=1e-16)
    def test_storage_matches_full_first_gradient_direction(self):
        R=uniform_jet(self.n,1.2,.7,.005,.07);B=ad.generator(np.full(self.n,1.2),.7,d_time=.07)
        alpha=self.rng.normal(size=len(self.g.edges))*.1
        M=reg.assemble_storage(self.g,R,B,numerical_eta=.005,kappa_s=.4,alpha=alpha)
        h=self.rng.normal(size=(self.n,2))*.1;hE=self.rng.normal(size=(self.n,2))*.1;direction=self.rng.normal(size=(self.n,2))
        def residual(v):
            H=ad.distribution(h[:,0],h[:,1],energy_L=hE[:,0],energy_T=hE[:,1],time_L=v[:,0],time_T=v[:,1])
            r,_,_,_=reg.kinetic_residual(self.g,R,H,B,numerical_eta=.005,kappa_s=.4,alpha=alpha)
            return ad.projections(r.value+r.first)
        np.testing.assert_allclose((M@direction.ravel()).reshape(self.n,2),residual(direction)-residual(np.zeros_like(h)),atol=1e-15)
        np.testing.assert_allclose((M-old.assemble_storage(self.g,R,B,kappa_s=.4,alpha=alpha)).toarray(),0,atol=1e-15)
    def test_zero_eta_reproduces_original_rhs_and_no_physical_source(self):
        R=uniform_jet(self.n,1.2,2.,0.,.07);B=ad.generator(np.full(self.n,1.2),2.,d_time=.07)
        h=np.column_stack((np.full(self.n,.4),np.zeros(self.n)));hE=h*.325
        args=dict(kappa_s=.4,projected_source=np.zeros_like(h),source_contract='explicit collisionless test')
        a=old.solve_distribution_rhs(self.g,R,h,hE,B,**args)
        b=reg.solve_distribution_rhs(self.g,R,h,hE,B,numerical_eta=0.,**args)
        np.testing.assert_allclose(a.distribution_rhs_per_s,b.distribution_rhs_per_s,atol=1e-15)
        np.testing.assert_array_equal(b.projected_source,0.)
    def test_nonequilibrium_both_modes_with_time_dependent_phi_has_zero_trace_regulator(self):
        # Independent R entries/derivatives are deliberately broad: no BCS
        # relation or solved kinetic equation can make this identity pass.
        def random():return self.rng.normal(size=(self.n,2,2))+1j*self.rng.normal(size=(self.n,2,2))
        R=ad.Jet(random(),random(),random(),random())
        h=self.rng.normal(size=(self.n,2))*.1;hE=self.rng.normal(size=(self.n,2));hdot=self.rng.normal(size=(self.n,2))
        H=ad.distribution(h[:,0],h[:,1],energy_L=hE[:,0],energy_T=hE[:,1],time_L=hdot[:,0],time_T=hdot[:,1])
        B=ad.generator(np.full(self.n,1.2),.7,d_time=.09,
            voltage=self.rng.normal(size=self.n)*.02,voltage_time=self.rng.normal(size=self.n)*.3)
        raw,_,_=ad.kinetic_residual(self.g,R,H,B,kappa=.4)
        corrected,_,_,sigma=reg.kinetic_residual(self.g,R,H,B,numerical_eta=.005,kappa_s=.4)
        np.testing.assert_allclose(ad.projections(sigma.value),0,atol=1e-15)
        np.testing.assert_allclose(ad.projections(sigma.first),0,atol=1e-15)
        np.testing.assert_allclose(ad.projections(corrected.value+corrected.first),ad.projections(raw.value+raw.first),atol=1e-15)

    def test_normalized_triangular_identity_independent_of_spectral_solution(self):
        R=uniform_jet(self.n,1.2,.7,.05,.07)
        R=ad.Jet.of(R.value,R.energy,R.time,ad.normalization_particular(R,.4))
        h=self.rng.normal(size=(self.n,2))*.1;hE=self.rng.normal(size=(self.n,2));hdot=self.rng.normal(size=(self.n,2))
        H=ad.distribution(h[:,0],h[:,1],energy_L=hE[:,0],energy_T=hE[:,1],time_L=hdot[:,0],time_T=hdot[:,1])
        B=ad.generator(np.full(self.n,1.2),.7,d_time=.07,
            voltage=self.rng.normal(size=self.n)*.02,voltage_time=self.rng.normal(size=self.n)*.3)
        Bs=ad.Jet.of(B.value+.05*ad.TAU3,B.energy,B.time,B.first)
        alpha=self.rng.normal(size=len(self.g.edges))*.2
        LK,_,K,_=reg.kinetic_residual(self.g,R,H,B,numerical_eta=.05,kappa_s=.4,alpha=alpha)
        LR,NR=ad.spectral_residual(self.g,R,Bs,alpha=alpha,kappa=.4);LA=ad.advanced(LR);A=ad.advanced(R)
        lhs=ad.star(R,LK,.4)+ad.star(LK,A,.4)
        rhs=(ad.star(LR,K,.4)+ad.star(K,LA,.4)).scale(-1)
        np.testing.assert_allclose(NR.value,0,atol=2e-15)
        np.testing.assert_allclose(NR.first,0,atol=2e-15)
        self.assertGreater(np.max(abs(LR.value)),.1)
        np.testing.assert_allclose(lhs.value,rhs.value,atol=2e-14)
        np.testing.assert_allclose(lhs.first,rhs.first,atol=2e-14)

    def test_guard_rejects_eta_in_physical_B(self):
        R=uniform_jet(self.n,1.2,2.,.03);h=np.column_stack((np.full(self.n,.4),np.zeros(self.n)))
        with self.assertRaisesRegex(ValueError,'eta is forbidden'):
            reg.solve_distribution_rhs(self.g,R,h,np.zeros_like(h),ad.generator(np.full(self.n,1.2),2.,eta=.03),
                numerical_eta=.03,kappa_s=.4,projected_source=np.zeros_like(h),source_contract='explicit')

if __name__=='__main__':unittest.main()
