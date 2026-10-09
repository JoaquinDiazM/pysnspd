"""Coupled physical operators, SI circuit feedback and honest rank diagnostics."""
import unittest
import numpy as np
from scipy.constants import hbar, elementary_charge as e
from pysnspd.experimental import adiabatic_matrix_usadel as ad
from pysnspd.experimental import retarded_spatial_usadel as ret
from pysnspd.experimental.thermal_spatial_usadel import rectangular_graph
from pysnspd.experimental.nambu_moments import NambuMoments
from pysnspd.experimental.joint_nambu_rows import JointNambuRows
from pysnspd.experimental.electrical_ports import ThesisCircuitParameters
from pysnspd.experimental.joint_evolution import (
    BulkCutLaw, JointSpectralSlice, JointSnapshotSystem, solve_affine_candidate)


def physical_small_system():
    g=rectangular_graph(3,2,2.,1.)
    n=g.n_nodes;left=np.array([0,1]);right=np.array([4,5]);fixed=np.r_[left,right]
    Eu=1.2e-22;d=np.full(n,1.2,dtype=complex);alpha=np.zeros(len(g.edges))
    energies=np.array([.3,.8,1.6,2.4])*Eu;weights=np.full(4,.4*Eu)
    bulk=BulkCutLaw(1.2,0.,0.,0.,"rigid-amplitude test bulk, q=0")
    slices=[];reference=[];thermal_flux=[]
    for energy in energies:
        a,b=ret.uniform_fields(d,.03-1j*energy/Eu)
        f,ft,rho=ret.fields(a,b);R0=ret.nambu_matrix(f,ft,rho)
        hl=np.full(n,np.tanh(energy/(2*.25*Eu)))
        h=np.column_stack((hl,np.zeros(n)))
        hE=np.column_stack(((1-hl**2)/(.5),np.zeros(n)))
        sl=JointSpectralSlice(g,d,alpha,R0,energy,.03*Eu,Eu,h,hE,
            fixed,bulk,projected_source=np.zeros((n,2)),
            source_contract="collisionless physical operator test; not a detector")
        slices.append(sl)
        _,flux,K=ad.kinetic_residual(g,ad.Jet.of(R0),ad.distribution(hl,np.zeros(n)),
            ad.generator(d,energy/Eu),alpha=alpha,kappa=hbar/(2*Eu))
        reference.append(K.value);thermal_flux.append(flux.value)
    moments=NambuMoments(energies,weights,np.full(n,3e-25),9e46,Eu,Eu/(2*e*608.))
    rows=JointNambuRows(moments,g.edges,left_nodes=left,right_nodes=right)
    circuit=ThesisCircuitParameters(Lk_ext_H=1.3e-7,V_bias_V=0.)
    return JointSnapshotSystem(g,d,alpha,slices,rows,K_thermal_zero=reference,
        thermal_gradient_J_per_gapbar=np.zeros(n,complex),
        thermal_edge_current_A=np.zeros(len(g.edges)),
        thermal_edge_flux_zero=thermal_flux,circuit=circuit,circuit_state=[0.,0.,0.],
        bulk_law=bulk,reference_contract="same instantaneous BCS gap and links")


class JointEvolutionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.system=physical_small_system()

    def test_zero_drive_is_one_stationary_system(self):
        s=self.system;v=s.evaluate(np.zeros(s.size))
        np.testing.assert_allclose(v.packed_rows,0.,atol=1e-12)
        self.assertLess(v.kinetic_residual_max,1e-12)
        self.assertLess(v.full_nambu_residual_max_per_mass,1e-10)
        np.testing.assert_array_equal(v.circuit_rate_SI_per_s,[0.,0.,0.])

    def test_all_fields_respond_in_the_same_affine_operator(self):
        s=self.system;rng=np.random.default_rng(510)
        a=rng.normal(size=s.size)*.03;b=rng.normal(size=s.size)*.02
        zero=s.evaluate(np.zeros(s.size));va=s.evaluate(a);vb=s.evaluate(b)
        both=s.evaluate(a+b)
        np.testing.assert_allclose(both.packed_rows,
            va.packed_rows+vb.packed_rows-zero.packed_rows,rtol=2e-10,atol=2e-11)
        self.assertGreater(np.linalg.norm(va.distribution_rate_per_s),1e9)
        self.assertGreater(np.linalg.norm(va.moments.gap_force_J_per_gapbar),1e-25)
        self.assertGreater(np.linalg.norm(va.moments.edge_current_A),1e-12)
        self.assertGreater(np.linalg.norm(va.circuit_rate_SI_per_s),1.)

    def test_ports_share_voltage_and_current_with_cm_and_josephson(self):
        s=self.system;x=np.zeros(s.size);x[-2]=.03;x[-1]=-.01
        v=s.evaluate(x);expected=s.circuit.rhs(s.y,.04*s.Vunit)
        np.testing.assert_allclose(v.circuit_rate_SI_per_s,expected,rtol=1e-14)
        np.testing.assert_allclose(v.gap_rate_per_s[s.fixed],
            -2j*e*v.phi_V[s.fixed]*s.gap[s.fixed]/hbar,rtol=1e-14)
        self.assertEqual(v.moments.left_equipotential_residual_V.max(),0.)
        self.assertEqual(v.moments.right_equipotential_residual_V.max(),0.)
        balance=s.circuit.power_balance(s.y,v.moments.port_voltage_V,
            derivative=v.circuit_rate_SI_per_s)
        self.assertLess(abs(balance.residual_W),1e-25)

    def test_assembled_matrix_predicts_independent_recomputed_residual(self):
        s=self.system;A,b=s.affine_matrix()
        x=np.random.default_rng(91).normal(size=s.size)*.02
        np.testing.assert_allclose(A@x+b,s.evaluate(x).packed_rows,
            rtol=2e-10,atol=2e-11)
        _,diagnostic=solve_affine_candidate(A,b)
        self.assertFalse(diagnostic["temporal_admission"])
        self.assertFalse(diagnostic["energy_admission"])
        self.assertTrue(diagnostic["no_regularization"])

    def test_common_voltage_is_an_exact_gauge_and_not_a_physical_mode(self):
        s=self.system;rng=np.random.default_rng(77);x=rng.normal(size=s.size)*.01
        one=s.evaluate(x);two=s.evaluate(x+.17*s.gauge_direction)
        np.testing.assert_allclose(one.packed_rows,two.packed_rows,rtol=1e-10,atol=1e-12)
        np.testing.assert_allclose(one.distribution_rate_per_s,two.distribution_rate_per_s,
            rtol=1e-10,atol=1e-3)
        np.testing.assert_allclose(one.circuit_rate_SI_per_s,two.circuit_rate_SI_per_s,
            rtol=1e-12,atol=1e-5)
        A,b=s.affine_matrix()
        np.testing.assert_allclose(A@s.gauge_direction,0.,atol=1e-11)
        _,info=solve_affine_candidate(A,b,gauge_direction=s.gauge_direction,
            gauge_reference_index=s.size-1)
        self.assertTrue(info["gauge_reduced"])
        self.assertTrue(info["no_physical_rows_removed"])
        self.assertEqual(info["physical_columns"],s.size-1)

    def test_no_extra_null_direction_is_repaired(self):
        A=np.diag([2.,1.,0.]);x,info=solve_affine_candidate(A,[1.,1.,.2])
        self.assertEqual(info["rank"],2)
        self.assertFalse(info["algebraically_solved_full_rank"])
        self.assertAlmostEqual(info["scaled_residual_max"],.2)
        np.testing.assert_allclose(x,[-.5,-1.,0.])

    def test_full_rank_candidate_is_rechecked_and_does_not_admit_time(self):
        A=np.array([[2.,1.],[.5,3.]])
        x,info=solve_affine_candidate(A,[.7,-.2])
        np.testing.assert_allclose(A@x+[.7,-.2],0.,atol=1e-15)
        self.assertTrue(info["algebraically_solved_full_rank"])
        self.assertFalse(info["temporal_admission"])


if __name__=="__main__":
    unittest.main()
