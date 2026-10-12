"""Physical circuit/power/Noether controls for the effective neutral ports."""
import unittest
import numpy as np
from scipy.constants import elementary_charge as e, hbar
from pysnspd.experimental.thermal_spatial_usadel import rectangular_graph
from pysnspd.experimental.electrical_ports import ThesisCircuitParameters
from pysnspd.experimental.effective_bcs_ports import (
    EffectiveNeutralPorts, freeze_exterior_inductance,
    bulk_boundary_load, boundary_material_power_W)


class EffectivePortsTests(unittest.TestCase):
    def setUp(self):
        self.graph=rectangular_graph(7,5,6.,2.)
        x=self.graph.coordinates_bar[:,0]
        self.left=np.flatnonzero(x==x.min())
        self.right=np.flatnonzero(x==x.max())
        self.ports=EffectiveNeutralPorts(self.graph,608.,self.left,self.right)
        self.current=21e-6

    def test_normal_strip_exact_resistance_complete_cuts_and_power(self):
        v=self.ports.solve(np.zeros(len(self.graph.edges)),self.current)
        self.assertAlmostEqual(v.Vdev_V/(608.*3.*self.current),1.,places=13)
        np.testing.assert_allclose(v.phi_V[self.left],v.Vdev_V,atol=0.)
        np.testing.assert_array_equal(v.phi_V[self.right],np.zeros(len(self.right)))
        self.assertLess(max(abs(v.interior_current_residual_A)),1e-18)
        self.assertLess(max(abs(v.port_current_residual_A)),1e-18)
        self.assertGreater(v.normal_joule_W,0.)
        self.assertTrue(np.all(v.nodal_joule_W>=0.))
        self.assertLess(abs(v.normal_joule_W-self.current*v.Vdev_V),1e-18)
        self.assertLess(abs(v.power_residual_W),1e-18)
        self.assertLess(abs(v.boundary_current_power_residual_W),1e-18)
        self.assertEqual(v.port_power_W,self.current*v.Vdev_V)
        self.assertLess(abs(v.decomposition_residual_W),1e-18)
        self.assertLess(abs(v.joule_allocation_residual_W),1e-18)
        # Vertex-dual half width at the two lateral corners.
        w=np.array([.5,1.,1.,1.,.5]);w/=sum(w)
        np.testing.assert_allclose(v.nodal_boundary_injection_A[self.left],w*self.current,rtol=2e-13,atol=1e-20)

    def test_forward_uniform_supercurrent_is_dc_with_zero_voltage_and_joule(self):
        tail,head=self.graph.edges.T
        dx=self.graph.coordinates_bar[head,0]-self.graph.coordinates_bar[tail,0]
        js=self.graph.conductance*dx*self.current/2.
        v=self.ports.solve(js,self.current)
        self.assertLess(max(abs(v.phi_V)),1e-16)
        self.assertLess(v.normal_joule_W,1e-28)
        cm=ThesisCircuitParameters(Lk_ext_H=105e-9,V_bias_V=1e4*self.current)
        c=self.ports.coupled(js,cm,[self.current,self.current,0.])
        self.assertLess(max(abs(c.circuit_rate_SI_per_s)),1e-8)
        self.assertEqual(c.baseline_subtracted_Vout_V,0.)
        self.assertLess(abs(c.circuit_power.residual_W),1e-18)

    def test_gauge_reference_does_not_change_drops_or_joule(self):
        rng=np.random.default_rng(61)
        js=rng.normal(size=len(self.graph.edges))*1e-7
        a=self.ports.solve(js,self.current)
        b=self.ports.solve(js,self.current,reference_potential_V=17.)
        np.testing.assert_array_equal(a.delta_phi_V,b.delta_phi_V)
        np.testing.assert_array_equal(a.I_total_A,b.I_total_A)
        self.assertEqual(a.normal_joule_W,b.normal_joule_W)
        self.assertEqual(a.boundary_power_W,b.boundary_power_W)
        np.testing.assert_allclose(b.phi_V-a.phi_V,17.,atol=1e-14)

    def test_same_energy_gauge_current_has_exact_Ward_and_phase_power(self):
        rng=np.random.default_rng(71);g=self.graph
        d=1.+.2*rng.normal(size=g.n_nodes)+.3j*rng.normal(size=g.n_nodes)
        alpha=.2*rng.normal(size=len(g.edges));P=np.exp(-1j*alpha)
        tail,head=g.edges.T;C=2.4e-22
        difference=d[tail]-P*d[head]
        gradient=np.zeros(g.n_nodes,complex)
        np.add.at(gradient,tail,2*C*g.conductance*difference)
        np.add.at(gradient,head,2*C*g.conductance*(d[head]-P.conj()*d[tail]))
        js=4*e/hbar*C*g.conductance*np.imag(d[tail].conj()*P*d[head])
        v=self.ports.solve(js,self.current)
        check=self.ports.gauge_power_check(d,gradient,js,v.phi_V)
        self.assertLess(check.noether_residual_max_J,1e-36)
        self.assertLess(abs(check.power_residual_W),1e-18)
        self.assertLess(abs(check.edge_superconducting_power_W-v.super_power_W),1e-18)
        self.assertLess(abs(v.boundary_power_W-v.normal_joule_W-check.field_phase_power_W),1e-18)

    def test_shared_radial_BCS_functional_current_and_port_power(self):
        from scipy.constants import Boltzmann as kB
        from pysnspd.experimental.effective_bcs_functional import EffectiveBCSFunctional,bcs_count_quadrature
        Eu=kB*8.65;D=5e-5;ell=np.sqrt(hbar*D/(2*Eu))
        xi,w=bcs_count_quadrature()
        f=EffectiveBCSFunctional(self.graph,energy_scale_J=Eu,N0_per_J_m3=9.153297821e46,
            ell0_m=ell,thickness_m=7e-9,count_nodes=xi,count_weights=w,
            gradient_law="bcs_radial",bath_temperature_K=.9,diffusion_m2_s=D)
        rng=np.random.default_rng(728)
        d=(.4+rng.random(self.graph.n_nodes))*np.exp(.2j*rng.normal(size=self.graph.n_nodes))
        p=f.fermi_dirac(d,np.full(self.graph.n_nodes,7.))
        alpha=.03*rng.normal(size=len(self.graph.edges))
        value=f.evaluate(d,p,alpha)
        potential=self.ports.solve(value.edge_current_A,self.current)
        check=self.ports.gauge_power_check(d,value.gap_gradient_J_per_gapbar,
            value.edge_current_A,potential.phi_V)
        self.assertLess(check.noether_residual_max_J,1e-36)
        self.assertLess(abs(check.power_residual_W),1e-17)
        self.assertLess(abs(potential.boundary_power_W-potential.normal_joule_W-check.field_phase_power_W),1e-17)

    def test_CM_keeps_all_three_states_baseline_and_rejects_reversal(self):
        js=np.zeros(len(self.graph.edges));cm=ThesisCircuitParameters(Lk_ext_H=110e-9,V_bias_V=.21)
        y=np.array([22e-6,20e-6,7e-5])
        c=self.ports.coupled(js,cm,y,baseline_Vout_V=3e-6)
        np.testing.assert_allclose(c.circuit_rate_SI_per_s,cm.rhs(y,c.electrical.Vdev_V),rtol=0.)
        self.assertAlmostEqual(c.baseline_subtracted_Vout_V,50*(y[0]-y[1])-3e-6,places=17)
        self.assertLess(abs(c.circuit_power.residual_W),1e-18)
        self.assertGreater(c.circuit_energy_J,0.)
        with self.assertRaisesRegex(ValueError,"reversal"):
            self.ports.coupled(js,cm,[22e-6,-1e-6,0.])

    def test_incomplete_cuts_and_inactive_disconnection_rejected(self):
        with self.assertRaisesRegex(ValueError,"every node"):
            EffectiveNeutralPorts(self.graph,608.,self.left[1:],self.right)
        from pysnspd.experimental.thermal_spatial_usadel import ThermalGraph
        g=ThermalGraph(self.graph.area_weights,self.graph.edges,np.zeros(len(self.graph.edges)),self.graph.coordinates_bar)
        with self.assertRaisesRegex(ValueError,"connect"):
            EffectiveNeutralPorts(g,608.,self.left,self.right)

    def test_frozen_inductance_partitions_wire_only(self):
        a=freeze_exterior_inductance(.002,320e-9)
        self.assertEqual(a.added_inductance_H,96e-9)
        self.assertAlmostEqual(a.unresolved_bulk_inductance_H,.002*(5e-6-320e-9),places=22)
        self.assertEqual(a.total_H,96e-9+a.unresolved_bulk_inductance_H)
        b=freeze_exterior_inductance(.002,5e-6)
        self.assertEqual(b.total_H,96e-9)
        with self.assertRaises(ValueError):
            freeze_exterior_inductance(.002,6e-6)

    def test_bulk_load_uses_radial_traction_and_declared_phase_current(self):
        g=self.graph;d=np.exp(.07j*g.coordinates_bar[:,0])*1.6
        bulk=d*1.01;gradient=np.zeros(g.n_nodes,complex)
        fixed=np.r_[self.left,self.right]
        gradient[fixed]=(3e-23+5e-24j)*bulk[fixed]/abs(bulk[fixed])
        load=bulk_boundary_load(self.ports,d,bulk,gradient,self.current)
        self.assertLess(max(abs(load.phase_torque_residual_J)),1e-36)
        np.testing.assert_allclose(load.radial_load_J_per_gapbar[fixed],
            3e-23*d[fixed]/abs(d[fixed]),rtol=3e-15)
        self.assertAlmostEqual(sum(load.nodal_target_injection_A[self.left]),self.current,places=18)
        self.assertAlmostEqual(sum(load.nodal_target_injection_A[self.right]),-self.current,places=18)
        np.testing.assert_array_equal(load.full_load_J_per_gapbar[self.ports.free],np.zeros(len(self.ports.free)))
        rate=1e10*(1+.1j)*d
        power=boundary_material_power_W(load.full_load_J_per_gapbar,rate)
        self.assertAlmostEqual(power,float(np.real(np.vdot(load.full_load_J_per_gapbar,rate))),places=25)


if __name__=="__main__":
    unittest.main()
