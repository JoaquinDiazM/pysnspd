"""Physical invariants of the DC interior-window preparation, without long solves."""
import unittest
from types import SimpleNamespace

import numpy as np

from pysnspd.experimental import thermal_spatial_usadel as thermal
from pysnspd.experimental.bulk_current_reference import physical_scales
from sandbox.stage5_prephoton.prepare_dc_reference import (
    bulk_state, fixed_inductance_partition, metrics_bulk,
)


class PrephotonDCReferenceTests(unittest.TestCase):
    def setUp(self):
        self.plan=dict(T_K=.9,Tc_K=8.65,matsubara_count=256,current_A=21.5e-6,
            width_m=80e-9,diffusion_m2_s=.5e-4,sheet_resistance_ohm=608.,
            active_length_m=5e-6,added_inductance_H=96e-9)
        self.scales=physical_scales(Tc_K=8.65,diffusion_m2_s=.5e-4,sheet_resistance_ohm=608.)

    def graph(self,length_m=240e-9,width_m=80e-9):
        ell=self.scales['ell0_m']
        original=thermal.rectangular_graph(9,4,length_m/ell,width_m/ell)
        x=original.coordinates_bar[:,0]
        ends=np.flatnonzero((x==x.min()) | (x==x.max()))
        return thermal.ThermalGraph(original.area_weights,original.edges,
            original.conductance,original.coordinates_bar,ends)

    def test_bulk_boundary_carries_the_same_suppressed_gap_and_spectrum(self):
        graph=self.graph()
        ref,fields=bulk_state(graph,self.plan)
        # A current-biased *interior* has no enforced recovery to unbiased BCS
        # at either cut. The same depressed amplitude holds on all nodes.
        self.assertLess(ref.gap_bar,1.7)
        np.testing.assert_allclose(abs(fields['delta_bar']),ref.gap_bar,atol=1e-15)
        np.testing.assert_allclose(abs(fields['u']),
            np.broadcast_to(ref.u[:,None],fields['u'].shape),atol=1e-15)
        fixed=graph.boundary_nodes
        eps=ref.epsilon_bar[:,None]
        # Check the actual depairing equation with complex spatial phase.
        residual=(eps*fields['u'][:,fixed]+ref.q_bar**2*fields['f'][:,fixed]
                  -fields['delta_bar'][None,fixed])
        np.testing.assert_allclose(residual,0.,atol=4e-14)
        self.assertGreater(np.max(abs(fields['u'][:,fixed]
            -fields['delta_bar'][None,fixed]/eps)),.1)

    def test_inductance_is_device_invariant_when_window_length_changes(self):
        short,long=self.graph(120e-9),self.graph(480e-9)
        ref,_=bulk_state(short,self.plan)
        a=fixed_inductance_partition(short,ref,self.plan)
        b=fixed_inductance_partition(long,ref,self.plan)
        self.assertAlmostEqual(a['total_equilibrium_inductance_H'],
                               b['total_equilibrium_inductance_H'],delta=1e-22)
        self.assertAlmostEqual(a['fixed_external_inductance_H']+
            a['resolved_equilibrium_inductance_H'],
            b['fixed_external_inductance_H']+b['resolved_equilibrium_inductance_H'],delta=1e-22)
        self.assertLess(b['fixed_external_inductance_H'],a['fixed_external_inductance_H'])
        self.assertGreater(a['fixed_external_inductance_H'],self.plan['added_inductance_H'])
        self.assertGreater(a['total_equilibrium_inductance_H'],96e-9)
        self.assertLess(a['total_equilibrium_inductance_H'],110e-9)

    def test_uniform_dc_flux_metrics_are_orientation_invariant(self):
        graph=self.graph()
        ref,fields=bulk_state(graph,self.plan)
        x=graph.coordinates_bar[:,0]
        tail,head=graph.edges.T
        # Exact finite-volume flux for a uniform sheet current: current density
        # times the transverse dual-face width, with each edge's orientation.
        current=ref.current_density_bar*graph.conductance*(x[head]-x[tail])
        measured=metrics_bulk(graph,fields['delta_bar'],SimpleNamespace(current_bar=current),
            ref,self.plan,self.scales['current_unit_A'])
        self.assertAlmostEqual(measured['current_reference_A'],self.plan['current_A'],delta=1e-14)
        self.assertLess(measured['crosscut_max_relative_error'],2e-14)
        self.assertLess(measured['gap_peak_to_peak_relative'],1e-14)
        reverse=thermal.ThermalGraph(graph.area_weights,graph.edges[:,::-1],
            graph.conductance,graph.coordinates_bar,graph.boundary_nodes)
        other=metrics_bulk(reverse,fields['delta_bar'],SimpleNamespace(current_bar=-current),
            ref,self.plan,self.scales['current_unit_A'])
        for key in ('current_left_A','current_right_A','current_reference_A',
                    'crosscut_max_relative_error'):
            self.assertAlmostEqual(measured[key],other[key],delta=1e-14)

    def test_metrics_detect_nonconservation_and_complete_current_loss(self):
        graph=self.graph();ref,fields=bulk_state(graph,self.plan)
        x=graph.coordinates_bar[:,0];tail,head=graph.edges.T
        current=ref.current_density_bar*graph.conductance*(x[head]-x[tail])
        # A local missing flux changes an internal crosscut, even when end
        # currents and their nominal bias agree.
        cut=(x[tail]>x.min()) & (x[head]<x.max()) & (x[head]>x[tail])
        defect=np.flatnonzero(cut)[0]
        broken=current.copy();broken[defect]=0.
        bad=metrics_bulk(graph,fields['delta_bar'],SimpleNamespace(current_bar=broken),
            ref,self.plan,self.scales['current_unit_A'])
        self.assertLess(bad['current_target_relative_error'],1e-12)
        self.assertGreater(bad['crosscut_max_relative_error'],.1)
        zero=metrics_bulk(graph,fields['delta_bar'],SimpleNamespace(current_bar=0*current),
            ref,self.plan,self.scales['current_unit_A'])
        self.assertEqual(zero['current_target_relative_error'],1.)
        self.assertEqual(zero['crosscut_max_relative_error'],0.)

    def test_material_geometry_and_inductance_partition_mismatch_rejected(self):
        with self.assertRaisesRegex(ValueError,'mesh width'):
            bulk_state(self.graph(width_m=100e-9),self.plan)
        ref,_=bulk_state(self.graph(),self.plan)
        with self.assertRaisesRegex(ValueError,'partition'):
            fixed_inductance_partition(self.graph(6e-6),ref,self.plan)
        with self.assertRaisesRegex(ValueError,'partition'):
            fixed_inductance_partition(self.graph(),ref,dict(self.plan,added_inductance_H=-1.))
        with self.assertRaisesRegex(ValueError,'nonzero'):
            bulk_state(self.graph(),dict(self.plan,current_A=0.))


if __name__=='__main__':
    unittest.main()
