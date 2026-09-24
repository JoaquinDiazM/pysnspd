"""Small stationary references; no device-size calculation or time evolution."""
from pathlib import Path
import tempfile
import unittest

import numpy as np
from pysnspd.experimental.thermal_spatial_usadel import ThermalGraph, rectangular_graph
from sandbox.stage4_core.biased_strip_reference import (
    SCHEMA, contact_partition, initial_state, solve_mode, solve_reference,
)


def graph():
    base=rectangular_graph(7,3,10.,2.)
    ends=np.flatnonzero((base.coordinates_bar[:,0]==0.) | (base.coordinates_bar[:,0]==10.))
    return ThermalGraph(base.area_weights,base.edges,base.conductance,base.coordinates_bar,ends)


def plan(bias):
    return dict(schema=SCHEMA,T_K=.5,Tc_K=1.,phase_bias=bias,matsubara_count=8,
        max_outer_iterations=80,max_newton_iterations=40,checkpoint_cadence=5,
        gap_rms_relative_tolerance=1e-5,spectral_tolerance=1e-8)


def evaluator(mesh,settings):
    def evaluate(d,epsilon,warm,sweep):
        solutions=[solve_mode(mesh,settings,d,e,None if warm is None else warm[n]) for n,e in enumerate(epsilon)]
        return solutions,[]
    return evaluate


class BiasedStripReferenceTests(unittest.TestCase):
    def test_natural_sides_and_declared_contact_phase(self):
        mesh=graph();settings=plan(8.)
        left,right,fixed=contact_partition(mesh)
        d,d0,eps=initial_state(mesh,settings)
        np.testing.assert_allclose(d[left],d0*np.exp(-4j))
        np.testing.assert_allclose(d[right],d0*np.exp(4j))
        # The initial unwrapped phase gradient is imposed across the full strip,
        # rather than selecting the shorter phase difference modulo 2*pi.
        edge=mesh.edges[np.abs(np.diff(mesh.coordinates_bar[mesh.edges,0],axis=1)[:,0])>0]
        np.testing.assert_allclose(np.angle(np.conj(d[edge[:,0]])*d[edge[:,1]]),8/6)
        with self.assertRaises(ValueError):
            contact_partition(rectangular_graph(7,3,10.,2.))

    def test_uniform_limit_converges_in_one_block_and_final_fields_are_consistent(self):
        mesh=graph();settings=plan(0.)
        with tempfile.TemporaryDirectory() as folder:
            output=Path(folder);(output/'checkpoints').mkdir()
            result=solve_reference(mesh,settings,evaluator(mesh,settings),output=output)
            self.assertEqual(result['status'],'FINITE_SUM_STATIONARY')
            self.assertEqual(result['completed_sweeps'],1)
            self.assertLess(result['metrics']['gap_mass_rms_relative'],1e-12)
            with np.load(output/'reference.npz') as saved:
                np.testing.assert_allclose(saved['f'],saved['u']*saved['g'])
                np.testing.assert_allclose(abs(saved['f'])**2+saved['g']**2,1.)
                np.testing.assert_allclose(saved['current_bar'],0.,atol=1e-13)
                self.assertIn('delta_bar',saved)
                self.assertIn('next_delta_bar',saved)

    def test_opposite_bias_reverses_current_and_material_budget_exit_is_not_admission(self):
        mesh=graph()
        values=[]
        for bias in (.7,-.7):
            settings=plan(bias)
            result=solve_reference(mesh,settings,evaluator(mesh,settings))
            self.assertEqual(result['status'],'FINITE_SUM_STATIONARY')
            metrics=result['metrics'];values.append(metrics)
            left=metrics['contact_reactions']['left'];right=metrics['contact_reactions']['right']
            self.assertLess(abs(left['outgoing_link_current_bar']+right['outgoing_link_current_bar']),1e-7)
            self.assertLess(abs(left['outgoing_link_current_bar']+left['total_phase_reaction_bar']),1e-9)
            self.assertLess(metrics['maximum_noether_residual'],1e-9)
            self.assertGreater(metrics['minimum_matsubara_g'],0.)
        self.assertAlmostEqual(values[0]['minimum_gap_relative'],values[1]['minimum_gap_relative'],places=12)
        self.assertAlmostEqual(values[0]['contact_reactions']['left']['outgoing_link_current_bar'],
            -values[1]['contact_reactions']['left']['outgoing_link_current_bar'],places=12)
        settings=plan(.7);settings['max_outer_iterations']=1
        incomplete=solve_reference(mesh,settings,evaluator(mesh,settings))
        self.assertEqual(incomplete['status'],'INCOMPLETE_MAXIMUM_ITERATIONS')
        self.assertFalse(incomplete['criterion_met'])


if __name__=='__main__':
    unittest.main()
