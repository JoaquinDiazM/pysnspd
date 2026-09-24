"""Admission, physical invariance and provenance at the DC campaign boundary."""
import copy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from pysnspd.experimental.thermal_spatial_usadel import rectangular_graph
from sandbox.stage4_core.coupled_campaign import sha
from sandbox.stage5_prephoton import run_prephoton_dc as workflow


class PrephotonCampaignTests(unittest.TestCase):
    def setUp(self):
        self.temporary=tempfile.TemporaryDirectory()
        self.directory=Path(self.temporary.name)
        self.reference=self.directory/'reference.npz'
        graph=rectangular_graph(3,3,2.,1.)
        x=graph.coordinates_bar[:,0]
        fixed=np.flatnonzero((x==x.min()) | (x==x.max()))
        gap=1.4*np.exp(.4j*x)
        eps=2*np.pi*(.9/8.65)*(np.arange(2)+.5)
        # A small serializable fixture, never claimed to be a spectral solve.
        np.savez_compressed(self.reference,area_weights=graph.area_weights,
            edges=graph.edges,conductance=graph.conductance,
            coordinates_bar=graph.coordinates_bar,boundary_nodes=fixed,
            delta_bar=gap,epsilon_bar=eps,u=gap[None,:]/eps[:,None],
            alpha=np.zeros(len(graph.edges)))
        self.summary=dict(admitted_for_dark_hold=True,reference_sha256=sha(self.reference),
            bulk=dict(gap_bulk_bar=1.4,current_target_A=21.5e-6),
            fixed_external_inductance_H=105.4e-9)
        self.write_summary()
        self.plan=dict(hold_template=dict(schema='pysnspd.stage5.prephoton_dc_hold.v1',
            T_K=.9,Tc_K=8.65,tau_ee_Tc_ps=6.,tau_ep_Tc_ps=24.7,
            duration_ps=10.,delta0_over_kBTc=1.764,sheet_resistance_ohm=608.,
            circuit=dict(R_bias_ohm=1e4,L_bias_H=1e-6,R_load_ohm=50.,C_couple_F=100e-12)),
            euler_safety=.4,maximum_dt_ps=.0001,
            source_hash_normalization='CRLF-to-LF',
            source_sha256={'pinned_fixture.py':'a'*64})
        self.item=dict(id='dc',reference='dc',step_multiplier=1.)

    def tearDown(self):
        self.temporary.cleanup()

    def write_summary(self):
        (self.directory/'summary.json').write_text(json.dumps(self.summary),encoding='utf8')

    def test_refining_time_changes_no_physical_parameter_or_inductance(self):
        original=copy.deepcopy(self.plan)
        with patch.object(workflow,'explicit_step_bound',return_value={'primary_step_ps':.00025}):
            coarse=workflow.hold_plan(self.plan,self.item,self.reference)
            fine=workflow.hold_plan(self.plan,dict(self.item,step_multiplier=.5),self.reference)
        self.assertEqual(coarse['dt_ps'],.0001)
        self.assertEqual(fine['dt_ps'],.00005)
        self.assertEqual(coarse['reference_sha256'],sha(self.reference))
        self.assertEqual(coarse['circuit']['Lk_ext_H'],105.4e-9)
        self.assertEqual(coarse['bias_current_A'],21.5e-6)
        self.assertEqual(coarse['duration_ps'],fine['duration_ps'])
        for result in (coarse,fine):
            result.pop('dt_ps');result.pop('observation_cadence')
        self.assertEqual(coarse,fine)
        self.assertEqual(self.plan,original)

    def test_step_estimate_receives_real_bulk_amplitude_and_honors_the_bound(self):
        observed={}
        def bound(graph,gap,epsilon,plan):
            observed.update(gap=gap.copy(),epsilon=epsilon.copy(),plan=plan)
            self.assertEqual(len(gap),graph.n_nodes)
            return {'primary_step_ps':.000025}
        with patch.object(workflow,'explicit_step_bound',side_effect=bound):
            result=workflow.hold_plan(self.plan,self.item,self.reference)
        np.testing.assert_array_equal(observed['gap'],np.full(9,1.4,complex))
        self.assertEqual(observed['plan']['euler_safety'],.4)
        self.assertEqual(result['dt_ps'],.000025)
        self.assertIn('half-step comparison',result['step_estimate']['scope'])

    def test_unadmitted_reference_is_rejected_before_any_numerical_estimate(self):
        self.summary['admitted_for_dark_hold']=False;self.write_summary()
        with patch.object(workflow,'explicit_step_bound') as estimator:
            with self.assertRaisesRegex(ValueError,'admitted, unchanged'):
                workflow.hold_plan(self.plan,self.item,self.reference)
            estimator.assert_not_called()

    def test_mutated_reference_is_rejected_before_loading_fields(self):
        with self.reference.open('ab') as stream:
            stream.write(b'modified after admission')
        with patch.object(workflow,'load_reference') as loader:
            with self.assertRaisesRegex(ValueError,'admitted, unchanged'):
                workflow.hold_plan(self.plan,self.item,self.reference)
            loader.assert_not_called()

    def test_workflow_pins_sources_and_reference_plan_contents(self):
        source=self.directory/'source.py';source.write_bytes(b'value = 1\r\n')
        reference_plan=self.directory/'reference_plan.json';reference_plan.write_text('{"current_A":2.15e-5}\n')
        plan=dict(schema=workflow.SCHEMA,source_hash_normalization='CRLF-to-LF',
            source_sha256={'source.py':hashlib.sha256(b'value = 1\n').hexdigest()},
            references=[dict(id='dc',plan='reference_plan.json',sha256=sha(reference_plan),estimated_seconds=10.)],
            holds=[dict(id='dc',reference='dc',step_multiplier=1.)],
            maximum_parallel_cases=2,maximum_workers_per_case=13)
        with patch.object(workflow,'ROOT',self.directory):
            workflow.validate_workflow(plan)
            source.write_bytes(b'value = 1\n')
            workflow.validate_workflow(plan)
            source.write_text('value = 2\n')
            with self.assertRaisesRegex(ValueError,'Pinned source changed'):
                workflow.validate_workflow(plan)
            source.write_text('value = 1\n')
            reference_plan.write_text('{"current_A":1.55e-5}\n')
            with self.assertRaisesRegex(ValueError,'Reference plan changed'):
                workflow.validate_workflow(plan)

    def test_analyzer_does_not_reuse_stale_admission_when_raw_cases_are_missing(self):
        from sandbox.stage5_prephoton import analyze_prephoton_dc as analysis
        plan=dict(references=[dict(id='dc')],holds=[dict(id='dc',reference='dc'),
                dict(id='dc_half',reference='dc')],comparisons={'time':['dc','dc_half']},
            acceptance=dict(gap_comparison_relative=.01,current_comparison_relative=.02,
                            normalized_dark_readout_difference=.001))
        (self.directory/'workflow_plan.json').write_text(json.dumps(plan),encoding='utf8')
        (self.directory/'analysis.json').write_text(json.dumps({'prephoton_dc_admitted':True}),encoding='utf8')
        with patch.object(analysis,'_plots',return_value=[]):
            result=analysis.analyze(self.directory)
        self.assertFalse(result['prephoton_dc_admitted'])
        self.assertFalse(result['references']['dc']['available'])
        self.assertFalse(result['holds']['dc_half']['available'])
        self.assertFalse(result['comparisons']['time']['admitted'])
        saved=json.loads((self.directory/'analysis.json').read_text(encoding='utf8'))
        self.assertFalse(saved['prephoton_dc_admitted'])
        self.assertTrue((self.directory/'analysis.md').is_file())


if __name__=='__main__':
    unittest.main()
