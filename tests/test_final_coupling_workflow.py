"""Small workflow/receipt tests; no physical queries or device solves."""
import contextlib
import hashlib
import io
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from sandbox.stage4_core import run_final_coupling as workflow
from sandbox.stage4_core.biased_coupled_response import (
    PhaseProgress, admitted_reference, durable_run,
)
from tests.test_coupled_campaign import resources


def make_plan(folder):
    settings=Path(folder)/'settings.json';settings.write_text('{}')
    digest=hashlib.sha256(settings.read_bytes()).hexdigest()
    def item(name,reference=None):
        row=dict(id=name,plan=str(settings),sha256=digest,estimated_seconds=1.)
        if reference is not None:
            row['reference']=reference
        return row
    value=dict(schema='pysnspd.stage4.final_coupling.v1',references=[item('bias2'),item('bias8')],
        responses=[item('control2','bias2'),item('base8','bias8'),item('refine8','bias8')])
    path=Path(folder)/'workflow.json';path.write_text(json.dumps(value))
    return value,path


class FinalCouplingWorkflowTests(unittest.TestCase):
    def test_dry_run_constructs_exact_routes_without_starting_children(self):
        with tempfile.TemporaryDirectory() as folder,contextlib.redirect_stdout(io.StringIO()) as text:
            plan,path=make_plan(folder);output=Path(folder)/'new'
            with patch('sys.argv',['runner','--plan',str(path),'--output-root',str(output)]), \
                 patch.object(workflow,'linux_resources',return_value=resources()), \
                 patch.object(workflow,'run_campaign',side_effect=AssertionError('dry-run executed a child')):
                self.assertEqual(workflow.main(),0)
            preview=json.loads(text.getvalue())
            self.assertEqual(preview['status'],'DRY_RUN')
            self.assertFalse(output.exists())
            self.assertEqual(preview['references']['cases'][0]['expected_outputs'],['reference.npz','manifest.json'])
            for row in preview['responses']:
                location=row['argv'][row['argv'].index('--reference')+1]
                self.assertTrue(location.endswith('reference.npz'))
                self.assertIn('manifest.json',row['expected_outputs'])
            self.assertLessEqual(preview['allocation']['total_process_cpu_ceiling'],28)

    def test_failed_reference_blocks_only_its_dependent_responses(self):
        with tempfile.TemporaryDirectory() as folder,contextlib.redirect_stdout(io.StringIO()):
            plan,path=make_plan(folder);output=Path(folder)/'new'
            stages=[]
            def run(campaign,root,**kw):
                stages.append(campaign)
                if len(stages)==1:
                    return dict(status='COMPLETED_WITH_FAILURES',cases=[
                        dict(id='bias2',status='SUCCEEDED'),dict(id='bias8',status='FAILED')])
                return dict(status='SUCCEEDED',cases=[dict(id='control2',status='SUCCEEDED')])
            with patch('sys.argv',['runner','--plan',str(path),'--output-root',str(output),'--execute']), \
                 patch.object(workflow,'linux_resources',return_value=resources()), \
                 patch.object(workflow,'run_campaign',side_effect=run):
                self.assertEqual(workflow.main(),2)
            self.assertEqual([case['id'] for case in stages[1]['cases']],['control2'])
            dependencies=json.loads((output/'dependencies.json').read_text())
            self.assertEqual(dependencies['blocked_responses'],['base8','refine8'])
            final=json.loads((output/'workflow_result.json').read_text())
            self.assertEqual(final['status'],'COMPLETED_WITH_INCOMPLETE_CASES')

    def test_reference_admission_requires_matching_final_stationary_receipt(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);path=root/'reference.npz';path.write_bytes(b'unchanged reference bytes')
            receipt=dict(status='FINITE_SUM_STATIONARY',latest_checkpoint=dict(path='reference.npz',
                sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
            (root/'manifest.json').write_text(json.dumps(receipt))
            self.assertEqual(admitted_reference(path)['status'],'FINITE_SUM_STATIONARY')
            receipt['status']='INCOMPLETE_MAXIMUM_ITERATIONS'
            (root/'manifest.json').write_text(json.dumps(receipt))
            with self.assertRaisesRegex(ValueError,'gap criterion'):
                admitted_reference(path)
            receipt['status']='FINITE_SUM_STATIONARY';receipt['latest_checkpoint']['sha256']='changed'
            (root/'manifest.json').write_text(json.dumps(receipt))
            with self.assertRaisesRegex(ValueError,'bytes differ'):
                admitted_reference(path)

    def test_eta_uses_energy_phase_cost_and_progress_is_throttled(self):
        now=[0.];records=[]
        meter=PhaseProgress({'thermal_tangent':256,'coupled_spectral_kinetic':4},
            emit=lambda line,**kw:records.append(json.loads(line)),clock=lambda:now[0])
        meter.begin('thermal_tangent')
        now[0]=1.;meter.completed(256)
        self.assertIsNone(records[-1]['eta_seconds'])
        meter.begin('coupled_spectral_kinetic')
        now[0]=11.;meter.completed()
        self.assertEqual(records[-1]['eta_seconds'],30.)
        count=len(records)
        now[0]=12.;meter.completed()
        self.assertEqual(len(records),count)
        now[0]=13.;meter.completed(2)
        self.assertEqual(records[-1]['eta_seconds'],0.)
        self.assertEqual(records[-1]['fraction'],1.)

    def test_postprocessing_failure_is_durable_and_restores_affinity(self):
        with tempfile.TemporaryDirectory() as folder,contextlib.redirect_stdout(io.StringIO()):
            root=Path(folder);(root/'integrated_moments.npz').write_bytes(b'completed evidence')
            manifest=dict(status='RUNNING',phase='postprocessing')
            with patch.object(os,'sched_getaffinity',create=True,return_value={2,3}), \
                 patch.object(os,'sched_setaffinity',create=True) as affinity:
                with self.assertRaisesRegex(ValueError,'singular postprocessing'):
                    with durable_run(root,manifest,{'coordinator_cpu':3}):
                        raise ValueError('singular postprocessing')
            self.assertEqual(affinity.call_args_list[-1].args,(0,{2,3}))
            final=json.loads((root/'manifest.json').read_text())
            self.assertEqual(final['status'],'FAILED')
            self.assertEqual(final['failure']['phase'],'postprocessing')
            self.assertEqual((root/'integrated_moments.npz').read_bytes(),b'completed evidence')


if __name__=='__main__':
    unittest.main()
