"""Runtime tests with synthetic queries; no spectral or physical campaign."""
from concurrent.futures import ThreadPoolExecutor
import contextlib
import io
import json
import os
from pathlib import Path
from types import SimpleNamespace
import tempfile
import time
import unittest
from unittest.mock import patch

import numpy as np
from sandbox.stage4_core import coupled_response as runner
from sandbox.stage4_core.coupled_campaign import allocate
from tests.test_coupled_campaign import plan as campaign_plan, resources


def plan():
    return dict(schema='pysnspd.stage4.coupled_response.v1', T_K=.9, Tc_K=8.65,
        matsubara_count=2, mode_count=1, mode_indices=[0], omegas=[.1],
        eta_over_gap=[.01], orders=[2], cutoffs=[8.], scope='Synthetic runtime test')


class CoupledResponseRuntimeTests(unittest.TestCase):
    def test_shared_admission_has_no_nested_ten_percent_reservation(self):
        allocation = allocate(campaign_plan(), resources())
        slot = allocation['slots'][0]
        inherited = dict(parent_pid=os.getppid(), cpus=slot['cpus'], workers=slot['workers'],
            shared_budget=allocation['shared_budget'], total_process_cpu_ceiling=28)
        visible = dict(resources(), logical_cpus=slot['cpus'])
        admitted = runner.admit_runtime(slot['workers'], visible, inherited)
        self.assertEqual(admitted['workers'], slot['workers'])
        self.assertEqual(admitted['coordinator_cpu'], slot['cpus'][-1])
        self.assertNotIn(admitted['coordinator_cpu'], admitted['worker_affinity_cpus'])
        for changed in (dict(visible, logical_cpus=slot['cpus'][:-1]),
                        dict(visible, cpu_quota_units=8), dict(visible, available_memory_bytes=2*1024**3)):
            with self.assertRaises(ValueError):
                runner.admit_runtime(slot['workers'], changed, inherited)
        with self.assertRaises(ValueError):
            runner.admit_runtime(28, resources())
        with self.assertRaises(ValueError):
            runner.admit_runtime(0, resources())

    def test_existing_output_refused_without_inspecting_resources(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); source = root/'plan.json'
            source.write_text(json.dumps(plan()))
            with patch('sys.argv', ['runner', '--plan', str(source), '--output', str(root)]), \
                 patch.object(runner, 'linux_resources', side_effect=AssertionError('must not inspect')):
                with self.assertRaises(FileExistsError):
                    runner.main()

    def test_final_manifest_records_failure_and_preserves_independent_success(self):
        graph = SimpleNamespace(n_nodes=3, coordinates_bar=np.zeros((3, 2)), area_weights=np.ones(3))
        def dummy(task):
            if task['kind'] == 'port':
                raise ValueError('declared dummy failure')
            return dict(kind='mode', status='DUMMY_SUCCESS')
        def pool(**kw):
            return ThreadPoolExecutor(max_workers=kw['max_workers'])
        with tempfile.TemporaryDirectory() as folder, contextlib.redirect_stdout(io.StringIO()):
            args = SimpleNamespace(output=Path(folder), workers=2)
            manifest = dict(status='RUNNING', sources={'dummy.py':'test'})
            with patch.object(runner, 'mesh_modes', return_value=(graph,np.ones(1),np.zeros((3,1)),np.zeros(3),1.,0.)), \
                 patch.object(runner, 'uniform_reference', return_value=(np.ones(3),None,None)), \
                 patch.object(runner, 'ProcessPoolExecutor', side_effect=pool), \
                 patch.object(runner, 'execute_task', side_effect=dummy), \
                 patch.object(runner, 'runtime_telemetry', return_value={'pid':os.getpid()}):
                code = runner.run_admitted(args,plan(),Path('unused'),{'worker_affinity_cpus':[0,1]},manifest,time.monotonic())
            result = json.loads((args.output/'results.json').read_text())
            final = json.loads((args.output/'manifest.json').read_text())
            self.assertEqual(code,1)
            self.assertEqual(final['status'],'COMPLETED_WITH_FAILURES')
            self.assertEqual(len(result['records']),1)
            self.assertEqual(len(result['failures']),1)
            self.assertEqual(final['outputs']['results.json'],runner.sha(args.output/'results.json'))
            self.assertEqual(final['sources'],{'dummy.py':'test'})
            self.assertTrue((args.output/'query_00000.failed.json').is_file())
            self.assertTrue((args.output/'query_00001.json').is_file())


if __name__ == '__main__':
    unittest.main()
