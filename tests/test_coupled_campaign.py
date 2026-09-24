"""Short dummy-process checks; no SNSPD solves and no background computation."""
import contextlib
import io
import json
import os
from pathlib import Path
import tempfile
import unittest

from sandbox.stage4_core.coupled_campaign import (
    SCHEMA, allocate, exact_argv, run_campaign, validate_plan,
)
from sandbox.stage4_core.parallel_runtime import THREAD_VARIABLES, linux_resources


DUMMY = '''import argparse,json,os,pathlib,time
p=argparse.ArgumentParser()
p.add_argument('--workers',type=int);p.add_argument('--output');p.add_argument('--sleep',type=float)
p.add_argument('--code',type=int,default=0);p.add_argument('--literal')
a=p.parse_args(); out=pathlib.Path(a.output)
assert not out.exists();out.mkdir()
print(json.dumps({'event':'DUMMY_START','fraction':0.25,'eta_seconds':a.sleep}),flush=True)
time.sleep(a.sleep)
(out/'summary.json').write_text(json.dumps({'workers':a.workers,'affinity':sorted(os.sched_getaffinity(0)),
 'literal':a.literal,'threads':{k:os.environ.get(k) for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS')}}))
print(json.dumps({'event':'DUMMY_END','fraction':1.0,'eta_seconds':0.0}),flush=True)
raise SystemExit(a.code)
'''


def case(name, *, delay=.02, code=0):
    return dict(id=name, estimated_seconds=1., expected_outputs=['summary.json'],
        argv=['{python}', '-c', DUMMY, '--workers', '{workers}', '--output', '{output}',
              '--sleep', str(delay), '--code', str(code), '--literal', 'a; $(false) literal spaces'])


def plan(cases=None):
    return dict(schema=SCHEMA, maximum_parallel_cases=2, heartbeat_seconds=.1,
                cases=cases or [case('a'), case('b')])


def resources(memory_gib=123, quota=None):
    return dict(logical_cpus=list(range(32)), physical_core_groups=[[i, i+16] for i in range(16)],
        available_memory_bytes=memory_gib*1024**3, cpu_quota_units=quota)


class CoupledCampaignTests(unittest.TestCase):
    def test_shared_budget_counts_each_case_coordinator_and_nested_workers(self):
        result = allocate(plan(), resources())
        self.assertEqual(result['total_process_cpu_ceiling'], 28)
        self.assertEqual(result['numerical_worker_ceiling'], 25)
        self.assertEqual(result['shared_budget']['physical_cores_reserved'], 2)
        sets = [set(slot['cpus']) for slot in result['slots']]
        self.assertFalse(sets[0] & sets[1])
        self.assertNotIn(result['shared_budget']['coordinator_cpu'], sets[0] | sets[1])
        self.assertLessEqual(result['estimated_memory_reservation_bytes'], .9*123*1024**3)
        tight = allocate(plan(), resources(memory_gib=8, quota=8))
        self.assertEqual(tight['parallel_cases'], 1)
        self.assertLessEqual(tight['total_process_cpu_ceiling'], .9*8)
        self.assertLessEqual(tight['estimated_memory_reservation_bytes'], .9*8*1024**3)

    def test_exact_argv_keeps_literals_and_plan_validation_prevents_output_collisions(self):
        argv = exact_argv(case('safe'), 3, Path('output with spaces'))
        self.assertEqual(argv[-1], 'a; $(false) literal spaces')
        self.assertIn('3', argv)
        self.assertEqual(argv[2], DUMMY)
        with self.assertRaises(ValueError):
            validate_plan(plan([case('same'), case('same')]))
        bad = case('bad'); bad['expected_outputs'] = ['../escape.json']
        with self.assertRaises(ValueError):
            validate_plan(plan([bad]))

    def test_existing_output_is_refused_before_resource_inspection(self):
        with tempfile.TemporaryDirectory() as folder:
            saved = Path(folder)/'keep.txt'; saved.write_text('untouched')
            with self.assertRaises(FileExistsError):
                run_campaign(plan(), folder)
            self.assertEqual(saved.read_text(), 'untouched')

    @unittest.skipUnless(hasattr(os, 'sched_setaffinity'), 'Linux affinity integration test')
    def test_parallel_failure_does_not_skip_siblings_and_success_is_durable(self):
        detected = linux_resources()
        try:
            allocation = allocate(plan(), detected)
        except ValueError as error:
            self.skipTest(str(error))
        if allocation['parallel_cases'] < 2:
            self.skipTest('The effective runtime budget cannot admit two dummy cases')
        previous = os.sched_getaffinity(0)
        with tempfile.TemporaryDirectory() as folder, contextlib.redirect_stdout(io.StringIO()):
            root = Path(folder)/'mixed'
            requested = plan([case('slow_ok', delay=.4), case('fast_failure', delay=.1, code=7),
                              case('later_ok', delay=.02)])
            result = run_campaign(requested, root, resources=detected)
            self.assertEqual(result['status'], 'COMPLETED_WITH_FAILURES')
            by_id = {row['id']: row for row in result['cases']}
            self.assertEqual(by_id['fast_failure']['returncode'], 7)
            self.assertEqual(by_id['slow_ok']['status'], 'SUCCEEDED')
            self.assertEqual(by_id['later_ok']['status'], 'SUCCEEDED')
            self.assertLess(by_id['fast_failure']['start_elapsed_seconds'], by_id['slow_ok']['end_elapsed_seconds'])
            self.assertLess(by_id['slow_ok']['start_elapsed_seconds'], by_id['fast_failure']['end_elapsed_seconds'])
            for record in result['cases']:
                self.assertEqual(record['attempts'], 1)
                payload = json.loads((Path(record['output'])/'summary.json').read_text())
                self.assertEqual(payload['affinity'], sorted(record['affinity_cpus']))
                self.assertEqual(payload['workers'], record['workers'])
                self.assertEqual(set(payload['threads'].values()), {'1'})
                self.assertEqual(payload['literal'], 'a; $(false) literal spaces')
                self.assertTrue((root/record['console_log']).is_file())
                self.assertTrue((root/'case_records'/f'{record["id"]}.json').is_file())
            events = [json.loads(line) for line in (root/'progress.jsonl').read_text().splitlines()]
            self.assertTrue(any(item['event'] == 'HEARTBEAT' and item['eta_seconds'] is not None for item in events))
            self.assertEqual(json.loads((root/'campaign.json').read_text())['status'], result['status'])
            success = run_campaign(plan([case('only')]), Path(folder)/'success', resources=detected)
            self.assertEqual(success['status'], 'SUCCEEDED')
            missing = case('missing'); missing['expected_outputs'] = ['absent.json']
            lost = run_campaign(plan([missing]), Path(folder)/'missing', resources=detected)
            self.assertEqual(lost['cases'][0]['status'], 'MISSING_OUTPUT')
            self.assertEqual(lost['cases'][0]['returncode'], 0)
        self.assertEqual(os.sched_getaffinity(0), previous)


if __name__ == '__main__':
    unittest.main()
