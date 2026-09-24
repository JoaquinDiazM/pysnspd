"""Exact-key scheduling and budgeting; no causal solver or long calculation."""
from concurrent.futures import ThreadPoolExecutor
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import time
import unittest

import numpy as np

from test_stage4_core_parameter import AnalyticCatalogue

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT/'sandbox/stage4_core'
sys.path.insert(0, str(SCRIPT))
import controls
import parallel_runtime as parallel


def fake_kernel(resolution, key):
    time.sleep(.001)
    return (np.array(key),), .001, 1


class TestStage4Parallel(unittest.TestCase):
    def resources(self, cores=16, memory_GiB=116):
        return dict(logical_cpus=list(range(2*cores)),
            physical_core_groups=[[i, i+cores] for i in range(cores)],
            available_memory_bytes=memory_GiB*1024**3, total_memory_bytes=123*1024**3)

    def test_geminga_budget_reserves_whole_cores_and_counts_coordinator(self):
        budget = parallel.resource_budget(self.resources())
        self.assertEqual(budget['workers'], 27)
        self.assertEqual(budget['total_processes'], 28)
        self.assertEqual(budget['physical_cores_used'], 14)
        self.assertEqual(set(range(32))-set(budget['selected_affinity_cpus']), {14, 15, 30, 31})
        self.assertNotIn(budget['coordinator_cpu'], budget['worker_affinity_cpus'])
        self.assertLessEqual(budget['estimated_memory_reservation_bytes'], budget['memory_limit_bytes'])

    def test_small_host_cannot_silently_exceed_fraction_or_memory(self):
        with self.assertRaisesRegex(ValueError, 'coordinator and worker'):
            parallel.resource_budget(self.resources(cores=1))
        with self.assertRaisesRegex(ValueError, 'RAM'):
            parallel.resource_budget(self.resources(memory_GiB=3))
        with self.assertRaises(ValueError):
            parallel.resource_budget(self.resources(), cpu_fraction=.95)
        low_memory = parallel.resource_budget(self.resources(memory_GiB=8))
        self.assertEqual(low_memory['workers'], 4)
        self.assertEqual(parallel.resource_budget(self.resources(), max_workers=4)['total_processes'], 5)

    def test_v2_limits_use_most_restrictive_visible_ancestor_and_headroom(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)/'cgroup'
            leaf = root/'team/job'
            leaf.mkdir(parents=True)
            proc = Path(temp)/'proc_cgroup'
            proc.write_text('0::/team/job\n')
            (root/'cpu.max').write_text('max 100000')
            (root/'team/cpu.max').write_text('400000 100000')
            (leaf/'cpu.max').write_text('800000 100000')
            (root/'memory.max').write_text('max')
            (root/'team/memory.max').write_text(str(20*1024**3))
            (root/'team/memory.current').write_text(str(15*1024**3))
            (leaf/'memory.max').write_text(str(10*1024**3))
            (leaf/'memory.current').write_text(str(2*1024**3))
            limits = parallel.cgroup_limits(proc, root)
            self.assertEqual(limits['cpu_quota_units'], 4.)
            self.assertEqual(limits['memory_headroom_bytes'], 5*1024**3)
            resources = self.resources()
            resources['cpu_quota_units'] = limits['cpu_quota_units']
            resources['available_memory_bytes'] = min(resources['available_memory_bytes'], limits['memory_headroom_bytes'])
            budget = parallel.resource_budget(resources)
            self.assertEqual(budget['logical_budget'], 3)
            self.assertLessEqual(budget['total_processes'], 3)
            self.assertEqual(budget['physical_cores_used'], 1)
            self.assertEqual(budget['physical_cores_reserved'], 15)
            self.assertLessEqual(budget['estimated_memory_reservation_bytes'], .9*5*1024**3)

    def test_v1_conventional_mounts_and_unlimited_cpu_are_supported(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)/'cgroup'
            cpu = root/'cpu,cpuacct/team/job'
            memory = root/'memory/team/job'
            cpu.mkdir(parents=True)
            memory.mkdir(parents=True)
            proc = Path(temp)/'proc_cgroup'
            proc.write_text('2:cpu,cpuacct:/team/job\n3:memory:/team/job\n')
            (cpu/'cpu.cfs_quota_us').write_text('-1')
            (cpu.parent/'cpu.cfs_quota_us').write_text('250000')
            (cpu.parent/'cpu.cfs_period_us').write_text('100000')
            (memory/'memory.limit_in_bytes').write_text(str(10*1024**3))
            (memory/'memory.usage_in_bytes').write_text(str(3*1024**3))
            limits = parallel.cgroup_limits(proc, root)
            self.assertEqual(limits['cpu_quota_units'], 2.5)
            self.assertEqual(limits['memory_headroom_bytes'], 7*1024**3)
            self.assertEqual(limits['cpu_records'][0]['status'], 'unlimited')

    def test_missing_cgroup_files_are_unavailable_and_readable_invalid_limits_fail(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)/'cgroup'
            root.mkdir()
            proc = Path(temp)/'proc_cgroup'
            proc.write_text('0::/\n')
            missing = parallel.cgroup_limits(proc, root)
            self.assertIsNone(missing['cpu_quota_units'])
            self.assertTrue(missing['unavailable'])
            (root/'cpu.max').write_text('10000 0')
            with self.assertRaisesRegex(ValueError, 'Invalid readable cgroup CPU'):
                parallel.cgroup_limits(proc, root)
            (root/'cpu.max').write_text('max 100000')
            (root/'memory.max').write_text('1000000')
            with self.assertRaisesRegex(ValueError, 'Invalid/incomplete finite cgroup memory'):
                parallel.cgroup_limits(proc, root)

    def test_cases_interleave_share_exact_keys_and_keep_neighbouring_floats_distinct(self):
        near = np.nextafter(.1, 1.)
        requests = [dict(id='first', resolution='r', keys=[(.1, .2), (.3, .4)]),
                    dict(id='second', resolution='r', keys=[(.1, .2), (near, .2)])]
        jobs = parallel.interleaved_jobs(requests)
        self.assertEqual(len(jobs), 3)
        self.assertEqual(jobs[0]['owners'], ['first', 'second'])
        self.assertNotEqual(jobs[0]['key'], jobs[-1]['key'])
        obtained = []
        with ThreadPoolExecutor(max_workers=2) as pool:
            parallel.prefetch_jobs(jobs, pool, 2, lambda job, *args: obtained.append(job), kernel=fake_kernel)
        self.assertEqual({job['key'] for job in obtained}, {job['key'] for job in jobs})

    def test_task_failure_is_propagated_without_retry(self):
        calls = []
        def fail(resolution, key):
            calls.append(key)
            raise RuntimeError('deliberate kernel failure')
        jobs = [dict(resolution='r', key=(.1, .2), owners=['one'])]
        with ThreadPoolExecutor(max_workers=1) as pool:
            with self.assertRaisesRegex(RuntimeError, 'deliberate kernel failure'):
                parallel.prefetch_jobs(jobs, pool, 1, lambda *args: None, kernel=fail)
        self.assertEqual(len(calls), 1)

    def test_traced_prefetch_exactly_covers_final_spatial_paths(self):
        plan = dict(material=dict(Tc_K=8.65, Tb_K=.9, width_m=80e-9, thickness_m=7e-9),
            geometry=dict(length_m=160e-9), mobility_pairs_ps=dict(memory=[.5, 2.47]))
        case = dict(id='light', phase='spatial', profile='suppressed', delta_reg=.05,
            elements_x=2, elements_y=1, degree=2, population='synthetic_fixed_count',
            finite_differences=False, boundary_response='fixed_prescribed_all_sides')
        source = AnalyticCatalogue()
        source.eta = 1e-8
        cat = controls.CachedDirect(source)
        keys = parallel.collect_exact_keys(plan, case, cat)
        self.assertGreater(len(keys), 15)
        self.assertEqual(cat.calls, 0)  # Collection never evaluates a real kernel.
        kernels = {key: source.energy_kernel(*key) for key in keys}
        cat.install_prefetched(kernels)
        with tempfile.TemporaryDirectory() as path:
            result = controls.spatial_case(plan, case, cat, lambda *args: None, Path(path))
        self.assertEqual(cat.calls, 0)
        self.assertGreater(cat.prefetch_hits, 0)
        serial = controls.CachedDirect(source)
        with tempfile.TemporaryDirectory() as path:
            expected = controls.spatial_case(plan, case, serial, lambda *args: None, Path(path))
        self.assertEqual(json.dumps(result, sort_keys=True), json.dumps(expected, sort_keys=True))
        self.assertGreater(serial.calls, 0)
        with self.assertRaisesRegex(RuntimeError, 'Unscheduled spectral key'):
            cat.energy_kernel(.123456789, .0123456789)

    def test_untraced_finite_differences_are_explicitly_rejected(self):
        with self.assertRaisesRegex(ValueError, 'finite_differences=false'):
            parallel.collect_exact_keys({}, dict(phase='spatial', finite_differences=True), None)


if __name__ == '__main__':
    unittest.main()
