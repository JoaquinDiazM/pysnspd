"""Scheduling and monitoring protections for the eight-case calibration campaign."""
import argparse
from contextlib import redirect_stdout
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('_eight_case_test', ROOT/'pipelines/15_calibration_basecases.py')
C = importlib.util.module_from_spec(spec)
spec.loader.exec_module(C)


class CalibrationCampaignTests(unittest.TestCase):
    def inventory(self):
        return dict(logical_cpus=list(range(32)), physical_core_groups=[[i, i+16] for i in range(16)],
                    available_memory_bytes=100*1024**3, total_memory_bytes=128*1024**3, cpu_quota_units=None)

    def test_all_eight_admitted_without_logical_overlap_or_reserved_core_use(self):
        result = C.allocation(self.inventory(), [])
        slots = result['slots']
        cpus = [cpu for slot in slots.values() for cpu in slot]
        self.assertEqual(len(slots), 8)
        self.assertEqual(len(cpus), len(set(cpus)))
        self.assertGreaterEqual(len(cpus)/32, .8)
        self.assertLessEqual((len(cpus)+len(result['coordinator_cpus']))/32, .9)
        reserved = {cpu for group in result['reserved_whole_core_groups'] for cpu in group}
        self.assertGreaterEqual(len(reserved), 4)
        self.assertFalse((set(cpus)|set(result['coordinator_cpus'])) & reserved)
        self.assertTrue(result['shared_SMT_cores'])

    def test_busy_or_quota_or_memory_constrained_host_refuses_all_before_start(self):
        for inventory, active in [(self.inventory(), [dict(pid=99, affinity_cpus=[0,16])]),
                                  (dict(self.inventory(), cpu_quota_units=24), []),
                                  (dict(self.inventory(), available_memory_bytes=40*1024**3), [])]:
            with self.assertRaises(ValueError):
                C.allocation(inventory, active)

    def test_complete_factorial_uses_one_common_physics_and_detection_horizon(self):
        cases = C.catalog()
        self.assertEqual({(c['sheet_resistance_ohm'], c['current_uA'], c['wavelength_nm']) for c in cases},
                         {(r,i,w) for r in (425,608) for i in (15.5,21.5) for w in (775,1550)})
        slots = C.allocation(self.inventory(), [])['slots']
        for c in cases:
            words = C.command(c, Path('/unused'), slots[c['name']], argparse.Namespace(ceiling_ps=100.))
            value = lambda flag: words[words.index(flag)+1]
            self.assertIn('--latency-only', words)
            self.assertEqual(value('--retained-fraction'), '0.9')
            self.assertEqual(value('--tau-ee-Tc-ps'), '5')
            self.assertEqual(value('--tau-escape-ps'), '20')
            self.assertEqual(float(value('--target-depairing-uA')), 38.8*80/120)
            self.assertEqual(value('--length-nm'), '300')
            self.assertEqual(value('--extra-time-ps'), '2')
            self.assertEqual(len(value('--cpu-list').split(',')), c['threads'])
            self.assertIn('13_dimensional_photon_run.py', ' '.join(words))

    def test_board_handles_silent_preparation_and_preserves_readable_error(self):
        with tempfile.TemporaryDirectory() as directory:
            board = C.Board(C.catalog(), Path(directory), 100., False)
            name = C.catalog()[0]['name']
            dest = Path(directory)/'runs'/name
            dest.mkdir(parents=True)
            (dest/'manifest.json').write_text(json.dumps(dict(status='PREPARATION')))
            board.ingest(name, 'Some non-JSON compilation message')
            output = io.StringIO()
            with redirect_stdout(output):
                board.render()
            self.assertIn('DC/source', output.getvalue())
            self.assertIn('Vout[mV]', output.getvalue())
            self.assertNotIn('"status"', output.getvalue())
            board.ingest(name, json.dumps(dict(status='RUNNING', time_ps=20., Vout_mV=.05, accepted_steps=500)))
            board.update(name, status='FAILED_EXIT', exit_code=1, last_message='See full raw traceback in case log')
            output = io.StringIO()
            with redirect_stdout(output):
                board.render(final=True)
            self.assertIn('FAILED', output.getvalue())
            self.assertIn('See full raw traceback', output.getvalue())


if __name__ == '__main__':
    unittest.main()
