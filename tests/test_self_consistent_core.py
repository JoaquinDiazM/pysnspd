"""Variational checks of the exact thermal gap block and accepted checkpoint."""
import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('self_consistent_core', ROOT/'sandbox/stage4_core/run_self_consistent_core.py')
core = importlib.util.module_from_spec(spec)
spec.loader.exec_module(core)
thermal = core.thermal


class TestSelfConsistentCore(unittest.TestCase):
    def setUp(self):
        self.graph = thermal.rectangular_graph(5, 5, 4., 4.)
        self.xy = self.graph.coordinates_bar.copy()
        self.xy[:, 0] -= 2
        self.radius = np.linalg.norm(self.xy, axis=1)
        self.d = (.7+.1*np.cos(self.xy[:, 0]))*np.exp(.2j*self.xy[:, 1])
        self.t, self.count = .5, 4

    def action(self, d, u):
        energy = np.dot(self.graph.area_weights, abs(d)**2)*np.log(self.t)
        for n, mode in enumerate(u):
            eps = 2*np.pi*self.t*(n+.5)
            energy += 2*np.pi*self.t*(thermal.spectral_energy_gradient(self.graph, d, eps, mode).energy
                +np.dot(self.graph.area_weights, abs(d)**2)/eps)
        return energy

    def test_gap_block_is_exact_minimum_of_full_action(self):
        u = [self.d/(2*np.pi*self.t*(n+.5)) for n in range(self.count)]
        f_sum = sum(mode/np.sqrt(1+abs(mode)**2) for mode in u)
        next_d, block = core.exact_gap_block(self.graph, self.d, f_sum, self.t, self.count, self.graph.boundary_nodes)
        observed = self.action(next_d, u)-self.action(self.d, u)
        self.assertAlmostEqual(observed, block['predicted_energy_change'], delta=2e-13)
        free = np.ones(len(self.d), bool); free[self.graph.boundary_nodes] = False
        gradient = 2*self.graph.area_weights*(block['coefficient']*next_d-2*np.pi*self.t*f_sum)
        np.testing.assert_allclose(gradient[free], 0., atol=1e-15)
        np.testing.assert_array_equal(next_d[~free], self.d[~free])
        self.assertLess(observed, 0.)

    def test_two_blocks_descend_and_checkpoint_keeps_solved_pair(self):
        solutions = [thermal.solve_frequency(self.graph, self.d, 2*np.pi*self.t*(n+.5), tol=1e-9)
            for n in range(self.count)]
        observation = thermal.evaluate_thermal(self.graph, self.d, self.t, solutions)
        next_d, block = core.exact_gap_block(self.graph, self.d, sum(s.f for s in solutions), self.t, self.count, self.graph.boundary_nodes)
        new = [thermal.solve_frequency(self.graph, next_d, solution.epsilon, initial_u=solution.u, tol=1e-9)
            for solution in solutions]
        self.assertLessEqual(thermal.evaluate_thermal(self.graph, next_d, self.t, new).energy,
            observation.energy+block['observed_energy_change']+1e-12)
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory); (out/'checkpoints').mkdir()
            state = dict(geometry=(self.graph, self.xy), d=self.d, solutions=solutions)
            record = core.save_checkpoint(out, 'control', 1, state, observation, next_d, {})
            with np.load(record['fields_path']) as fields:
                np.testing.assert_array_equal(fields['d'], self.d)
                np.testing.assert_array_equal(fields['next_d'], next_d)
                np.testing.assert_array_equal(fields['u'][0], solutions[0].u)
                self.assertGreater(np.max(abs(fields['d']-fields['next_d'])), 0.)
            metrics_only = core.save_checkpoint(out, 'control', 2, state, observation, next_d, {}, full_fields=False)
            self.assertFalse(metrics_only['full_checkpoint'])
            self.assertIsNone(metrics_only['fields_path'])
            self.assertTrue((out/'checkpoints/control_sweep0002.json').exists())
            self.assertFalse((out/'checkpoints/control_sweep0002.npz').exists())

    def test_rms_and_nonconvex_cutoff_rejection(self):
        next_d = self.d.copy(); free = np.ones(len(self.d), bool); free[self.graph.boundary_nodes] = False
        next_d[free] += .002
        metrics = core.residual_metrics(self.graph, self.radius, self.d, next_d, 2., 1.5)
        self.assertAlmostEqual(metrics['core_mass_rms_relative'], .001)
        self.assertAlmostEqual(metrics['global_free_maximum_relative'], .001)
        with self.assertRaises(ValueError):
            core.gap_coefficient(.001, 1)

    def test_storage_projection_counts_cadence_and_bounds_pilot(self):
        plan = dict(checkpoint_cadence=5, max_outer_iterations=150,
            cases=[dict(id='small', nx=5, ny=5, matsubara_count=4)])
        with tempfile.TemporaryDirectory() as directory:
            with patch.object(core.shutil, 'disk_usage', return_value=SimpleNamespace(free=10**9)):
                full = core.checkpoint_storage_budget(plan, Path(directory)/'new')
                pilot = core.checkpoint_storage_budget(plan, Path(directory)/'new', pilot_sweeps=2)
            self.assertEqual(full['cases']['small']['projected_full_checkpoints'], 32)
            self.assertEqual(pilot['cases']['small']['projected_full_checkpoints'], 2)
            self.assertEqual(full['maximum_recomputed_sweeps'], 4)
            with patch.object(core.shutil, 'disk_usage', return_value=SimpleNamespace(free=1)):
                with self.assertRaises(RuntimeError):
                    core.checkpoint_storage_budget(plan, Path(directory)/'new')


if __name__ == '__main__':
    unittest.main()
