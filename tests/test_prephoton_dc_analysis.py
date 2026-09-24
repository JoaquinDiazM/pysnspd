"""Admission regressions for missing data and near-zero DC observations."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from sandbox.stage5_prephoton.analyze_prephoton_dc import analyze


def write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data), encoding='utf8')


def fixture(root):
    names = ('base', 'mesh', 'long')
    plan = dict(references=[dict(id=name) for name in names],
        holds=[dict(id=name, reference=name) for name in names]+[dict(id='half', reference='base')],
        comparisons=dict(mesh=['base', 'mesh'], length=['base', 'long'], time=['base', 'half']),
        acceptance=dict(gap_comparison_relative=.01, current_comparison_relative=.02,
                        normalized_dark_readout_difference=.001),
        hold_template=dict(delta0_over_kBTc=1.764))
    write(root/'workflow_plan.json', plan)
    parameters = dict(T_K=.9, Tc_K=8.65, diffusion_m2_s=5e-5, sheet_resistance_ohm=608.,
                      width_m=80e-9, current_A=21.5e-6, added_inductance_H=96e-9)
    for name in names:
        directory = root/'references/cases'/name
        xy = np.array([[0., 0.], [0., 1.], [1., 0.], [1., 1.], [2., 0.], [2., 1.]])
        if name == 'long':
            xy[:, 0] *= 2
        gap = np.full(6, 1.5, complex)
        npz = dict(delta_bar=gap, coordinates_bar=xy, area_weights=np.ones(6),
                   boundary_nodes=np.array([0, 1, 4, 5]))
        resolved = 2e-10 if name == 'long' else 1e-10
        summary = dict(status='DC_REFERENCE_ADMITTED', stationary=True, admitted_for_dark_hold=True,
            plan=parameters, bulk=dict(gap_bulk_bar=1.5, current_target_A=21.5e-6,
                current_reference_A=21.5e-6), fixed_external_inductance_H=1e-7-resolved,
            resolved_equilibrium_inductance_H=resolved, total_equilibrium_inductance_H=1e-7)
        write(directory/'summary.json', summary)
        np.savez(directory/'reference.npz', **npz)
        for hold in [name]+(['half'] if name == 'base' else []):
            destination = root/'holds/cases'/hold
            circuit = dict(R_load_ohm=50., R_bias_ohm=1e4, L_bias_H=1e-6, C_couple_F=1e-10,
                           Lk_ext_H=summary['fixed_external_inductance_H'])
            times = [0., .5, 1.] if hold == 'half' else [0., 1.]
            history = [dict(time_ps=t, gap_drift_relative=1e-8*t, current_cut_mean_A=21.5e-6,
                            V_out_V=(1e-12 if hold == 'half' else 1e-20)*t,
                            terminal_voltage_V=5e-6) for t in times]
            write(destination/'summary.json', dict(status='COMPLETED', dc_hold_passed=True,
                checks=dict(stationary=True), duration_ps=1., dt_ps=.5 if hold == 'half' else 1.,
                initial_discrete_current_A=21.5e-6, requested_current_A=21.5e-6, circuit=circuit))
            write(destination/'history.json', history)
            np.savez(destination/'fields.npz', gap=np.asarray([gap for _ in times]), initial_gap=gap,
                     time_ps=times, area_weights=np.ones(6), coordinates_bar=xy)


class PrephotonAnalysisTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        fixture(self.root)
        self.no_plots = patch('sandbox.stage5_prephoton.analyze_prephoton_dc._plots', return_value=[])
        self.no_plots.start()

    def tearDown(self):
        self.no_plots.stop()
        self.temp.cleanup()

    def test_tiny_dark_voltage_uses_pulse_scale(self):
        result = analyze(self.root)
        self.assertTrue(result['prephoton_dc_admitted'])
        voltage = result['comparisons']['time']['measures']['trajectory_dark_voltage_difference_over_Rload_target']
        self.assertLess(voltage, 1e-8)
        self.assertTrue((self.root/'analysis.md').exists())

    def test_every_step_maxima_override_sparse_history(self):
        path = self.root/'holds/cases/half/summary.json'
        summary = json.loads(path.read_text())
        summary['maxima'] = dict(gap_drift_relative=.005, dark_readout_absolute_V=2e-8,
                                device_voltage_absolute_V=9e-6)
        summary['history_sampling'] = 'Stored observations only; maxima evaluate every accepted step'
        write(path, summary)
        record = analyze(self.root)['holds']['half']
        self.assertEqual(record['maximum_gap_drift_relative'], .005)
        self.assertEqual(record['maximum_absolute_dark_readout_V'], 2e-8)
        self.assertEqual(record['maximum_absolute_device_voltage_V'], 9e-6)
        self.assertTrue(all(value == 'SUMMARY_MAXIMUM_OVER_EVERY_ACCEPTED_STEP'
                            for value in record['maximum_sources'].values()))

    def test_sampled_fallback_is_explicit(self):
        record = analyze(self.root)['holds']['half']
        self.assertTrue(all(value == 'SAMPLED_HISTORY_ONLY_FALLBACK'
                            for value in record['maximum_sources'].values()))

    def test_missing_hold_cannot_be_admitted(self):
        (self.root/'holds/cases/half/summary.json').unlink()
        result = analyze(self.root)
        self.assertFalse(result['prephoton_dc_admitted'])
        self.assertTrue(result['holds']['base']['admitted'])
        self.assertFalse(result['holds']['half']['available'])

    def test_failed_manifest_overrules_stale_success(self):
        write(self.root/'holds/cases/half/manifest.json', dict(status='FAILED'))
        self.assertFalse(analyze(self.root)['prephoton_dc_admitted'])

    def test_equal_overlap_is_not_equal_horizon(self):
        path = self.root/'holds/cases/half/summary.json'
        summary = json.loads(path.read_text()); summary['duration_ps'] = 2.
        write(path, summary)
        result = analyze(self.root)
        self.assertFalse(result['prephoton_dc_admitted'])
        self.assertIn('horizon', result['holds']['half']['reason'])

    def test_nonfinite_field_is_rejected(self):
        path = self.root/'holds/cases/half/fields.npz'
        with np.load(path) as saved:
            arrays = dict(saved)
        arrays['gap'][1, 2] = np.nan
        np.savez(path, **arrays)
        result = analyze(self.root)
        self.assertFalse(result['prephoton_dc_admitted'])
        self.assertIn('Nonfinite', result['holds']['half']['reason'])

    def test_changed_physical_current_cannot_pass_mesh_comparison(self):
        path = self.root/'references/cases/mesh/summary.json'
        summary = json.loads(path.read_text()); summary['plan']['current_A'] *= 1.1
        write(path, summary)
        result = analyze(self.root)
        self.assertFalse(result['comparisons']['mesh']['admitted'])
        self.assertIn('current_A', result['comparisons']['mesh']['reason'])


if __name__ == '__main__':
    unittest.main()
