"""Adjudicate the frozen R2 contract from exact, hash-linked evidence.

No tolerance is selected here. Failed or missing evidence prevents admission.
The external certificate never changes the measured catalogue or admits NbN.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import platform
import sys

import numpy as np
import scipy

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
BASE = HERE.parent
sys.path.insert(0, str(ROOT))
from pysnspd.experimental.energy_catalog import OccupationEnergyCatalog, energy_at_count_batch, retarded_spectrum_batch
from causal_reference import vacuum_reference


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def metrics(value, reference):
    value, reference = np.asarray(value), np.asarray(reference)
    difference = abs(value-reference)
    return {'scaled': (difference/np.maximum(1, abs(reference))).tolist(),
            'relative': [float(difference[i]/abs(reference[i])) if reference[i] != 0 else None for i in range(3)]}


def main():
    paths = {
        'criteria': BASE/'acceptance_criteria.json', 'catalog': BASE/'catalogs/occupation_catalog.npz',
        'assessment': HERE/'final_assessment.json', 'spectral': HERE/'spectral_crosscheck.json',
        'shape': HERE/'shape_final.json', 'regressions': HERE/'support_regression_points.json',
        'interpolation': BASE/'candidates/ratio_final_17/interpolation.json',
        'exact_count_rule': BASE/'candidates/ratio_final_17/count_rule_check.json',
        'separated_controls': BASE/'quadrature_diagnostics.json',
        'causal_reference': HERE/'causal_reference.json', 'thermal_reference': HERE/'thermal_reference.json',
        'tail_reference': HERE/'tail_reference.json', 'material': BASE/'material_comparison.json',
    }
    missing = [name for name, path in paths.items() if not path.is_file()]
    if missing:
        raise RuntimeError('Cannot adjudicate: missing '+', '.join(missing))
    records = {name:json.loads(path.read_text(encoding='utf-8')) for name, path in paths.items() if name != 'catalog'}
    hashes = {name:sha(path) for name, path in paths.items()}
    criteria, assessment = records['criteria'], records['assessment']
    gates = criteria['numerical_gates']
    limits = [gates['energy_max_scaled_error'], gates['amplitude_force_max_scaled_error'], gates['gamma_response_max_scaled_error']]
    relative_limit = gates['nonzero_force_and_gamma_max_relative_error']
    def passes(measurement, fraction=1.):
        return (all(value <= limit*fraction for value, limit in zip(measurement['scaled'], limits))
                and all(measurement['relative'][i] is None or measurement['relative'][i] <= relative_limit*fraction for i in (1, 2)))
    decisions = []
    def gate(name, passed, details):
        decisions.append(dict(gate=name, status='PASS' if passed else 'FAIL', details=details))
    gate('frozen_criteria', hashes['criteria'] == 'ef8164fe9f7d7d27f4ecbba1a318fc8642ea3dcc9723dd2e1a9d9c1447dc6bfb', criteria['criteria_id'])
    for name in ('assessment', 'spectral', 'interpolation', 'separated_controls', 'causal_reference', 'thermal_reference', 'tail_reference'):
        gate(name+'_criteria_link', records[name].get('criteria_sha256') == hashes['criteria'], hashes['criteria'])
    for name in ('assessment', 'shape', 'interpolation', 'exact_count_rule'):
        gate(name+'_catalog_link', records[name].get('catalog_sha256') == hashes['catalog'], hashes['catalog'])
    source = ROOT/'pysnspd/experimental/energy_catalog.py'
    for name, key in (('assessment', 'implementation_sha256'), ('spectral', 'implementation_sha256'), ('shape', 'source_sha256'),
                      ('interpolation', 'source_sha256'), ('exact_count_rule', 'source_sha256')):
        gate(name+'_current_source_link', records[name].get(key) == sha(source), sha(source))
    gate('total_causal_values_and_API', assessment['status'] == 'PASS_TESTED_PHYSICAL_AND_API_GATES'
         and assessment['counts']['physical_cases'] == 234 and not assessment['failed_physical_cases'] and not assessment['exceptions'],
         dict(counts=assessment['counts'], worst_errors=assessment['summary']))
    gate('historical_regressions', assessment.get('support_regression_points_sha256') == hashes['regressions']
         and assessment['historical_support_regression_points'] == len(records['regressions']['points']),
         len(records['regressions']['points']))
    gate('dense_support', not assessment['dense_cell_support']['failures']
         and assessment['dense_cell_support']['minimum_energy'] > 0 and assessment['dense_cell_support']['minimum_adjacent_spacing'] > 0,
         assessment['dense_cell_support'])
    shape = records['shape']
    gate('exact_ratio_cell_minima_at_sampled_amplitudes', shape['status'] == 'PASS' and shape['crossing_intervals'] == 0
         and shape['endpoint_nonpositive_increments'] == 0 and shape['minimum_cubic_increment'] > 0,
         {key:shape[key] for key in ('amplitudes', 'checked_cells', 'minimum_cubic_increment', 'flagged_Bernstein_intervals', 'coverage_limit')})
    spectral = records['spectral']
    gate('causal_spectrum_and_count_inverse', spectral['status'] == 'PASS' and len(spectral['rows']) == 90 and not spectral['failures'],
         {key:max(row[key] for row in spectral['rows']) for key in
          ('max_unsquared_scaled_residual', 'max_normalization_scaled_residual', 'max_count_inverse_scaled_error')})
    interpolation = records['interpolation']
    gate('fixed_population_field_interpolation', interpolation['status'] == 'PASS_SAMPLED' and len(interpolation['rows']) == 210
         and not interpolation['failures'] and all(passes({'scaled':r['scaled_error'], 'relative':r['relative_error']}) for r in interpolation['rows']),
         interpolation['worst_errors'])

    # Complete the field-error separation for FD without changing p between the
    # catalogue and direct spectrum. Entropy is identical and cancels here.
    table = OccupationEnergyCatalog.load(paths['catalog'])
    thermal_interpolation = []
    for reference in records['thermal_reference']['rows']:
        a, g, temperature = reference['amplitude'], reference['gamma'], reference['temperature']
        e_table = table.energy_kernel(a, g)[0]
        p = np.exp(-e_table/temperature)/(1+np.exp(-e_table/temperature))
        actual = table.evaluate(a, g, p)
        energy = energy_at_count_batch(table.count_nodes, delta=a, gamma=g, eta=table.eta)
        c, s = retarded_spectrum_batch(energy, delta=a, gamma=g, eta=table.eta)
        direct = vacuum_reference(a, g)+4*np.sum(np.array([energy, s.imag/c.real, -s.real*s.imag/c.real])*(p*table.count_weights)[None,:], axis=1)
        thermal_interpolation.append(dict(amplitude=a, gamma=g, temperature=temperature, candidate=list(actual),
                                          direct_fixed_p_reference=direct.tolist(), errors=metrics(actual, direct)))
    gate('FD_fixed_population_field_interpolation', all(passes(row['errors']) for row in thermal_interpolation), thermal_interpolation)
    exact_rule = records['exact_count_rule']
    gate('exact_saved_count_rule', len(exact_rule['rows']) == 159 and not exact_rule.get('failures', [])
         and all(passes({'scaled':r['scaled_error'], 'relative':r['relative_error']}) for r in exact_rule['rows']), exact_rule['worst_errors'])

    controls = records['separated_controls']
    index = {(row['eta'], row['order'], row['cutoff'], row['amplitude'], row['gamma'], row['profile']):row for row in controls['rows']}
    summaries, reference_budget = [], []
    for eta in (1e-7, 3e-8, 1e-8):
        for order in (8, 12, 16, 24, 32):
            rows = [row for row in controls['rows'] if row['eta'] == eta and row['order'] == order and row['cutoff'] == 12.]
            gate(f'count_resolution_{eta}_{order}', len(rows) == 159, len(rows))
            observations = [{'scaled':row['same_eta_order32_cutoff12_scaled_difference'], 'relative':row['same_eta_order32_cutoff12_relative_difference']} for row in rows]
            summaries.append(dict(eta=eta, order=order, reference_order=32,
                                  max_scaled_change=[max(obs['scaled'][i] for obs in observations) for i in range(3)],
                                  max_relative_change=[max(obs['relative'][i] or 0 for obs in observations) for i in range(3)]))
            if order == 24:
                reference_budget.extend(observations)
    gate('finite_eta_quadrature_reference_budget', not controls['failures'] and all(passes(obs, .1) for obs in reference_budget),
         dict(comparison='order24 versus32 at each of three eta values; all159 cases', resolutions=summaries))
    cutoff_rows = []
    for cutoff in (6., 8., 12., 20.):
        selected = [row for row in controls['rows'] if row['eta'] == 1e-8 and row['order'] == 32 and row['cutoff'] == cutoff]
        gate(f'cutoff_evidence_{cutoff}', len(selected) == 159, len(selected))
        for row in selected:
            reference = index[(1e-8, 32, 20., row['amplitude'], row['gamma'], row['profile'])]
            cutoff_rows.append(dict(amplitude=row['amplitude'], gamma=row['gamma'], profile=row['profile'], cutoff=cutoff,
                                    errors=metrics(row['values'], reference['values'])))
    gate('selected_cutoff_12', all(passes(row['errors']) for row in cutoff_rows if row['cutoff'] == 12.),
         dict(compared_cutoffs=[6, 8, 12, 20], selected=12,
              max_scaled_change_12_to_20=[max(row['errors']['scaled'][i] for row in cutoff_rows if row['cutoff']==12.) for i in range(3)],
              independent_tail_max_absolute=max(abs(value) for row in records['tail_reference']['rows'] for value in row['causal_tail_to_infinity'])))
    ideal_index = {(row['amplitude'],row['gamma'],row['profile']):row['reference'] for row in records['causal_reference']['rows']}
    ideal_index.update({(row['amplitude'],row['gamma'],'FD_'+str(row['temperature'])):row['causal_reference']
                       for row in records['thermal_reference']['rows']})
    tail_budgets = []
    for row in records['tail_reference']['rows']:
        reference = ideal_index[(row['amplitude'],row['gamma'],row['profile'])]
        bounds = [abs(value)+row['error_estimate'] for value in row['causal_tail_to_infinity']]
        tail_budgets.append(dict(scaled=[bound/max(1,abs(ref)) for bound,ref in zip(bounds,reference)],
                                relative=[bound/abs(ref) if ref!=0 else None for bound,ref in zip(bounds,reference)]))
    gate('tail_to_infinity', len(tail_budgets)==318 and all(passes(obs, .1) for obs in tail_budgets),
         '318 independent adaptive tail integrals at cutoffs12 and20 meet the frozen10% reference-error budget for each observable; combined with same-eta cutoff comparison.')
    eta_summaries = []
    for eta in (1e-7, 3e-8, 1e-8):
        rows = [row for row in controls['rows'] if row['eta']==eta and row['order']==32 and row['cutoff']==12.]
        eta_summaries.append(dict(eta=eta, max_scaled_bias=[max(row['causal_scaled_error'][i] for row in rows) for i in range(3)],
                                 max_relative_bias=[max(row['causal_relative_error'][i] or 0 for row in rows) for i in range(3)]))
    selected_eta = [row for row in controls['rows'] if row['eta']==table.eta and row['order']==32 and row['cutoff']==12.]
    gate('selected_regulator_causal_limit', len(selected_eta)==159 and
         all(passes({'scaled':row['causal_scaled_error'], 'relative':row['causal_relative_error']}) for row in selected_eta),
         dict(selected_eta=table.eta, sequence=eta_summaries,
              interpretation='The coarser eta=1e-7 exceeds the Gamma gate; it is retained as convergence evidence, not an admitted regulator.'))
    gate('exact_rule_reference_link', exact_rule['reference_sha256'] == hashes['separated_controls'], hashes['separated_controls'])

    passed = all(row['status']=='PASS' for row in decisions)
    material = records['material']['experimental']
    certificate = dict(schema='pysnspd.stage1_r2.catalog_admission.v1', created_utc=datetime.now(timezone.utc).isoformat(),
                       numerical_status='PASS_SAMPLED_ELECTRONIC_DOMAIN' if passed else 'FAIL',
                       material_status=material['status'], combined_stage1_status='ELECTRONIC_CATALOGUE_ONLY; MATERIAL_NOT_ADMITTED',
                       next_permitted_validation='Synthetic one-or-two-cell validation within the tested contract' if passed else 'No advancement; resolve failed gates',
                       production_solver_connected=False,
                       physical_scope='Uniform electronic catalogue at positive amplitude plus separate exact normal point. Does not resolve the model0.4 core/principal-symbol issues, predict detector transients, or certify NbN phonon inputs.',
                       sampling_scope=dict(amplitude=[.08,1.5], physical_Gamma=[0,1.2], eta=table.eta,
                                           count_nodes=len(table.count_nodes), count_support=float(sum(table.count_weights)),
                                           populations='Five fixed profiles and nine thermal checks from frozen criteria;234 physical cases including off-node points.',
                                           nonuniform_limit='Amplitude remains sampled. No uniform error bound for arbitrary occupations or every interior field is asserted.'),
                       catalog_internal_status=table.vacuum.metadata.get('admission_status'),
                       external_certificate_policy='Admission attaches to the exact catalogue SHA and query-code SHA; the measured NPZ is not rewritten to change its builder status.',
                       evidence={name:dict(path=path.relative_to(ROOT).as_posix(),sha256=hashes[name]) for name,path in paths.items()},
                       query_code_sha256=sha(source), adjudicator_code_sha256=sha(Path(__file__)),
                       environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,platform=platform.platform()),
                       gates=decisions, material_reason='Published provenance is verified, but source units/basis remain unverified,361 DOS values are negative and4 duplicate abscissae conflict. No clipping or renormalization admits the file.',
                       final_total_errors=assessment['summary'])
    (BASE/'catalog_admission.json').write_text(json.dumps(certificate,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(numerical_status=certificate['numerical_status'],material_status=certificate['material_status'],
                         gates=len(decisions),failed_gates=[row for row in decisions if row['status']!='PASS']),indent=2))


if __name__=='__main__':
    main()
