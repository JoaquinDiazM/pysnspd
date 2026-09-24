"""Verify the final KWT delivery without tests, spectra or time evolution.

The same command works in Windows and Geminga. Exact bytes are required for
experimental sources, data and figures. Checkout newline differences are
accepted only for production files whose Git content is unchanged.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import runpy
import sys


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT/'docs/implementation/stage4/final_kwt_20260924'
OLD = ROOT/'docs/implementation/stage4/practical_time_review_20260924'
BASE = runpy.run_path(str(ROOT/'sandbox/stage4_core/verify_time_review.py'))
PREVIOUS = runpy.run_path(str(ROOT/'sandbox/stage4_core/verify_practical_time_review.py'))
read, sha, require = BASE['read'], BASE['sha'], BASE['require']
inherited_sources = PREVIOUS['inherited_sources']


def same_record(left, right):
    """Allow only numerical reduction roundoff across Windows/Linux BLAS."""
    if isinstance(left, bool) or isinstance(right, bool):
        return left is right
    if isinstance(left, (float, int)) and isinstance(right, (float, int)):
        return math.isclose(left, right, rel_tol=2e-12, abs_tol=1e-15)
    if isinstance(left, dict) and isinstance(right, dict):
        return left.keys() == right.keys() and all(same_record(left[key], right[key]) for key in left)
    if isinstance(left, list) and isinstance(right, list):
        return len(left) == len(right) and all(same_record(a, b) for a, b in zip(left, right))
    return left == right


def file_map(folder, mapping):
    for name, expected in mapping.items():
        path = Path(name.replace('\\', '/'))
        require(not path.is_absolute() and '..' not in path.parts, 'Unsafe relative file path')
        require(sha(folder/path) == expected, 'Changed recorded file: '+str(folder/path))
    return len(mapping)


def verify_budget(identity):
    resource, budget = identity['resources'], identity['budget']
    selected = set(budget['selected_affinity_cpus'])
    available = set(resource['logical_cpus'])
    worker_cpus = budget['worker_affinity_cpus']
    require(selected <= available, 'CPU affinity escaped the observed available set')
    require(len(worker_cpus) == budget['workers'] and len(set(worker_cpus)) == len(worker_cpus),
            'Worker CPU allocation is inconsistent')
    require(set(worker_cpus) <= selected and budget['coordinator_cpu'] in selected,
            'Worker/coordinator outside admitted affinity')
    require(budget['coordinator_cpu'] not in worker_cpus, 'Coordinator shares its logical CPU with a worker')
    require(budget['threads_per_process'] == 1, 'Numerical thread budget changed')
    require(budget['total_processes'] == budget['workers']+1, 'Nested process count is inconsistent')
    require(budget['total_processes'] <= .9*budget['effective_cpu_quota_units'], 'CPU budget exceeds 90%')
    require(budget['estimated_memory_reservation_bytes'] <= .9*resource['available_memory_bytes'],
            'Memory reservation exceeds 90%')
    whole_reserved = sum(not selected.intersection(group) for group in resource['physical_core_groups'])
    require(whole_reserved == budget['physical_cores_reserved'], 'Whole-core reservation differs')
    require(whole_reserved >= 2, 'Expected two complete reserved physical cores')
    return dict(workers=budget['workers'], coordinator=1, threads_per_process=1,
        assigned_logical_cpus=len(selected), entirely_reserved_physical_cores=whole_reserved,
        reservation_not_peak_RSS=True)


def verify_policy():
    policy = read(DATA/'acceptance_policy.json')
    analysis = read(OLD/'analysis.json')
    original = read(OLD/'raw/summary.json')
    require(sha(ROOT/policy['analysis_path']) == policy['analysis_sha256'], 'ETD2 source analysis changed')
    require(original['status'] == 'TEMPORAL_REFINEMENT_NOT_MET', 'Original ETD2 failure was rewritten')
    require(sum(row['admitted'] for row in original['refinement']['records']) == 95,
            'Original ETD2 certificate does not retain 95/96')
    require(not policy['original_certificate']['changed'], 'Original ETD2 gate changed')
    revised = policy['revised_policy']
    require(revised['relative_limit'] == .02 and revised['absolute_addition'] == 0.,
            'Different practical ETD2 policy')
    require(revised['posterior_policy_revision_declared'] and revised['same_channel_scale']
            and revised['no_remaining_tail_denominator'], 'Practical policy is not explicit')
    rows = analysis['comparisons']
    identical_times = analysis['temporal_refinement_scope']['bitwise_identical_observation_times_ps']
    admitted = [row['error_L2'] <= .02*row['initial_norm'] for row in rows]
    independent = [passed for row, passed in zip(rows, admitted) if row['time_ps'] not in identical_times]
    require(sum(admitted) == policy['passed_comparisons'] == len(rows) == 96,
            'Revised policy does not reproduce 96/96')
    require(sum(independent) == policy['passed_independent_comparisons'] == len(independent) == 24,
            'Revised policy does not reproduce 24/24 distinct-trajectory comparisons')
    require(not policy['tail_rerun_pending'] and not policy['precise_tail_relaxation_time_claimed'],
            'Old secondary tail is incorrectly left pending or overstated')
    require(not policy['stage4_complete'] and not policy['production_changed'], 'ETD2 scope overstated')
    return dict(original_certificate_passed=95, original_certificate_total=96,
        practical_policy_passed=96, practical_distinct_trajectory_passed=24,
        original_failure_preserved=True, secondary_tail_rerun_pending=False)


def verify_pilots():
    comparison = read(DATA/'pilot_comparison.json')
    file_map(DATA, comparison['inputs_sha256'])
    for folder_name, plan_name in [('pilot_initial', 'dual_kwt_plan.json'),
                                   ('pilot_predicted', 'dual_kwt_predicted_plan.json')]:
        folder = DATA/folder_name
        identity, plan, summary = [read(folder/name) for name in
                                  ('identity.json', 'executed_plan.json', 'summary.json')]
        require(sha(folder/'executed_plan.json') == identity['plan_sha256'] == sha(DATA/plan_name),
                'Pilot plan changed: '+folder_name)
        inherited_sources(plan['source_sha256'])
        require(summary['status'] == 'BOUNDED_STEP_PILOT_COMPLETE' and summary['covered_ps'] == .001,
                'Pilot status/horizon changed')
        require(summary['maximum_spectral_residual'] <= plan['spectral_tolerance'], 'Pilot spectral residual not admitted')
        verify_budget(identity)
    for name, item in comparison['source_verification'].items():
        source = ROOT/('sandbox/stage4_core' if name == 'dual_kwt_time.py' else 'tests')/name
        archive = DATA/'pilot_initial/executed_sources'/name
        require(sha(source) == sha(archive) == item['registered_sha256'] == item['current_sha256'],
                'Initial runner or test was overwritten')
    require(comparison['physics_solves_executed_by_analysis'] == 0, 'Pilot comparison executed physical work')
    return dict(pilots=2, completed_ps=.001, physical_model_and_time_integrator_unchanged=True)


def verify_thermal():
    folder = DATA/'thermal'
    summary, identity, plan, analysis, extraction = [read(folder/name) for name in
        ('summary.json', 'identity.json', 'executed_plan.json', 'analysis.json', 'extraction_receipt.json')]
    require(extraction['status'] == 'COMPLETED_RESULTS_EXTRACTED' and extraction['physical_solves'] == 0,
            'Incomplete or non-read-only thermal extraction')
    count = file_map(folder, extraction['source_hashes'])
    require(sha(folder/'executed_plan.json') == identity['plan_sha256'] == sha(DATA/'dual_kwt_predicted_plan.json'),
            'Full trajectory plan differs from the verified predictor plan')
    inherited_sources(plan['source_sha256'])
    require(summary['status'] == analysis['status'] == 'DUAL_KWT_THERMAL_COMPLETE', 'Full thermal trajectory incomplete')
    require(summary['completed_horizon_ps'] == analysis['completed_horizon_ps'] == 1., 'Wrong thermal horizon')
    require(summary['all_practical_criteria_met'] and summary['refinement_all_admitted'], 'Thermal margins not admitted')
    require(summary['maximum_spectral_residual'] <= plan['spectral_tolerance'] == 1e-7, 'Changed or exceeded spectral tolerance')
    require(plan['refinement_relative_limit'] == plan['energy_relative_limit'] == .02, 'Changed practical tolerance')
    require(summary['time_order'] == 1 and not summary['production_changed'] and not summary['photon'],
            'Thermal method/scope was overstated')
    require(not summary['stage4_complete'] and not summary['nonthermal_work_closed'], 'Thermal trajectory overstates full stage4')
    require(sha(ROOT/'sandbox/stage4_core/analyze_dual_kwt_results.py') == analysis['analysis_source_sha256'],
            'Thermal analysis source changed')
    analyzer = runpy.run_path(str(ROOT/'sandbox/stage4_core/analyze_dual_kwt_results.py'))
    reproduced = analyzer['analyze'](folder, ROOT)[0]
    for key, value in reproduced.items():
        require(same_record(analysis[key], value), 'Reproduced thermal analysis differs beyond reduction roundoff: '+key)
    require(analysis['admitted_comparisons'] == analysis['total_comparisons'] == 48
            and analysis['admitted_positive_time_comparisons'] == analysis['positive_time_comparisons'] == 40,
            'Incorrect thermal comparison counts')
    file_map(folder/'figures', analysis['figures_sha256'])
    qa = read(folder/'figure_qa.json')
    require(qa['status'] == 'PASS' and qa['figures_sha256'] == analysis['figures_sha256'], 'Thermal figure QA mismatch')
    require(qa['physics_solves'] == 0, 'Figure QA scope changed')
    remote_count = None
    remote = Path(extraction['source'])
    if sys.platform.startswith('linux') and remote.is_dir():
        file_map(remote, extraction['source_hashes'])
        remote_count = len(list(remote.glob('checkpoint_*.npz')))
        require(remote_count == extraction['checkpoints_retained_remotely'] == 581,
                'Remote accepted checkpoints are missing')
    return dict(copied_files_verified=count, accepted_steps=summary['accepted_steps'],
        remote_checkpoint_count=remote_count, horizon_ps=1., temporal_comparisons=48,
        positive_time_comparisons=40, maximum_spectral_residual=summary['maximum_spectral_residual'],
        data_reproduced_without_physics=True, resources=verify_budget(identity))


def verify_longitudinal():
    folder = DATA/'longitudinal'
    receipt = read(folder/'receipt.json')
    inherited_sources(receipt['source_sha256'])
    require(receipt['status'] == 'PASSED_CONTROL' and receipt['duration_ps'] == 1., 'Weak longitudinal control incomplete')
    require(receipt['nodes'] == 1712 and receipt['matsubara_count'] == 256, 'Different longitudinal graph/reference')
    require(1.8 < receipt['error_reduction_ratio'] < 2.2, 'Weak longitudinal Euler order inconsistent')
    for name in ('primary', 'refined'):
        result = receipt[name]
        require(result['metric_error_over_initial'] <= .02 and abs(result['integrated_balance_over_initial']) <= .02,
                'Weak longitudinal accuracy not admitted')
        require(0 <= result['minimum_occupation'] <= result['maximum_occupation'] <= 1,
                'Weak longitudinal Pauli support violated')
    require(receipt['full_operator_embedding_relative_error'] < 1e-10,
            'Longitudinal invariant mode does not reproduce its full spatial operator')
    require('Internal energy total' in receipt['excluded'] and 'Photon or hotbelt' in receipt['excluded'],
            'Weak diagnostic scope lost')
    require(not receipt['eta_is_physical_bath'], 'Diagnostic eta reinterpreted as a physical bath')
    return dict(status=receipt['status'], duration_ps=receipt['duration_ps'],
        weak_invariant_mode_only=True, internal_energy_total_validated=False)


def verify_tests():
    receipt = read(DATA/'unit_tests_receipt.json')
    require(receipt['status'] == 'PASS', 'Focused test receipt not passing')
    total = 0
    for run in receipt['runs']:
        require(run['status'] == 'PASS', 'A focused test run failed')
        inherited_sources(run['source_sha256'])
        total += run['tests']
    require(total == receipt['tests'] == 6, 'Missing focused test records')
    return dict(recorded_tests=total, tests_executed_now=0)


def verify_pdf(required=False):
    path = DATA/'report_qa.json'
    if not path.exists():
        require(not required, 'PDF QA required for completed delivery')
        return dict(status='NOT_YET_PRESENT')
    qa, context = read(path), read(DATA/'report_context.json')
    require(qa['status'] == qa['semantic_checks'] == 'PASS' and qa['all_pages_visually_inspected'], 'PDF QA incomplete')
    require(qa['pdf_sha256'] == context['pdf_sha256'] == sha(ROOT/qa['pdf_path']), 'PDF bytes differ from QA/context')
    require(qa['pages'] == context['pages'] == 3, 'Unexpected PDF page count')
    file_map(ROOT, context['source_sha256'])
    require(context['builder_sha256'] == sha(ROOT/'sandbox/stage4_core/build_final_kwt_report.py'), 'PDF builder changed')
    require(context['physics_solves_performed'] == 0 and not context['stage4_complete']
            and not context['production_changed'], 'PDF scope overstated')
    available_renderings = 0
    for row in qa['renderings']:
        require(row['visual_status'] == 'PASS', 'PDF page not visually admitted')
        path = ROOT/row['rendered_path']
        if path.exists():
            require(sha(path) == row['sha256'], 'Changed available page rendering')
            available_renderings += 1
    return dict(status='PASS', pages=qa['pages'], available_renderings_verified=available_renderings,
        pdf_sha256=qa['pdf_sha256'])


def manifest(required=False):
    path = DATA/'delivery_manifest.json'
    if not path.exists():
        require(not required, 'Delivery manifest is required')
        return dict(present=False, files_verified=0)
    data = read(path)
    rows = data['files']
    require(all(row['path'].replace('\\', '/') != path.relative_to(ROOT).as_posix() for row in rows),
            'Manifest cannot hash itself')
    require(not data.get('production_changed', False) and not data.get('stage4_complete', False),
            'Manifest overstates production/full-stage scope')
    count = BASE['records'](rows)
    return dict(present=True, files_verified=count)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--require-manifest', action='store_true')
    args = parser.parse_args()
    result = dict(verified=True, policy=verify_policy(), pilots=verify_pilots(),
        thermal=verify_thermal(), longitudinal=verify_longitudinal(), focused_tests=verify_tests(),
        pdf=verify_pdf(required=args.require_manifest), manifest=manifest(required=args.require_manifest),
        production=BASE['production_contract'](), physics_solves_executed_now=0,
        trajectories_executed_now=0, tests_executed_now=0, stage4_complete=False)
    external = Path('/home/jdiaz/GEMINGA_COMMANDS.md')
    if sys.platform.startswith('linux') and external.exists():
        require(sha(external) == sha(ROOT/'docs/GEMINGA_COMMANDS.md'), 'External command notebook differs')
        result['external_command_notebook_matches'] = True
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
