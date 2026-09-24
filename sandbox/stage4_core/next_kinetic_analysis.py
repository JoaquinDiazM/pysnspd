"""Hash-verified compact review of the frozen-spectrum Keldysh diagnostic."""
from pathlib import Path
import argparse
import hashlib
import json
import shutil

ROOT = Path(__file__).resolve().parents[2]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    identity = json.loads((args.raw/'identity.json').read_text())
    summary = json.loads((args.raw/'summary.json').read_text())
    plan_path = ROOT/'docs/implementation/stage4/self_consistent_review_20260924/next_kinetic_plan.json'
    if sha(plan_path) != identity['plan_sha256']:
        raise ValueError('Executed kinetic plan changed')
    for name, digest in {**identity['sources'], **identity['inputs']}.items():
        if sha(ROOT/name) != digest:
            raise ValueError('Executed kinetic source/input changed: '+name)
    if summary['status'] != 'FROZEN_KINETIC_RESPONSE_COMPLETE' or len(summary['records']) != 72 or len(summary['integrated']) != 12:
        raise ValueError('The complete registered frozen-spectrum diagnostic is required')
    if args.output.exists():
        raise FileExistsError('Fresh compact-evidence directory required')
    args.output.mkdir(parents=True)
    for name in ('identity.json', 'summary.json', 'progress.jsonl'):
        shutil.copyfile(args.raw/name, args.output/name)
    raw_integrity, representative = [], []
    for row in summary['records']+list(summary['integrated'].values()):
        path = args.raw/row['fields_path']
        if sha(path) != row['fields_sha256']:
            raise ValueError('Frozen diagnostic raw field changed: '+str(path))
        raw_integrity.append(dict(path=row['fields_path'], sha256=row['fields_sha256'], bytes=path.stat().st_size))
        if 'maximum_joint_ward_residual' in row or (row['energy_relative'] == .75 and row['eta_relative'] == .02):
            target = args.output/row['fields_path']
            target.parent.mkdir(exist_ok=True)
            shutil.copyfile(path, target)
            representative.append(dict(path=target.relative_to(args.output).as_posix(), sha256=sha(target)))
    probe_summary = {}
    for kind in ('radial', 'angular'):
        nonzero = [row for row in summary['records'] if row['energy_probe_weight'] > 0]
        maximum_flux = max(nonzero, key=lambda row:row['probes'][kind]['energy_flux_feedback_relative'])
        maximum_response = max(nonzero, key=lambda row:row['probes'][kind]['hT_to_hL_core_L2'])
        probe_summary[kind] = dict(maximum_energy_flux_change_relative=maximum_flux['probes'][kind]['energy_flux_feedback_relative'],
            largest_flux_change_id=maximum_flux['id'], largest_flux_change_metrics=maximum_flux['probes'][kind],
            maximum_hT_to_hL_core_L2=maximum_response['probes'][kind]['hT_to_hL_core_L2'],
            largest_response_id=maximum_response['id'],
            maximum_omitted_charge_residual_density_L2=max(row['probes'][kind]['omitted_charge_residual_density_L2'] for row in nonzero),
            maximum_corrected_charge_residual_density_L2=max(row['probes'][kind]['corrected_charge_residual_density_L2'] for row in nonzero))
    assessment = dict(schema='pysnspd.stage4.frozen_kinetic_review.v1',
        result='FROZEN_SPECTRUM_LINEAR_RESPONSE_COMPUTED_NOT_DYNAMIC_MODEL_ADMISSION',
        spectral_maps=72, mathematical_probe_responses=144, integrated_checks=12,
        runtime_seconds=summary['runtime_seconds'], verified_raw_files=len(raw_integrity),
        maximum_joint_ward_residual=max(row['maximum_joint_ward_residual'] for row in summary['integrated'].values()),
        maximum_uniform_energy_mode_free_residual=max(row['equilibrium_maximum_free_residual'] for row in summary['records']),
        exact_zero_force_increment=all(row['zero_increment_force_maximum'] == 0 for row in summary['records']),
        exact_zero_current_increment=all(row['zero_increment_current_maximum'] == 0 for row in summary['records']),
        probes=probe_summary,
        artificial_eta_leakage='Excluded from collisionless operator and reported separately; numerical eta is not a physical bath rate',
        interpretation='Frozen gap and potential: angular energy perturbations couple to a charge-mode response and modify energy flux. This does not alone invalidate the memory scalar closure with an evolving phase and electrochemical potential.',
        next_decisive_comparison='Compare spectral charge response and its moments with the electrochemical-potential and phase closure before choosing additional dynamic states.',
        mathematical_basis='Infinitesimal odd-compatible compact energy direction and radial/angular spatial directions; no finite Pauli-admitted distribution or photon deposit',
        energy_quadrature_scope='Finite 0..5 gap-reference trapezoid tests a shared algebraic identity, not convergence of a physical energy integral',
        not_admitted=['moving-gap energy dynamics', 'finite population trajectory', 'additional independent dynamic hT state', 'eta->0 closure', 'full thesis circuit coupling', 'photon', 'production promotion'],
        stage4_complete=False, production_changed=False)
    (args.output/'analysis.json').write_text(json.dumps(assessment, indent=2)+'\n')
    receipt = dict(remote_root='/home/jdiaz/scratch/pysnspd_stage4_frozen_kinetic_20260924',
        verified_raw_files=raw_integrity, raw_bytes=sum(row['bytes'] for row in raw_integrity),
        compact_fields=representative, analysis_script_sha256=sha(__file__),
        identity_sha256=sha(args.raw/'identity.json'), summary_sha256=sha(args.raw/'summary.json'))
    (args.output/'raw_integrity.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(assessment, indent=2))


if __name__ == '__main__':
    main()
