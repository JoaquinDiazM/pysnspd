"""Compact, hash-verified evidence from the registered retarded graph oracle."""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import numpy as np

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
    plan_path = ROOT/'docs/implementation/stage4/self_consistent_review_20260924/next_retarded_plan.json'
    plan = json.loads(plan_path.read_text())
    if sha(plan_path) != identity['plan_sha256']:
        raise ValueError('Plan mismatch')
    for name, digest in {**identity['sources'], **identity['inputs']}.items():
        if sha(ROOT/name) != digest:
            raise ValueError('Executed source/input changed: '+name)
    if summary['status'] != 'RETARDED_ORACLE_COMPLETE' or len(summary['records']) != 72:
        raise ValueError('This assessment requires the complete registered oracle')
    if args.output.exists():
        raise FileExistsError('Fresh compact-evidence directory required')
    args.output.mkdir(parents=True)
    for name in ('identity.json', 'summary.json', 'progress.jsonl'):
        shutil.copyfile(args.raw/name, args.output/name)
    arrays, certificates, records = {}, [], {}
    for row in summary['records']:
        path = args.raw/row['fields_path']
        if sha(path) != row['fields_sha256']:
            raise ValueError('Raw spectral fields hash mismatch: '+str(path))
        certificates.append(dict(path=row['fields_path'], sha256=row['fields_sha256'], bytes=path.stat().st_size))
        key = (row['case_id'], row['eta_relative'], row['energy_relative'])
        with np.load(path) as raw:
            arrays[key] = {name:raw[name].copy() for name in ('DOS', 'd', 'coordinates_bar', 'area_weights')}
            # Retain representative complex fields as well as all DOS maps.
            if row['eta_relative'] == .02 and row['energy_relative'] in (.5, 1.):
                np.savez_compressed(args.output/(row['id']+'.npz'), **{name:raw[name] for name in raw.files})
        records[key] = row
    for case in plan['cases']:
        selected = sorted((key for key in arrays if key[0] == case['id']), key=lambda key:(key[1], key[2]))
        first = arrays[selected[0]]
        np.savez_compressed(args.output/(case['id']+'_DOS_maps.npz'),
            eta_relative=np.array([key[1] for key in selected]),
            energy_relative=np.array([key[2] for key in selected]),
            DOS=np.array([arrays[key]['DOS'] for key in selected]),
            d=first['d'], coordinates_bar=first['coordinates_bar'], area_weights=first['area_weights'])
    mesh = []
    for eta in plan['etas_relative']:
        for energy in plan['energies_relative']:
            low = arrays['radial_65_N256', eta, energy]
            high = arrays['radial_129_N256', eta, energy]
            restricted = high['DOS'].reshape(129, 129)[::2, ::2].ravel()
            xy = low['coordinates_bar']; weights = low['area_weights']
            core = np.linalg.norm(xy, axis=1) <= plan['core_radius_ell0']
            norm = lambda value: float(np.sqrt(np.dot(weights[core], value[core]**2)))
            mesh.append(dict(eta_relative=eta, energy_relative=energy,
                core_DOS_L2_relative=norm(low['DOS']-restricted)/norm(restricted),
                core_DOS_maximum_absolute=float(np.max(abs(low['DOS'][core]-restricted[core]))),
                definition='65x65 DOS minus coincident 129x129 nodes, weighted in r<=4ell0; includes the independently relaxed gap difference'))
    contour = []
    for case in plan['cases']:
        for energy in plan['energies_relative']:
            low, high = arrays[case['id'], .02, energy], arrays[case['id'], .04, energy]
            core = np.linalg.norm(low['coordinates_bar'], axis=1) <= plan['core_radius_ell0']
            weight = low['area_weights'][core]
            difference = high['DOS'][core]-low['DOS'][core]
            contour.append(dict(case_id=case['id'], energy_relative=energy,
                core_DOS_L2_relative=float(np.sqrt(np.dot(weight, difference**2)/np.dot(weight, low['DOS'][core]**2))),
                core_mean_DOS_eta002=records[case['id'], .02, energy]['core_area_mean_DOS'],
                core_mean_DOS_eta004=records[case['id'], .04, energy]['core_area_mean_DOS'],
                definition='Numerical contour difference .04 vs .02, not a measured broadening or a rigorous eta->0 error bound'))
    all_steps = [step for row in summary['records'] for step in row['continuation']]
    assessment = dict(schema='pysnspd.stage4.retarded_oracle_review.v1',
        result='REGISTERED_CAUSAL_SPECTRAL_ORACLE_PASSED_NOT_KINETIC_ADMISSION',
        jobs=len(summary['records']), continuation_steps=len(all_steps), runtime_seconds=summary['runtime_seconds'],
        source_hashes_verified=True, raw_fields_hashes_verified=len(certificates),
        maximum_root_residual=max(step['root_residual'] for step in all_steps),
        maximum_normalization_residual=max(step['normalization_residual'] for step in all_steps),
        minimum_DOS_along_continuation=min(step['minimum_DOS'] for step in all_steps),
        maximum_radial_rms_residual=max(step['exterior']['maximum_rms_residual'] for step in all_steps),
        maximum_matsubara_action_bridge_error=max(row['bridge']['action_absolute_difference'] for row in summary['records']),
        maximum_matsubara_current_bridge_error=max(row['bridge']['current_maximum_difference'] for row in summary['records']),
        mesh_comparison=mesh, contour_comparison=contour,
        maximum_mesh_core_DOS_L2_relative=max(row['core_DOS_L2_relative'] for row in mesh),
        maximum_contour_core_DOS_L2_relative=max(row['core_DOS_L2_relative'] for row in contour),
        scope='Electronic local DOS Re(g); distinct from the phonon DOS input. Fixed thermal gaps and matched prescribed-tanh radial exterior.',
        accepted_for='Frozen-spectrum diagnostics and construction of the next nonequilibrium bridge',
        not_admitted=['eta->0 convergence', 'local quasiparticle counts/occupation closure', 'nonequilibrium energy force/current', 'full kinetic spatial dynamics', 'photon', 'production promotion'],
        stage4_complete=False, production_changed=False)
    (args.output/'analysis.json').write_text(json.dumps(assessment, indent=2)+'\n')
    receipt = dict(remote_root='/home/jdiaz/scratch/pysnspd_stage4_retarded_oracle_20260924',
        verified_downloaded_files=certificates, raw_bytes=sum(x['bytes'] for x in certificates),
        analysis_script_sha256=sha(__file__), identity_sha256=sha(args.raw/'identity.json'),
        summary_sha256=sha(args.raw/'summary.json'), compact_outputs={p.name:sha(p) for p in args.output.glob('*.npz')})
    (args.output/'raw_integrity.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({k:v for k,v in assessment.items() if k not in ('mesh_comparison', 'contour_comparison')}, indent=2))


if __name__ == '__main__':
    main()
