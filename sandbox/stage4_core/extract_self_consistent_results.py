"""Verify completed thermal-core checkpoints and export compact fields.

Read-only with respect to the completed simulation. No nonlinear solves or
time evolution. Full spectral arrays stay on Geminga, with verified hashes.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import multiprocessing as mp
import os
from pathlib import Path
import shutil
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
from parallel_runtime import initialize_affinity, limit_thread_environment, linux_resources, resource_budget
limit_thread_environment()
import numpy as np

PLAN = ROOT/'docs/implementation/stage4/spatial_energy_20260924/self_consistent_plan.json'


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for part in iter(lambda: stream.read(1024**2), b''):
            digest.update(part)
    return digest.hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf8'))


def write(path, value):
    with Path(path).open('x', encoding='utf8') as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write('\n')


def inspect_case(job):
    case, record, output, t, gap = job
    source = Path(record['fields_path'])
    if sha(source) != record['fields_sha256']:
        raise ValueError('Final checkpoint changed: '+str(source))
    with np.load(source) as data:
        fields = {name:data[name].copy() for name in data.files}
    if any(np.any(~np.isfinite(value)) for value in fields.values()):
        raise ValueError('Nonfinite checkpoint')
    d, u, f, g = (fields[name] for name in ('d', 'u', 'f', 'g'))
    count, nodes = f.shape
    if count != case['matsubara_count'] or nodes != case['nx']*case['ny']:
        raise ValueError('Unexpected final field dimensions')
    free = np.ones(nodes, bool)
    free[fields['boundary_nodes']] = False
    mass, edges, conductance = (fields[name] for name in ('area_weights', 'edges', 'conductance'))
    tail, head = edges.T
    coefficient = float(np.log(t)+np.sum(1/(np.arange(count)+.5)))
    fsum = np.sum(f, axis=0)
    target = 2*np.pi*t*fsum/coefficient
    predicted_next = target.copy()
    predicted_next[~free] = d[~free]
    reconstructed_force = 2*coefficient*d-4*np.pi*t*fsum
    reconstructed_energy = float(np.sum(mass*abs(d)**2)*np.log(t))
    reconstructed_current = np.zeros(len(edges))
    maximum_residual = 0.
    residual_by_n = []
    scale = mass*np.maximum(1., abs(d))
    for n in range(count):
        epsilon = 2*np.pi*t*(n+.5)
        mode_energy = np.sum(mass*(2*epsilon*abs(f[n])**2/(1+g[n])-2*np.real(np.conj(d)*f[n])))
        mode_energy += np.sum(conductance*(abs(f[n,tail]-f[n,head])**2+(g[n,tail]-g[n,head])**2))
        reconstructed_energy += 2*np.pi*t*(mode_energy+np.sum(mass*abs(d)**2)/epsilon)
        reconstructed_current += 4*np.pi*t*conductance*np.imag(np.conj(f[n,tail])*f[n,head])
        hz, hxy = mass*epsilon, mass*d
        np.add.at(hz, tail, conductance*g[n,head]); np.add.at(hz, head, conductance*g[n,tail])
        np.add.at(hxy, tail, conductance*f[n,head]); np.add.at(hxy, head, conductance*f[n,tail])
        value = float(np.max(abs((hz*u[n]-hxy)[free])/scale[free]))
        residual_by_n.append(value)
        maximum_residual = max(maximum_residual, value)
    core = free & (np.linalg.norm(fields['coordinates_bar'], axis=1) <= 4.)
    core_rms = float(np.sqrt(np.sum(mass[core]*abs((d-predicted_next)[core]/gap)**2)/np.sum(mass[core])))
    checks = dict(
        maximum_f_minus_g_u=float(np.max(abs(f-g*u))),
        maximum_normalization_error=float(np.max(abs(g*g+abs(f)**2-1))),
        maximum_positive_chart_error=float(np.max(abs(g-1/np.hypot(1., abs(u))))),
        maximum_next_d_error=float(np.max(abs(predicted_next-fields['next_d']))),
        maximum_force_error=float(np.max(abs(reconstructed_force-fields['force']))),
        maximum_current_error=float(np.max(abs(reconstructed_current-fields['current_bar']))),
        energy_error=abs(reconstructed_energy-record['metrics']['energy']),
        maximum_recomputed_spectral_residual=maximum_residual,
        core_mass_rms_relative=core_rms,
        core_mass_rms_record_error=abs(core_rms-record['metrics']['core_mass_rms_relative']),
        fixed_gap_block_change=float(np.max(abs(fields['next_d'][~free]-d[~free]))),
        coherent_pair='All saved d,u,f,g,force,current checked together; next_d remains a separate proposal')
    for name in ('maximum_f_minus_g_u', 'maximum_normalization_error', 'maximum_positive_chart_error',
                 'maximum_next_d_error', 'maximum_force_error', 'maximum_current_error',
                 'core_mass_rms_record_error', 'fixed_gap_block_change'):
        if checks[name] > 2e-11:
            raise ValueError('Incoherent fields: '+name)
    if checks['energy_error'] > 2e-10 or maximum_residual > 1e-6:
        raise ValueError('Energy or spectral stationarity failed')
    compact_names = ('d', 'next_d', 'force', 'current_bar', 'coordinates_bar', 'area_weights', 'edges', 'conductance', 'boundary_nodes')
    compact = {name:fields[name] for name in compact_names}
    compact.update(f_lowest=f[0], g_lowest=g[0], f_sum=fsum, spectral_residuals=np.asarray(residual_by_n))
    path = Path(output)/(case['id']+'.npz')
    np.savez_compressed(path, **compact)
    return dict(case_id=case['id'], full_checkpoint_path=str(source), full_checkpoint_sha256=sha(source),
        compact_file=path.name, compact_sha256=sha(path), compact_bytes=path.stat().st_size, checks=checks)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    if args.output.exists():
        raise FileExistsError('Choose a fresh extraction directory')
    identity, summary, plan = read(args.raw/'identity.json'), read(args.raw/'summary.json'), read(PLAN)
    if identity['plan_sha256'] != sha(PLAN):
        raise ValueError('Plan hash changed')
    sources = {}
    for name, expected in identity['sources'].items():
        actual = sha(ROOT/name)
        if actual != expected:
            raise ValueError('Executed source changed: '+name)
        sources[name] = actual
    if summary['status'] != 'FINITE_SUM_CORE_STATIONARY' or not summary['all_core_criteria_met']:
        raise ValueError('Run did not meet its declared core criterion')
    budget = resource_budget(linux_resources(), .9, .9, max_workers=4)
    args.output.mkdir(parents=True)
    (args.output/'histories').mkdir()
    for name in ('identity.json', 'summary.json', 'progress.jsonl'):
        shutil.copyfile(args.raw/name, args.output/name)
    shutil.copyfile(PLAN, args.output/'executed_plan.json')
    history = {}
    verified = []
    directories = [Path(identity['resumed_from']), args.raw]
    pilot_identity = read(directories[0]/'identity.json')
    for key in ('plan_sha256', 'sources', 'radial_reference_identity_sha256'):
        if pilot_identity[key] != identity[key]:
            raise ValueError('Pilot identity differs: '+key)
    for directory in directories:
        for path in sorted((directory/'checkpoints').glob('*.json')):
            record = read(path)
            key = (record['case_id'], record['completed_sweeps'])
            if key in history:
                raise ValueError('Duplicate completed sweep')
            history[key] = record
            shutil.copyfile(path, args.output/'histories'/path.name)
            if record['full_checkpoint']:
                actual = sha(record['fields_path'])
                if actual != record['fields_sha256']:
                    raise ValueError('Saved checkpoint hash differs')
                verified.append(dict(path=record['fields_path'], sha256=actual, bytes=Path(record['fields_path']).stat().st_size))
    initial = {case['id']:max(sweep for cid,sweep in history if cid==case['id'] and sweep<=2) for case in plan['cases']}
    final_records = {}
    for case in plan['cases']:
        cid = case['id']; last = summary['cases'][cid]
        if [sweep for cc,sweep in sorted(history) if cc==cid] != list(range(1,last['completed_sweeps']+1)):
            raise ValueError('Missing sweep history')
        if history[(cid,last['completed_sweeps'])] != last:
            raise ValueError('Final summary is not the final checkpoint')
        final_records[cid] = last
    jobs = [(case, final_records[case['id']], str(args.output), plan['T_K']/plan['Tc_K'], identity['gap_reference']) for case in plan['cases']]
    results = {}
    original = os.sched_getaffinity(0)
    try:
        os.sched_setaffinity(0, {budget['coordinator_cpu']})
        context = mp.get_context('spawn'); counter = context.Value('i', 0)
        with ProcessPoolExecutor(max_workers=budget['workers'], mp_context=context, initializer=initialize_affinity,
                initargs=(budget['worker_affinity_cpus'],counter)) as pool:
            futures = [pool.submit(inspect_case, job) for job in jobs]
            for future in as_completed(futures):
                result = future.result(); results[result['case_id']] = result
                done = len(results); elapsed = time.monotonic()-started
                print(json.dumps(dict(event='EXTRACTED', completed=done, total=len(jobs), elapsed_seconds=elapsed,
                    eta_seconds=elapsed/done*(len(jobs)-done), case=result['case_id'])), flush=True)
    finally:
        os.sched_setaffinity(0, original)
    receipt = dict(schema='pysnspd.stage4.self_consistent_extraction.v1', raw_directory=str(args.raw.resolve()),
        runtime_seconds=time.monotonic()-started, extraction_script_sha256=sha(__file__),
        plan_sha256=sha(PLAN), identity_sha256=sha(args.raw/'identity.json'), summary_sha256=sha(args.raw/'summary.json'),
        sources_matching_identity=sources, checkpoint_files_verified=verified,
        completed_sweep_records=len(history), final_cases=results, process_budget=budget,
        physical_time_steps=0, new_spectral_solves=0, production_changed=False)
    write(args.output/'extraction_receipt.json', receipt)
    print(json.dumps(dict(status='EXTRACTED_NO_SOLVES', runtime_seconds=receipt['runtime_seconds'],
        verified_checkpoints=len(verified), verified_sweeps=len(history))))


if __name__ == '__main__':
    main()
