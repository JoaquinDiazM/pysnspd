"""Frozen-spectrum charge-mode response with exact force/current increments.

Uses already computed retarded maps, not additional spectral solves. Prescribed
energy-mode directions are abstract linear probes. No finite population, moving
gap, electrostatic dynamics, detector time or photon is admitted by this test.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import multiprocessing as mp
import os
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
from parallel_runtime import initialize_affinity, limit_thread_environment, linux_resources, resource_budget
limit_thread_environment()
sys.path.insert(0, str(ROOT))
import numpy as np
from pysnspd.experimental import thermal_spatial_usadel as thermal
from pysnspd.experimental import frozen_kinetic_usadel as kinetic

DEFAULT_PLAN = ROOT/'docs/implementation/stage4/self_consistent_review_20260924/next_kinetic_plan.json'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    with Path(path).open('x', encoding='utf8') as stream:
        json.dump(value, stream, indent=2, allow_nan=False); stream.write('\n')


def worker(job):
    start = time.monotonic()
    plan, case, row = job['plan'], job['case'], job['record']
    with np.load(ROOT/case['fields_path']) as data:
        d, xy, mass = data['d'], data['coordinates_bar'], data['area_weights']
        graph = thermal.ThermalGraph(mass, data['edges'], data['conductance'], xy, data['boundary_nodes'])
    with np.load(Path(job['spectra_root'])/row['fields_path']) as raw:
        if not np.array_equal(d, raw['d']):
            raise ValueError('Gap changed between thermal core and retarded spectrum')
        operator = kinetic.assemble(graph, d, raw['g'], raw['f'], raw['f_tilde'])
    radius = np.linalg.norm(xy, axis=1)
    core = radius <= plan['probe_radius_ell0']
    free = np.ones(len(d), bool); free[graph.boundary_nodes] = False
    norm = lambda value: float(np.sqrt(np.dot(mass[core], abs(value[core])**2)))
    midpoint = (xy[graph.edges[:, 0]]+xy[graph.edges[:, 1]])/2
    active_edges = (graph.conductance > 0) & (np.linalg.norm(midpoint, axis=1) <= plan['probe_radius_ell0'])
    flux_norm = lambda value: float(np.sqrt(np.sum(abs(value[active_edges])**2/graph.conductance[active_edges])))
    equilibrium = operator.evaluate(np.column_stack((np.ones(len(d)), np.zeros(len(d)))), eta=row['z_real'])
    zero = operator.evaluate(np.zeros((len(d), 2)), eta=row['z_real'])
    records, fields = {}, {}
    xenergy = row['energy_relative']
    # Odd-compatible compact mathematical basis; no photon preparation or
    # physical population is inferred from its arbitrary normalization.
    energy_basis = xenergy*max(0., 1-xenergy*xenergy/25)**3
    for kind in plan['probes']:
        bump = np.maximum(0., 1-radius**2/plan['probe_radius_ell0']**2)**3
        hL = energy_basis*(bump if kind == 'radial' else xy[:, 0]/plan['probe_radius_ell0']*bump)
        original = np.column_stack((hL, np.zeros(len(d))))
        corrected = operator.charge_response(hL, graph.boundary_nodes)
        before = operator.evaluate(original, eta=row['z_real'])
        after = operator.evaluate(corrected, eta=row['z_real'])
        omitted = norm(before.residual[:, 1]/mass)
        remaining = norm(after.residual[:, 1]/mass)
        energy_before = norm(before.residual[:, 0]/mass)
        feedback = norm((after.residual[:, 0]-before.residual[:, 0])/mass)
        energy_flux_before = flux_norm(before.edge_flux[:, 0])
        energy_flux_feedback = flux_norm(after.edge_flux[:, 0]-before.edge_flux[:, 0])
        records[kind] = dict(omitted_charge_residual_density_L2=omitted,
            corrected_charge_residual_density_L2=remaining,
            maximum_free_charge_residual=float(np.max(abs(after.residual[free, 1]))),
            maximum_contact_charge_reaction=float(np.max(abs(after.residual[~free, 1]))),
            hT_maximum_absolute=float(np.max(abs(corrected[:, 1]))),
            hT_to_hL_core_L2=(norm(corrected[:, 1])/norm(hL) if norm(hL)>0 else None),
            energy_residual_feedback_relative=(feedback/energy_before if energy_before > 1e-14 else None),
            energy_flux_feedback_relative=(energy_flux_feedback/energy_flux_before if energy_flux_before > 1e-14 else None),
            energy_flux_feedback_norm=energy_flux_feedback,
            energy_flux_norm_definition='sqrt(sum(J_L^2/c)) on active edges with midpoint r<=4ell0',
            energy_residual_density_L2_after=norm(after.residual[:, 0]/mass),
            maximum_ward_residual=float(np.max(abs(after.ward_residual))),
            artificial_eta_energy_leakage_density_L2=norm(after.artificial_eta_leakage[:, 0]/mass),
            artificial_eta_charge_leakage_density_L2=norm(after.artificial_eta_leakage[:, 1]/mass))
        fields.update({kind+'_hL':hL, kind+'_hT':corrected[:, 1],
            kind+'_delta_gap_force':after.gap_force_increment,
            kind+'_delta_current':after.charge_current_increment,
            kind+'_charge_residual':after.residual[:, 1],
            kind+'_energy_residual':after.residual[:, 0],
            kind+'_energy_flux_before':before.edge_flux[:, 0],
            kind+'_energy_flux_after':after.edge_flux[:, 0],
            kind+'_omitted_charge_residual':before.residual[:, 1]})
    record = dict(id=row['id'], case_id=case['id'], energy_relative=row['energy_relative'],
        eta_relative=row['eta_relative'], energy_probe_weight=energy_basis, worker_seconds=time.monotonic()-start, worker_pid=os.getpid(),
        equilibrium_maximum_free_residual=float(np.max(abs(equilibrium.residual[free]))),
        equilibrium_artificial_eta_leakage_maximum=float(np.max(abs(equilibrium.artificial_eta_leakage))),
        zero_increment_force_maximum=float(np.max(abs(zero.gap_force_increment))),
        zero_increment_current_maximum=float(np.max(abs(zero.charge_current_increment))), probes=records)
    return record, fields


def integrated_checks(plan, records, stored, output):
    result = {}
    gap = plan['gap_reference_kBTc']
    for case in plan['cases']:
        with np.load(ROOT/case['fields_path']) as data:
            d, xy, edges = data['d'], data['coordinates_bar'], data['edges']
        for eta in sorted({r['eta_relative'] for r in records}):
            selected = sorted([r for r in records if r['case_id'] == case['id'] and r['eta_relative'] == eta], key=lambda r:r['energy_relative'])
            if len(selected) != len(plan['energies_relative']):
                continue
            energies = gap*np.array([r['energy_relative'] for r in selected])
            weights = np.zeros(len(energies))
            weights[:-1] += np.diff(energies)/2; weights[1:] += np.diff(energies)/2
            for kind in plan['probes']:
                force = sum(w*stored[r['id']][kind+'_delta_gap_force'] for w,r in zip(weights, selected))
                current = sum(w*stored[r['id']][kind+'_delta_current'] for w,r in zip(weights, selected))
                charge = sum(w*stored[r['id']][kind+'_charge_residual'] for w,r in zip(weights, selected))
                divergence = np.zeros(len(d))
                np.add.at(divergence, edges[:, 0], current); np.add.at(divergence, edges[:, 1], -current)
                ward = divergence+np.imag(np.conj(d)*force)-2*charge
                key = f"{case['id']}_eta{eta:.4f}_{kind}".replace('.', 'p')
                path = output/'integrated'/(key+'.npz')
                np.savez_compressed(path, coordinates_bar=xy, d=d, edges=edges,
                    delta_gap_force=force, delta_current=current, charge_residual=charge, ward_residual=ward,
                    energy_nodes_kBTc=energies, energy_weights_kBTc=weights)
                result[key] = dict(maximum_joint_ward_residual=float(np.max(abs(ward))),
                    maximum_gap_force_increment=float(np.max(abs(force))),
                    maximum_charge_current_increment=float(np.max(abs(current))),
                    fields_path=path.relative_to(output).as_posix(), fields_sha256=sha(path),
                    scope='Same finite-window trapezoid quadrature in force,current,residual; checks algebra, not physical energy-grid accuracy')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, default=DEFAULT_PLAN)
    parser.add_argument('--spectra-root', type=Path, required=True)
    parser.add_argument('--output-root', type=Path, required=True)
    parser.add_argument('--pilot', action='store_true')
    parser.add_argument('--max-workers', type=int)
    parser.add_argument('--execute', action='store_true')
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    if plan['schema'] != 'pysnspd.stage4.frozen_kinetic_response.v1' or set(plan['probes']) != {'radial', 'angular'}:
        raise ValueError('Unknown frozen-spectrum diagnostic plan')
    inputs = {}
    for case in plan['cases']:
        actual = sha(ROOT/case['fields_path'])
        if actual != case['fields_sha256']:
            raise ValueError('Completed thermal field changed')
        inputs[case['fields_path']] = actual
    summary_path = args.spectra_root/'summary.json'
    if sha(summary_path) != plan['spectral_summary_sha256']:
        raise ValueError('Frozen retarded campaign summary changed')
    summary = json.loads(summary_path.read_text())
    cases = {case['id']:case for case in plan['cases']}
    jobs, spectral_hashes = [], {}
    for record in summary['records']:
        if args.pilot and record['id'] not in plan['pilot_ids']:
            continue
        path = args.spectra_root/record['fields_path']
        if sha(path) != record['fields_sha256']:
            raise ValueError('Frozen retarded spectrum changed: '+record['id'])
        spectral_hashes[record['fields_path']] = record['fields_sha256']
        jobs.append(dict(plan=plan, case=cases[record['case_id']], record=record, spectra_root=str(args.spectra_root.resolve())))
    if not jobs:
        raise ValueError('No registered spectral maps selected')
    resources = linux_resources()
    budget = resource_budget(resources, .9, .9, min(len(jobs), args.max_workers or len(jobs)))
    paths = [Path(__file__), ROOT/'pysnspd/experimental/frozen_kinetic_usadel.py',
             ROOT/'pysnspd/experimental/thermal_spatial_usadel.py', HERE/'parallel_runtime.py']
    identity = dict(plan_sha256=sha(args.plan), sources={p.relative_to(ROOT).as_posix():sha(p) for p in paths},
        inputs=inputs, spectral_summary_sha256=sha(summary_path), spectral_fields=spectral_hashes,
        resources=resources, budget=budget, jobs=len(jobs), pilot=args.pilot, physical_time_steps=0,
        finite_population_state_admitted=False, production_changed=False)
    if not args.execute:
        print(json.dumps(dict(status='DRY_RUN', **identity), indent=2)); return
    if args.output_root.exists():
        raise FileExistsError('Choose a new output directory; no overwrite')
    args.output_root.mkdir(parents=True)
    for child in ('responses', 'integrated'):
        (args.output_root/child).mkdir()
    write_json(args.output_root/'identity.json', identity)
    started, records, stored = time.monotonic(), [], {}
    affinity = os.sched_getaffinity(0)
    log = (args.output_root/'progress.jsonl').open('x', encoding='utf8', buffering=1)
    def event(name, **values):
        line = json.dumps(dict(event=name, elapsed_seconds=time.monotonic()-started, **values))
        print(line, flush=True); log.write(line+'\n')
    try:
        event('START', jobs=len(jobs), budget=budget)
        os.sched_setaffinity(0, {budget['coordinator_cpu']})
        context = mp.get_context('spawn'); counter = context.Value('i', 0)
        with ProcessPoolExecutor(max_workers=budget['workers'], mp_context=context,
            initializer=initialize_affinity, initargs=(budget['worker_affinity_cpus'], counter)) as pool:
            futures = [pool.submit(worker, job) for job in jobs]
            try:
                for future in as_completed(futures):
                    record, fields = future.result()
                    path = args.output_root/'responses'/(record['id']+'.npz')
                    np.savez_compressed(path, **fields)
                    record.update(fields_path=path.relative_to(args.output_root).as_posix(), fields_sha256=sha(path))
                    records.append(record); stored[record['id']] = fields
                    write_json(path.with_suffix('.json'), record)
                    fraction = len(records)/len(jobs); elapsed = time.monotonic()-started
                    eta = elapsed*(1-fraction)/fraction
                    bar = '#'*int(24*fraction)+'-'*(24-int(24*fraction))
                    print(f'[{bar}] respuestas cineticas {len(records)}/{len(jobs)} | {elapsed/60:.1f} min | ETA {eta/60:.1f} min', flush=True)
                    event('PROGRESS', complete=len(records), total=len(jobs), eta_seconds=eta)
            except BaseException:
                for future in futures:
                    future.cancel()
                raise
        result = dict(status='FROZEN_KINETIC_PILOT_COMPLETE' if args.pilot else 'FROZEN_KINETIC_RESPONSE_COMPLETE',
            runtime_seconds=time.monotonic()-started, records=records,
            integrated={} if args.pilot else integrated_checks(plan, records, stored, args.output_root),
            stage4_complete=False, finite_population_state_admitted=False, full_energy_dynamics_admitted=False,
            eta_as_physical_relaxation=False, production_changed=False, physical_time_steps=0)
        result['runtime_seconds'] = time.monotonic()-started
        write_json(args.output_root/'summary.json', result)
        event('COMPLETE', status=result['status'], seconds=result['runtime_seconds'])
    except BaseException as exc:
        event('FAILED', error=type(exc).__name__, reason=str(exc))
        write_json(args.output_root/'failure.json', dict(error=type(exc).__name__, reason=str(exc)))
        raise
    finally:
        os.sched_setaffinity(0, affinity); log.close()


if __name__ == '__main__':
    main()
