"""Read completed archives and frozen criteria; no RHS or integration is called."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'sandbox/stage2_cells/recovery_20260921'),str(ROOT/'sandbox/stage2_cells')]
import assess_limited_time as contract
original=contract.original
FIELDS=('amplitudes','excitation_energy','electron_energy','phonon_energy',
        'escape','input','electron_to_phonon','condensate_heat','transport')


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def rel(path):return Path(path).resolve().relative_to(ROOT).as_posix()


def comparison(actual,reference):
    contract.compare(actual,reference)
    rows=[]
    for j,matches in original.shared_indices([actual],reference):
        rows.append(dict(time=float(reference['times'][j]),
            errors=original.checkpoint_errors(actual,matches[0],reference,j,FIELDS)))
    names=rows[0]['errors'];maxima={}
    for name in names:
        measured=[row['errors'][name]|{'time':row['time']} for row in rows]
        defined=all(value['comparison_defined'] for value in measured)
        maxima[name]=max(measured,key=lambda x:x['assessment_error']) if defined else {'comparison_defined':False}
    defined=all(v.get('comparison_defined') for v in maxima.values())
    worst=max(maxima,key=lambda name:maxima[name]['assessment_error']) if defined else None
    return dict(actual_path=rel(actual['path']),reference_path=rel(reference['path']),
        common_checkpoint_count=len(rows),common_times=[row['time'] for row in rows],
        maximum_error=maxima[worst]['assessment_error'] if worst else None,
        worst_observable=worst,per_observable_maxima=maxima,checkpoints=rows,
        status='DIAGNOSTIC_ONLY_NOT_ADMISSION')


def summary(run,criteria):
    record=run['record'];limit=criteria['coupled_trajectories']
    residuals=np.array([s['conserved']-record['initial']['conserved'] for s in record['snapshots']])
    return dict(path=rel(run['path']),json_sha256=sha(run['path']),
        trajectory_sha256=record['trajectory_sha256'],parameters=record['parameters'],
        runtime_seconds=record['runtime_seconds'],rhs_calls=record['rhs_calls'],
        energy_ledger_scaled_max=record['energy_ledger_scaled_max'],
        ledger_limit=limit['energy_ledger_scaled_max'],
        ledger_to_limit_ratio=record['energy_ledger_scaled_max']/limit['energy_ledger_scaled_max'],
        initial_conserved=record['initial']['conserved'],final_conserved=record['final']['conserved'],
        final_ledger_difference=float(residuals[-1]),
        maximum_ledger_time=float(run['times'][np.argmax(abs(residuals))]),
        instantaneous_residual_max=record['instantaneous_residual_max'],
        minimum_electron=record['minimum_electron'],maximum_electron=record['maximum_electron'],
        minimum_phonon=record['minimum_phonon'],invariant_checks=original.invariant_checks(record,criteria),
        integration=record.get('time_integration'),ledger_series=residuals.tolist(),times=run['times'].tolist())


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'docs/implementation/stage2/review_20260922/diagnosed_ssp.json')
    args=parser.parse_args()
    if args.output.exists():raise SystemExit('Existing diagnosis is preserved; choose a new output.')
    started=time.perf_counter();base=ROOT/'docs/implementation/stage2'
    criteria_path=base/'acceptance_criteria.json';criteria=json.loads(criteria_path.read_text())
    ssp=original.load_run(base/'review_20260922/raw/ssp/one_ssp_40.json')
    rk={n:original.load_run(base/f'review_20260922/raw/rk4/one_time_{n}.json') for n in (40,80,160,320)}
    short={n:original.load_run(base/f'recovery_20260921/one_ssp_short_{n}.json') for n in (10,20,40)}
    for run in [ssp,*rk.values(),*short.values()]:contract.verify_contract(run)
    ledger_short=[short[n]['record']['energy_ledger_scaled_max'] for n in (10,20,40)]
    short_ratios=[a/b for a,b in zip(ledger_short[:-1],ledger_short[1:])]
    short_orders=[float(np.log2(ratio)) for ratio in short_ratios]
    full_error=ssp['record']['energy_ledger_scaled_max']
    projections=[dict(steps=n,macro_dt=2/n,
        projected_ledger_from_third_order=full_error*(40/n)**3,
        projected_ledger_to_limit=full_error*(40/n)**3/criteria['coupled_trajectories']['energy_ledger_scaled_max'],
        runtime_linear_estimate_seconds=ssp['record']['runtime_seconds']*n/40,
        measured=False) for n in (80,160,320,640,1280)]
    comparisons={f'ssp40_vs_rk{n}':comparison(ssp,rk[n]) for n in (40,320)}
    adjacent={f'rk{a}_vs_rk{b}':comparison(rk[a],rk[b]) for a,b in ((40,80),(80,160),(160,320))}
    estimates={}
    for name in adjacent['rk160_vs_rk320']['per_observable_maxima']:
        earlier=adjacent['rk80_vs_rk160']['per_observable_maxima'][name]
        last=adjacent['rk160_vs_rk320']['per_observable_maxima'][name]
        a,b=earlier.get('absolute'),last.get('absolute')
        ratio=a/b if b and a is not None else None
        relative=last.get('relative')
        estimates[name]=dict(last_difference_relative=relative,last_difference_absolute=b,
            last_difference_time=last.get('time'),consecutive_absolute_difference_ratio=ratio,
            inferred_order=float(np.log2(ratio)) if ratio is not None and ratio>0 else None,
            estimated_rk320_remainder_relative=(relative/(ratio-1)
                if relative is not None and ratio is not None and ratio>1 else None),
            limitation='Heuristic only: maxima may occur at different checkpoints; asymptotic order and error sign are not certified.')
    final,reference=ssp['record']['final'],rk[320]['record']['final']
    total=lambda snap,key:float(np.asarray(snap[key]).sum())
    decomp=dict(vacuum_energy_difference=(total(final,'electron_energy')-total(final,'excitation_energy')
        -total(reference,'electron_energy')+total(reference,'excitation_energy')),
        qp_energy_difference=total(final,'excitation_energy')-total(reference,'excitation_energy'),
        phonon_energy_difference=total(final,'phonon_energy')-total(reference,'phonon_energy'),
        escaped_energy_difference=total(final,'escape')-total(reference,'escape'),
        input_energy_difference=total(final,'input')-total(reference,'input'),
        conserved_total_difference=final['conserved']-reference['conserved'])
    stats=ssp['record']['time_integration']['stats']
    inputs=[ssp['path'],*[r['path'] for r in rk.values()],*[r['path'] for r in short.values()]]
    input_hashes={rel(p):sha(p) for path in inputs for p in
        (path,path.with_suffix('.npz'),path.with_name(path.stem+'_initial.npz'))}
    chosen=[160,320,640,1280]
    result=dict(schema='pysnspd.stage2.ssp-failure-diagnosis.v1',status='NOT_CLOSED_NO_NEW_TRAJECTORY',
        scope='Postprocessing of complete recorded archives; no RHS evaluation or time integration.',
        criteria_sha256=sha(criteria_path),source_sha256=sha(__file__),
        assessor_sha256=sha(contract.__file__),input_hashes=input_hashes,
        ssp40=summary(ssp,criteria),rk4=[summary(rk[n],criteria) for n in rk],
        short_ssp=[summary(short[n],criteria) for n in short],
        short_ledger_halving_ratios=short_ratios,short_ledger_observed_orders=short_orders,
        diagnosis=dict(failed_gates=[k for k,v in original.invariant_checks(ssp['record'],criteria).items() if not v],
            limiter_active=stats['limited_event_evaluations']>0,
            limited_events=stats['limited_event_evaluations'],events_evaluated=stats['event_evaluations'],
            energy_weighted_flux_defect=stats['integrated_energy_weighted_flux_defect'],
            interpretation='The whole-duration SSP40 trajectory fails accumulated nonlinear energy accounting at the chosen time resolution. The continuous RHS balance and populations pass. With zero limited events, this failure cannot be attributed to event limiting; the observed short-time energy error decreases at approximately third order. This is numerical evidence, not evidence of a failed physical model.'),
        final_energy_difference_decomposition=decomp,comparisons=comparisons,
        adjacent_rk_comparisons=adjacent,rk320_error_floor_diagnostics=estimates,
        ledger_and_cost_projections=projections,
        proposed_next_diagnostic=dict(status='PROPOSAL_NOT_REGISTERED_NOT_LAUNCHED',
            one_cell_only=True,duration=2.,methods_unchanged=True,criteria_unchanged=True,
            tested_resolutions=[160,320,640],separate_step_halved_reference=1280,
            justification='80 steps is predicted to fail the unchanged energy gate. Use a separate SSP1280 reference to avoid mistaking the reused RK320 temporal-error floor for failed SSP refinement. A separately step-halved reference is explicitly permitted by the frozen criteria. Retain RK320 and short DOP853 as independent cross-checks, not an assumed exact long-time answer.',
            ordering='Run SSP160 first; stop immediately if invariant or population gates fail. Then SSP320,640,1280, followed by the unchanged three-level temporal assessment and a separately labelled RK320 cross-check. Do not launch two-cell/mesh phases until this diagnostic passes.',
            estimated_total_seconds=sum(ssp['record']['runtime_seconds']*n/40 for n in chosen),
            estimate_basis='Linear scaling from the measured 25.8-second SSP40 run on the same candidate mesh; a resource estimate, not a runtime guarantee.',
            required_before_any_restart='Freeze a new plan and output directory; preserve both failed batch histories. No solver, source, or historical plan is edited by this diagnosis.',
            scope_limit='Passing would validate this one-cell numerical case only; limiter-active fine-grid, two-cell, mesh and final model gates remain separate.'),
        runtime_seconds=time.perf_counter()-started)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x',encoding='utf-8') as stream:stream.write(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(status=result['status'],failed_gates=result['diagnosis']['failed_gates'],
        short_orders=short_orders,ssp40_vs_rk320_max=comparisons['ssp40_vs_rk320']['maximum_error'],
        worst_observable=comparisons['ssp40_vs_rk320']['worst_observable'],
        next_estimated_minutes=result['proposed_next_diagnostic']['estimated_total_seconds']/60,
        output=str(args.output)),indent=2))


if __name__=='__main__':main()
