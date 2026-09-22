"""Offline correctness guards for local RK4 bisection; no detector calculation."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from pysnspd.experimental.cell_validation import rk4_trajectory
from positive_rk4 import integrate, SupportViolation, IntegrationLimitError
import numpy as np


def validate_nonnegative(values):
    failed = np.flatnonzero(values < 0)
    if len(failed):
        i = int(failed[0])
        raise SupportViolation('negative toy population',
            details={'sector':'toy', 'index':i, 'value':float(values[i]),
                     'magnitude':float(-values[i]), 'lower_bound':0.})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit('Existing test evidence is preserved; choose a new output path.')
    started = time.perf_counter()
    checks, details = {}, {}

    # This includes nonautonomous forcing and a duration whose dt is not an
    # exact binary fraction, to check both arithmetic and the macro time grid.
    def smooth_rhs(t, y):
        return np.array([-.3*y[0]+.02*np.cos(t), .1*y[0]-.2*y[1], .07*y[1]])
    initial = np.array([.7,.2,0.])
    untouched = initial.copy()
    times0, states0 = rk4_trajectory(smooth_rhs, initial, .7, 13)
    times1, states1, stats = integrate(smooth_rhs, initial, .7, 13, validate_nonnegative)
    checks['no_rejection_is_bytewise_original_RK4'] = (
        times0.tobytes()==times1.tobytes() and states0.tobytes()==states1.tobytes()
        and not stats['rejected_attempts'] and stats['rhs_calls']==4*13
        and np.array_equal(initial, untouched))
    checks['uniform_macro_frame_and_step_preserved'] = (
        len(stats['accepted_substeps'])==13
        and all(row['depth']==0 and row['dt']==.7/13 for row in stats['accepted_substeps'])
        and [row['start_time'] for row in stats['accepted_substeps']]==times0[:-1].tolist())

    toy_initial = np.array([1.,0.])
    events=[]
    def exchange(t, y):
        return np.array([-4*y[0],4*y[0]])
    times, states, trial = integrate(exchange, toy_initial, 1., 1,
        validate_nonnegative, callback=events.append)
    checks['local_rejection_preserves_positive_state_and_linear_balance'] = (
        len(trial['rejected_attempts'])>0 and np.all(states>=0)
        and abs(states[-1].sum()-1.)<2e-15 and states[-1,0]>0
        and np.array_equal(toy_initial,np.array([1.,0.])))
    checks['rejected_and_accepted_substeps_have_complete_metadata'] = (
        all(set(('start_time','dt','depth','failed_stage','failed_stage_time','violation','rhs_calls')).issubset(row)
            and row['violation']['magnitude']>0 for row in trial['rejected_attempts'])
        and abs(sum(row['dt'] for row in trial['accepted_substeps'])-1.)<1e-15
        and sum(row['event']=='ACCEPT' for row in events)==len(trial['accepted_substeps'])
        and sum(row['event']=='REJECT' for row in events)==len(trial['rejected_attempts'])
        and events[0]['event']=='START' and events[-1]['event']=='STOP')
    details['positive_toy'] = trial

    # No error text matching: a catalogue/physics ValueError propagates exactly.
    marker=ValueError('unsupported catalogue field, not a support trial')
    def invalid_rhs(t,y):
        raise marker
    failed={}
    try:
        integrate(invalid_rhs, np.array([1.]), 1., 1, validate_nonnegative, stats=failed)
    except ValueError as exc:
        checks['non_support_RHS_error_propagates_without_retry'] = (
            exc is marker and failed['rhs_calls']==1 and not failed['rejected_attempts'])

    marker_validator=ValueError('invalid validator dimensions')
    def invalid_validator(y):
        raise marker_validator
    failed={}
    try:
        integrate(exchange, toy_initial, 1., 1, invalid_validator, stats=failed)
    except ValueError as exc:
        checks['non_support_validator_error_propagates_without_retry'] = (
            exc is marker_validator and failed['rhs_calls']==0 and not failed['rejected_attempts'])

    # Even a SupportViolation raised inside the RHS is not a classified result
    # of the supplied state validator and must not trigger a hidden retry.
    rhs_marker=SupportViolation('RHS-originated exception')
    def support_error_in_rhs(t,y):
        raise rhs_marker
    failed={}
    try:
        integrate(support_error_in_rhs, np.array([1.]), 1., 1, validate_nonnegative, stats=failed)
    except SupportViolation as exc:
        checks['only_validator_support_errors_trigger_bisection'] = (
            exc is rhs_marker and failed['rhs_calls']==1 and not failed['rejected_attempts'])

    failed={}
    try:
        integrate(exchange, toy_initial, 1., 1, validate_nonnegative,
                  max_rhs_calls=3, stats=failed)
    except IntegrationLimitError as exc:
        checks['RHS_budget_stops_before_excess_call_with_incomplete_evidence'] = (
            failed['rhs_calls']==3 and failed['completed_macro_steps']==0
            and exc.stats is failed and failed['status']=='INCOMPLETE_LIMIT')
    failed={}
    try:
        integrate(exchange, toy_initial, 1., 1, validate_nonnegative,
                  max_depth=0, stats=failed)
    except IntegrationLimitError as exc:
        checks['depth_budget_stops_without_repair'] = (
            failed['rhs_calls']==1 and failed['max_depth_reached']==0
            and len(failed['rejected_attempts'])==1 and not failed['accepted_substeps'])
    failed={}
    try:
        integrate(exchange, np.array([-1.,2.]), 1., 1, validate_nonnegative, stats=failed)
    except SupportViolation:
        checks['invalid_initial_state_never_repaired_or_integrated'] = (
            failed['rhs_calls']==0 and not failed['accepted_substeps'] and not failed['rejected_attempts'])

    # A validator that permits a deliberately large fourth-stage overshoot
    # but requires a positive endpoint exercises endpoint rejection separately.
    validation_calls=[0]
    def endpoint_test_validator(y):
        validation_calls[0]+=1
        # Initial, k1, k2, k3, k4, then endpoint of the first attempted step.
        if validation_calls[0]==6:
            raise SupportViolation('synthetic endpoint-only failure',{'magnitude':1.,'sector':'test'})
        validate_nonnegative(y)
    _, _, endpoint_stats=integrate(lambda t,y:-.1*y,np.array([1.]),1.,1,endpoint_test_validator)
    checks['endpoint_rejection_reintegrates_only_failed_interval'] = (
        endpoint_stats['rejected_attempts'][0]['failed_stage']=='endpoint'
        and len(endpoint_stats['accepted_substeps'])==2
        and all(row['dt']==.5 for row in endpoint_stats['accepted_substeps']))

    paths=[Path(__file__),Path(__file__).with_name('positive_rk4.py'),
           ROOT/'pysnspd/experimental/cell_validation.py']
    result=dict(schema='pysnspd.stage2.positive-rk4-offline-tests.v1',
        status='PASS' if len(checks)==11 and all(checks.values()) else 'FAIL',
        scope='Synthetic offline numerical guards only; no physical RHS or detector calculation.',
        checks=checks, details=details,
        source_hashes={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        runtime_seconds=time.perf_counter()-started)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x',encoding='utf-8') as stream:
        stream.write(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ('status','checks','source_hashes','runtime_seconds')},indent=2))
    if result['status']!='PASS':
        raise SystemExit(1)


if __name__=='__main__':
    main()
