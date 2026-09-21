"""Compare fresh coupled runs with a separate integrator/reference trajectory."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
ROOT = Path(__file__).resolve().parents[2]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_run(path):
    record = json.loads(path.read_text())
    if record.get('status') != 'COMPLETED_NOT_YET_ADJUDICATED':
        raise ValueError('time comparison requires completed trajectories')
    if sha(path.with_suffix('.npz')) != record.get('trajectory_sha256'):
        raise ValueError('trajectory NPZ does not match its recorded hash')
    initial_path = path.with_name(path.stem+'_initial.npz')
    if sha(initial_path) != record.get('initial_conditions_sha256'):
        raise ValueError('initial-condition NPZ does not match its recorded hash')
    with np.load(path.with_suffix('.npz')) as archive:
        values = {k: archive[k] for k in ('states','times','electron_count',
            'electron_weights','phonon_energies','phonon_capacities')}
    with np.load(initial_path) as archive:
        if not np.array_equal(archive['state'], values['states'][0]):
            raise ValueError('trajectory does not start at its declared initial state')
        for name in ('electron_count','electron_weights','phonon_energies','phonon_capacities'):
            if not np.array_equal(archive[name], values[name]):
                raise ValueError('trajectory grid changed from its declared initial grid')
    cells = len(record['final']['amplitudes'])
    block = 1+len(values['electron_weights'])+len(values['phonon_capacities'])
    if (values['states'].shape != (len(values['times']),cells*block+4*cells+1)
            or len(values['times']) != record['parameters']['steps']+1
            or len(record['snapshots']) != len(values['times'])
            or np.any(~np.isfinite(values['states']))
            or np.any(~np.isfinite(values['times']))
            or np.any(np.diff(values['times']) <= 0)
            or values['times'][0] != 0
            or values['times'][-1] != record['parameters']['duration']):
        raise ValueError('inconsistent trajectory shape, times or saved snapshots')
    values.update(record=record,path=path,cells=cells,block=block)
    return values


def require_same_problem(actual, reference):
    excluded={'steps','method','rtol','output'}
    ap={k:v for k,v in actual['record']['parameters'].items() if k not in excluded}
    rp={k:v for k,v in reference['record']['parameters'].items() if k not in excluded}
    if ap != rp:
        raise ValueError('time comparison requires identical physical and grid parameters')
    for name in ('electron_count','electron_weights','phonon_energies','phonon_capacities'):
        if not np.array_equal(actual[name],reference[name]):
            raise ValueError('time comparison requires identical physical grid axes and capacities')
    if not np.array_equal(actual['states'][0],reference['states'][0]):
        raise ValueError('time comparison requires identical initial populations and ledgers')
    for name in ('source_hashes','catalog_sha256','criteria_sha256'):
        if actual['record'].get(name) != reference['record'].get(name):
            raise ValueError('time comparison requires matching '+name)


def shared_indices(runs, reference):
    """Use exact common output times; never interpolate numerical populations."""
    tolerance=32*np.finfo(float).eps*max(1.,abs(reference['times'][-1]))
    indices=[]
    for j,time in enumerate(reference['times']):
        matches=[]
        for run in runs:
            k=int(np.argmin(abs(run['times']-time)))
            if abs(run['times'][k]-time)>tolerance:
                break
            matches.append(k)
        else:
            indices.append((j,matches))
    if len(indices)<3 or indices[0][0]!=0 or indices[-1][0]!=len(reference['times'])-1:
        raise ValueError('at least three common checkpoints including both endpoints are required')
    return indices


def measurement_error(difference, scale):
    floor=float(100*np.finfo(float).eps*max(1.,scale))
    degenerate=scale<=floor
    if degenerate:
        relative=None
        gate=0. if difference<=floor else None
        convention='Degenerate reference: only roundoff agreement is accepted; otherwise an independently declared activity scale is required'
    else:
        relative=difference/scale
        gate=0. if difference<=floor else relative
        convention='Difference normalized by the reference observable norm'
    return dict(absolute=difference,scale=scale,relative=relative,
        assessment_error=gate,zero_reference=degenerate,roundoff_limited=difference<=floor,
        comparison_defined=gate is not None,convention=convention)


def checkpoint_errors(run, k, reference, j, fields):
    actual_snapshot=run['record']['snapshots'][k]
    expected_snapshot=reference['record']['snapshots'][j]
    errors={}
    for name in fields:
        actual,expected=np.asarray(actual_snapshot[name]),np.asarray(expected_snapshot[name])
        difference=float(np.max(abs(actual-expected)))
        scale=float(np.max(abs(expected)))
        errors[name]=measurement_error(difference,scale)
    block=run['block'];ne=len(run['electron_weights'])
    for name,low,high,mass in (('electronic_population',1,1+ne,4*run['electron_weights']),
                              ('phonon_population',1+ne,block,run['phonon_capacities'])):
        actual=np.asarray([run['states'][k,i*block+low:i*block+high] for i in range(run['cells'])])
        expected=np.asarray([reference['states'][j,i*block+low:i*block+high] for i in range(run['cells'])])
        difference=float(np.sum(mass*abs(actual-expected)))
        scale=float(np.sum(mass*abs(expected)))
        errors[name]=measurement_error(difference,scale)
    return errors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('runs', nargs='+', type=Path)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    reference=load_run(args.reference)
    runs=[load_run(path) for path in args.runs]
    for run in runs:require_same_problem(run,reference)
    common=shared_indices(runs,reference)
    criteria_path = ROOT/'docs/implementation/stage2/acceptance_criteria.json'
    criteria = json.loads(criteria_path.read_text())
    if reference['record']['criteria_sha256']!=sha(criteria_path):
        raise ValueError('recorded acceptance criteria differ from the criteria being applied')
    fields = ('amplitudes','excitation_energy','electron_energy','phonon_energy',
              'escape','input','electron_to_phonon','condensate_heat','transport')
    cases = []
    for index,run in enumerate(runs):
        path,record=run['path'],run['record']
        checkpoints=[];errors={}
        for j,matches in common:
            time=float(reference['times'][j])
            measured=checkpoint_errors(run,matches[index],reference,j,fields)
            checkpoints.append(dict(time=time,errors=measured))
            for name,value in measured.items():
                score=lambda entry:float('inf') if entry['assessment_error'] is None else entry['assessment_error']
                if name not in errors or score(value)>score(errors[name]):
                    errors[name]=value|{'time_of_maximum':time}
        defined=all(v['comparison_defined'] for v in errors.values())
        cases.append(dict(path=path.as_posix(),sha256=sha(path),steps=record['parameters']['steps'],
                          errors=errors, max_error=max(v['assessment_error'] for v in errors.values()) if defined else None,
                          comparison_defined=defined,
                          checkpoints=checkpoints,final_errors=checkpoints[-1]['errors'],
                          ledger=record['energy_ledger_scaled_max'],
                          minimum_electron=record['minimum_electron'], maximum_electron=record['maximum_electron'],
                          minimum_phonon=record['minimum_phonon']))
    cases.sort(key=lambda v:v['steps'])
    if len({case['steps'] for case in cases}) != len(cases):
        raise ValueError('time convergence requires distinct step counts')
    ratios = []
    for previous, current in zip(cases[:-1],cases[1:]):
        ratios.append(previous['max_error']/current['max_error']
            if current['max_error'] not in (0.,None) and previous['max_error'] is not None else None)
    finest = cases[-1]
    tolerance = criteria['coupled_trajectories']['finest_time_observable_relative_error_max']
    defined=all(c['comparison_defined'] for c in cases)
    resolved = finest['max_error'] is None or finest['max_error'] > 100*np.finfo(float).eps
    convergence = defined and (not resolved or all(r is None or r >= 1.5 for r in ratios))
    ref=reference['record']
    reference_checks=dict(ledger=ref['energy_ledger_scaled_max']<=1e-7,
        instantaneous_balance=ref['instantaneous_residual_max']<=1e-10,
        populations=ref['minimum_electron']>=0 and ref['maximum_electron']<=1 and ref['minimum_phonon']>=0)
    passed = (defined and finest['max_error'] <= tolerance and convergence and len(cases)>=3
              and all(reference_checks.values())
              and all(c['ledger']<=1e-7 and c['minimum_electron']>=0 and c['maximum_electron']<=1
                      and c['minimum_phonon']>=0 for c in cases))
    result = dict(status='PASS' if passed else 'FAIL',criteria_sha256=sha(criteria_path),
                  reference_path=args.reference.as_posix(),reference_sha256=sha(args.reference),
                  reference_checks=reference_checks,
                  comparison='Maximum over shared physical checkpoints; endpoints retained separately; no population interpolation',
                  shared_checkpoint_times=[float(reference['times'][j]) for j,_ in common],
                  reviewer_sha256=sha(__file__),cases=cases, successive_max_error_reduction=ratios,
                  finest_max_error=finest['max_error'], tolerance=tolerance,
                  all_comparisons_defined=defined,
                  converged=convergence,roundoff_limited=not resolved)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:result[k] for k in ('status','finest_max_error','successive_max_error_reduction')},indent=2))
    if not passed: raise SystemExit(1)


if __name__=='__main__':
    main()
