"""Synthetic tests of the assessment tool, not physical trajectory evidence."""
from pathlib import Path
import contextlib
import io
import json
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/'sandbox/stage2_cells'))
import assess_time_refinement as assessor


def fixture(path,steps,error,method):
    times=np.linspace(0,1,steps+1)
    states=np.zeros((steps+1,10));states[:,0]=.6+.1*times+error*np.sin(np.pi*times)
    states[:,1:3]=.1;states[:,3:5]=.01
    axes=dict(electron_count=np.array([.1,1.]),electron_weights=np.ones(2)/2,
              phonon_energies=np.array([.2,.4]),phonon_capacities=np.ones(2))
    np.savez(path.with_suffix('.npz'),states=states,times=times,**axes)
    ip=path.with_name(path.stem+'_initial.npz')
    np.savez(ip,state=states[0],**axes)
    snapshots=[]
    for state in states:
        snapshots.append(dict(amplitudes=[float(state[0])],excitation_energy=[1.],
            electron_energy=[1.],phonon_energy=[1.],escape=[0.],input=[0.],
            electron_to_phonon=[0.],condensate_heat=[0.],transport=0.))
    record=dict(status='COMPLETED_NOT_YET_ADJUDICATED',parameters=dict(steps=steps,duration=1.,method=method,rtol=1e-10,output=str(path),heating=.01),
        final=snapshots[-1],snapshots=snapshots,source_hashes={'synthetic_fixture':'not_a_physical_certificate'},
        catalog_sha256='synthetic',criteria_sha256=assessor.sha(ROOT/'docs/implementation/stage2/acceptance_criteria.json'),
        trajectory_sha256=assessor.sha(path.with_suffix('.npz')),initial_conditions_sha256=assessor.sha(ip),
        energy_ledger_scaled_max=0.,instantaneous_residual_max=0.,minimum_electron=.1,maximum_electron=.1,minimum_phonon=.01)
    path.write_text(json.dumps(record))
    return path


def run(paths,reference,output):
    argv=sys.argv
    try:
        sys.argv=['assessor',*map(str,paths),'--reference',str(reference),'--output',str(output)]
        with contextlib.redirect_stdout(io.StringIO()):
            try:assessor.main()
            except SystemExit as exc:
                if exc.code!=1:raise
    finally:sys.argv=argv
    return json.loads(output.read_text())


def main():
    tmp=ROOT/'tmp/stage2_time_assessment_selfcheck';tmp.mkdir(parents=True,exist_ok=True)
    ref=fixture(tmp/'ref.json',4,0.,'dop853')
    convergent=[fixture(tmp/f'good_{steps}.json',steps,.0004*(4/steps)**4,'rk4') for steps in (4,8,16)]
    passed=run(convergent,ref,tmp/'good_result.json')
    bad=[fixture(tmp/f'bad_{steps}.json',steps,.01,'rk4') for steps in (4,8,16)]
    rejected=run(bad,ref,tmp/'bad_result.json')
    gates=dict(refined_time_error_passes=passed['status']=='PASS',
        interior_error_rejected=rejected['status']=='FAIL',
        equal_final_values_do_not_mask_interior=all(max(v['assessment_error'] for v in r['final_errors'].values())<1e-12 for r in rejected['cases']))
    gates['zero_reference_nonroundoff_difference_is_undefined']=assessor.measurement_error(1e-6,0.)['assessment_error'] is None
    gates['zero_reference_roundoff_difference_is_zero']=assessor.measurement_error(1e-16,0.)['assessment_error']==0.
    good=assessor.load_run(convergent[0]);reference=assessor.load_run(ref)
    good['record']['parameters']['heating']=.02
    try:assessor.require_same_problem(good,reference)
    except ValueError:gates['physical_parameter_mismatch_rejected']=True
    else:gates['physical_parameter_mismatch_rejected']=False
    good=assessor.load_run(convergent[0]);good['phonon_energies']=good['phonon_energies']+.001
    try:assessor.require_same_problem(good,reference)
    except ValueError:gates['equal_capacities_wrong_energy_axis_rejected']=True
    else:gates['equal_capacities_wrong_energy_axis_rejected']=False
    tampered=convergent[0].with_suffix('.npz')
    tampered.write_bytes(tampered.read_bytes()+b'changed')
    try:assessor.load_run(convergent[0])
    except ValueError:gates['NPZ_hash_tampering_rejected']=True
    else:gates['NPZ_hash_tampering_rejected']=False
    result=dict(status='PASS' if all(gates.values()) else 'FAIL',
        scope='Synthetic assessment-tool regressions only, no physical dynamics claim',gates=gates,
        good_finest_error=passed['finest_max_error'],good_reduction=passed['successive_max_error_reduction'],
        bad_interior_error=rejected['finest_max_error'],
        assessor_sha256=assessor.sha(ROOT/'sandbox/stage2_cells/assess_time_refinement.py'),
        selfcheck_sha256=assessor.sha(Path(__file__)))
    (Path(__file__).resolve().parent/'time_assessment_selfcheck.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))
    if not all(gates.values()):raise SystemExit(1)


if __name__=='__main__':main()
