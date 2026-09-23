"""Prepare or explicitly run one bounded guarded-wrapper contract smoke test.

This two-step trajectory tests the real output schema and its validator only.
It is not a temporal, grid, physical or stage2 admission measurement.
"""
from pathlib import Path
import argparse,hashlib,importlib.metadata,json,os,platform,subprocess,sys,time

ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import run_guarded_acceptance_batch as batch

REGISTRATION=ROOT/'docs/implementation/stage2/closure_prep_20260922/guarded_entrypoint_smoke_registration.json'
OUTPUT_ROOT=ROOT/'tmp/stage2_guarded_entrypoint_smoke_20260922'
OUTPUT=OUTPUT_ROOT/'one_guarded_two_steps.json'
GUARDED_SHA='9b77b91ff7b4e68786ea07a7e915f119d330f1c957ca22d2dd7b07028b00fad4'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()


def write_new(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x',encoding='utf-8') as stream:
        json.dump(value,stream,indent=2,allow_nan=False)
        stream.write('\n')


def registration():
    parameters=dict(case='one',steps=2,duration=.002,phonon_nodes=1025,
        infrared=.005,face_order=2,reaction_order=2,reaction_layout='resolved_panels',
        reaction_max_energy_panel=.125,reaction_outer_order=2,reaction_max_panel=.5,
        electron_refinement=1,reaction_method='projected',escape='inf',heating=0.,
        method=batch.METHOD,scenario='driven',rtol=1e-9)
    command=['python','-u',batch.ADAPTER]
    for name,value in parameters.items():
        command+=['--'+name.replace('_','-'),'rk4' if name=='method' else str(value)]
    command+=['--output',OUTPUT.relative_to(ROOT).as_posix()]
    paths=batch.required_contract_paths(dict(tasks=[]))
    paths.append(Path(__file__).relative_to(ROOT).as_posix())
    contracts={path:sha(ROOT/path) for path in sorted(set(paths))}
    if contracts[batch.INTEGRATOR]!=GUARDED_SHA:
        raise ValueError('This smoke registration is bound to the reviewed guarded source')
    return dict(schema='pysnspd.stage2.guarded-entrypoint-smoke.registration.v1',
        status='PREREGISTERED_NOT_EXECUTED',parameters=parameters,
        command_relative=command,source_hashes=contracts,
        output_directory=OUTPUT_ROOT.relative_to(ROOT).as_posix(),
        execution=dict(explicit_execute_required=True,foreground=True,threads=1,
            timeout_seconds=240,estimated_seconds_less_than=10,max_attempts=1),
        validation='Real run_guarded_coupled output, followed by run_guarded_acceptance_batch.validate_trajectory with complete source, parameter, NPZ, initial-state, grid, dependency and invariant checks.',
        final_stage2_admission=False,
        scope='Entrypoint/schema/validator smoke only; no time, grid or physical admission. Not an input to the full acceptance batch.')


def execute():
    expected=registration()
    if json.loads(REGISTRATION.read_text())!=expected:
        raise ValueError('Registration, parameters or sources changed')
    if OUTPUT_ROOT.exists():raise FileExistsError('Preserve previous attempt; no retry or overwrite')
    OUTPUT_ROOT.mkdir(parents=True)
    started=time.perf_counter()
    result=dict(schema='pysnspd.stage2.guarded-entrypoint-smoke.result.v1',
        status='INCOMPLETE',registration_sha256=sha(REGISTRATION),
        source_hashes=expected['source_hashes'],python=platform.python_version(),
        numpy=importlib.metadata.version('numpy'),scipy=importlib.metadata.version('scipy'),
        new_trajectories=0,new_macrosteps=0,physical_admission=False,stage2_admission=False)
    command=[sys.executable,*expected['command_relative'][1:]]
    environment=dict(os.environ)
    environment.update({key:'1' for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS')})
    result['actual_command']=command
    try:
        with (OUTPUT_ROOT/'smoke.console.log').open('x',encoding='utf-8') as log:
            completed=subprocess.run(command,cwd=ROOT,env=environment,stdout=log,
                stderr=subprocess.STDOUT,timeout=240,check=False)
        result['returncode']=completed.returncode
        if completed.returncode:raise RuntimeError('Guarded wrapper failed; inspect the preserved console log')
        result['new_trajectories']=1;result['new_macrosteps']=2
        result['validated_output_hashes']=batch.validate_trajectory(
            OUTPUT,expected['parameters'],expected['source_hashes'])
        if registration()!=expected:raise ValueError('A frozen source changed during the smoke test')
        result['status']='PASS_ENTRYPOINT_CONTRACT_ONLY'
    except Exception as exc:
        result['status']='FAIL_OR_INCOMPLETE';result['reason']=repr(exc)
        raise
    finally:
        result['runtime_seconds']=time.perf_counter()-started
        result['artifacts_sha256']={p.name:sha(p) for p in sorted(OUTPUT_ROOT.iterdir()) if p.is_file()}
        write_new(OUTPUT_ROOT/'smoke_result.json',result)
        print(json.dumps({key:result[key] for key in ('status','runtime_seconds','stage2_admission')},allow_nan=False),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    modes=parser.add_mutually_exclusive_group()
    modes.add_argument('--register',action='store_true')
    modes.add_argument('--execute',action='store_true')
    args=parser.parse_args()
    if args.register:
        write_new(REGISTRATION,registration())
        print('REGISTERED_NO_PHYSICS',REGISTRATION)
    elif args.execute:execute()
    else:
        print(json.dumps(registration(),indent=2))
        print('DRY RUN: no output directory, subprocess, integrator import or trajectory.')


if __name__=='__main__':main()
