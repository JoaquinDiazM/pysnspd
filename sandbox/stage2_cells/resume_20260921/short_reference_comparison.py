"""Three short, complete RK4 diagnostics against the user-run DOP853 reference.

Expected total below 90 seconds at the measured 0.17 seconds/RHS. This is
a local integration diagnostic on [0,0.1], not acceptance of the full stage.
The caller applies an external 240-second bound; no automatic retry occurs.
"""
from pathlib import Path
import hashlib,json,os,subprocess,sys,time

ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'docs/implementation/stage2/resume_20260921/short_time'

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    registration=OUT/'registration.json'
    if registration.exists():raise SystemExit('Preserve existing experiment; do not rerun automatically.')
    common=['--case','one','--duration','0.1','--escape','inf','--heating','0',
            '--method','rk4','--phonon-nodes','1025','--electron-refinement','1',
            '--reaction-layout','resolved_panels','--reaction-order','2']
    runs=[OUT/f'one_rk4_{n}.json' for n in (10,20,40)]
    registration.write_text(json.dumps(dict(
        purpose='Local temporal contrast against completed user-run DOP853, not full stage acceptance',
        common_arguments=common,steps=[10,20,40],
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        expected_seconds=90,permitted_external_bound_seconds=240,
        complete_independent_trajectories=True,tolerance_changes=False),indent=2)+'\n')
    env=os.environ|{k:'1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS')}
    started=time.perf_counter()
    for n,path in zip((10,20,40),runs):
        print(json.dumps(dict(event='START',steps=n,elapsed=time.perf_counter()-started)),flush=True)
        subprocess.run([sys.executable,'-u','sandbox/stage2_cells/run_coupled.py',
            *common,'--steps',str(n),'--output',str(path)],cwd=ROOT,env=env,check=True)
    print(json.dumps(dict(event='COMPLETED',elapsed=time.perf_counter()-started)),flush=True)

if __name__=='__main__':main()
