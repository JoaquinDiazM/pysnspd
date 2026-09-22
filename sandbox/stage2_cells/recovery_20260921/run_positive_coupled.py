"""Frozen cell equations with explicitly recorded support-rejecting RK4.

This adapter changes only time integration. It never edits populations or the
energy ledger. Records identify both the frozen physical runner and this adapter.
"""
from pathlib import Path
import argparse,hashlib,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT),str(ROOT/'sandbox/stage2_cells'),str(Path(__file__).parent)]
import run_coupled
import positive_rk4

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser(add_help=False)
    parser.add_argument('--integration-max-depth',type=int,default=20)
    parser.add_argument('--integration-max-rhs',type=int,default=20000)
    controls,remaining=parser.parse_known_args()
    if '--output' not in remaining:raise SystemExit('An explicit new --output is required.')
    output=Path(remaining[remaining.index('--output')+1])
    progress=output.with_suffix('.integration.jsonl')
    if any(p.exists() for p in (output,progress,output.with_suffix('.npz'),output.with_name(output.stem+'_initial.npz'))):
        raise SystemExit('Preserve existing outputs; choose a new path.')
    if '--method' in remaining and remaining[remaining.index('--method')+1]!='rk4':
        raise SystemExit('This adapter implements RK4 only.')
    output.parent.mkdir(parents=True,exist_ok=True)
    started=time.perf_counter();stats={};original=run_coupled.rk4_trajectory
    with progress.open('x',encoding='utf-8',buffering=1) as stream:
        def emit(event):
            row=dict(event,elapsed_seconds=time.perf_counter()-started)
            line=json.dumps(row,allow_nan=False)
            stream.write(line+'\n')
            if row.get('event')!='ACCEPT' or len(stats.get('accepted_substeps',[]))%25==0:
                print(line,flush=True)
        def replacement(rhs,initial,duration,steps):
            system=rhs.__self__
            def validate(state):
                if state.shape!=initial.shape or np.any(~np.isfinite(state)):
                    raise ValueError('Nonfinite or wrong-shaped state is not a support-rejection event.')
                for cell in range(system.cell_count):
                    offset=cell*system.block_size
                    for kind,lo,hi in (
                        ('electron',offset+1,offset+1+system.electron_size),
                        ('phonon',offset+1+system.electron_size,offset+system.block_size)):
                        values=state[lo:hi]
                        bad=np.flatnonzero((values<0)|((values>1)&(kind=='electron')))
                        if len(bad):
                            k=int(bad[np.argmax(np.maximum(-values[bad],values[bad]-1 if kind=='electron' else 0))])
                            raise positive_rk4.SupportViolation('Population outside physical support',
                                details=dict(cell=cell,sector=kind,index=k,value=float(values[k]),
                                    magnitude=float(max(-values[k],values[k]-1 if kind=='electron' else 0))))
            times,states,_=positive_rk4.integrate(rhs,initial,duration,steps,validate,
                max_depth=controls.integration_max_depth,max_rhs_calls=controls.integration_max_rhs,
                callback=emit,stats=stats)
            return times,states
        sys.argv=[sys.argv[0],*remaining]
        run_coupled.rk4_trajectory=replacement
        try:
            run_coupled.main()
        finally:
            run_coupled.rk4_trajectory=original
            if output.exists():
                record=json.loads(output.read_text())
                record['parameters']['method']='rk4_support_bisection'
                record['time_integration']=dict(method='Classical RK4; reject entire trial and bisect locally on support violation',
                    controls=vars(controls),stats=stats,clipping=False,energy_repair=False,
                    max_macro_step=record['parameters']['duration']/record['parameters']['steps'],
                    formal_time_admission='PENDING_FOR_REJECTION_ACTIVE_TRAJECTORIES')
                for source in (Path(__file__),Path(positive_rk4.__file__)):
                    record['source_hashes'][source.relative_to(ROOT).as_posix()]=sha(source)
                output.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n',encoding='utf-8')
                emit(dict(event='INTEGRATION_RECORD_SAVED',status=record['status'],
                    rhs_calls=stats.get('rhs_calls'),output=str(output)))

if __name__=='__main__':main()
