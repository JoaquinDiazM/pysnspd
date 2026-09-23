"""Recorded SSPRK3 experiment with a conservative subnormal-arithmetic guard."""
from pathlib import Path
import hashlib,json,sys,time
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT),str(ROOT/'sandbox/stage2_cells'),str(Path(__file__).parent)]
import run_coupled
import limited_ssp_guarded as guarded
import limited_ssp as legacy

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main():
    if '--output' not in sys.argv:raise SystemExit('Explicit new output required.')
    output=Path(sys.argv[sys.argv.index('--output')+1])
    progress=output.with_suffix('.integration.jsonl')
    if any(p.exists() for p in (output,progress,output.with_suffix('.npz'),output.with_name(output.stem+'_initial.npz'))):
        raise SystemExit('Existing evidence is preserved; choose a new output path.')
    if '--method' in sys.argv and sys.argv[sys.argv.index('--method')+1]!='rk4':
        raise SystemExit('Pass --method rk4 as the frozen runner entry; output declares SSPRK3 explicitly.')
    output.parent.mkdir(parents=True,exist_ok=True)
    stats={};started=time.perf_counter();original=run_coupled.rk4_trajectory
    with progress.open('x',encoding='utf-8',buffering=1) as stream:
        def emit(event):
            row=dict(event,elapsed_seconds=time.perf_counter()-started)
            line=json.dumps(row,allow_nan=False)
            stream.write(line+'\n');print(line,flush=True)
        def replacement(rhs,initial,duration,steps):
            times,states,_=guarded.integrate(rhs.__self__,initial,duration,steps,
                callback=emit,stats=stats)
            return times,states
        run_coupled.rk4_trajectory=replacement
        try:
            run_coupled.main()
        finally:
            run_coupled.rk4_trajectory=original
            if output.exists():
                record=json.loads(output.read_text())
                record['parameters']['method']='ssprk3_common_flux_guarded'
                record['time_integration']=dict(method='SSPRK3 with common conservative event-flux factors and subnormal arithmetic guard',
                    stats=stats,population_clipping=False,posthoc_energy_repair=False,
                    scope='Corrective numerical experiment; requires separate temporal and grid admission')
                for source in (Path(__file__),Path(guarded.__file__),Path(legacy.__file__)):
                    record['source_hashes'][source.relative_to(ROOT).as_posix()]=sha(source)
                output.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n',encoding='utf-8')
                emit(dict(event='INTEGRATION_RECORD_SAVED',status=record['status'],output=str(output)))

if __name__=='__main__':main()
