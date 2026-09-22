"""Inspect two complete first-step trials; never integrate a long trajectory."""
from pathlib import Path
import hashlib,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT),str(ROOT/'sandbox/stage2_cells')]
from run_coupled import setup

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    initial_path=ROOT/'tmp/stage2_validation_20260921_203419/two_electron_4_initial.npz'
    folder=ROOT/'docs/implementation/stage2/recovery_20260921'
    folder.mkdir(parents=True,exist_ok=True)
    output=folder/'first_step_diagnosis.json'
    if output.exists():raise RuntimeError('Do not overwrite a recorded diagnosis.')
    system,y,_=setup('two',electron_refinement=4)
    with np.load(initial_path) as saved:assert np.array_equal(y,saved['state'])
    def inspect(state):
        rows=[]
        for i in range(system.cell_count):
            start=i*system.block_size
            for kind,lower,upper,capacity in (
                ('electron',start+1,start+1+system.electron_size,4*system.catalog.count_weights),
                ('phonon',start+1+system.electron_size,start+system.block_size,system.phonons.capacities)):
                values=state[lower:upper];mask=(values<0)|((values>1)&(kind=='electron'))
                bad=np.flatnonzero(mask)
                rows.append(dict(cell=i,kind=kind,minimum=float(values.min()),maximum=float(values.max()),
                    invalid_count=len(bad),weighted_negative_mass=float(np.dot(capacity,np.maximum(-values,0))),
                    violations=[dict(index=int(k),value=float(values[k]),capacity=float(capacity[k]),
                        coordinate=float(system.catalog.count_nodes[k] if kind=='electron' else system.phonons.energies[k]))
                        for k in bad[np.argsort(values[bad])][:8]]))
        return dict(valid=not any(r['invalid_count'] for r in rows),rows=rows)
    def euler_bound(state,rhs):
        candidates=[]
        for i in range(system.cell_count):
            start=i*system.block_size
            for name,lo,hi in [('electron',start+1,start+1+system.electron_size),
                               ('phonon',start+1+system.electron_size,start+system.block_size)]:
                values,derivative=state[lo:hi],rhs[lo:hi]
                for mode,mask,ratio in [('lower',derivative<0,lambda:-values/derivative),
                    ('upper',(derivative>0)&(name=='electron'),lambda:(1-values)/derivative)]:
                    if mask.any():
                        with np.errstate(divide='ignore',invalid='ignore'):ratios=ratio()
                        k=np.flatnonzero(mask)[np.argmin(ratios[mask])]
                        candidates.append(dict(cell=i,kind=name,bound=mode,index=int(k),
                            limit=float(ratios[k]),value=float(values[k]),derivative=float(derivative[k])))
        return min(candidates,key=lambda r:r['limit']) if candidates else None
    started=time.perf_counter();trials=[]
    for steps in (160,320):
        dt=2/steps;ks=[];records=[];passed=True
        for stage in range(4):
            trial=y if stage==0 else y+dt*ks[-1]*(.5 if stage<3 else 1.)
            state_record=inspect(trial);state_record['stage']=stage+1
            records.append(state_record)
            if not state_record['valid']:passed=False;break
            derivative=system.rhs(dt*(0 if stage==0 else .5 if stage<3 else 1),trial)
            state_record['forward_euler_support_bound']=euler_bound(trial,derivative)
            ks.append(derivative)
        if passed:
            final=y+dt*(ks[0]+2*ks[1]+2*ks[2]+ks[3])/6
            records.append(dict(stage='accepted_endpoint',**inspect(final)))
            passed=records[-1]['valid']
        trials.append(dict(steps=steps,dt=dt,physical_first_step=passed,stages=records))
    result=dict(scope='First-step diagnostic only; no full trajectory or temporal admission',
        trials=trials,rhs_calls=system.rhs_calls,runtime_seconds=time.perf_counter()-started,
        initial_sha256=digest(initial_path),script_sha256=digest(Path(__file__)),
        physical_sources={p.relative_to(ROOT).as_posix():digest(p) for p in sorted((ROOT/'pysnspd/experimental').glob('*.py'))})
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
