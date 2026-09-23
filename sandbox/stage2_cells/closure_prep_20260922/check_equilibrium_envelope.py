"""Check actual equilibrium upper support by an exact causal spectral envelope.

No trajectory integration and no continuum quadrature. The proof and the
earlier independent quadrature comparison are retained separately. This entry
point recomputes gross in-support activity from the input trajectory itself.
"""
from pathlib import Path
import argparse,hashlib,json,sys,time
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT),str(ROOT/'sandbox/stage2_cells')]
import numpy as np
from assess_grid_refinement import read
from run_coupled import setup
from pysnspd.experimental.cell_validation import ElectronicCell
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('trajectories',nargs='+',type=Path);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise ValueError('Do not overwrite evidence')
    started=time.perf_counter();rows=[];sampling=[];inputs={}
    criteria=json.loads((ROOT/'docs/implementation/stage2/acceptance_criteria.json').read_text())
    limit=criteria['continuous_consistency']['upper_support']['relative_omitted_power_and_population_activity_max']
    minimum=criteria['catalogue_field_checks']['minimum_actual_trajectory_samples']
    for path in args.trajectories:
        loaded=read(path,ROOT);record=loaded['record'];par=record['parameters'];inputs[str(path)]=sha(path)
        if par['scenario']!='equilibrium':raise ValueError('This frozen application is the equilibrium case')
        system,_,debye=setup(par['case'],par['phonon_nodes'],par['infrared'],par['face_order'],float(par['escape']),par['heating'],par['scenario'],par['reaction_order'],par['reaction_method'],par['electron_refinement'],par['reaction_layout'],par['reaction_outer_order'],par['reaction_max_panel'],par['reaction_max_energy_panel'])
        times=loaded['arrays']['times'];states=loaded['arrays']['states']
        indices=np.unique(np.linspace(0,len(times)-1,11,dtype=int))
        if len(np.unique(times[indices]))<minimum:raise ValueError('Insufficient actual times')
        sampling.append(dict(trajectory=str(path),times_per_cell=len(indices),cells=system.cell_count))
        cache={};cut=debye.cutoff_energy_bar;lo=debye.infrared_cutoff_bar
        for index in indices:
            amplitudes,p,n=system.unpack(states[index])
            for i,a in enumerate(amplitudes):
                cell=ElectronicCell(system.catalog,a,system.gammas[i]);emax=float(cell.energies[-1]);emin=emax-cut
                if emin<=a:raise ValueError('The analytic bound requires minimum tail energy above amplitude')
                key=(float(a),system.gammas[i],p[i].tobytes(),n[i].tobytes())
                if key not in cache:
                    events=system.reaction_events(cell);lf,lb=events.log_activities(p[i],n[i])
                    gross=events.coefficients*(np.exp(lf)+np.exp(lb))
                    denom=dict(number=float(gross.sum()),power=float(gross@events.omega))
                    pmax=float(max(np.interp(emin,cell.energies,p[i]),np.max(p[i][cell.energies>=emin])))
                    nmax=float(np.max(n[i]));coherence=emin**2/(emin**2-a*a)
                    factor=system.rate_prefactor*coherence*pmax*nmax*debye.lambda_eph/cut**2
                    bounds=dict(number=factor*(cut**4-lo**4)/4,power=factor*(cut**5-lo**5)/5)
                    if any(v<=0 for v in denom.values()):raise ValueError('Undefined gross activity denominator')
                    cache[key]=dict(relative_upper_bounds={k:bounds[k]/denom[k] for k in bounds},
                        absolute_upper_bounds=bounds,gross_activity=denom,population_maxima=[pmax,nmax],
                        minimum_tail_energy=emin,coherence_upper_bound=coherence)
                rows.append(dict(trajectory=str(path),time=float(times[index]),cell=i,amplitude=float(a),gamma=system.gammas[i],**cache[key]))
    worst=max(v for row in rows for v in row['relative_upper_bounds'].values())
    proof=ROOT/'docs/implementation/stage2/closure_prep_20260922/EQUILIBRIUM_SUPPORT_ENVELOPE.md'
    result=dict(status='PASS' if worst<=limit else 'FAIL',worst_relative_upper_bound=worst,relative_limit=limit,
        sampling=sampling,samples=len(rows),rows=rows,inputs=inputs,source_sha256=sha(__file__),proof_sha256=sha(proof),
        criteria_sha256=sha(ROOT/'docs/implementation/stage2/acceptance_criteria.json'),runtime_seconds=time.perf_counter()-started,
        finite_eta_included_in_causal_bound=True,quadrature_refinement_required=False,new_trajectories=0,stage2_admission=False,
        scope='Absorption to unresolved high-energy states. No claim about arbitrary unrepresented incoming populations.')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2,allow_nan=False)
    print(json.dumps({k:result[k] for k in ('status','samples','worst_relative_upper_bound','runtime_seconds')}))
    if result['status']!='PASS':raise SystemExit(1)

if __name__=='__main__':main()
