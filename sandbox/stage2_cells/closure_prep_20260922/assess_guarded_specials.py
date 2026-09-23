"""Assess newly recorded guarded equilibrium and boundary trajectories."""
from pathlib import Path
import argparse,json,sys
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT),str(ROOT/'sandbox/stage2_cells'),str(Path(__file__).parent)]
from assess_guarded_time import verify_contract,sha

def assess_trajectory(path, task):
    import numpy as np
    from run_coupled import setup
    from pysnspd.experimental.cell_validation import ElectronicCell
    record = json.loads(path.read_text())
    if record['status'] != 'COMPLETED_NOT_YET_ADJUDICATED':
        raise ValueError('Incomplete trajectory')
    if sha(path.with_suffix('.npz')) != record['trajectory_sha256']:
        raise ValueError('Trajectory hash mismatch')
    par = record['parameters']
    for key, wanted in task['parameters'].items():
        actual = par[key]
        if key == 'method':
            if actual != 'ssprk3_common_flux_guarded':
                raise ValueError('Not the registered SSP method')
        elif str(actual) != str(wanted):
            raise ValueError(f'Parameter mismatch{key}: {actual} vs{wanted}')
    for key, digest in record['source_hashes'].items():
        if sha(ROOT/key) != digest:
            raise ValueError('Trajectory source mismatch:'+key)
    system, initial, _ = setup(par['case'], par['phonon_nodes'], par['infrared'],
        par['face_order'], float(par['escape']), par['heating'], par['scenario'],
        par['reaction_order'], par['reaction_method'], par['electron_refinement'],
        par['reaction_layout'], par['reaction_outer_order'], par['reaction_max_panel'],
        par['reaction_max_energy_panel'])
    with np.load(path.with_suffix('.npz')) as stored:
        times, states = stored['times'], stored['states']
        for key, expected in [('electron_count',system.catalog.count_nodes),
            ('electron_weights',system.catalog.count_weights),
            ('phonon_energies',system.phonons.energies),('phonon_capacities',system.phonons.capacities)]:
            if not np.array_equal(stored[key],expected):
                raise ValueError('Grid mismatch:'+key)
    if not np.array_equal(states[0],initial):
        raise ValueError('Actual initial state differs from reconstructed contract')
    gates = dict(ledger=record['energy_ledger_scaled_max']<=1e-7,
        instantaneous=record['instantaneous_residual_max']<=1e-10,
        electron_bounds=record['minimum_electron']>=-1e-13 and record['maximum_electron']<=1+1e-13,
        phonon_bounds=record['minimum_phonon']>=-1e-13,
        no_clipping=record['time_integration']['population_clipping'] is False,
        no_energy_repair=record['time_integration']['posthoc_energy_repair'] is False)
    row=dict(id=task['id'], trajectory=str(path.relative_to(ROOT)),
        trajectory_sha256=sha(path), archive_sha256=sha(path.with_suffix('.npz')),
        gates=gates, minimum_electron=record['minimum_electron'],
        maximum_electron=record['maximum_electron'], minimum_phonon=record['minimum_phonon'],
        energy_ledger_scaled_max=record['energy_ledger_scaled_max'],
        instantaneous_residual_max=record['instantaneous_residual_max'],
        runtime_seconds=record['runtime_seconds'], limiter=record['time_integration']['stats'])
    a0,p0,n0=system.unpack(initial)
    if par['scenario']=='equilibrium':
        rhs=system.rhs(0.,initial)
        drift_budget=par['duration']*1e-11+100*np.finfo(float).eps
        comparisons=[]
        for i in range(system.cell_count):
            start=i*system.block_size
            ep=slice(start+1,start+1+system.electron_size)
            ph=slice(ep.stop,start+system.block_size)
            ew=4*system.catalog.count_weights; pw=system.phonons.capacities
            es=max(1.,float(ew@p0[i])); ps=max(1.,float(pw@n0[i]))
            force=ElectronicCell(system.catalog,a0[i],system.gammas[i]).moments(p0[i])[1]
            observed=dict(cell=i,gamma=system.gammas[i],amplitude=float(a0[i]),
                force_absolute=float(abs(force)),
                amplitude_rhs_absolute=float(abs(rhs[start])),
                electronic_rhs_scaled_L1=float(ew@abs(rhs[ep])/es),
                phonon_rhs_scaled_L1=float(pw@abs(rhs[ph])/ps),
                amplitude_drift_max=float(np.max(abs(states[:,start]-a0[i]))),
                electronic_drift_scaled_L1=float(np.max(abs(states[:,ep]-p0[i])@ew)/es),
                phonon_drift_scaled_L1=float(np.max(abs(states[:,ph]-n0[i])@pw)/ps),
                electron_scale=es,phonon_scale=ps)
            comparisons.append(observed)
        gates['equilibrium_force']=all(v['force_absolute']<=1e-11 for v in comparisons)
        gates['equilibrium_full_rhs']=all(max(v[k] for k in ('amplitude_rhs_absolute','electronic_rhs_scaled_L1','phonon_rhs_scaled_L1'))<=1e-11 for v in comparisons)
        gates['equilibrium_full_trajectory']=all(max(v[k] for k in ('amplitude_drift_max','electronic_drift_scaled_L1','phonon_drift_scaled_L1'))<=drift_budget for v in comparisons)
        row['stationarity']=dict(drift_budget=drift_budget,cells=comparisons,
            all_checkpoint_count=len(times), ledger_absolute_max=float(np.max(abs(states[:,system.population_size:]))))
    else:
        _,pf,nf=system.unpack(states[-1])
        if par['scenario']=='phonon_vacuum':
            gates['exact_initial_phonon_vacuum']=bool(np.all(n0==0))
            gates['nontrivial_phonon_emission']=bool(np.dot(system.phonons.capacities,nf[0])>100*np.finfo(float).eps)
        else:
            gates['exact_initial_empty_and_full_electron_faces']=bool(np.any(p0==0) and np.any(p0==1) and np.all((p0==0)|(p0==1)))
            gates['nontrivial_interior_relaxation']=bool(np.dot(4*system.catalog.count_weights,abs(pf[0]-p0[0]))>100*np.finfo(float).eps)
    row['status']='PASS' if all(gates.values()) else 'FAIL'
    return row

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('runs',nargs=3,type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise ValueError('Preserve prior assessment')
    scenarios=[];rows=[]
    for path in args.runs:
        record=json.loads(path.read_text());verify_contract({'record':record})
        scenarios.append(record['parameters']['scenario'])
        rows.append(assess_trajectory(path,dict(id=path.stem,parameters={k:v for k,v in record['parameters'].items() if k!='output'})))
    if scenarios!=['equilibrium','phonon_vacuum','sparse_electrons']:
        raise ValueError('Exactly the registered equilibrium/vacuum/sparse scenarios, in that order, are required')
    result=dict(schema='pysnspd.stage2.guarded-special-assessment.v1',status='PASS' if all(row['status']=='PASS' for row in rows) else 'FAIL',cases=rows,reviewer_sha256=sha(__file__),stage2_admission=False,scope='Invariants, registered stationarity and nontrivial boundary relaxation; independent field/support checks are separate tasks in the same batch.')
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({k:result[k] for k in ('status','stage2_admission')}),flush=True)
    if result['status']!='PASS':raise SystemExit(1)

if __name__=='__main__':main()
