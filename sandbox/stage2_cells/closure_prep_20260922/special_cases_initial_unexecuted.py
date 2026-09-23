"""Registered short SSP equilibrium/boundary trajectories, with no kernel edits.

The caller applies ONE 240-second process-group timeout to the whole block.
--prepare only writes a preregistration; --execute preserves a fresh output root.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
CRITERIA = ROOT/'docs/implementation/stage2/acceptance_criteria.json'
FROZEN = ROOT/'docs/implementation/stage2/time_pass_20260922/two_cell_time_plan.json'
RUNNER = ROOT/'sandbox/stage2_cells/recovery_20260921/run_limited_coupled.py'
sys.path[:0] = [str(ROOT), str(ROOT/'sandbox/stage2_cells'),
               str(ROOT/'docs/implementation/stage2/review'),
               str(ROOT/'docs/implementation/stage1_r2/review')]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def registration():
    frozen = json.loads(FROZEN.read_text())
    template = dict(frozen['tasks'][0]['parameters'])
    template['method'] = 'rk4'  # Frozen wrapper entry; metadata becomes actual SSP.
    tasks = []
    for name, case, scenario, duration, steps, escape in (
        ('special_equilibrium_two', 'two', 'equilibrium', .1, 10, 15.),
        ('special_phonon_vacuum_one', 'one', 'phonon_vacuum', .2, 40, 'inf'),
        ('special_sparse_electrons_one', 'one', 'sparse_electrons', .2, 40, 'inf'),
    ):
        par = dict(template, case=case, scenario=scenario, duration=duration,
                   steps=steps, escape=escape, heating=0.)
        tasks.append(dict(id=name, parameters=par))
    sources = dict(frozen['source_hashes'])
    sources[Path(__file__).relative_to(ROOT).as_posix()] = sha(__file__)
    for ref in ('docs/implementation/stage1_r2/review/causal_reference.py',
                'docs/implementation/stage2/review/complementary_energy_reference.py'):
        sources[ref] = sha(ROOT/ref)
    return dict(schema='pysnspd.stage2.special-short-cases.registration.v1',
        status='PREREGISTERED_NOT_EXECUTED', tasks=tasks,
        source_hashes=sources, criteria_sha256=sha(CRITERIA),
        physical_kernels_changed=False, thresholds_changed=False,
        inherited_gates=dict(ledger_scaled_max=1e-7, instantaneous_scaled_max=1e-10,
            pauli_roundoff_max=1e-13, phonon_negative_roundoff_max=1e-13,
            equilibrium_weighted_rhs_scaled_max=1e-11,
            field_moment_scaled_or_relative_max=1e-3, derivative_identity_scaled_max=1e-7),
        stationarity=dict(
            norm='Per-cell electron L1 with4wx; phonon L1 with phonon capacities; amplitude absolute. Population scales=max(1,initial weighted count).',
            drift_budget='duration*1e-11 +100*machine_epsilon, after the declared normalization',
            amplitude_force_absolute_max=1e-11,
            note='A stationary implemented force is distinct from an independently referenced physical force. Near-zero physical force uses scale max(1,sum(abs(vacuum and occupation contributions))).'),
        field_reference=dict(
            method='Independent eta=0 real-angle reference at all native count nodes, same populations and weights; separate60-digit w-plane finite-eta checks at9 count nodes per actual finite-Gamma field.',
            samples='All11 distinct equilibrium times in both cells; zero-Gamma boundary fields use the existing analytic BCS assessor.',
            scope='The ideal reference comparison includes finite-eta bias. A sampled finite-eta spectral check alone is not a full moment or upper-support certificate.'),
        support_reference=dict(method='Independent high-energy fixed-point Usadel solve, original residual checked, integration in(E,Omega), external hole bounded by1. Orders16/32; actual populations, no extrapolation.',
            domain='Only energies greater than2*max(amplitude,Gamma,eta), which contains the registered unresolved upper absorption sector.',
            relative_bound_max=1e-3, quadrature_relative_difference_max=2e-4),
        limitations=['No time-refinement certificate for these short special trajectories.',
            'Upper support bounds outgoing absorption; no arbitrary populated external electron reservoir is assumed.',
            'No global stage2 or production admission follows from this block.'],
        execution=dict(expected_seconds=[100, 170], bounded_timeout_seconds=240,
            threads=1, account='jdiaz', no_retry_or_splitting=True,
            stop_on_first_failure=True))


def check_sources(plan):
    for path, expected in plan['source_hashes'].items():
        if sha(ROOT/path) != expected:
            raise ValueError('Source differs from registered contract: '+path)
    if sha(CRITERIA) != plan['criteria_sha256']:
        raise ValueError('Criteria changed')


def run_command(command, log):
    with Path(log).open('x', encoding='utf-8') as stream:
        completed = subprocess.run(command, cwd=ROOT, stdout=stream,
                                   stderr=subprocess.STDOUT, check=False)
    if completed.returncode:
        raise RuntimeError(f'Command failed({completed.returncode}); preserved log: {log}')


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
            if actual != 'ssprk3_common_flux_limited':
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


def equilibrium_fields(path, output):
    """Independent eta0 moments plus sparse finite-eta implicit spectral checks."""
    import numpy as np
    from causal_reference import parametric_state,t_at_count,vacuum_reference
    from complementary_energy_reference import reference,mp
    from pysnspd.experimental.energy_catalog import OccupationEnergyCatalog,retarded_spectrum
    from pysnspd.experimental.refined_cells import ComplementaryCountCatalog
    started=time.perf_counter();mp.mp.dps=60
    data=json.loads(path.read_text())
    with np.load(path.with_suffix('.npz')) as file:
        states,times=file['states'],file['times']
        x,w=file['electron_count'],file['electron_weights']
        nph=len(file['phonon_energies'])
    base=OccupationEnergyCatalog.load(ROOT/'docs/implementation/stage1_r2/catalogs/occupation_catalog.npz')
    cat=ComplementaryCountCatalog(base,x,w);block=1+len(x)+nph
    indices=np.unique(np.linspace(0,len(x)-1,9,dtype=int));rows=[];spectral=[];cache={}
    for time_value,state in zip(times,states):
        for i,g in enumerate((0.,.1)):
            a=float(state[i*block]);p=state[i*block+1:i*block+1+len(x)]
            if (a,g) not in cache:
                if g==0:
                    den=x*x+base.eta**2
                    e=x*np.sqrt(1+a*a/den)
                    da=x*a/np.sqrt(den*(den+a*a))
                    dg=-a*a*base.eta*x/(den*den+a*a*base.eta**2)
                    independent=np.array([e,da,dg])
                else:
                    independent=np.array([parametric_state(t_at_count(float(xx),a,g),a,g)[1:4] for xx in x]).T
                    actual=np.asarray(cat.energy_kernel(a,g))
                    c,s=retarded_spectrum(actual[0,indices],delta=a,gamma=g,eta=base.eta)
                    for node,seed in zip(indices,1j*a/s):
                        wanted,residual=reference(float(x[node]),a,g,base.eta,seed)
                        err=abs(actual[:,node]-wanted)/np.maximum(1.,abs(wanted))
                        spectral.append(dict(amplitude=a,gamma=g,count=float(x[node]),
                            scaled_errors=err.tolist(),reference_residual=residual))
                cache[(a,g)]=independent
            kernel=cache[(a,g)];vac=vacuum_reference(a,g)
            exc=4*np.sum(kernel*(w*p)[None,:],axis=1)
            wanted=vac+exc;actual=np.array(cat.evaluate(a,g,p))
            scales=np.maximum(1.,abs(vac)+abs(exc))
            errors=abs(actual-wanted)/scales
            h=1e-5
            fd=(cat.evaluate(a+h,g,p)[0]-cat.evaluate(a-h,g,p)[0])/(2*h)
            identity=abs(fd-actual[1])/max(1.,abs(actual[1]))
            rows.append(dict(time=float(time_value),cell=i,amplitude=a,gamma=g,
                actual=actual.tolist(),independent_reference=wanted.tolist(),
                scales=scales.tolist(),scaled_errors=errors.tolist(),derivative_identity_scaled=float(identity)))
    maxima=dict(moment_scaled_max=max(max(row['scaled_errors']) for row in rows),
        derivative_identity_scaled_max=max(row['derivative_identity_scaled'] for row in rows),
        finite_eta_kernel_scaled_max=max(max(row['scaled_errors']) for row in spectral))
    passed=(maxima['moment_scaled_max']<=1e-3 and maxima['derivative_identity_scaled_max']<=1e-7
            and maxima['finite_eta_kernel_scaled_max']<=1e-7)
    result=dict(status='PASS' if passed else 'FAIL',scope='Actual equilibrium fields only; all native-node moment comparisons against independent real-angle eta0 reference. Finite-eta spectral kernels additionally sampled with60-digit implicit reference; ideal comparison includes regulator bias.',
        sampling=dict(distinct_times_per_cell=len(np.unique(times)),cells=2,unique_fields=len(cache)),
        maxima=maxima,rows=rows,finite_eta_spectral_samples=spectral,
        source_sha256=sha(__file__),trajectory_sha256=sha(path),runtime_seconds=time.perf_counter()-started)
    save(output,result)
    if not passed:raise RuntimeError('Independent equilibrium field check failed')
    return {k:result[k] for k in ('status','scope','sampling','maxima','runtime_seconds')}


def equilibrium_support(path, output):
    """Bound missing absorption with an independent high-energy causal solve."""
    import numpy as np
    from run_coupled import setup
    from pysnspd.experimental.cell_validation import ElectronicCell
    started=time.perf_counter();data=json.loads(path.read_text());par=data['parameters']
    system,_,debye=setup(par['case'],par['phonon_nodes'],par['infrared'],par['face_order'],
        float(par['escape']),par['heating'],par['scenario'],par['reaction_order'],
        par['reaction_method'],par['electron_refinement'],par['reaction_layout'],
        par['reaction_outer_order'],par['reaction_max_panel'],par['reaction_max_energy_panel'])
    with np.load(path.with_suffix('.npz')) as file:
        times,states=file['times'],file['states']
    spectral_residuals=[]

    def spectrum(energy,a,g):
        if np.any(energy<=2*max(a,g,system.catalog.eta)):
            raise ValueError('Outside high-energy contraction domain')
        z=energy+1j*system.catalog.eta;u=z.copy()
        for iteration in range(40):
            w=np.sqrt(u*u-a*a)
            if np.any(w.real<=0) or np.any(w.imag<=0):
                raise ValueError('Noncausal independent high-energy branch')
            updated=z+1j*g*u/w
            change=np.max(abs(updated-u));u=updated
            if change<5e-15:break
        else:raise ArithmeticError('Independent high-energy solve did not converge')
        w=np.sqrt(u*u-a*a);c=u/w;s=1j*a/w
        residual=float(np.max(abs(a*c-(g*c-1j*z)*s)/np.maximum(1.,energy)))
        if residual>1e-11 or np.any(c.real<=0) or np.any(s.imag<0):
            raise ArithmeticError('Independent high-energy spectral residual/branch failed')
        spectral_residuals.append(residual)
        return c.real,s.imag

    def bound(a,g,emax,energies,p,n,order):
        cuts=[par['infrared'],.02,.05,.1,.2,.5,1.,1.5,2.,2.5,3.,4.]
        cuts=sorted(set(v for v in cuts if par['infrared']<=v<=4.))
        z,w=np.polynomial.legendre.leggauss(order)
        om=np.concatenate([lo+(z+1)*(hi-lo)/2 for lo,hi in zip(cuts[:-1],cuts[1:])])
        wo=np.concatenate([w*(hi-lo)/2 for lo,hi in zip(cuts[:-1],cuts[1:])])
        inside=emax-om[:,None]+om[:,None]*(z[None,:]+1)/2;outside=inside+om[:,None]
        if np.min(inside)<energies[0] or np.max(inside)>energies[-1]:
            raise ValueError('Tail would extrapolate electron populations')
        if np.min(om)<system.phonons.energies[0] or np.max(om)>system.phonons.energies[-1]:
            raise ValueError('Tail would extrapolate phonons')
        ni,ri=spectrum(inside,a,g);no,ro=spectrum(outside,a,g)
        coherence=ni*no-ri*ro
        if np.any(coherence<0):raise ValueError('Negative independent coherence')
        weighted=wo[:,None]*(om[:,None]*w[None,:]/2)*coherence*debye.alpha2F(om)[:,None]
        pp=np.interp(inside,energies,p);nn=np.interp(om,system.phonons.energies,n)
        rate=weighted*pp*nn[:,None]
        return np.array([np.sum(rate),np.sum(rate*om[:,None])])*system.rate_prefactor

    rows=[]
    for time_value,state in zip(times,states):
        aa,pp,nn=system.unpack(state)
        for i,(a,g) in enumerate(zip(aa,system.gammas)):
            cell=ElectronicCell(system.catalog,a,g);event=system.reaction_events(cell)
            lf,lb=event.log_activities(pp[i],nn[i]);gross=event.coefficients*(np.exp(lf)+np.exp(lb))
            activity=np.array([np.sum(gross),np.dot(event.omega,gross)])
            del event
            values=[bound(a,g,cell.energies[-1],cell.energies,pp[i],nn[i],order) for order in (16,32)]
            relative=values[-1]/np.maximum(activity,1e-300)
            quadrature=abs(values[-1]-values[-2])/np.maximum(abs(values[-1]),1e-300)
            rows.append(dict(time=float(time_value),cell=i,amplitude=float(a),gamma=g,
                orders=[16,32],bounds=[v.tolist() for v in values],gross_activity=activity.tolist(),
                relative_bounds=relative.tolist(),quadrature_relative_difference=quadrature.tolist()))
    worst=max(max(row['relative_bounds']) for row in rows)
    refinement=max(max(row['quadrature_relative_difference']) for row in rows)
    passed=worst<=1e-3 and refinement<=2e-4
    result=dict(status='PASS' if passed else 'FAIL',scope='Actual equilibrium populations, outgoing absorption into unrepresented higher states only; external holes<=1. No arbitrary incoming external electrons assumed. Finite-eta independent high-energy Usadel solve.',
        worst_relative_upper_bound=worst,maximum_quadrature_relative_difference=refinement,
        maximum_original_spectral_residual=max(spectral_residuals),
        distinct_times_per_cell=len(np.unique(times)),cells=2,rows=rows,
        source_sha256=sha(__file__),trajectory_sha256=sha(path),runtime_seconds=time.perf_counter()-started)
    save(output,result)
    if not passed:raise RuntimeError('Independent equilibrium upper-support check failed')
    return {k:result[k] for k in ('status','scope','worst_relative_upper_bound',
        'maximum_quadrature_relative_difference','maximum_original_spectral_residual','runtime_seconds')}


def execute(plan,path,output):
    check_sources(plan)
    if output.exists():raise ValueError('Existing output preserved; choose a fresh directory')
    output.mkdir(parents=True)
    started=time.perf_counter();report=dict(schema='pysnspd.stage2.special-cases.results.v1',
        status='INCOMPLETE',stage2_admission=False,plan_sha256=sha(path),
        host=os.uname().nodename if hasattr(os,'uname') else os.environ.get('COMPUTERNAME'),
        cases=[],postprocessing={},source_hashes=plan['source_hashes'])
    result_path=output/'special_assessment.json';save(result_path,report)
    try:
        paths=[]
        for task in plan['tasks']:
            target=output/(task['id']+'.json')
            command=[sys.executable,'-u',str(RUNNER)]
            for key,value in task['parameters'].items():
                command += ['--'+key.replace('_','-'),str(value)]
            command += ['--output',str(target)]
            print(json.dumps(dict(event='START_CASE',id=task['id'])),flush=True)
            run_command(command,target.with_suffix('.console.log'))
            result=assess_trajectory(target,task);report['cases'].append(result)
            report['runtime_seconds']=time.perf_counter()-started;save(result_path,report)
            print(json.dumps(dict(event='CASE_ASSESSED',id=task['id'],status=result['status'],runtime_seconds=result['runtime_seconds'])),flush=True)
            if result['status']!='PASS':raise RuntimeError('Case failed registered checks:'+task['id'])
            paths.append(target)
        report['postprocessing']['equilibrium_fields']=equilibrium_fields(paths[0],output/'special_equilibrium_fields.json')
        save(result_path,report)
        report['postprocessing']['equilibrium_support']=equilibrium_support(paths[0],output/'special_equilibrium_support.json')
        save(result_path,report)
        for check in ('fields','support'):
            target=output/f'special_boundary_{check}.json'
            run_command([sys.executable,str(ROOT/f'sandbox/stage2_cells/check_trajectory_{check}.py'),
                         *map(str,paths[1:]),'--output',str(target)],target.with_suffix('.console.log'))
            result=json.loads(target.read_text())
            report['postprocessing']['boundary_'+check]={k:result[k] for k in
                ('status','maxima','worst_relative_upper_bound','runtime_seconds') if k in result}
            if result['status']!='PASS':raise RuntimeError('Boundary postprocessing failed:'+check)
        check_sources(plan)
        report['status']='PASS_SHORT_SPECIAL_SCENARIOS'
        report['remaining']=['These cases are not a three-level temporal refinement or a global stage2 closure.']
    except Exception as exc:
        report['status']='FAILED_OR_INCOMPLETE';report['reason']=repr(exc)
        raise
    finally:
        report['runtime_seconds']=time.perf_counter()-started;save(result_path,report)
        print(json.dumps({k:report[k] for k in ('status','runtime_seconds')}),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True)
    parser.add_argument('--output-root',type=Path)
    modes=parser.add_mutually_exclusive_group(required=True)
    modes.add_argument('--prepare',action='store_true')
    modes.add_argument('--execute',action='store_true')
    args=parser.parse_args();plan_path=args.plan.resolve()
    if args.prepare:
        if plan_path.exists():raise ValueError('Registration already exists; preserve it')
        save(plan_path,registration());print('PREPARED_NO_CALCULATIONS',sha(plan_path));return
    if args.output_root is None:parser.error('--output-root required')
    execute(json.loads(plan_path.read_text()),plan_path,args.output_root.resolve())


if __name__=='__main__':main()
