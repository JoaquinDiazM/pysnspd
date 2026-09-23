"""Stage 3A registered static experiments with checkpoint-driven progress.

Default: describe the work, without calculating. --pilot runs only the registered
eight-cell timing/identity pilot. --execute runs the full registered campaign.
Existing outputs are never overwritten. --reuse may import an exactly compatible
completed pilot, including its source, input and output hashes.
"""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import shutil
import sys
import time
import traceback
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from pysnspd.experimental.energy_catalog import OccupationEnergyCatalog
from pysnspd.experimental.refined_cells import refined_count_catalog
from pysnspd.experimental.spatial_functional import PeriodicSpatialFunctional
from progress import Progress


def digest(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path,value):
    Path(path).write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')


def complex_record(values):
    values=np.asarray(values)
    return dict(real=values.real.tolist(),imag=values.imag.tolist())


def complex_array(value): return np.array(value['real'])+1j*np.array(value['imag'])


def frozen_sources(registration):
    paths=[
        'pysnspd/experimental/energy_catalog.py',
        'pysnspd/experimental/refined_cells.py',
        'pysnspd/experimental/spatial_functional.py',
        'sandbox/stage3_spatial/progress.py',
        'sandbox/stage3_spatial/run_static_batch.py',
        'docs/implementation/stage1_r2/catalogs/occupation_catalog.npz',
    ]
    values={p:digest(ROOT/p) for p in paths}
    values[str(Path(registration).resolve().relative_to(ROOT)).replace('\\','/')]=digest(registration)
    return values


def project(values,coordinates):
    # Explicit Fourier coefficients of the same physical modes on each grid.
    return np.array([np.mean(values*np.exp(-2j*np.pi*k*coordinates)) for k in (-2,-1,0,1,2)])


def uniform_oracle(model,catalog,case,p):
    """Independent helical link algebra; never calls the spatial evaluator."""
    a0=case['amplitude_base'];q0=case['q_bare_ell0'];h=model.h_bar
    angle=q0*h/2;a=a0*np.cos(angle);d=2*a0*np.sin(angle)/h
    s=a*a+.1**2;q=a*d/s;gamma=q*q/model.gap_ratio
    u,ua,ug=(catalog.evaluate(a,gamma,p) if np.any(p) else catalog.vacuum.evaluate(a,gamma))
    def directional(ad,dd):
        qd=((ad*d+a*dd)*s-a*d*2*a*ad)/s**2
        return ua*ad+ug*2*q*qd/model.gap_ratio+model.kappa*(2*d*dd-2*a*ad*q*q-2*a*a*q*qd)
    phase_derivative=directional(-a0*np.sin(angle)/2,a0*np.cos(angle)/h)
    amplitude_derivative=directional(np.cos(angle),2*np.sin(angle)/h)
    m=a0*a0/(a0*a0+.1**2);qc=m*q0
    uc,_,ugc=(catalog.evaluate(a0,qc*qc/model.gap_ratio,p) if np.any(p)
              else catalog.vacuum.evaluate(a0,qc*qc/model.gap_ratio))
    return dict(energy_bar=model.length_bar*(u+model.kappa*(d*d-a*a*q*q)),
        force_radial_density_bar=amplitude_derivative,link_current_bar=h*phase_derivative,
        continuum_energy_density=uc+model.kappa*a0*a0*(q0*q0-qc*qc),
        continuum_current_bar=2*m*qc*ugc/model.gap_ratio+2*model.kappa*a0*a0*(1-m*m)*q0)


def make_state(model,catalog,case):
    x=np.arange(model.cells)*model.h_bar
    normalized=x/model.length_bar
    q=case['q_bare_ell0']
    amplitude=case['amplitude_base']+case['amplitude_modulation']*np.cos(2*np.pi*normalized)
    phase=q*x+case['phase_modulation']*np.sin(2*np.pi*normalized)
    z=amplitude*np.exp(1j*phase)
    if case['occupation']=='vacuum':
        p=np.zeros_like(catalog.count_nodes)
    else:
        # One fixed preparation, independent of spatial mesh and case fields.
        ar=.9; qr=.1*ar*ar/(ar*ar+.1**2)
        energy=catalog.energy_kernel(ar,qr*qr/model.gap_ratio)[0]
        p=1/(np.exp(np.minimum(energy/.12,700.))+1)
        if case['occupation']=='nonthermal':
            p=p+.02*np.exp(-((catalog.count_nodes-.035)/.01)**2)
    assert np.all((p>=0)&(p<=1))
    return x,z,np.repeat(p[None,:],model.cells,axis=0),q*model.length_bar


def run_case(catalog,case,cells,registration,output,source_hashes,log_callback):
    started=time.perf_counter()
    geometry=registration['geometry']
    model=PeriodicSpatialFunctional(catalog,geometry['length_m'],
        geometry['width_m']*geometry['thickness_m'],cells)
    x,z,p,twist=make_state(model,catalog,case)
    direction=(.2+np.cos(2*np.pi*x/model.length_bar)+.7j*np.sin(4*np.pi*x/model.length_bar))*np.exp(1j*case['q_bare_ell0']*x)
    direction/=np.max(abs(direction))
    link_direction=.3+.7*np.cos(2*np.pi*(np.arange(cells)+.5)/cells)
    link_direction/=np.max(abs(link_direction))
    initial_path=output.with_name(output.stem+'_initial.npz')
    np.savez_compressed(initial_path,delta=z,p_links=p,cartesian_direction=direction,
        link_direction=link_direction,twist=twist,count_nodes=catalog.count_nodes)
    initial_hash=digest(initial_path)
    # Baseline, 24 energies, gauge, uniform, N symbols and one scoped refinement.
    progress=Progress(28+cells,f"{case['id']} / N={cells}",min_interval=2.,callback=log_callback)
    progress.start_task('energía, fuerza y corriente')
    callback=lambda i,n:progress.checkpoint(f'enlace {i+1}/{n}')
    base=model.evaluate(z,p,twist=twist,require_stability=False,on_cell=callback)
    progress.advance()
    predicted_force=float(np.sum(base.gradient_cartesian_bar[:,0]*direction.real+
                                  base.gradient_cartesian_bar[:,1]*direction.imag))
    predicted_current=float(np.dot(base.link_current_bar,link_direction))
    differences={}
    for kind,prediction in (('force',predicted_force),('current',predicted_current)):
        estimates=[]
        for step in registration['fd_steps']:
            energies=[]
            for multiplier in (-2,-1,1,2):
                progress.start_task(f'{kind}: h={step:g}, desplazamiento={multiplier:+d}')
                value=model.evaluate(z+multiplier*step*direction if kind=='force' else z,p,
                    twist=twist,link_phases=multiplier*step*link_direction if kind=='current' else None,
                    require_stability=False,on_cell=callback).energy_bar
                energies.append(value);progress.advance()
            estimates.append(float((energies[0]-8*energies[1]+8*energies[2]-energies[3])/(12*step)))
        tolerance=1e-8+1e-4*abs(estimates[-1])
        error=abs(estimates[-1]-prediction)
        uncertainty=abs(estimates[-1]-estimates[-2])
        differences[kind]=dict(predicted=prediction,fd_estimates=estimates,absolute_error=error,
            tolerance=tolerance,fd_uncertainty=uncertainty,
            status='INCONCLUSIVE_REFERENCE' if uncertainty>.2*tolerance else
                   'PASS' if error<=tolerance else 'FAIL_IDENTITY',
            passed=bool(error<=tolerance and uncertainty<=.2*tolerance))
    progress.start_task('fase global constante')
    rotated=model.evaluate(z*np.exp(.371j),p,twist=twist,require_stability=False,on_cell=callback)
    phase_error=abs(rotated.energy_bar-base.energy_bar)/max(model.length_bar,abs(base.energy_bar))
    gradient_complex=base.gradient_cartesian_bar[:,0]+1j*base.gradient_cartesian_bar[:,1]
    rotated_gradient=rotated.gradient_cartesian_bar[:,0]+1j*rotated.gradient_cartesian_bar[:,1]
    phase_force_error=float(np.max(abs(rotated_gradient-gradient_complex*np.exp(.371j)))/max(1.,np.max(abs(gradient_complex))))
    phase_current_error=float(np.max(abs(rotated.link_current_bar-base.link_current_bar))/max(1.,np.max(abs(base.link_current_bar))))
    phase_derivative=abs(float(np.real(np.vdot(gradient_complex,1j*z))))/max(1.,model.length_bar)
    progress.advance()
    progress.start_task('estado uniforme de referencia a la misma torsión')
    uniform_z=case['amplitude_base']*np.exp(1j*case['q_bare_ell0']*x)
    uniform=model.evaluate(uniform_z,p,twist=twist,require_stability=False,on_cell=callback)
    uniform_reference=uniform_oracle(model,catalog,case,p[0])
    ug_complex=uniform.gradient_cartesian_bar[:,0]+1j*uniform.gradient_cartesian_bar[:,1]
    urad=ug_complex*np.exp(-1j*case['q_bare_ell0']*x)/model.h_bar
    uniform_errors=dict(
        energy=abs(uniform.energy_bar-uniform_reference['energy_bar'])/max(model.length_bar,abs(uniform_reference['energy_bar'])),
        force=float(np.max(abs(urad-uniform_reference['force_radial_density_bar'])))/max(1.,abs(uniform_reference['force_radial_density_bar'])),
        current=float(np.max(abs(uniform.link_current_bar-uniform_reference['link_current_bar'])))/max(1.,abs(uniform_reference['link_current_bar'])))
    progress.advance()
    right=np.roll(z,-1).copy();right[-1]*=np.exp(1j*twist)
    midpoint=(z+right)/2;derivative=(right-z)/model.h_bar
    symbols=[]
    for i in range(cells):
        progress.start_task(f'D.36: enlace {i+1}/{cells}')
        symbol=model.principal_symbol(midpoint[i],derivative[i],p[i])
        symbols.append(dict(eigenvalues=symbol.eigenvalues.tolist(),uncertainty=symbol.uncertainty,
            stable=symbol.stable,gamma=symbol.gamma_bar,curvatures=list(symbol.gamma_curvatures)))
        progress.advance()
    progress.start_task('alcance espectral / contraste 630-1260')
    spectral=dict(status='NOT_RUN_PILOT_OR_COARSER_GRID',probes=[])
    if cells==32 and case['occupation']!='vacuum':
        finer_catalog=refined_count_catalog(catalog.base,refinement=2)
        finer=PeriodicSpatialFunctional(finer_catalog,geometry['length_m'],
            geometry['width_m']*geometry['thickness_m'],cells)
        _,_,finer_p,_=make_state(finer,finer_catalog,case)
        indices=sorted(set([int(np.argmin([s['eigenvalues'][0] for s in symbols])),
                            int(np.argmax(base.gamma_links_bar)),int(np.argmin(base.amplitude_links_bar))]))
        for i in indices:
            progress.checkpoint(f'contraste espectral en enlace {i+1}',force=True)
            fine=finer.principal_symbol(midpoint[i],derivative[i],finer_p[i])
            coarse=symbols[i]
            shift=abs(fine.eigenvalues[0]-coarse['eigenvalues'][0])
            margin=min(fine.eigenvalues[0],coarse['eigenvalues'][0])-fine.uncertainty-coarse['uncertainty']-shift
            spectral['probes'].append(dict(link=i,coarse_min=coarse['eigenvalues'][0],
                fine_min=float(fine.eigenvalues[0]),fine_uncertainty=fine.uncertainty,
                sampled_spectral_shift=float(shift),remaining_positive_margin=float(margin)))
        spectral['status']='PASS_SAMPLED_SPECTRAL_SIGN' if all(v['remaining_positive_margin']>0 for v in spectral['probes']) else 'UNRESOLVED_SAMPLED_SPECTRAL_SIGN'
    elif case['occupation']=='vacuum':
        spectral['status']='VACUUM_NO_OCCUPIED_SPECTRAL_QUADRATURE'
    progress.advance()
    density_gradient=(base.gradient_cartesian_bar[:,0]+1j*base.gradient_cartesian_bar[:,1])/model.h_bar
    reference_gradient=(uniform.gradient_cartesian_bar[:,0]+1j*uniform.gradient_cartesian_bar[:,1])/model.h_bar
    induced=(density_gradient-reference_gradient)*np.exp(-1j*case['q_bare_ell0']*x)
    current_induced=base.link_current_bar-uniform.link_current_bar
    analytic=None
    if case['id']=='weak_amplitude_vacuum':
        local_force=.5*(base.electronic_derivatives[:,1]+np.roll(base.electronic_derivatives[:,1],1))
        measured=density_gradient.real-local_force
        expected=2*model.kappa*case['amplitude_modulation']*(2*np.pi/model.length_bar)**2*np.cos(2*np.pi*x/model.length_bar)
        analytic=dict(relative_rms=float(np.linalg.norm(measured-expected)/np.linalg.norm(expected)),
                      measured=measured.tolist(),expected=expected.tolist())
    global_phase=dict(energy=phase_error,force_covariance=phase_force_error,
        current_invariance=phase_current_error,directional_energy=phase_derivative)
    structural_ok=(max(global_phase.values())<=1e-10 and max(uniform_errors.values())<=1e-4
        and all(s['stable'] for s in symbols) and spectral['status']!='UNRESOLVED_SAMPLED_SPECTRAL_SIGN'
        and not any(v['status']=='FAIL_IDENTITY' for v in differences.values()))
    passed=structural_ok and all(v['passed'] for v in differences.values())
    result=dict(schema='pysnspd.stage3.static-case.v1',case=case,cells=cells,
        status='PASS_STATIC_IDENTITIES_AND_SAMPLED_SYMBOL' if passed else
               'INCONCLUSIVE_REFERENCE' if structural_ok else 'FAIL_STATIC_CHECK',
        scope='Fixed populations, periodic 1D static functional; no time, physical terminals or circuit.',
        source_hashes=source_hashes,runtime_seconds=time.perf_counter()-started,
        initial_sha256=initial_hash,
        energy_bar=base.energy_bar,energy_J=base.energy_J,length_bar=model.length_bar,
        h_bar=model.h_bar,ell0_m=model.ell0_m,energy_scale_J=model.energy_scale_J,
        electron_nodes=len(catalog.count_nodes),derivative_checks=differences,
        global_phase_error=phase_error,global_phase_checks=global_phase,
        uniform_reference=uniform_reference,uniform_reference_errors=uniform_errors,
        principal_symbols=symbols,spectral_sign_comparison=spectral,
        symbol_scope='Longitudinal sign on the candidate spectral quadrature; FD uncertainty excludes spectral/regulator error.',
        excess_energy_density_bar=(base.energy_bar-uniform.energy_bar)/model.length_bar,
        induced_force_modes=complex_record(project(induced,x/model.length_bar)),
        induced_current_modes=complex_record(project(current_induced,(np.arange(cells)+.5)/cells)),
        analytic_gradient=analytic,
        observable_arrays=dict(x_over_ell0=x.tolist(),amplitude=abs(z).tolist(),
            phase=np.unwrap(np.angle(z)).tolist(),force_density=complex_record(density_gradient),
            current_A=base.current_A.tolist(),gamma=base.gamma_links_bar.tolist(),
            energy_density=base.energy_density_bar.tolist()))
    data_path=output.with_suffix('.npz')
    np.savez_compressed(data_path,delta=z,p_links=p,count_nodes=catalog.count_nodes,
        count_weights=catalog.count_weights,gradient=base.gradient_cartesian_bar,current_A=base.current_A)
    result['arrays_sha256']=digest(data_path)
    write_json(output,result)
    progress.finish(success=passed)
    return result


def negative_control(catalog):
    model=PeriodicSpatialFunctional(catalog,360e-9,120e-9*7e-9,8)
    amplitude=.6;gradient=1j*amplitude;rho=amplitude**2;s=rho+.1**2
    q=rho/s;gamma=q*q/model.gap_ratio
    symbol=model.principal_symbol(amplitude,gradient,np.zeros_like(catalog.count_nodes))
    fg=np.pi*amplitude/2-2*gamma/3
    fgg=-2/3
    tangent=2*model.kappa+rho/s**2*(2*fg/model.gap_ratio+
        4*q*q*fgg/model.gap_ratio**2-2*model.kappa*rho)
    independent=np.sort([2*model.kappa,tangent])
    error=float(np.max(abs(symbol.eigenvalues-independent)))
    return dict(status='PASS_EXPECTED_REJECTION' if not symbol.stable and tangent<0 and error<1e-8 else 'FAIL_NEGATIVE_CONTROL',
        amplitude_over_Delta0=.6,bare_q_times_ell0=1.,q_delta=q,
        measured_eigenvalues=symbol.eigenvalues.tolist(),reference_eigenvalues=independent.tolist(),
        error=error,uncertainty=symbol.uncertainty,admitted_for_evolution=False)


def assess_grids(results,registration,require_complete=False):
    expected={(c['id'],n) for c in registration['cases'] for n in registration['cell_counts']}
    keys=[(r['case']['id'],r['cells']) for r in results]
    if len(set(keys))!=len(keys):raise ValueError('Duplicate cases in spatial assessment')
    if not set(keys)<=expected:raise ValueError('Unexpected cases in spatial assessment')
    if require_complete and set(keys)!=expected:raise ValueError('Incomplete full spatial campaign')
    rows=[]
    for case in registration['cases']:
        values=sorted([r for r in results if r['case']['id']==case['id']],key=lambda r:r['cells'])
        if [r['cells'] for r in values]!=registration['cell_counts']:continue
        if case['id'].startswith('uniform'):
            rows.append(dict(case=case['id'],status='STATIC_IDENTITIES_ONLY',
                reason='Uniform induced responses vanish; no relative zero-over-zero grid admission.',
                energy_density=[r['energy_bar']/r['length_bar'] for r in values],
                mean_current_A=[float(np.mean(r['observable_arrays']['current_A'])) for r in values],
                continuum_reference=values[-1]['uniform_reference']))
            continue
        metrics={}
        for key in ('excess_energy_density_bar','induced_force_modes','induced_current_modes'):
            arrays=[np.atleast_1d(complex_array(r[key]) if isinstance(r[key],dict) else r[key]) for r in values]
            d1=float(np.linalg.norm(arrays[0]-arrays[1]));d2=float(np.linalg.norm(arrays[1]-arrays[2]))
            norm=float(np.linalg.norm(arrays[2]));floor=1e-8
            degenerate=norm<=floor
            ratio=d1/d2 if d2>0 else None
            estimate=(max(d2/3,d2/(ratio-1)) if ratio is not None and ratio>1 else
                      0. if d2==0 else None)
            # Observed-order extrapolation is a diagnostic estimate, not a rigorous bound.
            passed=(max(d1,d2)<=floor if degenerate else estimate is not None and d2<d1 and estimate<=.01*norm)
            metrics[key]=dict(coarse_medium_difference=d1,medium_fine_difference=d2,
                fine_norm=norm,richardson_conservative_estimate=estimate,
                relative_estimate=None if degenerate or estimate is None else estimate/norm,
                observed_order=None if d1<=0 or d2<=0 else float(np.log2(d1/d2)),
                degenerate=degenerate,passed=bool(passed),
                interpretation='ROUND_OFF_LIMITED_NO_RELATIVE_CERTIFICATE' if degenerate else 'CONDITIONAL_EXTRAPOLATION_ESTIMATE')
        if values[-1]['analytic_gradient'] is not None:
            errors=[v['analytic_gradient']['relative_rms'] for v in values]
            metrics['isolated_gradient']=dict(errors=errors,passed=bool(errors[-1]<=.01 and errors[2]<errors[1]<errors[0]))
        passed=all(v['passed'] for v in metrics.values())
        rows.append(dict(case=case['id'],status='PASS_REGISTERED_SPATIAL_RESPONSE' if passed else 'PRECISION_TARGET_NOT_MET',metrics=metrics))
    return rows


def verify_reuse(root,task,hashes,expected_case,expected_cells):
    path=Path(root)/(task+'.json')
    if not path.exists():return None
    record=json.loads(path.read_text(encoding='utf-8'))
    if record['case']!=expected_case or record['cells']!=expected_cells:
        raise ValueError('Case or mesh mismatch in reused task '+task)
    if record['source_hashes']!=hashes:raise ValueError('Source/input mismatch in reused case '+task)
    if record['status']!='PASS_STATIC_IDENTITIES_AND_SAMPLED_SYMBOL':raise ValueError('Cannot reuse failed/incomplete case '+task)
    if digest(path.with_suffix('.npz'))!=record['arrays_sha256']:raise ValueError('Corrupt reused arrays '+task)
    if digest(path.with_name(task+'_initial.npz'))!=record['initial_sha256']:raise ValueError('Corrupt initial state '+task)
    receipt=json.loads(path.with_suffix('.receipt.json').read_text())
    if receipt['json_sha256']!=digest(path):raise ValueError('Corrupt reused result '+task)
    return record


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--registration',default='docs/implementation/stage3/iteration_20260923/registration.json')
    parser.add_argument('--output-root')
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument('--pilot',action='store_true');mode.add_argument('--execute',action='store_true')
    parser.add_argument('--reuse',help='A verified previous pilot output directory')
    args=parser.parse_args()
    registration_path=(ROOT/args.registration).resolve()
    registration=json.loads(registration_path.read_text(encoding='utf-8'))
    for record in registration['inputs']:
        if digest(ROOT/record['path'])!=record['sha256']:raise ValueError('Registered input changed: '+record['path'])
    if digest(ROOT/registration['catalogue']['source'])!=registration['catalogue']['source_sha256']:
        raise ValueError('Registered catalogue changed')
    tasks=[(case,n) for case in registration['cases'] for n in registration['cell_counts']]
    if args.pilot:tasks=[(c,n) for c,n in tasks if c['id']=='weak_phase_thermal' and n==8]
    if not (args.pilot or args.execute):
        print(json.dumps(dict(mode='DESCRIPTION_ONLY',tasks=[f"{c['id']}_n{n}" for c,n in tasks],
            negative_control=True,output_required=True,progress='Flushed bars, elapsed time, approximate remaining time and per-link checkpoints.'),indent=2));return
    if not args.output_root:parser.error('--output-root is required for execution')
    out=Path(args.output_root).resolve()
    if out.exists():raise SystemExit('Output directory already exists. Choose a fresh one; use --reuse for validated prior work.')
    out.mkdir(parents=True)
    sources=frozen_sources(registration_path)
    write_json(out/'manifest.json',dict(mode='PILOT' if args.pilot else 'FULL_STATIC_CAMPAIGN',
        source_hashes=sources,host=platform.node(),python=platform.python_version(),numpy=np.__version__,
        tasks=[f"{c['id']}_n{n}" for c,n in tasks]))
    reused={}
    for case,cells in tasks:
        task=f"{case['id']}_n{cells}"
        value=verify_reuse(args.reuse,task,sources,case,cells) if args.reuse else None
        if value is not None:reused[task]=value
    # Relative cost weights distinguish vacuum algebra from occupied spectral work.
    # Seed only verified prior results; their historical runtimes do not count as
    # freshly executed work when the terminal estimate updates.
    def weight(case,cells):return cells*(.005 if case['occupation']=='vacuum' else 1.)
    reuse_weight=sum(weight(c,n) for c,n in tasks if f"{c['id']}_n{n}" in reused)
    pilot_seconds=next((v['runtime_seconds'] for v in reused.values()
                        if v['case']['id']=='weak_phase_thermal' and v['cells']==8),None)
    with (out/'progress.jsonl').open('w',encoding='utf-8',buffering=1) as log:
        def log_progress(row): log.write(json.dumps(row,allow_nan=False)+'\n');log.flush()
        batch=Progress(len(tasks),'Etapa 3A',min_interval=2.,callback=log_progress,
            total_weight=sum(weight(c,n) for c,n in tasks),completed=len(reused),completed_weight=reuse_weight)
        started=time.perf_counter();results=[]
        try:
            catalog=refined_count_catalog(OccupationEnergyCatalog.load(ROOT/'docs/implementation/stage1_r2/catalogs/occupation_catalog.npz'))
            batch.note('Control negativo D.36; no avance de casos hasta completar sus comprobaciones')
            negative=negative_control(catalog);write_json(out/'negative_control.json',negative)
            if negative['status']!='PASS_EXPECTED_REJECTION':raise RuntimeError('D.36 negative control failed')
            for case,cells in tasks:
                task=f"{case['id']}_n{cells}"
                if task in reused:
                    for suffix in ('.json','.npz','.receipt.json','_initial.npz'):
                        shutil.copy2(Path(args.reuse)/(task+suffix),out/(task+suffix))
                    result=reused[task];batch.note('Reutilización verificada; no se recalcula '+task)
                else:
                    estimate=None if pilot_seconds is None else pilot_seconds*weight(case,cells)/8
                    batch.start_task(task,estimated_seconds=estimate,weight=weight(case,cells))
                    def case_progress(row):
                        log_progress(row)
                        batch.checkpoint('estimación del lote; detalle del caso en la línea anterior')
                    result=run_case(catalog,case,cells,registration,out/(task+'.json'),sources,case_progress)
                    write_json(out/(task+'.receipt.json'),dict(json_sha256=digest(out/(task+'.json')),
                        arrays_sha256=digest(out/(task+'.npz')),initial_sha256=result['initial_sha256'],source_hashes=sources))
                    batch.advance()
                results.append(result)
                if result['status']=='FAIL_STATIC_CHECK':
                    raise RuntimeError('Static identity or stability failure in '+task+'; result preserved')
            grids=assess_grids(results,registration,require_complete=not args.pilot)
            passed=(all(r['status']!='PRECISION_TARGET_NOT_MET' for r in grids) and
                    all(r['status']=='PASS_STATIC_IDENTITIES_AND_SAMPLED_SYMBOL' for r in results))
            summary=dict(status=('PILOT_PASS_FULL_CAMPAIGN_PENDING' if args.pilot and passed else
                'PILOT_REFERENCE_INCONCLUSIVE' if args.pilot else
                'STATIC_CAMPAIGN_COMPLETE' if passed else 'STATIC_CAMPAIGN_PRECISION_TARGET_NOT_MET'),
                runtime_seconds=time.perf_counter()-started,completed_cases=len(results),
                negative_control=negative,grid_comparisons=grids,
                stage3_closed=False,production_promotion=False,new_trajectories=0,
                circuit_and_physical_boundaries_implemented=False,source_hashes=sources)
            write_json(out/'summary.json',summary);batch.finish(success=passed)
            print(json.dumps({k:summary[k] for k in ('status','runtime_seconds','completed_cases','stage3_closed')}),flush=True)
            if not passed:raise SystemExit(2)
        except BaseException as exc:
            if not isinstance(exc,SystemExit):
                write_json(out/'failure.json',dict(exception=type(exc).__name__,reason=str(exc),
                    runtime_seconds=time.perf_counter()-started,traceback=traceback.format_exc()))
                batch.finish(success=False)
            raise


if __name__=='__main__':main()
