"""Stage 3A nodal fourth-order static experiments with explicit validation scopes.

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
from pysnspd.experimental.spatial_nodal import PeriodicNodalSpatialFunctional
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
        'pysnspd/experimental/spatial_nodal.py',
        'sandbox/stage3_spatial/progress.py',
        'sandbox/stage3_spatial/run_nodal_batch.py',
        'docs/implementation/stage1_r2/catalogs/occupation_catalog.npz',
    ]
    values={p:digest(ROOT/p) for p in paths}
    values[str(Path(registration).resolve().relative_to(ROOT)).replace('\\','/')]=digest(registration)
    return values


def project(values,coordinates):
    # Explicit Fourier coefficients of the same physical modes on each grid.
    return np.array([np.mean(values*np.exp(-2j*np.pi*k*coordinates)) for k in (-2,-1,0,1,2)])


def uniform_oracle(model,catalog,case,p):
    """Independent D4/L4 helical algebra; never calls the spatial evaluator.

    Q is the first-derivative symbol; K is the fourth-order negative-Laplacian
    symbol. K is not Q**2. Differentiation holds nodal populations fixed.
    """
    a0=case['amplitude_base'];q0=case['q_bare_ell0'];h=model.h_bar
    angle=q0*h;a=a0;s=a*a+.1**2;m=a*a/s
    Q=(8*np.sin(angle)-np.sin(2*angle))/(6*h)
    K=(30-32*np.cos(angle)+2*np.cos(2*angle))/(12*h*h)
    q=m*Q;gamma=q*q/model.gap_ratio
    u,ua,ug=(catalog.evaluate(a,gamma,p) if np.any(p) else catalog.vacuum.evaluate(a,gamma))
    dq_da=2*a*.1**2*Q/s**2
    B=2*q*ug/model.gap_ratio-2*model.kappa*a*a*q
    amplitude_derivative=ua+B*dq_da+2*model.kappa*a*(K-q*q)
    phase_derivative=(B*m*(8*np.cos(angle)-2*np.cos(2*angle))/(6*h)
        +model.kappa*a*a*(32*np.sin(angle)-4*np.sin(2*angle))/(12*h*h))
    m=a0*a0/(a0*a0+.1**2);qc=m*q0
    uc,_,ugc=(catalog.evaluate(a0,qc*qc/model.gap_ratio,p) if np.any(p)
              else catalog.vacuum.evaluate(a0,qc*qc/model.gap_ratio))
    return dict(energy_bar=model.length_bar*(u+model.kappa*a*a*(K-q*q)),
        force_radial_density_bar=amplitude_derivative,link_current_bar=h*phase_derivative,
        nodal_amplitude_bar=a,first_derivative_symbol_Q=Q,negative_laplacian_symbol_K=K,
        q_delta_bar=q,gamma_bar=gamma,
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


def fd_required(cells,registration,*,pilot=False):
    return bool(pilot or cells in registration['fd_cells'])


def derivative_metadata(cells,registration,checks,*,pilot=False):
    if not fd_required(cells,registration,pilot=pilot):
        if checks:raise ValueError('Skipped FD scope cannot contain fabricated checks')
        return dict(status='NOT_REPEATED_UNIT_AND_N16_REFERENCES',performed=False,
            scope='No directional finite differences on this grid; structural checks below are separate. Derivative evidence is delegated to unit tests and registered N16 references, whose completion must be reported independently.',
            registered_fd_cells=list(registration['fd_cells']))
    if set(checks)!={'force','current'}:
        raise ValueError('Executed FD scope requires both force and current checks')
    status=('FAIL_IDENTITY' if any(v['status']=='FAIL_IDENTITY' for v in checks.values()) else
            'INCONCLUSIVE_REFERENCE' if any(not v['passed'] for v in checks.values()) else 'PASS')
    return dict(status=status,performed=True,registered_fd_cells=list(registration['fd_cells']),
                scope='Pilot N8 derivative evidence only' if pilot else 'Registered N16 directional references',
                **checks)


def accepted_case(record):
    metadata=record.get('derivative_checks',{})
    if record['status']=='PASS_STATIC_IDENTITIES_AND_SAMPLED_SYMBOL':
        return metadata.get('performed') is True and metadata.get('status')=='PASS'
    if record['status']=='PASS_STATIC_STRUCTURE_WITHOUT_REPEATED_FD':
        return (metadata.get('performed') is False
                and metadata.get('status')=='NOT_REPEATED_UNIT_AND_N16_REFERENCES')
    return False


def spectral_probe_record(link,coarse,fine,coarse_moments,fine_moments):
    """Report both stiffness branches and moments, not just the minimum.

    The radial eigenvalue may stay exactly 2*kappa while the phase branch
    changes under spectral refinement. Its unchanged minimum is therefore
    insufficient to describe the actual spectral sensitivity.
    """
    coarse_eigen=np.asarray(coarse['eigenvalues'],float)
    fine_eigen=np.asarray(fine.eigenvalues,float)
    coarse_matrix=np.asarray(coarse['matrix'],float)
    fine_matrix=np.asarray(fine.matrix,float)
    coarse_moments=np.asarray(coarse_moments,float)
    fine_moments=np.asarray(fine_moments,float)
    shift=float(abs(fine_eigen[0]-coarse_eigen[0]))
    margin=float(min(fine_eigen[0],coarse_eigen[0])-fine.uncertainty-coarse['uncertainty']-shift)
    matrix_difference=float(np.linalg.norm(fine_matrix-coarse_matrix,2))
    return dict(link=int(link),coarse_min=float(coarse_eigen[0]),fine_min=float(fine_eigen[0]),
        coarse_uncertainty=float(coarse['uncertainty']),fine_uncertainty=float(fine.uncertainty),
        coarse_eigenvalues=coarse_eigen.tolist(),fine_eigenvalues=fine_eigen.tolist(),
        eigenvalue_absolute_differences=abs(fine_eigen-coarse_eigen).tolist(),
        eigenvalue_scaled_differences=(abs(fine_eigen-coarse_eigen)/np.maximum(1.,abs(fine_eigen))).tolist(),
        coarse_matrix=coarse_matrix.tolist(),fine_matrix=fine_matrix.tolist(),
        matrix_operator_norm_difference=matrix_difference,
        matrix_scaled_operator_norm_difference=matrix_difference/max(1.,float(np.linalg.norm(fine_matrix,2))),
        electronic_moments=dict(names=['u','u_a','u_Gamma'],
            coarse=coarse_moments.tolist(),fine=fine_moments.tolist(),
            absolute_differences=abs(fine_moments-coarse_moments).tolist(),
            scaled_differences=(abs(fine_moments-coarse_moments)/np.maximum(1.,abs(fine_moments))).tolist()),
        scale_definition='Absolute difference / max(1, absolute fine value); matrices use spectral operator norm.',
        sampled_spectral_shift=shift,remaining_positive_margin=margin)


def run_case(catalog,case,cells,registration,output,source_hashes,log_callback,*,pilot=False):
    started=time.perf_counter()
    geometry=registration['geometry']
    model=PeriodicNodalSpatialFunctional(catalog,geometry['length_m'],
        geometry['width_m']*geometry['thickness_m'],cells)
    x,z,p,twist=make_state(model,catalog,case)
    direction=(.2+np.cos(2*np.pi*x/model.length_bar)+.7j*np.sin(4*np.pi*x/model.length_bar))*np.exp(1j*case['q_bare_ell0']*x)
    direction/=np.max(abs(direction))
    link_direction=.3+.7*np.cos(2*np.pi*(np.arange(cells)+.5)/cells)
    link_direction/=np.max(abs(link_direction))
    initial_path=output.with_name(output.stem+'_initial.npz')
    np.savez_compressed(initial_path,delta=z,p_nodes=p,cartesian_direction=direction,
        link_direction=link_direction,twist=twist,count_nodes=catalog.count_nodes)
    initial_hash=digest(initial_path)
    perform_fd=fd_required(cells,registration,pilot=pilot)
    # Baseline, optional directional energies, gauge, uniform, N symbols, spectral.
    progress=Progress(4+cells+(8*len(registration['fd_steps']) if perform_fd else 0),
        f"{case['id']} / N={cells}",min_interval=2.,callback=log_callback)
    progress.start_task('energía, fuerza y corriente')
    callback=lambda i,n:progress.checkpoint(f'nodo {i+1}/{n}')
    base=model.evaluate(z,p,twist=twist,require_stability=False,on_cell=callback)
    progress.advance()
    predicted_force=float(np.sum(base.gradient_cartesian_bar[:,0]*direction.real+
                                  base.gradient_cartesian_bar[:,1]*direction.imag))
    predicted_current=float(np.dot(base.link_current_bar,link_direction))
    differences={}
    for kind,prediction in ((('force',predicted_force),('current',predicted_current)) if perform_fd else ()):
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
        tolerance=(registration['thresholds']['derivative_absolute']
                   +registration['thresholds']['derivative_relative']*abs(estimates[-1]))
        error=abs(estimates[-1]-prediction)
        uncertainty=abs(estimates[-1]-estimates[-2])
        differences[kind]=dict(predicted=prediction,fd_estimates=estimates,absolute_error=error,
            tolerance=tolerance,fd_uncertainty=uncertainty,
            status='INCONCLUSIVE_REFERENCE' if uncertainty>registration['thresholds']['reference_budget_fraction']*tolerance else
                   'PASS' if error<=tolerance else 'FAIL_IDENTITY',
            passed=bool(error<=tolerance and uncertainty<=registration['thresholds']['reference_budget_fraction']*tolerance))
    if not perform_fd:
        progress.note('FD no repetidas en esta malla; evidencia derivativa de unidades y N16 se informa por separado')
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
    derivative=model.covariant_derivative(z,twist=twist)
    symbols=[]
    for i in range(cells):
        progress.start_task(f'D.36: nodo {i+1}/{cells}')
        symbol=model.principal_symbol(z[i],derivative[i],p[i])
        symbols.append(dict(matrix=symbol.matrix.tolist(),eigenvalues=symbol.eigenvalues.tolist(),
            uncertainty=symbol.uncertainty,stable=symbol.stable,gamma=symbol.gamma_bar,
            curvatures=list(symbol.gamma_curvatures),fd_steps=list(symbol.fd_steps)))
        progress.advance()
    progress.start_task('alcance espectral / contraste 630-1260')
    spectral=dict(status='NOT_RUN_PILOT_OR_COARSER_GRID',probes=[])
    if cells==32 and case['occupation']!='vacuum':
        finer_catalog=refined_count_catalog(catalog.base,refinement=2)
        finer=PeriodicNodalSpatialFunctional(finer_catalog,geometry['length_m'],
            geometry['width_m']*geometry['thickness_m'],cells)
        _,_,finer_p,_=make_state(finer,finer_catalog,case)
        indices=sorted(set([int(np.argmin([s['eigenvalues'][0] for s in symbols])),
                            int(np.argmax(base.gamma_nodes_bar)),int(np.argmin(base.amplitude_nodes_bar))]))
        for i in indices:
            progress.checkpoint(f'contraste espectral en nodo {i+1}',force=True)
            fine=finer.principal_symbol(z[i],derivative[i],finer_p[i])
            fine_moments=finer.local_density(z[i],derivative[i],finer_p[i]).electronic_derivatives
            spectral['probes'].append(spectral_probe_record(i,symbols[i],fine,
                base.electronic_derivatives[i],fine_moments))
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
        local_force=base.electronic_derivatives[:,1]
        measured=density_gradient.real-local_force
        expected=2*model.kappa*case['amplitude_modulation']*(2*np.pi/model.length_bar)**2*np.cos(2*np.pi*x/model.length_bar)
        analytic=dict(relative_rms=float(np.linalg.norm(measured-expected)/np.linalg.norm(expected)),
                      measured=measured.tolist(),expected=expected.tolist())
    global_phase=dict(energy=phase_error,force_covariance=phase_force_error,
        current_invariance=phase_current_error,directional_energy=phase_derivative)
    structural_ok=(max(global_phase.values())<=registration['thresholds']['global_phase_scaled_error']
        and max(uniform_errors.values())<=registration['thresholds'].get('uniform_oracle_scaled_error',1e-4)
        and all(s['stable'] for s in symbols) and spectral['status']!='UNRESOLVED_SAMPLED_SPECTRAL_SIGN'
        and not any(v['status']=='FAIL_IDENTITY' for v in differences.values()))
    derivative_checks=derivative_metadata(cells,registration,differences,pilot=pilot)
    passed=structural_ok and (not perform_fd or derivative_checks['status']=='PASS')
    result=dict(schema='pysnspd.stage3.nodal-static-case.v1',case=case,cells=cells,
        status=('PASS_STATIC_IDENTITIES_AND_SAMPLED_SYMBOL' if perform_fd else
                'PASS_STATIC_STRUCTURE_WITHOUT_REPEATED_FD') if passed else
               'INCONCLUSIVE_REFERENCE' if structural_ok else 'FAIL_STATIC_CHECK',
        scope='Fixed nodal populations, periodic 1D fourth-order static functional; no time, physical terminals or circuit.',
        declared_spatial_order=registration['declared_spatial_order'],
        source_hashes=source_hashes,runtime_seconds=time.perf_counter()-started,
        initial_sha256=initial_hash,
        energy_bar=base.energy_bar,energy_J=base.energy_J,length_bar=model.length_bar,
        h_bar=model.h_bar,ell0_m=model.ell0_m,energy_scale_J=model.energy_scale_J,
        electron_nodes=len(catalog.count_nodes),derivative_checks=derivative_checks,
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
            current_A=base.current_A.tolist(),gamma=base.gamma_nodes_bar.tolist(),
            energy_density=base.energy_density_bar.tolist()))
    data_path=output.with_suffix('.npz')
    np.savez_compressed(data_path,delta=z,p_nodes=p,count_nodes=catalog.count_nodes,
        count_weights=catalog.count_weights,gradient=base.gradient_cartesian_bar,current_A=base.current_A)
    result['arrays_sha256']=digest(data_path)
    write_json(output,result)
    progress.finish(success=passed)
    return result


def negative_control(catalog):
    model=PeriodicNodalSpatialFunctional(catalog,360e-9,120e-9*7e-9,8)
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
            denominator=registration.get('richardson_denominator_by_observable',{}).get(
                key,registration['richardson_denominator'])
            arrays=[np.atleast_1d(complex_array(r[key]) if isinstance(r[key],dict) else r[key]) for r in values]
            d1=float(np.linalg.norm(arrays[0]-arrays[1]));d2=float(np.linalg.norm(arrays[1]-arrays[2]))
            norm=float(np.linalg.norm(arrays[2]))
            floor=registration['thresholds']['spatial_response_absolute_floor']
            target=registration['thresholds']['nonlinear_mesh_response_target_relative']
            # A vanishing real-field current is an analytic zero, not a
            # relative-convergence result. No other small response receives
            # this exception merely because its measured value is small.
            analytic_zero=(key=='induced_current_modes' and
                case.get('q_bare_ell0')==0. and case.get('phase_modulation')==0.)
            near_zero=norm<=floor
            difference_unresolved=min(d1,d2)<=floor
            degenerate=near_zero or difference_unresolved
            ratio=d1/d2 if d2>0 else None
            estimate=None
            relative_certificate=False
            absolute_control_passed=None
            if analytic_zero:
                absolute_control_passed=bool(max(float(np.linalg.norm(a)) for a in arrays)<=floor)
                passed=absolute_control_passed
                status='PASS_ANALYTIC_ZERO_ABSOLUTE' if passed else 'FAIL_ANALYTIC_ZERO_ABSOLUTE'
            elif degenerate:
                # None deliberately means no relative verdict, rather than
                # promoting rounded differences to a convergence PASS.
                passed=None
                status='ROUND_OFF_LIMITED'
            elif ratio is None or ratio<=1:
                passed=None
                status='INCONCLUSIVE_SPATIAL_ESTIMATE'
            else:
                # Conditional asymptotic estimate, never a rigorous bound.
                estimate=max(d2/denominator,d2/(ratio-1))
                passed=bool(estimate<=target*norm)
                relative_certificate=passed
                status='PASS_RELATIVE_ESTIMATE' if passed else 'FAIL_PRECISION_TARGET'
            metrics[key]=dict(coarse_medium_difference=d1,medium_fine_difference=d2,
                fine_norm=norm,richardson_conservative_estimate=estimate,
                richardson_denominator=denominator,
                relative_estimate=None if estimate is None else estimate/norm,
                observed_order=None if d1<=0 or d2<=0 else float(np.log2(d1/d2)),
                near_zero_response=near_zero,difference_unresolved=difference_unresolved,
                degenerate=degenerate,passed=passed,status=status,
                relative_certificate=relative_certificate,analytic_zero=analytic_zero,
                absolute_control_passed=absolute_control_passed,
                interpretation=('ANALYTIC_ZERO_ABSOLUTE_CONTROL_NO_RELATIVE_CERTIFICATE' if analytic_zero else
                    'ROUND_OFF_LIMITED_NO_RELATIVE_CERTIFICATE' if degenerate else
                    'CONDITIONAL_EXTRAPOLATION_ESTIMATE'))
        if values[-1]['analytic_gradient'] is not None:
            errors=[v['analytic_gradient']['relative_rms'] for v in values]
            gradient_passed=bool(errors[-1]<=registration['thresholds']['gradient_finest_relative_error']
                                 and errors[2]<errors[1]<errors[0])
            metrics['isolated_gradient']=dict(errors=errors,passed=gradient_passed,
                status='PASS_ANALYTIC_GRADIENT' if gradient_passed else 'FAIL_ANALYTIC_GRADIENT')
        failed=any(v['passed'] is False for v in metrics.values())
        unresolved=[key for key,value in metrics.items() if value['passed'] is None]
        status=('PRECISION_TARGET_NOT_MET' if failed else 'SPATIAL_RESPONSE_SCOPE_LIMITED' if unresolved
                else 'PASS_REGISTERED_SPATIAL_RESPONSE')
        rows.append(dict(case=case['id'],status=status,metrics=metrics,
            relative_certificate_observables=[key for key,value in metrics.items() if value.get('relative_certificate')],
            absolute_control_observables=[key for key,value in metrics.items() if value.get('absolute_control_passed')],
            unresolved_observables=unresolved))
    return rows


def cost_weight(case,cells,registration,*,pilot=False,reused=None):
    """Relative work estimate, not elapsed time or a physical calculation."""
    reused_fd=(reused is not None and reused.get('derivative_checks',{}).get('performed') is True)
    evaluations=32 if fd_required(cells,registration,pilot=pilot) or reused_fd else 8
    return cells*evaluations*(.005 if case['occupation']=='vacuum' else 1.)


def verify_reuse(root,task,hashes,expected_case,expected_cells):
    path=Path(root)/(task+'.json')
    if not path.exists():return None
    record=json.loads(path.read_text(encoding='utf-8'))
    if record.get('schema')!='pysnspd.stage3.nodal-static-case.v1':
        raise ValueError('Only this nodal representation can be reused: '+task)
    if record['case']!=expected_case or record['cells']!=expected_cells:
        raise ValueError('Case or mesh mismatch in reused task '+task)
    if record['source_hashes']!=hashes:raise ValueError('Source/input mismatch in reused case '+task)
    if not accepted_case(record):raise ValueError('Cannot reuse failed/incomplete or incorrectly scoped case '+task)
    if digest(path.with_suffix('.npz'))!=record['arrays_sha256']:raise ValueError('Corrupt reused arrays '+task)
    if digest(path.with_name(task+'_initial.npz'))!=record['initial_sha256']:raise ValueError('Corrupt initial state '+task)
    receipt=json.loads(path.with_suffix('.receipt.json').read_text())
    if receipt['json_sha256']!=digest(path):raise ValueError('Corrupt reused result '+task)
    return record


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--registration',default='docs/implementation/stage3/nodal_20260923/registration.json')
    parser.add_argument('--output-root')
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument('--pilot',action='store_true');mode.add_argument('--execute',action='store_true')
    parser.add_argument('--reuse',help='A verified previous pilot output directory')
    args=parser.parse_args()
    registration_path=(ROOT/args.registration).resolve()
    registration=json.loads(registration_path.read_text(encoding='utf-8'))
    if (registration['declared_spatial_order']!=4
            or registration['richardson_denominator']!=15
            or registration['fd_cells']!=[16]):
        raise ValueError('This runner requires the registered fourth-order nodal/FD scope')
    for record in registration['inputs']:
        if digest(ROOT/record['path'])!=record['sha256']:raise ValueError('Registered input changed: '+record['path'])
    if digest(ROOT/registration['catalogue']['source'])!=registration['catalogue']['source_sha256']:
        raise ValueError('Registered catalogue changed')
    tasks=[(case,n) for case in registration['cases'] for n in registration['cell_counts']]
    if args.pilot:tasks=[(c,n) for c,n in tasks if c['id']=='weak_phase_thermal' and n==8]
    if not (args.pilot or args.execute):
        print(json.dumps(dict(mode='DESCRIPTION_ONLY',tasks=[f"{c['id']}_n{n}" for c,n in tasks],
            negative_control=True,output_required=True,declared_spatial_order=4,
            fd_cells=registration['fd_cells'],pilot_runs_directional_fd=True,
            progress='Flushed bars, elapsed time, approximate remaining time and per-node checkpoints.'),indent=2));return
    if not args.output_root:parser.error('--output-root is required for execution')
    out=Path(args.output_root).resolve()
    if out.exists():raise SystemExit('Output directory already exists. Choose a fresh one; use --reuse for validated prior work.')
    out.mkdir(parents=True)
    sources=frozen_sources(registration_path)
    write_json(out/'manifest.json',dict(mode='PILOT' if args.pilot else 'FULL_STATIC_CAMPAIGN',
        source_hashes=sources,host=platform.node(),python=platform.python_version(),numpy=np.__version__,
        declared_spatial_order=4,fd_cells=registration['fd_cells'],pilot_runs_directional_fd=args.pilot,
        tasks=[f"{c['id']}_n{n}" for c,n in tasks]))
    reused={}
    for case,cells in tasks:
        task=f"{case['id']}_n{cells}"
        value=verify_reuse(args.reuse,task,sources,case,cells) if args.reuse else None
        if value is not None and fd_required(cells,registration,pilot=args.pilot) and not value['derivative_checks']['performed']:
            raise ValueError('Reused case lacks required directional FD evidence: '+task)
        if value is not None:reused[task]=value
    # Relative cost weights distinguish vacuum algebra from occupied spectral work.
    # Seed only verified prior results; their historical runtimes do not count as
    # freshly executed work when the terminal estimate updates.
    def weight(case,cells):
        return cost_weight(case,cells,registration,pilot=args.pilot,
                           reused=reused.get(f"{case['id']}_n{cells}"))
    reuse_weight=sum(weight(c,n) for c,n in tasks if f"{c['id']}_n{n}" in reused)
    pilot_seconds=next((v['runtime_seconds'] for v in reused.values()
                        if v['case']['id']=='weak_phase_thermal' and v['cells']==8
                        and v['derivative_checks']['performed']),34.)
    with (out/'progress.jsonl').open('w',encoding='utf-8',buffering=1) as log:
        def log_progress(row): log.write(json.dumps(row,allow_nan=False)+'\n');log.flush()
        batch=Progress(len(tasks),'Etapa 3A nodal D4',min_interval=2.,callback=log_progress,
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
                    estimate=pilot_seconds*weight(case,cells)/(8*32)
                    batch.start_task(task,estimated_seconds=estimate,weight=weight(case,cells))
                    def case_progress(row):
                        log_progress(row)
                        batch.checkpoint('estimación del lote; detalle del caso en la línea anterior')
                    result=run_case(catalog,case,cells,registration,out/(task+'.json'),sources,case_progress,pilot=args.pilot)
                    write_json(out/(task+'.receipt.json'),dict(json_sha256=digest(out/(task+'.json')),
                        arrays_sha256=digest(out/(task+'.npz')),initial_sha256=result['initial_sha256'],source_hashes=sources))
                    batch.advance()
                    if case['id']=='weak_phase_thermal' and cells==8 and result['derivative_checks']['performed']:
                        pilot_seconds=result['runtime_seconds']
                results.append(result)
                if result['status']=='FAIL_STATIC_CHECK':
                    raise RuntimeError('Static identity or stability failure in '+task+'; result preserved')
            grids=assess_grids(results,registration,require_complete=not args.pilot)
            passed=(all(r['status']!='PRECISION_TARGET_NOT_MET' for r in grids) and
                    all(accepted_case(r) for r in results))
            scope_limited=any(r['status']=='SPATIAL_RESPONSE_SCOPE_LIMITED' for r in grids)
            summary=dict(status=('PILOT_PASS_FULL_CAMPAIGN_PENDING' if args.pilot and passed else
                'PILOT_REFERENCE_INCONCLUSIVE' if args.pilot else
                'STATIC_CAMPAIGN_COMPLETE_SCOPE_LIMITED' if passed and scope_limited else
                'STATIC_CAMPAIGN_COMPLETE' if passed else 'STATIC_CAMPAIGN_PRECISION_TARGET_NOT_MET'),
                runtime_seconds=time.perf_counter()-started,completed_cases=len(results),
                declared_spatial_order=4,registered_fd_cells=registration['fd_cells'],
                cases_with_directional_fd=[f"{r['case']['id']}_n{r['cells']}" for r in results if r['derivative_checks']['performed']],
                cases_without_repeated_fd=[f"{r['case']['id']}_n{r['cells']}" for r in results if not r['derivative_checks']['performed']],
                negative_control=negative,grid_comparisons=grids,
                spatial_precision_fully_resolved=bool(not args.pilot and passed and not scope_limited),
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
