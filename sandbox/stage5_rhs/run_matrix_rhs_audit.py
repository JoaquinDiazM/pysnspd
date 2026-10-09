"""User-run first-gradient two-mode audit on an admitted actual DC dual graph.

Snapshot directions are declared actuators, not photons or time trajectories.
No neutrality, heat source or material self-energy is fabricated to close the
observed energy-rate defect. Spectral queries alone may take >5 minutes.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import replace
import hashlib
import json
import multiprocessing as mp
import os
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from sandbox.stage4_core.parallel_runtime import limit_thread_environment, linux_resources, resource_budget, initialize_affinity
limit_thread_environment()
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.constants import Boltzmann, elementary_charge, hbar
from scipy.sparse.linalg import splu
from pysnspd.experimental import thermal_spatial_usadel as thermal
from pysnspd.experimental import retarded_spatial_usadel as retarded
from pysnspd.experimental import adiabatic_matrix_usadel as ad
from pysnspd.experimental import frozen_kinetic_usadel as frozen
from pysnspd.experimental.matrix_kinetic_rhs import solve_distribution_rhs


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path,value):
    with Path(path).open('x',encoding='utf8') as stream:
        json.dump(value,stream,indent=2,allow_nan=False);stream.write('\n')


def load_inputs(plan):
    ref=Path(plan['reference']);observer=Path(plan['thermal_observer'])
    summary=json.loads((ref.parent/'summary.json').read_text())
    observed=json.loads((observer/'summary.json').read_text())
    manifest=json.loads((observer/'manifest.json').read_text())
    if not summary.get('admitted_for_dark_hold') or sha(ref)!=summary['reference_sha256']:
        raise ValueError('Unchanged admitted DC reference required')
    if (manifest.get('status')!='COMPLETED' or manifest['reference_sha256']!=sha(ref)
            or observed['count']<summary['plan']['matsubara_count']):
        raise ValueError('Matching completed finite-action observer with at least DC mode count required')
    with np.load(ref) as z:
        graph=thermal.ThermalGraph(z['area_weights'],z['edges'],z['conductance'],
            z['coordinates_bar'],z['boundary_nodes'])
        base=z['delta_bar'].copy();alpha=z['alpha'].copy();ep=z['epsilon_bar'][0]
    with np.load(observer/'fields.npz') as z:
        Eu=Boltzmann*summary['plan']['Tc_K'];gap=z['delta_J']/Eu
        direction=z['audit_direction_gap_J_per_s']/Eu
        Ugradient=z['internal_gradient_J_per_gapbar'].copy()
        Fgradient=z['free_gradient_J_per_gapbar'].copy()
        profile=z['probe_profile'].copy()
        area=z['area_m2'].copy()
    with np.load(observer/'spectral_receipts/mode_00000.npz') as z:
        u0=z['u'].copy()
        if abs(z['epsilon_bar']-ep)>1e-12:
            raise ValueError('Lowest Matsubara anchor changed')
    if any(a.shape!=(graph.n_nodes,) for a in (gap,direction,Ugradient,Fgradient,profile,area)):
        raise ValueError('Observer fields do not match admitted graph')
    if observed['completed_modes']!=observed['count'] or manifest['count']!=observed['count']:
        raise ValueError('Incomplete observer spectral count')
    if observed['depletion']!=manifest['depletion']:
        raise ValueError('Observer metadata disagrees on the probed state')
    ell0_sq=hbar*summary['plan']['diffusion_m2_s']/(2*Eu)
    x=graph.coordinates_bar[:,0]*np.sqrt(ell0_sq)
    offset=(x-.5*(x.min()+x.max()))/(manifest['probe_span_nm']*1e-9)
    expected_profile=np.where(abs(offset)<.5,.5+.5*np.cos(2*np.pi*offset),0.)
    expected_profile[graph.boundary_nodes]=0.
    expected_gap=base*(1-manifest['depletion']*expected_profile)
    expected_direction=-manifest['direction_rate_per_s']*base*expected_profile
    if (not np.allclose(profile,expected_profile,rtol=1e-12,atol=1e-14)
            or not np.allclose(gap,expected_gap,rtol=1e-12,atol=1e-14)
            or not np.allclose(direction,expected_direction,rtol=1e-12,atol=1e-4)
            or not np.allclose(area,graph.area_weights*ell0_sq,rtol=1e-12,atol=0)):
        raise ValueError('Observer fields do not match its declared reference/state/direction/geometry')
    for gradient,key in ((Ugradient,'direction_internal_energy_rate_W'),(Fgradient,'direction_free_energy_rate_W')):
        rate=float(np.real(np.vdot(gradient,direction)))
        if not np.isclose(rate,observed[key],rtol=1e-12,atol=1e-27):
            raise ValueError('Observer gradient contraction disagrees with its completed receipt')
    if np.any(direction[graph.boundary_nodes]!=0) or np.any(gap[graph.boundary_nodes]!=base[graph.boundary_nodes]):
        raise ValueError('This audit requires fixed gap/spectral cuts; moving cuts need their work')
    return dict(graph=graph,base=base,gap=gap,direction=direction,Ugradient=Ugradient,
        Fgradient=Fgradient,profile=profile,area=area,u0=u0,anchor=ep,alpha=alpha,
        q=summary['q_bar'],summary=summary,observed=observed,Eu=Eu)


def outer_angle(z,d,q,initial):
    theta=complex(initial)
    for _ in range(60):
        residual=z*np.sin(theta)-d*np.cos(theta)+q*q*np.sin(theta)*np.cos(theta)
        if abs(residual)<2e-12:return theta
        denominator=z*np.cos(theta)+d*np.sin(theta)+q*q*np.cos(2*theta)
        update=-residual/denominator
        for reduction in range(30):
            trial=theta+2.**(-reduction)*update
            r=z*np.sin(trial)-d*np.cos(trial)+q*q*np.sin(trial)*np.cos(trial)
            if abs(r)<abs(residual):theta=trial;break
        else:raise RuntimeError('Intrinsic current-carrying retarded bulk root failed; no BCS replacement')
    raise RuntimeError('Intrinsic bulk continuation iteration budget exhausted')


def spectral_tangent(graph,d,z,root,alpha,direction,*,theta,bulk_d,q,phase):
    # Existing retarded.residual_jacobian expects complex storage even for a
    # real physical gap: its hm accumulation otherwise discards Im(f_tilde).
    d=np.asarray(d,complex);direction=np.asarray(direction,complex)
    value,J=retarded.residual_jacobian(graph,d,z,root.a,root.b,alpha)
    fixed=graph.boundary_nodes;free=np.ones(graph.n_nodes,bool);free[fixed]=False
    freecomp=np.flatnonzero(np.repeat(free,2));fixedcomp=np.flatnonzero(np.repeat(~free,2))
    factor=splu(J[freecomp][:,freecomp].tocsc())
    denominator=z*np.cos(theta)+bulk_d*np.sin(theta)+q*q*np.cos(2*theta)
    thetaE=1j*np.sin(theta)/denominator
    spectralE=.5*thetaE/(np.cos(theta/2)**2)
    knownE=np.zeros((graph.n_nodes,2),complex)
    knownE[fixed,0]=spectralE*phase;knownE[fixed,1]=spectralE*np.conj(phase)
    sources=(2j*graph.area_weights[:,None]*np.column_stack((root.a,root.b)),
        -graph.area_weights[:,None]*np.column_stack((direction.conj()*root.a**2-direction,
                                                     direction*root.b**2-direction.conj())))
    tangents=[];residuals=[]
    for source,known in zip(sources,(knownE,np.zeros_like(knownE))):
        raw=known.ravel().copy()
        raw[freecomp]=factor.solve((source.ravel()-J@raw)[freecomp])
        error=(J@raw-source.ravel()).reshape(graph.n_nodes,2)
        absolute=float(np.max(abs(error[free]),initial=0))
        relative=float(np.max(abs(error[free])/(graph.area_weights[free,None]*
            np.maximum(1.,np.max(abs(source[free]/graph.area_weights[free,None]),axis=1))[:,None]),initial=0))
        residuals.append((absolute,relative))
        aD,bD=raw.reshape(graph.n_nodes,2).T
        denominator=1+root.a*root.b;product=root.b*aD+root.a*bD
        df=2*aD/denominator-2*root.a*product/denominator**2
        dft=2*bD/denominator-2*root.b*product/denominator**2
        dg=-2*product/denominator**2
        matrix=np.empty((graph.n_nodes,2,2),complex)
        matrix[:,0,0],matrix[:,0,1],matrix[:,1,0],matrix[:,1,1]=dg,df,dft,-dg
        tangents.append(matrix)
    R=np.empty((graph.n_nodes,2,2),complex)
    R[:,0,0],R[:,0,1],R[:,1,0],R[:,1,1]=root.g,root.f,root.f_tilde,-root.g
    return ad.Jet.of(R,*tangents),dict(energy_residual_max=residuals[0][0],
        energy_residual_relative=residuals[0][1],time_residual_max_per_s=residuals[1][0],
        time_residual_relative=residuals[1][1])


def gauge_case(inputs,R,B,h,hE,kappa,alpha):
    g=inputs['graph'];xy=g.coordinates_bar
    x=(xy[:,0]-xy[:,0].min())/(xy[:,0].max()-xy[:,0].min())
    y=(xy[:,1]-xy[:,1].min())/(xy[:,1].max()-xy[:,1].min())
    chi=.2*np.sin(np.pi*x)*np.sin(2*np.pi*y);chi[g.boundary_nodes]=0
    U=ad.links(-chi)
    Rg,Bg=ad.transport(U,R,kappa),ad.transport(U,B,kappa)
    tail,head=g.edges.T;alphaG=alpha+chi[head]-chi[tail]
    answer=solve_distribution_rhs(g,Rg,h,hE,Bg,kappa_s=kappa,
        projected_source=np.zeros_like(h),source_contract='collisionless prescribed snapshot; static gauge equivalent',
        alpha=alphaG,fixed_nodes=g.boundary_nodes,fixed_rhs_per_s=np.zeros((len(g.boundary_nodes),2)))
    return answer


def query(payload):
    plan,energy,eta=payload['plan'],payload['energy'],payload['eta']
    w=load_inputs(plan);g=w['graph'];d=w['gap'];direction=w['direction'];alpha=w['alpha']
    fixed=g.boundary_nodes;phase=w['base'][fixed]/abs(w['base'][fixed])
    bulk_d=float(np.mean(abs(w['base'][fixed])))
    if not np.allclose(abs(w['base'][fixed]),bulk_d,rtol=1e-12,atol=0):
        raise ValueError('Homogeneous admitted bulk-cut amplitudes required')
    u0=w['u0'];gg=1/np.hypot(1.,abs(u0));ff=gg*u0
    a,b=ff/(1+gg),np.conj(ff)/(1+gg)
    theta=np.arctan(float(np.mean(abs(u0[fixed]))));target=eta-1j*energy
    for fraction in np.linspace(0,1,plan['continuation_steps']+1)[1:]:
        z=w['anchor']+(target-w['anchor'])*fraction
        theta=outer_angle(z,bulk_d,w['q'],theta)
        s=np.tan(theta/2)
        root=retarded.solve(g,d,z,alpha,initial_a=a,initial_b=b,fixed_nodes=fixed,
            fixed_a=s*phase,fixed_b=s*np.conj(phase),tolerance=plan['spectral_tolerance'],
            max_iterations=100,causal_tolerance=1e-8)
        a,b=root.a,root.b
    R,tangent_error=spectral_tangent(g,d,target,root,alpha,direction,
        theta=theta,bulk_d=bulk_d,q=w['q'],phase=phase)
    kappa=hbar/(2*w['Eu'])
    Bs=ad.generator(d,energy,d_time=direction,eta=eta)
    R,metrics=ad.solve_first_spectral(g,R,Bs,alpha=alpha,kappa=kappa,
        fixed_nodes=fixed,fixed_first=np.zeros((len(fixed),2,2)),
        tolerance=plan['first_gradient_tolerance'],maxiter=plan['first_gradient_maxiter'])
    B=ad.generator(d,energy,d_time=direction)
    t=w['summary']['plan']['T_K']/w['summary']['plan']['Tc_K']
    h0=np.tanh(energy/(2*t));h0E=(1-h0*h0)/(2*t)
    center=plan['population_band_center_J']/w['Eu'];width=plan['population_band_width_J']/w['Eu']
    band=np.exp(-((energy-center)/width)**2);bandE=-2*(energy-center)/width**2*band
    profile=w['profile'];amp=plan['energy_mode_increment']
    xy=g.coordinates_bar;y=(xy[:,1]-.5*(xy[:,1].min()+xy[:,1].max()))/(xy[:,1].max()-xy[:,1].min())
    ht=plan['charge_mode_increment']*profile*y
    h=np.column_stack((h0-amp*profile*band,ht*band))
    hE=np.column_stack((h0E-amp*profile*bandE,ht*bandE))
    response=solve_distribution_rhs(g,R,h,hE,B,kappa_s=kappa,
        projected_source=np.zeros_like(h),source_contract='collisionless nonthermal audit direction; no physical rates inferred',
        alpha=alpha,fixed_nodes=fixed,fixed_rhs_per_s=np.zeros((len(fixed),2)))
    gauged=gauge_case(w,R,B,h,hE,kappa,alpha)
    op=frozen.assemble(g,d,root.g,root.f,root.f_tilde,alpha)
    frozen_obs=op.evaluate(h,eta=eta)
    increment_obs=op.evaluate(np.column_stack((h[:,0]-h0,h[:,1])),eta=eta)
    N0sheet=1/(2*elementary_charge**2*w['summary']['plan']['diffusion_m2_s']*w['summary']['plan']['sheet_resistance_ohm'])
    EJ=energy*w['Eu'];weightless_measure=2*N0sheet*EJ*w['area']
    # Coefficients below are densities per joule; quadrature weights multiply in coordinator.
    moving=-float(np.dot(weightless_measure,R.time[:,0,0].real*(h[:,0]-h0)))
    population=-float(np.dot(weightless_measure,root.g.real*response.distribution_rhs_per_s[:,0]))
    energy_density=-float(np.dot(weightless_measure,root.g.real*(h[:,0]-h0)))
    force_increment_bar=increment_obs.gap_force_increment
    U0=w['observed']['energy_scale_J']
    force_work_per_J=U0/w['Eu']*float(np.real(np.vdot(force_increment_bar,direction)))
    flux=response.edge_flux_zero+response.edge_flux_first
    divergence=np.zeros((g.n_nodes,2));tail,head=g.edges.T
    np.add.at(divergence,tail,flux);np.add.at(divergence,head,-flux)
    boundary_into_per_J=2*U0/(kappa*w['Eu'])*energy*float(np.sum(divergence[fixed,0]))
    spectral,norm=ad.spectral_residual(g,R,Bs,alpha=alpha,kappa=kappa)
    free=np.ones(g.n_nodes,bool);free[fixed]=False
    record=dict(energy_J=EJ,energy_meV=EJ/elementary_charge*1e3,eta_J=eta*w['Eu'],
        rho_min=float(np.min(root.g.real)),rho_max=float(np.max(root.g.real)),
        root_residual=root.residual,spectral_tangent=tangent_error,
        first_gradient=metrics,first_over_zero=float(np.linalg.norm(R.first)/np.linalg.norm(R.value)),
        projected_rhs_residual=response.free_projected_max,linear_storage_residual=response.linear_solve_max,
        full_nambu_residual_per_mass=float(np.max(abs(response.full_nambu_residual[free]/g.area_weights[free,None,None]),initial=0)),
        frozen_noether_max=float(np.max(abs(frozen_obs.ward_residual),initial=0)),
        numerical_eta_leakage_max=float(np.max(abs(frozen_obs.artificial_eta_leakage),initial=0)),
        storage_nonzeros=response.storage_matrix_s.nnz,
        gauge_hdot_max_difference_per_s=float(np.max(abs(gauged.distribution_rhs_per_s-response.distribution_rhs_per_s),initial=0)),
        gauge_flux_max_difference=float(np.max(abs(gauged.edge_flux_zero+gauged.edge_flux_first-flux),initial=0)),
        max_abs_hLdot_per_s=float(np.max(abs(response.distribution_rhs_per_s[:,0]))),
        max_abs_hTdot_per_s=float(np.max(abs(response.distribution_rhs_per_s[:,1]))),
        pauli_margin=response.pauli_margin,
        moving_dos_rate_density_W_per_J=moving,population_rate_density_W_per_J=population,
        excitation_energy_density_dimensionless=energy_density,
        candidate_force_work_density_W_per_J=force_work_per_J,
        boundary_energy_flux_into_density_W_per_J=boundary_into_per_J,
        physical_time_steps=0,photon=False,neutrality_solved=False,finite_balance_admitted=False)
    arrays=dict(rho=root.g.real,rho_dot_per_s=R.time[:,0,0].real,h=h,h_energy=hE,
        hdot_per_s=response.distribution_rhs_per_s,R0=R.value,RE=R.energy,Rt=R.time,R1=R.first,
        edge_flux_zero=response.edge_flux_zero,edge_flux_first=response.edge_flux_first,
        projected_rhs_residual=response.projected_residual,boundary_reaction=response.boundary_reaction,
        full_nambu_residual=response.full_nambu_residual,
        storage_data=response.storage_matrix_s.data,storage_indices=response.storage_matrix_s.indices,
        storage_indptr=response.storage_matrix_s.indptr,storage_shape=response.storage_matrix_s.shape)
    return record,arrays


def energy_rule(plan,Eu):
    panels=np.asarray(plan['energy_panel_edges_J'],float)/Eu
    if panels.ndim!=1 or len(panels)<2 or panels[0]<0 or np.any(np.diff(panels)<=0):
        raise ValueError('Ordered nonnegative physical energy panels required')
    if type(plan['energy_quadrature_order']) is not int or plan['energy_quadrature_order']<1:
        raise ValueError('Positive integer physical-energy quadrature order required')
    x,w=leggauss(plan['energy_quadrature_order'])
    energy=((panels[:-1,None]+panels[1:,None])/2+np.diff(panels)[:,None]*x/2).ravel()
    weight=(np.diff(panels)[:,None]*w/2*Eu).ravel()
    return energy,weight


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--execute',action='store_true')
    args=parser.parse_args();plan=json.loads(args.plan.read_text())
    if plan.get('schema')!='pysnspd.stage5.matrix_rhs_snapshot_audit.v1':raise ValueError('Unknown audit plan')
    if (np.any(~np.isfinite(plan['eta_values_J'])) or np.any(np.asarray(plan['eta_values_J'])<=0)
            or type(plan['continuation_steps']) is not int or plan['continuation_steps']<1
            or plan['population_band_width_J']<=0):
        raise ValueError('Positive explicit eta, continuation count and population energy width required')
    inputs=load_inputs(plan);Eu=inputs['Eu'];energy,weights=energy_rule(plan,Eu)
    jobs=[dict(plan=plan,energy=float(E),eta=float(eta)/Eu,index=k,eta_index=j)
        for j,eta in enumerate(plan['eta_values_J']) for k,E in enumerate(energy)]
    resources=linux_resources();budget=resource_budget(resources,.9,.9,min(args.workers,len(jobs)))
    source_files=(Path(__file__),ROOT/'pysnspd/experimental/matrix_kinetic_rhs.py',
        ROOT/'pysnspd/experimental/adiabatic_matrix_usadel.py',ROOT/'pysnspd/experimental/retarded_spatial_usadel.py')
    identity=dict(status='DECLARED_MATRIX_RHS_SNAPSHOT_AUDIT' if args.execute else 'DRY_RUN_NO_SPECTRAL_PHYSICS',
        plan_sha256=sha(args.plan),plan=plan,
        sources={str(p.relative_to(ROOT)):sha(p) for p in source_files},
        reference_sha256=sha(plan['reference']),observer_summary_sha256=sha(Path(plan['thermal_observer'])/'summary.json'),
        observer_fields_sha256=sha(Path(plan['thermal_observer'])/'fields.npz'),
        observer_lowest_mode_sha256=sha(Path(plan['thermal_observer'])/'spectral_receipts/mode_00000.npz'),
        resources=resources,budget=budget,jobs=len(jobs),physical_time_steps=0,photon=False,
        output_exists=args.output.exists(),scope='actual DC dual mesh; prescribed finite nonthermal snapshot/direction, no heat/neutrality repair')
    if not args.execute:print(json.dumps(identity,indent=2));return 0
    if args.output.exists():raise FileExistsError('Use a fresh output directory; no overwrite or resume')
    args.output.mkdir(parents=True);(args.output/'queries').mkdir()
    write(args.output/'identity.json',identity)
    original=os.sched_getaffinity(0);os.sched_setaffinity(0,{budget['coordinator_cpu']})
    started=time.monotonic();records=[];context=mp.get_context('spawn');counter=context.Value('i',0)
    try:
        with ProcessPoolExecutor(budget['workers'],mp_context=context,initializer=initialize_affinity,
                initargs=(budget['worker_affinity_cpus'],counter)) as pool:
            pending={pool.submit(query,job):job for job in jobs}
            for future in as_completed(pending):
                job=pending[future];record,arrays=future.result();record.update(index=job['index'],eta_index=job['eta_index'])
                name=f"eta{job['eta_index']:02d}_E{job['index']:05d}"
                np.savez_compressed(args.output/'queries'/f'{name}.npz',**arrays)
                write(args.output/'queries'/f'{name}.json',record);records.append(record)
                elapsed=time.monotonic()-started
                print(json.dumps(dict(event='MATRIX_RHS_QUERY',completed=len(records),total=len(jobs),
                    elapsed_seconds=elapsed,eta_seconds=elapsed*(len(jobs)-len(records))/len(records))),flush=True)
        summaries=[]
        for j,eta in enumerate(plan['eta_values_J']):
            rows=sorted((r for r in records if r['eta_index']==j),key=lambda r:r['index'])
            sums={k:float(np.dot(weights,[r[k] for r in rows])) for k in
                ('moving_dos_rate_density_W_per_J','population_rate_density_W_per_J',
                 'candidate_force_work_density_W_per_J','boundary_energy_flux_into_density_W_per_J')}
            thermal_rate=float(np.real(np.vdot(inputs['Ugradient'],inputs['direction'])))
            thermal_force_rate=float(np.real(np.vdot(inputs['Fgradient'],inputs['direction'])))
            stored_rate=thermal_rate+sums['moving_dos_rate_density_W_per_J']+sums['population_rate_density_W_per_J']
            candidate_work=thermal_force_rate+sums['candidate_force_work_density_W_per_J']
            boundary=sums['boundary_energy_flux_into_density_W_per_J']
            summaries.append(dict(eta_J=eta,thermal_internal_rate_W=thermal_rate,
                distribution_increment_energy_J=float(np.dot(weights,[r['excitation_energy_density_dimensionless'] for r in rows])),
                moving_dos_rate_W=sums['moving_dos_rate_density_W_per_J'],
                population_rate_W=sums['population_rate_density_W_per_J'],electronic_rate_W=stored_rate,
                candidate_gap_actuator_work_W=candidate_work,boundary_energy_flux_into_W=boundary,
                independent_candidate_work_defect_W=stored_rate-candidate_work-boundary,
                maxima={k:max(r[k] for r in rows) for k in ('projected_rhs_residual','frozen_noether_max',
                    'gauge_hdot_max_difference_per_s','first_over_zero','full_nambu_residual_per_mass')},
                interpretation='Defect measured independently; candidate actuator work is F plus frozen nonthermal force, not an admitted finite-energy closure'))
        write(args.output/'summary.json',dict(status='MATRIX_RHS_SNAPSHOT_AUDIT_COMPLETED_NOT_DETECTOR_ADMISSION',
            summaries=summaries,records=records,physical_time_steps=0,photon=False,
            neutrality_solved=False,physical_collisions=False,full_recovery=False,
            finite_balance_admitted=False,elapsed_seconds=time.monotonic()-started))
    except BaseException as error:
        write(args.output/'failure.json',dict(reason=str(error),completed_queries=len(records),
            no_retry=True,no_repair=True));raise
    finally:os.sched_setaffinity(0,original)
    return 0


if __name__=='__main__':raise SystemExit(main())
