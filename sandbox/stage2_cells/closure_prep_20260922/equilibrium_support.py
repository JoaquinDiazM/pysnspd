"""Independent high-energy absorption bound for the finite-Gamma equilibrium.

Extends only the reference domain, not the dynamics. The same real-angle cubic
as the gapless reference is inverted well ABOVE the gapped threshold; original
Usadel residuals and finite-eta bias are checked independently at every node.
"""
from pathlib import Path
import argparse,hashlib,json,sys,time
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT),str(ROOT/'sandbox/stage2_cells'),str(ROOT/'docs/implementation/stage2/review')]
import numpy as np
from continuous_reactions import real_angle,segments
from run_coupled import setup
from pysnspd.experimental.cell_validation import ElectronicCell
from pysnspd.experimental.energy_catalog import retarded_spectrum
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()

def spectrum(energy,a,g):
    e=np.asarray(energy,float)
    if np.any(e<=max(2*a,4.)):
        raise ValueError('This independent extension is restricted to the high-energy tail')
    if g==0:
        den=np.sqrt(e*e-a*a)
        return e/den,a/den,0.
    lo=np.zeros_like(e);hi=np.log1p((e/g)**2)+4
    t=np.log1p((e/g)**2)
    for _ in range(50):
        actual,derivative,*_=real_angle(t,a,g)
        relative=abs(actual-e)/e
        if float(relative.max())<3e-13:break
        lo=np.where(actual<e,t,lo);hi=np.where(actual>=e,t,hi)
        proposal=t-(actual-e)/derivative
        update=np.where((proposal>lo)&(proposal<hi),proposal,(lo+hi)/2)
        t=np.where(relative<3e-13,t,update)
    else:raise ArithmeticError('Independent high-energy inversion did not converge')
    _,_,s,y,v,k=real_angle(t,a,g)
    c=np.sqrt(k)*y-1j*s*v;f=s*y+1j*np.sqrt(k)*v
    residual=float(np.max(abs(a*c-(g*c-1j*e)*f)/np.maximum(1,e)))
    if residual>1e-11:raise ArithmeticError('Independent Usadel equation residual failed')
    return c.real,f.imag,residual

def bound(a,g,emax,population,phonons,alpha,eta,order):
    cuts=[.005,.02,.05,.1,.2,.5,1.,1.5,2.,2.5,3.,4.]
    om,wo=segments(cuts,order,square=True);z,w=np.polynomial.legendre.leggauss(order)
    ei=emax-om[:,None]+om[:,None]*(z[None,:]+1)/2;eo=ei+om[:,None]
    ni,ri,erri=spectrum(ei,a,g);no,ro,erro=spectrum(eo,a,g)
    coherence=ni*no-ri*ro
    ci,fi=retarded_spectrum(ei,delta=a,gamma=g,eta=eta)
    co,fo=retarded_spectrum(eo,delta=a,gamma=g,eta=eta)
    finite=ci.real*co.real-fi.imag*fo.imag
    if np.any(coherence<=0) or np.any(~np.isfinite(finite)) or np.any(finite<=0):
        raise ValueError('Nonpositive or invalid tail coherence')
    regulator_bias=float(np.max(abs(finite-coherence)/coherence))
    if regulator_bias>1e-7:raise ValueError('Finite-eta coherence differs beyond preregistered budget')
    p=population(ei);n=phonons(om)
    if np.any((p<0)|(p>1)) or np.any(n<0):raise ValueError('Invalid tail populations')
    integrand=wo[:,None]*om[:,None]*w[None,:]/2*coherence*alpha(om)[:,None]*p*n[:,None]
    return dict(number=float(integrand.sum()),power=float(np.sum(integrand*om[:,None])),
                spectral_residual=max(erri,erro),eta_relative_bias=regulator_bias,
                minimum_sampled_energy=float(ei.min()))

def contract(trajectory):
    r=json.loads(trajectory.read_text())
    paths=[Path(__file__),ROOT/'docs/implementation/stage2/review/continuous_reactions.py',
           ROOT/'docs/implementation/stage1_r2/review/causal_reference.py',
           ROOT/'docs/implementation/stage1_r2/review/gamma_zero_reference.py',
           ROOT/'docs/implementation/stage2/acceptance_criteria.json']
    return dict(schema='pysnspd.stage2.equilibrium-support-registration.v1',
        trajectory_sha256=sha(trajectory),archive_sha256=sha(trajectory.with_suffix('.npz')),
        source_hashes={**r['source_hashes'],**{p.relative_to(ROOT).as_posix():sha(p) for p in paths}},
        quadrature_orders=[16,32],minimum_times_per_cell=10,upper_activity_relative_limit=1e-3,
        reference_refinement_relative_limit=2e-4,causal_residual_scaled_limit=1e-11,
        finite_eta_coherence_relative_budget=1e-7,
        scope='Absorption into unresolved high-energy states; external hole bounded by1. No bound on arbitrary incoming unrepresented populations.',
        reference='Ideal real-angle Usadel solution in E>max(2Delta,4). Original equation checked at all quadrature nodes. Finite-eta difference separately measured.',
        execution=dict(expected_seconds=[5,30],timeout_seconds=240,new_trajectories=0))

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('trajectory',type=Path)
    p.add_argument('--registration',type=Path,required=True);p.add_argument('--output',type=Path)
    p.add_argument('--prepare',action='store_true');args=p.parse_args()
    current=contract(args.trajectory)
    if args.prepare:
        with args.registration.open('x',encoding='utf-8') as f:json.dump(current,f,indent=2)
        print('PREREGISTERED_NO_CALCULATION');return
    if json.loads(args.registration.read_text())!=current:raise ValueError('Contract changed')
    if args.output is None or args.output.exists():raise ValueError('A fresh output is required')
    started=time.perf_counter();record=json.loads(args.trajectory.read_text());par=record['parameters']
    if par['scenario']!='equilibrium':raise ValueError('Equilibrium scenario required')
    system,_,debye=setup(par['case'],par['phonon_nodes'],par['infrared'],par['face_order'],float(par['escape']),par['heating'],par['scenario'],par['reaction_order'],par['reaction_method'],par['electron_refinement'],par['reaction_layout'],par['reaction_outer_order'],par['reaction_max_panel'],par['reaction_max_energy_panel'])
    with np.load(args.trajectory.with_suffix('.npz'),allow_pickle=False) as data:
        times,states=data['times'],data['states']
    if len(np.unique(times))<10:raise ValueError('At least10 distinct times required')
    rows=[];cache={}
    for ti,state in zip(times,states):
        amplitudes,populations,phonon=system.unpack(state)
        for i,a in enumerate(amplitudes):
            g=system.gammas[i];cell=ElectronicCell(system.catalog,a,g)
            key=(float(a),g,populations[i].tobytes(),phonon[i].tobytes())
            if key not in cache:
                network=system.reaction_events(cell);lf,lb=network.log_activities(populations[i],phonon[i])
                gross=network.coefficients*(np.exp(lf)+np.exp(lb))
                denominators=dict(number=float(gross.sum()),power=float(gross@network.omega))
                refinements=[bound(a,g,cell.energies[-1],lambda e:np.interp(e,cell.energies,populations[i]),
                    lambda om:np.interp(om,system.phonons.energies,phonon[i]),debye.alpha2F,system.catalog.eta,order) for order in (16,32)]
                relative={k:system.rate_prefactor*refinements[-1][k]/denominators[k] for k in ('number','power')}
                difference={k:abs(refinements[-1][k]-refinements[0][k])/abs(refinements[-1][k]) for k in ('number','power')}
                cache[key]=dict(relative_upper_bounds=relative,reference_relative_changes=difference,
                    reference_orders=[16,32],bounds=refinements,gross_activity=denominators)
            rows.append(dict(time=float(ti),cell=i,amplitude=float(a),gamma=g,**cache[key]))
    maximum=max(v for row in rows for v in row['relative_upper_bounds'].values())
    refinement=max(v for row in rows for v in row['reference_relative_changes'].values())
    status='PASS' if maximum<=1e-3 and refinement<=2e-4 else 'FAIL'
    result=dict(status=status,registration_sha256=sha(args.registration),source_sha256=sha(__file__),
        trajectory_sha256=sha(args.trajectory),samples=len(rows),distinct_times_per_cell=len(np.unique(times)),
        worst_relative_upper_bound=maximum,worst_reference_relative_change=refinement,
        scope=current['scope'],reference=current['reference'],rows=rows,runtime_seconds=time.perf_counter()-started,
        new_trajectories=0,stage2_admission=False)
    with args.output.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2,allow_nan=False)
    print(json.dumps({k:result[k] for k in ('status','samples','worst_relative_upper_bound','worst_reference_relative_change','runtime_seconds')}))
    if status!='PASS':raise SystemExit(1)

if __name__=='__main__':main()
