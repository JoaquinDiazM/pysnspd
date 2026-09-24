"""Post-process radial spectral results and independently check candidate variation."""
from pathlib import Path
import hashlib,json,sys
import numpy as np
from scipy.special import zeta,digamma
from scipy.integrate import trapezoid
from numpy.polynomial.legendre import leggauss
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];sys.path.insert(0,str(HERE))
import radial_usadel_reference as ref
OUT=ROOT/'docs/implementation/stage4/followup_20260923'
RAW=OUT/'raw/stage4_radial_reference_20260923'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()


def uniform_free(d,q,t,n=600):
    """Thermal free energy, not integrated from a previously computed force."""
    eps=2*np.pi*t*(np.arange(n)[:,None]+.5);gamma=q*q
    lo=np.zeros((n,len(d)));hi=d[None,:]/eps
    for _ in range(48):
        mid=(lo+hi)/2;positive=eps*mid+gamma*mid/np.sqrt(1+mid*mid)>d
        hi=np.where(positive,mid,hi);lo=np.where(positive,lo,mid)
    tangent=(lo+hi)/2;s=tangent/np.sqrt(1+tangent*tangent);c=1/np.sqrt(1+tangent*tangent)
    den=eps+gamma;tail=lambda p:zeta(p,n+.5+gamma/(2*np.pi*t))/(2*np.pi*t)**p
    lin=np.log(t)+digamma(.5+gamma/(2*np.pi*t))-digamma(.5)
    correction=d*d/den+2*eps*s*s/(1+c)-2*d*s+gamma*s*s
    return d*d*lin+2*np.pi*t*(np.sum(correction,axis=0)+d**4/4*(tail(3)-gamma*tail(4)))


def gateaux(t,gap):
    x,w=leggauss(28);panels=[.001,.01,.03,.06,.12,.25,.5,1.,2.,4.,8.]
    r=np.concatenate([(a+b)/2+(b-a)/2*x for a,b in zip(panels[:-1],panels[1:])]);weights=np.concatenate([(b-a)/2*w for a,b in zip(panels[:-1],panels[1:])])*2*np.pi*r
    d,dp,_,_=ref.profile(r,gap);checks=[]
    for delta_fraction,a,b in [(.05,.01,.5),(.1,.01,.5),(.2,.01,.5),(.1,.25,4.)]:
        inside=(r>a)&(r<b);h=np.zeros_like(r);hp=np.zeros_like(r);normal=((b-a)/2)**4
        h[inside]=(r[inside]-a)**2*(b-r[inside])**2/normal
        hp[inside]=2*(r[inside]-a)*(b-r[inside])*(a+b-2*r[inside])/normal
        delta=delta_fraction*ref.GAP_ZERO
        def energy(step):
            dd=d+step*h;ddp=dp+step*hp;q=dd*dd/((dd*dd+delta*delta)*r)
            return float(weights@(uniform_free(dd,q,t)+ref.K0*(ddp*ddp+dd*dd/(r*r)-dd*dd*q*q)))
        analytic=float(weights@(ref.candidate_force(r,gap,delta_fraction,t)['force']*h))
        estimates=[(energy(step)-energy(-step))/(2*step) for step in (2e-5,1e-5)]
        checks.append(dict(delta_fraction=delta_fraction,compact_variation_support=[a,b],
            predicted_energy_derivative=analytic,centered_energy_derivatives=estimates,
            absolute_errors=[abs(value-analytic) for value in estimates],
            relative_errors=[abs(value-analytic)/abs(analytic) for value in estimates]))
    return checks


def interval_integral(r,values,a,b):
    inside=(r>a)&(r<b);rr=np.r_[a,r[inside],b]
    yy=np.r_[np.interp(a,r,values),values[inside],np.interp(b,r,values)]
    return float(trapezoid(yy,rr))


def main():
    summary=json.loads((RAW/'summary.json').read_text());identity=json.loads((RAW/'identity.json').read_text())
    checkfiles=[dict(path=row['arrays'],matches=sha(RAW/row['arrays'])==row['sha256']) for row in summary['task_records']]
    checkfiles += [dict(path=case['fields_file'],matches=sha(RAW/case['fields_file'])==case['fields_sha256']) for case in summary['cases'].values()]
    allfiles=all(row['matches'] for row in checkfiles)
    f=np.load(RAW/'R8_N256.npz');r=f['radius'];raw=f['force_raw'];corrected=f['force_with_leading_tail_estimate'];mask=r<=4
    t=.9/8.65;gap=identity['gap_kBTc'];ell=identity['ell0_nm'];d,dp,dpp,lap=ref.profile(r,gap)
    n=512;eps=2*np.pi*t*(np.arange(n)[:,None]+.5);en=np.hypot(eps,d)
    x0=2*d*np.log(t)+4*np.pi*t*(np.sum(d**3/(eps*en*(eps+en)),axis=0)+sum(coef*zeta(power,n+.5)/(2*np.pi*t)**power for power,coef in ((3,d**3/2),(5,-3*d**5/8),(7,5*d**7/16))))
    ordinary=x0-2*ref.K0*lap
    l2=lambda v:np.sqrt(trapezoid(2*np.pi*r[mask]*v[mask]**2,r[mask]))
    candidates={}
    for fraction in (.05,.1,.2):
        candidate=f[f'candidate_force_delta_{fraction:g}'];peak=int(np.argmax(abs(candidate[mask])))
        task=np.load(RAW/'tasks'/f'R8_delta{fraction:g}.npz');q=task['q_delta'];uf=task['uniform_force'];moment=task['uniform_q_moment'];delta=fraction*ref.GAP_ZERO
        parts=dict(uniform_depaired_force=uf,complex_gradient_completion=-2*ref.K0*lap,
            phase_gradient_subtraction=-2*ref.K0*d*q*q,
            amplitude_dependence_of_q_delta=(moment-2*ref.K0*d*d*q)*2*d*delta*delta/((d*d+delta*delta)**2*r))
        difference=candidate-corrected;density=2*np.pi*r*difference*difference
        total=interval_integral(r,density,r[0],4.)
        zones=[]
        for a,b in [(r[0],.25),(.25,.5),(.5,1.),(1.,4.),(.5,4.)]:
            integral=interval_integral(r,density,a,b);den=interval_integral(r,2*np.pi*r*corrected*corrected,a,b)
            zones.append(dict(radial_interval_ell0=[a,b],radial_interval_nm=[a*ell,b*ell],squared_error_fraction_of_total=integral/total,relative_force_L2_error=math_sqrt(integral/den)))
        opposite=np.flatnonzero((candidate*corrected<0)&mask)
        splits=np.split(opposite,np.flatnonzero(np.diff(opposite)>1)+1) if len(opposite) else []
        candidates[str(fraction)]=dict(raw_relative_L2_error=float(l2(candidate-raw)/l2(raw)),
            leading_tail_relative_L2_error=float(l2(difference)/l2(corrected)),
            peak=dict(radius_ell0=float(r[peak]),radius_nm=float(r[peak]*ell),candidate=float(candidate[peak]),
                reference_raw=float(raw[peak]),reference_with_leading_tail=float(corrected[peak]),
                ordinary_BCS_plus_K0=float(ordinary[peak]),uniform_flow_embedding_correction=float((candidate-ordinary)[peak]),
                explicit_components={name:float(values[peak]) for name,values in parts.items()}),
            component_sum_max_error=float(max(abs(sum(parts.values())-candidate))),
            flow_embedding_correction_L2=float(l2(candidate-ordinary)),
            decomposition_cross_term=float(2*trapezoid(2*np.pi*r[mask]*(candidate-ordinary)[mask]*(ordinary-corrected)[mask],r[mask])),
            error_by_radial_zone=zones,
            sampled_opposite_sign_intervals_ell0=[[float(r[group[0]]),float(r[group[-1]])] for group in splits if len(group)])
    checks=gateaux(t,gap)
    result=dict(schema='pysnspd.stage4.radial_reference_analysis.v1',
        integrity=dict(summary_sha256=sha(RAW/'summary.json'),identity_sha256=sha(RAW/'identity.json'),
            registered_outputs_checked=len(checkfiles),all_registered_output_hashes_match=allfiles,
            plan_hash_matches=sha(OUT/'radial_reference_plan.json')==identity['plan_sha256'],
            current_source_hashes_match={name:sha(ROOT/name)==digest for name,digest in identity['sources'].items()}),
        measured_runtime_seconds=summary['runtime_seconds'],reference_cutoff_and_outer_comparisons=summary['comparisons'],
        reference_solver_diagnostics=dict(maximum_BVP_rms_residual=max(row['diagnostics'].get('maximum_rms_residual',0.) for row in summary['task_records']),
            maximum_BVP_boundary_residual=max(max(abs(v) for v in row['diagnostics']['boundary_residual']) for row in summary['task_records'] if row['kind']=='mode')),
        candidate_force_variation_checks=checks,
        max_candidate_Gateaux_relative_error_final_step=max(row['relative_errors'][-1] for row in checks),
        ordinary_BCS_plus_K0=dict(meaning='Counterfactual local BCS free energy plus ordinary Cartesian K0 gradient; removing the flow embedding, not an adopted model',
            force_L2=float(l2(ordinary)),relative_L2_error_against_leading_tail_reference=float(l2(ordinary-corrected)/l2(corrected))),
        candidates=candidates,
        diagnosis=dict(candidate_variation_formula_supported=True,physical_core_admitted=False,stage4_complete=False,
            neither_mesh_refinement_nor_positive_mobility_rescaling_can_fix_observed_force=True,
            full_core_model_change_applied=False),
        limits=['Prescribed thermal radial vortex, not a self-consistent vortex, barrier or strip geometry',
            'Leading-tail curve is an asymptotic estimate; cutoff and boundary comparisons are observed sensitivity, not rigorous uncertainty bounds',
            'Mean free path is unknown; subnanometer candidate peaks should not be interpreted as measured microscopic forces',
            'Mismatch remains about27% outside one ell0 and36% outside half ell0 for delta=.10, so the diagnosis is not solely a subnanometer peak',
            'This thermal free-energy comparison does not by itself construct the nonthermal kinetic closure or calibrate KWT'],
        recommendation='Freeze the current local core completion as a failed reference for low-temperature radial-core fidelity. Before advancing physical core transients, introduce a thermal spatial-response reference backend and an energy-level alternative/completion that reproduces the same-ensemble spatial response, deriving its force and current consistently. Keep the existing solver as a controlled baseline; do not tune delta or KWT to hide the discrepancy, and do not promote the ordinary BCS+K0 counterfactual.',
        analysis_script_sha256=sha(Path(__file__)))
    output=OUT/'radial_reference_analysis.json'
    if output.exists():raise FileExistsError(output)
    output.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf8')
    print(json.dumps(dict(Gateaux=checks,max_relative=result['max_candidate_Gateaux_relative_error_final_step'],integrity=result['integrity']),indent=2))


def math_sqrt(value):return float(np.sqrt(value))


if __name__=='__main__':main()
