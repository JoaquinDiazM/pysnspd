"""High-precision independent w-plane reference for complementary E(x).

Production solves complex spectral roots then inverts their count. This check
instead solves two real equations in w, at60 decimal digits, and differentiates
that implicit system. The production solution is only an initial Newton seed.
mpmath is a temporary reviewer dependency, not a production dependency.
"""
from pathlib import Path
import json
import sys
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[4]
sys.path[:0]=[str(ROOT),str(ROOT/'tmp/stage2_review_deps')]
import mpmath as mp
from pysnspd.experimental.energy_catalog import OccupationEnergyCatalog,retarded_spectrum
from pysnspd.experimental.refined_cells import ComplementaryCountCatalog,refined_count_catalog


def sha(path):
    import hashlib
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def reference(x,amplitude,gamma,eta,seed):
    x,a,g,h=map(lambda v:mp.mpf(str(v)),(x,amplitude,gamma,eta))
    def fields(wr,wi,aa,gg):
        w=mp.mpc(wr,wi)
        u=mp.sqrt(w*w+aa*aa)
        z=u*(1-1j*gg/w)
        count=mp.re(w-1j*gg*aa*aa/(2*w*w))
        return count,mp.im(z),mp.re(z)
    def equations(wr,wi):
        v=fields(wr,wi,a,g)
        return (v[0]-x)/x,(v[1]-h)/h
    wr,wi=mp.findroot(equations,(mp.mpf(float(seed.real)),mp.mpf(float(seed.imag))),tol=mp.mpf('1e-48'),maxsteps=35)
    if wr<=0 or wi<=0:raise ArithmeticError('noncausal w-plane root')
    jac=mp.matrix([[mp.diff(lambda q:fields(q,wi,a,g)[i],wr),
                    mp.diff(lambda q:fields(wr,q,a,g)[i],wi)] for i in (0,1)])
    grad=mp.matrix([[mp.diff(lambda q:fields(q,wi,a,g)[2],wr),
                     mp.diff(lambda q:fields(wr,q,a,g)[2],wi)]])
    values=[fields(wr,wi,a,g)[2]]
    for field in ('amplitude','gamma'):
        if field=='amplitude':
            partial=[mp.diff(lambda q:fields(wr,wi,q,g)[i],a) for i in (0,1,2)]
        else:partial=[mp.diff(lambda q:fields(wr,wi,a,q)[i],g) for i in (0,1,2)]
        direction=mp.lu_solve(jac,-mp.matrix(partial[:2]))
        values.append(partial[2]+(grad*direction)[0])
    residual=max(abs(v) for v in equations(wr,wi))
    return np.array([float(v) for v in values]),float(residual)


def main():
    mp.mp.dps=60;started=time.perf_counter();out=Path(__file__).resolve().parent
    cp=ROOT/'docs/implementation/stage1_r2/catalogs/occupation_catalog.npz'
    base=OccupationEnergyCatalog.load(cp)
    points=np.array([2.113248654051871e-9,1e-7,1e-5,.001,.1,1.,11.99])
    sample=ComplementaryCountCatalog(base,points,np.ones(len(points)))
    fields=[(.72,0.),(.72,1e-12),(.72,1e-8),(.72,.2),(.72,.719999),
            (.72,.72),(.72,.720001),(.35,.65),(.08,1.2),(1.5,.002)]
    rows=[];failures=[]
    for a,g in fields:
        try:
            actual=np.asarray(sample.energy_kernel(a,g))
            c,s=retarded_spectrum(actual[0],delta=a,gamma=g,eta=base.eta)
            seeds=1j*a/s
            for i,x in enumerate(points):
                wanted,residual=reference(x,a,g,base.eta,seeds[i])
                error=abs(actual[:,i]-wanted)/np.maximum(1,abs(wanted))
                original=a*c[i]-(g*c[i]+base.eta-1j*actual[0,i])*s[i]
                rows.append(dict(amplitude=a,gamma=g,count=float(x),actual=actual[:,i].tolist(),
                    independent_reference=wanted.tolist(),scaled_errors=error.tolist(),
                    reference_equation_relative_residual=residual,
                    original_spectral_scaled_residual=float(abs(original)/max(1,a,g,actual[0,i]))))
        except Exception as exc:
            failures.append(dict(amplitude=a,gamma=g,error=repr(exc)))
    mesh=[]
    for level in (1,2,4):
        cat=refined_count_catalog(base,level)
        energy,_,_=cat.energy_kernel(.72,0.)
        mesh.append(dict(refinement=level,nodes=len(cat.count_nodes),count_min=float(cat.count_nodes[0]),
            count_max=float(cat.count_nodes[-1]),total_count_weight=float(np.sum(cat.count_weights)),
            quadrature_count_moment_error=float(abs(np.dot(cat.count_weights,cat.count_nodes)-72)),
            minimum_energy_increment=float(np.min(np.diff(energy))),
            arrays_readonly=not cat.count_nodes.flags.writeable and not cat.count_weights.flags.writeable))
    maxima=dict(energy=max(r['scaled_errors'][0] for r in rows),
                amplitude_derivative=max(r['scaled_errors'][1] for r in rows),
                gamma_derivative=max(r['scaled_errors'][2] for r in rows),
                spectral=max(r['original_spectral_scaled_residual'] for r in rows))
    gates=dict(no_failures=not failures,energy=maxima['energy']<=1e-5,
               derivatives=max(maxima['amplitude_derivative'],maxima['gamma_derivative'])<=1e-7,
               causal=maxima['spectral']<=1e-8,positive_ordered_mesh=all(r['minimum_energy_increment']>0 for r in mesh))
    data=dict(schema='pysnspd.stage2.independent_complementary_energy.v1',
        status='PASS' if all(gates.values()) else 'FAIL',gates=gates,maxima=maxima,
        method='Solve x=Re[w-iGammaDelta^2/(2w^2)] and eta=Im[sqrt(w^2+Delta^2)*(1-iGamma/w)] at60 digits; implicit Jacobian differentiation in(wr,wi)',
        derivative_order=['energy','fixed_count_amplitude_derivative','fixed_count_Gamma_derivative'],
        rows=rows,failures=failures,meshes=mesh,mpmath_version=mp.__version__,decimal_digits=60,
        catalog_sha256=sha(cp),wrapper_sha256=sha(ROOT/'pysnspd/experimental/refined_cells.py'),
        source_sha256=sha(Path(__file__)),criteria_sha256=sha(ROOT/'docs/implementation/stage2/acceptance_criteria.json'),
        amendment_sha256=sha(out/'complementary_grid_amendment.json'),runtime_seconds=time.perf_counter()-started)
    (out/'complementary_energy_reference.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=data['status'],maxima=maxima,failures=failures,runtime_seconds=data['runtime_seconds'])))
    if not all(gates.values()):raise SystemExit(1)


if __name__=='__main__':main()
