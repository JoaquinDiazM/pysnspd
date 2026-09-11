"""Independent 70-digit check of Hermite arithmetic, not an observable-current test.

The real-spectrum cases compare d/dGamma log(E_i-E_{i-1}) at the midpoint
of a single Gamma interval. These are slopes at individual quadrature nodes.
They are NOT errors of the integrated supercurrent, kinetic rates or a detector
transient. The compensated test retains the high-precision endpoint difference
before converting it to double precision; it diagnoses representation loss.

Review-only dependencies: NumPy, SciPy and mpmath 1.3.0. No constructor or
production module is changed. Results are written beside this script.
"""
from __future__ import annotations
import json,sys,time,hashlib
from pathlib import Path
import numpy as np
import mpmath as mp
ROOT=next(parent for parent in Path(__file__).resolve().parents
          if (parent/'pysnspd/experimental/energy_catalog.py').is_file())
sys.path.insert(0,str(ROOT))
from pysnspd.experimental.energy_catalog import energy_at_count_batch,retarded_spectrum_batch
mp.mp.dps=70

def naive(v0,v1,m0,m1,h,t=.5):
    return ((6*t*t-6*t)*v0/h+(3*t*t-4*t+1)*m0+
            (-6*t*t+6*t)*v1/h+(3*t*t-2*t)*m1)

def stable(dv,m0,m1,h,t=.5):
    return 6*t*(1-t)*(dv/h)+(1-4*t+3*t*t)*m0+(-2*t+3*t*t)*m1

def mp_node(x,gamma,a,eta):
    e0=x*mp.sqrt(1+a*a/(x*x+eta*eta))
    if gamma==0:
        u=mp.mpc(e0,eta)
    else:
        def equations(ur,ui):
            u=mp.mpc(ur,ui);w=mp.sqrt(u*u-a*a)
            z=u-1j*gamma*u/w
            count=mp.re(w-1j*gamma*a*a/(2*w*w))
            return mp.im(z)-eta,count-x
        ur,ui=mp.findroot(equations,(e0,eta),tol=mp.mpf('1e-62'),maxsteps=40)
        u=mp.mpc(ur,ui)
    w=mp.sqrt(u*u-a*a);z=u-1j*gamma*u/w
    c=u/w;s=1j*a/w
    return mp.re(z),-mp.re(s)*mp.im(s)/mp.re(c)

def main():
    start=time.perf_counter()
    linear=[]
    for slope in ('0.123','-1e-8','10'):
        v0=mp.mpf('-10');m=mp.mpf(slope);h=mp.mpf('2.45e-14');v1=v0+h*m
        f0,f1,fm,fh=map(float,(v0,v1,m,h))
        linear.append(dict(slope=float(m),rounded_endpoints_equal=f0==f1,
                           naive=naive(f0,f1,fm,fm,fh),
                           stable_rounded_endpoints=stable(f1-f0,fm,fm,fh),
                           stable_compensated_difference=stable(float(v1-v0),fm,fm,fh)))
    a,eta=mp.mpf('.72'),mp.mpf('1e-8')
    xs=[mp.mpf(s) for s in ('1e-9','1e-8','1e-7','1e-6','1e-5','1e-4','.001','.01','.1','1')]
    node0=[mp_node(x,mp.mpf(0),a,eta) for x in xs]
    e0=[v[0] for v in node0];s0=[v[1] for v in node0]
    records=[]
    for hs in ('2.45e-14','1e-12','1e-10'):
        h=mp.mpf(hs);node1=[mp_node(x,h,a,eta) for x in xs]
        e1=[v[0] for v in node1];s1=[v[1] for v in node1]
        batch0=energy_at_count_batch(np.array(xs,float),delta=float(a),gamma=0,eta=float(eta))
        batch1=energy_at_count_batch(np.array(xs,float),delta=float(a),gamma=float(h),eta=float(eta))
        bc0,bs0=retarded_spectrum_batch(batch0,delta=float(a),gamma=0,eta=float(eta))
        bc1,bs1=retarded_spectrum_batch(batch1,delta=float(a),gamma=float(h),eta=float(eta))
        bm0=-bs0.real*bs0.imag/bc0.real;bm1=-bs1.real*bs1.imag/bc1.real
        binc0=np.diff(batch0,prepend=0);binc1=np.diff(batch1,prepend=0)
        bd0=np.diff(bm0,prepend=0)/binc0;bd1=np.diff(bm1,prepend=0)/binc1
        for i,x in enumerate(xs):
            inc0=e0[i]-(e0[i-1] if i else 0);inc1=e1[i]-(e1[i-1] if i else 0)
            v0,v1=mp.log(inc0),mp.log(inc1)
            m0=(s0[i]-(s0[i-1] if i else 0))/inc0
            m1=(s1[i]-(s1[i-1] if i else 0))/inc1
            reference=stable(v1-v0,m0,m1,h,mp.mpf('.5'))
            f0,f1,fm0,fm1,fh=map(float,(v0,v1,m0,m1,h))
            n=naive(f0,f1,fm0,fm1,fh);s=stable(f1-f0,fm0,fm1,fh)
            comp=stable(float(v1-v0),fm0,fm1,fh)
            actual=stable(np.log(binc1[i])-np.log(binc0[i]),bd0[i],bd1[i],fh)
            denom=max(1,abs(float(reference)))
            records.append(dict(x=float(x),h=float(h),true_log_change=float(v1-v0),
                                reference_Hermite_derivative=float(reference),
                                rounded_log_endpoints_equal=f0==f1,
                                naive_scaled_error=abs(n-float(reference))/denom,
                                stable_rounded_scaled_error=abs(s-float(reference))/denom,
                                compensated_scaled_error=abs(comp-float(reference))/denom,
                                actual_batch_stable_scaled_error=abs(actual-float(reference))/denom,
                                absolute_batch_energy_error=float(abs(mp.mpf(float(batch1[i]))-e1[i]))))
    result=dict(precision_decimal_digits=mp.mp.dps,mpmath_version=mp.__version__,
                quantity="d/dGamma log(E_i-E_(i-1)) at individual count-quadrature nodes",
                error_normalization="absolute derivative error divided by max(1, abs(reference derivative))",
                interpretation="Representation/arithmetic control only; not an integrated-current or transient error.",
                compensated_reference="High-precision delta(log increment) stored separately before conversion to double precision.",
                runtime_seconds=time.perf_counter()-start,
                module_sha256=hashlib.sha256((ROOT/'pysnspd/experimental/energy_catalog.py').read_bytes()).hexdigest(),
                linear_cases=linear,real_spectrum_cases=records,
                max_errors={key:max(r[key] for r in records) for key in
                            ('naive_scaled_error','stable_rounded_scaled_error','compensated_scaled_error','actual_batch_stable_scaled_error')})
    path=Path(__file__).with_name('hermite_precision.json')
    path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='real_spectrum_cases'},indent=2))
if __name__=='__main__':main()
