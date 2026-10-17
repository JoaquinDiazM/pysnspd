"""Finite-Q dirty-limit equilibrium free-energy correction and derivatives.

Dimensionless phi=(Omega_Q-Omega_0)/(N0*Eu^2*volume).
Inputs rho=|Delta/Eu|^2, theta=T/Tc, Q=pair-breaking energy/Eu.
Outputs: phi,rho,Q,theta,rhorho,rhoQ,QQ,rhotheta,Qtheta,thetatheta.
Native nonequilibrium kinetics are not supplied or claimed by this helper.
The off-shell Matsubara variational construction follows Clem--Kogan (2012),
https://arxiv.org/html/1207.6421 ; the independent SciPy reference is in tests.
Launchers set NUMBA_CACHE_DIR under scratch before import, retaining JIT cache
reuse without creating numerical cache files in the source checkout.
"""
from __future__ import annotations
import math
import numpy as np
try:
    from numba import njit, prange
    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False
    # Keep the numerical algorithm and its explicit failure contract identical.
    # Performance is optional; the scalar/array loops remain usable in Python.
    def njit(*args, **kwargs):
        if args and callable(args[0]):
            return args[0]
        return lambda function: function
    prange = range
_PI=math.pi
_EPS=np.finfo(float).eps
_DIGAMMA_HALF=-1.9635100260214234794
_BERNOULLI=np.array([1/6,-1/30,1/42,-1/30,5/66,-691/2730],float)
_BFAC=np.array([1/12,-1/720,1/30240,-1/1209600,1/47900160,-691/1307674368000],float)
_ZETA=np.array([1.6449340668482264,1.2020569031595942,1.0823232337111382,
1.0369277551433700,1.0173430619844492,1.0083492773819228,1.0040773561979443,
1.0020083928260822,1.0009945751278180,1.0004941886041194,1.0002460865533080,
1.0001227133475785],float)
# m,j,p,coefficient, shifted-Hurwitz argument?
_TERMS=np.array([
[2,0,3,1/4,1],[2,1,4,-1/4,1],[2,0,3,-1/4,0],
[3,0,5,1/8,0],[3,0,5,-1/8,1],[3,1,6,3/8,1],[3,2,7,-1/4,1],
[4,0,7,5/64,1],[4,0,7,-5/64,0],[4,1,8,-29/64,1],
[4,2,9,48/64,1],[4,3,10,-24/64,1],
[5,0,9,7/128,0],[5,0,9,-7/128,1],[5,1,10,65/128,1],
[5,2,11,-190/128,1],[5,3,12,220/128,1],[5,4,13,-88/128,1]],float)

@njit(cache=True)
def _digamma_data(q):
    """D(q),psi1,psi2,q*psi1,q^2*psi2; no overflow at large q."""
    if q<.01:
        D=0.;p1=0.;p2=0.
        for j in range(12):
            k=j+1
            a=(1. if k%2 else -1.)*(2.**(k+1)-1)*_ZETA[j]
            D+=a*q**k
            p1+=k*a*q**(k-1)
            if k>=2:p2+=k*(k-1)*a*q**(k-2)
        return D,p1,p2,q*p1,q*q*p2
    z=.5+q;y=z;rd=0.;r1=0.;r2=0.;rs1=0.;rs2=0.
    while y<16.:
        inv=1/y
        rd-=inv;r1+=inv*inv;r2-=2*inv**3
        rs1+=q*inv*inv;rs2-=2*(q*inv)**2*inv
        y+=1.
    inv=1/y;inv2=inv*inv
    d=math.log(y)-.5*inv
    p1=inv+.5*inv2;p2=-inv2-inv2*inv
    s1=(q/y)*(1+.5*inv);s2=-(q/y)**2*(1+inv)
    power=inv2
    for j in range(6):
        k=j+1;b=_BERNOULLI[j]
        d-=b/(2*k)*power
        p1+=b*power*inv
        p2-=(2*k+1)*b*power*inv2
        s1+=(q/y)*b*power
        s2-=(q/y)**2*(2*k+1)*b*power
        power*=inv2
    return d+rd-_DIGAMMA_HALF,p1+r1,p2+r2,s1+rs1,s2+rs2

@njit(cache=True)
def _tail_table(t,Q,N):
    """K[j,p]=Q**j * 2*pi*t sum_{n>=N}(omega_n+Q)**(-p).
    Euler-Maclaurin at z>=64; scaled representation avoids Q**j overflow.
    """
    tab=np.zeros((7,16))
    W=Q+2*_PI*t*(N+.5)
    scale=2*_PI*t/W
    rq=Q/W
    for p in range(2,16):
        Z=1/(p-1)+.5*scale
        rising=1.;power=1.
        for k in range(1,7):
            # rising (p)_(2k-1)
            if k==1:rising=p
            else:rising*=(p+2*k-3)*(p+2*k-2)
            power*=scale*scale
            Z+=_BFAC[k-1]*rising*power
        for j in range(7):
            tab[j,p]=rq**j*W**(j+1-p)*Z
    return tab

@njit(cache=True)
def _v_root(rho,omega,Q):
    if Q==0:return 1/omega
    v=1/(omega+Q)
    for iteration in range(64):
        c=1/math.sqrt(1+rho*v*v)
        val=omega*v+Q*v*c-1
        if abs(val)<=4*_EPS:return v
        H=omega+Q*c**3
        new=v-val/H
        if abs(new-v)<=4*_EPS*abs(new):
            return new
        v=new
    # A parallel worker must never escape before initializing its output.
    return np.nan

@njit(cache=True)
def _partials_one(rho,t,Q,min_terms=64,tail_scale=32.):
    if rho<0 or t<=0 or Q<0 or not math.isfinite(rho+t+Q):
        return np.full(10,np.nan)
    required=tail_scale*math.sqrt(rho)/(2*_PI*t)
    # Fail explicitly; never clamp an unresolved Matsubara sum.
    if not math.isfinite(required) or required>100000:
        return np.full(10,np.nan)
    N=max(min_terms,int(math.ceil(required)))
    if N>100000:return np.full(10,np.nan)
    q=Q/(2*_PI*t)
    D,p1,p2,s1,s2=_digamma_data(q)
    out=np.zeros(10)
    out[0]=rho*D;out[1]=D;out[2]=rho*p1/(2*_PI*t);out[3]=-rho*s1/t
    out[4]=0.;out[5]=p1/(2*_PI*t);out[6]=rho*p2/(2*_PI*t)**2
    out[7]=-s1/t
    out[8]=-rho/(2*_PI*t*t)*(p1+(s2/q if q else 0.))
    out[9]=rho/(t*t)*(2*s1+s2)
    for n in range(N):
        omega=2*_PI*t*(n+.5);W=omega+Q
        E=math.sqrt(omega*omega+rho);r0=1/E;c0=omega/E
        v=_v_root(rho,omega,Q)
        if not math.isfinite(v) or v<=0:
            return np.full(10,np.nan)
        c=1/math.sqrt(1+rho*v*v);r=v*c;H=omega+Q*c**3
        dr=r-r0
        cdiff=rho*dr*(r+r0)/(c+c0)
        B=E*rho*dr*dr*(1+rho*(r+r0)**2/(c+c0)**2)+Q*rho*r*r
        phi=2*_PI*t*B
        pr=2*_PI*t*(r0-r)
        pQ=2*_PI*t*rho*r*r
        pt=2*_PI*(B+2*omega*cdiff)
        prr=_PI*t*(omega*r**3/H-r0**3)
        prQ=2*_PI*t*r*c**3/H
        pQQ=-4*_PI*t*rho*r*r*c**3/H
        prt=2*_PI*((r0-r)+omega*(-omega*r0**3+r*c*c/H))
        pQt=2*_PI*rho*r*r*(1-2*omega*c*c/H)
        ptt=2*_PI/t*(4*omega*cdiff-2*omega*omega*rho*(r*r*c/H-r0**3))
        if Q==0:
            phi=0.;pr=0.;pt=0.;prr=0.;prt=0.;ptt=0.
        # Subtract the finite linear-rho contribution already summed exactly.
        lin=1/omega-1/W
        out[0]+=phi-2*_PI*t*rho*lin
        out[1]+=pr-2*_PI*t*lin
        out[2]+=pQ-2*_PI*t*rho/W**2
        out[3]+=pt+2*_PI*rho*Q/W**2
        out[4]+=prr
        out[5]+=prQ-2*_PI*t/W**2
        out[6]+=pQQ+4*_PI*t*rho/W**3
        out[7]+=prt+2*_PI*Q/W**2
        out[8]+=pQt-2*_PI*rho*(Q-omega)/W**3
        out[9]+=ptt-4*_PI*rho*Q*omega/(t*W**3)
    shifted=_tail_table(t,Q,N);base=_tail_table(t,0.,N)
    for row in _TERMS:
        m=int(row[0]);j=int(row[1]);p=int(row[2]);coef=row[3];shift=int(row[4])
        tab=shifted if shift else base
        K=tab[j,p]
        Kt=((1-p)*K+p*tab[j+1,p+1])/t
        Ktt=(p*(p-1)*K-2*p*p*tab[j+1,p+1]+p*(p+1)*tab[j+2,p+2])/(t*t)
        if shift:
            KQ=-p*tab[j,p+1]
            KQQ=p*(p+1)*tab[j,p+2]
            KQt=(p*p*tab[j,p+1]-p*(p+1)*tab[j+1,p+2])/t
            if j>=1:
                KQ+=j*tab[j-1,p]
                KQQ-=2*j*p*tab[j-1,p+1]
                KQt+=j*((1-p)*tab[j-1,p]+p*tab[j,p+1])/t
            if j>=2:KQQ+=j*(j-1)*tab[j-2,p]
        else:KQ=0.;KQQ=0.;KQt=0.
        a=coef*rho**m;ar=coef*m*rho**(m-1);arr=coef*m*(m-1)*rho**(m-2)
        out[0]+=a*K;out[1]+=ar*K;out[2]+=a*KQ;out[3]+=a*Kt
        out[4]+=arr*K;out[5]+=ar*KQ;out[6]+=a*KQQ
        out[7]+=ar*Kt;out[8]+=a*KQt;out[9]+=a*Ktt
    # Exact reference identity: phi(rho,theta,0)=0 for every rho/theta.
    # Keep Q derivatives: they contain the nonzero dirty-limit stiffness.
    if Q==0:
        for j in (0,1,3,4,7,9):out[j]=0.
    return out

@njit(cache=True,parallel=True)
def _partials_array(rho,t,Q,min_terms,tail_scale):
    out=np.full((10,len(rho)),np.nan)
    for i in prange(len(rho)):
        val=_partials_one(rho[i],t[i],Q[i],min_terms,tail_scale)
        for j in range(10):out[j,i]=val[j]
    return out

def finite_depairing_partials(rho,theta,Q,*,min_terms=64,tail_scale=32.):
    """Ten broadcast derivatives; worker failure is a declared Python error."""
    a,b,c=np.broadcast_arrays(np.asarray(rho,float),np.asarray(theta,float),np.asarray(Q,float))
    if (np.any(~np.isfinite(a)) or np.any(~np.isfinite(b)) or np.any(~np.isfinite(c))
            or np.any(a<0) or np.any(b<=0) or np.any(c<0)):
        raise ValueError("Finite rho>=0, theta>0, Q>=0 required")
    if (not np.isfinite(min_terms) or int(min_terms)!=min_terms or min_terms<64
            or not np.isfinite(tail_scale) or tail_scale<16):
        raise ValueError("Integer explicit modes>=64 and finite controlled tail scale>=16 required")
    values=_partials_array(a.ravel(),b.ravel(),c.ravel(),int(min_terms),float(tail_scale))
    bad=np.flatnonzero(~np.all(np.isfinite(values),axis=0))
    if len(bad):
        i=int(bad[0])
        raise ArithmeticError("Finite-Q Matsubara evaluation failed or declared cold resolution exceeded"
            +": flat_node="+str(i)+", rho="+str(float(a.ravel()[i]))
            +", theta="+str(float(b.ravel()[i]))+", Q="+str(float(c.ravel()[i])))
    return tuple(row.reshape(a.shape)for row in values)
