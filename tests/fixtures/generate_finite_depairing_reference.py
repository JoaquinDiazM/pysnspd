"""Regenerate the finite-depairing mathematical reference with SciPy.

Independent implementation: monotone native v root, variational-envelope
derivatives, 4096 explicit Matsubara terms and a rho**5 Hurwitz tail.
No pysnspd functions or production summation helper are imported.

omega_n=(2*n+1)*pi*t; c=(1+rho*v**2)**(-1/2).
v solves omega*v+Q*v*c=1. B is the positive free-energy difference
E*rho*(r-r0)**2*(1+rho*(r+r0)**2/(c+c0)**2)+Q*rho*r**2.
The exactly summed rho term uses scipy.special.digamma. Higher tail
coefficients use scipy.special.zeta, with analytic envelope derivatives.

This validates the specified equilibrium variational mathematics. It does
not establish admission of a nonequilibrium SNSPD constitutive model.
"""
import numpy as np
from scipy.special import digamma, polygamma, zeta
FIELDS=("phi","rho","Q","t","rhorho","rhoQ","QQ","rhot","Qt","tt")
TERMS=((2,.25,0,3,False),(2,-.25,1,4,False),(2,-.25,0,3,True),
 (3,.125,0,5,True),(3,-.125,0,5,False),(3,.375,1,6,False),(3,-.25,2,7,False),
 (4,5/64,0,7,False),(4,-5/64,0,7,True),(4,-29/64,1,8,False),(4,.75,2,9,False),(4,-.375,3,10,False),
 (5,7/128,0,9,True),(5,-7/128,0,9,False),(5,65/128,1,10,False),(5,-95/64,2,11,False),(5,55/32,3,12,False),(5,-11/16,4,13,False))
def linear(rho,t,Q):
    q=Q/(2*np.pi*t);z=.5+q
    D=digamma(z)-digamma(.5)
    p1,p2=polygamma(1,z),polygamma(2,z)
    DQ=p1/(2*np.pi*t);DQQ=p2/(2*np.pi*t)**2
    Dt=-q*p1/t;DQt=-(p1+q*p2)/(2*np.pi*t*t);Dtt=(2*q*p1+q*q*p2)/(t*t)
    return np.array([rho*D,D,rho*DQ,rho*Dt,0.,DQ,rho*DQQ,Dt,rho*DQt,rho*Dtt])
def tail(rho,t,Q,N):
    C=2*np.pi*t;out=np.zeros(10)
    cache={}
    def hs(p,reference):
        key=(p,reference)
        if key in cache:return cache[key]
        shift=0 if reference else Q
        h=[C**(1-k)*zeta(k,N+.5+shift/C) for k in (p,p+1,p+2)]
        H,Hp1,Hp2=h
        Ht=((1-p)*H+p*shift*Hp1)/t
        Htt=(p*(p-1)*H-2*p*p*shift*Hp1+p*(p+1)*shift*shift*Hp2)/(t*t)
        if reference:HQ=HQQ=HQt=0.
        else:HQ=-p*Hp1;HQQ=p*(p+1)*Hp2;HQt=(p*p*Hp1-p*(p+1)*Q*Hp2)/t
        cache[key]=(H,HQ,Ht,HQQ,HQt,Htt);return cache[key]
    for m,c,j,p,reference in TERMS:
        H,HQ,Ht,HQQ,HQt,Htt=hs(p,reference)
        a0=rho**m;a1=m*rho**(m-1);a2=m*(m-1)*rho**(m-2)
        q0=Q**j;q1=j*Q**(j-1) if j else 0.;q2=j*(j-1)*Q**(j-2) if j>1 else 0.
        out+=c*np.array([a0*q0*H,a1*q0*H,a0*(q1*H+q0*HQ),a0*q0*Ht,a2*q0*H,
            a1*(q1*H+q0*HQ),a0*(q2*H+2*q1*HQ+q0*HQQ),a1*q0*Ht,a0*(q1*Ht+q0*HQt),a0*q0*Htt])
    return out
def finite(rho,t,Q,N):
    w=(2*np.arange(N)+1)*np.pi*t
    v=1/(w+Q)
    for _ in range(32):
        c=1/np.sqrt(1+rho*v*v);H=w+Q*c*c*c
        dv=(1-w*v-Q*v*c)/H
        v+=dv
        if np.max(abs(dv)/(1+abs(v)))<2*np.finfo(float).eps:break
    c=1/np.sqrt(1+rho*v*v);H=w+Q*c*c*c;r=v*c
    E=np.sqrt(w*w+rho);r0=1/E;c0=w/E
    diff=r-r0
    B=E*rho*diff*diff*(1+rho*(r+r0)**2/(c+c0)**2)+Q*rho*r*r
    Br=r0-r;Brr=.5*(w*r**3/H-r0**3);BQ=rho*r*r
    BrQ=r*c**3/H;BQQ=-2*rho*r*r*c**3/H
    Bw=-2*rho*(r0-r)*(r0+r)/(c+c0)
    Brw=-w*r0**3+r*c*c/H;BQw=-2*rho*r*r*c*c/H
    Bww=-2*rho*(r*r*c/H-r0**3)
    # Subtract the linearly resummed term before performing the thermal sums.
    W=w+Q;A=Q/(w*W)
    Aw=-(Q*(w+W))/(w*w*W*W)
    Aww=2*Q*(W*W+W*w+w*w)/(w**3*W**3)
    B-=rho*A;Br-=A;BQ-=rho/W**2;BrQ-=1/W**2;BQQ+=2*rho/W**3
    Bw-=rho*Aw;Brw-=Aw;BQw+=2*rho/W**3;Bww-=rho*Aww
    return np.array([2*np.pi*t*np.sum(B),2*np.pi*t*np.sum(Br),2*np.pi*t*np.sum(BQ),
      2*np.pi*np.sum(B+w*Bw),2*np.pi*t*np.sum(Brr),2*np.pi*t*np.sum(BrQ),2*np.pi*t*np.sum(BQQ),
      2*np.pi*np.sum(Br+w*Brw),2*np.pi*np.sum(BQ+w*BQw),(2*np.pi/t)*np.sum(2*w*Bw+w*w*Bww)])
def phi(rho,t,Q,N):
    return linear(rho,t,Q)+finite(rho,t,Q,N)+tail(rho,t,Q,N)


REFERENCE_CASES=(
    (0.,.05,0.),
    (0.,.05,1.),
    (1e-12,.05,1e6),
    (4.,.05,.01),
    (1.,.104,1.),
    (.1,1.,100.),
    (4.,5.,1e6),
    (1.,.104,0.),
)
def generate(output):
    import json,hashlib
    from pathlib import Path
    count=4096
    cases=[dict(rho=rho,theta=theta,Q=Q,values=phi(rho,theta,Q,count).tolist())
        for rho,theta,Q in REFERENCE_CASES]
    data=dict(schema="pysnspd.finite_depairing.reference.v1",fields=list(FIELDS),
        reference=dict(generator=Path(__file__).name,
            generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            explicit_matsubara_terms=count,tail="rho^5 all-Q Hurwitz-zeta expansion",
            linear_term="rho*(digamma(1/2+Q/(2*pi*theta))-digamma(1/2))",
            root="omega*v+Q*v/sqrt(1+rho*v^2)=1",
            derivation="Stationary Matsubara variational envelope and its Hessian",
            scope="Mathematical summation oracle; no nonequilibrium physical admission"),
        cases=cases)
    Path(output).write_text(json.dumps(data,indent=2)+"\n")
if __name__=="__main__":
    import argparse
    from pathlib import Path
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=Path(__file__).with_name("finite_depairing_reference.json"))
    args=parser.parse_args();generate(args.output);print(args.output)
