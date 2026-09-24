"""Matrix-free exponential actions with a measured full-space ODE defect.

Arnoldi vectors span a numerical approximation space, not a physical truncation
to prescribed profiles. The integrated defect is a local backward-error measure;
it is not a rigorous forward-error bound for an arbitrary nonnormal operator.
"""
from dataclasses import dataclass
import numpy as np
from scipy.linalg import expm


@dataclass(frozen=True)
class ExponentialStep:
    state: np.ndarray
    derivative: np.ndarray
    dimension: int
    integrated_defect: float
    quadrature_difference: float
    endpoint_defect: float
    tolerance: float
    accepted: bool
    action_count: int


def exponential_step(action,state,dt,*,rtol=1e-3,atol=1e-7,max_dimension=64,active_size=None):
    x=np.asarray(state,float)
    if x.ndim!=1 or np.any(~np.isfinite(x)) or not dt>0 or not 0<rtol<1 or atol<=0 or max_dimension<2:
        raise ValueError('Finite real state, positive interval/tolerances and dimension>=2 required')
    active_size=len(x) if active_size is None else int(active_size)
    if not 0<active_size<=len(x):raise ValueError('Invalid physical coordinate count')
    beta=float(np.linalg.norm(x));maximum=min(max_dimension,len(x))
    if beta==0:
        return ExponentialStep(x.copy(),x.copy(),0,0.,0.,0.,atol,True,0)
    basis=np.zeros((len(x),maximum+1));H=np.zeros((maximum+1,maximum));basis[:,0]=x/beta
    calls=0
    def defect_integral(matrix,size,coefficient,order):
        nodes,weights=np.polynomial.legendre.leggauss(order)
        return dt/2*sum(w*abs(coefficient*expm(dt*(z+1)/2*matrix)[size-1,0]) for z,w in zip(nodes,weights))
    for j in range(maximum):
        v=np.asarray(action(basis[:,j]),float);calls+=1
        if v.shape!=x.shape or np.any(~np.isfinite(v)):raise RuntimeError('Invalid full-space operator action')
        original=float(np.linalg.norm(v))
        for _ in range(2):
            coefficients=basis[:,:j+1].T@v
            H[:j+1,j]+=coefficients;v-=basis[:,:j+1]@coefficients
        H[j+1,j]=np.linalg.norm(v)
        breakdown=H[j+1,j]<=64*np.finfo(float).eps*max(1.,original)
        if not breakdown:basis[:,j+1]=v/H[j+1,j]
        size=j+1
        if not(breakdown or size==maximum or size>=8 and size%8==0):continue
        small=H[:size,:size];coefficients=beta*expm(dt*small)[:,0]
        result=basis[:,:size]@coefficients
        derivative=basis[:,:size]@(small@coefficients)
        factor=beta*H[size,size-1]
        integral=0. if breakdown else defect_integral(small,size,factor,16)
        difference=0. if breakdown else abs(integral-defect_integral(small,size,factor,8))
        endpoint=float(abs(factor*expm(dt*small)[size-1,0]))
        # Exclude the exact constant affine coordinate from physical accuracy.
        tolerance=atol+rtol*max(float(np.linalg.norm(x[:active_size])),float(np.linalg.norm(result[:active_size])))
        accepted=bool(integral+difference<=tolerance)
        if accepted or size==maximum or breakdown:
            # Independent full-space evaluation, including discarded Arnoldi
            # directions; the analytic defect above is not the only check.
            endpoint=float(np.linalg.norm(np.asarray(action(result),float)-derivative));calls+=1
            return ExponentialStep(result,derivative,size,float(integral),float(difference),endpoint,tolerance,accepted,calls)
    raise RuntimeError('Unreachable Arnoldi termination')


def affine_action(action,forcing):
    """Exact augmentation of xdot=Jx+b; the last coordinate is constant one."""
    b=np.asarray(forcing,float)
    def augmented(value):
        v=np.asarray(value,float)
        if v.shape!=(len(b)+1,):raise ValueError('Invalid augmented coordinate')
        return np.r_[np.asarray(action(v[:-1]))+v[-1]*b,0.]
    return augmented


def advance(action,state,interval,*,rtol=1e-3,atol=1e-7,max_dimension=64,maximum_steps=2048,progress=None,active_size=None):
    """Adaptive numerical restarts; rejected steps never change physical state."""
    x=np.asarray(state,float).copy();elapsed=0.;step=float(interval);history=[];attempts=0
    if not interval>0:raise ValueError('Positive observation interval required')
    while elapsed<interval:
        attempts+=1
        if attempts>maximum_steps:raise RuntimeError('Arnoldi restart budget exhausted; no silent accuracy reduction')
        step=min(step,interval-elapsed)
        local_atol=atol*step/interval
        result=exponential_step(action,x,step,rtol=rtol*step/interval,atol=local_atol,max_dimension=max_dimension,active_size=active_size)
        row=dict(start=elapsed,step=step,dimension=result.dimension,integrated_defect=result.integrated_defect,
            quadrature_difference=result.quadrature_difference,endpoint_defect=result.endpoint_defect,
            tolerance=result.tolerance,accepted=result.accepted,operator_actions=result.action_count)
        history.append(row)
        if progress is not None:progress(row,result.state if result.accepted else None)
        if result.accepted:
            x=result.state;elapsed+=step
            step=min(interval-elapsed,step*1.5)
        else:
            step/=2
            if step<=np.finfo(float).eps*max(1.,interval):raise RuntimeError('Time step underflow; no accuracy fallback')
    return x,history
