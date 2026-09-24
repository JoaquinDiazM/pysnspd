"""Arnoldi actions of h*phi_k(hJ) on physical coordinates only.

phi_1(z)=(exp(z)-1)/z and phi_2(z)=(exp(z)-1-z)/z**2, with
their continuous values at zero. Small dense augmented exponentials evaluate
these functions without dividing by a singular J. No auxiliary clock or
constant becomes a numerically evolved physical state.

The integrated full-space defect is a backward-error diagnostic. It is not
a forward-error bound for a nonnormal operator; trajectory refinement remains
necessary. This module changes no constitutive law or physical time scale.
"""
from dataclasses import dataclass
import numpy as np
from scipy.linalg import expm


@dataclass(frozen=True)
class PhiAction:
    value: np.ndarray
    dimension: int
    integrated_defect: float
    quadrature_difference: float
    endpoint_defect: float
    tolerance: float
    accepted: bool
    action_count: int


def _small_phi_column(matrix, time, order):
    """Return phi_order(time*matrix)e_1, including singular matrices."""
    n=len(matrix)
    augmented=np.zeros((n+order,n+order))
    augmented[:n,:n]=time*matrix
    augmented[0,n]=1.
    for index in range(n,n+order-1):augmented[index,index+1]=1.
    return expm(augmented)[:n,n+order-1]


def phi_action(action, vector, interval, *, order=1, rtol=1e-3,
               atol=1e-8, max_dimension=64):
    """Approximate interval*phi_order(interval*J)vector for order 1 or 2.

    `action` receives only the physical real coordinates of a vector. The
    returned defect is integrated over the polynomially forced ODE defining
    phi_k, with the final scaling to interval*phi_k. No rejected approximation
    is silently promoted; the caller must check `accepted`.
    """
    raw=np.asarray(vector)
    if np.iscomplexobj(raw):raise ValueError('Use real Cartesian physical coordinates')
    b=np.asarray(raw,float)
    if (b.ndim!=1 or not len(b) or np.any(~np.isfinite(b)) or
        not np.isfinite(interval) or interval<=0 or order not in (1,2) or
        not np.isfinite(rtol) or not 0<rtol<1 or
        not np.isfinite(atol) or atol<=0 or type(max_dimension) is not int or max_dimension<1):
        raise ValueError('Finite nonempty vector, positive interval/tolerances, order 1/2 and positive dimension required')
    h=float(interval);beta=float(np.linalg.norm(b));limit=min(max_dimension,len(b))
    if beta==0:return PhiAction(np.zeros_like(b),0,0.,0.,0.,float(atol),True,0)
    V=np.zeros((len(b),limit+1));H=np.zeros((limit+1,limit));V[:,0]=b/beta
    calls=0
    for j in range(limit):
        raw_action=np.asarray(action(V[:,j].copy()));calls+=1
        if np.iscomplexobj(raw_action):raise RuntimeError('Operator action must use real Cartesian coordinates')
        w=np.asarray(raw_action,float)
        if w.shape!=b.shape or np.any(~np.isfinite(w)):
            raise RuntimeError('Invalid full-space operator action')
        original_norm=float(np.linalg.norm(w))
        for _ in range(2):
            projection=V[:,:j+1].T@w
            H[:j+1,j]+=projection;w-=V[:,:j+1]@projection
        H[j+1,j]=np.linalg.norm(w)
        breakdown=bool(H[j+1,j]<=64*np.finfo(float).eps*max(1.,original_norm))
        if not breakdown:V[:,j+1]=w/H[j+1,j]
        n=j+1
        if not(breakdown or n==limit or n>=8 and n%8==0):continue
        small=H[:n,:n];column=_small_phi_column(small,h,order)
        result=beta*h*(V[:,:n]@column)
        coefficient=beta*float(H[n,n-1])
        def integrated(quadrature_order):
            nodes,weights=np.polynomial.legendre.leggauss(quadrature_order)
            values=[]
            for z in nodes:
                s=h*(z+1)/2
                values.append(abs(coefficient*s**order*_small_phi_column(small,s,order)[-1]))
            return float(h**(1-order)*h/2*np.dot(weights,values))
        defect=0. if breakdown else integrated(16)
        variation=0. if breakdown else abs(defect-integrated(8))
        endpoint=abs(coefficient*h*float(column[-1]))
        tolerance=float(atol+rtol*np.linalg.norm(result))
        accepted=bool(defect+variation<=tolerance)
        if not np.all(np.isfinite(result)) or not np.isfinite(defect+variation+endpoint):
            raise RuntimeError('Nonfinite projected exponential or defect')
        if accepted or n==limit or breakdown:
            return PhiAction(result,n,defect,variation,float(endpoint),tolerance,accepted,calls)
    raise RuntimeError('Unreachable Arnoldi termination')


def frozen_affine_step(action, state, rhs_at_state, interval, **options):
    """Exact affine form x_next=x+h*phi_1(hJ)F(x), up to Arnoldi error.

    This form preserves a known constant forcing algebraically. It does not
    introduce the drifting constant coordinate of a truncated augmented state.
    """
    x=np.asarray(state,float)
    increment=phi_action(action,rhs_at_state,interval,order=1,**options)
    if x.shape!=increment.value.shape:raise ValueError('State and RHS dimensions differ')
    return x+increment.value,increment


def etd2_trial(action, rhs, state, interval, *, rhs_at_state=None, **options):
    """Exponential trapezoidal ETD2 with one fixed, full-space linear part.

    F(x)=J*x+N(x); the first-order predictor is exact for affine F. The
    correction h*phi_2(hJ)[F(p)-F(x)-J(p-x)] restores second order for the
    nonlinear remainder. No inversion of J or auxiliary constant is used.
    A caller must reject an unsuccessful action or an excessive correction;
    the correction is an embedded step-control indicator, not a rigorous
    forward-error bound. The caller also controls the physical admissibility
    of the predictor and corrected state through its RHS evaluation.
    """
    raw=np.asarray(state)
    if np.iscomplexobj(raw):raise ValueError('Use real Cartesian physical coordinates')
    x=np.asarray(raw,float)
    def checked(value):
        value=np.asarray(value)
        if np.iscomplexobj(value) or value.shape!=x.shape or np.any(~np.isfinite(value)):
            raise ValueError('A finite real full-state RHS is required')
        return value.astype(float)
    f=checked(rhs(x) if rhs_at_state is None else rhs_at_state)
    first=phi_action(action,f,interval,order=1,**options)
    if not first.accepted:return None,dict(accepted=False,first=first,second=None)
    predictor=x+first.value
    remainder=checked(rhs(predictor))-f-checked(action(first.value))
    second=phi_action(action,remainder,interval,order=2,**options)
    result=predictor+second.value if second.accepted else None
    return result,dict(accepted=second.accepted,first=first,second=second,
        predictor=predictor,correction=second.value)
