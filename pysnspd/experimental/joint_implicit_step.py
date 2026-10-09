"""Constrained backward-Euler step for a complete, explicitly supplied DAE.

For F(y, ydot)=0 this solves F(y_new, D(y_new-y_old)/dt)=0, retaining
all supplied rows. D is the declared differential-state mask; algebraic
variables have zero derivative. The caller must evaluate the instantaneous
spectra, distributions, collision kernels and thermal reference at y_new.
A frozen snapshot mobility matrix plus a fitted diagonal is not this DAE.

Only discretization is provided here: no new gap friction, bath, thermal
closure, or numerical spectral broadening is introduced. A result is a
nonlinear solver result, never a detector, energy or temporal admission.
initial_residual_max records F(y_old,0), a stationary residual rather than
an initial-consistency certificate: a valid differential state can have a
nonzero physical derivative. Backward Euler does not establish that an
algebraic initial projection is physically valid. The full first-gradient model must still pass its own
small-gradient and step-refinement checks.

Pauli-safe use: represent electronic modes by u_plus=hL+hT and
u_minus=hL-hT, with bounds [-1,1]; use n_ph>=0 for each phonon mode.
Bounds constrain trial states via a fraction-to-boundary line search. No
state, population, rate, row or residual is clipped or repaired.
"""
from __future__ import annotations
from dataclasses import dataclass
import warnings
import numpy as np
from scipy.linalg import LinAlgWarning, solve
from scipy.sparse import diags, issparse
from scipy.sparse.linalg import MatrixRankWarning, splu


class ImplicitStepRejected(RuntimeError):
    """No accepted discrete state; preserve the previous state and diagnostics."""


@dataclass(frozen=True)
class ImplicitStepResult:
    state: np.ndarray
    derivative_per_s: np.ndarray
    scaled_residual: np.ndarray
    iterations: int
    residual_evaluations: int
    jacobian_evaluations: int
    initial_residual_max: float
    accepted_residual_max: float
    minimum_bound_margin_scaled: float
    maximum_increment_scaled: float
    minimum_line_search_fraction: float
    no_clipping_or_regularization: bool = True
    energy_admission: bool = False
    temporal_admission: bool = False
    scope: str = 'complete supplied DAE, one backward-Euler nonlinear solve only'

    def receipt(self):
        """Small JSON-safe convergence record; physical state is stored separately."""
        finite=lambda value:float(value) if np.isfinite(value) else None
        return dict(status='DISCRETE_DAE_ROWS_CONVERGED_PENDING_PHYSICAL_ADMISSION',
            iterations=self.iterations,residual_evaluations=self.residual_evaluations,
            jacobian_evaluations=self.jacobian_evaluations,
            initial_residual_max=self.initial_residual_max,
            initial_residual_convention='Scaled F(y_old,0); stationary test, not DAE consistency',
            accepted_residual_max=self.accepted_residual_max,
            minimum_bound_margin_scaled=finite(self.minimum_bound_margin_scaled),
            maximum_increment_scaled=self.maximum_increment_scaled,
            minimum_line_search_fraction=self.minimum_line_search_fraction,
            no_clipping_or_regularization=self.no_clipping_or_regularization,
            energy_admission=self.energy_admission,temporal_admission=self.temporal_admission,
            rank_certified=False,scope=self.scope)


def _vector(value, size, name, *, allow_infinite=False):
    raw=np.asarray(value)
    if raw.shape!=(size,) or np.iscomplexobj(raw) or raw.dtype.kind=='b':
        raise ValueError(name+' must be a real vector of the full state size')
    out=np.array(raw,float,copy=True)
    if np.any(np.isnan(out)) or (not allow_infinite and np.any(~np.isfinite(out))):
        raise ValueError(name+' contains nonfinite entries')
    return out


def _matrix(value,size,name):
    matrix=value.tocsc() if issparse(value) else np.asarray(value)
    if matrix.shape!=(size,size) or np.iscomplexobj(matrix.data if issparse(matrix) else matrix):
        raise ValueError(name+' must be a square real matrix retaining every DAE row')
    if np.any(~np.isfinite(matrix.data if issparse(matrix) else matrix)):
        raise ValueError(name+' contains nonfinite entries')
    return matrix


def _feasible_fraction(x, direction, lower, upper):
    limits=np.full_like(x,np.inf)
    np.divide(upper-x,direction,out=limits,where=direction>0)
    np.divide(lower-x,direction,out=limits,where=direction<0)
    return min(1.,float(np.min(limits,initial=np.inf)))


def _linear_direction(matrix,rhs,tolerance):
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('error',LinAlgWarning)
            warnings.simplefilter('error',MatrixRankWarning)
            step=splu(matrix.tocsc()).solve(rhs) if issparse(matrix) else solve(matrix,rhs,assume_a='gen')
    except (RuntimeError,ValueError,np.linalg.LinAlgError,LinAlgWarning,MatrixRankWarning) as error:
        raise ImplicitStepRejected('Discrete DAE Jacobian is singular/ill-conditioned; no artificial mobility') from error
    error=np.asarray(matrix@step-rhs).ravel()
    if np.any(~np.isfinite(step)) or np.max(abs(error),initial=0)>tolerance*max(1.,float(np.max(abs(rhs),initial=0))):
        raise ImplicitStepRejected('Discrete Newton linear solve fails its independent row residual')
    return step


def solve_backward_euler(previous,dt_s,*,residual,jacobian,differential_mask,
        variable_scales,row_scales,lower_bounds,upper_bounds,initial_guess=None,
        residual_tolerance=1e-9,linear_tolerance=1e-10,max_iterations=30,
        max_line_search=40,boundary_fraction=.995):
    """Solve all supplied physical DAE rows simultaneously without clipping.

    residual(state,derivative) returns one REAL physical row per state entry.
    jacobian(state,derivative) returns (dF/dstate,dF/dderivative), evaluated
    at the same endpoint. Explicit positive SI variable/row scales prevent
    mixing amperes, volts, gap amplitudes and populations in raw norms.
    A gauge must be fixed consistently by the caller before this square solve.
    Derivatives of algebraic variables are identically zero. An initial guess
    must already satisfy support; no hidden interior nudging is performed.

    The method is first order, not a time-accuracy certificate. Any rejected
    solve raises ImplicitStepRejected and returns no modified state.
    """
    old=np.asarray(previous)
    if old.ndim!=1 or not len(old):raise ValueError('One full nonempty previous state required')
    n=len(old);old=_vector(old,n,'previous')
    if not np.isfinite(dt_s) or dt_s<=0:raise ValueError('Positive finite physical dt_s required')
    if not callable(residual) or not callable(jacobian):raise ValueError('Complete residual and Jacobian callables required')
    mask=np.asarray(differential_mask)
    if mask.shape!=(n,) or mask.dtype.kind!='b':raise ValueError('Explicit Boolean differential_mask required')
    scale=_vector(variable_scales,n,'variable_scales');rows=_vector(row_scales,n,'row_scales')
    if np.any(scale<=0) or np.any(rows<=0):raise ValueError('Positive explicit SI scales required')
    low=_vector(lower_bounds,n,'lower_bounds',allow_infinite=True)
    high=_vector(upper_bounds,n,'upper_bounds',allow_infinite=True)
    if np.any(low>=high):raise ValueError('Each state needs distinct ordered bounds')
    if np.any(old<low) or np.any(old>high):raise ValueError('Previous state violates support; no clipping')
    x=old.copy() if initial_guess is None else _vector(initial_guess,n,'initial_guess')
    if np.any(x<low) or np.any(x>high):raise ValueError('Initial guess violates support; no clipping')
    if not 0<boundary_fraction<1:raise ValueError('boundary_fraction must lie strictly between zero and one')
    if (not np.isfinite(residual_tolerance) or residual_tolerance<=0 or
            not np.isfinite(linear_tolerance) or linear_tolerance<=0 or
            type(max_iterations) is not int or max_iterations<1 or
            type(max_line_search) is not int or max_line_search<1):
        raise ValueError('Finite positive numerical limits required')
    evaluations=0;jacobians=0;minimum_fraction=1.
    def rate(z):return mask*(z-old)/dt_s
    def evaluate(z,zdot=None):
        nonlocal evaluations
        evaluations+=1
        return _vector(residual(z,rate(z) if zdot is None else zdot),n,'physical DAE residual')/rows
    initial=float(np.max(abs(evaluate(old,np.zeros(n))),initial=0))
    f=evaluate(x)
    for iteration in range(max_iterations+1):
        error=float(np.max(abs(f),initial=0))
        if error<=residual_tolerance:
            margin=np.minimum((x-low)/scale,(high-x)/scale)
            for value in (x,f):value.setflags(write=False)
            derivative=rate(x);derivative.setflags(write=False)
            return ImplicitStepResult(x,derivative,f,iteration,evaluations,jacobians,
                initial,error,float(np.min(margin,initial=np.inf)),
                float(np.max(abs((x-old)/scale),initial=0)),minimum_fraction)
        if iteration==max_iterations:break
        fy,fdot=jacobian(x,rate(x));jacobians+=1
        fy=_matrix(fy,n,'dF/dstate');fdot=_matrix(fdot,n,'dF/dderivative')
        if issparse(fy) or issparse(fdot):
            from scipy.sparse import csc_matrix
            J=diags(1/rows)@(csc_matrix(fy)+csc_matrix(fdot)@diags(mask.astype(float)/dt_s))@diags(scale)
        else:J=(fy+fdot*mask[None,:]/dt_s)*scale[None,:]/rows[:,None]
        delta=_linear_direction(J,-f,linear_tolerance)*scale
        feasible=_feasible_fraction(x,delta,low,high)
        if feasible<=0:
            raise ImplicitStepRejected('Newton direction leaves exhausted support; no projection or clipping')
        alpha=1. if feasible>=1. else boundary_fraction*feasible
        merit=float(np.dot(f,f));accepted=False
        for _ in range(max_line_search):
            trial=x+alpha*delta
            # Rounding beyond a bound is a rejected trial, never a clip.
            if np.all(trial>=low) and np.all(trial<=high):
                trial_f=evaluate(trial)
                if float(np.dot(trial_f,trial_f))<=(1-1e-4*alpha)*merit:
                    x,f=trial,trial_f;accepted=True;minimum_fraction=min(minimum_fraction,alpha);break
            alpha*=.5
        if not accepted:raise ImplicitStepRejected('No feasible residual-decreasing discrete step')
    raise ImplicitStepRejected('Full DAE residual did not converge at declared tolerance')
