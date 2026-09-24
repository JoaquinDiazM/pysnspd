"""Same thermal Newton/Armijo method, evaluating energy increments stably.

The alpha=0 graph action, stationarity residual, Jacobian, tolerance and
boundary conditions are identical to thermal_spatial_usadel. Only the Armijo
comparison avoids subtracting nearly equal extensive energy values.
"""
import warnings
import numpy as np
from scipy.sparse.linalg import MatrixRankWarning,spsolve
from . import thermal_spatial_usadel as thermal
from .thermal_snapshot import spectral_difference


def solve_frequency(graph,d,epsilon,*,fixed_nodes=None,fixed_u=None,initial_u=None,tol=1e-8,
                    max_iterations=60,max_backtracks=30):
    d=thermal._complex_vector(d,graph.n_nodes,'d');alpha=np.zeros(len(graph.edges))
    if not np.isfinite(epsilon) or epsilon<=0 or not np.isfinite(tol) or tol<=0:
        raise ValueError('Positive finite epsilon and tolerance required')
    if type(max_iterations) is not int or max_iterations<1 or type(max_backtracks) is not int or max_backtracks<1:
        raise ValueError('Positive integer iteration limits required')
    if (fixed_nodes is None)!=(fixed_u is None):raise ValueError('Supply fixed_nodes and fixed_u together')
    fixed=np.array([],int) if fixed_nodes is None else thermal._fixed_indices(fixed_nodes,graph.n_nodes)
    u=np.array(d/epsilon if initial_u is None else thermal._complex_vector(initial_u,graph.n_nodes,'initial_u'),copy=True)
    if len(fixed):u[fixed]=thermal._complex_vector(fixed_u,len(fixed),'fixed_u')
    free=np.ones(graph.n_nodes,bool);free[fixed]=False;components=np.flatnonzero(np.repeat(free,2))
    scale=graph.area_weights*np.maximum(1.,abs(d));backtracks=[]
    for iteration in range(max_iterations+1):
        evaluation,jacobian=thermal.spectral_residual_jacobian(graph,d,epsilon,u,alpha)
        residual=float(np.max(abs(evaluation.residual[free])/scale[free])) if np.any(free) else 0.
        gradient_residual=float(np.max(abs(evaluation.gradient_u[free])/scale[free])) if np.any(free) else 0.
        if residual<=tol:
            return thermal.SpectralSolution(float(epsilon),u.copy(),evaluation.f,evaluation.g,evaluation.energy,
                residual,gradient_residual,iteration,tuple(backtracks),fixed.copy(),thermal._signature(graph,d,alpha))
        if iteration==max_iterations:raise RuntimeError(f'Stable-energy Newton iteration limit; residual={residual:.6g}')
        rhs=np.column_stack((evaluation.residual.real,evaluation.residual.imag)).ravel()
        with warnings.catch_warnings():
            warnings.simplefilter('error',MatrixRankWarning)
            direction_free=spsolve(jacobian[components][:,components],-rhs[components])
        if np.any(~np.isfinite(direction_free)):raise RuntimeError('Nonfinite spectral Newton direction')
        direction=np.zeros(2*graph.n_nodes);direction[components]=direction_free
        dz=direction[::2]+1j*direction[1::2];slope=float(np.real(np.vdot(evaluation.gradient_u,dz)))
        if not np.isfinite(slope) or slope>=0:raise RuntimeError(f'Newton direction is not energy descent; slope={slope:.6g}')
        for reduction in range(max_backtracks):
            step=2.**(-reduction);trial=u+step*dz
            difference=spectral_difference(graph,d,u,d,trial,epsilon)['renormalized_energy_difference']
            if difference<=1e-4*step*slope:
                u=trial;backtracks.append(reduction);break
        else:raise RuntimeError('Stable-energy Newton Armijo failed; no fallback or tolerance relaxation')
    raise AssertionError('Unreachable Newton state')
