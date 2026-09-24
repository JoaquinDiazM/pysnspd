"""Analytic continuation of the finite thermal Usadel graph action.

The two stereographic fields a,b are independent complex unknowns. Conjugating
b during a retarded solve is incorrect: b=a.conj() only on the real positive
Matsubara axis. z=eta-i*E with eta>0 specifies the retarded half-plane in the
Matsubara convention. This module solves spectral roots, never minimizes the
real part of a complex action. It contains no distribution/kinetic closure.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import MatrixRankWarning, spsolve
import warnings

from .thermal_spatial_usadel import ThermalGraph, _alpha, _complex_vector, _fixed_indices


def fields(a, b):
    a, b = np.asarray(a, complex), np.asarray(b, complex)
    if a.shape != b.shape or np.any(~np.isfinite(a)) or np.any(~np.isfinite(b)):
        raise ValueError('Two matching finite complex stereographic fields required')
    denominator = 1+a*b
    if np.any(abs(denominator) < 1e-14):
        raise FloatingPointError('Stereographic chart approached a pole; no clipping')
    return 2*a/denominator, 2*b/denominator, (1-a*b)/denominator


def uniform_fields(d, z):
    """Exact causal homogeneous BCS anchor, with positive-real square root."""
    d = np.asarray(d, complex)
    z = complex(z)
    if not np.isfinite(z) or z.real <= 0 or np.any(~np.isfinite(d)):
        raise ValueError('Finite d and z in the open right half-plane required')
    root = np.sqrt(z*z+abs(d)**2)
    root = np.where(root.real < 0, -root, root)
    return d/(z+root), np.conj(d)/(z+root)


@dataclass(frozen=True)
class RetardedEvaluation:
    action: complex
    residual: np.ndarray
    f: np.ndarray
    f_tilde: np.ndarray
    g: np.ndarray
    link_derivative_alpha: np.ndarray
    normalization_residual: float


def evaluate(graph: ThermalGraph, d, z, a, b, alpha=None):
    """Complex analytic action and its stationary root equations.

    The d-dependent counterterm is omitted: it does not affect spectral roots.
    d is a physical prescribed gap and its conjugate is a fixed coefficient;
    a,b remain independent in every expression and derivative.
    """
    size = graph.n_nodes
    d = _complex_vector(d, size, 'd')
    a, b = _complex_vector(a, size, 'a'), _complex_vector(b, size, 'b')
    z = complex(z)
    if not np.isfinite(z) or z.real <= 0:
        raise ValueError('z must lie in the open right half-plane')
    phase = np.exp(-1j*_alpha(graph, alpha))
    f, ft, g = fields(a, b)
    tail, head = graph.edges.T
    c, m = graph.conductance, graph.area_weights
    df, dft = f[tail]-phase*f[head], ft[tail]-np.conj(phase)*ft[head]
    action = np.dot(m, 2*z*(1-g)-np.conj(d)*f-d*ft)
    action += np.dot(c, df*dft+(g[tail]-g[head])**2)
    hz, hp, hm = m*z, m*d, m*np.conj(d)
    np.add.at(hz, tail, c*g[head]); np.add.at(hz, head, c*g[tail])
    np.add.at(hp, tail, c*phase*f[head]); np.add.at(hp, head, c*np.conj(phase)*f[tail])
    np.add.at(hm, tail, c*np.conj(phase)*ft[head]); np.add.at(hm, head, c*phase*ft[tail])
    residual = np.column_stack((hm*a*a+2*hz*a-hp, hp*b*b+2*hz*b-hm))
    derivative = 1j*c*(phase*f[head]*ft[tail]-np.conj(phase)*ft[head]*f[tail])
    return RetardedEvaluation(complex(action), residual, f, ft, g, derivative,
        float(np.max(abs(g*g+f*ft-1))))


def residual_jacobian(graph, d, z, a, b, alpha=None):
    """Exact sparse complex Jacobian; no numerical differencing or conjugation."""
    size = graph.n_nodes
    a, b = _complex_vector(a, size, 'a'), _complex_vector(b, size, 'b')
    observation = evaluate(graph, d, z, a, b, alpha)
    phase = np.exp(-1j*_alpha(graph, alpha))
    tail, head = graph.edges.T
    c, m = graph.conductance, graph.area_weights
    hz = m*complex(z)
    hp, hm = m*np.asarray(d, complex), m*np.conj(d)
    np.add.at(hz, tail, c*observation.g[head]); np.add.at(hz, head, c*observation.g[tail])
    np.add.at(hp, tail, c*phase*observation.f[head]); np.add.at(hp, head, c*np.conj(phase)*observation.f[tail])
    np.add.at(hm, tail, c*np.conj(phase)*observation.f_tilde[head]); np.add.at(hm, head, c*phase*observation.f_tilde[tail])
    nodes = np.arange(size)
    rows = [2*nodes, 2*nodes+1]
    cols = [2*nodes, 2*nodes+1]
    values = [2*hm*a+2*hz, 2*hp*b+2*hz]
    denominator2 = (1+a*b)**2
    fa, fb = 2/denominator2, -2*a*a/denominator2
    fta, ftb = -2*b*b/denominator2, 2/denominator2
    ga, gb = -2*b/denominator2, -2*a/denominator2
    for i, j, transport in ((tail, head, phase), (head, tail, np.conj(phase))):
        for component, ff, ftf, gg in ((0, fa, fta, ga), (1, fb, ftb, gb)):
            rows += [2*i, 2*i+1]
            cols += [2*j+component, 2*j+component]
            values += [c*(a[i]**2*np.conj(transport)*ftf[j]+2*a[i]*gg[j]-transport*ff[j]),
                c*(b[i]**2*transport*ff[j]+2*b[i]*gg[j]-np.conj(transport)*ftf[j])]
    jacobian = coo_matrix((np.concatenate(values), (np.concatenate(rows), np.concatenate(cols))),
        shape=(2*size, 2*size)).tocsr()
    return observation, jacobian


@dataclass(frozen=True)
class RetardedSolution:
    z: complex
    a: np.ndarray
    b: np.ndarray
    f: np.ndarray
    f_tilde: np.ndarray
    g: np.ndarray
    residual: float
    normalization_residual: float
    iterations: int
    backtracks: tuple[int, ...]
    minimum_DOS: float


def solve(graph, d, z, alpha=None, *, initial_a=None, initial_b=None,
          fixed_nodes=None, fixed_a=None, fixed_b=None, tolerance=1e-9,
          max_iterations=60, max_backtracks=24, causal_tolerance=1e-8):
    """Damped complex Newton root solve; branch selection needs continuation.

    Every accepted Newton update decreases the scaled root residual. Such
    numerical updates are not physical time steps or energy minimization.
    Positive DOS is checked, not clipped, at each returned stationary solution.
    """
    size = graph.n_nodes
    d = _complex_vector(d, size, 'd')
    if not np.isfinite(tolerance) or tolerance <= 0 or not np.isfinite(causal_tolerance) or causal_tolerance < 0:
        raise ValueError('Positive finite tolerance and nonnegative causal tolerance required')
    if type(max_iterations) is not int or max_iterations < 1 or type(max_backtracks) is not int or max_backtracks < 1:
        raise ValueError('Positive integer solver budgets required')
    if (initial_a is None) != (initial_b is None):
        raise ValueError('Both initial stereographic fields are required together')
    if initial_a is None:
        a, b = uniform_fields(d, z)
    else:
        a, b = _complex_vector(initial_a, size, 'a').copy(), _complex_vector(initial_b, size, 'b').copy()
    if (fixed_nodes is None) != (fixed_a is None) or (fixed_nodes is None) != (fixed_b is None):
        raise ValueError('Fixed nodes and both fixed stereographic fields are required together')
    fixed = np.array([], int) if fixed_nodes is None else _fixed_indices(fixed_nodes, size)
    if len(fixed):
        a[fixed], b[fixed] = _complex_vector(fixed_a, len(fixed), 'fixed_a'), _complex_vector(fixed_b, len(fixed), 'fixed_b')
    free = np.ones(size, bool); free[fixed] = False
    components = np.flatnonzero(np.repeat(free, 2))
    scale = graph.area_weights*np.maximum(1., abs(d))
    reductions = []
    for iteration in range(max_iterations+1):
        observation, jacobian = residual_jacobian(graph, d, z, a, b, alpha)
        scaled = observation.residual[free]/scale[free, None]
        residual = float(np.max(abs(scaled))) if len(scaled) else 0.
        if residual <= tolerance:
            minimum = float(np.min(observation.g.real))
            if minimum < -causal_tolerance:
                raise RuntimeError(f'Negative DOS on stationary spectral branch: {minimum}; no clipping')
            return RetardedSolution(complex(z), a.copy(), b.copy(), observation.f,
                observation.f_tilde, observation.g, residual, observation.normalization_residual,
                iteration, tuple(reductions), minimum)
        if iteration == max_iterations:
            raise RuntimeError(f'Retarded Newton iteration budget exhausted; residual={residual:.6g}')
        with warnings.catch_warnings():
            warnings.simplefilter('error', MatrixRankWarning)
            update = spsolve(jacobian[components][:, components], -observation.residual.ravel()[components])
        if np.any(~np.isfinite(update)):
            raise RuntimeError('Nonfinite retarded Newton update')
        direction = np.zeros((size, 2), complex)
        direction.ravel()[components] = update
        merit = float(np.vdot(scaled.ravel(), scaled.ravel()).real)
        for reduction in range(max_backtracks):
            fraction = 2.**(-reduction)
            aa, bb = a+fraction*direction[:, 0], b+fraction*direction[:, 1]
            try:
                trial = evaluate(graph, d, z, aa, bb, alpha).residual[free]/scale[free, None]
            except FloatingPointError:
                # An invalid trial coordinate is rejected, never clipped. A
                # shorter step follows the same Newton direction and action.
                continue
            trial_merit = float(np.vdot(trial.ravel(), trial.ravel()).real)
            if trial_merit <= (1-1e-4*fraction)*merit:
                a, b = aa, bb
                reductions.append(reduction)
                break
        else:
            raise RuntimeError('Retarded root merit line search failed; no alternate root or clipping')
    raise AssertionError('Unreachable retarded root state')
