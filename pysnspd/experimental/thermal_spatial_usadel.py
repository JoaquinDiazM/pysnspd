"""Experimental graph discretization of the thermal spatial Usadel energy.

Units: d=Delta/(kB Tc), epsilon=omega/(kB Tc), X=x/ell0 and
ell0=sqrt(hbar D/(2 kB Tc)). The area weights integrate dX dY. This is a
finite positive-Matsubara sum, with no tail, kinetic evolution or K0 addition.
U_ij=exp(-i alpha_ij) transports the head spectral field to the tail.
The physical oriented current is minus the derivative with respect to alpha.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import warnings

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import MatrixRankWarning, spsolve


@dataclass(frozen=True)
class ThermalGraph:
    area_weights: np.ndarray
    edges: np.ndarray
    conductance: np.ndarray
    coordinates_bar: np.ndarray | None = None
    boundary_nodes: np.ndarray | None = None

    def __post_init__(self):
        mass = np.array(self.area_weights, dtype=float, copy=True)
        raw_edges = np.asarray(self.edges)
        edges = np.array(raw_edges, dtype=int, copy=True)
        conductance = np.array(self.conductance, dtype=float, copy=True)
        if mass.ndim != 1 or not len(mass) or np.any(~np.isfinite(mass)) or np.any(mass <= 0):
            raise ValueError('Positive finite node area weights required')
        if (raw_edges.dtype.kind not in 'iu' or edges.ndim != 2 or edges.shape[1] != 2
                or np.any(edges < 0) or np.any(edges >= len(mass)) or np.any(edges[:, 0] == edges[:, 1])):
            raise ValueError('Oriented graph edges must connect distinct valid nodes')
        if conductance.shape != (len(edges),) or np.any(~np.isfinite(conductance)) or np.any(conductance < 0):
            raise ValueError('Finite nonnegative edge conductance required; zero edges are inactive')
        for name, value in (('area_weights', mass), ('edges', edges), ('conductance', conductance)):
            value.setflags(write=False)
            object.__setattr__(self, name, value)
        if self.coordinates_bar is not None:
            coordinates = np.array(self.coordinates_bar, dtype=float, copy=True)
            if coordinates.shape != (len(mass), 2) or np.any(~np.isfinite(coordinates)):
                raise ValueError('One finite planar coordinate per node required')
            coordinates.setflags(write=False)
            object.__setattr__(self, 'coordinates_bar', coordinates)
        if self.boundary_nodes is not None:
            nodes = _fixed_indices(self.boundary_nodes, len(mass))
            nodes.setflags(write=False)
            object.__setattr__(self, 'boundary_nodes', nodes)

    @property
    def n_nodes(self):
        return len(self.area_weights)

    @property
    def signature(self):
        digest = hashlib.sha256()
        for array in (self.area_weights, self.edges, self.conductance):
            digest.update(str(array.shape).encode())
            digest.update(array.tobytes())
        return digest.hexdigest()


def rectangular_graph(nx, ny, length_bar, width_bar):
    """Vertex-dual orthogonal rectangle, including half boundary measures."""
    if type(nx) is not int or type(ny) is not int or min(nx, ny) < 2:
        raise ValueError('At least two nodes per Cartesian direction required')
    if not np.isfinite(length_bar) or not np.isfinite(width_bar) or min(length_bar, width_bar) <= 0:
        raise ValueError('Positive finite rectangle dimensions required')
    x, y = np.linspace(0., length_bar, nx), np.linspace(-width_bar/2, width_bar/2, ny)
    hx, hy = length_bar/(nx-1), width_bar/(ny-1)
    wx, wy = np.full(nx, hx), np.full(ny, hy)
    wx[[0, -1]] /= 2
    wy[[0, -1]] /= 2
    ids = np.arange(nx*ny).reshape(nx, ny)
    edges, weights = [], []
    for i in range(nx-1):
        for j in range(ny):
            edges.append((ids[i, j], ids[i+1, j]))
            weights.append(wy[j]/hx)
    for i in range(nx):
        for j in range(ny-1):
            edges.append((ids[i, j], ids[i, j+1]))
            weights.append(wx[i]/hy)
    boundary = np.unique(np.r_[ids[0], ids[-1], ids[:, 0], ids[:, -1]])
    return ThermalGraph(np.kron(wx, wy), np.asarray(edges), np.asarray(weights),
        np.column_stack((np.repeat(x, ny), np.tile(y, nx))), boundary)


def _complex_vector(value, size, name):
    raw = np.asarray(value)
    if raw.shape == (size, 2) and not np.iscomplexobj(raw):
        result = raw[:, 0]+1j*raw[:, 1]
    else:
        result = np.asarray(value, dtype=complex)
    if result.shape != (size,) or np.any(~np.isfinite(result)):
        raise ValueError(name+' must contain one finite complex value per node')
    return result


def _alpha(graph, alpha):
    raw = np.zeros(len(graph.edges)) if alpha is None else np.asarray(alpha)
    if np.iscomplexobj(raw) or raw.shape != (len(graph.edges),) or np.any(~np.isfinite(raw)):
        raise ValueError('One finite real alpha per edge required')
    return np.asarray(raw, dtype=float)


def _fixed_indices(value, size):
    raw = np.asarray(value)
    if raw.dtype.kind == 'b' and raw.shape == (size,):
        nodes = np.flatnonzero(raw)
    elif raw.dtype.kind in 'iu' and raw.ndim == 1:
        nodes = np.asarray(raw, dtype=int)
    else:
        raise ValueError('Fixed nodes must be a Boolean mask or integer indices')
    if np.any(nodes < 0) or np.any(nodes >= size) or len(np.unique(nodes)) != len(nodes):
        raise ValueError('Fixed nodes must be unique valid indices')
    return nodes


def _fields(u):
    # hypot avoids overflow of |u|^2 when diagnosing an unsuccessful iterate.
    g = 1/np.hypot(1., abs(u))
    return g*u, g


@dataclass(frozen=True)
class SpectralEvaluation:
    energy: float
    gradient_u: np.ndarray
    residual: np.ndarray
    f: np.ndarray
    g: np.ndarray
    link_derivative_alpha: np.ndarray


def spectral_energy_gradient(graph, d, epsilon, u, alpha=None):
    """One spectral energy, its full real gradient in complex notation and R.

    dE=Re(vdot(gradient_u,du)); the Matsubara prefactor 2*pi*t and the
    d-dependent renormalization counterterm are not part of this one mode.
    R_i=Hz_i*u_i-Hxy_i has the same stationary points as gradient_u.
    """
    d = _complex_vector(d, graph.n_nodes, 'd')
    u = _complex_vector(u, graph.n_nodes, 'u')
    if not np.isfinite(epsilon) or epsilon <= 0:
        raise ValueError('Positive finite Matsubara epsilon required')
    phase = np.exp(-1j*_alpha(graph, alpha))
    f, g = _fields(u)
    tail, head = graph.edges.T
    c, m = graph.conductance, graph.area_weights
    difference = f[tail]-phase*f[head]
    one_minus_g = abs(f)**2/(1+g)
    energy = np.dot(m, 2*epsilon*one_minus_g-2*np.real(np.conj(d)*f))
    energy += np.dot(c, abs(difference)**2+(g[tail]-g[head])**2)
    hz, hxy = m*epsilon, m*d
    np.add.at(hz, tail, c*g[head])
    np.add.at(hz, head, c*g[tail])
    np.add.at(hxy, tail, c*phase*f[head])
    np.add.at(hxy, head, c*np.conj(phase)*f[tail])
    residual = hz*u-hxy
    projection = np.real(np.conj(u)*residual)
    gradient = 2*g*(residual-g*g*u*projection)
    derivative = -2*c*np.imag(np.conj(f[tail])*phase*f[head])
    return SpectralEvaluation(float(energy), gradient, residual, f, g, derivative)


def spectral_residual_jacobian(graph, d, epsilon, u, alpha=None):
    """Exact sparse Jacobian of R in interleaved (Re u_i, Im u_i) order."""
    d = _complex_vector(d, graph.n_nodes, 'd')
    u = _complex_vector(u, graph.n_nodes, 'u')
    evaluated = spectral_energy_gradient(graph, d, epsilon, u, alpha)
    phase = np.exp(-1j*_alpha(graph, alpha))
    xy = np.column_stack((u.real, u.imag))
    hz = graph.area_weights*epsilon
    tail, head = graph.edges.T
    np.add.at(hz, tail, graph.conductance*evaluated.g[head])
    np.add.at(hz, head, graph.conductance*evaluated.g[tail])
    rows = list(np.arange(2*graph.n_nodes))
    cols = rows.copy()
    values = list(np.repeat(hz, 2))
    for i, j, weight, transport in zip(tail, head, graph.conductance, phase):
        if weight == 0:
            continue
        for target, source, rotation in ((i, j, transport), (j, i, np.conj(transport))):
            rot = np.array([[rotation.real, -rotation.imag], [rotation.imag, rotation.real]])
            gs = evaluated.g[source]
            block = weight*(-gs*rot+gs**3*np.outer(rot@xy[source]-xy[target], xy[source]))
            for a in range(2):
                for b in range(2):
                    rows.append(2*target+a)
                    cols.append(2*source+b)
                    values.append(block[a, b])
    jacobian = coo_matrix((values, (rows, cols)), shape=(2*graph.n_nodes,)*2).tocsr()
    return evaluated, jacobian


def _signature(graph, d, alpha):
    digest = hashlib.sha256(graph.signature.encode())
    digest.update(np.asarray(d, dtype=complex).tobytes())
    digest.update(np.asarray(alpha, dtype=float).tobytes())
    return digest.hexdigest()


@dataclass(frozen=True)
class SpectralSolution:
    epsilon: float
    u: np.ndarray
    f: np.ndarray
    g: np.ndarray
    energy: float
    residual: float
    gradient_residual: float
    iterations: int
    line_search_steps: tuple[int, ...]
    fixed_nodes: np.ndarray
    source_signature: str
    converged: bool = True


def solve_frequency(graph, d, epsilon, alpha=None, *, fixed_nodes=None,
                    fixed_u=None, initial_u=None, tol=1e-8,
                    max_iterations=60, max_backtracks=30):
    """Sparse damped Newton; errors stop the solve without another algorithm.

    fixed_u are imposed spectral boundary values, independent of a varied d.
    Free boundaries are the natural zero-flux boundary of the graph energy.
    A returned solution is spectrally stationary for the prescribed d; the
    condensate is not assumed self-consistent or dynamically admitted.
    """
    d = _complex_vector(d, graph.n_nodes, 'd')
    alpha = _alpha(graph, alpha)
    if not np.isfinite(epsilon) or epsilon <= 0 or not np.isfinite(tol) or tol <= 0:
        raise ValueError('Positive finite epsilon and tolerance required')
    if type(max_iterations) is not int or max_iterations < 1 or type(max_backtracks) is not int or max_backtracks < 1:
        raise ValueError('Positive integer iteration limits required')
    if (fixed_nodes is None) != (fixed_u is None):
        raise ValueError('Supply fixed_nodes and fixed_u together')
    fixed = np.array([], dtype=int) if fixed_nodes is None else _fixed_indices(fixed_nodes, graph.n_nodes)
    u = np.array(d/epsilon if initial_u is None else _complex_vector(initial_u, graph.n_nodes, 'initial_u'), copy=True)
    if len(fixed):
        u[fixed] = _complex_vector(fixed_u, len(fixed), 'fixed_u')
    free = np.ones(graph.n_nodes, dtype=bool)
    free[fixed] = False
    components = np.flatnonzero(np.repeat(free, 2))
    scale = graph.area_weights*np.maximum(1., abs(d))
    backtracks = []
    for iteration in range(max_iterations+1):
        evaluation, jacobian = spectral_residual_jacobian(graph, d, epsilon, u, alpha)
        residual = float(np.max(abs(evaluation.residual[free])/scale[free])) if np.any(free) else 0.
        grad_residual = float(np.max(abs(evaluation.gradient_u[free])/scale[free])) if np.any(free) else 0.
        if residual <= tol:
            return SpectralSolution(float(epsilon), u.copy(), evaluation.f, evaluation.g,
                evaluation.energy, residual, grad_residual, iteration, tuple(backtracks),
                fixed.copy(), _signature(graph, d, alpha))
        if iteration == max_iterations:
            raise RuntimeError(f'Spectral Newton iteration limit; residual={residual:.6g}')
        rhs = np.column_stack((evaluation.residual.real, evaluation.residual.imag)).ravel()
        with warnings.catch_warnings():
            warnings.simplefilter('error', MatrixRankWarning)
            direction_free = spsolve(jacobian[components][:, components], -rhs[components])
        if np.any(~np.isfinite(direction_free)):
            raise RuntimeError('Nonfinite spectral Newton direction')
        direction = np.zeros(2*graph.n_nodes)
        direction[components] = direction_free
        dz = direction[::2]+1j*direction[1::2]
        slope = float(np.real(np.vdot(evaluation.gradient_u, dz)))
        if not np.isfinite(slope) or slope >= 0:
            raise RuntimeError(f'Spectral Newton direction is not energy-descent; slope={slope:.6g}')
        for reduction in range(max_backtracks):
            step = 2.**(-reduction)
            trial = u+step*dz
            trial_evaluation = spectral_energy_gradient(graph, d, epsilon, trial, alpha)
            if trial_evaluation.energy <= evaluation.energy+1e-4*step*slope:
                u = trial
                backtracks.append(reduction)
                break
        else:
            raise RuntimeError('Spectral Newton energy line search failed; no fallback')
    raise AssertionError('Unreachable spectral iteration state')


@dataclass(frozen=True)
class ThermalEvaluation:
    energy: float
    gap_gradient: np.ndarray
    link_derivative_alpha: np.ndarray
    current_bar: np.ndarray
    noether_residual: np.ndarray
    modes: int
    maximum_spectral_residual: float
    scope: str


def evaluate_thermal(graph, d, t, solutions, alpha=None):
    """On-shell finite-sum energy/force/current, with fixed spectral boundaries.

    gap_gradient is the integrated real Cartesian gradient represented as a
    complex vector. Divide by area_weights for its volume density. current_bar
    is oriented tail-to-head and equals -dF/dalpha. SI current needs the common
    physical energy factor times 2e/hbar. There is no separate K0 term.
    """
    d = _complex_vector(d, graph.n_nodes, 'd')
    alpha = _alpha(graph, alpha)
    if not np.isfinite(t) or t <= 0:
        raise ValueError('Positive finite T/Tc required')
    solutions = tuple(solutions)
    if not solutions:
        raise ValueError('At least one positive Matsubara solution required')
    signature = _signature(graph, d, alpha)
    energy = float(np.dot(graph.area_weights, abs(d)**2)*np.log(t))
    force = 2*graph.area_weights*d*np.log(t)
    derivative = np.zeros(len(graph.edges))
    max_residual = 0.
    spectral_boundary_torque = np.zeros(graph.n_nodes)
    for n, solution in enumerate(solutions):
        epsilon = 2*np.pi*t*(n+.5)
        if (not solution.converged or solution.source_signature != signature
                or not np.isclose(solution.epsilon, epsilon, rtol=1e-13, atol=0)):
            raise ValueError('Solutions must belong to these fields/links/graph and consecutive Matsubara frequencies')
        value = spectral_energy_gradient(graph, d, epsilon, solution.u, alpha)
        energy += 2*np.pi*t*(value.energy+np.dot(graph.area_weights, abs(d)**2)/epsilon)
        force += 4*np.pi*t*graph.area_weights*(d/epsilon-value.f)
        derivative += 2*np.pi*t*value.link_derivative_alpha
        max_residual = max(max_residual, solution.residual)
        spectral_boundary_torque += 2*np.pi*t*np.imag(np.conj(solution.u)*value.gradient_u)
    current = -derivative
    divergence = np.zeros(graph.n_nodes)
    np.add.at(divergence, graph.edges[:, 0], current)
    np.add.at(divergence, graph.edges[:, 1], -current)
    # For fixed spectral boundary values their conjugate phase work remains in
    # the identity; the interior spectral residual vanishes at stationarity.
    noether = np.imag(np.conj(d)*force)+divergence+spectral_boundary_torque
    return ThermalEvaluation(float(energy), force, derivative, current, noether,
        len(solutions), max_residual,
        'Thermal finite Matsubara sum; imposed d and spectral boundaries; no tail, KWT, kinetic evolution or physical-core admission')
