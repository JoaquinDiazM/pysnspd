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
from scipy.sparse.linalg import MatrixRankWarning, spsolve, splu
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
    d = _complex_vector(d, size, 'd')
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


def _direction_vector(value, size, name, *, real=False):
    """A scalar zero/direction or one finite direction per node/edge."""
    raw = np.asarray(value)
    if real and np.iscomplexobj(raw):
        raise ValueError(name+' must be real')
    result = np.asarray(value, float if real else complex)
    if result.ndim == 0:
        result = np.full(size, result, dtype=result.dtype)
    if result.shape != (size,) or np.any(~np.isfinite(result)):
        raise ValueError(name+' must be a finite scalar or one direction per entry')
    return result


def nambu_matrix(f, f_tilde, g):
    """Traceless R=[[g,f],[f_tilde,-g]], without conjugating f_tilde."""
    f, ft, g = np.broadcast_arrays(f, f_tilde, g)
    result = np.empty(f.shape+(2, 2), complex)
    result[..., 0, 0], result[..., 0, 1] = g, f
    result[..., 1, 0], result[..., 1, 1] = ft, -g
    return result


def field_direction(a, b, da, db):
    """Exact Cartesian derivative of the stereographic fields and R0."""
    a, b, da, db = np.broadcast_arrays(
        np.asarray(a, complex), np.asarray(b, complex), da, db)
    if not all(np.all(np.isfinite(x)) for x in (a, b, da, db)):
        raise ValueError('Finite stereographic fields and directions required')
    denominator = 1+a*b
    if np.any(abs(denominator) < 1e-14):
        raise FloatingPointError('Stereographic chart approached a pole; no clipping')
    product = b*da+a*db
    df = 2*da/denominator-2*a*product/denominator**2
    dft = 2*db/denominator-2*b*product/denominator**2
    dg = -2*product/denominator**2
    return df, dft, dg, nambu_matrix(df, dft, dg)


def current_carrying_cut_tangent(theta, z, gap_bar, q_bar, phase, *,
                                 gap_direction=0., q_direction=0.,
                                 phase_direction=0., z_direction=0.):
    """Derivative of the SAME homogeneous current-carrying cut law.

    The supplied causal theta satisfies
      z*sin(theta)-gap_bar*cos(theta)+q_bar**2*sin(theta)*cos(theta)=0.
    phase is a real angle in radians, not exp(i*phase). Returned (da,db)
    preserve the independent retarded partners a=tan(theta/2)*exp(i*phase),
    b=tan(theta/2)*exp(-i*phase). A physical gap-magnitude/q/phase time
    direction is supplied per second; the result is then per second.
    For d/d(E/Eu), use z_direction=-1j with fixed numerical eta. This
    function does not select a root, alter eta or determine a boundary law.
    """
    theta, z = np.asarray(theta, complex), np.asarray(z, complex)
    real_inputs = (gap_bar, q_bar, phase, gap_direction, q_direction, phase_direction)
    if any(np.iscomplexobj(np.asarray(x)) for x in real_inputs):
        raise ValueError('Gap magnitude, q, phase and their directions must be real')
    theta, z, d, q, phase, dd, dq, dp, dz = np.broadcast_arrays(
        theta, z, *[np.asarray(x, float) for x in real_inputs],
        np.asarray(z_direction, complex))
    if not all(np.all(np.isfinite(x)) for x in (theta, z, d, q, phase, dd, dq, dp, dz)):
        raise ValueError('Finite current-carrying cut data required')
    if np.any(d < 0) or np.any(z.real <= 0):
        raise ValueError('Nonnegative gap and causal open-half-plane z required')
    sine, cosine = np.sin(theta), np.cos(theta)
    equation = z*sine-d*cosine+q*q*sine*cosine
    scale = np.maximum(1., abs(z)+d+q*q)
    if np.any(abs(equation) > 1e-8*scale):
        raise ValueError('Supplied theta does not solve the declared current-carrying cut law')
    denominator = z*cosine+d*sine+q*q*np.cos(2*theta)
    if np.any(abs(denominator) <= 32*np.finfo(float).eps*scale):
        raise FloatingPointError('Singular current-carrying cut tangent; no regularization')
    direction = (dd*cosine-dz*sine-2*q*dq*sine*cosine)/denominator
    chart = np.tan(theta/2)
    chart_direction = .5*direction/(np.cos(theta/2)**2)
    da = np.exp(1j*phase)*(chart_direction+1j*chart*dp)
    db = np.exp(-1j*phase)*(chart_direction-1j*chart*dp)
    if np.any(~np.isfinite(da)) or np.any(~np.isfinite(db)):
        raise FloatingPointError('Nonfinite current-carrying cut tangent')
    return da, db


@dataclass(frozen=True)
class RetardedDirection:
    """Full-node spectral derivative, including prescribed-cut reactions.

    Units follow the supplied direction: energy derivative is per E/Eu,
    time derivative is per second; residual_derivative has graph residual
    units times those same derivative units. No boundary reaction is dropped.
    """
    da: np.ndarray
    db: np.ndarray
    df: np.ndarray
    df_tilde: np.ndarray
    dg: np.ndarray
    matrix: np.ndarray
    residual_derivative: np.ndarray
    fixed_nodes: np.ndarray
    boundary_reaction: np.ndarray
    free_equation_residual_max: float
    free_equation_residual_relative: float
    fixed_law_residual_max: float
    normalization_derivative_max: float


class RetardedTangent:
    """Reuse one exact sparse Jacobian for moving gap/link/cut directions.

    Solve J_ff dx_f=-F_parameter,f-J_fb dx_b with explicitly prescribed
    cut dx_b. Default zero cut directions retain the old fixed-cut controls;
    they do not assert that physical superconducting cuts remain stationary.
    Gap directions are physical complex Cartesian variations, so d* varies
    as conjugate(dd); a,b remain independent complex spectral coordinates.
    """
    def __init__(self, graph, d, solution, alpha=None, *, fixed_nodes=(),
                 spectral_tolerance=1e-7):
        self.graph = graph
        self.d = _complex_vector(d, graph.n_nodes, 'd').copy()
        self.alpha = _alpha(graph, alpha).copy()
        self.z = complex(solution.z)
        self.a = _complex_vector(solution.a, graph.n_nodes, 'a').copy()
        self.b = _complex_vector(solution.b, graph.n_nodes, 'b').copy()
        if not np.isfinite(spectral_tolerance) or spectral_tolerance <= 0:
            raise ValueError('Positive finite spectral tolerance required')
        raw_fixed = np.asarray(fixed_nodes)
        self.fixed_nodes = (np.empty(0, int) if raw_fixed.size == 0
                            else _fixed_indices(fixed_nodes, graph.n_nodes))
        self.free = np.ones(graph.n_nodes, bool)
        self.free[self.fixed_nodes] = False
        self.components = np.flatnonzero(np.repeat(self.free, 2))
        self.fixed_components = np.flatnonzero(np.repeat(~self.free, 2))
        self.value, self._jacobian = residual_jacobian(
            graph, self.d, self.z, self.a, self.b, self.alpha)
        for storage in (self._jacobian.data, self._jacobian.indices, self._jacobian.indptr):
            storage.setflags(write=False)
        scale = graph.area_weights*np.maximum(1., abs(self.d))
        self.reference_spectral_residual = float(np.max(
            abs(self.value.residual[self.free])/scale[self.free, None], initial=0.))
        if self.reference_spectral_residual > spectral_tolerance:
            raise ValueError('Supplied spectrum is not stationary on the free nodes')
        self.matrix = nambu_matrix(self.value.f, self.value.f_tilde, self.value.g)
        self.factor = (splu(self.jacobian[self.components][:, self.components].tocsc())
                       if len(self.components) else None)

    @property
    def jacobian(self):
        """Read-only stationary CSR Jacobian; use .copy() for assembly edits."""
        return self._jacobian

    def _fixed_coordinates(self, fixed_a_direction, fixed_b_direction,
                           fixed_matrix_direction):
        count = len(self.fixed_nodes)
        if fixed_matrix_direction is not None:
            if fixed_a_direction is not None or fixed_b_direction is not None:
                raise ValueError('Supply fixed matrix OR both stereographic directions')
            raw = np.asarray(fixed_matrix_direction, complex)
            if raw.shape != (count, 2, 2) or np.any(~np.isfinite(raw)):
                raise ValueError('One finite 2x2 R0 derivative per prescribed node required')
            fixed = self.fixed_nodes
            f, ft, g = self.value.f[fixed], self.value.f_tilde[fixed], self.value.g[fixed]
            denominator = 1+g
            if np.any(abs(denominator) < 1e-14):
                raise FloatingPointError('Prescribed derivative approached chart pole; no clipping')
            dg, df, dft = raw[:, 0, 0], raw[:, 0, 1], raw[:, 1, 0]
            da = (df*denominator-f*dg)/denominator**2
            db = (dft*denominator-ft*dg)/denominator**2
            reconstructed = field_direction(self.a[fixed], self.b[fixed], da, db)[3]
            scale = np.maximum(1., np.max(abs(raw), axis=(-2, -1)))
            if np.any(np.max(abs(reconstructed-raw), axis=(-2, -1)) > 1e-9*scale):
                raise ValueError('Prescribed R0 derivative violates traceless normalization tangent')
            return da, db
        if (fixed_a_direction is None) != (fixed_b_direction is None):
            raise ValueError('Both prescribed stereographic directions are required together')
        if fixed_a_direction is None:
            return np.zeros(count, complex), np.zeros(count, complex)
        da = _direction_vector(fixed_a_direction, count, 'fixed_a_direction')
        db = _direction_vector(fixed_b_direction, count, 'fixed_b_direction')
        return da, db

    def apply(self, gap_direction=0., *, z_direction=0., alpha_direction=0.,
              fixed_a_direction=None, fixed_b_direction=None,
              fixed_matrix_direction=None):
        """Exact tangent of the same stationary graph, with moving cuts.

        Use z_direction=-1j for RE=dR0/d(E/Eu). Use gap_direction=d(Delta/Eu)/dt
        and alpha_direction=dalpha/dt in s^-1 for Rt. No term from a constant
        numerical eta is differentiated into a physical relaxation or source.
        fixed_matrix_direction supplies a normalized traceless R0 derivative;
        alternatively prescribe both da and db. Values refer to fixed_nodes
        in their original order. Boundary residuals are explicit reactions.
        """
        n = self.graph.n_nodes
        dd = _direction_vector(gap_direction, n, 'gap_direction')
        dalpha = _direction_vector(alpha_direction, len(self.graph.edges),
                                   'alpha_direction', real=True)
        dz = np.asarray(z_direction, complex)
        if dz.ndim or not np.isfinite(dz):
            raise ValueError('Finite scalar z_direction required')
        m, c = self.graph.area_weights, self.graph.conductance
        a, b = self.a, self.b
        parameter = m[:, None]*np.column_stack((
            np.conj(dd)*a*a-dd+2*dz*a, dd*b*b-np.conj(dd)+2*dz*b))
        tail, head = self.graph.edges.T
        phase = np.exp(-1j*self.alpha)
        for i, j, transport, sign in ((tail, head, phase, 1),
                                      (head, tail, np.conj(phase), -1)):
            dhp = -1j*sign*c*dalpha*transport*self.value.f[j]
            dhm = 1j*sign*c*dalpha*np.conj(transport)*self.value.f_tilde[j]
            np.add.at(parameter[:, 0], i, dhm*a[i]*a[i]-dhp)
            np.add.at(parameter[:, 1], i, dhp*b[i]*b[i]-dhm)
        fixed_da, fixed_db = self._fixed_coordinates(
            fixed_a_direction, fixed_b_direction, fixed_matrix_direction)
        direction = np.zeros((n, 2), complex)
        direction[self.fixed_nodes, 0] = fixed_da
        direction[self.fixed_nodes, 1] = fixed_db
        flat = direction.ravel()
        if self.factor is not None:
            flat[self.components] = self.factor.solve(
                (-parameter.ravel()-self.jacobian@flat)[self.components])
        if np.any(~np.isfinite(direction)):
            raise FloatingPointError('Nonfinite retarded tangent; no repair applied')
        residual = (self.jacobian@flat+parameter.ravel()).reshape(n, 2)
        df, dft, dg, matrix = field_direction(a, b, direction[:, 0], direction[:, 1])
        norm = self.matrix@matrix+matrix@self.matrix
        free_scale = m*np.maximum(1., np.max(abs(parameter/m[:, None]), axis=1))
        fixed_error = np.max(abs(direction[self.fixed_nodes]-
            np.column_stack((fixed_da, fixed_db))), initial=0.)
        return RetardedDirection(direction[:, 0].copy(), direction[:, 1].copy(),
            df, dft, dg, matrix, residual, self.fixed_nodes.copy(),
            residual[self.fixed_nodes].copy(),
            float(np.max(abs(residual[self.free]), initial=0.)),
            float(np.max(abs(residual[self.free])/free_scale[self.free, None], initial=0.)),
            float(fixed_error), float(np.max(abs(norm), initial=0.)))
