"""Two-mode kinetic response with exact E+/- shifts on the Usadel graph.

Convention exp(-i*omega*t), E+/-=E+/-hbar*omega/(2*energy_unit).
All response amplitudes and projections are complex.  This is an explicit
linear-response building block, not a self-consistent detector closure.
Static link phases are allowed; a perturbation of vector potential requires
its own link vertex and is not silently inferred here.
"""
from __future__ import annotations

from dataclasses import dataclass
import warnings
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import MatrixRankWarning, spsolve

I2 = np.eye(2, dtype=complex)
TAU3 = np.diag([1., -1.]).astype(complex)


def project(matrix):
    """Complex harmonic L/T projections. Taking real here would lose phase."""
    matrix = np.asarray(matrix, complex)
    return np.stack(((matrix[..., 0, 0]+matrix[..., 1, 1])/4,
                     (matrix[..., 0, 0]-matrix[..., 1, 1])/4), axis=-1)


def gap_cartesian_moments(K, area_weights):
    """Two Fourier amplitudes Gx,Gy; conjugation would reverse frequency."""
    K = np.asarray(K, complex)
    mass = np.asarray(area_weights, float)
    return np.column_stack((.5j*mass*(K[:, 0, 1]+K[:, 1, 0]),
                             .5*mass*(K[:, 0, 1]-K[:, 1, 0])))


def gap_vertex(delta_x, delta_y, delta_voltage=0.):
    """delta B=delta D+i*delta_v*I, with independent complex x,y amplitudes."""
    x, y, v = np.broadcast_arrays(delta_x, delta_y, delta_voltage)
    result = np.zeros(x.shape+(2, 2), complex)
    result[..., 0, 1], result[..., 1, 0] = x+1j*y, x-1j*y
    result[..., 0, 0] = result[..., 1, 1] = 1j*v
    return result


def _matrices(value, n, name):
    array = np.asarray(value, complex)
    try:
        array = np.broadcast_to(array, (n, 2, 2)).copy()
    except ValueError as error:
        raise ValueError(f'{name}: expected one or {n} full 2x2 matrices') from error
    if np.any(~np.isfinite(array)):
        raise ValueError(f'{name}: finite matrices required')
    return array


def _scalars(value, n, name):
    array = np.broadcast_to(np.asarray(value, complex), (n,)).copy()
    if np.any(~np.isfinite(array)):
        raise ValueError(f'{name}: finite values required')
    return array


@dataclass(frozen=True)
class HarmonicKineticObservation:
    residual: np.ndarray
    edge_flux: np.ndarray
    delta_K: np.ndarray
    local_residual: np.ndarray
    gap_cartesian: np.ndarray
    charge_current: np.ndarray


@dataclass(frozen=True)
class HarmonicKineticOperator:
    graph: object
    K_basis: np.ndarray
    K_source: np.ndarray
    pair_blocks: np.ndarray
    edge_tail_blocks: np.ndarray
    edge_head_blocks: np.ndarray
    local_source: np.ndarray
    edge_source: np.ndarray
    source: np.ndarray
    matrix: object

    def evaluate(self, delta_h):
        h = np.asarray(delta_h, complex)
        if h.shape != (self.graph.n_nodes, 2) or np.any(~np.isfinite(h)):
            raise ValueError('One finite complex (delta_hL,delta_hT) row per node required')
        tail, head = self.graph.edges.T
        local = self.local_source+np.einsum('noc,nc->no', self.pair_blocks, h)
        flux = self.edge_source+np.einsum('eoc,ec->eo', self.edge_tail_blocks, h[tail])
        flux += np.einsum('eoc,ec->eo', self.edge_head_blocks, h[head])
        residual = local.copy()
        np.add.at(residual, tail, flux); np.add.at(residual, head, -flux)
        K = self.K_source+np.einsum('nc,ncab->nab', h, self.K_basis)
        return HarmonicKineticObservation(residual, flux, K, local,
            gap_cartesian_moments(K, self.graph.area_weights), 2*flux[:, 1])

    def solve(self, fixed_nodes, fixed_h=0., *, forcing=None, residual_tolerance=1e-8):
        """Solve matrix*h+source=forcing with explicitly fixed distributions.

        A nonzero forcing is a declared external kinetic source, not inferred
        heat. Its moments must be included by the caller in its energy ledger.
        """
        n = self.graph.n_nodes
        fixed = np.asarray(fixed_nodes, int)
        if fixed.ndim != 1 or not len(fixed) or np.any(fixed < 0) or np.any(fixed >= n) or len(np.unique(fixed)) != len(fixed):
            raise ValueError('Explicit unique valid distribution contacts required')
        h = np.zeros((n, 2), complex)
        known = np.asarray(fixed_h, complex)
        if known.shape == (n, 2):
            known = known[fixed]
        h[fixed] = np.broadcast_to(known, (len(fixed), 2))
        if np.any(~np.isfinite(h)):
            raise ValueError('Finite fixed distributions required')
        rhs_external = np.zeros((n, 2), complex) if forcing is None else np.asarray(forcing, complex)
        if rhs_external.shape != (n, 2) or np.any(~np.isfinite(rhs_external)):
            raise ValueError('Finite forcing in each projected kinetic equation required')
        free = np.ones((n, 2), bool); free[fixed] = False; free = free.ravel()
        rhs = (rhs_external-self.source).ravel()-self.matrix@h.ravel()
        with warnings.catch_warnings():
            warnings.simplefilter('error', MatrixRankWarning)
            h.ravel()[free] = spsolve(self.matrix[free][:, free], rhs[free])
        if np.any(~np.isfinite(h)):
            raise RuntimeError('Nonfinite harmonic kinetic response; no regularization added')
        result = self.evaluate(h)
        error = (result.residual-rhs_external).ravel()[free]
        scale = max(float(np.max(abs(rhs[free]), initial=0)), 1.)
        maximum = float(np.max(abs(error), initial=0))
        if maximum > residual_tolerance*scale:
            raise RuntimeError(f'Harmonic kinetic solve residual {maximum} exceeds {residual_tolerance*scale}')
        return h, result, {'maximum_free_residual': maximum, 'residual_scale': scale,
            'free_variables': int(np.count_nonzero(free)),
            'fixed_nodes': int(len(fixed)), 'time_convention': 'exp(-i omega t)'}


def assemble(graph, R_plus, R_minus, A_plus, A_minus, B_plus, B_minus,
             delta_R, delta_A, delta_B, h_plus, h_minus, *, alpha=None):
    """Assemble exact finite-frequency kinetic response for given spectra.

    B_plus/minus must be the *physical* generators D-iE_plus/minus*tau3,
    without numerical eta. The spectral backgrounds and responses may use
    causal eta, whose effect must then be assessed as regulator dependence.
    All source terms, including moving gap and electric potential, are kept.
    """
    n = graph.n_nodes
    Rp, Rm, Ap, Am, Bp, Bm, dR, dA, dB = (_matrices(value, n, name) for value, name in (
        (R_plus, 'R_plus'), (R_minus, 'R_minus'), (A_plus, 'A_plus'), (A_minus, 'A_minus'),
        (B_plus, 'B_plus'), (B_minus, 'B_minus'), (delta_R, 'delta_R'), (delta_A, 'delta_A'), (delta_B, 'delta_B')))
    hp, hm = _scalars(h_plus, n, 'h_plus'), _scalars(h_minus, n, 'h_minus')
    # Physical B has purely imaginary diagonals for real E and background phi.
    # Reject an accidentally imported real eta*tau3 from the spectral problem.
    if np.max(abs(np.real(Bp[:, (0, 1), (0, 1)]))) > 1e-12 or np.max(abs(np.real(Bm[:, (0, 1), (0, 1)]))) > 1e-12:
        raise ValueError('Physical kinetic B must exclude the real numerical eta*tau3 regulator')
    H_basis = np.stack((I2, TAU3))
    K_basis = Rp[:, None]@H_basis-H_basis@Am[:, None]
    K_source = dR*hm[:, None, None]-hp[:, None, None]*dA
    Kp, Km = hp[:, None, None]*(Rp-Ap), hm[:, None, None]*(Rm-Am)
    tail, head = graph.edges.T
    alpha = np.zeros(len(tail)) if alpha is None else np.asarray(alpha, float)
    if alpha.shape != (len(tail),) or np.any(~np.isfinite(alpha)):
        raise ValueError('One finite static link phase per edge required')
    transport = np.ones((len(tail), 2, 2), complex)
    transport[:, 0, 1], transport[:, 1, 0] = np.exp(-1j*alpha), np.exp(1j*alpha)
    def at_head(value):
        return value[head]*transport
    Rpj, Amj, Kpj, Kmj, dRj, dAj, Ksj = (at_head(value) for value in (Rp, Am, Kp, Km, dR, dA, K_source))
    Kbj = K_basis[head]*transport[:, None]
    coefficient = graph.conductance[:, None, None, None]/2
    tail_matrices = coefficient*(K_basis[tail]@Amj[:, None]-Rpj[:, None]@K_basis[tail])
    head_matrices = coefficient*(Rp[tail, None]@Kbj-Kbj@Am[tail, None])
    tail_blocks, head_blocks = project(tail_matrices).transpose(0, 2, 1), project(head_matrices).transpose(0, 2, 1)
    mass = graph.area_weights[:, None, None, None]/2
    local_blocks = project(-mass*(Bp[:, None]@K_basis-K_basis@Bm[:, None])).transpose(0, 2, 1)
    source_flux_matrix = graph.conductance[:, None, None]/2*(
        Rp[tail]@Ksj+dR[tail]@Kmj+Kp[tail]@dAj+K_source[tail]@Amj
        -Rpj@K_source[tail]-dRj@Km[tail]-Kpj@dA[tail]-Ksj@Am[tail])
    edge_source = project(source_flux_matrix)
    local_source = project(-graph.area_weights[:, None, None]/2*(
        Bp@K_source+dB@Km-K_source@Bm-Kp@dB))
    source = local_source.copy()
    np.add.at(source, tail, edge_source); np.add.at(source, head, -edge_source)
    rows, cols, values = [], [], []
    nodes = np.arange(n)
    for output in range(2):
        for component in range(2):
            rows += [2*nodes+output, 2*tail+output, 2*tail+output, 2*head+output, 2*head+output]
            cols += [2*nodes+component, 2*tail+component, 2*head+component, 2*tail+component, 2*head+component]
            values += [local_blocks[:, output, component], tail_blocks[:, output, component],
                head_blocks[:, output, component], -tail_blocks[:, output, component], -head_blocks[:, output, component]]
    matrix = coo_matrix((np.concatenate(values), (np.concatenate(rows), np.concatenate(cols))), shape=(2*n, 2*n)).tocsr()
    return HarmonicKineticOperator(graph, K_basis, K_source, local_blocks,
        tail_blocks, head_blocks, local_source, edge_source, source, matrix)


def modal_coefficients(eigenvalue, R_plus, R_minus, A_plus, A_minus,
                       B_plus, B_minus, delta_R, delta_A, delta_B,
                       h_plus, h_minus):
    """Exact scalar graph-mode restriction for spatially uniform backgrounds.

    eigenvalue is that of positive L v=lambda M v. All delta fields share v;
    this helper returns residual per nodal mass and per mode amplitude. It is
    not valid for a spatially varying or biased spectral background by merely
    replacing its coefficients with averages.
    """
    lam = float(eigenvalue)
    if not np.isfinite(lam) or lam < 0:
        raise ValueError('Finite nonnegative positive-Laplacian eigenvalue required')
    Rp, Rm, Ap, Am, Bp, Bm, dR, dA, dB = (_matrices(value, 1, name)[0] for value, name in (
        (R_plus, 'R_plus'), (R_minus, 'R_minus'), (A_plus, 'A_plus'), (A_minus, 'A_minus'),
        (B_plus, 'B_plus'), (B_minus, 'B_minus'), (delta_R, 'delta_R'), (delta_A, 'delta_A'), (delta_B, 'delta_B')))
    hp, hm = complex(h_plus), complex(h_minus)
    Kp, Km = hp*(Rp-Ap), hm*(Rm-Am)
    K_source = dR*hm-hp*dA
    K_basis = np.stack([Rp@H-H@Am for H in (I2, TAU3)])
    def residual(K, dr, da, db):
        divergence = lam/2*(dr@Km+K@Am-Rp@K-Kp@da)
        local = -.5*(Bp@K+db@Km-K@Bm-Kp@db)
        return project(divergence+local), project(divergence)
    zero = np.zeros((2, 2), complex)
    columns = [residual(K, zero, zero, zero) for K in K_basis]
    source, source_divergence = residual(K_source, dR, dA, dB)
    return {'matrix': np.column_stack([column[0] for column in columns]),
            'source': source, 'K_basis': K_basis, 'K_source': K_source,
            'divergence_matrix': np.column_stack([column[1] for column in columns]),
            'divergence_source': source_divergence}
