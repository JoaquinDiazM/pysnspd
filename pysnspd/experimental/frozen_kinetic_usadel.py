"""Static Keldysh response on a frozen spatial retarded Usadel spectrum.

This is a two-component linear diagnostic, not a population time integrator.
hL and hT are real distribution directions. The pair-conversion commutator is
retained. Numerical eta is excluded from collisions and its fictitious bath
leakage is returned separately, never interpreted as a physical relaxation.
"""
from __future__ import annotations

from dataclasses import dataclass
import warnings
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import MatrixRankWarning, spsolve

from .thermal_spatial_usadel import _alpha, _complex_vector, _fixed_indices

TAU3 = np.diag([1., -1.]).astype(complex)


def _project(matrix):
    """Real charge/energy projections: Tr(I/τ3 times matrix)/4."""
    return np.stack(((matrix[..., 0, 0]+matrix[..., 1, 1]).real/4,
                     (matrix[..., 0, 0]-matrix[..., 1, 1]).real/4), axis=-1)


def _real_fields(h, size):
    raw = np.asarray(h)
    if np.iscomplexobj(raw) or raw.shape != (size, 2) or np.any(~np.isfinite(raw)):
        raise ValueError('One finite real (hL,hT) row per node required')
    return np.asarray(raw, float)


@dataclass(frozen=True)
class KineticObservation:
    """pair_conversion is the residual contribution, minus the RHS source."""
    residual: np.ndarray
    edge_flux: np.ndarray
    pair_conversion: np.ndarray
    gap_force_increment: np.ndarray
    charge_current_increment: np.ndarray
    ward_residual: np.ndarray
    artificial_eta_leakage: np.ndarray


@dataclass(frozen=True)
class FrozenKineticOperator:
    graph: object
    d: np.ndarray
    R: np.ndarray
    A: np.ndarray
    K_basis: np.ndarray
    pair_blocks: np.ndarray
    edge_tail_blocks: np.ndarray
    edge_head_blocks: np.ndarray
    matrix: object

    def evaluate(self, directions, *, eta=0.):
        h = _real_fields(directions, self.graph.n_nodes)
        if not np.isfinite(eta) or eta < 0:
            raise ValueError('Nonnegative numerical eta required')
        tail, head = self.graph.edges.T
        flux = np.einsum('eoc,ec->eo', self.edge_tail_blocks, h[tail])
        flux += np.einsum('eoc,ec->eo', self.edge_head_blocks, h[head])
        pair = np.einsum('noc,nc->no', self.pair_blocks, h)
        residual = pair.copy()
        np.add.at(residual, tail, flux); np.add.at(residual, head, -flux)
        K = np.einsum('nc,ncab->nab', h, self.K_basis)
        force = 1j*self.graph.area_weights*K[:, 0, 1]
        current = 2*flux[:, 1]
        divergence = np.zeros(self.graph.n_nodes)
        np.add.at(divergence, tail, current); np.add.at(divergence, head, -current)
        ward = divergence+np.imag(np.conj(self.d)*force)-2*residual[:, 1]
        leakage = self.graph.area_weights[:, None]*eta*self.R[:, 0, 0].real[:, None]*h
        return KineticObservation(residual, flux, pair, force, current, ward, leakage)

    def charge_response(self, delta_hL, fixed_nodes):
        hL = np.asarray(delta_hL)
        if np.iscomplexobj(hL) or hL.shape != (self.graph.n_nodes,) or np.any(~np.isfinite(hL)):
            raise ValueError('One finite real prescribed energy-mode direction per node required')
        fixed = _fixed_indices(fixed_nodes, self.graph.n_nodes)
        if not len(fixed):
            raise ValueError('Explicit fixed charge-mode contacts required; no gauge or electrostatic closure is inferred')
        free = np.ones(self.graph.n_nodes, bool); free[fixed] = False
        if np.any(hL[fixed] != 0):
            raise ValueError('This diagnostic requires zero distribution increments at prescribed contacts')
        charge = np.arange(1, 2*self.graph.n_nodes, 2)
        energy = charge-1
        Ltt = self.matrix[charge][:, charge]
        rhs = -np.asarray(self.matrix[charge][:, energy]@hL)
        hT = np.zeros(self.graph.n_nodes)
        with warnings.catch_warnings():
            warnings.simplefilter('error', MatrixRankWarning)
            hT[free] = spsolve(Ltt[free][:, free], rhs[free])
        if np.any(~np.isfinite(hT)):
            raise RuntimeError('Nonfinite frozen-spectrum charge response; no regularizing leakage added')
        return np.column_stack((hL, hT))


def assemble(graph, d, g, f, f_tilde, alpha=None):
    size = graph.n_nodes
    d, g, f, ft = (_complex_vector(value, size, name) for value, name in
                   ((d, 'd'), (g, 'g'), (f, 'f'), (f_tilde, 'f_tilde')))
    if np.max(abs(g*g+f*ft-1)) > 1e-7 or np.min(g.real) < -1e-8:
        raise ValueError('Normalized causal frozen spectral fields required')
    R = np.empty((size, 2, 2), complex)
    R[:, 0, 0], R[:, 0, 1], R[:, 1, 0], R[:, 1, 1] = g, f, ft, -g
    A = -TAU3@np.conj(R).transpose(0, 2, 1)@TAU3
    K = np.stack((R-A, R@TAU3-TAU3@A), axis=1)
    gap = np.zeros_like(R); gap[:, 0, 1], gap[:, 1, 0] = d, np.conj(d)
    pair_matrix = -.5*graph.area_weights[:, None, None, None]*(gap[:, None]@K-K@gap[:, None])
    pair_blocks = _project(pair_matrix).transpose(0, 2, 1)
    tail, head = graph.edges.T
    phase = np.exp(-1j*_alpha(graph, alpha))
    transport = np.ones((len(tail), 2, 2), complex)
    transport[:, 0, 1], transport[:, 1, 0] = phase, np.conj(phase)
    Rj, Aj, Kj = R[head]*transport, A[head]*transport, K[head]*transport[:, None]
    c = graph.conductance[:, None, None, None]/2
    tail_matrices = c*(K[tail]@Aj[:, None]-Rj[:, None]@K[tail])
    head_matrices = c*(R[tail, None]@Kj-Kj@A[tail, None])
    tail_blocks = _project(tail_matrices).transpose(0, 2, 1)
    head_blocks = _project(head_matrices).transpose(0, 2, 1)
    rows, cols, values = [], [], []
    nodes = np.arange(size)
    for output in range(2):
        for component in range(2):
            rows += [2*nodes+output, 2*tail+output, 2*tail+output, 2*head+output, 2*head+output]
            cols += [2*nodes+component, 2*tail+component, 2*head+component, 2*tail+component, 2*head+component]
            values += [pair_blocks[:, output, component], tail_blocks[:, output, component],
                head_blocks[:, output, component], -tail_blocks[:, output, component], -head_blocks[:, output, component]]
    matrix = coo_matrix((np.concatenate(values), (np.concatenate(rows), np.concatenate(cols))),
                        shape=(2*size, 2*size)).tocsr()
    return FrozenKineticOperator(graph, d.copy(), R, A, K, pair_blocks, tail_blocks, head_blocks, matrix)
