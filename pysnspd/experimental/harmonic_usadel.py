"""Finite-frequency complex Nambu response on a stationary Usadel graph.

The caller supplies stationary branches at E+Omega/2 and E-Omega/2, where
Omega=hbar*omega/energy_unit. No temporal-gradient expansion, real projection,
distribution closure or retarded-to-advanced response conjugation is made.
Static Peierls links transport all matrix entries in the same gauge.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu

from .thermal_spatial_usadel import _alpha, _fixed_indices

I2 = np.eye(2, dtype=complex)
TAU3 = np.diag([1., -1.]).astype(complex)


def _matrices(value, size, name):
    raw = np.asarray(value, complex)
    try:
        result = np.broadcast_to(raw, (size, 2, 2)).copy()
    except ValueError as error:
        raise ValueError(name+' must contain one complex 2x2 matrix per node') from error
    if np.any(~np.isfinite(result)):
        raise ValueError(name+' must be finite')
    return result


def _transport(matrices, phase):
    result = np.array(matrices, complex, copy=True)
    result[..., 0, 1] *= phase
    result[..., 1, 0] *= np.conj(phase)
    return result


def stationary_residual(graph, R, B, *, alpha=None):
    """Full stationary matrix equation; fixed-node equations are not discarded here."""
    R, B = (_matrices(value, graph.n_nodes, name) for value, name in ((R, 'R'), (B, 'B')))
    tail, head = graph.edges.T
    phase = np.exp(-1j*_alpha(graph, alpha))
    other = _transport(R[head], phase)
    flux = graph.conductance[:, None, None]/2*(R[tail]@other-other@R[tail])
    result = -graph.area_weights[:, None, None]/2*(B@R-R@B)
    np.add.at(result, tail, flux)
    np.add.at(result, head, -_transport(flux, np.conj(phase)))
    return result


def response_residual(graph, R_plus, R_minus, B_plus, B_minus, delta_R,
                      delta_B, *, alpha=None):
    """Return full four-entry spectral residual, normalization and edge flux."""
    n = graph.n_nodes
    Rp, Rm, Bp, Bm, X, dB = (_matrices(value, n, name) for value, name in
        ((R_plus, 'R_plus'), (R_minus, 'R_minus'), (B_plus, 'B_plus'),
         (B_minus, 'B_minus'), (delta_R, 'delta_R'), (delta_B, 'delta_B')))
    tail, head = graph.edges.T
    phase = np.exp(-1j*_alpha(graph, alpha))
    Rpj, Rmj, Xj = (_transport(field[head], phase) for field in (Rp, Rm, X))
    flux = graph.conductance[:, None, None]/2*(Rp[tail]@Xj+X[tail]@Rmj-
                                              Rpj@X[tail]-Xj@Rm[tail])
    residual = -graph.area_weights[:, None, None]/2*(Bp@X-X@Bm+dB@Rm-Rp@dB)
    np.add.at(residual, tail, flux)
    np.add.at(residual, head, -_transport(flux, np.conj(phase)))
    return residual, Rp@X+X@Rm, flux


@dataclass(frozen=True)
class HarmonicUsadelResponse:
    delta_R: np.ndarray
    spectral_residual: np.ndarray
    normalization_residual: np.ndarray
    edge_flux: np.ndarray
    metrics: dict


class HarmonicUsadelFactor:
    """Reusable sparse LU on the two-dimensional mixed-normalization tangent.

    `R_plus` and `R_minus` can be retarded or advanced backgrounds, but must be
    supplied and solved consistently on the same branch. For advanced response
    use its own factor with B_A=D+(-eta-iE)*tau3. The identity component of
    delta_R is retained: finite-frequency normalization need not imply zero
    trace. A dagger at the same harmonic frequency is not used to infer delta_A.

    A fixed response is an explicit boundary input to solve(); there is no
    inferred contact voltage, reservoir spectrum or physical forcing.
    """
    def __init__(self, graph, R_plus, R_minus, B_plus, B_minus, *, alpha=None,
                 fixed_nodes=(), background_tolerance=1e-7,
                 normalization_tolerance=1e-8, rank_tolerance=1e-8):
        self.graph = graph
        self.R_plus, self.R_minus, self.B_plus, self.B_minus = (
            _matrices(value, graph.n_nodes, name) for value, name in
            ((R_plus, 'R_plus'), (R_minus, 'R_minus'), (B_plus, 'B_plus'), (B_minus, 'B_minus')))
        self.alpha = _alpha(graph, alpha).copy()
        fixed_array = np.asarray(fixed_nodes)
        if fixed_array.size == 0:
            fixed_array = np.array([], dtype=int)
        self.fixed_nodes = _fixed_indices(fixed_array, graph.n_nodes)
        self.free = np.ones(graph.n_nodes, bool); self.free[self.fixed_nodes] = False
        for value in (background_tolerance, normalization_tolerance, rank_tolerance):
            if not np.isfinite(value) or value <= 0:
                raise ValueError('Positive finite background/normalization/rank tolerances required')
        self.normalization_tolerance = float(normalization_tolerance)
        self.background_metrics = {}
        for name, R, B in (('plus', self.R_plus, self.B_plus), ('minus', self.R_minus, self.B_minus)):
            norm = float(np.max(abs(R@R-I2), initial=0.))
            trace = float(np.max(abs(np.trace(R, axis1=1, axis2=2)), initial=0.))
            residual = stationary_residual(graph, R, B, alpha=self.alpha)
            scale = graph.area_weights*np.maximum(1., np.max(abs(B), axis=(1, 2)))
            scaled = float(np.max(abs(residual[self.free])/scale[self.free, None, None], initial=0.))
            self.background_metrics[name] = dict(normalization_max=norm, trace_max=trace,
                                                 full_stationary_scaled_max=scaled)
            if norm > normalization_tolerance or trace > normalization_tolerance:
                raise ValueError('Background '+name+' is not a normalized traceless Nambu involution')
            if scaled > background_tolerance:
                raise ValueError('Background '+name+' is not stationary on the free nodes')
        units = np.eye(4, dtype=complex).reshape(4, 2, 2)
        basis = []
        for Rp, Rm in zip(self.R_plus, self.R_minus):
            constraint = (Rp@units+units@Rm).reshape(4, 4).T
            _, singular, vh = np.linalg.svd(constraint)
            threshold = rank_tolerance*max(1., singular[0])
            if np.count_nonzero(singular > threshold) != 2:
                raise ValueError('Mixed normalization must have exactly two independent constraints')
            basis.append(vh.conj().T[:, 2:])
        self.basis = np.asarray(basis)
        rows, cols, values = [], [], []

        def add_block(i, j, matrix_values):
            block = self.basis[i].conj().T@np.asarray(matrix_values).reshape(2, 4).T
            rr, cc = np.indices((2, 2))
            rows.extend((2*i+rr).ravel()); cols.extend((2*j+cc).ravel()); values.extend(block.ravel())

        for i in range(graph.n_nodes):
            q = self.basis[i].T.reshape(2, 2, 2)
            local = -graph.area_weights[i]/2*(self.B_plus[i]@q-q@self.B_minus[i])
            add_block(i, i, local)
        for (i, j), c, a in zip(graph.edges, graph.conductance, self.alpha):
            phase = np.exp(-1j*a)
            qi = self.basis[i].T.reshape(2, 2, 2)
            qj = self.basis[j].T.reshape(2, 2, 2)
            Rp_j = _transport(self.R_plus[j], phase)
            Rm_j = _transport(self.R_minus[j], phase)
            qj_transported = _transport(qj, phase)
            from_i = c/2*(qi@Rm_j-Rp_j@qi)
            from_j = c/2*(self.R_plus[i]@qj_transported-qj_transported@self.R_minus[i])
            add_block(i, i, from_i); add_block(i, j, from_j)
            add_block(j, i, -_transport(from_i, np.conj(phase)))
            add_block(j, j, -_transport(from_j, np.conj(phase)))
        self.matrix = coo_matrix((values, (rows, cols)),
            shape=(2*graph.n_nodes, 2*graph.n_nodes)).tocsc()
        self.free_components = np.flatnonzero(np.repeat(self.free, 2))
        self.factor = splu(self.matrix[self.free_components][:, self.free_components]) if len(self.free_components) else None

    def solve(self, delta_B, fixed_response=None, *, residual_tolerance=1e-7):
        """Solve one complex forcing; verify all matrix entries after projection.

        The residual is scaled by area*max(1,|B+|,|B-|)*max(1,|delta_R|,|delta_B|)
        per node. Absolute maxima are returned as well. Boundary residuals are
        reported separately because fixed contacts supply their own sources.
        """
        n = self.graph.n_nodes
        dB = _matrices(delta_B, n, 'delta_B')
        if not np.isfinite(residual_tolerance) or residual_tolerance <= 0:
            raise ValueError('Positive finite residual_tolerance required')
        fixed = np.zeros((n, 2, 2), complex)
        if len(self.fixed_nodes):
            if fixed_response is None:
                raise ValueError('Fixed contact responses must be supplied explicitly')
            raw = np.asarray(fixed_response, complex)
            if raw.shape == (n, 2, 2):
                raw = raw[self.fixed_nodes]
            fixed[self.fixed_nodes] = _matrices(raw, len(self.fixed_nodes), 'fixed_response')
            mismatch = self.R_plus@fixed+fixed@self.R_minus
            if np.max(abs(mismatch), initial=0.) > self.normalization_tolerance*max(1., np.max(abs(fixed), initial=0.)):
                raise ValueError('Fixed response violates mixed-frequency normalization')
        elif fixed_response is not None and np.asarray(fixed_response).size:
            raise ValueError('No fixed nodes were declared')
        initial, _, _ = response_residual(self.graph, self.R_plus, self.R_minus,
            self.B_plus, self.B_minus, fixed, dB, alpha=self.alpha)
        projected = np.einsum('nki,nk->ni', self.basis.conj(), initial.reshape(n, 4)).ravel()
        coordinates = np.zeros(2*n, complex)
        if self.factor is not None:
            coordinates[self.free_components] = self.factor.solve(-projected[self.free_components])
        delta_R = fixed+np.einsum('nki,ni->nk', self.basis, coordinates.reshape(n, 2)).reshape(n, 2, 2)
        if np.any(~np.isfinite(delta_R)):
            raise RuntimeError('Nonfinite harmonic response')
        residual, normalization, flux = response_residual(self.graph, self.R_plus, self.R_minus,
            self.B_plus, self.B_minus, delta_R, dB, alpha=self.alpha)
        matrix_scale = np.maximum(1., np.maximum(np.max(abs(self.B_plus), axis=(1, 2)),
                                                  np.max(abs(self.B_minus), axis=(1, 2))))
        response_scale = np.maximum(1., np.maximum(np.max(abs(delta_R), axis=(1, 2)),
                                                    np.max(abs(dB), axis=(1, 2))))
        scale = self.graph.area_weights*matrix_scale*response_scale
        full = float(np.max(abs(residual[self.free])/scale[self.free, None, None], initial=0.))
        normalization_max = float(np.max(abs(normalization), initial=0.))
        linear = self.matrix@coordinates+projected
        metrics = dict(full_spectral_scaled_max=full,
            full_spectral_absolute_max=float(np.max(abs(residual[self.free]), initial=0.)),
            fixed_boundary_residual_max=float(np.max(abs(residual[self.fixed_nodes]), initial=0.)),
            normalization_absolute_max=normalization_max,
            projected_linear_absolute_max=float(np.max(abs(linear[self.free_components]), initial=0.)),
            residual_tolerance=float(residual_tolerance),
            maximum_identity_response=float(np.max(abs(np.trace(delta_R, axis1=1, axis2=2)/2), initial=0.)),
            full_four_entry_residual_checked=True)
        if full > residual_tolerance:
            raise RuntimeError(f'Full harmonic spectral residual failed after projected solve: {full:.6g}')
        if normalization_max > self.normalization_tolerance*max(1., np.max(abs(delta_R), initial=0.)):
            raise RuntimeError('Full mixed-frequency normalization failed')
        return HarmonicUsadelResponse(delta_R, residual, normalization, flux, metrics)
