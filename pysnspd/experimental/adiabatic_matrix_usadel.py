"""First temporal-gradient Nambu algebra on the admitted Usadel graph.

This module supplies residuals, not a claimed closed detector integrator.
All 2x2 entries (including the identity component) are retained.  A temporal
coefficient kappa=hbar/(2*energy_unit) is explicit; kappa=1 means time is
measured in that unit.  Numerical eta is excluded from the kinetic local
term unless a caller explicitly supplies physical self-energies elsewhere.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import lsmr

I2 = np.eye(2, dtype=complex)
TAU3 = np.diag([1., -1.]).astype(complex)


@dataclass(frozen=True)
class Jet:
    """Zeroth value, its E/t derivatives, and independent first correction."""
    value: np.ndarray
    energy: np.ndarray
    time: np.ndarray
    first: np.ndarray

    @classmethod
    def of(cls, value, energy=None, time=None, first=None):
        value = np.asarray(value, complex)
        if value.shape[-2:] != (2, 2) or np.any(~np.isfinite(value)):
            raise ValueError('Finite full 2x2 matrices required')
        fields = [value]
        for field in (energy, time, first):
            result = np.zeros_like(value) if field is None else np.broadcast_to(np.asarray(field, complex), value.shape).copy()
            if np.any(~np.isfinite(result)):
                raise ValueError('Finite jet entries required')
            fields.append(result)
        return cls(*fields)

    def __add__(self, other):
        return Jet(*(getattr(self, f)+getattr(other, f) for f in ('value', 'energy', 'time', 'first')))

    def __sub__(self, other):
        return Jet(*(getattr(self, f)-getattr(other, f) for f in ('value', 'energy', 'time', 'first')))

    def scale(self, scalar):
        scalar = np.asarray(scalar)
        if scalar.ndim:
            scalar = scalar[..., None, None]
        return Jet(*(scalar*getattr(self, f) for f in ('value', 'energy', 'time', 'first')))

    def take(self, index):
        return Jet(*(getattr(self, f)[index] for f in ('value', 'energy', 'time', 'first')))

    def dagger(self):
        return Jet(*(np.swapaxes(np.conj(getattr(self, f)), -1, -2) for f in ('value', 'energy', 'time', 'first')))


def star(a: Jet, b: Jet, kappa=1.):
    """Truncated Moyal product; no product of first corrections is retained."""
    return Jet(a.value@b.value,
               a.energy@b.value+a.value@b.energy,
               a.time@b.value+a.value@b.time,
               a.first@b.value+a.value@b.first
               +1j*kappa*(a.energy@b.time-a.time@b.energy))


def commutator(a, b, kappa=1.):
    return star(a, b, kappa)-star(b, a, kappa)


def advanced(retarded):
    adj = retarded.dagger()
    return Jet(*(-TAU3@getattr(adj, field)@TAU3 for field in ('value', 'energy', 'time', 'first')))


def keldysh(retarded, distribution, kappa=1.):
    return star(retarded, distribution, kappa)-star(distribution, advanced(retarded), kappa)


def distribution(hL, hT, *, energy_L=0., energy_T=0., time_L=0., time_T=0.):
    def matrix(left, right):
        left, right = np.broadcast_arrays(left, right)
        return left[..., None, None]*I2+right[..., None, None]*TAU3
    return Jet.of(matrix(hL, hT), matrix(energy_L, energy_T), matrix(time_L, time_T))


def generator(d, energy, *, d_time=0., voltage=0., voltage_time=0., eta=0.):
    """B=D-iE*tau3+i*v*I, v=e*phi/energy_unit; eta only for R/A."""
    d, dt, v, vt = np.broadcast_arrays(np.asarray(d, complex), d_time, voltage, voltage_time)
    gap = np.zeros(d.shape+(2, 2), complex)
    gap[..., 0, 1], gap[..., 1, 0] = d, d.conj()
    gap_time = np.zeros_like(gap)
    gap_time[..., 0, 1], gap_time[..., 1, 0] = dt, dt.conj()
    return Jet.of(gap+(eta-1j*energy)*TAU3+1j*v[..., None, None]*I2,
                  -1j*TAU3, gap_time+1j*vt[..., None, None]*I2)


def links(alpha, alpha_time=None):
    alpha = np.asarray(alpha, float)
    alpha_time = np.zeros_like(alpha) if alpha_time is None else np.broadcast_to(alpha_time, alpha.shape)
    value = np.zeros(alpha.shape+(2, 2), complex)
    value[..., 0, 0], value[..., 1, 1] = np.exp(-.5j*alpha), np.exp(.5j*alpha)
    return Jet.of(value, time=-.5j*alpha_time[..., None, None]*(TAU3@value))


def transport(link, field, kappa=1.):
    return star(star(link, field, kappa), link.dagger(), kappa)


def _scatter(graph, flux, link, kappa):
    tail, head = graph.edges.T
    backwards = transport(link.dagger(), flux, kappa).scale(-1)
    result = []
    for field in ('value', 'energy', 'time', 'first'):
        array = np.zeros((graph.n_nodes, 2, 2), complex)
        np.add.at(array, tail, getattr(flux, field))
        np.add.at(array, head, getattr(backwards, field))
        result.append(array)
    return Jet(*result)


def spectral_residual(graph, retarded, B, *, alpha=None, alpha_time=None, kappa=1.):
    """Full matrix graph residual and star normalization minus identity."""
    tail, head = graph.edges.T
    link = links(np.zeros(len(tail)) if alpha is None else alpha, alpha_time)
    other = transport(link, retarded.take(head), kappa)
    flux = commutator(retarded.take(tail), other, kappa).scale(graph.conductance/2)
    residual = _scatter(graph, flux, link, kappa)-commutator(B, retarded, kappa).scale(graph.area_weights/2)
    normalization = star(retarded, retarded, kappa)-Jet.of(np.broadcast_to(I2, retarded.value.shape))
    return residual, normalization


def kinetic_residual(graph, retarded, h, B_without_eta, *, alpha=None, alpha_time=None, kappa=1.):
    """Matrix kinetic residual with no inferred collision or relaxation bath.

    Supply the common physical B, excluding numerical eta.  The caller must
    add a conserving collision self-energy, neutrality and gap self-consistency
    before treating this residual as a complete time-evolution system.
    """
    if np.any(np.abs(np.trace(B_without_eta.energy, axis1=-2, axis2=-1)) > 1e-12):
        raise ValueError('Unexpected scalar energy dependence in physical generator')
    tail, head = graph.edges.T
    link = links(np.zeros(len(tail)) if alpha is None else alpha, alpha_time)
    R, A = retarded, advanced(retarded)
    K = keldysh(R, h, kappa)
    Rj, Aj, Kj = (transport(link, x.take(head), kappa) for x in (R, A, K))
    flux = (star(R.take(tail), Kj, kappa)+star(K.take(tail), Aj, kappa)
            -star(Rj, K.take(tail), kappa)-star(Kj, A.take(tail), kappa)).scale(graph.conductance/2)
    local = commutator(B_without_eta, K, kappa).scale(graph.area_weights/2)
    return _scatter(graph, flux, link, kappa)-local, flux, K


def projections(matrix):
    """Real L/T trace projections, identical to FrozenKineticOperator."""
    return np.stack(((matrix[..., 0, 0]+matrix[..., 1, 1]).real/4,
                     (matrix[..., 0, 0]-matrix[..., 1, 1]).real/4), axis=-1)


def normalization_particular(retarded, kappa=1.):
    """Required first correction commuting with R0; usually includes I2."""
    source = -1j*kappa*(retarded.energy@retarded.time-retarded.time@retarded.energy)
    return .5*retarded.value@source


def solve_first_spectral(graph, retarded, B, *, alpha=None, alpha_time=None,
                         fixed_nodes=(), fixed_first=None, kappa=1., tolerance=1e-10,
                         maxiter=None):
    """Constrained sparse first-order R1 solve, retaining all four entries.

    Stacks linearized spectral residual and normalization; LSMR resolves the
    consistent overdetermined equations without a fictitious damping term.
    Exact residuals are returned: no least-squares status is an acceptance.
    Energy/time derivatives must already be derivatives of the same R0.
    """
    n = graph.n_nodes
    zero = replace(retarded, first=np.zeros_like(retarded.value))
    base, norm = spectral_residual(graph, zero, B, alpha=alpha, alpha_time=alpha_time, kappa=kappa)
    R, b = zero.value, B.value
    basis = np.eye(4, dtype=complex).reshape(4, 2, 2)
    row, col, val = [], [], []
    def block(i, j, values, offset=0):
        values = np.asarray(values).reshape(4, 4).T
        rr, cc = np.indices((4, 4))
        row.extend((offset+4*i+rr).ravel()); col.extend((4*j+cc).ravel()); val.extend(values.ravel())
    for i in range(n):
        block(i, i, -.5*(b[i]@basis-basis@b[i]))
        block(i, i, R[i]@basis+basis@R[i], offset=4*n)
    alpha = np.zeros(len(graph.edges)) if alpha is None else np.asarray(alpha)
    p = links(alpha).value
    for edge, (i, j) in enumerate(graph.edges):
        P = p[edge]; Pd = P.conj().T
        Rj, Ri = P@R[j]@Pd, Pd@R[i]@P
        pb, pdb = P@basis@Pd, Pd@basis@P
        c = graph.conductance[edge]/2
        block(i, i, c/graph.area_weights[i]*(basis@Rj-Rj@basis))
        block(i, j, c/graph.area_weights[i]*(R[i]@pb-pb@R[i]))
        block(j, j, c/graph.area_weights[j]*(basis@Ri-Ri@basis))
        block(j, i, c/graph.area_weights[j]*(R[j]@pdb-pdb@R[j]))
    matrix = coo_matrix((np.asarray(val), (np.asarray(row), np.asarray(col))), shape=(8*n, 4*n)).tocsr()
    rhs = -np.concatenate(((base.first/graph.area_weights[:, None, None]).ravel(), norm.first.ravel()))
    fixed = np.asarray(fixed_nodes, int)
    if np.any(fixed < 0) or np.any(fixed >= n) or len(np.unique(fixed)) != len(fixed):
        raise ValueError('Unique valid fixed nodes required')
    known = np.zeros((n, 2, 2), complex)
    if len(fixed):
        known[fixed] = 0 if fixed_first is None else np.broadcast_to(fixed_first, (len(fixed), 2, 2))
    columns = np.ones((n, 4), bool); columns[fixed] = False; columns = columns.ravel()
    rows = np.ones((2, n, 4), bool); rows[0, fixed] = False; rows[1, fixed] = False; rows = rows.ravel()
    rhs = (rhs-matrix@known.ravel())[rows]
    system = matrix[rows][:, columns]
    answer = lsmr(system, rhs, atol=tolerance, btol=tolerance, maxiter=maxiter)
    solved = known.ravel(); solved[columns] = answer[0]; solved = solved.reshape(n, 2, 2)
    result = replace(retarded, first=solved)
    residual, normalization = spectral_residual(graph, result, B, alpha=alpha, alpha_time=alpha_time, kappa=kappa)
    free = np.ones(n, bool); free[fixed] = False
    metrics = {'lsmr_stop': int(answer[1]), 'iterations': int(answer[2]),
               'spectral_first_max': float(np.max(abs(residual.first[free]/graph.area_weights[free, None, None]), initial=0)),
               'normalization_first_max': float(np.max(abs(normalization.first[free]), initial=0)),
               'spectral_zero_max': float(np.max(abs(residual.value[free]/graph.area_weights[free, None, None]), initial=0)),
               'normalization_zero_max': float(np.max(abs(normalization.value), initial=0)),
               'identity_first_max': float(np.max(abs(np.trace(solved, axis1=-2, axis2=-1)/2), initial=0))}
    return result, metrics
