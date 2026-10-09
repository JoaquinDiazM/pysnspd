"""Two-mode temporal RHS of the first-gradient spatial Usadel residual.

This solves distributions on a SUPPLIED moving spectrum; it does not close
the gap, electrostatic neutrality, collision self-energy or detector circuit.
Time derivatives are per second, energy derivatives per dimensionless E/Eu,
and kappa_s=hbar/(2*Eu). Both modes and nonlocal temporal storage are retained.
No inverse-DOS floor, eta collision bath, least-squares repair or clipping.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import warnings
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import MatrixRankWarning, spsolve

from . import adiabatic_matrix_usadel as ad


def _real_fields(value, size, name):
    raw = np.asarray(value)
    if np.iscomplexobj(raw) or raw.shape != (size, 2) or np.any(~np.isfinite(raw)):
        raise ValueError(name+' must be finite real (hL,hT) rows per node')
    return np.asarray(raw, float)


def _positive(value, name):
    raw = np.asarray(value)
    if (raw.ndim or np.iscomplexobj(raw) or raw.dtype.kind == 'b'
            or not np.isfinite(raw) or raw <= 0):
        raise ValueError(name+' must be a positive real scalar')
    return float(raw)


def _stationary_jet(jet):
    """Coefficients of hdot use R0, RE only: R1*hdot is second order."""
    return replace(jet, time=np.zeros_like(jet.value), first=np.zeros_like(jet.value))


def _flux(Ri, Rj, Ki, Kj, Ai, Aj, coefficient, kappa):
    return (ad.star(Ri, Kj, kappa)+ad.star(Ki, Aj, kappa)
            -ad.star(Rj, Ki, kappa)-ad.star(Kj, Ai, kappa)).scale(coefficient/2)


def assemble_storage(graph, retarded, B_without_eta, *, alpha=None, kappa_s):
    """Sparse exact coefficient M of hdot in Pi(R0+R1), at each energy.

    The 2x2 node/edge blocks are differentiated analytically through the same
    Jet products as the residual. Work is linear in nodes+edges; no full-graph
    impulse per unknown and no assumption M=-kappa*mass*rho*I are made.
    M has units seconds because hdot is supplied in s^-1.
    """
    kappa = _positive(kappa_s, 'kappa_s')
    n = graph.n_nodes
    if retarded.value.shape != (n, 2, 2) or B_without_eta.value.shape != (n, 2, 2):
        raise ValueError('One Nambu matrix per graph node required')
    R, B = _stationary_jet(retarded), _stationary_jet(B_without_eta)
    A = ad.advanced(R)
    tail, head = graph.edges.T
    alpha = np.zeros(len(tail)) if alpha is None else np.asarray(alpha)
    if np.iscomplexobj(alpha) or alpha.shape != tail.shape or np.any(~np.isfinite(alpha)):
        raise ValueError('One finite real link phase per graph edge required')
    P = ad.links(alpha)
    Ri, Ai = R.take(tail), A.take(tail)
    Rj, Aj = ad.transport(P, R.take(head), kappa), ad.transport(P, A.take(head), kappa)
    node_blocks = np.zeros((n, 2, 2))
    tail_blocks = np.zeros((len(tail), 2, 2))
    head_blocks = np.zeros_like(tail_blocks)
    for component, basis in enumerate((ad.I2, ad.TAU3)):
        H = ad.Jet.of(np.zeros_like(R.value), time=np.broadcast_to(basis, R.value.shape))
        K = ad.keldysh(R, H, kappa)
        local = ad.commutator(B, K, kappa).scale(-graph.area_weights/2)
        node_blocks[:, :, component] = ad.projections(local.first)
        Ki = K.take(tail)
        Kj = ad.transport(P, K.take(head), kappa)
        zero = ad.Jet.of(np.zeros_like(Ki.value))
        tail_flux = _flux(Ri, Rj, Ki, zero, Ai, Aj, graph.conductance, kappa)
        head_flux = _flux(Ri, Rj, zero, Kj, Ai, Aj, graph.conductance, kappa)
        tail_blocks[:, :, component] = ad.projections(tail_flux.first)
        head_blocks[:, :, component] = ad.projections(head_flux.first)
    rows, columns, values = [], [], []
    nodes = np.arange(n)
    for output in range(2):
        for component in range(2):
            rows.extend((2*nodes+output, 2*tail+output, 2*tail+output,
                         2*head+output, 2*head+output))
            columns.extend((2*nodes+component, 2*tail+component, 2*head+component,
                            2*tail+component, 2*head+component))
            values.extend((node_blocks[:, output, component], tail_blocks[:, output, component],
                           head_blocks[:, output, component], -tail_blocks[:, output, component],
                           -head_blocks[:, output, component]))
    matrix = coo_matrix((np.concatenate(values),
        (np.concatenate(rows), np.concatenate(columns))), shape=(2*n, 2*n)).tocsr()
    matrix.eliminate_zeros()
    return matrix


@dataclass(frozen=True)
class MatrixKineticRHS:
    distribution_rhs_per_s: np.ndarray
    storage_matrix_s: object
    frozen_residual: np.ndarray
    known_temporal_residual: np.ndarray
    projected_source: np.ndarray
    projected_residual: np.ndarray
    boundary_reaction: np.ndarray
    edge_flux_zero: np.ndarray
    edge_flux_first: np.ndarray
    full_nambu_residual: np.ndarray
    free_projected_max: float
    linear_solve_max: float
    pauli_margin: float
    kappa_s: float
    source_contract: str
    scope: str = 'two-mode kinetic RHS on supplied first-gradient spectrum; no detector admission'


def solve_distribution_rhs(graph, retarded, h, h_energy, B_without_eta, *,
                           kappa_s, projected_source, source_contract,
                           alpha=None, alpha_time_per_s=None,
                           fixed_nodes=(), fixed_rhs_per_s=None,
                           require_pauli_support=True):
    """Solve M*hdot=C-Pi(R0)-Pi(R1 at hdot=0), retaining both modes.

    projected_source is the SUPPLIED two trace projections of the physical
    collision/source residual, integrated per graph node. Zero is explicit.
    This module does not convert fitted rates or power to self-energies.
    source_contract must state the approximation and signs. In the normal
    limit C=-kappa*m*hdot_collision. Eta in a retarded regulator never becomes
    part of B or C. Prescribed-node reactions are reported, not discarded.
    The full matrix residual is diagnostic only: a two-trace collision input
    does not determine all anomalous collision entries or gap self-consistency.
    """
    kappa = _positive(kappa_s, 'kappa_s')
    n = graph.n_nodes
    h = _real_fields(h, n, 'h')
    hE = _real_fields(h_energy, n, 'h_energy')
    source = _real_fields(projected_source, n, 'projected_source')
    if not isinstance(source_contract, str) or not source_contract.strip():
        raise ValueError('Explicit physical source/self-energy contract required')
    support = 1-np.max(abs(np.column_stack((h[:, 0]+h[:, 1], h[:, 0]-h[:, 1]))))
    if require_pauli_support and support < -32*np.finfo(float).eps:
        raise ValueError('hL +/- hT leave Pauli support; no clipping applied')
    # B must explicitly exclude the usual numerical eta*tau3 regulator.
    if np.any(abs((B_without_eta.value[:, 0, 0]-B_without_eta.value[:, 1, 1]).real) > 1e-12):
        raise ValueError('Numerical eta is forbidden in physical kinetic B')
    M = assemble_storage(graph, retarded, B_without_eta, alpha=alpha, kappa_s=kappa)
    H0 = ad.distribution(h[:, 0], h[:, 1], energy_L=hE[:, 0], energy_T=hE[:, 1])
    base, _, _ = ad.kinetic_residual(graph, retarded, H0, B_without_eta,
        alpha=alpha, alpha_time=alpha_time_per_s, kappa=kappa)
    frozen, temporal = ad.projections(base.value), ad.projections(base.first)
    rhs = (source-frozen-temporal).ravel()
    fixed = np.asarray(fixed_nodes)
    if fixed.size == 0:
        fixed = np.empty(0, int)
    if (fixed.dtype.kind not in 'iu' or fixed.ndim != 1 or np.any(fixed < 0)
            or np.any(fixed >= n) or len(np.unique(fixed)) != len(fixed)):
        raise ValueError('Unique integer prescribed nodes required')
    fixed = fixed.astype(int, copy=False)
    value = np.zeros((n, 2))
    if len(fixed):
        if fixed_rhs_per_s is None:
            raise ValueError('Explicit time law of prescribed distributions required')
        value[fixed] = _real_fields(fixed_rhs_per_s, len(fixed), 'fixed_rhs_per_s')
    elif fixed_rhs_per_s is not None and np.asarray(fixed_rhs_per_s).size:
        raise ValueError('Boundary RHS supplied without prescribed nodes')
    free = np.ones(n, bool)
    free[fixed] = False
    components = np.flatnonzero(np.repeat(free, 2))
    vector = value.ravel()
    if len(components):
        with warnings.catch_warnings():
            warnings.simplefilter('error', MatrixRankWarning)
            try:
                vector[components] = spsolve(M[components][:, components],
                    (rhs-M@vector)[components])
            except MatrixRankWarning as error:
                raise FloatingPointError('Singular physical temporal storage; no DOS floor or regularization') from error
    if np.any(~np.isfinite(vector)):
        raise FloatingPointError('Nonfinite kinetic RHS; no repair applied')
    H = ad.distribution(h[:, 0], h[:, 1], energy_L=hE[:, 0], energy_T=hE[:, 1],
                        time_L=value[:, 0], time_T=value[:, 1])
    residual, flux, _ = ad.kinetic_residual(graph, retarded, H, B_without_eta,
        alpha=alpha, alpha_time=alpha_time_per_s, kappa=kappa)
    projected = ad.projections(residual.value+residual.first)-source
    solve_residual = M@vector-rhs
    return MatrixKineticRHS(value.copy(), M, frozen, temporal, source.copy(), projected,
        projected[fixed].copy(), ad.projections(flux.value), ad.projections(flux.first),
        residual.value+residual.first, float(np.max(abs(projected[free]), initial=0)),
        float(np.max(abs(solve_residual[components]), initial=0)), float(support), kappa,
        source_contract)
