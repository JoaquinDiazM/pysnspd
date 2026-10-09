"""Common SI moment rows for a simultaneous first-gradient detector RHS.

This module integrates supplied full Nambu fields; it never solves a spectrum,
adds Ohmic current, assumes a DOS floor, changes populations, or supplies a gap
friction.  Gap self-consistency is G0+G1=0 with Gamma_res=0.  The rank of that
candidate is a property of the assembled physical system, not guaranteed here.
Interior current divergence is returned for a Ward/continuity diagnostic.  It
is not silently appended to the kinetic-T and phase-gap equations, which can
contain the same conservation identity.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.constants import elementary_charge
from scipy.sparse import coo_matrix, csr_matrix, diags, issparse
from .nambu_moments import NambuMoments


def _real(value, shape, name):
    raw = np.asarray(value)
    if (raw.shape != shape or np.iscomplexobj(raw) or raw.dtype.kind == 'b'
            or np.any(~np.isfinite(raw))):
        raise ValueError(name+' must be finite real with shape '+str(shape))
    return np.asarray(raw, float)


def _indices(value, size, name, *, nonempty=False):
    raw = np.asarray(value)
    if raw.size == 0 and not nonempty:
        return np.empty(0, int)
    if (raw.ndim != 1 or raw.dtype.kind not in 'iu' or not raw.size
            or np.any(raw < 0) or np.any(raw >= size)
            or len(np.unique(raw)) != len(raw)):
        raise ValueError(name+' must contain unique valid integer indices')
    return raw.astype(int, copy=True)


@dataclass(frozen=True)
class JointRowScales:
    neutrality_V: np.ndarray
    gap_force_J_per_gapbar: np.ndarray
    current_A: float


@dataclass(frozen=True)
class JointNambuRowValue:
    neutrality_V: np.ndarray
    gap_force_J_per_gapbar: np.ndarray
    edge_current_A: np.ndarray
    nodal_current_divergence_A: np.ndarray
    interior_current_divergence_A: np.ndarray
    left_current_A: float
    right_current_A: float
    port_current_residual_A: np.ndarray
    port_voltage_V: float
    left_equipotential_residual_V: np.ndarray
    right_equipotential_residual_V: np.ndarray
    thermal_real_axis_discrepancy_A: np.ndarray
    reference_contract: str
    scope: str = 'common moment rows only; rank and constitutive closure not presumed'

    def pack(self, *, neutrality_nodes, force_nodes, scales,
             include_port_currents=True, include_interior_current=False):
        """Explicit row selection; complex force is interleaved Re/Im.

        Returned numbers are dimensionless only after division by the supplied
        physical SI scales. No energy, charge, or current row is discarded by
        a numerical norm. Extra continuity rows require an explicit caller
        choice because conservation can make them algebraically dependent.
        """
        n = len(self.neutrality_V)
        neutral = _indices(neutrality_nodes, n, 'neutrality_nodes')
        force = _indices(force_nodes, n, 'force_nodes')
        charge_scale = _real(scales.neutrality_V, (n,), 'neutrality scale')
        pair_scale = _real(scales.gap_force_J_per_gapbar, (n,), 'gap-force scale')
        current_scale = _real(scales.current_A, (), 'current scale').item()
        if np.any(charge_scale <= 0) or np.any(pair_scale <= 0) or current_scale <= 0:
            raise ValueError('All physical row scales must be positive')
        pairs = self.gap_force_J_per_gapbar[force]/pair_scale[force]
        values = [self.neutrality_V[neutral]/charge_scale[neutral],
                  np.column_stack((pairs.real, pairs.imag)).ravel()]
        if include_port_currents:
            values.append(self.port_current_residual_A/current_scale)
        if include_interior_current:
            values.append(self.interior_current_divergence_A/current_scale)
        return np.concatenate(values)


class JointNambuRows:
    """Same-spectrum neutrality, pair force and thermally completed current.

    Edges point tail->head; the incidence has +1 at tail and -1 at head.
    IL=sum_left(B I), IR=-sum_right(B I), both positive along the bias.
    The cut voltage is mean(phi_left)-mean(phi_right); deviations of all
    individual cut nodes are retained, so an unequipotential cut is visible.
    N0 is the normal DOS PER SPIN, energy is physical J, Eu=kB*Tc, and gapbar
    is Delta/Eu. The thermal reference must be at the instantaneous same
    gap, links and gauge. It is not the old DC state held after a photon.
    """
    def __init__(self, moments: NambuMoments, edges, *, left_nodes, right_nodes):
        if not isinstance(moments, NambuMoments):
            raise TypeError('A NambuMoments instance fixes the SI convention')
        n = len(moments.volume)
        raw = np.asarray(edges)
        if (raw.ndim != 2 or raw.shape[1:] != (2,) or raw.dtype.kind not in 'iu'
                or np.any(raw < 0) or np.any(raw >= n)
                or np.any(raw[:, 0] == raw[:, 1])):
            raise ValueError('Oriented edges must connect distinct valid nodes')
        self.edges = raw.astype(int, copy=True)
        self.left = _indices(left_nodes, n, 'left_nodes', nonempty=True)
        self.right = _indices(right_nodes, n, 'right_nodes', nonempty=True)
        if np.intersect1d(self.left, self.right).size:
            raise ValueError('Longitudinal cuts must be disjoint')
        self.interior = np.setdiff1d(np.arange(n), np.r_[self.left, self.right])
        self.moments = moments
        m = len(self.edges)
        ids = np.arange(m)
        self.incidence = coo_matrix((np.r_[np.ones(m), -np.ones(m)],
            (np.r_[self.edges[:, 0], self.edges[:, 1]], np.r_[ids, ids])),
            shape=(n, m)).tocsr()
        self.scales = JointRowScales(np.full(n, moments.Eu/elementary_charge),
            moments.N0*moments.Eu**2*moments.volume, 2*moments.I0)

    def evaluate(self, *, K_zero, K_first, K_thermal_zero, phi_V,
                 thermal_gradient_J_per_gapbar, thermal_edge_current_A,
                 edge_flux_zero, edge_flux_first, thermal_edge_flux_zero,
                 imposed_current_A, reference_contract):
        """Evaluate common rows, including the independent Matsubara baseline.

        The prescribed Is is an instantaneous circuit state, not a gap-rate
        constraint fabricated from the current defect. Both cuts are reported
        separately; their disagreement must not be hidden by averaging.
        """
        obs = self.moments
        current = _real(imposed_current_A, (), 'imposed_current_A').item()
        phi = _real(phi_V, obs.volume.shape, 'phi_V')
        flux_shape = (len(obs.energy), len(self.edges), 2, 2)
        if np.shape(edge_flux_zero) != flux_shape or np.shape(edge_flux_first) != flux_shape:
            raise ValueError('One full matrix flux per energy and declared edge required')
        values = obs.evaluate(K_zero=K_zero, K_first=K_first,
            K_thermal_zero=K_thermal_zero, phi_V=phi,
            thermal_gradient_J_per_gapbar=thermal_gradient_J_per_gapbar,
            edge_flux_zero=edge_flux_zero, edge_flux_first=edge_flux_first,
            reference_contract=reference_contract)
        completed = obs.complete_current(thermal_edge_current_A=thermal_edge_current_A,
            edge_flux_zero=edge_flux_zero, thermal_edge_flux_zero=thermal_edge_flux_zero,
            edge_flux_first=edge_flux_first, reference_contract=reference_contract)
        divergence = np.asarray(self.incidence@completed.edge_current_A).ravel()
        left = float(np.sum(divergence[self.left]))
        right = -float(np.sum(divergence[self.right]))
        phi_left, phi_right = float(np.mean(phi[self.left])), float(np.mean(phi[self.right]))
        return JointNambuRowValue(values.neutrality_residual_V,
            values.gap_force_J_per_gapbar, completed.edge_current_A, divergence,
            divergence[self.interior], left, right, np.array([left-current, right-current]),
            phi_left-phi_right, phi[self.left]-phi_left, phi[self.right]-phi_right,
            completed.thermal_quadrature_discrepancy_A, reference_contract)

    def linear_direction(self, *, K_first_direction, edge_flux_first_direction,
                  phi_direction_V, reference_contract, K_zero_direction=None,
                  K_thermal_zero_direction=None, edge_flux_zero_direction=None,
                  thermal_edge_flux_zero_direction=None,
                  thermal_gradient_direction_J_per_gapbar=None,
                  thermal_edge_current_direction_A=None, imposed_current_direction_A=0.):
        """Exact affine/Jacobian moment action for arbitrary supplied directions.

        For a rate closure at a fixed state, only K1, J1 and phi vary: the
        default K0 and all independent thermal directions are zero. If a
        Newton method varies the instantaneous state, its K0/thermal tangent
        must be supplied explicitly, not assumed equal to a DC reference.
        This integrates directions; it does not derive a spectrum tangent.
        """
        obs = self.moments
        zero = np.zeros(obs.shape, complex)
        flux = np.zeros((len(obs.energy), len(self.edges), 2, 2), complex)
        return self.evaluate(
            K_zero=zero if K_zero_direction is None else K_zero_direction,
            K_first=K_first_direction,
            K_thermal_zero=zero if K_thermal_zero_direction is None else K_thermal_zero_direction,
            phi_V=phi_direction_V,
            thermal_gradient_J_per_gapbar=np.zeros(len(obs.volume), complex)
                if thermal_gradient_direction_J_per_gapbar is None else thermal_gradient_direction_J_per_gapbar,
            thermal_edge_current_A=np.zeros(len(self.edges))
                if thermal_edge_current_direction_A is None else thermal_edge_current_direction_A,
            edge_flux_zero=flux if edge_flux_zero_direction is None else edge_flux_zero_direction,
            edge_flux_first=edge_flux_first_direction,
            thermal_edge_flux_zero=flux if thermal_edge_flux_zero_direction is None
                else thermal_edge_flux_zero_direction,
            imposed_current_A=imposed_current_direction_A, reference_contract=reference_contract)

    direction = linear_direction


@dataclass(frozen=True)
class GaugeReduction:
    matrix: object
    kept_columns: np.ndarray
    reference_index: int
    gauge_direction: np.ndarray
    scaled_ward_residual: float
    scaled_ward_bound: float

    def expand(self, reduced):
        value = _real(reduced, self.kept_columns.shape, 'reduced unknowns')
        result = np.zeros(len(self.gauge_direction))
        result[self.kept_columns] = value
        return result

    def representative(self, full):
        value = _real(full, self.gauge_direction.shape, 'full unknowns')
        return value-self.gauge_direction*(value[self.reference_index]
            /self.gauge_direction[self.reference_index])


def checked_gauge_reduction(matrix, gauge_direction, *, reference_index,
                            row_scales, column_scales, tolerance=1e-10):
    """Remove ONE pure gauge column only after checking its full Ward vector.

    The caller supplies the actual simultaneous gauge direction, including
    condensate rates, distribution/cut tangents and potential. Setting phi_R
    alone is not a test of gauge covariance. All physical rows are retained;
    the caller decides whether a dependent conservation row is diagnostic.
    A physical nullspace beyond this gauge is never regularized here.
    """
    A = csr_matrix(matrix) if issparse(matrix) else np.asarray(matrix)
    if len(A.shape) != 2 or np.iscomplexobj(A.data if issparse(A) else A):
        raise ValueError('A finite real physical matrix is required')
    if np.any(~np.isfinite(A.data if issparse(A) else A)):
        raise ValueError('A finite real physical matrix is required')
    m, n = A.shape
    g = _real(gauge_direction, (n,), 'gauge_direction')
    rows = _real(row_scales, (m,), 'row_scales')
    columns = _real(column_scales, (n,), 'column_scales')
    tol = _real(tolerance, (), 'tolerance').item()
    if np.any(rows <= 0) or np.any(columns <= 0) or tol <= 0:
        raise ValueError('Positive physical scales and tolerance required')
    if (isinstance(reference_index, (bool, np.bool_)) or not isinstance(reference_index, (int, np.integer))
            or not 0 <= reference_index < n or g[reference_index] == 0):
        raise ValueError('Reference index must have a nonzero gauge component')
    scaled = diags(1/rows)@A@diags(columns) if issparse(A) else A*columns[None, :]/rows[:, None]
    q = g/columns
    ward = np.asarray(scaled@q).ravel()
    bound = np.asarray(abs(scaled)@abs(q)).ravel()
    residual = float(np.max(abs(ward), initial=0.))
    scale = float(np.max(bound, initial=0.))
    if residual > tol*scale:
        raise FloatingPointError('Gauge Ward identity failed; do not freeze potential or remove a physical equation')
    kept = np.delete(np.arange(n), reference_index)
    return GaugeReduction(A[:, kept], kept, int(reference_index), g.copy(), residual, scale)
