"""Effective quasineutral electrical ports for a conserving gap model.

The normal network is the inherited sheet conductance on the SAME graph.
Each complete longitudinal cut is equipotential; its injected current is an
unknown distribution, while its sum is the imposed detector current Is.
Supercurrent must be the exact gauge derivative of the shared energy:
  I_e = -(2e/hbar) partial U / partial alpha_e.
No separately sampled continuum current or dynamic lumped inductance is used.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.constants import elementary_charge as e, hbar
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import splu

from .electrical_ports import ThesisCircuitParameters, CircuitOutputs, CircuitPowerBalance
from .thermal_spatial_usadel import ThermalGraph


def _real(value, shape, name):
    raw = np.asarray(value)
    if raw.shape != shape or np.iscomplexobj(raw) or raw.dtype.kind == "b" or np.any(~np.isfinite(raw)):
        raise ValueError(name + " requires finite real values of shape " + str(shape))
    return np.asarray(raw, float)


def _scalar(value, name):
    return float(_real(value, (), name))


def _nodes(value, n, name):
    raw = np.asarray(value)
    if (raw.ndim != 1 or raw.dtype.kind not in "iu" or not raw.size
            or np.any(raw < 0) or np.any(raw >= n) or len(np.unique(raw)) != raw.size):
        raise ValueError(name + " requires distinct valid integer nodes")
    return np.sort(np.asarray(raw, int))


@dataclass(frozen=True)
class EffectivePotential:
    phi_V: np.ndarray
    delta_phi_V: np.ndarray
    I_normal_A: np.ndarray
    I_super_A: np.ndarray
    I_total_A: np.ndarray
    Vdev_V: float
    left_current_A: float
    right_current_A: float
    port_current_residual_A: np.ndarray
    interior_current_residual_A: np.ndarray
    nodal_boundary_injection_A: np.ndarray
    nodal_joule_W: np.ndarray
    face_joule_W: np.ndarray
    normal_joule_W: float
    super_power_W: float
    boundary_power_W: float
    face_power_W: float
    power_residual_W: float
    decomposition_residual_W: float
    joule_allocation_residual_W: float
    reference_potential_V: float

    @property
    def imposed_port_power_W(self):
        # The imposed current is reconstructed independently of the reaction.
        return float((self.left_current_A-self.port_current_residual_A[0])*self.Vdev_V)

    @property
    def port_power_W(self):
        return self.imposed_port_power_W

    @property
    def boundary_current_power_residual_W(self):
        return self.boundary_power_W-self.imposed_port_power_W


@dataclass(frozen=True)
class EffectiveCircuitPort:
    electrical: EffectivePotential
    circuit_rate_SI_per_s: np.ndarray
    circuit_outputs: CircuitOutputs
    circuit_power: CircuitPowerBalance
    circuit_energy_J: float
    baseline_subtracted_Vout_V: float
    global_current_direction_ok: bool


@dataclass(frozen=True)
class EffectiveGaugePower:
    noether_residual_J: np.ndarray
    noether_residual_max_J: float
    field_phase_power_W: float
    edge_superconducting_power_W: float
    power_residual_W: float


@dataclass(frozen=True)
class FrozenExteriorInductance:
    total_H: float
    added_inductance_H: float
    unresolved_bulk_inductance_H: float
    bulk_differential_inductance_per_length_H_per_m: float
    total_active_length_m: float
    resolved_length_m: float
    unresolved_length_m: float
    scope: str = "Frozen once from the chosen DC branch; never updated during a pulse"


def freeze_exterior_inductance(bulk_differential_inductance_per_length_H_per_m,
        resolved_length_m, *, total_active_length_m=5e-6, added_inductance_H=96e-9):
    """Partition only the 5 um wire; the added 96 nH is never reduced."""
    per_length = _scalar(bulk_differential_inductance_per_length_H_per_m, "bulk differential inductance")
    resolved = _scalar(resolved_length_m, "resolved_length_m")
    total = _scalar(total_active_length_m, "total_active_length_m")
    added = _scalar(added_inductance_H, "added_inductance_H")
    if per_length <= 0 or total <= 0 or added <= 0 or resolved <= 0 or resolved > total:
        raise ValueError("Positive inductance and 0 < resolved_length <= total_active_length required")
    length = total-resolved
    bulk = per_length*length
    return FrozenExteriorInductance(added+bulk, added, bulk, per_length, total, resolved, length)


class EffectiveNeutralPorts:
    """Cached SI potential solve with two complete equipotential strip cuts.

    The quasineutral approximation eliminates fast charge storage. It preserves
    the full external three-state thesis circuit. Longitudinal cut reactions
    are returned individually: uniform per-node injection is never imposed.
    A fixed normal conductance is an explicit effective-model assumption.
    """
    def __init__(self, graph: ThermalGraph, sheet_resistance_ohm, left_nodes, right_nodes,
            *, longitudinal_axis=0, validate_full_width=True):
        if not isinstance(graph, ThermalGraph):
            raise TypeError("A ThermalGraph fixes geometry and edge orientation")
        self.graph = graph
        self.sheet_resistance_ohm = _scalar(sheet_resistance_ohm, "sheet_resistance_ohm")
        if self.sheet_resistance_ohm <= 0:
            raise ValueError("Positive sheet resistance required")
        n = graph.n_nodes
        self.left = _nodes(left_nodes, n, "left_nodes")
        self.right = _nodes(right_nodes, n, "right_nodes")
        if np.intersect1d(self.left, self.right).size:
            raise ValueError("The two strip cuts must be disjoint")
        if longitudinal_axis not in (0, 1):
            raise ValueError("longitudinal_axis must be 0 or 1")
        if validate_full_width:
            if graph.coordinates_bar is None:
                raise ValueError("Full-width cut validation needs graph coordinates")
            coordinate = graph.coordinates_bar[:, longitudinal_axis]
            lo, hi = float(np.min(coordinate)), float(np.max(coordinate))
            tolerance = 256*np.finfo(float).eps*max(1., abs(lo), abs(hi), hi-lo)
            if hi-lo <= tolerance:
                raise ValueError("A nonzero longitudinal length is required")
            expected_left = np.flatnonzero(abs(coordinate-lo) <= tolerance)
            expected_right = np.flatnonzero(abs(coordinate-hi) <= tolerance)
            if not np.array_equal(self.left, expected_left) or not np.array_equal(self.right, expected_right):
                raise ValueError("Electrical cuts must contain every node at each longitudinal extremum")
        self.free = np.setdiff1d(np.arange(n), np.r_[self.left, self.right])
        labels = np.empty(n, int)
        labels[self.free] = np.arange(len(self.free))
        labels[self.left] = len(self.free)
        labels[self.right] = len(self.free)+1
        self.labels = labels
        self.node_to_group = labels
        self.group_count = len(self.free)+2
        self.prolongation = coo_matrix((np.ones(n), (np.arange(n), labels)), shape=(n, self.group_count)).tocsr()
        self.left_group, self.right_group = len(self.free), len(self.free)+1
        ng = len(self.free)+2
        tail, head = graph.edges.T
        ids = np.arange(len(tail))
        self.incidence = coo_matrix((np.r_[np.ones(len(ids)), -np.ones(len(ids))],
            (np.r_[tail, head], np.r_[ids, ids])), shape=(n, len(ids))).tocsr()
        self.reduced_incidence = coo_matrix((np.r_[np.ones(len(ids)), -np.ones(len(ids))],
            (np.r_[labels[tail], labels[head]], np.r_[ids, ids])), shape=(ng, len(ids))).tocsr()
        self.reduced_incidence.eliminate_zeros()
        self.conductance_S = np.asarray(graph.conductance/self.sheet_resistance_ohm, float)
        self.laplacian = (self.reduced_incidence.multiply(self.conductance_S) @ self.reduced_incidence.T).tocsr()
        # Exposed grouped SI operator for a coupled implicit block solve.
        self.normal_laplacian = self.laplacian
        self.reference_group = self.right_group
        self.group_injection_per_A = np.zeros(ng)
        self.group_injection_per_A[self.left_group] = 1.
        self.group_injection_per_A[self.right_group] = -1.
        count, _ = connected_components(self.laplacian, directed=False)
        if count != 1:
            raise ValueError("Positive-conductance electrical graph must connect both cuts and every interior node")
        self.unknown = np.arange(ng) != self.right_group
        self.factor = splu(self.laplacian[self.unknown][:, self.unknown].tocsc())

    def solve(self, edge_supercurrent_A, Is_A, *, reference_potential_V=0.):
        """I=I_super+G*(phi_tail-phi_head), B I=0 in the interior."""
        supercurrent = _real(edge_supercurrent_A, (len(self.graph.edges),), "edge_supercurrent_A")
        current = _scalar(Is_A, "Is_A")
        reference = _scalar(reference_potential_V, "reference_potential_V")
        source = np.zeros(self.laplacian.shape[0])
        source[self.left_group], source[self.right_group] = current, -current
        rhs = source-self.reduced_incidence@supercurrent
        relative_group = np.zeros_like(source)
        relative_group[self.unknown] = self.factor.solve(rhs[self.unknown])
        relative = relative_group[self.labels]
        phi = relative+reference
        tail, head = self.graph.edges.T
        drop = relative[tail]-relative[head]
        normal = self.conductance_S*drop
        total = supercurrent+normal
        divergence = np.asarray(self.incidence@total).ravel()
        left_current, right_current = float(np.sum(divergence[self.left])), -float(np.sum(divergence[self.right]))
        reactions = np.zeros(self.graph.n_nodes)
        reactions[self.left], reactions[self.right] = divergence[self.left], divergence[self.right]
        face_joule = self.conductance_S*drop*drop
        nodal_joule = np.zeros(self.graph.n_nodes)
        np.add.at(nodal_joule, tail, .5*face_joule)
        np.add.at(nodal_joule, head, .5*face_joule)
        vdev = float(relative_group[self.left_group])
        # This uses independently reconstructed cut reactions, not Is*Vdev.
        boundary_power = float(np.dot(relative, reactions))
        face_power = float(np.dot(drop, total))
        joule, super_power = float(np.sum(face_joule)), float(np.dot(drop, supercurrent))
        if not all(np.all(np.isfinite(v)) for v in (phi, normal, total, divergence, nodal_joule)):
            raise FloatingPointError("Nonfinite potential/current/power; no correction applied")
        return EffectivePotential(phi, drop, normal, supercurrent.copy(), total, vdev,
            left_current, right_current, np.array([left_current-current, right_current-current]),
            divergence[self.free], reactions, nodal_joule, face_joule, joule,
            super_power, boundary_power, face_power, boundary_power-face_power,
            face_power-joule-super_power, float(np.sum(nodal_joule)-joule), reference)

    def coupled(self, edge_supercurrent_A, circuit: ThesisCircuitParameters, state,
            *, baseline_Vout_V=0., reference_potential_V=0., require_forward_current=True):
        """Same instantaneous Is, device voltage and complete inherited CM RHS."""
        raw = state.as_array() if hasattr(state, "as_array") else state
        y = _real(raw, (3,), "circuit state (Ib,Is,vc)")
        direction_ok = bool(y[0] > 0 and y[1] > 0)
        if require_forward_current and not direction_ok:
            raise ValueError("Global bias/detector current reversal is outside the admitted observation window")
        electrical = self.solve(edge_supercurrent_A, y[1], reference_potential_V=reference_potential_V)
        derivative = circuit.rhs(y, electrical.Vdev_V)
        outputs = circuit.outputs(y, electrical.Vdev_V)
        return EffectiveCircuitPort(electrical, derivative, outputs,
            circuit.power_balance(y, electrical.Vdev_V, derivative=derivative),
            circuit.energy(y), outputs.Vout_V-_scalar(baseline_Vout_V, "baseline_Vout_V"), direction_ok)

    def gauge_power_check(self, delta_bar, full_gradient_J_per_gapbar, edge_supercurrent_A,
            phi_V):
        """Verify SI Ward identity using a full real Cartesian energy gradient.

        For a Wirtinger force F returned by an energy functional supply 2*F.
        This independent identity is not enforced or repaired by the solver.
        """
        n = self.graph.n_nodes
        delta = np.asarray(delta_bar, complex)
        gradient = np.asarray(full_gradient_J_per_gapbar, complex)
        phi = _real(phi_V, (n,), "phi_V")
        current = _real(edge_supercurrent_A, (len(self.graph.edges),), "edge_supercurrent_A")
        if delta.shape != (n,) or gradient.shape != (n,) or np.any(~np.isfinite(delta)) or np.any(~np.isfinite(gradient)):
            raise ValueError("Finite gap and full energy gradient per graph node required")
        phase_gradient = np.imag(delta.conj()*gradient)
        divergence = np.asarray(self.incidence@current).ravel()
        noether = phase_gradient+hbar/(2*e)*divergence
        phase_velocity = -2j*e/hbar*phi*delta
        field_power = float(np.real(np.vdot(gradient, phase_velocity)))
        tail, head = self.graph.edges.T
        edge_power = float(np.dot(current, phi[tail]-phi[head]))
        return EffectiveGaugePower(noether, float(np.max(abs(noether), initial=0.)),
            field_power, edge_power, field_power-edge_power)


@dataclass(frozen=True)
class EffectiveBoundaryLoad:
    full_load_J_per_gapbar: np.ndarray
    radial_load_J_per_gapbar: np.ndarray
    phase_load_J_per_gapbar: np.ndarray
    nodal_target_injection_A: np.ndarray
    left_width_weights: np.ndarray
    right_width_weights: np.ndarray
    phase_torque_residual_J: np.ndarray
    scope: str = "Natural bulk traction; its signed material work is external, never heat"


def boundary_material_power_W(load_J_per_gapbar, covariant_gap_rate_per_s):
    """Independent conjugate work Re sum(conj(load) D_t delta), in watts."""
    load = np.asarray(load_J_per_gapbar, complex)
    rate = np.asarray(covariant_gap_rate_per_s, complex)
    if load.ndim != 1 or rate.shape != load.shape or np.any(~np.isfinite(load)) or np.any(~np.isfinite(rate)):
        raise ValueError("Finite matching boundary loads and covariant gap rates required")
    return float(np.real(np.vdot(load, rate)))


def bulk_boundary_load(ports: EffectiveNeutralPorts, delta_bar, bulk_delta_bar,
        bulk_full_gradient_J_per_gapbar, Is_A, *, left_width_weights=None,
        right_width_weights=None):
    """Continue the shared bulk functional without fixing an interior profile.

    Radial traction is taken from the SAME discrete functional on the chosen
    uniform current-carrying bulk branch. Phase traction prescribes only each
    cut's total current, using declared width weights. Default weights are the
    cut dual areas, proportional to width for the admitted straight strip.
    A caller evolving the cut material phase must retain the resulting signed
    boundary material work in addition to electrical phi-drop work.
    """
    if not isinstance(ports, EffectiveNeutralPorts):
        raise TypeError("EffectiveNeutralPorts fixes the two full-width cuts")
    n = ports.graph.n_nodes
    delta, bulk, gradient = [np.asarray(v, complex) for v in
        (delta_bar, bulk_delta_bar, bulk_full_gradient_J_per_gapbar)]
    if any(v.shape != (n,) or np.any(~np.isfinite(v)) for v in (delta, bulk, gradient)):
        raise ValueError("Finite actual/bulk gap and full bulk gradient per node required")
    fixed = np.r_[ports.left, ports.right]
    if np.any(abs(delta[fixed]) == 0) or np.any(abs(bulk[fixed]) == 0):
        raise ValueError("Superconducting continuation requires nonzero cut gap")
    current = _scalar(Is_A, "Is_A")
    def weights(value, nodes, name):
        if value is None:
            area = ports.graph.area_weights[nodes]
            return np.asarray(area/np.sum(area))
        w = _real(value, (len(nodes),), name)
        if np.any(w <= 0) or abs(float(np.sum(w))-1.) > 64*np.finfo(float).eps*len(w):
            raise ValueError("Declared positive cut weights must sum to one; no normalization applied")
        return w
    wl = weights(left_width_weights, ports.left, "left_width_weights")
    wr = weights(right_width_weights, ports.right, "right_width_weights")
    injection = np.zeros(n)
    injection[ports.left], injection[ports.right] = current*wl, -current*wr
    direction = delta[fixed]/abs(delta[fixed])
    bulk_direction = bulk[fixed]/abs(bulk[fixed])
    radial = np.zeros(n, complex)
    radial[fixed] = direction*np.real(bulk_direction.conj()*gradient[fixed])
    torque = -hbar/(2*e)*injection
    phase = np.zeros(n, complex)
    phase[fixed] = 1j*direction*torque[fixed]/abs(delta[fixed])
    load = radial+phase
    check = np.imag(delta.conj()*load)-torque
    return EffectiveBoundaryLoad(load, radial, phase, injection, wl, wr, check)
