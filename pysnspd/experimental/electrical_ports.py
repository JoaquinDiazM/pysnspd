"""Electrical port algebra and the external thesis circuit, in SI units.

This opt-in module does not evolve a superconducting field or select its normal
conductivity. The caller supplies every face conductance G=sigma*S/length and
the superconducting face currents. Edge (i,j) points from i to j; the incidence
matrix has +1 at i and -1 at j, so B I equals the externally injected current.

The circuit implements CM.4 and CM.7--CM.10 in
docs/modelo_v0_4/actualizaciones/circuito_memoria_20260922.md. Its external
inductance is mandatory: the historical total of 10 nH is NOT an external
inductance default. Prescribed-voltage/resistance equilibria below are circuit
oracles, not equilibria or an admission of a coupled detector model.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve


def _real_array(value, name):
    raw = np.asarray(value)
    if np.iscomplexobj(raw):
        raise ValueError(f"{name} must be real")
    result = np.array(raw, dtype=float, copy=True)
    if np.any(~np.isfinite(result)):
        raise ValueError(f"{name} must be finite")
    return result


def _scalar(value, name):
    result = _real_array(value, name)
    if result.ndim != 0:
        raise ValueError(f"{name} must be a scalar")
    return float(result)


def _node(value, size, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise ValueError(f"{name} must be an integer node index")
    if not 0 <= value < size:
        raise ValueError(f"{name} is outside the graph")
    return int(value)


@dataclass(frozen=True)
class PotentialSolution:
    """Passive voltage and powers; residuals are measurements, not admissions.

    With injections only +Is at left and -Is at right, port_power_W is
    Is*Vdev_V. With additional injected ports, it instead includes all of them.
    super_power_W is signed work, while normal_joule_W is dissipated power.
    No absolute-value voltage convention or population/energy repair is used.
    """
    phi_V: np.ndarray
    delta_phi_V: np.ndarray
    I_normal_A: np.ndarray
    I_total_A: np.ndarray
    Vdev_V: float
    current_residual_A: np.ndarray
    injection_imbalance_A: float
    port_power_W: float
    face_power_W: float
    normal_joule_W: float
    super_power_W: float
    power_residual_W: float
    decomposition_residual_W: float
    reference_node: int
    reference_potential_V: float


def solve_potential(edges, conductance_S, supercurrent_A, injection_A, *,
                    left_node, right_node, reference_node=None,
                    reference_potential_V=0.):
    """Solve B G B.T phi = injection - B I_super on a connected graph.

    ``injection_A`` defines the number of nodes and is positive into a node
    from the external circuit. ``conductance_S`` and ``supercurrent_A`` contain
    one entry per directed edge. Parallel edges are allowed, self edges are
    rejected. The reference defaults to the right terminal with phi_right=0.
    A different constant reference tests electrical gauge invariance.

    Input current balance is checked to floating point roundoff relative to
    the input currents; no current is changed to enforce it. The returned
    nodal residual includes the reference equation omitted in the linear solve.
    """
    injection = _real_array(injection_A, "injection_A")
    if injection.ndim != 1 or len(injection) < 2:
        raise ValueError("injection_A must contain at least two nodes")
    nodes = len(injection)
    raw_edges = np.asarray(edges)
    if (raw_edges.ndim != 2 or raw_edges.shape[1] != 2 or len(raw_edges) == 0
            or raw_edges.dtype.kind not in "iu"):
        raise ValueError("edges must be a nonempty (faces,2) integer array")
    if np.any(raw_edges < 0) or np.any(raw_edges >= nodes):
        raise ValueError("edge index is outside the graph")
    edge = raw_edges.astype(np.intp, copy=True)
    tail, head = edge.T
    if np.any(tail == head):
        raise ValueError("self edges are not electrical faces")
    conductance = _real_array(conductance_S, "conductance_S")
    supercurrent = _real_array(supercurrent_A, "supercurrent_A")
    if conductance.shape != (len(edge),) or supercurrent.shape != (len(edge),):
        raise ValueError("one conductance and superconducting current per edge required")
    if np.any(conductance <= 0):
        raise ValueError("conductance_S must be strictly positive")
    left = _node(left_node, nodes, "left_node")
    right = _node(right_node, nodes, "right_node")
    if left == right:
        raise ValueError("the two device terminals must be distinct")
    reference = right if reference_node is None else _node(reference_node, nodes, "reference_node")
    reference_value = _scalar(reference_potential_V, "reference_potential_V")
    # Connectivity is topological and is checked independently of solver warnings.
    neighbours = [[] for _ in range(nodes)]
    for i, j in edge:
        neighbours[i].append(j)
        neighbours[j].append(i)
    reached = {reference}
    pending = [reference]
    while pending:
        for neighbour in neighbours[pending.pop()]:
            if neighbour not in reached:
                reached.add(neighbour)
                pending.append(neighbour)
    if len(reached) != nodes:
        raise ValueError("electrical graph must be connected")
    imbalance = math.fsum(injection)
    scale = math.fsum(abs(injection))
    if abs(imbalance) > 64*np.finfo(float).eps*scale:
        raise ValueError("external current injections do not balance")
    faces = np.arange(len(edge))
    B = coo_matrix((np.r_[np.ones(len(edge)), -np.ones(len(edge))],
                    (np.r_[tail, head], np.r_[faces, faces])),
                   shape=(nodes, len(edge))).tocsr()
    laplacian = (B.multiply(conductance) @ B.T).tocsr()
    rhs = injection - B @ supercurrent
    unknown = np.arange(nodes) != reference
    relative_phi = np.zeros(nodes)
    relative_phi[unknown] = spsolve(laplacian[unknown][:, unknown], rhs[unknown])
    phi = relative_phi + reference_value
    # Use the relative solution for differences to avoid subtracting a large
    # irrelevant reference potential from neighbouring potentials.
    drops = relative_phi[tail] - relative_phi[head]
    normal = conductance*drops
    total = supercurrent + normal
    residual = np.asarray(B @ total - injection)
    port_power = float(np.dot(injection, phi))
    face_power = float(np.dot(drops, total))
    joule = float(np.dot(conductance, drops*drops))
    super_power = float(np.dot(drops, supercurrent))
    vdev = float(relative_phi[left] - relative_phi[right])
    if (not all(np.all(np.isfinite(v)) for v in (phi, drops, normal, total, residual))
            or not all(np.isfinite(v) for v in (port_power, face_power, joule, super_power, vdev))):
        raise FloatingPointError("nonfinite electrical solution or power")
    return PotentialSolution(phi, drops, normal, total, vdev, residual, imbalance,
        port_power, face_power, joule, super_power, port_power-face_power,
        face_power-joule-super_power, reference, reference_value)


@dataclass(frozen=True)
class ThesisCircuitState:
    """Circuit state CM.1: bias current, total detector current, capacitor voltage."""
    Ib_A: float
    Is_A: float
    vc_V: float

    def as_array(self):
        return _real_array([self.Ib_A, self.Is_A, self.vc_V], "circuit state")


def _state(y):
    result = y.as_array() if isinstance(y, ThesisCircuitState) else _real_array(y, "circuit state")
    if result.shape != (3,):
        raise ValueError("circuit state must have shape (3,) in order (Ib,Is,vc)")
    return result


@dataclass(frozen=True)
class CircuitOutputs:
    Ib_A: float
    Is_A: float
    I_RF_A: float
    vc_V: float
    Vd_V: float
    Vdev_V: float
    Vout_V: float


@dataclass(frozen=True)
class CircuitPowerBalance:
    stored_energy_rate_W: float
    device_power_W: float
    bias_joule_W: float
    load_joule_W: float
    source_power_W: float
    residual_W: float


@dataclass(frozen=True)
class ThesisCircuitParameters:
    """Constant external network of CM.4; all parameters are SI.

    Lk_ext_H must describe only material outside the resolved energy domain.
    Its positivity is necessary, not proof of a correct inductance partition.
    Other defaults are the documented historical circuit values. A circuit
    state may be a length-three array or ``ThesisCircuitState``.
    """
    Lk_ext_H: float
    R_bias_ohm: float = 1.e4
    L_bias_H: float = 1.e-6
    R_load_ohm: float = 50.
    C_couple_F: float = 100.e-12
    V_bias_V: float = .300

    def __post_init__(self):
        for name in ("Lk_ext_H", "R_bias_ohm", "L_bias_H", "R_load_ohm", "C_couple_F", "V_bias_V"):
            value = _scalar(getattr(self, name), name)
            if name != "V_bias_V" and value <= 0:
                raise ValueError(f"{name} must be strictly positive")
            object.__setattr__(self, name, value)

    def rhs(self, y, Vdev_V):
        """CM.4 derivatives (A/s, A/s, V/s) for a supplied passive port voltage."""
        Ib, Is, vc = _state(y)
        voltage = _scalar(Vdev_V, "Vdev_V")
        irf = Ib-Is
        return np.array([(self.V_bias_V-self.R_bias_ohm*Ib-vc-self.R_load_ohm*irf)/self.L_bias_H,
                         (vc+self.R_load_ohm*irf-voltage)/self.Lk_ext_H,
                         irf/self.C_couple_F])

    def energy(self, y):
        """External inductive and capacitive storage CM.7, in joules."""
        Ib, Is, vc = _state(y)
        return float(.5*(self.L_bias_H*Ib*Ib+self.Lk_ext_H*Is*Is+self.C_couple_F*vc*vc))

    def outputs(self, y, Vdev_V):
        Ib, Is, vc = _state(y)
        voltage = _scalar(Vdev_V, "Vdev_V")
        irf = float(Ib-Is)
        vout = self.R_load_ohm*irf
        return CircuitOutputs(float(Ib), float(Is), irf, float(vc), float(vc+vout), voltage, vout)

    def power_balance(self, y, Vdev_V, *, derivative=None):
        """CM.8 instantaneous residual, optionally for a supplied state derivative.

        This is not an integrated energy ledger or a coupled time-admission test.
        The device power is Is*Vdev, not the bias current times the readout voltage.
        """
        state = _state(y)
        voltage = _scalar(Vdev_V, "Vdev_V")
        dy = self.rhs(state, voltage) if derivative is None else _state(derivative)
        Ib, Is, vc = state
        stored = float(np.dot([self.L_bias_H*Ib, self.Lk_ext_H*Is, self.C_couple_F*vc], dy))
        device = float(Is*voltage)
        bias = float(self.R_bias_ohm*Ib*Ib)
        load = float(self.R_load_ohm*(Ib-Is)**2)
        source = float(self.V_bias_V*Ib)
        return CircuitPowerBalance(stored, device, bias, load, source,
                                   stored+device+bias+load-source)

    def stationary_for_prescribed_voltage(self, Vdev_V):
        """CM.10 circuit oracle for an externally prescribed constant voltage."""
        voltage = _scalar(Vdev_V, "Vdev_V")
        current = (self.V_bias_V-voltage)/self.R_bias_ohm
        return np.array([current, current, voltage])

    def stationary_for_prescribed_resistance(self, Rdev_ohm):
        """Circuit-only oracle with Vdev=Rdev*Is; no film is replaced implicitly."""
        resistance = _scalar(Rdev_ohm, "Rdev_ohm")
        if resistance < 0:
            raise ValueError("prescribed device resistance must be nonnegative")
        current = self.V_bias_V/(self.R_bias_ohm+resistance)
        return np.array([current, current, resistance*current])
