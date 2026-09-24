"""Bulk end ports for the isothermal graph KWT/DC preparation diagnostic.

Each complete longitudinal end is equipotential. Contracting its nodes leaves
the spatial mesh and all spectral fields untouched; only the electrical solve
uses the quotient graph. Parallel faces retain their original conductances.
The normal Ohmic law is the inherited thermal closure, not a nonequilibrium
kinetic constitutive law.
"""
from dataclasses import dataclass

import numpy as np
from scipy.constants import Boltzmann, elementary_charge

from .electrical_ports import solve_potential
from .thermal_weak_response import ThermalKWTNormal


@dataclass(frozen=True)
class DCVoltage:
    phi_V: np.ndarray
    potential_v: np.ndarray
    normal_current_bar: np.ndarray
    Vdev_V: float
    port_power_W: float
    normal_joule_W: float
    current_residual_A: np.ndarray
    left_nodes: np.ndarray
    right_nodes: np.ndarray


class DCCurrentPorts:
    """Passive orientation: injection +I at x_min, -I at x_max; V=phiL-phiR."""
    def __init__(self, graph, *, Tc_K, sheet_resistance_ohm):
        if any(not np.isfinite(v) or v <= 0 for v in (Tc_K, sheet_resistance_ohm)):
            raise ValueError('Positive Tc and sheet resistance required')
        self.graph = graph
        self.current_scale_A = Boltzmann*Tc_K/(2*elementary_charge*sheet_resistance_ohm)
        self.voltage_scale_V = Boltzmann*Tc_K/(2*elementary_charge)
        x = graph.coordinates_bar[:, 0]
        self.left = np.flatnonzero(np.isclose(x, x.min(), rtol=0, atol=1e-10))
        self.right = np.flatnonzero(np.isclose(x, x.max(), rtol=0, atol=1e-10))
        if not np.array_equal(np.sort(graph.boundary_nodes), np.sort(np.r_[self.left, self.right])):
            raise ValueError('Declared fixed spectral/gap nodes must be the two longitudinal ends')
        interior = np.setdiff1d(np.arange(graph.n_nodes), graph.boundary_nodes)
        mapping = np.empty(graph.n_nodes, dtype=int)
        mapping[self.left] = 0; mapping[self.right] = 1
        mapping[interior] = np.arange(2, 2+len(interior))
        self.mapping = mapping
        quotient_edges = mapping[graph.edges]
        self.keep = ((quotient_edges[:, 0] != quotient_edges[:, 1]) & (graph.conductance > 0))
        self.edges = quotient_edges[self.keep]
        self.conductance_S = graph.conductance[self.keep]/sheet_resistance_ohm
        self.n_nodes = len(interior)+2

    def solve(self, supercurrent_bar, imposed_current_A):
        current = np.asarray(supercurrent_bar, float)
        if current.shape != (len(self.graph.edges),) or np.any(~np.isfinite(current)):
            raise ValueError('One finite superconducting current per original mesh edge required')
        if not np.isfinite(imposed_current_A):
            raise ValueError('Finite imposed circuit current required')
        injection = np.zeros(self.n_nodes)
        injection[0] = imposed_current_A; injection[1] = -imposed_current_A
        result = solve_potential(self.edges, self.conductance_S,
            current[self.keep]*self.current_scale_A, injection, left_node=0, right_node=1)
        phi = result.phi_V[self.mapping]
        potential = phi/self.voltage_scale_V
        tail, head = self.graph.edges.T
        normal = self.graph.conductance*(potential[tail]-potential[head])
        return DCVoltage(phi, potential, normal, result.Vdev_V,
            result.port_power_W, result.normal_joule_W, result.current_residual_A,
            self.left, self.right)


class CurrentDrivenThermalKWT(ThermalKWTNormal):
    """Reuse the KWT mobility, replacing only its port-current potential solve."""
    def __init__(self, graph, d, *, ports, imposed_current_A, **kwargs):
        super().__init__(graph, d, **kwargs)
        self.ports = ports
        self.imposed_current_A = float(imposed_current_A)
        self.last_voltage = None
        self.last_current = None

    def potential(self, current):
        if self.last_current is None or not np.array_equal(current, self.last_current):
            self.last_voltage = self.ports.solve(current, self.imposed_current_A)
            self.last_current = np.asarray(current).copy()
        return self.last_voltage.potential_v


def advance_bulk_contact_phases(gap, potential_v, fixed_nodes, *, dt_ps, tD_ps):
    """Josephson gauge velocity: dtheta/d(t/tD)=-v/2, with v=2e phi/(kBTc)."""
    result = np.asarray(gap, complex).copy()
    if not np.isfinite(dt_ps) or dt_ps <= 0 or not np.isfinite(tD_ps) or tD_ps <= 0:
        raise ValueError('Positive finite physical and material time scales required')
    result[fixed_nodes] *= np.exp(-.5j*np.asarray(potential_v)[fixed_nodes]*dt_ps/tD_ps)
    return result
