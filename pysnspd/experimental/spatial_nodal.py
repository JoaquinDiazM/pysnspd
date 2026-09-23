"""Fourth-order nodal spatial energy, without changing the archived 3A scheme.

Let C transport a nearest neighbour with exp(i*link_phase). With X=x/ell0,
    D4 = [8(C-C*)-(C**2-C* **2)]/(12h),
    L4 = [30I-16(C+C*)+(C**2+C* **2)]/(12h**2).
Here C* means the adjoint, L4 is positive, and D4 is skew-adjoint. The energy
is h*sum[u(|z|,q_delta**2/r,p)-kappa*|z|**2*q_delta**2]
plus kappa*h*Re(z* L4 z). All forces and currents differentiate this energy.
The local amplitude is |z_i|, independent of a change in the gauge links.

The compatible L4 term is essential: D4 alone has a checkerboard null mode.
For t=sin(kh/2)**2, L4-D4*D4 has symbol 16*t**3*(2+t)/(9h**2)>=0.
Together with rho*q_delta**2<=|D4z|**2 this makes the integrated gradient
remainder nonnegative. Unitary covariant shifts preserve this inequality.

Populations are nodal and fixed during variations. This is a different
spatial representation from midpoint p_links, not an implicit remapping of
those populations. Physical units, regularization, catalogue support and the
continuous local D.36 diagnostic are inherited without modification.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from .energy_catalog import E_CHARGE_C, HBAR_J_S
from .spatial_functional import (
    PeriodicSpatialFunctional, PrincipalSymbol, SpatialAdmissibilityError)


@dataclass(frozen=True)
class NodalSpatialEvaluation:
    energy_bar: float
    energy_J: float
    gradient_cartesian_bar: np.ndarray
    wirtinger_force_bar: np.ndarray
    amplitude_gradient_bar: np.ndarray | None
    link_current_bar: np.ndarray
    current_A: np.ndarray
    current_density_A_m2: np.ndarray
    amplitude_nodes_bar: np.ndarray
    gamma_nodes_bar: np.ndarray
    q_delta_nodes_bar: np.ndarray
    # Nodal bookkeeping density; only its sum is the integrated energy.
    energy_density_bar: np.ndarray
    electronic_derivatives: np.ndarray
    principal_symbols: tuple[PrincipalSymbol, ...] | None
    covariant_derivative_bar: np.ndarray
    gradient_remainder_bar: float


class PeriodicNodalSpatialFunctional(PeriodicSpatialFunctional):
    """One catalogue query per node with fourth-order covariant operators.

    The nearest-neighbour current includes contributions from every length-two
    path through that link. No effective next-neighbour current replaces it.
    ``require_stability=False`` skips local D.36 checks, not domain checks.
    Passing D.36 remains a local continuous-symbol check; it does not certify
    every discrete perturbation, resolution or nonlinear trajectory.
    """

    scheme = "nodal_D4_L4_covariant_energy_v1"

    def _state_and_links(self, delta_bar, twist, link_phases):
        z = np.asarray(delta_bar, dtype=complex)
        if z.shape != (self.cells,) or np.any(~np.isfinite(z)):
            raise ValueError("delta_bar must contain one finite complex value per node")
        if not np.isfinite(twist):
            raise ValueError("twist must be finite")
        if link_phases is None:
            phases = np.zeros(self.cells)
        else:
            raw = np.asarray(link_phases)
            if np.iscomplexobj(raw):
                raise ValueError("link_phases must be real")
            phases = np.array(raw, dtype=float, copy=True)
        if phases.shape != (self.cells,) or np.any(~np.isfinite(phases)):
            raise ValueError("link_phases must have one finite real value per link")
        phases[-1] += twist
        return z, np.exp(1j*phases)

    @staticmethod
    def _shift(z, links):
        return links*np.roll(z, -1)

    @staticmethod
    def _adjoint_shift(z, links):
        return np.roll(np.conj(links)*z, 1)

    def _operators(self, z, links):
        forward = self._shift(z, links)
        backward = self._adjoint_shift(z, links)
        forward2 = self._shift(forward, links)
        backward2 = self._adjoint_shift(backward, links)
        derivative = (8*(forward-backward)-(forward2-backward2))/(12*self.h_bar)
        positive_laplacian = (30*z-16*(forward+backward)+forward2+backward2)/(12*self.h_bar**2)
        return derivative, positive_laplacian

    def covariant_derivative(self, delta_bar, *, twist=0., link_phases=None):
        """Return D4 Delta in the inherited dimensionless coordinate X."""
        z, links = self._state_and_links(delta_bar, twist, link_phases)
        return self._operators(z, links)[0]

    def positive_laplacian(self, delta_bar, *, twist=0., link_phases=None):
        """Return L4 Delta, the positive approximation to -D_X**2 Delta."""
        z, links = self._state_and_links(delta_bar, twist, link_phases)
        return self._operators(z, links)[1]

    def _link_current(self, z, gd, links):
        """Differentiate D4 and L4, depositing each path on its original links."""
        current = np.zeros(self.cells)
        path = np.ones(self.cells, dtype=complex)
        # D4=sum c_j(C**j-C* **j); energy stiffness=sum a_j|C**j z-z|²/h².
        for hop, c, a in ((1, 2/(3*self.h_bar), 4/3),
                          (2, -1/(12*self.h_bar), -1/12)):
            path *= np.roll(links, -(hop-1))
            transported = path*np.roll(z, -hop)
            endpoint_gd = np.roll(gd, -hop)
            derivative_work = self.h_bar*c*np.real(1j*(
                np.conj(gd)*transported+np.conj(endpoint_gd)*np.conj(path)*z))
            stiffness_work = 2*self.kappa*a/self.h_bar*np.imag(np.conj(z)*transported)
            path_current = derivative_work+stiffness_work
            for edge_in_path in range(hop):
                current += np.roll(path_current, edge_in_path)
        return current

    def evaluate(self, delta_bar, p_nodes, *, twist: float = 0., link_phases=None,
                 require_stability: bool = True,
                 on_node: Callable[[int, int], None] | None = None,
                 on_cell: Callable[[int, int], None] | None = None) -> NodalSpatialEvaluation:
        """Evaluate the nodal energy, Cartesian force and nearest-link current.

        Twist and current signs match the archived module. ``p_nodes`` has
        shape (cells,count_states). A callback runs after each node, with its
        zero-based index and total node count. ``on_cell`` is an exclusive
        compatibility alias; supplying both callbacks is an error.
        """
        z, links = self._state_and_links(delta_bar, twist, link_phases)
        raw = np.asarray(p_nodes)
        if np.iscomplexobj(raw):
            raise ValueError("p_nodes must be real")
        p = np.asarray(raw, dtype=float)
        if (p.shape != (self.cells, len(self.catalog.count_nodes))
                or np.any(~np.isfinite(p)) or np.any((p < 0) | (p > 1))):
            raise ValueError("p_nodes must have shape (cells,count_states) and lie in [0,1]")
        if on_node is not None and on_cell is not None:
            raise ValueError("supply on_node or its on_cell alias, not both")
        callback = on_node if on_node is not None else on_cell
        d, lap = self._operators(z, links)
        rows = []
        symbols = [] if require_stability else None
        cache = {}
        for i in range(self.cells):
            rows.append(self._local(z[i], d[i], p[i], cache))
            if require_stability:
                symbol = self._symbol(z[i], d[i], p[i], cache)
                symbols.append(symbol)
                if not symbol.stable:
                    raise SpatialAdmissibilityError(i, symbol)
            if callback is not None:
                callback(i, self.cells)
        electronic = np.asarray([row.electronic_derivatives for row in rows])
        q = np.array([row.q_delta_bar for row in rows])
        gamma = np.array([row.gamma_bar for row in rows])
        amplitude = abs(z)
        rho = amplitude**2
        s = rho+self.delta_regularizer_bar**2
        # Only the local subtraction -kappa*rho*q² belongs to these partials.
        # The positive L4 quadratic is differentiated separately, exactly.
        gz = np.array([row.gradient_delta_bar for row in rows])
        v = (2*q*electronic[:, 2]/self.gap_ratio-2*self.kappa*rho*q)/s
        gd = 1j*v*z
        adjoint_d_gd = -self._operators(gd, links)[0]
        gradient = self.h_bar*(gz+adjoint_d_gd+2*self.kappa*lap)
        link_current = self._link_current(z, gd, links)
        gradient_density = self.kappa*(np.real(np.conj(z)*lap)-rho*q*q)
        density = electronic[:, 0]+gradient_density
        energy = float(self.h_bar*np.sum(density))
        radial = (np.real(np.conj(z/amplitude)*gradient) if np.all(amplitude > 0) else None)
        current_A = (2*E_CHARGE_C/HBAR_J_S)*self.energy_scale_J*link_current
        return NodalSpatialEvaluation(
            energy, self.energy_scale_J*energy,
            np.column_stack((gradient.real, gradient.imag)), gradient/(2*self.h_bar),
            radial, link_current, current_A, current_A/self.cross_section_m2,
            amplitude, gamma, q, density, electronic,
            tuple(symbols) if symbols is not None else None, d,
            float(self.h_bar*np.sum(gradient_density)))
