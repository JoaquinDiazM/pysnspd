"""Experimental periodic 1D discretization of model 0.4, D.8 and D.36.

Coordinates are X=x/ell0, z=Delta/Delta0, ell0=sqrt(hbar*D/(2*kB*Tc)).
The density unit is N0*Delta0**2, hence kappa=pi/4 and
Gamma/Delta0=q_delta**2/(Delta0/(kB*Tc)). Populations live at link
quadrature points in state-count coordinates and stay FIXED in variations.

A link transports the right node with U=exp(+i*link_phase); physically
link_phase=-(2e/hbar)*integral(A dx). A boundary twist is added to the last
link. Midpoint quadrature and a covariant difference define ONE energy.
Both Cartesian forces and link current below are derivatives of that energy;
the current is not an independently sampled continuum-current formula.

Only the static functional and a local longitudinal principal-symbol check
are provided. No KWT time step, potential solve, reservoir or circuit is
activated. Numerical sign checks do not replace spectral/mesh convergence.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from .energy_catalog import BCS_GAP_RATIO, E_CHARGE_C, HBAR_J_S, K_B_J_K


@dataclass(frozen=True)
class LocalEnergyDensity:
    energy_bar: float
    gradient_delta_bar: complex
    gradient_derivative_bar: complex
    amplitude_bar: float
    gamma_bar: float
    q_delta_bar: float
    electronic_derivatives: tuple[float, float, float]


@dataclass(frozen=True)
class PrincipalSymbol:
    """D.36 at fixed local z and p, not a nodal energy Hessian.

    ``uncertainty`` estimates only differentiation/roundoff error of this
    catalogue's local symbol, in eigenvalue units. It is not a bound on the
    physical regulator or count-quadrature error. Those require separate
    resolution comparisons. A positive sign is accepted only with a margin
    above this reported numerical uncertainty.
    """
    matrix: np.ndarray
    eigenvalues: np.ndarray
    uncertainty: float
    stable: bool
    amplitude_bar: float
    gamma_bar: float
    fd_steps: tuple[float, ...]
    gamma_curvatures: tuple[float, ...]


@dataclass(frozen=True)
class SpatialEvaluation:
    energy_bar: float
    energy_J: float
    # Full real gradient of the normalized, spatially integrated energy.
    gradient_cartesian_bar: np.ndarray
    # Density force dU/dDelta*: half the real gradient, divided by h_bar.
    wirtinger_force_bar: np.ndarray
    # Fixed-phase radial gradient of the integrated energy. Undefined at zero.
    amplitude_gradient_bar: np.ndarray | None
    link_current_bar: np.ndarray
    current_A: np.ndarray
    current_density_A_m2: np.ndarray
    amplitude_links_bar: np.ndarray
    gamma_links_bar: np.ndarray
    q_delta_links_bar: np.ndarray
    energy_density_bar: np.ndarray
    electronic_derivatives: np.ndarray
    principal_symbols: tuple[PrincipalSymbol, ...] | None


class SpatialAdmissibilityError(ValueError):
    """At least one D.36 sign is negative or not numerically resolved."""

    def __init__(self, link: int, symbol: PrincipalSymbol):
        self.link = link
        self.symbol = symbol
        super().__init__(
            f"D.36 rejected link {link}: lambda_min={symbol.eigenvalues[0]:.8g}, "
            f"numerical uncertainty={symbol.uncertainty:.8g}")


class PeriodicSpatialFunctional:
    """A uniform 1D mesh with periodic amplitude and an optional phase twist.

    ``p_links`` must have shape (cells, count_states); no implicit remapping or
    thermalization occurs. The catalogue's own field support is checked at
    midpoint quadrature points. Its exact normal point is separate from any
    unsupported small positive amplitude; no vortex-core extrapolation occurs.

    ``evaluate(..., require_stability=False)`` omits the D.36 calculation and
    makes NO spatial-admission claim. This is useful for differentiation tests
    and explicit negative controls. It never disables catalogue-domain checks.
    """

    delta_regularizer_bar = .1
    kappa = np.pi/4

    def __init__(self, catalog, length_m: float, cross_section_m2: float,
                 cells: int, Tc_K: float | None = None):
        if type(cells) is not int or cells < 3:
            raise ValueError("cells must be an integer at least three")
        for value, name in ((length_m, "length_m"),
                            (cross_section_m2, "cross_section_m2")):
            if not np.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive")
        vacuum = catalog.vacuum
        if Tc_K is None:
            Tc_K = vacuum.metadata.get("Tc_K", vacuum.delta0_J/(BCS_GAP_RATIO*K_B_J_K))
        if not np.isfinite(Tc_K) or Tc_K <= 0:
            raise ValueError("Tc_K must be finite and positive")
        ratio = vacuum.delta0_J/(K_B_J_K*Tc_K)
        if not np.isclose(ratio, BCS_GAP_RATIO, rtol=1e-10, atol=0.):
            raise ValueError("Tc_K and catalogue Delta0 violate the stated BCS scale convention")
        self.catalog = catalog
        self.length_m = float(length_m)
        self.cross_section_m2 = float(cross_section_m2)
        self.cells = cells
        self.Tc_K = float(Tc_K)
        self.gap_ratio = float(ratio)
        self.ell0_m = float(np.sqrt(HBAR_J_S*vacuum.D_m2_s/(2*K_B_J_K*Tc_K)))
        self.length_bar = self.length_m/self.ell0_m
        self.h_bar = self.length_bar/cells
        self.density_scale_J_m3 = vacuum.N0_per_J_m3*vacuum.delta0_J**2
        self.energy_scale_J = self.density_scale_J_m3*self.cross_section_m2*self.ell0_m

    def _population(self, occupation):
        p = np.asarray(occupation, dtype=float)
        if (p.shape != self.catalog.count_nodes.shape or np.any(~np.isfinite(p))
                or np.any((p < 0) | (p > 1))):
            raise ValueError("occupation must match count nodes and lie in [0,1]")
        return p

    def _potential(self, amplitude, gamma, p, cache):
        # Cache is local to one evaluation, has exact float keys, and stores
        # kernels independent of p. No parameter rounding or persistent stale
        # population/cross-catalogue cache is permitted.
        key = (float(amplitude), float(gamma))
        if key not in cache:
            vacuum = np.asarray(self.catalog.vacuum.evaluate(*key), float)
            cache[key] = [vacuum, None]
        vacuum, kernels = cache[key]
        if not np.any(p):
            return tuple(float(v) for v in vacuum)
        if kernels is None:
            kernels = tuple(np.asarray(k) for k in self.catalog.energy_kernel(*key))
            cache[key][1] = kernels
        weighted = self.catalog.count_weights*p
        return tuple(float(v+4*np.dot(k, weighted)) for v, k in zip(vacuum, kernels))

    @staticmethod
    def _complex(value, name):
        value = complex(value)
        if not np.isfinite(value):
            raise ValueError(f"{name} must be finite")
        return value

    def local_density(self, delta_mid_bar: complex, gradient_bar: complex,
                      occupation) -> LocalEnergyDensity:
        """Local D.8 energy and two Cartesian partial gradients.

        A complex gradient g represents (d/dRe, d/dIm), so the differential
        is Re(conj(g)*dz). ``gradient_derivative_bar`` is partial e/partial
        (dDelta/dX), holding Delta and p fixed.
        """
        z = self._complex(delta_mid_bar, "delta_mid_bar")
        d = self._complex(gradient_bar, "gradient_bar")
        return self._local(z, d, self._population(occupation), {})

    def _local(self, z, d, p, cache):
        amplitude = abs(z)
        rho = amplitude**2
        s = rho+self.delta_regularizer_bar**2
        q = float(np.imag(np.conj(z)*d)/s)
        gamma = q*q/self.gap_ratio
        u, fa, fg = self._potential(amplitude, gamma, p, cache)
        # B = partial e/partial q before differentiating q=P/(rho+delta^2).
        b = 2*q*fg/self.gap_ratio-2*self.kappa*rho*q
        v = b/s
        radial = fa*z/amplitude if amplitude else 0j
        gz = radial-2*self.kappa*q*q*z-1j*v*d-2*v*q*z
        gd = 2*self.kappa*d+1j*v*z
        value = u+self.kappa*(abs(d)**2-rho*q*q)
        return LocalEnergyDensity(float(value), complex(gz), complex(gd),
                                  float(amplitude), float(gamma), q, (u, fa, fg))

    def principal_symbol(self, delta_mid_bar: complex, gradient_bar: complex,
                         occupation, *, gamma_step: float = 2e-4) -> PrincipalSymbol:
        """Longitudinal D.36 using an analytic chain rule and refined FD in Gamma.

        Only u_GammaGamma needs a numerical derivative. Three nested stencils
        differentiate u_Gamma, with fixed p and one-sided stencils when needed
        by the catalogue's Gamma support. At q=0 its coefficient is exactly
        zero, so no unnecessary differencing across Gamma=0 is attempted.
        The maximum change between consecutive stencils is amplified by eight
        and propagated into the eigenvalue uncertainty; roundoff is included.
        """
        z = self._complex(delta_mid_bar, "delta_mid_bar")
        d = self._complex(gradient_bar, "gradient_bar")
        p = self._population(occupation)
        if not np.isfinite(gamma_step) or gamma_step <= 0:
            raise ValueError("gamma_step must be finite and positive")
        return self._symbol(z, d, p, {}, gamma_step)

    def _symbol(self, z, d, p, cache, gamma_step=2e-4):
        local = self._local(z, d, p, cache)
        amplitude, gamma, q = local.amplitude_bar, local.gamma_bar, local.q_delta_bar
        rho = amplitude**2
        s = rho+self.delta_regularizer_bar**2
        fg = local.electronic_derivatives[2]
        steps, estimates, roundoff = (), (), 0.
        if q != 0 and rho != 0:
            lo, hi = self.catalog.vacuum.gamma_axis[[0, -1]]
            step = min(float(gamma_step)*max(1., gamma), float((hi-lo)/8))
            if gamma-lo >= step and hi-gamma >= step:
                def stencil(h):
                    values = [self._potential(amplitude, gamma+j*h, p, cache)[2]
                              for j in (-1, 1)]
                    return (values[1]-values[0])/(2*h), sum(abs(v) for v in values)/h
            else:
                sign = 1. if hi-gamma >= gamma-lo else -1.
                step = min(step, float((hi-gamma if sign > 0 else gamma-lo)/3))
                def stencil(h):
                    values = [self._potential(amplitude, gamma+sign*j*h, p, cache)[2]
                              for j in (0, 1, 2)]
                    return sign*(-3*values[0]+4*values[1]-values[2])/(2*h), sum(
                        abs(c*v) for c, v in zip((3, 4, 1), values))/h
            if step <= 0 or gamma+step == gamma:
                raise ValueError("Gamma support cannot resolve a principal-symbol stencil")
            steps = (step, step/2, step/4)
            computed = [stencil(h) for h in steps]
            estimates = tuple(float(item[0]) for item in computed)
            roundoff = 32*np.finfo(float).eps*max(item[1] for item in computed)
            fgg = estimates[-1]
            fgg_error = 8*max(abs(estimates[1]-estimates[0]),
                            abs(estimates[2]-estimates[1]))+roundoff
        else:
            fgg = fgg_error = 0.
        eqq = 2*fg/self.gap_ratio+4*q*q*fgg/self.gap_ratio**2-2*self.kappa*rho
        jz = np.array([-z.imag, z.real])
        matrix = 2*self.kappa*np.eye(2)+(eqq/s**2)*np.outer(jz, jz)
        eigenvalues = np.linalg.eigvalsh(matrix)
        uncertainty = (rho/s**2)*(4*q*q/self.gap_ratio**2)*fgg_error
        uncertainty += 64*np.finfo(float).eps*max(1., np.linalg.norm(matrix, 2))
        finite = np.all(np.isfinite(matrix)) and np.isfinite(uncertainty)
        return PrincipalSymbol(matrix, eigenvalues, float(uncertainty),
                               bool(finite and eigenvalues[0] > uncertainty),
                               amplitude, gamma, steps, estimates)

    def evaluate(self, delta_bar, p_links, *, twist: float = 0., link_phases=None,
                 require_stability: bool = True,
                 on_cell: Callable[[int, int], None] | None = None) -> SpatialEvaluation:
        """Evaluate static energy, forces and current, optionally checking D.36.

        Twist convention: Delta_N=exp(i*twist)*Delta_0. A helical field can be
        entered as Delta_i=a*exp(i*q*X_i), twist=q*L_bar and zero link phases.
        Equivalently use a constant field, link_phases=twist/cells and twist=0;
        do not supply the twist twice. Current is positive for the conjugate
        positive link phase. ``on_cell(index,total)`` is called after each link
        (zero-based index), including its optional D.36 check.
        """
        z = np.asarray(delta_bar, dtype=complex)
        p = np.asarray(p_links, dtype=float)
        if z.shape != (self.cells,) or np.any(~np.isfinite(z)):
            raise ValueError("delta_bar must contain one finite complex value per node")
        if (p.shape != (self.cells, len(self.catalog.count_nodes))
                or np.any(~np.isfinite(p)) or np.any((p < 0) | (p > 1))):
            raise ValueError("p_links must have shape (cells,count_states) and lie in [0,1]")
        if not np.isfinite(twist):
            raise ValueError("twist must be finite")
        phases = (np.zeros(self.cells) if link_phases is None else
                  np.array(link_phases, dtype=float, copy=True))
        if phases.shape != (self.cells,) or np.any(~np.isfinite(phases)):
            raise ValueError("link_phases must have one finite real value per link")
        phases[-1] += twist
        transport = np.exp(1j*phases)
        right = transport*np.roll(z, -1)
        midpoint = (z+right)/2
        derivative = (right-z)/self.h_bar
        densities = []
        symbols = [] if require_stability else None
        cache = {}
        for i in range(self.cells):
            densities.append(self._local(midpoint[i], derivative[i], p[i], cache))
            if require_stability:
                symbol = self._symbol(midpoint[i], derivative[i], p[i], cache)
                symbols.append(symbol)
                if not symbol.stable:
                    raise SpatialAdmissibilityError(i, symbol)
            if on_cell is not None:
                on_cell(i, self.cells)
        gz = np.array([row.gradient_delta_bar for row in densities])
        gd = np.array([row.gradient_derivative_bar for row in densities])
        left_gradient = self.h_bar*gz/2-gd
        right_gradient = self.h_bar*gz/2+gd
        gradient = left_gradient+np.roll(np.conj(transport)*right_gradient, 1)
        link_current = np.imag(np.conj(right)*right_gradient)
        density = np.array([row.energy_bar for row in densities])
        energy = float(self.h_bar*np.sum(density))
        radial = (np.real(np.conj(z/abs(z))*gradient) if np.all(abs(z) > 0) else None)
        current_A = (2*E_CHARGE_C/HBAR_J_S)*self.energy_scale_J*link_current
        return SpatialEvaluation(
            energy, self.energy_scale_J*energy,
            np.column_stack((gradient.real, gradient.imag)), gradient/(2*self.h_bar),
            radial, link_current, current_A, current_A/self.cross_section_m2,
            np.array([row.amplitude_bar for row in densities]),
            np.array([row.gamma_bar for row in densities]),
            np.array([row.q_delta_bar for row in densities]), density,
            np.array([row.electronic_derivatives for row in densities]),
            tuple(symbols) if symbols is not None else None)
