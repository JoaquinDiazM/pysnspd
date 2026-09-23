"""Thermal superconducting reservoir branch from D.26--D.27 and CM.6.

This opt-in module uses the caller's electronic catalogue and its actual count
quadrature. theta=kB*Tb/Delta0 is explicit; theta=0 is the vacuum limit. Energy
and first derivatives share the same Fermi grand potential, including the
regularized uniform-gradient remainder. The derivative of an FD occupation
is NOT added to a fixed-occupation electronic force: it cancels against the
entropy variation when differentiating the thermal free energy.

Branch roots are confined to a supplied amplitude bracket. Local stability
requires positive F_aa and positive reduced dJ/dq, with finite-difference
uncertainty reported. Checking a finite continuation path does not establish
global uniqueness or stability between its samples. No hidden mesh search,
current clipping, or extrapolation selects a different root.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from itertools import product

import numpy as np
from scipy.optimize import brentq
from scipy.special import expit

from .energy_catalog import BCS_GAP_RATIO, E_CHARGE_C, HBAR_J_S, K_B_J_K


@dataclass(frozen=True)
class ReservoirEvaluation:
    free_energy_bar: float
    amplitude_force_bar: float
    current_bar: float
    amplitude_bar: float
    q_bare_bar: float
    gamma_bar: float
    occupation: np.ndarray
    electronic_moments_bar: tuple[float, float, float]
    current_A: float


@dataclass(frozen=True)
class ReservoirBranch(ReservoirEvaluation):
    amplitude_curvature: float
    differential_current_bar: float
    amplitude_curvature_uncertainty: float
    differential_current_uncertainty: float
    thermal_hessian: np.ndarray
    hessian_uncertainties: np.ndarray
    maxwell_disagreement: float
    derivative_steps: tuple[tuple[float, float], ...]
    amplitude_bracket: tuple[float, float]
    checked_qs: tuple[float, ...]
    stable: bool
    L_res_diff_H: float
    bath_theta: float
    resolved_length_m: float
    current_scale_A: float
    ell0_m: float


@dataclass(frozen=True)
class InductancePartition:
    total_reference_H: float
    resolved_differential_H: float
    exterior_fixed_H: float
    reference_current_A: float
    reference_q_bare_bar: float
    reference_amplitude_bar: float
    bath_theta: float
    resolved_length_m: float
    scope: str = "Uniform reference branch on the declared resolved length; exterior value fixed before dynamics."


class ReservoirBranchError(ValueError):
    """The requested bracket or locally stable superconducting branch failed."""


class SuperconductingReservoir:
    """Uniform thermal branch with an explicitly identified resolved length.

    The finite-difference uncertainty describes numerical differentiation of
    this catalogue only. It does not include count-grid, regulator, material,
    geometry, or nonuniform-field uncertainty. The supplied resolved length
    is the length between the energy-domain ports, not an inferred hotspot
    length. No circuit is integrated by this class.
    """

    kappa = np.pi/4
    delta_regularizer_bar = .1

    def __init__(self, catalog, bath_theta: float, *, cross_section_m2: float,
                 resolved_length_m: float, Tc_K: float | None = None,
                 derivative_step: float = 2e-4):
        if not np.isfinite(bath_theta) or bath_theta < 0:
            raise ValueError("bath_theta must be finite and nonnegative")
        for value, name in ((cross_section_m2, "cross_section_m2"),
                            (resolved_length_m, "resolved_length_m"),
                            (derivative_step, "derivative_step")):
            if not np.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive")
        vacuum = catalog.vacuum
        if Tc_K is None:
            Tc_K = vacuum.metadata.get("Tc_K", vacuum.delta0_J/(BCS_GAP_RATIO*K_B_J_K))
        if not np.isfinite(Tc_K) or Tc_K <= 0:
            raise ValueError("Tc_K must be finite and positive")
        ratio = vacuum.delta0_J/(K_B_J_K*Tc_K)
        if not np.isclose(ratio, BCS_GAP_RATIO, rtol=1e-10, atol=0.):
            raise ValueError("Tc_K and catalogue Delta0 violate the stated BCS convention")
        self.catalog = catalog
        self.bath_theta = float(bath_theta)
        self.cross_section_m2 = float(cross_section_m2)
        self.resolved_length_m = float(resolved_length_m)
        self.Tc_K = float(Tc_K)
        self.gap_ratio = float(ratio)
        self.derivative_step = float(derivative_step)
        self.ell0_m = float(np.sqrt(HBAR_J_S*vacuum.D_m2_s/(2*K_B_J_K*Tc_K)))
        self.length_bar = resolved_length_m/self.ell0_m
        self.density_scale_J_m3 = vacuum.N0_per_J_m3*vacuum.delta0_J**2
        self.current_scale_A = (2*E_CHARGE_C/HBAR_J_S)*self.density_scale_J_m3*cross_section_m2*self.ell0_m

    def evaluate(self, amplitude_bar: float, q_bare_bar: float) -> ReservoirEvaluation:
        """Return F, F_a and F_q at fixed bath temperature in normalized units."""
        a, q = float(amplitude_bar), float(q_bare_bar)
        if not np.isfinite(a) or a < 0 or not np.isfinite(q):
            raise ValueError("amplitude must be nonnegative and both fields finite")
        rho, d2 = a*a, self.delta_regularizer_bar**2
        s = rho+d2
        m = rho/s
        ma = 2*a*d2/s**2
        q_delta = m*q
        gamma = q_delta*q_delta/self.gap_ratio
        vacuum = self.catalog.vacuum.evaluate(a, gamma)  # Domain checked; never clip.
        if self.bath_theta == 0:
            p = np.zeros_like(self.catalog.count_nodes)
            free_electronic, ua, ug = vacuum
            electronic_energy = vacuum[0]
        else:
            energy, ea, eg = self.catalog.energy_kernel(a, gamma)
            theta = self.bath_theta
            p = expit(-energy/theta)
            weights = self.catalog.count_weights
            free_electronic = vacuum[0]-4*theta*np.dot(weights, np.logaddexp(0., -energy/theta))
            electronic_energy = vacuum[0]+4*np.dot(weights*energy, p)
            ua = vacuum[1]+4*np.dot(weights*ea, p)
            ug = vacuum[2]+4*np.dot(weights*eg, p)
        remainder = self.kappa*rho*q*q*(1-m*m)
        gamma_a = 2*q_delta*ma*q/self.gap_ratio
        gamma_q = 2*q_delta*m/self.gap_ratio
        remainder_a = self.kappa*q*q*(2*a*(1-m*m)-2*rho*m*ma)
        remainder_q = 2*self.kappa*rho*q*(1-m*m)
        fa = ua+ug*gamma_a+remainder_a
        fq = ug*gamma_q+remainder_q
        p = np.asarray(p, float)
        p.setflags(write=False)
        return ReservoirEvaluation(float(free_electronic+remainder), float(fa), float(fq),
            a, q, float(gamma), p, (float(electronic_energy), float(ua), float(ug)),
            float(self.current_scale_A*fq))

    def _bracket(self, amplitude_bracket):
        bracket = np.asarray(amplitude_bracket, float)
        if (bracket.shape != (2,) or np.any(~np.isfinite(bracket))
                or bracket[0] <= 0 or bracket[0] >= bracket[1]):
            raise ValueError("amplitude_bracket needs two ordered, positive, finite bounds")
        lo, hi = self.catalog.vacuum.delta_axis[[0, -1]]
        if bracket[0] < lo or bracket[1] > hi:
            raise ValueError("amplitude_bracket lies outside catalogue support")
        return float(bracket[0]), float(bracket[1])

    def _cached(self, a, q, cache):
        key = (float(a), float(q))
        if key not in cache:
            cache[key] = self.evaluate(*key)
        return cache[key]

    def _amplitude_root(self, q, bracket, cache):
        lo, hi = bracket
        f_lo = self._cached(lo, q, cache).amplitude_force_bar
        f_hi = self._cached(hi, q, cache).amplitude_force_bar
        if f_lo*f_hi > 0:
            raise ReservoirBranchError("amplitude force does not bracket a root; no alternate-root search performed")
        root = brentq(lambda a: self._cached(a, q, cache).amplitude_force_bar,
                      lo, hi, xtol=2e-11, rtol=2e-13)
        return self._cached(root, q, cache)

    def _stability(self, point, bracket, cache):
        a, q = point.amplitude_bar, point.q_bare_bar
        a_lo, a_hi = self.catalog.vacuum.delta_axis[[0, -1]]
        ha = min(self.derivative_step*max(1., a), (a-a_lo)/3, (a_hi-a)/3)
        m = a*a/(a*a+self.delta_regularizer_bar**2)
        gamma_max = self.catalog.vacuum.gamma_axis[-1]
        q_max = np.sqrt(self.gap_ratio*gamma_max)/m
        hq = min(self.derivative_step*max(1., abs(q)), (q_max-abs(q))/3)
        if min(ha, hq) <= 0 or a+ha == a or q+hq == q:
            raise ReservoirBranchError("catalogue support cannot resolve local branch stability")
        steps, estimates, roundoff = [], [], []
        for factor in (1., .5, .25):
            da, dq = ha*factor, hq*factor
            ap, am = self._cached(a+da, q, cache), self._cached(a-da, q, cache)
            qp, qm = self._cached(a, q+dq, cache), self._cached(a, q-dq, cache)
            hessian = np.array([
                [(ap.amplitude_force_bar-am.amplitude_force_bar)/(2*da),
                 (qp.amplitude_force_bar-qm.amplitude_force_bar)/(2*dq)],
                [(ap.current_bar-am.current_bar)/(2*da), (qp.current_bar-qm.current_bar)/(2*dq)]])
            magnitudes = np.array([
                [(abs(ap.amplitude_force_bar)+abs(am.amplitude_force_bar))/da,
                 (abs(qp.amplitude_force_bar)+abs(qm.amplitude_force_bar))/dq],
                [(abs(ap.current_bar)+abs(am.current_bar))/da,
                 (abs(qp.current_bar)+abs(qm.current_bar))/dq]])
            steps.append((float(da), float(dq)))
            estimates.append(hessian)
            roundoff.append(32*np.finfo(float).eps*np.maximum(1., magnitudes))
        hessian = estimates[-1]
        error = 8*np.maximum(abs(estimates[1]-estimates[0]), abs(estimates[2]-estimates[1]))
        error += np.maximum.reduce(roundoff)
        aa, aq, qa, qq = hessian[0, 0], hessian[0, 1], hessian[1, 0], hessian[1, 1]
        if not np.all(np.isfinite(hessian)) or aa <= error[0, 0]:
            raise ReservoirBranchError("amplitude curvature is nonpositive or unresolved")
        slope = qq-qa*aq/aa  # Implicit differentiation of F_a(a(q),q)=0.
        possible = [c*b/den for den, b, c in product(
            (aa-error[0, 0], aa+error[0, 0]),
            (aq-error[0, 1], aq+error[0, 1]),
            (qa-error[1, 0], qa+error[1, 0]))]
        lower = qq-error[1, 1]-max(possible)
        upper = qq+error[1, 1]-min(possible)
        slope_error = max(abs(slope-lower), abs(upper-slope))
        if not np.isfinite(slope) or slope <= slope_error:
            raise ReservoirBranchError("reduced differential current is nonpositive or unresolved")
        inductance = (HBAR_J_S/(2*E_CHARGE_C))*self.length_bar/(self.current_scale_A*slope)
        return ReservoirBranch(**point.__dict__, amplitude_curvature=float(aa),
            differential_current_bar=float(slope), amplitude_curvature_uncertainty=float(error[0, 0]),
            differential_current_uncertainty=float(slope_error), thermal_hessian=hessian,
            hessian_uncertainties=error, maxwell_disagreement=float(abs(aq-qa)),
            derivative_steps=tuple(steps), amplitude_bracket=bracket, checked_qs=(q,), stable=True,
            L_res_diff_H=float(inductance), bath_theta=self.bath_theta,
            resolved_length_m=self.resolved_length_m, current_scale_A=self.current_scale_A,
            ell0_m=self.ell0_m)

    def solve_at_q(self, q_bare_bar: float, amplitude_bracket=(.8, 1.02), *,
                   continuation_qs=()) -> ReservoirBranch:
        """Select a root in the explicit weak-branch bracket and check stability.

        q=0 is always checked as the low-current anchor. Optional intermediate
        q values must be ordered along the path from zero to the requested q;
        they are checked exactly as supplied. No unrequested grid is searched.
        Returned checked_qs states the finite extent of this branch evidence.
        """
        q = float(q_bare_bar)
        if not np.isfinite(q):
            raise ValueError("q_bare_bar must be finite")
        bracket = self._bracket(amplitude_bracket)
        intermediate = np.asarray(continuation_qs, float)
        if intermediate.ndim != 1 or np.any(~np.isfinite(intermediate)):
            raise ValueError("continuation_qs must be a finite ordered path")
        if len(intermediate):
            if q == 0:
                raise ValueError("a zero-current request has no intermediate continuation path")
            fractions = intermediate/q
            if np.any((fractions <= 0) | (fractions >= 1)) or np.any(np.diff(fractions) <= 0):
                raise ValueError("continuation_qs must lie strictly in order between zero and q")
        path = (0., *map(float, intermediate), q) if q != 0 else (0.,)
        cache = {}
        branch = None
        for value in path:
            point = self._amplitude_root(value, bracket, cache)
            branch = self._stability(point, bracket, cache)
        return replace(branch, checked_qs=tuple(path))

    def solve_at_current(self, current_A: float, *, q_bracket,
                         amplitude_bracket=(.8, 1.02)) -> ReservoirBranch:
        """Invert current on explicitly bracketed, sampled stable weak states.

        Endpoints, every scalar-solver trial, and the resulting state are
        checked. This is not a global monotonicity proof for a broad bracket;
        the caller must register a weak branch and any additional continuation
        samples it needs. Trial checks add work but do not hide a mesh sweep.
        """
        target = float(current_A)
        bounds = np.asarray(q_bracket, float)
        if (not np.isfinite(target) or bounds.shape != (2,) or np.any(~np.isfinite(bounds))
                or bounds[0] >= bounds[1]):
            raise ValueError("finite current and an ordered finite q_bracket are required")
        bracket = self._bracket(amplitude_bracket)
        endpoints = [self.solve_at_q(q, bracket) for q in bounds]
        residuals = [state.current_A-target for state in endpoints]
        if residuals[0]*residuals[1] > 0:
            raise ReservoirBranchError("current is outside the supplied stable bracket; no clipping performed")
        cache = {}
        checked = {0., *map(float, bounds)}
        def residual(q):
            point = self._amplitude_root(q, bracket, cache)
            self._stability(point, bracket, cache)
            checked.add(float(q))
            return point.current_A-target
        root = brentq(residual, *bounds, xtol=2e-11, rtol=2e-13)
        result = self.solve_at_q(root, bracket)
        return replace(result, checked_qs=tuple(sorted(checked | set(result.checked_qs))))

    def inductance_partition(self, reference: ReservoirBranch, *,
                             total_reference_H: float = 10e-9) -> InductancePartition:
        """CM.6 at one identified reference state; never a dynamic subtraction."""
        if not np.isfinite(total_reference_H) or total_reference_H <= 0:
            raise ValueError("total_reference_H must be finite and positive")
        matching = (reference.bath_theta == self.bath_theta
                    and reference.resolved_length_m == self.resolved_length_m
                    and reference.current_scale_A == self.current_scale_A
                    and reference.ell0_m == self.ell0_m)
        if not matching or not reference.stable:
            raise ValueError("reference branch does not belong to these reservoir scales")
        resolved = float(reference.L_res_diff_H)
        exterior = float(total_reference_H-resolved)
        if not np.isfinite(resolved) or resolved <= 0 or exterior <= 0:
            raise ReservoirBranchError("CM.6 needs a positive resolved and exterior inductance; no repair performed")
        return InductancePartition(float(total_reference_H), resolved, exterior,
            reference.current_A, reference.q_bare_bar, reference.amplitude_bar,
            self.bath_theta, self.resolved_length_m)
