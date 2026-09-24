"""Homogeneous dirty-superconductor DC reservoirs for an interior strip.

The conventions are those of :mod:`thermal_spatial_usadel`: d=Delta/(kB Tc),
X=x/ell0, ell0=sqrt(hbar D/(2 kB Tc)), and q=d(theta)/dX at zero magnetic field.
For each positive Matsubara frequency, the homogeneous Usadel equation is

    epsilon*u + q**2*u/sqrt(1+u**2) = d,  f=u/sqrt(1+u**2).

It is the theta variation of the SAME continuum energy density
2*epsilon*(1-g)-2*d*f+q**2*f**2. The gap uses the same finite, renormalized
Matsubara sum as the graph solver. No fitted depairing factor or volume source
is added. Reservoir spectra therefore include current depairing: u != d/epsilon.

The ascending I(q) branch is the homogeneous current-biased branch. This is
not a certificate against vortex entry, transverse instability, or heating.
An irregular finite graph approximates this continuum state; the actual graph
must be relaxed with these reservoirs and its stationarity/uniformity measured.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.constants import Boltzmann, elementary_charge, hbar
from scipy.optimize import brentq, minimize_scalar


def _parameters(temperature_ratio, matsubara_count):
    t = float(temperature_ratio)
    if not np.isfinite(t) or not 0 < t < 1:
        raise ValueError('Require 0 < temperature_ratio < 1')
    if type(matsubara_count) is not int or matsubara_count < 1:
        raise ValueError('Require a positive integer Matsubara count')
    epsilon = 2*np.pi*t*(np.arange(matsubara_count)+.5)
    # Otherwise the finite-cutoff renormalized gap energy is unbounded below.
    coefficient = np.log(t)+2*np.pi*t*np.sum(1/epsilon)
    if coefficient <= 0:
        raise ValueError('Finite Matsubara cutoff is too small at this temperature')
    return t, epsilon


def homogeneous_spectrum(gap_bar, epsilon_bar, depairing_bar):
    """Return positive-Matsubara (u,f,g) for a specified gap and depairing.

    ``depairing_bar=q_bar**2`` in the continuum. Accepting it explicitly also
    permits an independent check of the exact Cartesian graph symbol
    2*(1-cos(q*h))/h**2. It does not change a graph operator in production.
    The scalar spectral equation is strictly increasing in u; safeguarded
    Newton uses its physical bracket [0,d/epsilon].
    """
    d, a = float(gap_bar), float(depairing_bar)
    ep = np.asarray(epsilon_bar, dtype=float)
    if (not np.isfinite(d) or d < 0 or not np.isfinite(a) or a < 0
            or ep.ndim != 1 or not len(ep) or np.any(~np.isfinite(ep))
            or np.any(ep <= 0)):
        raise ValueError('Nonnegative gap/depairing and positive finite frequencies required')
    if a == 0 or d == 0:
        u = d/ep
    else:
        lower, upper = np.zeros_like(ep), d/ep
        u = d/(ep+a)
        for _ in range(100):
            g = 1/np.hypot(1., u)
            residual = ep*u+a*u*g-d
            if np.max(abs(residual)) <= 8*np.finfo(float).eps*max(1., d):
                break
            lower = np.where(residual < 0, u, lower)
            upper = np.where(residual >= 0, u, upper)
            newton = u-residual/(ep+a*g**3)
            u = np.where((newton > lower) & (newton < upper),
                         newton, .5*(lower+upper))
        else:
            raise RuntimeError('Monotone homogeneous spectral solve did not converge')
    g = 1/np.hypot(1., u)
    return u, u*g, g


@dataclass(frozen=True)
class BulkCurrentReference:
    temperature_ratio: float
    matsubara_count: int
    q_bar: float
    gap_bar: float
    epsilon_bar: np.ndarray
    u: np.ndarray
    f: np.ndarray
    g: np.ndarray
    current_density_bar: float
    gap_equation_residual: float

    def total_current_A(self, *, width_m, Tc_K, diffusion_m2_s,
                        sheet_resistance_ohm):
        """Integrated strip current; graph link unit I0=kB*Tc/(2*e*Rsheet).

        The continuum dimensionless sheet-current density is
        jbar=4*pi*t*q*sum(f_n**2). Thus I=I0*(width/ell0)*jbar.
        Film thickness is already represented by the sheet resistance.
        """
        scale = physical_scales(Tc_K=Tc_K, diffusion_m2_s=diffusion_m2_s,
                                sheet_resistance_ohm=sheet_resistance_ohm)
        if not np.isfinite(width_m) or width_m <= 0:
            raise ValueError('Positive finite strip width required')
        return float(scale['current_unit_A']*width_m/scale['ell0_m']
                     *self.current_density_bar)


def physical_scales(*, Tc_K, diffusion_m2_s, sheet_resistance_ohm):
    values = (Tc_K, diffusion_m2_s, sheet_resistance_ohm)
    if any(not np.isfinite(value) or value <= 0 for value in values):
        raise ValueError('Positive finite material parameters required')
    E0 = Boltzmann*Tc_K
    return dict(energy_unit_J=E0, ell0_m=float(np.sqrt(hbar*diffusion_m2_s/(2*E0))),
                current_unit_A=float(E0/(2*elementary_charge*sheet_resistance_ohm)))


def current_density_derivative(reference):
    """Analytic d(jbar)/d(qbar), including the self-consistent gap change.

    Differentiate the spectral equation and the undivided finite-sum gap
    equation together. A normal state has zero supercurrent derivative.
    """
    if reference.gap_bar == 0:
        return 0.
    t, q, ep, f, g = (reference.temperature_ratio, reference.q_bar,
                     reference.epsilon_bar, reference.f, reference.g)
    coefficient = np.log(t)+2*np.pi*t*np.sum(1/ep)
    response = g**3/(ep+q*q*g**3)
    radial_stiffness = coefficient-2*np.pi*t*np.sum(response)
    if radial_stiffness <= 0:
        raise ValueError('Positive homogeneous radial stiffness required')
    gap_prime = -4*np.pi*t*q*np.sum(response*f)/radial_stiffness
    f_prime = response*(gap_prime-2*q*f)
    return float(4*np.pi*t*(np.sum(f*f)+2*q*np.sum(f*f_prime)))


def differential_inductance_per_length_H_m(reference, *, width_m, Tc_K,
                                           diffusion_m2_s, sheet_resistance_ohm):
    """Bulk differential kinetic inductance at this fixed DC operating point.

    L'= [hbar/(2e)]/[ell0*dI/dqbar]. It includes the current-induced gap
    suppression but is a CONSTANT evaluated once, not a dynamic circuit state.
    Only the ascending homogeneous branch admits a positive value. This does
    not supply the inductance of a different-width series wire; evaluate that
    wire at its own width/current or use its measured fixed inductance.
    """
    scale = physical_scales(Tc_K=Tc_K, diffusion_m2_s=diffusion_m2_s,
                            sheet_resistance_ohm=sheet_resistance_ohm)
    if not np.isfinite(width_m) or width_m <= 0:
        raise ValueError('Positive finite strip width required')
    derivative = current_density_derivative(reference)
    if derivative <= 0:
        raise ValueError('Differential inductance requires the ascending current branch')
    return float(hbar/(2*elementary_charge)
                 /(scale['current_unit_A']*width_m*derivative))


def solve_bulk_reference(temperature_ratio, matsubara_count, q_bar):
    """Self-consistent continuum bulk gap and spectrum at imposed q.

    Returns the normal solution when depairing exceeds the finite-sum
    superconducting endpoint. The endpoint is not the smaller current maximum.
    Arrays are read-only so a saved reservoir cannot be modified by a warm start.
    """
    t, ep = _parameters(temperature_ratio, matsubara_count)
    q = float(q_bar)
    if not np.isfinite(q):
        raise ValueError('Finite gauge-invariant phase gradient required')
    a = q*q
    if not np.isfinite(a):
        raise ValueError('Finite squared phase gradient required')

    def gap_equation(d):
        ratio = 1/(ep+a) if d == 0 else homogeneous_spectrum(d, ep, a)[1]/d
        return float(np.log(t)+2*np.pi*t*np.sum(1/ep-ratio))

    normal_coefficient = gap_equation(0.)
    if normal_coefficient >= 0:
        gap = 0.
    else:
        upper = 4.
        while gap_equation(upper) <= 0:
            upper *= 2
            if upper > 1e6:
                raise RuntimeError('No finite-sum superconducting gap bracket found')
        gap = brentq(gap_equation, 0., upper, xtol=2e-13, rtol=2e-14)
    u, f, g = homogeneous_spectrum(gap, ep, a)
    current = float(4*np.pi*t*q*np.sum(f*f))
    # The undivided equation is also zero for a normal solution.
    residual = float(gap*gap_equation(gap))
    for array in (ep, u, f, g):
        array.setflags(write=False)
    return BulkCurrentReference(t, matsubara_count, q, float(gap), ep, u, f, g,
                                current, residual)


def normal_endpoint_q(temperature_ratio, matsubara_count):
    """q where the superconducting gap merges with the normal branch."""
    t, ep = _parameters(temperature_ratio, matsubara_count)

    def equation(q):
        return float(np.log(t)+2*np.pi*t*np.sum(1/ep-1/(ep+q*q)))

    upper = 1.
    while equation(upper) < 0:
        upper *= 2
    return float(brentq(equation, 0., upper, xtol=2e-13))


def find_depairing_current(temperature_ratio, matsubara_count):
    """Return the positive maximum of the homogeneous finite-sum I(q)."""
    q_end = normal_endpoint_q(temperature_ratio, matsubara_count)
    result = minimize_scalar(
        lambda q: -solve_bulk_reference(temperature_ratio, matsubara_count, q).current_density_bar,
        bounds=(0., q_end), method='bounded', options={'xatol': 2e-10})
    if not result.success:
        raise RuntimeError('Homogeneous depairing-current search failed')
    return solve_bulk_reference(temperature_ratio, matsubara_count, float(result.x))


def solve_bulk_current(temperature_ratio, matsubara_count, target_current_A, *,
                       width_m, Tc_K, diffusion_m2_s, sheet_resistance_ohm):
    """Select the ascending homogeneous branch for a prescribed DC current.

    Current reversal reverses q while preserving the gap. Targets beyond the
    homogeneous depairing maximum are rejected, never clipped to that maximum.
    """
    target = float(target_current_A)
    if not np.isfinite(target):
        raise ValueError('Finite target current required')
    material = dict(width_m=width_m, Tc_K=Tc_K, diffusion_m2_s=diffusion_m2_s,
                    sheet_resistance_ohm=sheet_resistance_ohm)
    critical = find_depairing_current(temperature_ratio, matsubara_count)
    critical_A = critical.total_current_A(**material)
    if abs(target) > critical_A:
        raise ValueError(f'Target |I|={abs(target):.9g} A exceeds homogeneous '
                         f'depairing current {critical_A:.9g} A')
    if target == 0:
        return solve_bulk_reference(temperature_ratio, matsubara_count, 0.)
    q = brentq(lambda value: solve_bulk_reference(
        temperature_ratio, matsubara_count, value).total_current_A(**material)-abs(target),
        0., critical.q_bar, xtol=2e-13, rtol=2e-14)
    return solve_bulk_reference(temperature_ratio, matsubara_count, np.sign(target)*q)


def intrinsic_boundary_fields(graph, reference, *, nodes=None, axis=0,
                              phase_offset=0., coordinate_origin_bar=0.):
    """Bulk fields sampled on this section, including matched reservoir spectra.

    With ``nodes=None`` the result also supplies a complete initial state.
    For end boundaries pass their indices; lateral boundaries remain the natural
    no-flux boundaries of the graph. No BCS zero-current contact is introduced.
    At zero magnetic field alpha=0. In a co-moving gauge with real gap use
    alpha_ij=-q*(X_j-X_i), and remove the common phase from d,u,f.
    """
    coordinates = graph.coordinates_bar
    if coordinates is None or axis not in (0, 1):
        raise ValueError('Planar graph coordinates and axis 0 or 1 required')
    if not np.isfinite(phase_offset) or not np.isfinite(coordinate_origin_bar):
        raise ValueError('Finite phase origin required')
    if nodes is None:
        selected = np.arange(graph.n_nodes)
    else:
        selected = np.asarray(nodes)
        if (selected.dtype.kind not in 'iu' or selected.ndim != 1
                or np.any(selected < 0) or np.any(selected >= graph.n_nodes)
                or len(np.unique(selected)) != len(selected)):
            raise ValueError('Unique valid integer node indices required')
    theta = phase_offset+reference.q_bar*(coordinates[selected, axis]-coordinate_origin_bar)
    phase = np.exp(1j*theta)
    return dict(delta_bar=reference.gap_bar*phase,
                u=reference.u[:, None]*phase[None, :],
                f=reference.f[:, None]*phase[None, :],
                g=np.broadcast_to(reference.g[:, None], (len(reference.g), len(selected))).copy())
