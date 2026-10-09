"""Internal-energy diagnostic of the SAME finite thermal Usadel graph action.

U_N = F_N - t*dF_N/dt, t=T/Tc. Gap, links, graph, material scales and the
NUMBER of Matsubara modes are held fixed, not the physical energy cutoff.
The action is relative to its normal state and contains no normal-metal
background energy, ultraviolet tail or fitted heat. This module is not a
finite-photon energy closure and does not evolve distributions or fields.

Prescribed spectral boundaries can move with temperature, even at fixed gap.
Their temperature derivative must be supplied explicitly, and their reaction
work is retained. The free-node derivative uses the same spectral Jacobian.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Real
import math
import numpy as np
from scipy.sparse.linalg import spsolve

from . import thermal_spatial_usadel as thermal


@dataclass(frozen=True)
class ThermalInternalEnergy:
    free_energy_bar: float
    temperature_derivative_bar: float
    internal_energy_bar: float
    explicit_temperature_derivative_bar: float
    boundary_reaction_temperature_derivative_bar: float
    free_stationarity_temperature_derivative_bar: float
    temperature_tangents: tuple[np.ndarray, ...]
    maximum_free_tangent_residual: float
    maximum_spectral_residual: float
    temperature_ratio: float
    modes: int
    maximum_matsubara_energy_bar: float
    energy_scale_J: float
    boundary_contract: str
    reference: str = 'same finite thermal action minus its normal state; no normal background'
    cutoff_contract: str = 'fixed mode count; physical cutoff moves with temperature'
    scope: str = 'independent finite-action energy diagnostic; no tail or dynamic heat admission'

    @property
    def free_energy_J(self):
        return self.energy_scale_J*self.free_energy_bar

    @property
    def internal_energy_J(self):
        return self.energy_scale_J*self.internal_energy_bar


def _temperature_tangent(graph, d, solution, alpha, temperature_ratio, fixed_derivative):
    """At fixed gap: J du/dt = -m*u*(epsilon/t), with explicit fixed rows."""
    value, jacobian = thermal.spectral_residual_jacobian(
        graph, d, solution.epsilon, solution.u, alpha)
    fixed = thermal._fixed_indices(solution.fixed_nodes, graph.n_nodes)
    free = np.ones(graph.n_nodes, bool)
    free[fixed] = False
    free_components = np.flatnonzero(np.repeat(free, 2))
    fixed_components = np.flatnonzero(np.repeat(~free, 2))
    boundary = np.asarray(fixed_derivative, complex)
    if boundary.shape != fixed.shape or np.any(~np.isfinite(boundary)):
        raise ValueError('one finite complex du/d(T/Tc) per prescribed spectral boundary required')
    # The fixed values use exactly solution.fixed_nodes order.
    full_boundary = np.zeros(graph.n_nodes, complex)
    full_boundary[fixed] = boundary
    tangent_real = np.column_stack((full_boundary.real, full_boundary.imag)).ravel()
    epsilon_rate = solution.epsilon/temperature_ratio
    source = -graph.area_weights*solution.u*epsilon_rate
    source_real = np.column_stack((source.real, source.imag)).ravel()
    if len(free_components):
        rhs = source_real[free_components]-jacobian[free_components][:, fixed_components]@tangent_real[fixed_components]
        tangent_real[free_components] = spsolve(jacobian[free_components][:, free_components], rhs)
    if np.any(~np.isfinite(tangent_real)):
        raise FloatingPointError('nonfinite thermal temperature tangent; no regularization added')
    tangent = tangent_real[0::2]+1j*tangent_real[1::2]
    residual = jacobian@tangent_real-source_real
    maximum = float(np.max(abs(residual[free_components]), initial=0.))
    reaction = float(np.real(np.vdot(value.gradient_u[fixed], tangent[fixed])))
    stationarity = float(np.real(np.vdot(value.gradient_u[free], tangent[free])))
    tangent.setflags(write=False)
    return value, tangent, reaction, stationarity, maximum


def evaluate_internal_energy(graph, d, temperature_ratio, solutions, alpha=None, *,
                             energy_scale_J, boundary_contract,
                             boundary_temperature_derivatives=None):
    """Observe F_N and U_N off the gap equation, on stationary spectral roots.

    boundary_temperature_derivatives[n] has one entry per
    solutions[n].fixed_nodes, in that exact order. It is du_b/d(T/Tc),
    at fixed prescribed gap and links. Examples are zero for externally fixed
    spectral values, or -u_b/t for the BCS law u_b=d_b/epsilon_n(t).
    There is no inferred contact law. boundary_contract documents the law.
    A graph with no fixed spectral nodes may omit this derivative sequence.

    energy_scale_J is the caller's fixed N0*(kB*Tc)^2*film_thickness*ell0^2.
    The returned derivative is with respect to t=T/Tc, not kelvin or time.
    U_N is relative to the same finite-action normal state. At fixed bath
    temperature the missing normal background is gap independent, but this
    fact does not certify the cutoff/tail or charge-sector energy.
    """
    scale_error = 'an explicit finite positive real physical energy scale in joules is required'
    if isinstance(energy_scale_J, (bool, np.bool_)) or not isinstance(energy_scale_J, Real):
        raise ValueError(scale_error)
    try:
        physical_scale = float(energy_scale_J)
    except (OverflowError, ValueError) as error:
        raise ValueError(scale_error) from error
    if not math.isfinite(physical_scale) or physical_scale <= 0:
        raise ValueError(scale_error)
    if not isinstance(boundary_contract, str) or not boundary_contract.strip():
        raise ValueError('an explicit spectral boundary temperature contract is required')
    solutions = tuple(solutions)
    # This checks source signatures, consecutive frequencies and stationarity receipts.
    reference = thermal.evaluate_thermal(graph, d, temperature_ratio, solutions, alpha)
    gap = thermal._complex_vector(d, graph.n_nodes, 'd')
    links = thermal._alpha(graph, alpha)
    t = float(temperature_ratio)
    fixed_sets = tuple(thermal._fixed_indices(s.fixed_nodes, graph.n_nodes) for s in solutions)
    if boundary_temperature_derivatives is None:
        if any(len(fixed) for fixed in fixed_sets):
            raise ValueError('prescribed spectral nodes require explicit temperature derivatives; zero is not inferred')
        boundaries = tuple(np.empty(0, complex) for _ in solutions)
    else:
        boundaries = tuple(boundary_temperature_derivatives)
        if len(boundaries) != len(solutions):
            raise ValueError('one boundary temperature derivative vector per Matsubara mode required')
    mass_gap_squared = float(np.dot(graph.area_weights, abs(gap)**2))
    weight = 2*np.pi*t
    # The explicit counterterm weight/epsilon is temperature independent
    # at fixed mode count; it is retained in F but cancels from its derivative.
    explicit_derivative = mass_gap_squared/t
    boundary_derivative = 0.
    free_derivative = 0.
    tangents = []
    tangent_residual = 0.
    for solution, boundary in zip(solutions, boundaries):
        value, tangent, reaction, stationarity, residual = _temperature_tangent(
            graph, gap, solution, links, t, boundary)
        # Envelope derivative with fixed spectral fields, followed by the
        # explicit free-stationarity and prescribed-boundary reactions.
        partial_epsilon = float(2*np.dot(graph.area_weights, abs(value.f)**2/(1+value.g)))
        explicit_derivative += weight/t*(value.energy+solution.epsilon*partial_epsilon)
        boundary_derivative += weight*reaction
        free_derivative += weight*stationarity
        tangents.append(tangent)
        tangent_residual = max(tangent_residual, residual)
    derivative = explicit_derivative+boundary_derivative+free_derivative
    internal = reference.energy-t*derivative
    if not np.all(np.isfinite((derivative, internal))):
        raise FloatingPointError('nonfinite finite-action internal energy')
    return ThermalInternalEnergy(reference.energy, float(derivative), float(internal),
        float(explicit_derivative), float(boundary_derivative), float(free_derivative),
        tuple(tangents), tangent_residual, reference.maximum_spectral_residual,
        t, reference.modes, float(solutions[-1].epsilon), physical_scale,
        boundary_contract)
