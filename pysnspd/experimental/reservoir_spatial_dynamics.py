"""Instantaneous fixed-radius, prescribed-current reservoir loads for KWT.

The phase load is prescribed from an independently selected reservoir current,
never chosen to cancel a measured discrete phase gradient. Radial reactions
enforce the explicitly specified zero radius velocity inside the KWT solve;
velocities are not overwritten after solving. Their work is retained.

This supplies the field-side, fixed-current snapshot of D.27. It does not
implement a moving radius a_b(Is(t)), kinetic reservoir exchange, a normal
boundary trace discretization, or a time integrator. The graph potential must
use the same injected currents as these phase loads. Interior-edge normal
currents are not automatically external-boundary normal fluxes.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .energy_catalog import E_CHARGE_C, HBAR_J_S
from .spatial_dynamics import kwt_spatial_response, SpatialKWTResponse


@dataclass(frozen=True)
class FixedRadiusCurrentLoad:
    load_cartesian_bar: np.ndarray
    radial_load_cartesian_bar: np.ndarray
    phase_load_cartesian_bar: np.ndarray
    terminal_nodes: np.ndarray
    terminal_injection_A: np.ndarray
    phase_load_bar: np.ndarray
    radial_reaction_bar: np.ndarray
    fixed_radii_bar: np.ndarray
    current_scale_A: float


@dataclass(frozen=True)
class ReservoirKWTResponse:
    kwt: SpatialKWTResponse
    load: FixedRadiusCurrentLoad
    terminal_radius_rates_bar: np.ndarray
    radial_reaction_work_rate_bar: float
    phase_reservoir_work_rate_bar: float
    reservoir_work_rate_bar: float
    decomposition_residual_bar: float


def fixed_radius_current_load(delta_bar, gradient_cartesian_bar, terminal_nodes,
                              terminal_injection_A, fixed_radii_bar, *, energy_scale_J):
    """Integrated load with phase component Ltheta=-I_in/current_scale.

    With graph B_tail=+1,B_head=-1, positive left injection gives the negative
    left outward phase load. At fixed radius a>0 the tangential real-gradient
    component is Ltheta/a. The radial component equals the radial energy
    gradient and is the reaction enforcing da/dtau=0, not a prescribed phase
    value. Both independent terminal phases remain dynamical.
    """
    z = np.asarray(delta_bar, complex)
    raw_g, raw_i, raw_r = map(np.asarray, (gradient_cartesian_bar, terminal_injection_A, fixed_radii_bar))
    if any(np.iscomplexobj(v) for v in (raw_g, raw_i, raw_r)):
        raise ValueError("gradients, injections and radii must be real")
    gradient = np.asarray(raw_g, float)
    injection = np.asarray(raw_i, float)
    radii = np.asarray(raw_r, float)
    nodes = np.asarray(terminal_nodes)
    if (z.ndim != 1 or gradient.shape != (len(z), 2) or np.any(~np.isfinite(z))
            or np.any(~np.isfinite(gradient))):
        raise ValueError("finite independent fields and full Cartesian gradients required")
    if (nodes.shape != (2,) or nodes.dtype.kind not in "iu" or nodes[0] == nodes[1]
            or np.any(nodes < 0) or np.any(nodes >= len(z))):
        raise ValueError("two distinct integer terminal DOFs required")
    if (injection.shape != (2,) or radii.shape != (2,) or np.any(~np.isfinite(injection))
            or np.any(~np.isfinite(radii)) or np.any(radii <= 0)):
        raise ValueError("two finite injections and positive fixed radii required")
    if abs(np.sum(injection)) > 64*np.finfo(float).eps*np.sum(abs(injection)):
        raise ValueError("terminal current injections must balance")
    if not np.isfinite(energy_scale_J) or energy_scale_J <= 0:
        raise ValueError("positive finite energy scale required")
    actual_radius = abs(z[nodes])
    if np.any(abs(actual_radius-radii) > 64*np.finfo(float).eps*np.maximum(actual_radius, radii)):
        raise ValueError("snapshot does not satisfy the prescribed radius; no projection performed")
    current_scale = 2*E_CHARGE_C/HBAR_J_S*float(energy_scale_J)
    phase_load = -injection/current_scale
    direction = z[nodes]/actual_radius
    complex_gradient = gradient[:, 0]+1j*gradient[:, 1]
    radial_reaction = np.real(np.conj(direction)*complex_gradient[nodes])
    radial = np.zeros(len(z), complex)
    tangent = np.zeros(len(z), complex)
    radial[nodes] = direction*radial_reaction
    tangent[nodes] = 1j*direction*phase_load/actual_radius
    cartesian = lambda value: np.column_stack((value.real, value.imag))
    return FixedRadiusCurrentLoad(cartesian(radial+tangent), cartesian(radial), cartesian(tangent),
        nodes.astype(np.intp, copy=True), injection.copy(), phase_load, radial_reaction,
        radii.copy(), current_scale)


def fixed_radius_reservoir_kwt(delta_bar, gradient_cartesian_bar, quadrature_mass_bar,
                               mobility, temperatures_bar, phi_V, *, quadrature_to_node=None,
                               terminal_nodes, terminal_injection_A, fixed_radii_bar,
                               energy_scale_J):
    """KWT with explicit fixed-radius reactions and signed reservoir work."""
    load = fixed_radius_current_load(delta_bar, gradient_cartesian_bar, terminal_nodes,
                                    terminal_injection_A, fixed_radii_bar, energy_scale_J=energy_scale_J)
    response = kwt_spatial_response(delta_bar, gradient_cartesian_bar, quadrature_mass_bar,
        mobility, temperatures_bar, phi_V, quadrature_to_node=quadrature_to_node,
        boundary_load_bar=load.load_cartesian_bar)
    z = response.delta_bar
    nodes = load.terminal_nodes
    radius_rates = np.real(np.conj(z[nodes])/abs(z[nodes])*response.material_velocity_bar[nodes])
    complex_load = lambda values: values[:, 0]+1j*values[:, 1]
    radial_work = float(np.real(np.vdot(complex_load(load.radial_load_cartesian_bar), response.material_velocity_bar)))
    phase_work = float(np.real(np.vdot(complex_load(load.phase_load_cartesian_bar), response.material_velocity_bar)))
    return ReservoirKWTResponse(response, load, radius_rates, radial_work, phase_work,
        response.boundary_work_rate_bar, response.boundary_work_rate_bar-radial_work-phase_work)
