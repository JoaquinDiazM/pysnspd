"""Energy-consistent algebra for an experimental spatial time RHS.

These helpers implement no time integrator, potential discretization or kinetic
collision kernel. The caller supplies the common-energy gradient, its exact
gauge-link currents, instantaneous electronic spectra, and geometric masses.
In particular Gamma(q) work is already in that gradient. Adding a separate
E_Gamma*Gamma_dot term to the field work, or another spectral drift to p(x),
would count it twice (D.20--21).

KWT is assembled on identified quadrature nodes. Several quadrature records
may share a field DOF, with distinct spectra/temperatures and physical masses.
Their positive temporal matrices are added before solving for the one field
velocity. This supports a one-hot identification map, not a general spatial
interpolation/prolongation matrix. No averaging of equivalent temperatures is
used. A general interpolation requires a global assembled positive solve.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .energy_catalog import E_CHARGE_C, HBAR_J_S


def _real(value, name):
    raw = np.asarray(value)
    if np.iscomplexobj(raw):
        raise ValueError(f"{name} must be real")
    result = np.asarray(raw, float)
    if np.any(~np.isfinite(result)):
        raise ValueError(f"{name} must be finite")
    return result


def _mass(value):
    mass = _real(value, "quadrature_mass_bar")
    if mass.ndim != 1 or not len(mass) or np.any(mass <= 0):
        raise ValueError("quadrature masses must be a nonempty positive vector")
    return mass


def _complex_vector(value, name):
    result = np.asarray(value, complex)
    if result.ndim != 1 or not len(result) or np.any(~np.isfinite(result)):
        raise ValueError(f"{name} must be a nonempty finite complex vector")
    return result


def _gradient(value, size, name):
    result = _real(value, name)
    if result.shape != (size, 2):
        raise ValueError(f"{name} must contain the two integrated Cartesian components per DOF")
    return result[:, 0]+1j*result[:, 1]


def _positive_scalar(value, name):
    result = _real(value, name)
    if result.ndim != 0 or result <= 0:
        raise ValueError(f"{name} must be a positive scalar")
    return float(result)


@dataclass(frozen=True)
class SpatialKWTResponse:
    material_velocity_bar: np.ndarray
    field_velocity_bar: np.ndarray
    heat_density_bar: np.ndarray
    quadrature_mass_bar: np.ndarray
    quadrature_to_node: np.ndarray
    node_mass_bar: np.ndarray
    radial_temporal_weight: np.ndarray
    tangential_temporal_weight: np.ndarray
    condensate_heat_rate_bar: float
    field_energy_rate_bar: float
    boundary_work_rate_bar: float
    electrical_phase_work_rate_bar: float
    identity_residual_bar: float
    delta_bar: np.ndarray
    time_scale_s: float


def equivalent_temperatures(cells, occupations):
    """D.11 on the supplied INSTANTANEOUS spectra; no catalogue reconstruction."""
    if len(cells) != len(occupations) or not len(cells):
        raise ValueError("one occupation row per instantaneous electronic cell required")
    return np.array([cell.equivalent_temperature(p) for cell, p in zip(cells, occupations)])


def kwt_spatial_response(delta_bar, gradient_cartesian_bar, quadrature_mass_bar,
                         mobility, temperatures_bar, phi_V, *,
                         quadrature_to_node=None, boundary_load_bar=None):
    """D.12--13 with a full real integrated gradient and covariant time derivative.

    G is dUbar/d(Re z,Im z), not the half-sized Wirtinger force. The returned
    material velocity is D_tau z; field_velocity is dz/dtau with
    D_tau z=dz/dtau+i*(2e*t_ref/hbar)*phi*z. Positive boundary_load denotes an
    externally applied conjugate load: KWT acts on G-load, while its work is
    Re(load* D_tau z), as in D.31. Overwriting boundary velocities afterwards
    invalidates this identity unless the constraint work is supplied.
    """
    z = _complex_vector(delta_bar, "delta_bar")
    size = len(z)
    gradient = _gradient(gradient_cartesian_bar, size, "gradient_cartesian_bar")
    load = (np.zeros(size, complex) if boundary_load_bar is None else
            _gradient(boundary_load_bar, size, "boundary_load_bar"))
    mass = _mass(quadrature_mass_bar)
    temperature = _real(temperatures_bar, "temperatures_bar")
    phi = _real(phi_V, "phi_V")
    if temperature.shape != mass.shape or np.any(temperature < 0) or phi.shape != (size,):
        raise ValueError("one nonnegative temperature per quadrature and one potential per field DOF required")
    if quadrature_to_node is None:
        if len(mass) != size:
            raise ValueError("different quadrature/DOF counts require an explicit identification map")
        mapping = np.arange(size)
    else:
        mapping = np.asarray(quadrature_to_node)
        if (mapping.shape != mass.shape or mapping.dtype.kind not in "iu"
                or np.any(mapping < 0) or np.any(mapping >= size)):
            raise ValueError("quadrature_to_node must identify exactly one valid DOF per quadrature")
        mapping = mapping.astype(np.intp, copy=True)
    node_mass = np.bincount(mapping, weights=mass, minlength=size)
    if np.any(node_mass <= 0):
        raise ValueError("every field DOF needs positive assembled quadrature mass")
    radial, tangential = np.empty(len(mass)), np.empty(len(mass))
    for j, node in enumerate(mapping):
        coefficient = mobility.coefficients(abs(z[node]), temperature[j])
        denominator = 2*coefficient["Abar"]*coefficient["tau0bar"]
        radial[j] = denominator*coefficient["R"]
        tangential[j] = denominator/coefficient["R"]
    radial_sum = np.bincount(mapping, weights=mass*radial, minlength=size)
    tangential_sum = np.bincount(mapping, weights=mass*tangential, minlength=size)
    direction = np.ones(size, complex)
    nonzero = abs(z) > 0
    direction[nonzero] = z[nonzero]/abs(z[nonzero])
    projected_force = np.conj(direction)*(gradient-load)
    vr = -projected_force.real/radial_sum
    vt = -projected_force.imag/tangential_sum
    material = direction*(vr+1j*vt)
    heat = radial*vr[mapping]**2+tangential*vt[mapping]**2
    gauge_rate = -1j*(2*E_CHARGE_C*mobility.scales.t_ref_s/HBAR_J_S)*phi*z
    velocity = material+gauge_rate
    heat_rate = float(np.dot(mass, heat))
    field_rate = float(np.real(np.vdot(gradient, velocity)))
    boundary_work = float(np.real(np.vdot(load, material)))
    phase_work = float(np.real(np.vdot(gradient, gauge_rate)))
    if not np.all(np.isfinite(velocity)) or not np.all(np.isfinite(heat)):
        raise FloatingPointError("nonfinite spatial KWT response")
    return SpatialKWTResponse(material, velocity, heat, mass.copy(), mapping,
        node_mass, radial_sum, tangential_sum, heat_rate, field_rate,
        boundary_work, phase_work, field_rate+heat_rate-boundary_work-phase_work,
        z.copy(), mobility.scales.t_ref_s)


@dataclass(frozen=True)
class GaugePowerCheck:
    noether_residual_bar: np.ndarray
    field_phase_power_W: float
    edge_superconducting_power_W: float
    power_residual_W: float


def gauge_power_check(delta_bar, gradient_cartesian_bar, edges, link_current_bar,
                      phi_V, *, energy_scale_J, time_scale_s):
    """Check the required g_theta+B*J=0 and its electrical power identity.

    Edge (tail,head) uses B_tail=+1, B_head=-1. Currents must be exact gauge-link
    derivatives of this energy; a separately sampled continuum js need not
    satisfy this discrete identity. No current or potential is repaired here.
    """
    z = _complex_vector(delta_bar, "delta_bar")
    gradient = _gradient(gradient_cartesian_bar, len(z), "gradient_cartesian_bar")
    edge = np.asarray(edges)
    current = _real(link_current_bar, "link_current_bar")
    phi = _real(phi_V, "phi_V")
    energy_scale = _positive_scalar(energy_scale_J, "energy_scale_J")
    time_scale = _positive_scalar(time_scale_s, "time_scale_s")
    if (edge.ndim != 2 or edge.shape[1] != 2 or edge.dtype.kind not in "iu"
            or np.any(edge < 0) or np.any(edge >= len(z))
            or current.shape != (len(edge),) or phi.shape != (len(z),)):
        raise ValueError("edges, currents and potentials must refer to the same field DOFs")
    divergence = np.zeros(len(z))
    np.add.at(divergence, edge[:, 0], current)
    np.add.at(divergence, edge[:, 1], -current)
    phase_gradient = np.imag(np.conj(z)*gradient)
    gauge_velocity = -1j*(2*E_CHARGE_C*time_scale/HBAR_J_S)*phi*z
    field = energy_scale/time_scale*np.real(np.vdot(gradient, gauge_velocity))
    edge_power = (2*E_CHARGE_C/HBAR_J_S)*energy_scale*np.dot(current, phi[edge[:, 0]]-phi[edge[:, 1]])
    return GaugePowerCheck(phase_gradient+divergence, float(field), float(edge_power), float(field-edge_power))


@dataclass(frozen=True)
class NormalHeatDistribution:
    face_power_W: np.ndarray
    quadrature_power_W: np.ndarray
    total_power_W: float
    allocation_residual_W: float


def normal_heat_distribution(conductance_S, delta_phi_V, allocation_q_face):
    """D.18 Joule power using an explicit geometric, conservative allocation.

    Each column assigns one electrical face's G*(delta_phi)**2 to quadrature
    volumes. Columns must be nonnegative and sum to one; they are not silently
    normalized. Conserving the total does not certify the allocation's local
    spatial accuracy. Use voltage drops from the potential solver, avoiding
    loss of precision from an irrelevant large potential reference.
    """
    conductance = _real(conductance_S, "conductance_S")
    drop = _real(delta_phi_V, "delta_phi_V")
    allocation = _real(allocation_q_face, "allocation_q_face")
    if (conductance.ndim != 1 or not len(conductance) or np.any(conductance <= 0)
            or drop.shape != conductance.shape or allocation.ndim != 2
            or not allocation.shape[0] or allocation.shape[1] != len(conductance)
            or np.any(allocation < 0)):
        raise ValueError("positive conductances and a nonnegative quadrature-by-face allocation required")
    sums = np.sum(allocation, axis=0)
    if np.any(abs(sums-1) > 64*np.finfo(float).eps*max(1, allocation.shape[0])):
        raise ValueError("each face heat allocation must sum to one; no normalization applied")
    power = conductance*drop*drop
    q_power = allocation@power
    if np.any(~np.isfinite(power)) or np.any(~np.isfinite(q_power)):
        raise FloatingPointError("nonfinite normal heating power")
    return NormalHeatDistribution(power, q_power, float(np.sum(power)), float(np.sum(q_power)-np.sum(power)))


@dataclass(frozen=True)
class SpatialHeatDeposition:
    electron_rhs: tuple[np.ndarray, ...]
    source_density_bar: np.ndarray
    deposited_density_rate_bar: np.ndarray
    condensate_rate_bar: float
    normal_rate_bar: float
    deposited_rate_bar: float
    moment_residual_bar: float


def deposit_spatial_heat(cells, occupations, kwt_response, normal_heat_power_W, *,
                         energy_scale_J, scales):
    """Deposit QDelta plus normal Joule power exactly once, using D.18 unchanged."""
    mass = kwt_response.quadrature_mass_bar
    normal = _real(normal_heat_power_W, "normal_heat_power_W")
    scale = _positive_scalar(energy_scale_J, "energy_scale_J")
    if (len(cells) != len(mass) or len(occupations) != len(mass)
            or normal.shape != mass.shape or np.any(normal < 0)):
        raise ValueError("one electronic cell, population and nonnegative normal power per quadrature required")
    if scales.t_ref_s != kwt_response.time_scale_s:
        raise ValueError("KWT and heat source must use the same reference time")
    density = kwt_response.heat_density_bar+normal*scales.t_ref_s/(scale*mass)
    rhs, moments = [], []
    for i, (cell, p) in enumerate(zip(cells, occupations)):
        vacuum = cell.catalog.vacuum
        if vacuum.delta0_J != scales.delta0_J or vacuum.N0_per_J_m3 != scales.N0_per_J_m3:
            raise ValueError("electronic cells and dynamical scales must use the same energy/DOS units")
        amplitude = abs(kwt_response.delta_bar[kwt_response.quadrature_to_node[i]])
        if not np.isclose(cell.amplitude, amplitude, rtol=2e-13, atol=0.):
            raise ValueError("heating cell amplitude does not match the instantaneous field")
        value = cell.heating(p, density[i], scales.bath_temperature_bar)
        rhs.append(value)
        moments.append(4*np.dot(cell.weights*cell.energies, value))
    moments = np.asarray(moments)
    normal_rate = float(np.sum(normal)*scales.t_ref_s/scale)
    deposited = float(np.dot(mass, moments))
    return SpatialHeatDeposition(tuple(rhs), density, moments,
        kwt_response.condensate_heat_rate_bar, normal_rate, deposited,
        deposited-kwt_response.condensate_heat_rate_bar-normal_rate)


@dataclass(frozen=True)
class SpatialEnergyRate:
    field_rate_bar: float
    explicit_link_work_bar: float
    electron_population_rate_bar: float
    phonon_rate_bar: float
    total_rate_bar: float


def spatial_energy_rate(gradient_cartesian_bar, field_velocity_bar, cells,
                        electron_rhs, quadrature_mass_bar, *,
                        phonon_density_rates_bar=None, link_current_bar=None,
                        link_phase_velocity=None):
    """Chain rule for the SAME instantaneous energy, including moving Gamma.

    Gamma work is in G*dDelta/dtau. The population term is only 4*w*E*dp/dtau
    at fixed count. Optional alpha_dot represents explicitly time-dependent
    gauge links; static A=0 uses neither optional link argument.
    """
    velocity = _complex_vector(field_velocity_bar, "field_velocity_bar")
    gradient = _gradient(gradient_cartesian_bar, len(velocity), "gradient_cartesian_bar")
    mass = _mass(quadrature_mass_bar)
    if len(cells) != len(mass) or len(electron_rhs) != len(mass):
        raise ValueError("one electronic RHS and instantaneous spectrum per quadrature required")
    population = 0.
    for weight, cell, dp in zip(mass, cells, electron_rhs):
        dp = _real(dp, "electron_rhs")
        if dp.shape != cell.energies.shape:
            raise ValueError("electronic RHS must match its instantaneous spectrum")
        population += weight*4*np.dot(cell.weights*cell.energies, dp)
    if phonon_density_rates_bar is None:
        phonon = 0.
    else:
        rates = _real(phonon_density_rates_bar, "phonon_density_rates_bar")
        if rates.shape != mass.shape:
            raise ValueError("one phonon energy-density rate per quadrature required")
        phonon = float(np.dot(mass, rates))
    link = 0.
    if (link_current_bar is None) != (link_phase_velocity is None):
        raise ValueError("both explicit link current and link phase velocity are required")
    if link_current_bar is not None:
        current = _real(link_current_bar, "link_current_bar")
        rate = _real(link_phase_velocity, "link_phase_velocity")
        if current.ndim != 1 or rate.shape != current.shape:
            raise ValueError("one velocity per explicit gauge link required")
        link = float(np.dot(current, rate))
    field = float(np.real(np.vdot(gradient, velocity)))
    return SpatialEnergyRate(field, link, float(population), phonon, float(field+link+population+phonon))


@dataclass(frozen=True)
class SpatialCircuitPowerBalance:
    domain_energy_rate_W: float
    combined_energy_rate_W: float
    expected_external_rate_W: float
    device_port_residual_W: float
    circuit_residual_W: float
    residual_W: float


def cm9_power_balance(domain_rate_bar, *, energy_scale_J, time_scale_s,
                      circuit_balance, escape_power_W=0., reservoir_power_into_W=0.):
    """CM.9, retaining signed escape and reservoir work and no duplicate QDelta.

    ``reservoir_power_into_W`` contains the condensate traction work plus the
    incoming quasiparticle energy flux. Circuit resistor losses remain outside
    the film; they are not added to D.18. Escape can be negative for a phonon
    population below the bath. This reports a residual, never repairs one.
    """
    scale = _positive_scalar(energy_scale_J, "energy_scale_J")
    time = _positive_scalar(time_scale_s, "time_scale_s")
    values = [domain_rate_bar, escape_power_W, reservoir_power_into_W,
              circuit_balance.stored_energy_rate_W, circuit_balance.device_power_W,
              circuit_balance.source_power_W, circuit_balance.bias_joule_W,
              circuit_balance.load_joule_W, circuit_balance.residual_W]
    if np.asarray(values).shape != (9,) or np.any(~np.isfinite(values)):
        raise ValueError("finite scalar domain, boundary and circuit powers required")
    domain = float(domain_rate_bar*scale/time)
    combined = domain+circuit_balance.stored_energy_rate_W
    expected = (circuit_balance.source_power_W-circuit_balance.bias_joule_W
                -circuit_balance.load_joule_W-escape_power_W+reservoir_power_into_W)
    port = domain-circuit_balance.device_power_W+escape_power_W-reservoir_power_into_W
    return SpatialCircuitPowerBalance(domain, float(combined), float(expected), float(port),
                                      float(circuit_balance.residual_W), float(combined-expected))
