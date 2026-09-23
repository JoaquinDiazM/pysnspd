"""Fixed FD reservoir face for the existing native electronic transport events.

The active cell and bath occupations meet at the SAME event energy, through
NativeTransportEvents' geometric Pauli activities and barycentric deposition.
They are not subtracted at equal state-count indices when spectra differ.
The reservoir distribution is fixed; its opposite energy/number ledger is
recorded, rather than evolving an invented finite reservoir volume.

Geometry sets D*t_ref*A_face/(V_cell*distance), the finite-volume face factor.
Both sides of the internal event pair use the active-cell normalization. Only
the active RHS evolves, and multiplying its moments by V_cell gives the
absolute bath exchange. This wrapper does not support unequal finite cells,
different material energy scales, charge/potential transport, or phase work.
It provides an RHS, not a positivity guarantee for an arbitrary time step.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .cell_transport import NativeTransportEvents
from .cell_validation import ElectronicCell


@dataclass(frozen=True)
class ReservoirFaceGeometry:
    cell_volume_m3: float
    face_area_m2: float
    center_to_reservoir_m: float
    time_scale_s: float

    def __post_init__(self):
        for name in ("cell_volume_m3", "face_area_m2", "center_to_reservoir_m", "time_scale_s"):
            value = float(getattr(self, name))
            if not np.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive")
            object.__setattr__(self, name, value)

    def transport_coefficient(self, diffusivity_m2_s: float) -> float:
        """Rate in tau=t/t_ref; distance must be chosen explicitly by the caller."""
        diffusion = float(diffusivity_m2_s)
        denominator = self.cell_volume_m3*self.center_to_reservoir_m
        if not np.isfinite(diffusion) or diffusion <= 0 or not np.isfinite(denominator) or denominator <= 0:
            raise ValueError("diffusivity and the geometric denominator must be finite and positive")
        coefficient = diffusion*self.time_scale_s*self.face_area_m2/denominator
        if not np.isfinite(coefficient) or coefficient <= 0:
            raise ValueError("the supplied geometry gives an unresolved transport coefficient")
        return float(coefficient)


@dataclass(frozen=True)
class ReservoirTransportExchange:
    cell_rhs: np.ndarray
    event_flux_into_cell: np.ndarray
    number_rate_bar_into_cell: float
    energy_rate_bar_into_cell: float
    reservoir_number_rate_bar: float
    reservoir_energy_rate_bar: float
    cell_quasiparticle_rate_per_s: float
    cell_power_W: float
    reservoir_quasiparticle_rate_per_s: float
    reservoir_power_W: float
    event_number_moment_residual: float
    event_energy_moment_residual: float
    shared_energy_support_bar: tuple[float, float]


@dataclass(frozen=True)
class ElectronicReservoirTransport:
    """A fixed-field electronic boundary, with signs positive INTO the active cell.

    ``reservoir_cell`` is a spectral view, not a finite bath control volume.
    D.26's stable thermal branch can supply its amplitude/Gamma; that selection
    is the caller's responsibility. The fixed FD occupation is evaluated on
    the reservoir's own energies at the explicit bath_theta=kBT/Delta0.

    Equal catalogue units and count quadratures are required because the
    existing native paired-RHS interface uses equal-length population arrays.
    This is an explicit current implementation domain, not a physical claim
    that distinct spectra have identical occupations at equal count labels.
    """
    cell: ElectronicCell
    reservoir_cell: ElectronicCell
    bath_theta: float
    geometry: ReservoirFaceGeometry
    gauss_order: int = 2

    def __post_init__(self):
        theta = float(self.bath_theta)
        if not np.isfinite(theta) or theta < 0:
            raise ValueError("bath_theta must be finite and nonnegative")
        if not isinstance(self.geometry, ReservoirFaceGeometry):
            raise TypeError("geometry must explicitly describe the active cell and reservoir face")
        if type(self.gauss_order) is not int or self.gauss_order < 1:
            raise ValueError("gauss_order must be a positive integer")
        left, right = self.reservoir_cell.catalog, self.cell.catalog
        for name in ("delta0_J", "N0_per_J_m3", "D_m2_s"):
            if getattr(left.vacuum, name) != getattr(right.vacuum, name):
                raise ValueError(f"reservoir and cell must share {name}; unit/material conversion is not implemented")
        if left.eta != right.eta:
            raise ValueError("reservoir and cell must share the same causal regulator")
        if (not np.array_equal(left.count_nodes, right.count_nodes)
                or not np.array_equal(left.count_weights, right.count_weights)):
            raise ValueError("the existing paired transport RHS requires the same native count quadrature")
        coefficient = self.geometry.transport_coefficient(right.vacuum.D_m2_s)
        # Logical left/right here means bath/active, irrespective of which
        # physical terminal this face occupies. Electrical port signs are separate.
        events = NativeTransportEvents(self.reservoir_cell, self.cell,
                                       diffusion_over_length_squared=coefficient,
                                       gauss_order=self.gauss_order)
        p = self.reservoir_cell.fermi_dirac(theta)
        p.setflags(write=False)
        object.__setattr__(self, "bath_theta", theta)
        object.__setattr__(self, "events", events)
        object.__setattr__(self, "reservoir_population", p)
        object.__setattr__(self, "coefficient_per_reference_time", coefficient)

    def evaluate(self, cell_occupation) -> ReservoirTransportExchange:
        """RHS and conservative opposite ledger, with no bath update or clipping.

        Normalized number is in N0*Delta0*V_cell and normalized energy in
        N0*Delta0**2*V_cell, each per reference time. SI moments are absolute
        quasiparticles/s and watts. Event-versus-native moment residuals are
        exposed as diagnostics; neither the RHS nor populations are repaired.
        Only the explicitly intersected native energy support participates.
        """
        p = self.cell.validate(cell_occupation)
        populations = (self.reservoir_population, p)
        paired_rhs = self.events.rhs(populations)
        flux = self.events.rates(populations)
        rhs = np.asarray(paired_rhs[1], float)
        number = float(4*np.dot(self.cell.weights, rhs))
        energy = float(4*np.dot(self.cell.weights*self.cell.energies, rhs))
        event_number = float(np.sum(flux))
        event_energy = float(np.dot(self.events.energies, flux))
        if np.any(~np.isfinite(rhs)) or np.any(~np.isfinite(flux)):
            raise FloatingPointError("reservoir transport produced a nonfinite native RHS")
        vacuum = self.cell.catalog.vacuum
        number_scale = vacuum.N0_per_J_m3*vacuum.delta0_J*self.geometry.cell_volume_m3/self.geometry.time_scale_s
        energy_scale = number_scale*vacuum.delta0_J
        number_si, energy_si = number*number_scale, energy*energy_scale
        return ReservoirTransportExchange(rhs, np.asarray(flux, float), number, energy,
            -number, -energy, float(number_si), float(energy_si), -float(number_si), -float(energy_si),
            number-event_number, energy-event_energy, self.events.shared_support)
