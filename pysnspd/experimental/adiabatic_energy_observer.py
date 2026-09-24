"""Independent SI energy observations for a candidate adiabatic kinetic model.

No state is evolved, projected or repaired here. The caller supplies the
off-shell thermal INTERNAL energy and its directional derivative; a thermal
free energy is not an interchangeable input. The electronic correction is the
neutral, electron-hole-symmetric spectral moment from EC.2 in the stage4
energy/circuit derivation. A charge/gauge completion must be supplied separately
when the dynamical approximation requires one; this module does not invent it.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np


def _real(value, name):
    raw = np.asarray(value)
    if np.iscomplexobj(raw):
        raise ValueError(name + ' must be real')
    out = np.asarray(raw, float)
    if np.any(~np.isfinite(out)):
        raise ValueError(name + ' must be finite')
    return out


def _scalar(value, name):
    out = _real(value, name)
    if out.ndim:
        raise ValueError(name + ' must be scalar')
    return float(out)


@dataclass(frozen=True)
class ElectronicEnergy:
    thermal_internal_J: float
    distribution_increment_J: float
    total_J: float


@dataclass(frozen=True)
class ElectronicEnergyRate:
    thermal_internal_rate_W: float
    moving_dos_rate_W: float
    distribution_rate_W: float
    total_rate_W: float


class AdiabaticElectronicEnergy:
    """Ueq - 2 N0 integral(V E rho [hL-h0] dE), with energies in joules.

    Arrays rho, hL and their derivatives have shape (energies, spatial nodes).
    N0 is the SINGLE-SPIN normal DOS in J^-1 m^-3. Energy weights are positive
    quadrature measures in J, volumes in m^3. h0 is a fixed thermal reference
    at these energy nodes, not recomputed from a fitted instantaneous T.
    There is no numerical DOS floor and no Pauli projection. Signed
    distribution directions can be observed without pretending to be finite
    admitted occupations; state-support admission belongs to the integrator.
    """
    def __init__(self, energy_J, weights_J, volume_m3, N0_per_J_m3, h0):
        self.energies = _real(energy_J, 'energy_J').copy()
        self.weights = _real(weights_J, 'weights_J').copy()
        self.volume = _real(volume_m3, 'volume_m3').copy()
        self.N0 = _scalar(N0_per_J_m3, 'N0_per_J_m3')
        if (self.energies.ndim != 1 or not self.energies.size
                or self.weights.shape != self.energies.shape
                or self.volume.ndim != 1 or not self.volume.size
                or np.any(self.energies < 0) or np.any(np.diff(self.energies) <= 0)
                or np.any(self.weights <= 0) or np.any(self.volume <= 0)
                or self.N0 <= 0):
            raise ValueError('Ordered nonnegative energies and positive measures/N0 required')
        reference = _real(h0, 'h0')
        if reference.shape == self.energies.shape:
            reference = np.broadcast_to(reference[:, None], (len(self.energies), len(self.volume)))
        self.shape = (len(self.energies), len(self.volume))
        if reference.shape != self.shape:
            raise ValueError('h0 must have one entry per energy or per energy/node')
        self.h0 = reference.copy()
        self.measure = 2*self.N0*(self.weights*self.energies)[:, None]*self.volume[None, :]
        for value in (self.energies, self.weights, self.volume, self.h0, self.measure):
            value.setflags(write=False)

    def _field(self, value, name, *, dos=False):
        result = _real(value, name)
        if result.shape != self.shape:
            raise ValueError(name + ' must have shape (energies, nodes)')
        if dos and np.any(result < 0):
            raise ValueError('Negative DOS is not admitted or clipped')
        return result

    def energy(self, thermal_internal_J, rho, hL):
        thermal = _scalar(thermal_internal_J, 'thermal_internal_J')
        dos = self._field(rho, 'rho', dos=True)
        distribution = self._field(hL, 'hL')
        increment = -float(np.sum(self.measure*dos*(distribution-self.h0)))
        return ElectronicEnergy(thermal, increment, thermal+increment)

    def rate(self, thermal_internal_rate_W, rho, hL, rho_dot_per_s, hL_dot_per_s):
        """Chain rule on a FIXED physical-energy quadrature, not count labels.

        The caller's hLdot must include the kinetic level drift required by
        its own equations. The rhoDot term cannot substitute for that drift.
        The stored h0, weights, energies and volumes remain fixed.
        """
        thermal = _scalar(thermal_internal_rate_W, 'thermal_internal_rate_W')
        dos = self._field(rho, 'rho', dos=True)
        distribution = self._field(hL, 'hL')
        dos_dot = self._field(rho_dot_per_s, 'rho_dot_per_s')
        h_dot = self._field(hL_dot_per_s, 'hL_dot_per_s')
        moving = -float(np.sum(self.measure*dos_dot*(distribution-self.h0)))
        population = -float(np.sum(self.measure*dos*h_dot))
        return ElectronicEnergyRate(thermal, moving, population, thermal+moving+population)


def phonon_energy_J(energy_J, weights_J, density_per_J_m3, occupation, volume_m3):
    """Positive phonon occupation energy, with no zero-point contribution.

    Density can be one curve or one curve per spatial node. The same volume
    ownership must be used at a 2D/1D interface; it is not repaired here.
    """
    energy = _real(energy_J, 'phonon energy')
    weights = _real(weights_J, 'phonon weights')
    volume = _real(volume_m3, 'phonon volume')
    if (energy.ndim != 1 or not energy.size or weights.shape != energy.shape
            or volume.ndim != 1 or not volume.size or np.any(energy < 0)
            or np.any(np.diff(energy) <= 0) or np.any(weights <= 0) or np.any(volume <= 0)):
        raise ValueError('Positive phonon measures and ordered nonnegative energies required')
    shape = (len(energy), len(volume))
    density = _real(density_per_J_m3, 'phonon density')
    if density.shape == energy.shape:
        density = np.broadcast_to(density[:, None], shape)
    population = _real(occupation, 'phonon occupation')
    if density.shape != shape or population.shape != shape or np.any(density < 0) or np.any(population < 0):
        raise ValueError('Nonnegative phonon density/occupation on every energy and node required')
    return float(np.sum((energy*weights)[:, None]*density*population*volume[None, :]))


@dataclass(frozen=True)
class GlobalPower:
    """SI powers for CM.9; positive signs are stated in each field name.

    KWT heating, internal e-ph exchange and device Is*Vdev are excluded from
    this external balance: each cancels between two included subsystems.
    """
    source_into_W: float
    bias_resistor_out_W: float
    load_resistor_out_W: float
    phonon_escape_out_W: float = 0.
    contact_heat_out_W: float = 0.
    reservoir_work_into_W: float = 0.

    def net_into_W(self):
        values = {name: _scalar(getattr(self, name), name) for name in self.__dataclass_fields__}
        # Contact heat and escape can be signed because a bath can inject heat.
        if values['bias_resistor_out_W'] < 0 or values['load_resistor_out_W'] < 0:
            raise ValueError('Passive resistor dissipation cannot be negative')
        return (values['source_into_W']-values['bias_resistor_out_W']
                -values['load_resistor_out_W']-values['phonon_escape_out_W']
                -values['contact_heat_out_W']+values['reservoir_work_into_W'])


def instantaneous_balance_residual_W(electronic_rate_W, phonon_rate_W, circuit_rate_W, powers):
    """Independent stored-energy derivative minus supplied external power."""
    if not isinstance(powers, GlobalPower):
        raise TypeError('powers must be GlobalPower')
    rate = sum(_scalar(value, name) for value, name in
               ((electronic_rate_W, 'electronic_rate_W'), (phonon_rate_W, 'phonon_rate_W'),
                (circuit_rate_W, 'circuit_rate_W')))
    return rate-powers.net_into_W()


def integrated_balance(times_s, total_energy_J, net_external_into_W):
    """Accepted-step energy difference minus trapezoidal external work.

    Returns cumulative work and residual arrays. It never defines or changes
    energy from the supplied powers, and does not infer an acceptance scale.
    """
    times = _real(times_s, 'times_s')
    energy = _real(total_energy_J, 'total_energy_J')
    power = _real(net_external_into_W, 'net_external_into_W')
    if (times.ndim != 1 or len(times) < 2 or energy.shape != times.shape
            or power.shape != times.shape or np.any(np.diff(times) <= 0)):
        raise ValueError('Matching histories at strictly increasing accepted times required')
    work = np.r_[0., np.cumsum(.5*np.diff(times)*(power[1:]+power[:-1]))]
    return work, energy-energy[0]-work
