"""Small electronic-cell experiments using the admitted R2 catalogue.

This module has no production imports or hooks. Energies are normalized by
Delta0, energy densities by N0*Delta0**2, and time by an explicitly chosen
synthetic reference time. It implements D.11, the BGK part of D.17 and D.18;
it does not supply electron-phonon rates or a calibrated KWT mobility. The
closure runner tests a homogeneous condensate with declared synthetic mobility.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.optimize import brentq
from scipy.special import expit, log_expit, xlogy

from .energy_catalog import OccupationEnergyCatalog, retarded_spectrum


def occupation_array(occupation: np.ndarray, shape: tuple[int, ...]) -> np.ndarray:
    """Validate Pauli support without clipping even a small violation."""
    p = np.asarray(occupation, dtype=float)
    if p.shape != shape or np.any(~np.isfinite(p)) or np.any((p < 0) | (p > 1)):
        raise ValueError("occupation must have the declared shape and lie in [0,1]")
    return p


def _positive(value: float, name: str, *, allow_zero: bool = False) -> float:
    value = float(value)
    if not np.isfinite(value) or value < 0 or (value == 0 and not allow_zero):
        raise ValueError(f"{name} must be finite and {'nonnegative' if allow_zero else 'positive'}")
    return value


@dataclass(frozen=True)
class ElectronicCell:
    """A fixed-field view of the real catalogue and its count quadrature."""

    catalog: OccupationEnergyCatalog
    amplitude: float
    gamma: float

    def __post_init__(self) -> None:
        kernels = self.catalog.energy_kernel(self.amplitude, self.gamma)
        for name, values in zip(("energies", "amplitude_kernel", "gamma_kernel"), kernels):
            values = np.array(values, copy=True)
            values.setflags(write=False)
            object.__setattr__(self, name, values)

    @property
    def weights(self) -> np.ndarray:
        return self.catalog.count_weights

    def validate(self, occupation: np.ndarray) -> np.ndarray:
        return occupation_array(occupation, self.energies.shape)

    def excitation_energy(self, occupation: np.ndarray) -> float:
        return float(4 * np.dot(self.weights * self.energies, self.validate(occupation)))

    def moments(self, occupation: np.ndarray) -> tuple[float, float, float]:
        return self.catalog.evaluate(self.amplitude, self.gamma, self.validate(occupation))

    def quasiparticle_count(self, occupation: np.ndarray) -> float:
        return float(4 * np.dot(self.weights, self.validate(occupation)))

    def entropy(self, occupation: np.ndarray) -> float:
        p = self.validate(occupation)
        return float(-4 * np.dot(self.weights, xlogy(p, p) + xlogy(1-p, 1-p)))

    def fermi_dirac(self, temperature: float) -> np.ndarray:
        temperature = _positive(temperature, "temperature", allow_zero=True)
        return np.zeros_like(self.energies) if temperature == 0 else expit(-self.energies / temperature)

    def equivalent_temperature(self, occupation: np.ndarray) -> float:
        """Invert D.11 on this finite count rule; reject an unattainable energy."""
        target = self.excitation_energy(occupation)
        if target == 0:
            return 0.0
        maximum = float(2 * np.dot(self.weights, self.energies))
        if target >= maximum:
            raise ValueError("energy has no finite positive-temperature FD state on this support")

        def residual(temperature: float) -> float:
            if temperature == 0:
                return -target
            return float(4 * np.dot(self.weights * self.energies,
                                   expit(-self.energies / temperature)) - target)

        upper = max(1.0, float(np.max(self.energies)))
        for _ in range(64):
            if residual(upper) > 0:
                return float(brentq(residual, 0.0, upper,
                                    xtol=1e-14, rtol=2e-14))
            upper *= 2
        raise ValueError("temperature inversion could not bracket the finite-table energy")

    def bgk(self, occupation: np.ndarray, tau: float) -> tuple[np.ndarray, float]:
        """Energy-conserving D.17 relaxation; quasiparticle count is not constrained."""
        p = self.validate(occupation)
        tau = _positive(tau, "tau")
        temperature = self.equivalent_temperature(p)
        return (self.fermi_dirac(temperature) - p) / tau, temperature

    def heating(self, occupation: np.ndarray, power: float,
                bath_temperature: float) -> np.ndarray:
        """D.18 heat source, normalized on the same E(x), including p=0."""
        p = self.validate(occupation)
        power = _positive(power, "power", allow_zero=True)
        bath_temperature = _positive(bath_temperature, "bath_temperature")
        if power == 0:
            return np.zeros_like(p)
        if np.all(p == 1):
            raise ValueError("a saturated count support cannot receive positive power")
        temperature = max(self.equivalent_temperature(p), bath_temperature)
        available = p < 1
        log_shape = np.full(p.shape, -np.inf)
        log_shape[available] = (np.log(self.energies[available])
                               + log_expit(-self.energies[available] / temperature)
                               + np.log1p(-p[available]))
        shape = np.exp(log_shape - np.max(log_shape))
        normalization = 4 * np.dot(self.weights * self.energies, shape)
        if not np.isfinite(normalization) or normalization <= 0:
            raise FloatingPointError("heat source has no resolved positive capacity")
        return power * shape / normalization

    def condensate_relaxation(self, occupation: np.ndarray,
                              mobility: float) -> tuple[float, float]:
        """Homogeneous amplitude relaxation and its positive deposited power.

        Mobility is a supplied synthetic constant in the closure tests. This
        helper does not identify it with the calibrated KWT law of D.12.
        """
        mobility = _positive(mobility, "mobility")
        force = self.moments(occupation)[1]
        return -mobility*force, mobility*force*force


def rk4_trajectory(rhs: Callable[[float, np.ndarray], np.ndarray], initial: np.ndarray,
                   duration: float, steps: int) -> tuple[np.ndarray, np.ndarray]:
    """Classical RK4 without clipping, projection or an energy correction."""
    duration = _positive(duration, "duration")
    if not isinstance(steps, (int, np.integer)) or steps < 1:
        raise ValueError("steps must be a positive integer")
    initial = np.asarray(initial, float)
    if np.any(~np.isfinite(initial)):
        raise ValueError("initial state must be finite")
    times = np.linspace(0, duration, steps + 1)
    states = np.empty((steps + 1, *initial.shape))
    states[0] = initial
    dt = duration / steps
    for i, time in enumerate(times[:-1]):
        y = states[i]
        k1 = rhs(float(time), y)
        k2 = rhs(float(time + dt/2), y + dt*k1/2)
        k3 = rhs(float(time + dt/2), y + dt*k2/2)
        k4 = rhs(float(time + dt), y + dt*k3)
        states[i+1] = y + dt*(k1 + 2*k2 + 2*k3 + k4)/6
        if np.any(~np.isfinite(states[i+1])):
            raise FloatingPointError("time integrator produced a nonfinite state")
    return times, states


def linear_count_moments(count: np.ndarray, energy: np.ndarray,
                         shared_energy: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Integrate hats and E*hats on a positive, piecewise-linear E(x).

    No endpoint extrapolation is performed. On each overlapping energy segment
    dx/dE is constant; three-point Gauss integrates the degree-two products
    exactly. The result is a new stated quadrature, not the R2 native Gauss rule.
    """
    x, e, face = (np.asarray(values, float) for values in (count, energy, shared_energy))
    if (x.ndim != 1 or e.shape != x.shape or face.ndim != 1 or len(x) < 2 or len(face) < 2
            or np.any(~np.isfinite(x)) or np.any(~np.isfinite(e)) or np.any(~np.isfinite(face))
            or np.any(np.diff(x) <= 0) or np.any(np.diff(e) <= 0)
            or np.any(np.diff(face) <= 0) or x[0] < 0 or e[0] < 0 or face[0] <= 0):
        raise ValueError("ordered finite count, energy and positive shared-energy axes required")
    moments = np.zeros((2, len(face)))
    nodes, weights = np.polynomial.legendre.leggauss(3)
    for i in range(len(e)-1):
        density = (x[i+1]-x[i])/(e[i+1]-e[i])
        first = max(0, int(np.searchsorted(face, e[i], side="right")-1))
        last = min(len(face)-2, int(np.searchsorted(face, e[i+1], side="left")))
        for j in range(first, last+1):
            low, high = max(e[i], face[j]), min(e[i+1], face[j+1])
            if high <= low:
                continue
            points = low+(nodes+1)*(high-low)/2
            measure = weights*(high-low)*density/2
            right = (points-face[j])/(face[j+1]-face[j])
            for k, hat in ((j, 1-right), (j+1, right)):
                moments[0, k] += np.dot(measure, hat)
                moments[1, k] += np.dot(measure, points*hat)
    return moments[0], moments[1]


@dataclass(frozen=True)
class EnergyFacePair:
    """Energy-conservative two-cell interface with an explicitly lumped capacity.

    The shared primary samples are f(E_j), not equal-index p(x) in different
    spectra. H_j=integral E*hat_j dx and C_j=H_j/E_j make the reconstructed
    electronic energy exactly conservative and give a positive two-state
    generator at each E_j. Quasiparticle count uses M_j=integral hat_j dx;
    its lumping error must be measured and converged. This preparatory scheme
    is not yet the final discretization of all moments of D.17. When one cell
    has a higher upper energy, the last degree of freedom represents its
    private tail as a constant population, disconnected from the interface.
    """

    left: ElectronicCell
    right: ElectronicCell
    shared_energy: np.ndarray
    diffusion_over_length_squared: float = .2

    def __post_init__(self) -> None:
        face = np.asarray(self.shared_energy, float)
        coefficient = _positive(self.diffusion_over_length_squared, "diffusion_over_length_squared")
        if (face.ndim != 1 or len(face) < 3 or np.any(~np.isfinite(face))
                or face[0] <= 0 or np.any(np.diff(face) <= 0)):
            raise ValueError("shared energy grid must be finite, positive and increasing")
        if self.left.catalog.eta != self.right.catalog.eta:
            raise ValueError("the interface requires the same verified spectral regulator")
        if (face[0] > min(self.left.energies[0], self.right.energies[0])
                or face[-1] < max(self.left.energies[-1], self.right.energies[-1])):
            raise ValueError("shared representation must cover every native state without extrapolation")
        common_upper = min(self.left.energies[-1], self.right.energies[-1])
        boundary = np.flatnonzero(face == common_upper)
        if len(boundary) != 1 or np.count_nonzero(face > common_upper) > 1:
            raise ValueError("represent the common upper endpoint explicitly and at most one private-tail degree of freedom")
        common_size = int(boundary[0])+1
        common_face = face[:common_size]
        count_moments, energy_moments, diffusion = [], [], []
        for cell in (self.left, self.right):
            # E(0)=0 is the exact finite-eta endpoint; no state is added there.
            count = np.r_[0., cell.catalog.count_nodes]
            energy = np.r_[0., cell.energies]
            mass, energetic = linear_count_moments(count, energy, common_face)
            if common_size < len(face):
                tail_mass = 0.; tail_energy = 0.
                for xa, xb, ea, eb in zip(count[:-1], count[1:], energy[:-1], energy[1:]):
                    low = max(ea, common_upper)
                    if eb > low:
                        density = (xb-xa)/(eb-ea)
                        tail_mass += density*(eb-low)
                        tail_energy += density*(eb*eb-low*low)/2
                mass = np.r_[mass, tail_mass]
                energetic = np.r_[energetic, tail_energy]
            c, s = retarded_spectrum(face, delta=cell.amplitude, gamma=cell.gamma,
                                     eta=cell.catalog.eta)
            longitudinal = c.real*c.real - s.imag*s.imag
            if np.any(~np.isfinite(longitudinal)) or np.any(longitudinal < 0):
                raise FloatingPointError("longitudinal spectral diffusion is not nonnegative")
            count_moments.append(mass)
            energy_moments.append(energetic)
            diffusion.append(longitudinal)
        mass, energetic, diffusion = map(np.asarray, (count_moments, energy_moments, diffusion))
        capacities = energetic/face[None, :]
        active = (np.all(capacities > 0, axis=0) & (np.sum(diffusion, axis=0) > 0)
                  & (face <= common_upper))
        # Integrate the interface hats only over the common interval. The
        # boundary receives its left half-weight, and the private tail zero.
        face_weights = np.zeros_like(face)
        face_weights[:common_size] = np.r_[np.diff(common_face)[0]/2,
                                          (common_face[2:]-common_face[:-2])/2,
                                          np.diff(common_face)[-1]/2]
        conductance = np.zeros_like(face)
        conductance[active] = (coefficient*face_weights[active]
                               *2*diffusion[0, active]*diffusion[1, active]
                               /(diffusion[0, active]+diffusion[1, active]))
        for name, values in (("shared_energy", face), ("count_moments", mass),
                             ("energy_moments", energetic), ("capacities", capacities),
                             ("longitudinal_diffusion", diffusion), ("conductance", conductance),
                             ("active", active)):
            values = np.array(values, copy=True)
            values.setflags(write=False)
            object.__setattr__(self, name, values)
        object.__setattr__(self, "common_upper_energy", float(common_upper))
        object.__setattr__(self, "common_size", common_size)

    def validate(self, population: np.ndarray) -> np.ndarray:
        return occupation_array(population, self.capacities.shape)

    def rhs(self, population: np.ndarray) -> np.ndarray:
        f = self.validate(population)
        flux = self.conductance*(f[0]-f[1])
        rhs = np.zeros_like(f)
        rhs[0, self.active] = -flux[self.active]/self.capacities[0, self.active]
        rhs[1, self.active] = flux[self.active]/self.capacities[1, self.active]
        return rhs

    def exact(self, population: np.ndarray, time: float) -> np.ndarray:
        """Independent exponential solution of each positive two-cell exchange."""
        f = self.validate(population)
        time = _positive(time, "time", allow_zero=True)
        result = f.copy()
        cl, cr = self.capacities[:, self.active]
        # Convex form avoids cancelling a cold population out of a much larger
        # mean at t=0. expm1 also retains a small, representable transfer.
        transfer = -np.expm1(-self.conductance[self.active]*(1/cl+1/cr)*time)
        left_weight = cr/(cl+cr)*transfer
        right_weight = cl/(cl+cr)*transfer
        result[0, self.active] = (1-left_weight)*f[0, self.active]+left_weight*f[1, self.active]
        result[1, self.active] = (1-right_weight)*f[1, self.active]+right_weight*f[0, self.active]
        return result

    def reconstructed_energy(self, population: np.ndarray) -> np.ndarray:
        return 4*np.sum(self.energy_moments*self.validate(population), axis=1)

    def quasiparticle_count(self, population: np.ndarray) -> float:
        return float(4*np.sum(self.count_moments*self.validate(population)))

    def native_populations(self, population: np.ndarray) -> np.ndarray:
        f = self.validate(population)
        if any(cell.energies[0] < self.shared_energy[0] or cell.energies[-1] > self.shared_energy[-1]
               for cell in (self.left, self.right)):
            raise ValueError("native population reconstruction would extrapolate the represented support")
        populations = []
        for cell, values in zip((self.left, self.right), f):
            common = cell.energies <= self.common_upper_energy
            result = np.empty_like(cell.energies)
            result[common] = np.interp(cell.energies[common], self.shared_energy[:self.common_size],
                                       values[:self.common_size])
            if np.any(~common):
                result[~common] = values[-1]  # the explicitly retained constant private-tail state
            populations.append(result)
        return np.asarray(populations)

    def native_energy(self, population: np.ndarray) -> np.ndarray:
        return np.asarray([cell.excitation_energy(p) for cell, p in
                           zip((self.left, self.right), self.native_populations(population))])
