"""Conservative discrete electron--phonon events for experimental cell tests.

Energies use Delta0, time uses a declared t_ref, and particle densities use
N0*Delta0. Electronic capacities are 4*w_i. Phonon capacities are
g_ph/N0*d(Omega/Delta0); consequently both sectors share energy N0*Delta0**2.

Barycentric phonon stoichiometry conserves energy exactly. Its weighted Bose
activities preserve detailed balance, but approximate the continuous kernel
at finite phonon resolution. Fractional powers are not Lipschitz at n=0;
physical boundary signs do not promise positivity for arbitrary time steps.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable
import numpy as np

from .energy_catalog import retarded_spectrum


def _positive_vector(values, name: str, *, allow_zero: bool = False):
    array = np.asarray(values, dtype=float)
    if (array.ndim != 1 or len(array) == 0 or np.any(~np.isfinite(array))
            or np.any(array < 0 if allow_zero else array <= 0)):
        raise ValueError(f"{name} must be a finite {'nonnegative' if allow_zero else 'positive'} vector")
    array = array.copy()
    array.flags.writeable = False
    return array


def _readonly(values):
    result = np.asarray(values).copy()
    result.flags.writeable = False
    return result


def _log_nonnegative(values):
    result = np.full_like(values, -np.inf, dtype=float)
    np.log(values, out=result, where=values > 0)
    return result


def _signed_exponential_difference(log_forward, log_reverse):
    """exp(a)-exp(b), retaining detailed-balance cancellation and exact zeros."""
    a, b = np.asarray(log_forward), np.asarray(log_reverse)
    result = np.zeros_like(a)
    positive = a > b
    negative = b > a
    result[positive] = np.exp(a[positive])*(-np.expm1(b[positive]-a[positive]))
    result[negative] = -np.exp(b[negative])*(-np.expm1(a[negative]-b[negative]))
    return result


@dataclass(frozen=True)
class PhononGrid:
    """Fixed positive-energy modes and an explicitly compact coupling function.

    ``coupling_support`` declares where alpha2F may be nonzero; the function is
    evaluated only there. A sharp Debye upper edge is allowed and is not
    silently smoothed. Energies must bracket every event with nonzero coupling.
    No zero-energy mode with divergent Bose occupation is inserted.
    """
    energies: np.ndarray
    capacities: np.ndarray
    alpha2F: Callable[[np.ndarray], np.ndarray]
    coupling_support: tuple[float, float]
    label: str

    def __post_init__(self):
        energies = _positive_vector(self.energies, "phonon energies")
        capacities = _positive_vector(self.capacities, "phonon capacities")
        if len(energies) < 2 or np.any(np.diff(energies) <= 0) or energies.shape != capacities.shape:
            raise ValueError("phonon energies must increase and match capacities")
        support = np.asarray(self.coupling_support, dtype=float)
        if support.shape != (2,) or np.any(~np.isfinite(support)) or not (0 <= support[0] < support[1]):
            raise ValueError("coupling_support must be an explicit finite nonnegative interval")
        if not callable(self.alpha2F) or not isinstance(self.label, str) or not self.label.strip():
            raise ValueError("a coupling function and provenance label are required")
        object.__setattr__(self, "energies", energies)
        object.__setattr__(self, "capacities", capacities)
        object.__setattr__(self, "coupling_support", tuple(float(x) for x in support))

    def coupling(self, omega):
        omega = np.asarray(omega, dtype=float)
        if np.any(~np.isfinite(omega)) or np.any(omega < 0):
            raise ValueError("event energies must be finite and nonnegative")
        inside = (omega >= self.coupling_support[0]) & (omega <= self.coupling_support[1])
        result = np.zeros_like(omega)
        if np.any(inside):
            values = np.asarray(self.alpha2F(omega[inside]), dtype=float)
            if values.shape != omega[inside].shape or np.any(~np.isfinite(values)) or np.any(values < 0):
                raise ValueError("alpha2F must return finite nonnegative values with the query shape")
            result[inside] = values
        return result


class ElectronPhononEvents:
    """One shared list of S/R events updates both sectors with opposite energy.

    Scattering enumerates i<j. Recombination enumerates i<=j, uses half the
    off-diagonal weight for i=j, and removes two particles there. This is the
    unordered quadrature of D.15/B.7, not a second recombination convention.
    """

    def __init__(self, electron_energies, electron_weights, coherence_ratio,
                 phonons: PhononGrid, *, rate_prefactor: float, label: str):
        energy = _positive_vector(electron_energies, "electron energies")
        weights = _positive_vector(electron_weights, "electron weights")
        ratio = np.asarray(coherence_ratio, dtype=float)
        if (energy.shape != weights.shape or ratio.shape != energy.shape
                or np.any(np.diff(energy) <= 0) or np.any(~np.isfinite(ratio))
                or np.any(np.abs(ratio) > 1)):
            raise ValueError("electron arrays must match, energies increase, and |R2/N1| must not exceed one")
        if not isinstance(phonons, PhononGrid):
            raise ValueError("phonons must be an explicit PhononGrid")
        if not np.isfinite(rate_prefactor) or rate_prefactor <= 0:
            raise ValueError("rate_prefactor=8*pi*Delta0*t_ref/hbar must be finite and positive")
        if not isinstance(label, str) or not label.strip():
            raise ValueError("an electronic spectrum provenance label is required")
        si, sj = np.triu_indices(len(energy), 1)
        ri, rj = np.triu_indices(len(energy), 0)
        index_i = np.r_[si, ri]
        index_j = np.r_[sj, rj]
        recombination = np.r_[np.zeros(len(si), dtype=bool), np.ones(len(ri), dtype=bool)]
        omega = np.r_[energy[sj]-energy[si], energy[ri]+energy[rj]]
        coherence = np.r_[1-ratio[si]*ratio[sj], 1+ratio[ri]*ratio[rj]]
        if np.any(coherence < 0):
            raise ValueError("negative coherence factor: incompatible spectral branch")
        coupling = phonons.coupling(omega)
        coefficients = rate_prefactor*weights[index_i]*weights[index_j]*coherence*coupling
        coefficients[recombination & (index_i == index_j)] *= .5
        kept = coefficients > 0
        uncovered = kept & ((omega < phonons.energies[0]) | (omega > phonons.energies[-1]))
        if np.any(uncovered):
            rejected = omega[uncovered]
            raise ValueError(f"{len(rejected)} coupled events lack a phonon bracket: energy range [{rejected.min()}, {rejected.max()}]")
        if np.any(~np.isfinite(coefficients)):
            raise ValueError("event coefficients are not representable as finite floats")
        self.metadata = {
            "schema": "pysnspd.electron_phonon_events.v1", "spectrum_label": label,
            "phonon_label": phonons.label, "rate_prefactor": float(rate_prefactor),
            "rate_prefactor_convention": "8*pi*Delta0*t_ref/hbar", "electron_capacity": "4*w_i",
            "phonon_capacity": "(g_ph/N0)*d(Omega/Delta0)",
            "energy_unit": "N0*Delta0^2", "coupling_support": list(phonons.coupling_support),
            "candidate_scattering": int(len(si)), "candidate_recombination": int(len(ri)),
            "active_scattering": int(np.count_nonzero(kept & ~recombination)),
            "active_recombination": int(np.count_nonzero(kept & recombination)),
            "zero_coupling_events": int(np.count_nonzero(coupling == 0)),
            "zero_coherence_events": int(np.count_nonzero(coherence == 0)),
            "outside_compact_coupling_events": int(np.count_nonzero((omega < phonons.coupling_support[0]) | (omega > phonons.coupling_support[1]))),
            "phonon_interpolation": "barycentric stoichiometry and geometric Bose activities; no occupation clipping",
        }
        self.electron_energies = energy
        self.electron_weights = weights
        self.electron_capacities = _readonly(4*weights)
        self.coherence_ratio = _readonly(ratio)
        self.phonons = phonons
        for name, values in (("index_i", index_i), ("index_j", index_j),
                             ("recombination", recombination), ("omega", omega),
                             ("coefficients", coefficients)):
            setattr(self, name, _readonly(values[kept]))
        upper = np.searchsorted(phonons.energies, self.omega, side="left")
        lower = np.maximum(0, upper-1)
        span = phonons.energies[upper]-phonons.energies[lower]
        beta_upper = np.divide(self.omega-phonons.energies[lower], span,
                               out=np.zeros_like(span), where=span > 0)
        self.phonon_lower, self.phonon_upper = _readonly(lower), _readonly(upper)
        self.beta_upper, self.beta_lower = _readonly(beta_upper), _readonly(1-beta_upper)
        if np.any(beta_upper < 0) or np.any(beta_upper > 1):
            raise ValueError("barycentric event weights are outside their physical interval")

    @classmethod
    def from_arrays(cls, electron_energies, electron_weights, coherence_ratio,
                    phonons: PhononGrid, *, rate_prefactor: float, label: str):
        return cls(electron_energies, electron_weights, coherence_ratio, phonons,
                   rate_prefactor=rate_prefactor, label=label)

    @classmethod
    def from_cell(cls, cell, phonons: PhononGrid, *, rate_prefactor: float):
        c, s = retarded_spectrum(cell.energies, delta=cell.amplitude, gamma=cell.gamma,
                                eta=cell.catalog.eta)
        if np.any(c.real <= 0):
            raise ValueError("catalogue states must have positive causal DOS")
        return cls(cell.energies, cell.weights, s.imag/c.real, phonons,
                   rate_prefactor=rate_prefactor,
                   label=f"R2 fixed-count catalogue: delta={cell.amplitude}, Gamma={cell.gamma}, eta={cell.catalog.eta}")

    def validate(self, p, n):
        p, n = np.asarray(p, dtype=float), np.asarray(n, dtype=float)
        if (p.shape != self.electron_energies.shape or np.any(~np.isfinite(p))
                or np.any(p < 0) or np.any(p > 1)):
            raise ValueError("electronic occupations must be finite and in [0,1]")
        if n.shape != self.phonons.energies.shape or np.any(~np.isfinite(n)) or np.any(n < 0):
            raise ValueError("phonon occupations must be finite and nonnegative")
        return p, n

    def _barycentric_log(self, nodal_logs):
        lower = np.zeros_like(self.omega)
        upper = np.zeros_like(self.omega)
        np.multiply(self.beta_lower, nodal_logs[self.phonon_lower], out=lower, where=self.beta_lower > 0)
        np.multiply(self.beta_upper, nodal_logs[self.phonon_upper], out=upper, where=self.beta_upper > 0)
        return lower+upper

    def log_activities(self, p, n):
        p, n = self.validate(p, n)
        logp, loghole = _log_nonnegative(p), _log_nonnegative(1-p)
        emission = self._barycentric_log(np.log1p(n))
        absorption = self._barycentric_log(_log_nonnegative(n))
        i, j, rec = self.index_i, self.index_j, self.recombination
        forward = np.where(rec, logp[i]+logp[j], loghole[i]+logp[j])+emission
        reverse = np.where(rec, loghole[i]+loghole[j], logp[i]+loghole[j])+absorption
        return forward, reverse

    def rates(self, p, n):
        forward, reverse = self.log_activities(p, n)
        result = self.coefficients*_signed_exponential_difference(forward, reverse)
        if np.any(~np.isfinite(result)):
            raise FloatingPointError("event rates overflowed; no repair was applied")
        return result

    def rhs_from_rates(self, rates):
        rates = np.asarray(rates, dtype=float)
        if rates.shape != self.omega.shape or np.any(~np.isfinite(rates)):
            raise ValueError("one finite signed rate is required per event")
        electronic = np.bincount(self.index_i, weights=np.where(self.recombination, -rates, rates),
                                  minlength=len(self.electron_energies))
        electronic += np.bincount(self.index_j, weights=-rates, minlength=len(self.electron_energies))
        phononic = np.bincount(self.phonon_lower, weights=self.beta_lower*rates,
                               minlength=len(self.phonons.energies))
        phononic += np.bincount(self.phonon_upper, weights=self.beta_upper*rates,
                                minlength=len(self.phonons.energies))
        dp, dn = electronic/self.electron_capacities, phononic/self.phonons.capacities
        if np.any(~np.isfinite(dp)) or np.any(~np.isfinite(dn)):
            raise FloatingPointError("occupation derivatives overflowed; no repair was applied")
        return dp, dn

    def rhs(self, p, n):
        return self.rhs_from_rates(self.rates(p, n))

    def energy(self, p, n):
        p, n = self.validate(p, n)
        return float(np.dot(self.electron_capacities*self.electron_energies, p)
                     +np.dot(self.phonons.capacities*self.phonons.energies, n))

    def energy_rate(self, dp, dn):
        return float(np.dot(self.electron_capacities*self.electron_energies, dp)
                     +np.dot(self.phonons.capacities*self.phonons.energies, dn))
