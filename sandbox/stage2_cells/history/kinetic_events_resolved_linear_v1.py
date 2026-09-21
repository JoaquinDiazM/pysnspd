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
from numpy.polynomial.legendre import leggauss

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
                 phonons: PhononGrid, *, rate_prefactor: float, label: str,
                 _admission_only: bool = False):
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
        self.electron_energies = energy
        self.electron_weights = weights
        self.electron_capacities = _readonly(4*weights)
        self.coherence_ratio = _readonly(ratio)
        self.phonons = phonons
        if _admission_only:
            self.metadata = {"spectrum_label": label, "scope": "array admission only; no quadratic native-pair list constructed"}
            return
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


class ProjectedElectronPhononEvents(ElectronPhononEvents):
    """Resolve coupling cutoffs while retaining the native R2 population vector.

    Positive nested quadrature is performed in the reconstructed count
    coordinate. Each virtual electron contributes barycentric stoichiometry
    to two native energies. Geometric Pauli activities use those same weights,
    so both particle/energy moments and FD/BE balance are retained. This is a
    convergent discretization, not an assertion that the finite reconstruction
    equals the continuous D.14 kernel. No population or energy is projected
    after a step. At mixed empty/full endpoints an activity can vanish; the
    original unnormalized products are deliberately retained.
    """

    @classmethod
    def from_arrays(cls, electron_energies, electron_weights, count_nodes, coherence_ratio,
                    phonons: PhononGrid, *, rate_prefactor: float,
                    quadrature_order: int = 64, label: str,
                    integration_layout: str = "global", outer_order: int = 2,
                    max_count_panel: float = .5, max_energy_panel: float = .125):
        return cls(electron_energies, electron_weights, count_nodes, coherence_ratio, phonons,
                   rate_prefactor=rate_prefactor, quadrature_order=quadrature_order, label=label,
                   integration_layout=integration_layout, outer_order=outer_order, max_count_panel=max_count_panel,
                   max_energy_panel=max_energy_panel)

    @classmethod
    def from_cell(cls, cell, phonons: PhononGrid, *, rate_prefactor: float,
                  quadrature_order: int = 64, integration_layout: str = "global",
                  outer_order: int = 2, max_count_panel: float = .5, max_energy_panel: float = .125):
        c, s = retarded_spectrum(cell.energies, delta=cell.amplitude, gamma=cell.gamma,
                                eta=cell.catalog.eta)
        if np.any(c.real <= 0):
            raise ValueError("catalogue states must have positive causal DOS")
        return cls(cell.energies, cell.weights, cell.catalog.count_nodes, s.imag/c.real,
                   phonons, rate_prefactor=rate_prefactor, quadrature_order=quadrature_order,
                   label=f"R2 fixed-count catalogue: delta={cell.amplitude}, Gamma={cell.gamma}, eta={cell.catalog.eta}",
                   integration_layout=integration_layout, outer_order=outer_order, max_count_panel=max_count_panel,
                   max_energy_panel=max_energy_panel)

    def __init__(self, electron_energies, electron_weights, count_nodes, coherence_ratio,
                 phonons: PhononGrid, *, rate_prefactor: float, quadrature_order: int,
                 label: str, ultraviolet_partition: tuple[float, float] | None = None,
                 integration_layout: str = "global", outer_order: int = 2,
                 max_count_panel: float = .5, max_energy_panel: float = .125):
        if type(quadrature_order) is not int or quadrature_order < 2:
            raise ValueError("quadrature_order must be an integer at least two")
        if (integration_layout not in ("global", "native_intervals", "energy_panels", "resolved_panels") or type(outer_order) is not int
                or outer_order < 1 or not np.isfinite(max_count_panel) or max_count_panel <= 0
                or not np.isfinite(max_energy_panel) or max_energy_panel <= 0):
            raise ValueError("invalid positive count-panel quadrature configuration")
        # Reuse admission and capacities, then replace only the event quadrature.
        super().__init__(electron_energies, electron_weights, coherence_ratio, phonons,
                         rate_prefactor=rate_prefactor, label=label,
                         _admission_only=ultraviolet_partition is None)
        counts = _positive_vector(count_nodes, "native count nodes")
        if counts.shape != self.electron_energies.shape or np.any(np.diff(counts) <= 0):
            raise ValueError("native count nodes must increase and match the electron grid")
        self.count_nodes = counts
        node, weight = leggauss(quadrature_order)
        outer_node, outer_weight_rule = leggauss(outer_order)
        energy = self.electron_energies
        xmin, xmax, emin, emax = counts[0], counts[-1], energy[0], energy[-1]
        omega_low, omega_high = phonons.coupling_support
        channels = []
        for rec in (False, True):
            # For recombination integrate only Ei<=Ej, avoiding double counting.
            outer_max_energy = min(emax, omega_high/2) if rec else emax
            if outer_max_energy <= emin:
                continue
            outer_max_count = np.interp(outer_max_energy, energy, counts)
            if integration_layout != "global":
                bounds = np.r_[counts[counts < outer_max_count], outer_max_count]
                outer = (bounds[:-1, None]+np.diff(bounds)[:, None]*(outer_node[None, :]+1)/2).ravel()
                outer_weight = (np.diff(bounds)[:, None]*outer_weight_rule[None, :]/2).ravel()
            else:
                outer = (outer_max_count-xmin)*(node+1)/2+xmin
                outer_weight = (outer_max_count-xmin)*weight/2
            ei = np.interp(outer, counts, energy)
            if rec:
                inner_low_energy = np.maximum(ei, omega_low-ei)
                inner_high_energy = np.minimum(emax, omega_high-ei)
            else:
                inner_low_energy = ei+omega_low
                inner_high_energy = np.minimum(emax, ei+omega_high)
            valid = inner_high_energy > inner_low_energy
            if not np.any(valid):
                continue
            ei, outer, outer_weight = ei[valid], outer[valid], outer_weight[valid]
            low = np.interp(inner_low_energy[valid], energy, counts)
            high = np.interp(inner_high_energy[valid], energy, counts)
            if integration_layout != "global":
                if integration_layout == "resolved_panels":
                    first_inside = np.searchsorted(counts, low, side="right")
                    last_inside = np.searchsorted(counts, high, side="left")
                    count_pieces = last_inside-first_inside+1
                    count_owner = np.repeat(np.arange(len(low)), count_pieces)
                    first_offset = np.r_[0, np.cumsum(count_pieces)[:-1]]
                    number = np.arange(len(count_owner))-np.repeat(first_offset, count_pieces)
                    left_index = np.maximum(first_inside[count_owner]+number-1, 0)
                    right_index = np.minimum(first_inside[count_owner]+number, len(counts)-1)
                    count_start = np.where(number == 0, low[count_owner], counts[left_index])
                    count_end = np.where(number == count_pieces[count_owner]-1,
                                         high[count_owner], counts[right_index])
                    e_start = np.interp(count_start, counts, energy)
                    e_end = np.interp(count_end, counts, energy)
                    omega_start = e_start+ei[count_owner] if rec else e_start-ei[count_owner]
                    omega_end = e_end+ei[count_owner] if rec else e_end-ei[count_owner]
                    first_energy_panel = np.floor(omega_start/max_energy_panel).astype(int)
                    energy_pieces = np.ceil(omega_end/max_energy_panel).astype(int)-first_energy_panel
                    segment = np.repeat(np.arange(len(count_owner)), energy_pieces)
                    first_offset = np.r_[0, np.cumsum(energy_pieces)[:-1]]
                    number = np.arange(len(segment))-np.repeat(first_offset, energy_pieces)
                    energy_left = np.maximum(omega_start[segment],
                                              (first_energy_panel[segment]+number)*max_energy_panel)
                    energy_right = np.minimum(omega_end[segment],
                                               (first_energy_panel[segment]+number+1)*max_energy_panel)
                    # Each count piece lies within one linear E(x) segment.
                    # Intersection with fixed Omega panels introduces no extrapolation.
                    density = (count_end[segment]-count_start[segment])/(e_end[segment]-e_start[segment])
                    starts = count_start[segment]+(energy_left-omega_start[segment])*density
                    span = (energy_right-energy_left)*density
                    owner = count_owner[segment]
                elif integration_layout == "energy_panels":
                    omega_lower = inner_low_energy[valid]+ei if rec else inner_low_energy[valid]-ei
                    omega_upper = inner_high_energy[valid]+ei if rec else inner_high_energy[valid]-ei
                    first_panel = np.floor(omega_lower/max_energy_panel).astype(int)
                    panel_count = np.ceil(omega_upper/max_energy_panel).astype(int)-first_panel
                else:
                    panel_count = np.ceil((high-low)/max_count_panel).astype(int)
                if integration_layout != "resolved_panels":
                    owner = np.repeat(np.arange(len(low)), panel_count)
                    first = np.r_[0, np.cumsum(panel_count)[:-1]]
                    number = np.arange(len(owner))-np.repeat(first, panel_count)
                if integration_layout == "energy_panels":
                    omega_start = np.maximum(omega_lower[owner], (first_panel[owner]+number)*max_energy_panel)
                    omega_end = np.minimum(omega_upper[owner], (first_panel[owner]+number+1)*max_energy_panel)
                    energy_start = omega_start-ei[owner] if rec else omega_start+ei[owner]
                    energy_end = omega_end-ei[owner] if rec else omega_end+ei[owner]
                    starts = np.interp(energy_start, energy, counts)
                    ends = np.interp(energy_end, energy, counts)
                    span = ends-starts
                elif integration_layout == "native_intervals":
                    widths = (high-low)/panel_count
                    starts = low[owner]+number*widths[owner]
                    span = widths[owner]
                inner = starts[:, None]+span[:, None]*(node[None, :]+1)/2
                inner_weight = span[:, None]*weight[None, :]/2
                e_i = np.broadcast_to(ei[owner, None], inner.shape).ravel()
                r_i = np.broadcast_to(np.interp(outer, counts, self.coherence_ratio)[owner, None], inner.shape).ravel()
                pair_weight = (outer_weight[owner, None]*inner_weight).ravel()
            else:
                inner = low[:, None]+(high-low)[:, None]*(node[None, :]+1)/2
                inner_weight = (high-low)[:, None]*weight[None, :]/2
                e_i = np.broadcast_to(ei[:, None], inner.shape).ravel()
                r_i = np.broadcast_to(np.interp(outer, counts, self.coherence_ratio)[:, None], inner.shape).ravel()
                pair_weight = (outer_weight[:, None]*inner_weight).ravel()
            e_j = np.interp(inner.ravel(), counts, energy)
            r_j = np.interp(inner.ravel(), counts, self.coherence_ratio)
            omega = e_i+e_j if rec else e_j-e_i
            coefficient = (rate_prefactor*pair_weight
                           *(1+r_i*r_j if rec else 1-r_i*r_j)*phonons.coupling(omega))
            channels.append((e_i, e_j, omega, coefficient, np.full(len(omega), rec, dtype=bool)))
        if channels:
            e_i, e_j, omega, coefficient, rec = [np.concatenate([channel[i] for channel in channels]) for i in range(5)]
        else:
            e_i = e_j = omega = coefficient = np.array([], dtype=float)
            rec = np.array([], dtype=bool)
        if ultraviolet_partition is not None:
            lower, upper = ultraviolet_partition
            if not (np.isfinite(lower) and np.isfinite(upper) and 0 <= lower < upper <= omega_high):
                raise ValueError("UV partition must be an increasing interval within the coupling support")

            def fraction(values):
                result = np.zeros_like(values)
                result[values >= upper] = 1.
                interior = (values > lower) & (values < upper)
                t = (values[interior]-lower)/(upper-lower)
                # Evaluate the complementary half by symmetry to retain [0,1].
                left = t <= .5
                u = np.where(left, t, 1-t)
                smooth = u**3*(10+u*(-15+6*u))
                result[interior] = np.where(left, smooth, 1-smooth)
                return result

            native_coefficient = self.coefficients*(1-fraction(self.omega))
            coefficient = np.r_[coefficient*fraction(omega), native_coefficient]
            e_i = np.r_[e_i, self.electron_energies[self.index_i]]
            e_j = np.r_[e_j, self.electron_energies[self.index_j]]
            omega = np.r_[omega, self.omega]
            rec = np.r_[rec, self.recombination]
        if np.any(coefficient < 0) or np.any(~np.isfinite(coefficient)):
            raise ValueError("projected event coefficients must be finite and nonnegative")
        active = coefficient > 0
        e_i, e_j, omega, coefficient, rec = [array[active] for array in (e_i, e_j, omega, coefficient, rec)]
        if np.any(omega < phonons.energies[0]) or np.any(omega > phonons.energies[-1]):
            raise ValueError("projected coupled event lacks a phonon bracket")
        self.target_energy_i, self.target_energy_j = _readonly(e_i), _readonly(e_j)
        self.omega, self.coefficients, self.recombination = _readonly(omega), _readonly(coefficient), _readonly(rec)
        self.phonon_lower, self.phonon_upper, self.beta_lower, self.beta_upper = self._map_energy(omega, phonons.energies)
        self.electron_lower_i, self.electron_upper_i, self.electron_beta_lower_i, self.electron_beta_upper_i = self._map_energy(e_i, energy)
        self.electron_lower_j, self.electron_upper_j, self.electron_beta_lower_j, self.electron_beta_upper_j = self._map_energy(e_j, energy)
        original_admission = self.metadata
        self.metadata = {
            "schema": "pysnspd.projected_electron_phonon_events.v1",
            "native_pair_admission": original_admission,
            "spectrum_label": label, "phonon_label": phonons.label,
            "rate_prefactor": float(rate_prefactor), "rate_prefactor_convention": "8*pi*Delta0*t_ref/hbar",
            "electron_capacity": "4*w_i", "phonon_capacity": "(g_ph/N0)*d(Omega/Delta0)",
            "energy_unit": "N0*Delta0^2", "coupling_support": list(phonons.coupling_support),
            "quadrature_order": quadrature_order, "native_count_support": [float(xmin), float(xmax)],
            "integration_layout": integration_layout, "outer_order": outer_order,
            "max_count_panel": max_count_panel,
            "max_energy_panel": max_energy_panel,
            "active_scattering": int(np.count_nonzero(~rec)), "active_recombination": int(np.count_nonzero(rec)),
            "recombination_counting": "ordered continuous domain Ei<=Ej; the diagonal has zero integration measure",
            "electron_reconstruction": "piecewise-linear E(x) and R2/N1 on native count nodes; no extrapolation",
            "electron_stoichiometry": "two positive native-energy barycentric weights per virtual electron",
            "electron_activities": "unnormalized geometric Pauli activities using the same stoichiometric weights",
            "support_limitation": "The small intervals outside first/last native count nodes are excluded explicitly; tail/error diagnostics are mandatory",
            "ultraviolet_partition": list(ultraviolet_partition) if ultraviolet_partition is not None else None,
            "partition_policy": "alpha*(1-chi) native pair quadrature plus alpha*chi nested quadrature; chi is a C2 quintic smooth step" if ultraviolet_partition is not None else None,
        }
        if hasattr(self, "index_i"):
            del self.index_i, self.index_j  # native pair indices do not describe projected events

    @staticmethod
    def _map_energy(target, grid):
        if np.any(target < grid[0]) or np.any(target > grid[-1]):
            raise ValueError("virtual energy outside native support; no extrapolation")
        upper = np.searchsorted(grid, target, side="left")
        lower = np.maximum(0, upper-1)
        span = grid[upper]-grid[lower]
        beta = np.divide(target-grid[lower], span, out=np.zeros_like(target), where=span > 0)
        return _readonly(lower), _readonly(upper), _readonly(1-beta), _readonly(beta)

    def _electron_activity_log(self, nodal_logs, suffix):
        lower = getattr(self, "electron_lower_"+suffix)
        upper = getattr(self, "electron_upper_"+suffix)
        beta_lower = getattr(self, "electron_beta_lower_"+suffix)
        beta_upper = getattr(self, "electron_beta_upper_"+suffix)
        first, second = np.zeros_like(self.omega), np.zeros_like(self.omega)
        np.multiply(beta_lower, nodal_logs[lower], out=first, where=beta_lower > 0)
        np.multiply(beta_upper, nodal_logs[upper], out=second, where=beta_upper > 0)
        return first+second

    def log_activities(self, p, n):
        p, n = self.validate(p, n)
        logp, loghole = _log_nonnegative(p), _log_nonnegative(1-p)
        pi, pj = self._electron_activity_log(logp, "i"), self._electron_activity_log(logp, "j")
        hi, hj = self._electron_activity_log(loghole, "i"), self._electron_activity_log(loghole, "j")
        emission = self._barycentric_log(np.log1p(n))
        absorption = self._barycentric_log(_log_nonnegative(n))
        return (np.where(self.recombination, pi+pj, hi+pj)+emission,
                np.where(self.recombination, hi+hj, pi+hj)+absorption)

    def rhs_from_rates(self, rates):
        rates = np.asarray(rates, dtype=float)
        if rates.shape != self.omega.shape or np.any(~np.isfinite(rates)):
            raise ValueError("one finite signed rate is required per projected event")
        electronic = np.zeros_like(self.electron_energies)
        for suffix, signed in (("i", np.where(self.recombination, -rates, rates)), ("j", -rates)):
            electronic += np.bincount(getattr(self, "electron_lower_"+suffix),
                weights=getattr(self, "electron_beta_lower_"+suffix)*signed, minlength=len(electronic))
            electronic += np.bincount(getattr(self, "electron_upper_"+suffix),
                weights=getattr(self, "electron_beta_upper_"+suffix)*signed, minlength=len(electronic))
        phononic = np.bincount(self.phonon_lower, weights=self.beta_lower*rates,
                               minlength=len(self.phonons.energies))
        phononic += np.bincount(self.phonon_upper, weights=self.beta_upper*rates,
                                minlength=len(self.phonons.energies))
        dp, dn = electronic/self.electron_capacities, phononic/self.phonons.capacities
        if np.any(~np.isfinite(dp)) or np.any(~np.isfinite(dn)):
            raise FloatingPointError("projected occupation derivatives overflowed; no repair was applied")
        return dp, dn


class HybridElectronPhononEvents(ProjectedElectronPhononEvents):
    """Use native Gauss pairs in the bulk and nested quadrature at the UV cut.

    A smooth partition of unity changes the quadrature, not alpha2F. Retaining
    native pairs in the bulk avoids imposing a coarse reconstructed population
    on a narrow nonthermal band. Both positive pieces preserve the same event
    identities. The default UV interval is [0.75,0.95]*Omega_max and is recorded.
    """

    def __init__(self, electron_energies, electron_weights, count_nodes, coherence_ratio,
                 phonons: PhononGrid, *, rate_prefactor: float, quadrature_order: int,
                 label: str, ultraviolet_partition: tuple[float, float] | None = None,
                 integration_layout: str = "global", outer_order: int = 2,
                 max_count_panel: float = .5, max_energy_panel: float = .125):
        if ultraviolet_partition is None:
            cutoff = phonons.coupling_support[1]
            ultraviolet_partition = (.75*cutoff, .95*cutoff)
        super().__init__(electron_energies, electron_weights, count_nodes, coherence_ratio, phonons,
                         rate_prefactor=rate_prefactor, quadrature_order=quadrature_order,
                         label=label, ultraviolet_partition=ultraviolet_partition,
                         integration_layout=integration_layout, outer_order=outer_order, max_count_panel=max_count_panel,
                         max_energy_panel=max_energy_panel)
        self.metadata["schema"] = "pysnspd.hybrid_electron_phonon_events.v1"
