"""Fixed-energy transport events on the native R2 occupation quadrature.

This experimental weak discretization shares electron number AND native energy
between two cells. At a face energy, the participating occupation is represented
by its log odds on two adjacent native states. The same barycentric coefficients
distribute the event. This preserves Fermi equilibrium and the physical boundary
faces, but needs convergence tests against D.17; it is not a new physical kernel.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from .cell_validation import ElectronicCell, occupation_array
from .energy_catalog import retarded_spectrum


def _geometric_log(values, indices, weights):
    selected = values[indices]
    terms = np.zeros_like(weights)
    active = weights > 0
    with np.errstate(divide="ignore"):
        terms[active] = weights[active]*np.log(selected[active])
    return np.sum(terms, axis=1)


def _difference_of_exponentials(forward, backward):
    """Signed net rate without cancellation; exact zero factors remain zero."""
    result = np.zeros_like(forward)
    positive = np.isfinite(forward) & (forward >= backward)
    negative = np.isfinite(backward) & (backward > forward)
    result[positive] = np.exp(forward[positive])*(-np.expm1(backward[positive]-forward[positive]))
    result[negative] = np.exp(backward[negative])*np.expm1(forward[negative]-backward[negative])
    return result


def _brackets(energies, query):
    upper = np.searchsorted(energies, query, side="right")
    # Query points are interior Gauss points of the shared support.
    if np.any(upper < 1) or np.any(upper >= len(energies)):
        raise ValueError("face energy has no two native bracketing states")
    lower = upper-1
    right = (query-energies[lower])/(energies[upper]-energies[lower])
    weights = np.stack((1-right, right), axis=1)
    if np.any(weights < 0) or np.any(weights > 1):
        raise ValueError("transport would extrapolate a population")
    return np.stack((lower, upper), axis=1), weights


@dataclass(frozen=True)
class NativeTransportEvents:
    left: ElectronicCell
    right: ElectronicCell
    diffusion_over_length_squared: float = .03
    gauss_order: int = 2

    def __post_init__(self):
        coefficient = float(self.diffusion_over_length_squared)
        if not np.isfinite(coefficient) or coefficient < 0:
            raise ValueError("diffusion must be finite and nonnegative")
        if not isinstance(self.gauss_order, (int, np.integer)) or self.gauss_order < 1:
            raise ValueError("gauss_order must be a positive integer")
        if self.left.catalog.eta != self.right.catalog.eta:
            raise ValueError("both cells need the same causal regulator")
        low = max(self.left.energies[0], self.right.energies[0])
        high = min(self.left.energies[-1], self.right.energies[-1])
        if high <= low:
            raise ValueError("the electronic supports do not overlap")
        native = np.r_[self.left.energies, self.right.energies]
        breaks = np.unique(np.r_[low, native[(native > low) & (native < high)], high])
        nodes, weights = np.polynomial.legendre.leggauss(self.gauss_order)
        widths = np.diff(breaks)
        energy = (breaks[:-1, None]+widths[:, None]*(nodes+1)/2).ravel()
        measures = (widths[:, None]*weights/2).ravel()
        indices, barycentric, longitudinal = [], [], []
        for cell in (self.left, self.right):
            ix, b = _brackets(cell.energies, energy)
            c, s = retarded_spectrum(energy, delta=cell.amplitude,
                                     gamma=cell.gamma, eta=cell.catalog.eta)
            dl = c.real*c.real-s.imag*s.imag
            if np.any(~np.isfinite(dl)) or np.any(dl < 0):
                raise FloatingPointError("nonphysical longitudinal spectral diffusion")
            indices.append(ix); barycentric.append(b); longitudinal.append(dl)
        dl = np.asarray(longitudinal)
        denominator = np.sum(dl, axis=0)
        harmonic = np.divide(2*dl[0]*dl[1], denominator,
                             out=np.zeros_like(denominator), where=denominator > 0)
        conductance = 4*coefficient*measures*harmonic
        for name, value in (("energies", energy), ("measures", measures),
                            ("indices", np.asarray(indices)),
                            ("barycentric", np.asarray(barycentric)),
                            ("conductance", conductance)):
            value = np.array(value, copy=True)
            value.setflags(write=False)
            object.__setattr__(self, name, value)
        object.__setattr__(self, "shared_support", (float(low), float(high)))

    def rates(self, populations):
        left = self.left.validate(populations[0])
        right = self.right.validate(populations[1])
        la = _geometric_log(left, self.indices[0], self.barycentric[0])
        lb = _geometric_log(1-left, self.indices[0], self.barycentric[0])
        ra = _geometric_log(right, self.indices[1], self.barycentric[1])
        rb = _geometric_log(1-right, self.indices[1], self.barycentric[1])
        return self.conductance*_difference_of_exponentials(la+rb, ra+lb)

    def rhs(self, populations):
        flux = self.rates(populations)
        result = []
        for side, cell in enumerate((self.left, self.right)):
            count_rate = np.bincount(self.indices[side].ravel(),
                                    weights=(self.barycentric[side]*flux[:, None]).ravel(),
                                    minlength=len(cell.energies))
            result.append((2*side-1)*count_rate/(4*cell.weights))
        return np.asarray(result)

    def transferred_power(self, populations):
        return float(np.dot(self.energies, self.rates(populations)))

    def event_energy_residual(self):
        represented = [np.sum(cell.energies[self.indices[i]]*self.barycentric[i], axis=1)
                       for i, cell in enumerate((self.left, self.right))]
        return represented[1]-represented[0]
