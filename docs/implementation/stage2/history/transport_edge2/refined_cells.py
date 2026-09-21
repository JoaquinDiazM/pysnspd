"""Complementary occupation meshes with a common causal electronic energy.

The archived R2 catalogue remains immutable. Its vacuum, units, field domain
and regulator are reused, while the count quadrature is explicitly replaced.
The excitation energy is obtained directly from the same causal state count;
its derivatives follow implicit differentiation at FIXED count, not separate
interpolation of a force. This opt-in numerical refinement is not production.
"""
from dataclasses import dataclass
import numpy as np
from .energy_catalog import (OccupationEnergyCatalog, energy_at_count_batch,
                             retarded_spectrum, segmented_count_quadrature)


@dataclass(frozen=True)
class ComplementaryCountCatalog:
    base: OccupationEnergyCatalog
    count_nodes: np.ndarray
    count_weights: np.ndarray

    def __post_init__(self):
        x, w = (np.array(a, dtype=float, copy=True) for a in
                (self.count_nodes, self.count_weights))
        if (x.ndim != 1 or len(x) < 2 or w.shape != x.shape
                or np.any(~np.isfinite(x)) or np.any(~np.isfinite(w))
                or np.any(x <= 0) or np.any(w <= 0) or np.any(np.diff(x) <= 0)):
            raise ValueError('positive ordered count nodes and positive finite weights required')
        for key, value in (('count_nodes', x), ('count_weights', w)):
            value.setflags(write=False)
            object.__setattr__(self, key, value)

    @property
    def vacuum(self):
        return self.base.vacuum

    @property
    def eta(self):
        return self.base.eta

    def energy_kernel(self, delta, gamma):
        self.vacuum.evaluate(delta, gamma)  # Preserve the admitted field domain.
        x = self.count_nodes
        if delta == 0:
            return x.copy(), np.zeros_like(x), np.zeros_like(x)
        energy = energy_at_count_batch(x, delta=delta, gamma=gamma, eta=self.eta)
        c, s = retarded_spectrum(energy, delta=delta, gamma=gamma, eta=self.eta)
        derivative_delta = s.imag/c.real
        derivative_gamma = -s.real*derivative_delta
        if gamma == 0:
            # Algebraic derivative avoids subtractive loss next to the BCS edge.
            derivative_delta = x*delta/np.sqrt((x*x+self.eta*self.eta)
                                               *(x*x+self.eta*self.eta+delta*delta))
        values = (energy, derivative_delta, derivative_gamma)
        if (np.any(~np.isfinite(values)) or np.any(energy <= 0)
                or np.any(np.diff(energy) <= 0)):
            raise FloatingPointError('direct causal energy lost finite monotone support')
        return values

    def evaluate(self, delta, gamma, occupation):
        p = np.asarray(occupation, float)
        if (p.shape != self.count_nodes.shape or np.any(~np.isfinite(p))
                or np.any((p < 0) | (p > 1))):
            raise ValueError('occupation must match complementary count nodes and lie in [0,1]')
        return tuple(float(u+4*np.dot(self.count_weights*k, p)) for u, k in
                     zip(self.vacuum.evaluate(delta, gamma), self.energy_kernel(delta, gamma)))


def refined_count_catalog(base, refinement=1, *, bulk_spacing=.02, order=2):
    """Explicit nested interval refinement; no interpolation of an old p vector.

    Smooth experimental initial distributions are sampled on each new mesh.
    Gauss nodes themselves are not nested. The near-edge count intervals resolve
    the finite regulator; bulk and tail intervals bound interpolation distances.
    """
    if type(refinement) is not int or refinement < 1:
        raise ValueError('refinement must be a positive integer')
    if not np.isfinite(bulk_spacing) or bulk_spacing <= 0:
        raise ValueError('bulk spacing must be positive and finite')
    edge = np.r_[0., np.geomspace(1e-8, .01, 13)]
    pieces = [edge]
    for low, high, step in ((.01, 3., bulk_spacing), (3., 6., .05), (6., 12., .15)):
        intervals = int(np.ceil((high-low)/step))*refinement
        pieces.append(np.linspace(low, high, intervals+1)[1:])
    breaks = np.concatenate(pieces)
    x, w = segmented_count_quadrature(breaks, order)
    return ComplementaryCountCatalog(base, x, w)
