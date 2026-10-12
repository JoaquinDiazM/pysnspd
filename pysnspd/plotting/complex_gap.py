"""Saved-data gap visualization: complex spatial P1, then magnitude.

The raster is a display sampling of the saved finite-element field. It is
neither a new evolution mesh nor a reconstruction of unsaved time dynamics.
"""
from __future__ import annotations

import numpy as np
from matplotlib.tri import LinearTriInterpolator


def complex_gap_at(triangulation, delta, x, y):
    """Evaluate |P1 Delta|, masking points outside the saved triangulation."""
    values = np.asarray(delta)
    if (values.shape != triangulation.x.shape or values.ndim != 1
            or not np.all(np.isfinite(values))):
        raise ValueError("One finite complex gap value per saved node required")
    xx, yy = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float))
    if not np.all(np.isfinite(xx)) or not np.all(np.isfinite(yy)):
        raise ValueError("Finite physical sampling coordinates required")
    real = LinearTriInterpolator(triangulation, values.real)(xx, yy)
    imag = LinearTriInterpolator(triangulation, values.imag)(xx, yy)
    return np.ma.abs(real + 1j * imag)


def complex_gap_raster(triangulation, delta, *, bounds=None, samples=601):
    """Sample |P1 Delta| on a physical rectangle; at most samples per axis.

    Coordinates inherit the caller's length units. Raster values are evaluated
    before magnitude; no interpolation of nodal magnitudes is performed.
    """
    if isinstance(samples, bool) or int(samples) != samples or samples < 3:
        raise ValueError("An integer raster resolution of at least three required")
    if bounds is None:
        bounds = (triangulation.x.min(), triangulation.x.max(),
                  triangulation.y.min(), triangulation.y.max())
    limits = np.asarray(bounds, float)
    if (limits.shape != (4,) or not np.all(np.isfinite(limits))
            or limits[1] <= limits[0] or limits[3] <= limits[2]):
        raise ValueError("Finite ordered physical raster bounds required")
    lengths = limits[[1, 3]] - limits[[0, 2]]
    counts = np.maximum(2, np.ceil((int(samples) - 1) * lengths / lengths.max()).astype(int) + 1)
    x = np.linspace(limits[0], limits[1], counts[0])
    y = np.linspace(limits[2], limits[3], counts[1])
    xx, yy = np.meshgrid(x, y)
    return x, y, complex_gap_at(triangulation, delta, xx, yy)


def draw_complex_gap(ax, triangulation, delta, *, unit_J, vmin, vmax,
                     cmap="viridis", bounds=None, samples=601):
    """Render sampled |P1 Delta| with nearest cells, never Gouraud |Delta_node|."""
    unit = float(unit_J)
    if not np.isfinite(unit) or unit <= 0:
        raise ValueError("A finite positive energy unit required")
    x, y, gap = complex_gap_raster(triangulation, delta, bounds=bounds, samples=samples)
    return ax.pcolormesh(x, y, gap / unit, shading="nearest", cmap=cmap,
                        vmin=vmin, vmax=vmax, rasterized=True)
