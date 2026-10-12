"""Physical gap magnitude in saved time; candidate analysis only.

Reconstruct complex P1 in space at EACH saved frame, take its magnitude,
then interpolate that scalar observable in time. No evolution or gate.
"""
from __future__ import annotations
import numpy as np


def spatial_magnitude_at(snap, sampler, time_s):
    """Return (1-f)|P1 Delta_j|+f|P1 Delta_{j+1}| on a fixed ROI.

    Invariant under arbitrary global phase rotations of each saved frame.
    This is a saved-data interpolation, not a reconstruction of unsaved
    dynamics. It is never used for phase winding or electrical events.
    """
    times = np.asarray(snap["time_s"], dtype=float)
    delta = np.asarray(snap["delta_J"])
    nodes = np.asarray(sampler["nodes"], dtype=int)
    bary = np.asarray(sampler["barycentric"], dtype=float)
    t = float(time_s)
    if (times.ndim != 1 or len(times) < 1 or not np.all(np.isfinite(times))
            or np.any(np.diff(times) <= 0)
            or delta.ndim != 2 or len(delta) != len(times)
            or not np.all(np.isfinite(delta)) or not np.isfinite(t)
            or nodes.ndim != 2 or nodes.shape != bary.shape
            or not np.all(np.isfinite(bary))
            or np.any(nodes < 0) or np.any(nodes >= delta.shape[1])):
        raise ValueError("Finite, ordered saved frames and a valid P1 sampler required")
    if t < times[0] or t > times[-1]:
        raise ValueError("No temporal extrapolation is permitted")

    def magnitude(index):
        return np.abs(np.sum(delta[index][nodes] * bary, axis=1))

    j = int(np.searchsorted(times, t))
    if j == 0:
        return magnitude(0)
    if j == len(times):
        return magnitude(-1)
    if t == times[j]:
        return magnitude(j)
    fraction = (t - times[j - 1]) / (times[j] - times[j - 1])
    return ((1 - fraction) * magnitude(j - 1)
            + fraction * magnitude(j))


def field_L2_saved_magnitude(reference, candidate, reference_sampler,
                             candidate_sampler, grid_s, dc_gap_scale_J):
    """Magnitude-field L2/DC at an explicit shared analysis grid; no gate."""
    scale = float(dc_gap_scale_J)
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("A finite positive DC gap scale is required")
    values = []
    for t in np.asarray(grid_s, dtype=float):
        a = spatial_magnitude_at(reference, reference_sampler, t)
        b = spatial_magnitude_at(candidate, candidate_sampler, t)
        if a.shape != b.shape:
            raise ValueError("Both samplers must use the same physical ROI points")
        values.append(float(np.sqrt(np.mean((a - b) ** 2)) / scale))
    return np.asarray(values)
