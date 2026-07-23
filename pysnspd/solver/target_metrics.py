"""Stationary-state target diagnostics for the SS gTDGL smoke runs.

The routines in this module are deliberately lightweight: they do not decide
new physics, they only build a smooth metallic-contact seed and quantify the
three checks needed before moving to photon dynamics:

1. bulk gauge-fixed physical stationarity of the phase gradient and
   electric-potential gradient, excluding the normal-contact conversion region,
2. a contact-healing length of order a few physical coherence lengths,
3. finite-volume current continuity.

For the present no-screening A=0 backend the physical phase-gradient diagnostic
is the edge superfluid momentum ``Q = grad(arg(Delta))`` stored by the current
adapter.  This is invariant under constant phase shifts.  The electrostatic
diagnostic is the edge gradient of ``phi``; it is invariant under constant
potential offsets.  Metallic contacts are allowed to have conversion fields, so
the stationarity gate is evaluated on a bulk-edge mask away from the contacts.  In a fully electromagnetic gauge treatment these would be
replaced by the gauge-covariant combinations involving A.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from pysnspd.solver.diagnostics import current_residual, max_current_residual
from pysnspd.gtdgl.material import GTDGLMaterial

def _snapshot_2d(value: Any) -> np.ndarray:
    arr = np.asarray(value, dtype=float)
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    return arr


def _edge_phi_gradient_snapshots(history: dict[str, np.ndarray]) -> np.ndarray:
    direct = history.get("edge_phi_gradient_snapshot_V_m")
    if direct is not None:
        arr = _snapshot_2d(direct)
        if arr.shape[0] >= 2:
            return arr

    phi = _snapshot_2d(history.get("phi_snapshot_V", []))
    edge_i = np.asarray(history.get("edge_i", []), dtype=np.int64).reshape(-1)
    edge_j = np.asarray(history.get("edge_j", []), dtype=np.int64).reshape(-1)
    edge_length = np.asarray(history.get("edge_length_m", []), dtype=float).reshape(-1)
    if phi.shape[0] < 2 or edge_i.size == 0 or edge_j.size != edge_i.size or edge_length.size != edge_i.size:
        return np.empty((0, 0), dtype=float)
    if int(np.max(edge_i, initial=-1)) >= phi.shape[1] or int(np.max(edge_j, initial=-1)) >= phi.shape[1]:
        return np.empty((0, 0), dtype=float)
    length = np.maximum(edge_length, 1.0e-300)
    return (phi[:, edge_j] - phi[:, edge_i]) / length[None, :]


def _active_phase_edges(
    history: dict[str, np.ndarray],
    *,
    edge_active_threshold: float,
    bulk_exclusion_xi: float = 0.0,
) -> np.ndarray | None:
    explicit = history.get("stationarity_active_edge_mask")
    if explicit is not None:
        mask = np.asarray(explicit, dtype=bool).reshape(-1)
    else:
        psi_r = _snapshot_2d(history.get("psi_snapshot_real_J", []))
        psi_i = _snapshot_2d(history.get("psi_snapshot_imag_J", []))
        edge_i = np.asarray(history.get("edge_i", []), dtype=np.int64).reshape(-1)
        edge_j = np.asarray(history.get("edge_j", []), dtype=np.int64).reshape(-1)
        if psi_r.shape[0] < 1 or psi_i.shape != psi_r.shape or edge_i.size == 0 or edge_j.size != edge_i.size:
            return None
        if int(np.max(edge_i, initial=-1)) >= psi_r.shape[1] or int(np.max(edge_j, initial=-1)) >= psi_r.shape[1]:
            return None
        psi = psi_r[-1] + 1j * psi_i[-1]
        amp_edge = 0.5 * (np.abs(psi[edge_i]) + np.abs(psi[edge_j]))
        finite_amp = amp_edge[np.isfinite(amp_edge)]
        if finite_amp.size == 0:
            return None
        bulk = float(np.nanpercentile(finite_amp, 90.0))
        threshold = float(np.clip(edge_active_threshold, 0.0, 0.95)) * max(bulk, 1.0e-300)
        mask = amp_edge >= threshold

    terminal_edges = history.get("normal_terminal_edge_mask")
    if terminal_edges is not None:
        term = np.asarray(terminal_edges, dtype=bool).reshape(-1)
        if term.size == mask.size:
            mask = mask & ~term

    d_contact = np.asarray(history.get("edge_distance_from_contact_m", []), dtype=float).reshape(-1)
    xi_hist = np.asarray(history.get("stationarity_xi_m", []), dtype=float).reshape(-1)
    if (
        d_contact.size == mask.size
        and xi_hist.size
        and np.isfinite(xi_hist[0])
        and float(bulk_exclusion_xi) > 0.0
    ):
        min_distance = float(bulk_exclusion_xi) * max(float(xi_hist[0]), 1.0e-300)
        mask = mask & np.isfinite(d_contact) & (d_contact >= min_distance)

    if not np.any(mask):
        return None
    return mask


def _edge_field_change_metrics(field: np.ndarray, mask: np.ndarray | None, *, abs_tol: float) -> dict[str, float | int]:
    arr = _snapshot_2d(field)
    total = int(arr.shape[1]) if arr.ndim == 2 else 0
    if arr.ndim != 2 or arr.shape[0] < 2 or total == 0:
        return {
            "rel_change": float("nan"),
            "abs_change": float("nan"),
            "rms_final": float("nan"),
            "active_fraction": 0.0,
            "active_count": 0,
            "total_count": total,
        }
    if mask is None:
        active = np.ones(total, dtype=bool)
    else:
        active = np.asarray(mask, dtype=bool).reshape(-1)
        if active.size != total:
            active = np.ones(total, dtype=bool)
    finite = np.isfinite(arr[-1]) & np.isfinite(arr[-2]) & active
    if np.count_nonzero(finite) == 0:
        return {
            "rel_change": float("nan"),
            "abs_change": float("nan"),
            "rms_final": float("nan"),
            "active_fraction": 0.0,
            "active_count": 0,
            "total_count": total,
        }
    final = arr[-1, finite]
    prev = arr[-2, finite]
    diff = final - prev
    rms_final = float(np.sqrt(np.nanmean(final * final)))
    rms_diff = float(np.sqrt(np.nanmean(diff * diff)))
    scale = max(rms_final, float(abs_tol), 1.0e-300)
    return {
        "rel_change": rms_diff / scale,
        "abs_change": rms_diff,
        "rms_final": rms_final,
        "active_fraction": float(np.count_nonzero(finite) / max(total, 1)),
        "active_count": int(np.count_nonzero(finite)),
        "total_count": total,
    }


def _transverse_delta_profiles(
    *,
    x_m: np.ndarray,
    delta_over_delta0: np.ndarray,
    xi_m: float,
    bulk_exclusion_xi: float,
) -> tuple[np.ndarray, float]:
    x = np.asarray(x_m, dtype=float).reshape(-1)
    values = np.asarray(delta_over_delta0, dtype=float)
    if values.ndim != 2 or values.shape[1] != x.size or x.size == 0:
        return np.empty((0, 0), dtype=float), float("nan")
    xmin = float(np.nanmin(x))
    xmax = float(np.nanmax(x))
    length = max(xmax - xmin, 1.0e-300)
    exclusion_m = (
        max(0.0, float(bulk_exclusion_xi)) * float(xi_m)
        if np.isfinite(xi_m) and xi_m > 0.0
        else 0.1 * length
    )
    lo = xmin + min(exclusion_m, 0.4 * length)
    hi = xmax - min(exclusion_m, 0.4 * length)
    bulk_nodes = np.isfinite(x) & (x >= lo) & (x <= hi)
    if np.count_nonzero(bulk_nodes) < 16:
        lo = xmin + 0.1 * length
        hi = xmax - 0.1 * length
        bulk_nodes = np.isfinite(x) & (x >= lo) & (x <= hi)
    bulk_length = max(hi - lo, 1.0e-300)
    if np.isfinite(xi_m) and xi_m > 0.0:
        n_bins = int(np.clip(np.ceil(2.0 * bulk_length / xi_m), 16, 128))
    else:
        n_bins = 64
    edges = np.linspace(lo, hi, n_bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    bin_index = np.clip(np.digitize(x, edges) - 1, 0, n_bins - 1)
    profiles = np.full((values.shape[0], n_bins), np.nan, dtype=float)
    for bin_number in range(n_bins):
        mask = bulk_nodes & (bin_index == bin_number)
        if np.any(mask):
            profiles[:, bin_number] = np.nanmedian(values[:, mask], axis=1)
    for row_index in range(profiles.shape[0]):
        finite = np.isfinite(profiles[row_index])
        if np.count_nonzero(finite) >= 2:
            profiles[row_index, ~finite] = np.interp(
                centers[~finite],
                centers[finite],
                profiles[row_index, finite],
            )
    valid_columns = np.all(np.isfinite(profiles), axis=0)
    profiles = profiles[:, valid_columns]
    bin_width_m = bulk_length / max(n_bins, 1)
    return profiles, float(bin_width_m)


def _count_true_runs(mask: np.ndarray, *, minimum_bins: int) -> int:
    active = np.asarray(mask, dtype=bool).reshape(-1)
    if active.size == 0:
        return 0
    padded = np.concatenate(([False], active, [False])).astype(np.int8)
    transitions = np.diff(padded)
    starts = np.flatnonzero(transitions == 1)
    stops = np.flatnonzero(transitions == -1)
    widths = stops - starts
    return int(np.count_nonzero(widths >= max(1, int(minimum_bins))))


def _tail_scalar_envelope_metrics(values: np.ndarray, *, absolute_scale: float) -> tuple[float, float]:
    array = np.asarray(values, dtype=float).reshape(-1)
    array = array[np.isfinite(array)]
    if array.size < 4:
        return float("nan"), float("nan")
    median = float(np.nanmedian(array))
    scale = max(abs(median), abs(float(absolute_scale)), 1.0e-300)
    span = float((np.nanpercentile(array, 95.0) - np.nanpercentile(array, 5.0)) / scale)
    midpoint = max(1, array.size // 2)
    first = float(np.nanmedian(array[:midpoint]))
    second = float(np.nanmedian(array[midpoint:]))
    drift = abs(second - first) / scale
    return span, float(drift)


def _first_binned_crossing_distance(
    distance: np.ndarray,
    values: np.ndarray,
    threshold: float,
    *,
    bin_width_m: float,
    max_distance: float,
) -> float:
    d = np.asarray(distance, dtype=float).reshape(-1)
    v = np.asarray(values, dtype=float).reshape(-1)
    finite = np.isfinite(d) & np.isfinite(v) & (d >= 0.0) & (d <= max_distance)
    if np.count_nonzero(finite) == 0:
        return float("nan")
    d = d[finite]
    v = v[finite]
    order = np.argsort(d)
    d = d[order]
    v = v[order]
    nbins = max(2, int(np.ceil(max_distance / max(bin_width_m, 1.0e-300))))
    edges = np.linspace(0.0, max_distance, nbins + 1)
    for lo, hi in zip(edges[:-1], edges[1:]):
        mask = (d >= lo) & (d < hi if hi < max_distance else d <= hi)
        if np.count_nonzero(mask) < 1:
            continue
        if float(np.nanmedian(v[mask])) >= threshold:
            return float(0.5 * (lo + hi))
    return float("nan")
