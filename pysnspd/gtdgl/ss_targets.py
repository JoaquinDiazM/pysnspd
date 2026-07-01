"""Stationary-state target diagnostics for the SS gTDGL smoke runs.

The routines in this module are deliberately lightweight: they do not decide
new physics, they only build a smooth metallic-contact seed and quantify the
three checks needed before moving to photon dynamics:

1. gauge-fixed physical stationarity of the phase gradient and electric-potential
   gradient,
2. a contact-healing length of order a few physical coherence lengths,
3. finite-volume current continuity.

For the present no-screening A=0 backend the physical phase-gradient diagnostic
is the edge superfluid momentum ``Q = grad(arg(Delta))`` stored by the current
adapter.  This is invariant under constant phase shifts.  The electrostatic
diagnostic is the edge gradient of ``phi``; it is invariant under constant
potential offsets.  In a fully electromagnetic gauge treatment these would be
replaced by the gauge-covariant combinations involving A.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from pysnspd.gtdgl.diagnostics import current_residual, max_current_residual
from pysnspd.gtdgl.material import GTDGLMaterial


@dataclass(frozen=True)
class ProximitySeedDiagnostics:
    """Summary of the smooth metallic-terminal seed envelope."""

    enabled: bool
    healing_target_xi: float
    target_bulk_fraction: float
    xi_m: float
    envelope_length_m: float
    target_recovery_length_m: float
    min_envelope: float
    median_envelope: float
    max_envelope: float

    def as_dict(self) -> dict[str, Any]:
        return {
            "enabled": bool(self.enabled),
            "healing_target_xi": float(self.healing_target_xi),
            "target_bulk_fraction": float(self.target_bulk_fraction),
            "xi_m": float(self.xi_m),
            "envelope_length_m": float(self.envelope_length_m),
            "target_recovery_length_m": float(self.target_recovery_length_m),
            "min_envelope": float(self.min_envelope),
            "median_envelope": float(self.median_envelope),
            "max_envelope": float(self.max_envelope),
        }


@dataclass(frozen=True)
class ContactRecoveryDiagnostics:
    """Final |Delta| recovery length measured from left/right contacts."""

    bulk_delta_over_delta0: float
    threshold_fraction: float
    xi_m: float
    left_recovery_length_m: float
    right_recovery_length_m: float
    left_recovery_length_xi: float
    right_recovery_length_xi: float
    mean_recovery_length_xi: float
    min_allowed_xi: float
    max_allowed_xi: float
    passes: bool
    reason: str

    def as_dict(self) -> dict[str, Any]:
        return {
            "bulk_delta_over_delta0": float(self.bulk_delta_over_delta0),
            "threshold_fraction": float(self.threshold_fraction),
            "xi_m": float(self.xi_m),
            "left_recovery_length_m": float(self.left_recovery_length_m),
            "right_recovery_length_m": float(self.right_recovery_length_m),
            "left_recovery_length_xi": float(self.left_recovery_length_xi),
            "right_recovery_length_xi": float(self.right_recovery_length_xi),
            "mean_recovery_length_xi": float(self.mean_recovery_length_xi),
            "min_allowed_xi": float(self.min_allowed_xi),
            "max_allowed_xi": float(self.max_allowed_xi),
            "passes": bool(self.passes),
            "reason": str(self.reason),
        }


@dataclass(frozen=True)
class StationarityDiagnostics:
    """Gauge-fixed stationarity of edge phase-gradient and phi-gradient fields."""

    phase_gradient_rel_change: float
    phi_gradient_rel_change: float
    phase_gradient_abs_change_m_inv: float
    phi_gradient_abs_change_V_m: float
    phase_gradient_rms_final_m_inv: float
    phi_gradient_rms_final_V_m: float
    active_edge_fraction: float
    active_edge_count: int
    total_edge_count: int
    eta_R_final: float
    eta_R_window_max: float
    tolerance_phase_gradient_rel: float
    tolerance_phi_gradient_rel: float
    tolerance_phase_gradient_abs_m_inv: float
    tolerance_phi_gradient_abs_V_m: float
    edge_active_threshold_over_bulk: float
    passes: bool
    reason: str

    def as_dict(self) -> dict[str, Any]:
        return {
            "diagnostic": "gauge_fixed_edge_gradients_A_eq_0_v1",
            "phase_gradient_rel_change": float(self.phase_gradient_rel_change),
            "phi_gradient_rel_change": float(self.phi_gradient_rel_change),
            "phase_gradient_abs_change_m_inv": float(self.phase_gradient_abs_change_m_inv),
            "phi_gradient_abs_change_V_m": float(self.phi_gradient_abs_change_V_m),
            "phase_gradient_rms_final_m_inv": float(self.phase_gradient_rms_final_m_inv),
            "phi_gradient_rms_final_V_m": float(self.phi_gradient_rms_final_V_m),
            "active_edge_fraction": float(self.active_edge_fraction),
            "active_edge_count": int(self.active_edge_count),
            "total_edge_count": int(self.total_edge_count),
            "eta_R_final_info_only": float(self.eta_R_final),
            "eta_R_window_max_info_only": float(self.eta_R_window_max),
            "tolerance_phase_gradient_rel": float(self.tolerance_phase_gradient_rel),
            "tolerance_phi_gradient_rel": float(self.tolerance_phi_gradient_rel),
            "tolerance_phase_gradient_abs_m_inv": float(self.tolerance_phase_gradient_abs_m_inv),
            "tolerance_phi_gradient_abs_V_m": float(self.tolerance_phi_gradient_abs_V_m),
            "edge_active_threshold_over_bulk": float(self.edge_active_threshold_over_bulk),
            "passes": bool(self.passes),
            "reason": str(self.reason),
        }


@dataclass(frozen=True)
class ContinuityDiagnostics:
    """Current-continuity metrics in the final stationary state."""

    rms_residual_rel: float
    max_residual_rel: float
    native_poisson_residual_rel: float
    tolerance_rms: float
    tolerance_max: float
    tolerance_poisson: float
    passes: bool
    reason: str

    def as_dict(self) -> dict[str, Any]:
        return {
            "rms_residual_rel": float(self.rms_residual_rel),
            "max_residual_rel": float(self.max_residual_rel),
            "native_poisson_residual_rel": float(self.native_poisson_residual_rel),
            "tolerance_rms": float(self.tolerance_rms),
            "tolerance_max": float(self.tolerance_max),
            "tolerance_poisson": float(self.tolerance_poisson),
            "passes": bool(self.passes),
            "reason": str(self.reason),
        }


def physical_xi_m(material: GTDGLMaterial, Te_K: np.ndarray | float) -> float:
    """Return the median Appendix-B coherence length in meters."""
    xi2 = np.asarray(material.xi_mod_squared_m2(Te_K), dtype=float)
    xi = float(np.sqrt(np.nanmedian(np.maximum(xi2, 1.0e-300))))
    if not np.isfinite(xi) or xi <= 0.0:
        raise ValueError(f"Invalid physical coherence length xi={xi!r}.")
    return xi


def apply_terminal_proximity_seed(
    psi_dimensionless: np.ndarray,
    *,
    nodes_m: np.ndarray,
    material: GTDGLMaterial,
    Te_K: np.ndarray,
    healing_target_xi: float | None,
    target_bulk_fraction: float = 0.95,
    terminal_value: complex | float = 0.0,
) -> tuple[np.ndarray, ProximitySeedDiagnostics]:
    """Apply a smooth |Delta| ramp near metallic terminals to the initial seed.

    The ramp is only an initial condition; the actual boundary condition remains
    the solver's terminal Dirichlet value.  The profile is chosen so that a
    tanh envelope reaches ``target_bulk_fraction`` at approximately
    ``healing_target_xi`` physical coherence lengths from either terminal.
    """

    psi = np.asarray(psi_dimensionless, dtype=np.complex128).copy()
    nodes = np.asarray(nodes_m, dtype=float)
    Te = np.asarray(Te_K, dtype=float)
    if healing_target_xi is None or float(healing_target_xi) <= 0.0:
        diag = ProximitySeedDiagnostics(
            enabled=False,
            healing_target_xi=0.0,
            target_bulk_fraction=float(target_bulk_fraction),
            xi_m=physical_xi_m(material, Te),
            envelope_length_m=0.0,
            target_recovery_length_m=0.0,
            min_envelope=1.0,
            median_envelope=1.0,
            max_envelope=1.0,
        )
        return psi, diag

    frac = float(np.clip(target_bulk_fraction, 0.50, 0.999999))
    xi = physical_xi_m(material, Te)
    target_length = float(healing_target_xi) * xi
    envelope_length = target_length / max(float(np.arctanh(frac)), 1.0e-300)
    xmin = float(np.nanmin(nodes[:, 0]))
    xmax = float(np.nanmax(nodes[:, 0]))
    d_left = np.maximum(nodes[:, 0] - xmin, 0.0)
    d_right = np.maximum(xmax - nodes[:, 0], 0.0)
    d_contact = np.minimum(d_left, d_right)
    envelope = np.tanh(d_contact / max(envelope_length, 1.0e-300))
    envelope = np.clip(envelope, 0.0, 1.0)

    # Preserve phase while imposing a smooth amplitude ramp.  Boundary nodes are
    # still clamped exactly by TDGLSolver.apply_terminal_psi().
    phase = np.exp(1j * np.angle(psi))
    amp = np.abs(psi)
    term = complex(terminal_value)
    if abs(term) > 0.0:
        psi = phase * (abs(term) + envelope * np.maximum(amp - abs(term), 0.0))
    else:
        psi = phase * amp * envelope

    diag = ProximitySeedDiagnostics(
        enabled=True,
        healing_target_xi=float(healing_target_xi),
        target_bulk_fraction=frac,
        xi_m=xi,
        envelope_length_m=float(envelope_length),
        target_recovery_length_m=float(target_length),
        min_envelope=float(np.nanmin(envelope)),
        median_envelope=float(np.nanmedian(envelope)),
        max_envelope=float(np.nanmax(envelope)),
    )
    return psi, diag


def contact_recovery_diagnostics(
    *,
    psi_dimensionless: np.ndarray,
    nodes_m: np.ndarray,
    material: GTDGLMaterial,
    Te_K: np.ndarray,
    threshold_fraction: float = 0.95,
    min_allowed_xi: float = 1.5,
    max_allowed_xi: float = 4.0,
    bin_width_m: float | None = None,
) -> ContactRecoveryDiagnostics:
    """Measure where |Delta| reaches a fraction of the bulk value from contacts."""

    psi = np.asarray(psi_dimensionless, dtype=np.complex128).reshape(-1)
    nodes = np.asarray(nodes_m, dtype=float)
    amp = np.abs(psi)
    xi = physical_xi_m(material, Te_K)
    xmin = float(np.nanmin(nodes[:, 0]))
    xmax = float(np.nanmax(nodes[:, 0]))
    length = max(xmax - xmin, 1.0e-300)
    d_left = np.maximum(nodes[:, 0] - xmin, 0.0)
    d_right = np.maximum(xmax - nodes[:, 0], 0.0)
    d_contact = np.minimum(d_left, d_right)

    # Estimate bulk from the central region.  Fall back to the top decile if the
    # geometry is too short for the requested exclusion band.
    bulk_mask = d_contact >= min(0.35 * length, 4.0 * xi)
    if np.count_nonzero(bulk_mask) < max(10, amp.size // 20):
        cutoff = np.nanpercentile(amp, 90.0)
        bulk_vals = amp[amp >= cutoff]
    else:
        bulk_vals = amp[bulk_mask]
    bulk = float(np.nanmedian(bulk_vals)) if np.size(bulk_vals) else float(np.nanmedian(amp))
    bulk = max(bulk, 1.0e-300)
    threshold = float(np.clip(threshold_fraction, 0.50, 0.999999)) * bulk
    bw = float(bin_width_m) if bin_width_m is not None and bin_width_m > 0.0 else max(0.25 * xi, 1.0e-12)

    left_len = _first_binned_crossing_distance(d_left, amp, threshold, bin_width_m=bw, max_distance=0.5 * length)
    right_len = _first_binned_crossing_distance(d_right, amp, threshold, bin_width_m=bw, max_distance=0.5 * length)
    left_xi = left_len / xi
    right_xi = right_len / xi
    mean_xi = 0.5 * (left_xi + right_xi)
    finite = np.isfinite(left_xi) and np.isfinite(right_xi)
    passes = bool(finite and min_allowed_xi <= left_xi <= max_allowed_xi and min_allowed_xi <= right_xi <= max_allowed_xi)
    if passes:
        reason = "contact recovery length is within the requested xi window"
    elif not finite:
        reason = "failed to find a contact recovery crossing"
    else:
        reason = "contact recovery length is outside the requested xi window"

    return ContactRecoveryDiagnostics(
        bulk_delta_over_delta0=bulk,
        threshold_fraction=float(threshold_fraction),
        xi_m=xi,
        left_recovery_length_m=float(left_len),
        right_recovery_length_m=float(right_len),
        left_recovery_length_xi=float(left_xi),
        right_recovery_length_xi=float(right_xi),
        mean_recovery_length_xi=float(mean_xi),
        min_allowed_xi=float(min_allowed_xi),
        max_allowed_xi=float(max_allowed_xi),
        passes=passes,
        reason=reason,
    )


def stationarity_diagnostics(
    *,
    history: dict[str, np.ndarray],
    material: GTDGLMaterial,
    phase_gradient_rel_tol: float = 1.0e-4,
    phi_gradient_rel_tol: float = 1.0e-4,
    phase_gradient_abs_tol_m_inv: float = 1.0e3,
    phi_gradient_abs_tol_V_m: float = 1.0e2,
    edge_active_threshold: float = 0.05,
    eta_window: int = 20,
    # Deprecated aliases kept so older tests/calls fail softly instead of
    # changing semantics silently.  They are ignored for the pass/fail gate.
    delta_rel_tol: float | None = None,
    phi_rel_tol: float | None = None,
    eta_tol: float | None = None,
) -> StationarityDiagnostics:
    """Evaluate gauge-fixed temporal stationarity from final snapshots.

    The previous smoke diagnostic compared ``Delta`` and ``phi`` themselves.
    That is too gauge-sensitive for the present objective: a global phase shift
    or a constant electrostatic offset should not make a stationary state fail.

    This diagnostic therefore compares edge fields between the last two stored
    snapshots:

    * ``Q_edge = grad(arg(Delta))`` in m^-1, taken from the current adapter.
    * ``grad(phi)_edge`` in V/m, built from node phi snapshots and edge lengths.

    Edges whose final |Delta| is close to zero are excluded because the phase of
    Delta is undefined there.  This matters for metallic contacts where terminal
    sites are intentionally clamped to |Delta| = 0.
    """

    del material, delta_rel_tol, phi_rel_tol, eta_tol

    q_edge = _snapshot_2d(
        history.get("edge_phase_gradient_snapshot_m_inv", history.get("edge_Q_snapshot_m_inv", []))
    )
    phi_grad = _edge_phi_gradient_snapshots(history)
    active = _active_phase_edges(history, edge_active_threshold=edge_active_threshold)

    q_diag = _edge_field_change_metrics(q_edge, active, abs_tol=float(phase_gradient_abs_tol_m_inv))
    phi_diag = _edge_field_change_metrics(phi_grad, active, abs_tol=float(phi_gradient_abs_tol_V_m))

    eta = np.asarray(history.get("eta_R", []), dtype=float).reshape(-1)
    eta_final = float(eta[-1]) if eta.size else float("nan")
    w = max(1, int(eta_window))
    eta_window_max = float(np.nanmax(eta[-w:])) if eta.size else float("nan")

    q_pass = bool(
        np.isfinite(q_diag["rel_change"])
        and (
            q_diag["rel_change"] <= float(phase_gradient_rel_tol)
            or q_diag["abs_change"] <= float(phase_gradient_abs_tol_m_inv)
        )
    )
    phi_pass = bool(
        np.isfinite(phi_diag["rel_change"])
        and (
            phi_diag["rel_change"] <= float(phi_gradient_rel_tol)
            or phi_diag["abs_change"] <= float(phi_gradient_abs_tol_V_m)
        )
    )
    passes = bool(q_pass and phi_pass)
    if passes:
        reason = "phase-gradient and phi-gradient are stationary within requested tolerances"
    else:
        reason = "gauge-fixed gradient stationarity tolerances were not all satisfied"

    return StationarityDiagnostics(
        phase_gradient_rel_change=float(q_diag["rel_change"]),
        phi_gradient_rel_change=float(phi_diag["rel_change"]),
        phase_gradient_abs_change_m_inv=float(q_diag["abs_change"]),
        phi_gradient_abs_change_V_m=float(phi_diag["abs_change"]),
        phase_gradient_rms_final_m_inv=float(q_diag["rms_final"]),
        phi_gradient_rms_final_V_m=float(phi_diag["rms_final"]),
        active_edge_fraction=float(q_diag["active_fraction"]),
        active_edge_count=int(q_diag["active_count"]),
        total_edge_count=int(q_diag["total_count"]),
        eta_R_final=eta_final,
        eta_R_window_max=eta_window_max,
        tolerance_phase_gradient_rel=float(phase_gradient_rel_tol),
        tolerance_phi_gradient_rel=float(phi_gradient_rel_tol),
        tolerance_phase_gradient_abs_m_inv=float(phase_gradient_abs_tol_m_inv),
        tolerance_phi_gradient_abs_V_m=float(phi_gradient_abs_tol_V_m),
        edge_active_threshold_over_bulk=float(edge_active_threshold),
        passes=passes,
        reason=reason,
    )


def continuity_diagnostics(
    *,
    currents,
    mesh,
    material: GTDGLMaterial,
    target_current_A: float,
    history: dict[str, np.ndarray],
    rms_tol: float = 1.0e-6,
    max_tol: float = 1.0e-3,
    poisson_tol: float = 1.0e-9,
) -> ContinuityDiagnostics:
    """Evaluate total-current continuity using SI and native Poisson metrics."""

    rms_rel = float(current_residual(currents, mesh, material, target_current_A))
    max_rel = float(max_current_residual(currents, mesh, material, target_current_A))
    pois = np.asarray(history.get("pytdgl_like_poisson_residual_rel", []), dtype=float)
    poisson_rel = float(pois[-1]) if pois.size else float("nan")
    passes = bool(
        np.isfinite(rms_rel)
        and np.isfinite(max_rel)
        and np.isfinite(poisson_rel)
        and rms_rel <= float(rms_tol)
        and max_rel <= float(max_tol)
        and poisson_rel <= float(poisson_tol)
    )
    if passes:
        reason = "current continuity is within requested tolerances"
    else:
        reason = "current-continuity tolerances were not all satisfied"

    return ContinuityDiagnostics(
        rms_residual_rel=rms_rel,
        max_residual_rel=max_rel,
        native_poisson_residual_rel=poisson_rel,
        tolerance_rms=float(rms_tol),
        tolerance_max=float(max_tol),
        tolerance_poisson=float(poisson_tol),
        passes=passes,
        reason=reason,
    )


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


def _active_phase_edges(history: dict[str, np.ndarray], *, edge_active_threshold: float) -> np.ndarray | None:
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
