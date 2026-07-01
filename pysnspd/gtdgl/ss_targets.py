"""Stationary-state target diagnostics for the SS gTDGL smoke runs.

The routines in this module are deliberately lightweight: they do not decide
new physics, they only build a smooth metallic-contact seed and quantify the
three checks needed before moving to photon dynamics:

1. temporal stationarity of Delta and phi,
2. a contact-healing length of order a few physical coherence lengths,
3. finite-volume current continuity.
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
    """Temporal stationarity measured from the last two stored snapshots."""

    delta_rel_change: float
    phi_rel_change: float
    delta_abs_change_over_delta0: float
    phi_abs_change_V: float
    eta_R_final: float
    eta_R_window_max: float
    tolerance_delta_rel: float
    tolerance_phi_rel: float
    tolerance_eta: float
    passes: bool
    reason: str

    def as_dict(self) -> dict[str, Any]:
        return {
            "delta_rel_change": float(self.delta_rel_change),
            "phi_rel_change": float(self.phi_rel_change),
            "delta_abs_change_over_delta0": float(self.delta_abs_change_over_delta0),
            "phi_abs_change_V": float(self.phi_abs_change_V),
            "eta_R_final": float(self.eta_R_final),
            "eta_R_window_max": float(self.eta_R_window_max),
            "tolerance_delta_rel": float(self.tolerance_delta_rel),
            "tolerance_phi_rel": float(self.tolerance_phi_rel),
            "tolerance_eta": float(self.tolerance_eta),
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
    delta_rel_tol: float = 1.0e-4,
    phi_rel_tol: float = 1.0e-4,
    eta_tol: float = 1.0e-5,
    eta_window: int = 20,
) -> StationarityDiagnostics:
    """Evaluate temporal stationarity from final snapshots and eta history."""

    psi_r = np.asarray(history.get("psi_snapshot_real_J", []), dtype=float)
    psi_i = np.asarray(history.get("psi_snapshot_imag_J", []), dtype=float)
    phi = np.asarray(history.get("phi_snapshot_V", []), dtype=float)
    delta0 = max(float(material.delta0_J), 1.0e-300)

    if psi_r.ndim == 2 and psi_i.shape == psi_r.shape and psi_r.shape[0] >= 2:
        psi = psi_r + 1j * psi_i
        dpsi = psi[-1] - psi[-2]
        psi_scale = max(float(np.sqrt(np.nanmean(np.abs(psi[-1]) ** 2))), delta0)
        delta_abs_change_over_delta0 = float(np.nanmax(np.abs(dpsi)) / delta0)
        delta_rel_change = float(np.sqrt(np.nanmean(np.abs(dpsi) ** 2)) / max(psi_scale, 1.0e-300))
    else:
        delta_abs_change_over_delta0 = float("nan")
        delta_rel_change = float("nan")

    if phi.ndim == 2 and phi.shape[0] >= 2:
        dphi = phi[-1] - phi[-2]
        phi_abs_change = float(np.nanmax(np.abs(dphi)))
        phi_scale = max(float(np.ptp(phi[-1])), 1.0e-12)
        phi_rel_change = phi_abs_change / phi_scale
    else:
        phi_abs_change = float("nan")
        phi_rel_change = float("nan")

    eta = np.asarray(history.get("eta_R", []), dtype=float).reshape(-1)
    eta_final = float(eta[-1]) if eta.size else float("nan")
    w = max(1, int(eta_window))
    eta_window_max = float(np.nanmax(eta[-w:])) if eta.size else float("nan")
    passes = bool(
        np.isfinite(delta_rel_change)
        and np.isfinite(phi_rel_change)
        and np.isfinite(eta_window_max)
        and delta_rel_change <= float(delta_rel_tol)
        and phi_rel_change <= float(phi_rel_tol)
        and eta_window_max <= float(eta_tol)
    )
    if passes:
        reason = "Delta and phi are stationary within requested tolerances"
    else:
        reason = "stationarity tolerances were not all satisfied"

    return StationarityDiagnostics(
        delta_rel_change=delta_rel_change,
        phi_rel_change=phi_rel_change,
        delta_abs_change_over_delta0=delta_abs_change_over_delta0,
        phi_abs_change_V=phi_abs_change,
        eta_R_final=eta_final,
        eta_R_window_max=eta_window_max,
        tolerance_delta_rel=float(delta_rel_tol),
        tolerance_phi_rel=float(phi_rel_tol),
        tolerance_eta=float(eta_tol),
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
