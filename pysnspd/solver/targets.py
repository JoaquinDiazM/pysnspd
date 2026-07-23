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
    bulk_exclusion_xi: float
    bulk_exclusion_length_m: float
    bulk_edge_fraction: float
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
            "diagnostic": "bulk_gauge_fixed_edge_gradients_A_eq_0_v1",
            "phase_gradient_rel_change": float(self.phase_gradient_rel_change),
            "phi_gradient_rel_change": float(self.phi_gradient_rel_change),
            "phase_gradient_abs_change_m_inv": float(self.phase_gradient_abs_change_m_inv),
            "phi_gradient_abs_change_V_m": float(self.phi_gradient_abs_change_V_m),
            "phase_gradient_rms_final_m_inv": float(self.phase_gradient_rms_final_m_inv),
            "phi_gradient_rms_final_V_m": float(self.phi_gradient_rms_final_V_m),
            "active_edge_fraction": float(self.active_edge_fraction),
            "active_edge_count": int(self.active_edge_count),
            "total_edge_count": int(self.total_edge_count),
            "bulk_exclusion_xi": float(self.bulk_exclusion_xi),
            "bulk_exclusion_length_m": float(self.bulk_exclusion_length_m),
            "bulk_edge_fraction": float(self.bulk_edge_fraction),
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
class DynamicStationarityDiagnostics:
    """Tail-envelope and condensate-morphology test for a dynamic SS state."""

    passes: bool
    reason: str
    sufficient_history: bool
    morphology_regime: str
    tail_snapshot_count: int
    tail_duration_ps: float
    suppressed_band_counts: tuple[int, ...]
    suppressed_band_count_final: int
    psl_count_final: int
    topology_count_stable: bool
    new_suppressed_band_in_tail: bool
    profile_relative_fluctuation: float
    profile_relative_drift: float
    voltage_relative_span: float
    voltage_relative_drift: float
    mean_delta_final_over_delta0: float
    normal_like_fraction_final: float
    tolerance_profile_relative: float
    tolerance_voltage_relative: float
    psl_threshold_over_delta0: float
    normal_like_fraction_threshold: float
    minimum_tail_duration_ps: float

    def as_dict(self) -> dict[str, Any]:
        return {
            "diagnostic": "dynamic_ss_tail_envelope_and_x_profile_v1",
            "passes": bool(self.passes),
            "reason": str(self.reason),
            "sufficient_history": bool(self.sufficient_history),
            "morphology_regime": str(self.morphology_regime),
            "tail_snapshot_count": int(self.tail_snapshot_count),
            "tail_duration_ps": float(self.tail_duration_ps),
            "suppressed_band_counts": [int(value) for value in self.suppressed_band_counts],
            "suppressed_band_count_final": int(self.suppressed_band_count_final),
            "psl_count_final": int(self.psl_count_final),
            "topology_count_stable": bool(self.topology_count_stable),
            "new_suppressed_band_in_tail": bool(self.new_suppressed_band_in_tail),
            "profile_relative_fluctuation": float(self.profile_relative_fluctuation),
            "profile_relative_drift": float(self.profile_relative_drift),
            "voltage_relative_span": float(self.voltage_relative_span),
            "voltage_relative_drift": float(self.voltage_relative_drift),
            "mean_delta_final_over_delta0": float(self.mean_delta_final_over_delta0),
            "normal_like_fraction_final": float(self.normal_like_fraction_final),
            "tolerance_profile_relative": float(self.tolerance_profile_relative),
            "tolerance_voltage_relative": float(self.tolerance_voltage_relative),
            "psl_threshold_over_delta0": float(self.psl_threshold_over_delta0),
            "normal_like_fraction_threshold": float(self.normal_like_fraction_threshold),
            "minimum_tail_duration_ps": float(self.minimum_tail_duration_ps),
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
    bulk_exclusion_xi: float = 4.0,
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

    The pass/fail gate is additionally restricted to a bulk-edge mask if the
    history contains ``edge_distance_from_contact_m`` and ``stationarity_xi_m``.
    By default this excludes edges within ``4 xi`` of either metallic contact,
    where a physical conversion field and normal current are expected.
    """

    del material, delta_rel_tol, phi_rel_tol, eta_tol

    q_edge = _snapshot_2d(
        history.get("edge_phase_gradient_snapshot_m_inv", history.get("edge_Q_snapshot_m_inv", []))
    )
    phi_grad = _edge_phi_gradient_snapshots(history)
    active = _active_phase_edges(
        history,
        edge_active_threshold=edge_active_threshold,
        bulk_exclusion_xi=bulk_exclusion_xi,
    )

    q_diag = _edge_field_change_metrics(q_edge, active, abs_tol=float(phase_gradient_abs_tol_m_inv))
    phi_diag = _edge_field_change_metrics(phi_grad, active, abs_tol=float(phi_gradient_abs_tol_V_m))

    eta = np.asarray(history.get("eta_R", []), dtype=float).reshape(-1)
    eta_final = float(eta[-1]) if eta.size else float("nan")
    w = max(1, int(eta_window))
    eta_window_max = float(np.nanmax(eta[-w:])) if eta.size else float("nan")

    xi_hist = np.asarray(history.get("stationarity_xi_m", []), dtype=float).reshape(-1)
    xi_m = float(xi_hist[0]) if xi_hist.size and np.isfinite(xi_hist[0]) else float("nan")
    bulk_excl = max(0.0, float(bulk_exclusion_xi))
    bulk_excl_m = bulk_excl * xi_m if np.isfinite(xi_m) else float("nan")
    bulk_frac = float(q_diag["active_fraction"])

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
        bulk_exclusion_xi=bulk_excl,
        bulk_exclusion_length_m=float(bulk_excl_m),
        bulk_edge_fraction=bulk_frac,
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


def dynamic_stationarity_diagnostics(
    *,
    history: dict[str, np.ndarray],
    nodes_m: np.ndarray,
    delta0_J: float,
    tail_snapshots: int = 4,
    minimum_tail_duration_ps: float = 2.0,
    profile_relative_tolerance: float = 5.0e-2,
    voltage_relative_tolerance: float = 5.0e-2,
    voltage_absolute_scale_V: float = 1.0e-4,
    psl_threshold_over_delta0: float = 0.75,
    normal_like_fraction_threshold: float = 0.85,
    minimum_band_width_xi: float = 0.25,
    bulk_exclusion_xi: float = 4.0,
) -> DynamicStationarityDiagnostics:
    """Identify a time-periodic or weakly oscillating dynamic SS attractor.

    Strict stationarity compares the last two gauge-fixed fields point by point
    and therefore rejects translating vortices or a persistent phase-slip
    cycle.  This complementary diagnostic asks whether the late-time envelope
    has stopped evolving:

    * the transverse-median ``|Delta|(x)`` profile has small fluctuation and
      drift over several snapshots;
    * the number of cross-width suppressed bands is unchanged throughout the
      tail, so no new phase-slip line appears there;
    * the terminal-voltage oscillation envelope and drift are small.

    Local moving vortex cores affect only a minority of nodes in an x bin and
    are intentionally suppressed by the transverse median.  A globally
    collapsed condensate is reported separately as ``normal_like`` rather than
    being misidentified as one phase-slip line.
    """

    nodes = np.asarray(nodes_m, dtype=float)
    x_m = nodes[:, 0] if nodes.ndim == 2 and nodes.shape[1] >= 1 else np.array([], dtype=float)
    snap_t_s = np.asarray(
        history.get("snapshot_t_s", history.get("delta_snapshot_t_s", [])),
        dtype=float,
    ).reshape(-1)
    psi_r = _snapshot_2d(history.get("psi_snapshot_real_J", []))
    psi_i = _snapshot_2d(history.get("psi_snapshot_imag_J", []))
    if psi_r.shape == psi_i.shape and psi_r.shape[1:] == (x_m.size,):
        delta_over_delta0 = np.hypot(psi_r, psi_i) / max(float(delta0_J), 1.0e-300)
    else:
        delta_meV = _snapshot_2d(history.get("delta_snapshot_meV", []))
        delta0_meV = float(delta0_J) / 1.602176634e-22
        delta_over_delta0 = delta_meV / max(delta0_meV, 1.0e-300)

    n_available = min(snap_t_s.size, delta_over_delta0.shape[0])
    requested_tail = max(3, int(tail_snapshots))
    n_tail = min(requested_tail, n_available)
    if n_tail > 0:
        tail_t_s = snap_t_s[n_available - n_tail : n_available]
        tail_delta = delta_over_delta0[n_available - n_tail : n_available]
    else:
        tail_t_s = np.array([], dtype=float)
        tail_delta = np.empty((0, x_m.size), dtype=float)
    tail_duration_ps = (
        float((tail_t_s[-1] - tail_t_s[0]) / 1.0e-12) if tail_t_s.size >= 2 else 0.0
    )

    xi_values = np.asarray(history.get("stationarity_xi_m", []), dtype=float).reshape(-1)
    xi_m = float(xi_values[0]) if xi_values.size and np.isfinite(xi_values[0]) else float("nan")
    profiles, bin_width_m = _transverse_delta_profiles(
        x_m=x_m,
        delta_over_delta0=tail_delta,
        xi_m=xi_m,
        bulk_exclusion_xi=float(bulk_exclusion_xi),
    )
    minimum_width_m = (
        float(minimum_band_width_xi) * xi_m if np.isfinite(xi_m) and xi_m > 0.0 else 0.0
    )
    minimum_bins = max(1, int(np.ceil(minimum_width_m / max(bin_width_m, 1.0e-300))))
    band_counts = tuple(
        _count_true_runs(profile < float(psl_threshold_over_delta0), minimum_bins=minimum_bins)
        for profile in profiles
    )
    topology_stable = bool(len(band_counts) >= 3 and len(set(band_counts)) == 1)
    new_band = bool(
        len(band_counts) >= 2
        and any(current > previous for previous, current in zip(band_counts[:-1], band_counts[1:]))
    )

    if profiles.size:
        tail_mean_profile = np.nanmean(profiles, axis=0)
        profile_scale = max(float(np.sqrt(np.nanmean(tail_mean_profile**2))), 5.0e-2)
        deviations = np.sqrt(np.nanmean((profiles - tail_mean_profile[None, :]) ** 2, axis=1))
        profile_fluctuation = float(np.nanmax(deviations) / profile_scale)
        profile_drift = float(
            np.sqrt(np.nanmean((profiles[-1] - profiles[0]) ** 2)) / profile_scale
        )
        final_profile = profiles[-1]
        normal_fraction = float(np.mean(final_profile < float(psl_threshold_over_delta0)))
        mean_delta_final = float(np.nanmean(final_profile))
    else:
        profile_fluctuation = float("nan")
        profile_drift = float("nan")
        normal_fraction = float("nan")
        mean_delta_final = float("nan")

    normal_like = bool(
        np.isfinite(normal_fraction)
        and normal_fraction >= float(normal_like_fraction_threshold)
    )
    final_band_count = int(band_counts[-1]) if band_counts else 0
    if normal_like:
        morphology_regime = "normal_like"
    elif final_band_count > 0:
        morphology_regime = "superconducting_with_suppressed_bands"
    else:
        morphology_regime = "superconducting_without_suppressed_bands"
    psl_count_final = 0 if normal_like else final_band_count

    history_t_s = np.asarray(history.get("t_s", []), dtype=float).reshape(-1)
    terminal_voltage = np.asarray(history.get("terminal_voltage_V", []), dtype=float).reshape(-1)
    if history_t_s.size and terminal_voltage.size:
        n_voltage = min(history_t_s.size, terminal_voltage.size)
        history_t_s = history_t_s[-n_voltage:]
        terminal_voltage = terminal_voltage[-n_voltage:]
        tail_start_s = float(tail_t_s[0]) if tail_t_s.size else float(history_t_s[-1])
        voltage_tail = terminal_voltage[history_t_s >= tail_start_s]
    else:
        voltage_tail = np.array([], dtype=float)
    voltage_span, voltage_drift = _tail_scalar_envelope_metrics(
        voltage_tail,
        absolute_scale=float(voltage_absolute_scale_V),
    )

    sufficient_history = bool(
        n_tail >= 3
        and profiles.shape[0] >= 3
        and tail_duration_ps >= float(minimum_tail_duration_ps)
        and voltage_tail.size >= 4
    )
    finite_metrics = bool(
        np.isfinite(profile_fluctuation)
        and np.isfinite(profile_drift)
        and np.isfinite(voltage_span)
        and np.isfinite(voltage_drift)
    )
    passes = bool(
        sufficient_history
        and finite_metrics
        and topology_stable
        and not new_band
        and profile_fluctuation <= float(profile_relative_tolerance)
        and profile_drift <= float(profile_relative_tolerance)
        and voltage_span <= float(voltage_relative_tolerance)
        and voltage_drift <= float(voltage_relative_tolerance)
    )
    if passes:
        reason = "late-time voltage envelope and condensate morphology are dynamically stationary"
    elif not sufficient_history:
        reason = "insufficient late-time duration or snapshots for dynamic-stationarity classification"
    elif not topology_stable or new_band:
        reason = "suppressed-band topology still changes in the late-time snapshot tail"
    else:
        reason = "late-time morphology or voltage envelope exceeds dynamic-stationarity tolerances"

    return DynamicStationarityDiagnostics(
        passes=passes,
        reason=reason,
        sufficient_history=sufficient_history,
        morphology_regime=morphology_regime,
        tail_snapshot_count=int(n_tail),
        tail_duration_ps=tail_duration_ps,
        suppressed_band_counts=band_counts,
        suppressed_band_count_final=final_band_count,
        psl_count_final=psl_count_final,
        topology_count_stable=topology_stable,
        new_suppressed_band_in_tail=new_band,
        profile_relative_fluctuation=profile_fluctuation,
        profile_relative_drift=profile_drift,
        voltage_relative_span=voltage_span,
        voltage_relative_drift=voltage_drift,
        mean_delta_final_over_delta0=mean_delta_final,
        normal_like_fraction_final=normal_fraction,
        tolerance_profile_relative=float(profile_relative_tolerance),
        tolerance_voltage_relative=float(voltage_relative_tolerance),
        psl_threshold_over_delta0=float(psl_threshold_over_delta0),
        normal_like_fraction_threshold=float(normal_like_fraction_threshold),
        minimum_tail_duration_ps=float(minimum_tail_duration_ps),
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

from pysnspd.solver.target_metrics import (
    _active_phase_edges,
    _count_true_runs,
    _edge_field_change_metrics,
    _edge_phi_gradient_snapshots,
    _first_binned_crossing_distance,
    _snapshot_2d,
    _tail_scalar_envelope_metrics,
    _transverse_delta_profiles,
)
