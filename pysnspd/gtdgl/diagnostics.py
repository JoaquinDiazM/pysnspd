"""Scalar diagnostics and seed/material helpers for OE7 relaxation."""
from __future__ import annotations

import numpy as np

from pysnspd.gtdgl.material import GTDGLMaterial
from pysnspd.gtdgl.state import CurrentFields


def current_residual(
    currents: CurrentFields,
    mesh,
    material: GTDGLMaterial | None = None,
    target_current_A: float | None = None,
) -> float:
    """Dimensionless RMS residual of div(j_tot).

    New notebook-order runs pass ``material`` and ``target_current_A`` so the
    scale is the imposed average current density divided by the mesh spacing.
    Older smoke tests called this with only ``(currents, mesh)``; in that case
    we fall back to the reconstructed total-current RMS scale.
    """
    div = np.asarray(currents.node_div_jtot_A_m3, dtype=float)
    h = float(getattr(mesh, "target_spacing_m", getattr(mesh, "xi_mesh_m", 1.0e-9)))
    if material is not None and target_current_A is not None:
        jscale = abs(target_current_density_A_m2(material, float(target_current_A)))
    else:
        jscale = float(
            np.sqrt(
                np.nanmean(
                    currents.node_jtot_x_A_m2**2 + currents.node_jtot_y_A_m2**2
                )
            )
        )
    scale = max(jscale / max(h, 1.0e-300), 1.0)
    return float(np.sqrt(np.nanmean(div * div)) / scale)


def max_current_residual(
    currents: CurrentFields,
    mesh,
    material: GTDGLMaterial | None = None,
    target_current_A: float | None = None,
) -> float:
    div = np.asarray(currents.node_div_jtot_A_m3, dtype=float)
    h = float(getattr(mesh, "target_spacing_m", getattr(mesh, "xi_mesh_m", 1.0e-9)))
    if material is not None and target_current_A is not None:
        jscale = abs(target_current_density_A_m2(material, float(target_current_A)))
    else:
        jscale = float(
            np.sqrt(
                np.nanmean(
                    currents.node_jtot_x_A_m2**2 + currents.node_jtot_y_A_m2**2
                )
            )
        )
    scale = max(jscale / max(h, 1.0e-300), 1.0)
    return float(np.nanmax(np.abs(div)) / scale)


def normal_current_fraction_rms(currents: CurrentFields) -> float:
    """RMS normal-current fraction relative to total-current edge scale."""
    num = float(np.sqrt(np.nanmean(currents.edge_jn_A_m2**2))) if currents.edge_jn_A_m2.size else 0.0
    den = float(np.sqrt(np.nanmean(currents.edge_jtot_A_m2**2))) if currents.edge_jtot_A_m2.size else 0.0
    return num / max(den, 1.0e-300)


def current_density_maxima_A_m2(currents: CurrentFields) -> tuple[float, float]:
    """Return max |j_n| and max |j_tot| from edge fields."""
    jn_max = float(np.nanmax(np.abs(currents.edge_jn_A_m2))) if currents.edge_jn_A_m2.size else 0.0
    jt_max = float(np.nanmax(np.abs(currents.edge_jtot_A_m2))) if currents.edge_jtot_A_m2.size else 0.0
    return jn_max, jt_max


def normal_current_fraction_max(currents: CurrentFields) -> float:
    jn_max, jt_max = current_density_maxima_A_m2(currents)
    return jn_max / max(jt_max, 1.0e-300)


def seed_target_current_A(seed) -> float:
    """Extract imposed transport current from an OE6 seed-like object."""
    for name in ("I_bias_A", "target_current_A", "current_A"):
        if hasattr(seed, name):
            value = float(getattr(seed, name))
            if np.isfinite(value):
                return value
    metadata = getattr(seed, "metadata", None)
    if isinstance(metadata, dict):
        for name in ("I_bias_A", "target_current_A", "current_A"):
            if name in metadata:
                value = float(metadata[name])
                if np.isfinite(value):
                    return value
    return 0.0


def seed_q_bias_m_inv(seed, *, target_current_A: float | None = None) -> float:
    """Extract seed phase-gradient q."""
    for name in ("q_bias_m_inv", "target_q_m_inv", "q_m_inv"):
        if hasattr(seed, name):
            value = float(getattr(seed, name))
            if np.isfinite(value):
                return value
    metadata = getattr(seed, "metadata", None)
    if isinstance(metadata, dict):
        for name in ("q_bias_m_inv", "target_q_m_inv", "q_m_inv"):
            if name in metadata:
                value = float(metadata[name])
                if np.isfinite(value):
                    return value
    if target_current_A is not None and abs(float(target_current_A)) <= 0.0:
        return 0.0
    return 0.0


def seed_delta_bias_J(seed, *, fallback: float) -> float:
    """Extract stationary terminal amplitude from the OE6 seed if available."""
    for name in ("delta_bias_J", "target_delta_J", "Delta_bias_J"):
        if hasattr(seed, name):
            value = float(getattr(seed, name))
            if np.isfinite(value) and value >= 0.0:
                return value
    metadata = getattr(seed, "metadata", None)
    if isinstance(metadata, dict):
        for name in ("delta_bias_J", "target_delta_J", "Delta_bias_J"):
            if name in metadata:
                value = float(metadata[name])
                if np.isfinite(value) and value >= 0.0:
                    return value
    if hasattr(seed, "node_delta_J"):
        arr = np.asarray(seed.node_delta_J, dtype=float)
        finite = arr[np.isfinite(arr)]
        if finite.size:
            return float(np.nanmedian(finite))
    return float(fallback)


def target_current_density_A_m2(material: GTDGLMaterial, target_current_A: float) -> float:
    return float(target_current_A) / max(material.width_m * material.thickness_m, 1.0e-300)


def suggest_next_dt(
    *,
    dt_s: float,
    max_amp2_change_rel: float,
    retries: int,
    adaptive: bool,
    target: float,
    shrink_factor: float,
    grow_factor: float,
    dt_min_s: float,
    dt_max_s: float,
) -> float:
    """Notebook-style adaptive dt rule for accepted SS steps.

    Hard failures are handled before a step is accepted.  This function only
    chooses the next tentative step from the diagnostics of the accepted step.
    """
    if not adaptive:
        return float(dt_s)
    if retries > 0 or max_amp2_change_rel > 0.75 * target:
        return max(float(dt_s) * float(shrink_factor), float(dt_min_s))
    if max_amp2_change_rel < 0.20 * target:
        return min(float(dt_s) * float(grow_factor), float(dt_max_s))
    return min(float(dt_s), float(dt_max_s))


def stationary_trial_rejection_reason(
    *,
    amp2_change_rel: float,
    edge_pairbreaking_max: float,
    edge_js_over_javg: float,
    edge_jtot_over_javg: float,
    eta_reject: float,
    max_pairbreaking_accept: float,
    max_js_over_javg_accept: float,
    max_jtot_over_javg_accept: float,
) -> str | None:
    """Return a reason for rejecting a trial SS step, or ``None``.

    OE7 is a stationary relaxation toward a superconducting branch.  A trial
    that already contains a local depairing spike or a huge edge-current spike
    is treated like a failed adaptive step, not as a physical event.  This is
    intentionally stricter than the future photon dynamics, where such events
    may be physical and should be handled by the photon solver instead.
    """
    vals = (
        float(amp2_change_rel),
        float(edge_pairbreaking_max),
        float(edge_js_over_javg),
        float(edge_jtot_over_javg),
    )
    if not all(np.isfinite(v) for v in vals):
        return "nonfinite_trial_diagnostic"
    if amp2_change_rel > float(eta_reject):
        return f"amp2_change_rel>{eta_reject:.3e}"
    if edge_pairbreaking_max > float(max_pairbreaking_accept):
        return f"pairbreaking>{max_pairbreaking_accept:.3g}"
    if edge_js_over_javg > float(max_js_over_javg_accept):
        return f"js_over_javg>{max_js_over_javg_accept:.3g}"
    if edge_jtot_over_javg > float(max_jtot_over_javg_accept):
        return f"jtot_over_javg>{max_jtot_over_javg_accept:.3g}"
    return None

