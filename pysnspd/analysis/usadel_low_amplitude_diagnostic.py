"""Temporary diagnostics for the low-amplitude Usadel current ambiguity.

This module intentionally does not participate in PRE generation or in either
solver.  It lets the D1 plotting pipeline compare three constitutive closures
using an existing strict PRE current table:

``current``
    Production interpolation of ``j_s(T_e, |Delta|, q)``.
``stiffness``
    Diagnostic reconstruction of ``kappa=j_s/(|Delta|^2 q)`` interpolated on
    ``(T_e, |Delta|^2, |q|)``, with the exact ``|Delta| -> 0`` anchor.
``direct``
    Direct Matsubara evaluation at every requested point.

The stiffness path is a candidate evaluated in isolation, not a production
implementation.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from pysnspd.gtdgl.usadel_current import interpolate_strict_usadel_current_table
from pysnspd.usadel.calibration import matsubara_energy_axis_J, solve_matsubara_s_values
from pysnspd.usadel.parameters import E_CHARGE_C, HBAR_J_S, K_B_J_K


@dataclass(frozen=True)
class DiagnosticCurrentCatalog:
    """Strict PRE current table plus the material data needed for recomputation."""

    Te_axis_K: np.ndarray
    delta_axis_J: np.ndarray
    q_axis_m_inv: np.ndarray
    js_A_m2: np.ndarray
    stiffness_A_m_inv_J_inv2: np.ndarray
    D_m2_s: float
    sigma_n_S_m: float
    Tc_K: float
    delta0_J: float
    n_matsubara: int
    metadata: dict[str, Any]

    @property
    def first_positive_delta_J(self) -> float:
        positive = self.delta_axis_J[self.delta_axis_J > 0.0]
        if positive.size == 0:
            raise ValueError("The PRE delta axis has no positive entry.")
        return float(positive[0])


@dataclass(frozen=True)
class ConstitutiveCurves:
    """Current curves and compact error/slope metrics for one ``(T_e, q)``."""

    Te_K: float
    q_m_inv: float
    delta_J: np.ndarray
    current_A_m2: np.ndarray
    stiffness_A_m2: np.ndarray
    direct_A_m2: np.ndarray
    current_exponent: np.ndarray
    stiffness_exponent: np.ndarray
    direct_exponent: np.ndarray
    metrics: dict[str, float]


@dataclass(frozen=True)
class NotchDiagnostic:
    """One-dimensional notch proxy for the Allmaras phase-source singularity."""

    x_m: np.ndarray
    delta_J: np.ndarray
    current_A_m2: np.ndarray
    stiffness_A_m2: np.ndarray
    direct_A_m2: np.ndarray
    gl_A_m2: np.ndarray
    source_current_A_m3_J_inv: np.ndarray
    source_stiffness_A_m3_J_inv: np.ndarray
    source_direct_A_m3_J_inv: np.ndarray
    metrics: dict[str, float]


def load_diagnostic_current_catalog(npz_path: str | Path) -> DiagnosticCurrentCatalog:
    """Load a strict current table and construct the candidate stiffness table."""

    source = Path(npz_path)
    with np.load(source, allow_pickle=True) as data:
        table = _first_array(data, ("js_A_m2", "j_s_A_m2", "js_T_delta_q_A_m2"))
        Te_axis = _first_array(data, ("Te_axis_K", "T_axis_K", "temperature_axis_K"))
        delta_axis = _first_array(data, ("delta_axis_J", "delta_values_J", "Delta_axis_J"))
        q_axis = _first_array(data, ("q_axis_m_inv", "q_values_m_inv", "q_grid_m_inv"))
        metadata = dict(data["metadata"].item()) if "metadata" in data.files else {}
        n_matsubara = _scalar_from_npz_or_metadata(
            data,
            "js_table_n_matsubara",
            metadata,
            "n_matsubara_configured",
        )

    Te = _strict_axis(Te_axis, "Te_axis_K", positive=True)
    delta = _strict_axis(delta_axis, "delta_axis_J", nonnegative=True)
    q = _strict_axis(q_axis, "q_axis_m_inv", nonnegative=True)
    js = np.asarray(table, dtype=float)
    expected = (Te.size, delta.size, q.size)
    if js.shape != expected:
        raise ValueError(f"Expected strict js[Te,delta,q] shape {expected}, got {js.shape}.")
    if np.any(~np.isfinite(js)):
        raise ValueError("The PRE current table contains non-finite values.")
    if delta.size < 2 or delta[0] != 0.0 or q.size < 2 or q[0] != 0.0:
        raise ValueError("The diagnostic requires delta and q axes anchored exactly at zero.")

    D = _positive_metadata_float(metadata, "D_m2_s")
    sigma = _positive_metadata_float(metadata, "sigma_n_S_m")
    Tc = _positive_metadata_float(metadata, "Tc_K")
    delta0 = _positive_metadata_float(metadata, "delta0_J")
    n_m = int(n_matsubara)
    if n_m <= 0:
        raise ValueError("A positive Matsubara cutoff is required.")

    stiffness = build_candidate_stiffness_table(
        Te_axis_K=Te,
        delta_axis_J=delta,
        q_axis_m_inv=q,
        js_A_m2=js,
        D_m2_s=D,
        sigma_n_S_m=sigma,
        n_matsubara=n_m,
    )
    return DiagnosticCurrentCatalog(
        Te_axis_K=Te,
        delta_axis_J=delta,
        q_axis_m_inv=q,
        js_A_m2=js,
        stiffness_A_m_inv_J_inv2=stiffness,
        D_m2_s=D,
        sigma_n_S_m=sigma,
        Tc_K=Tc,
        delta0_J=delta0,
        n_matsubara=n_m,
        metadata=metadata,
    )


def build_candidate_stiffness_table(
    *,
    Te_axis_K: np.ndarray,
    delta_axis_J: np.ndarray,
    q_axis_m_inv: np.ndarray,
    js_A_m2: np.ndarray,
    D_m2_s: float,
    sigma_n_S_m: float,
    n_matsubara: int,
) -> np.ndarray:
    """Re-express an existing table as ``kappa(T, delta^2, |q|)``.

    Positive ``delta``/``q`` nodes are an algebraic change of variables.  The
    zero-amplitude row is filled from the exact Matsubara expansion.  The
    zero-q column uses its analytic limit, avoiding a numerical division by q.
    """

    Te = np.asarray(Te_axis_K, dtype=float)
    delta = np.asarray(delta_axis_J, dtype=float)
    q = np.asarray(q_axis_m_inv, dtype=float)
    js = np.asarray(js_A_m2, dtype=float)
    expected = (Te.size, delta.size, q.size)
    if js.shape != expected:
        raise ValueError(f"Expected js shape {expected}, got {js.shape}.")

    stiffness = np.empty_like(js)
    denom = delta[None, :, None] ** 2 * q[None, None, :]
    np.divide(js, denom, out=stiffness, where=denom > 0.0)

    gamma = 0.5 * HBAR_J_S * float(D_m2_s) * q * q
    for iT, temperature in enumerate(Te):
        eps = matsubara_energy_axis_J(T_K=float(temperature), n_matsubara=int(n_matsubara))
        prefactor = 2.0 * np.pi * K_B_J_K * float(temperature) * float(sigma_n_S_m) / E_CHARGE_C

        # Exact Delta -> 0 coefficient for every q node.
        stiffness[iT, 0, :] = prefactor * np.sum(
            1.0 / (eps[None, :] + gamma[:, None]) ** 2,
            axis=1,
        )

        # At Gamma=q=0, s=Delta/sqrt(eps^2+Delta^2) exactly.
        stiffness[iT, :, 0] = prefactor * np.sum(
            1.0 / (eps[None, :] ** 2 + delta[:, None] ** 2),
            axis=1,
        )

    if np.any(~np.isfinite(stiffness)) or np.any(stiffness <= 0.0):
        raise ValueError("Candidate stiffness table is not finite and strictly positive.")
    return stiffness


def compare_constitutive_curves(
    catalog: DiagnosticCurrentCatalog,
    *,
    Te_K: float,
    q_m_inv: float,
    delta_J: np.ndarray,
) -> ConstitutiveCurves:
    """Evaluate current interpolation, candidate stiffness, and direct Usadel."""

    delta = np.asarray(delta_J, dtype=float).reshape(-1)
    if np.any(delta <= 0.0):
        raise ValueError("Constitutive diagnostic points must have delta > 0.")
    q = float(q_m_inv)
    T = float(Te_K)
    if q <= 0.0:
        raise ValueError("The constitutive comparison requires q > 0.")

    old = interpolate_strict_usadel_current_table(
        table=catalog.js_A_m2,
        Te_axis_K=catalog.Te_axis_K,
        delta_axis_J=catalog.delta_axis_J,
        q_axis_m_inv=catalog.q_axis_m_inv,
        q_edge_m_inv=np.full(delta.shape, q),
        delta_edge_J=delta,
        Te_edge_K=np.full(delta.shape, T),
    )
    candidate_kappa = _trilinear(
        catalog.stiffness_A_m_inv_J_inv2,
        catalog.Te_axis_K,
        catalog.delta_axis_J**2,
        catalog.q_axis_m_inv,
        np.full(delta.shape, T),
        delta**2,
        np.full(delta.shape, q),
    )
    candidate = candidate_kappa * delta * delta * q
    direct = direct_matsubara_current_density(
        delta_J=delta,
        Te_K=T,
        q_m_inv=q,
        D_m2_s=catalog.D_m2_s,
        sigma_n_S_m=catalog.sigma_n_S_m,
        n_matsubara=catalog.n_matsubara,
    )

    p_old = local_log_exponent(delta, old)
    p_candidate = local_log_exponent(delta, candidate)
    p_direct = local_log_exponent(delta, direct)
    first = catalog.first_positive_delta_J
    fit = (delta >= 0.01 * first) & (delta <= 0.25 * first)
    metrics = {
        "Te_K": T,
        "q_m_inv": q,
        "q_um_inv": q * 1.0e-6,
        "first_positive_delta_over_delta0": first / catalog.delta0_J,
        "current_low_amplitude_exponent": _median_finite(p_old[fit]),
        "stiffness_low_amplitude_exponent": _median_finite(p_candidate[fit]),
        "direct_low_amplitude_exponent": _median_finite(p_direct[fit]),
        "current_max_relative_error_below_first_node": _max_relative_error(old[delta <= first], direct[delta <= first]),
        "stiffness_max_relative_error_below_first_node": _max_relative_error(
            candidate[delta <= first],
            direct[delta <= first],
        ),
    }
    return ConstitutiveCurves(
        Te_K=T,
        q_m_inv=q,
        delta_J=delta,
        current_A_m2=old,
        stiffness_A_m2=candidate,
        direct_A_m2=direct,
        current_exponent=p_old,
        stiffness_exponent=p_candidate,
        direct_exponent=p_direct,
        metrics=metrics,
    )


def direct_matsubara_current_density(
    *,
    delta_J: np.ndarray,
    Te_K: float,
    q_m_inv: float,
    D_m2_s: float,
    sigma_n_S_m: float,
    n_matsubara: int,
) -> np.ndarray:
    """Direct dirty-limit Matsubara current for an arbitrary delta vector."""

    delta = np.asarray(delta_J, dtype=float).reshape(-1)
    T = float(Te_K)
    q = float(q_m_inv)
    eps = matsubara_energy_axis_J(T_K=T, n_matsubara=int(n_matsubara))
    gamma = 0.5 * HBAR_J_S * float(D_m2_s) * q * q
    prefactor = 2.0 * np.pi * K_B_J_K * T * float(sigma_n_S_m) / E_CHARGE_C
    out = np.empty_like(delta)
    for index, amplitude in enumerate(delta):
        s = solve_matsubara_s_values(delta_J=float(amplitude), gamma_J=gamma, eps_n_J=eps)
        out[index] = prefactor * q * float(np.sum(s * s))
    return out


def build_notch_diagnostic(
    catalog: DiagnosticCurrentCatalog,
    *,
    Te_K: float,
    q_m_inv: float,
    n_points: int = 801,
    half_width_nm: float = 36.0,
    notch_sigma_nm: float = 8.0,
    minimum_delta_fraction_of_first_node: float = 1.0e-4,
) -> NotchDiagnostic:
    """Build a smooth 1D notch and the phase-source factor ``d(jUs-jGL)/dx/R``."""

    x = np.linspace(-float(half_width_nm), float(half_width_nm), int(n_points)) * 1.0e-9
    first = catalog.first_positive_delta_J
    floor = float(minimum_delta_fraction_of_first_node) * first
    outside = 1.25 * first
    sigma = float(notch_sigma_nm) * 1.0e-9
    delta = floor + (outside - floor) * (1.0 - np.exp(-0.5 * (x / sigma) ** 2))
    q = float(q_m_inv)
    T = float(Te_K)
    current = interpolate_strict_usadel_current_table(
        table=catalog.js_A_m2,
        Te_axis_K=catalog.Te_axis_K,
        delta_axis_J=catalog.delta_axis_J,
        q_axis_m_inv=catalog.q_axis_m_inv,
        q_edge_m_inv=np.full(delta.shape, q),
        delta_edge_J=delta,
        Te_edge_K=np.full(delta.shape, T),
    )
    candidate_kappa = _trilinear(
        catalog.stiffness_A_m_inv_J_inv2,
        catalog.Te_axis_K,
        catalog.delta_axis_J**2,
        catalog.q_axis_m_inv,
        np.full(delta.shape, T),
        delta**2,
        np.full(delta.shape, q),
    )
    candidate = candidate_kappa * delta * delta * q
    direct = direct_matsubara_current_density(
        delta_J=delta,
        Te_K=T,
        q_m_inv=q,
        D_m2_s=catalog.D_m2_s,
        sigma_n_S_m=catalog.sigma_n_S_m,
        n_matsubara=catalog.n_matsubara,
    )

    kappa_gl = (
        np.pi
        * catalog.sigma_n_S_m
        / (4.0 * E_CHARGE_C * K_B_J_K * catalog.Tc_K)
    )
    gl = kappa_gl * delta * delta * float(q_m_inv)
    source_old = np.gradient(current - gl, x, edge_order=2) / delta
    source_candidate = np.gradient(candidate - gl, x, edge_order=2) / delta
    source_direct = np.gradient(direct - gl, x, edge_order=2) / delta
    exact_scale = max(float(np.max(np.abs(source_direct))), np.finfo(float).tiny)
    metrics = {
        "Te_K": float(Te_K),
        "q_m_inv": float(q_m_inv),
        "q_um_inv": float(q_m_inv) * 1.0e-6,
        "notch_min_delta_over_delta0": float(np.min(delta) / catalog.delta0_J),
        "notch_outer_delta_over_delta0": float(np.max(delta) / catalog.delta0_J),
        "current_source_peak_over_direct_peak": float(np.max(np.abs(source_old)) / exact_scale),
        "stiffness_source_peak_over_direct_peak": float(np.max(np.abs(source_candidate)) / exact_scale),
        "current_source_rms_relative_to_direct": _relative_rms(source_old, source_direct),
        "stiffness_source_rms_relative_to_direct": _relative_rms(source_candidate, source_direct),
    }
    return NotchDiagnostic(
        x_m=x,
        delta_J=delta,
        current_A_m2=current,
        stiffness_A_m2=candidate,
        direct_A_m2=direct,
        gl_A_m2=gl,
        source_current_A_m3_J_inv=source_old,
        source_stiffness_A_m3_J_inv=source_candidate,
        source_direct_A_m3_J_inv=source_direct,
        metrics=metrics,
    )


def local_log_exponent(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Return ``d log(|y|)/d log(x)`` with finite positive safeguards."""

    xx = np.asarray(x, dtype=float)
    yy = np.abs(np.asarray(y, dtype=float))
    tiny = np.finfo(float).tiny
    return np.gradient(np.log(np.maximum(yy, tiny)), np.log(np.maximum(xx, tiny)), edge_order=2)


def _trilinear(
    table: np.ndarray,
    axis0: np.ndarray,
    axis1: np.ndarray,
    axis2: np.ndarray,
    x0: np.ndarray,
    x1: np.ndarray,
    x2: np.ndarray,
) -> np.ndarray:
    i0, j0, w0 = _bracket(axis0, x0)
    i1, j1, w1 = _bracket(axis1, x1)
    i2, j2, w2 = _bracket(axis2, x2)
    c00 = table[i0, i1, i2] * (1.0 - w2) + table[i0, i1, j2] * w2
    c01 = table[i0, j1, i2] * (1.0 - w2) + table[i0, j1, j2] * w2
    c10 = table[j0, i1, i2] * (1.0 - w2) + table[j0, i1, j2] * w2
    c11 = table[j0, j1, i2] * (1.0 - w2) + table[j0, j1, j2] * w2
    c0 = c00 * (1.0 - w1) + c01 * w1
    c1 = c10 * (1.0 - w1) + c11 * w1
    return c0 * (1.0 - w0) + c1 * w0


def _bracket(axis: np.ndarray, points: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    values = np.asarray(axis, dtype=float)
    x = np.clip(np.asarray(points, dtype=float), values[0], values[-1])
    if values.size == 1:
        zero = np.zeros(x.shape, dtype=np.int64)
        return zero, zero, np.zeros(x.shape, dtype=float)
    upper = np.clip(np.searchsorted(values, x, side="right"), 1, values.size - 1)
    lower = upper - 1
    weight = (x - values[lower]) / np.maximum(values[upper] - values[lower], 1.0e-300)
    return lower, upper, np.clip(weight, 0.0, 1.0)


def _first_array(data: Any, names: tuple[str, ...]) -> np.ndarray:
    for name in names:
        if name in data.files:
            return np.asarray(data[name], dtype=float)
    raise ValueError(f"None of the required PRE arrays is present: {', '.join(names)}.")


def _strict_axis(
    values: np.ndarray,
    name: str,
    *,
    positive: bool = False,
    nonnegative: bool = False,
) -> np.ndarray:
    axis = np.asarray(values, dtype=float).reshape(-1)
    if axis.size == 0 or np.any(~np.isfinite(axis)) or np.any(np.diff(axis) <= 0.0):
        raise ValueError(f"{name} must be finite and strictly increasing.")
    if positive and np.any(axis <= 0.0):
        raise ValueError(f"{name} must be positive.")
    if nonnegative and np.any(axis < 0.0):
        raise ValueError(f"{name} must be nonnegative.")
    return axis


def _positive_metadata_float(metadata: dict[str, Any], key: str) -> float:
    value = float(metadata.get(key, np.nan))
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError(f"PRE metadata requires positive finite {key}.")
    return value


def _scalar_from_npz_or_metadata(data: Any, npz_key: str, metadata: dict[str, Any], metadata_key: str) -> int:
    if npz_key in data.files:
        return int(np.asarray(data[npz_key]).reshape(()).item())
    return int(metadata.get(metadata_key, 0))


def _median_finite(values: np.ndarray) -> float:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    return float(np.median(finite)) if finite.size else float("nan")


def _max_relative_error(candidate: np.ndarray, reference: np.ndarray) -> float:
    ref = np.asarray(reference, dtype=float)
    scale = np.maximum(np.abs(ref), np.finfo(float).tiny)
    return float(np.max(np.abs(np.asarray(candidate, dtype=float) - ref) / scale))


def _relative_rms(candidate: np.ndarray, reference: np.ndarray) -> float:
    cand = np.asarray(candidate, dtype=float)
    ref = np.asarray(reference, dtype=float)
    denom = max(float(np.sqrt(np.mean(ref * ref))), np.finfo(float).tiny)
    return float(np.sqrt(np.mean((cand - ref) ** 2)) / denom)


__all__ = [
    "ConstitutiveCurves",
    "DiagnosticCurrentCatalog",
    "NotchDiagnostic",
    "build_candidate_stiffness_table",
    "build_notch_diagnostic",
    "compare_constitutive_curves",
    "direct_matsubara_current_density",
    "load_diagnostic_current_catalog",
    "local_log_exponent",
]
