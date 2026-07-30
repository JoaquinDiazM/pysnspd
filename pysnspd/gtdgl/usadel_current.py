"""Regular gauge-invariant Usadel edge current for the gTDGL backend.

The production SS and photon paths require the PRE stiffness table

    js_stiffness_A_per_m_J2[Te, |Delta|^2, |q|].

It is combined with the regular edge pair flow

    P_ij = Im(conj(Delta_i) U_ij Delta_j) / ell_ij

to form ``j_ij = kappa_ij P_ij``.  This expression is gauge invariant, remains
regular at suppressed amplitudes, and is exactly zero when either endpoint has
zero order-parameter amplitude.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from pysnspd.gtdgl.material import GTDGLMaterial
from pysnspd.mesh.operators import (
    FVOperators,
    divergence_from_edge_scalar,
    edge_average,
    edge_scalar_to_node_vector_least_squares,
)


@dataclass(frozen=True)
class UsadelSupercurrentDiagnostics:
    """Edge/node fields from the production stiffness closure."""

    available: bool
    backend: str
    reason: str
    edge_q_m_inv: np.ndarray
    edge_delta_J: np.ndarray
    edge_delta2_J2: np.ndarray
    edge_Te_K: np.ndarray
    edge_pair_flow_J2_m_inv: np.ndarray
    edge_stiffness_A_per_m_J2: np.ndarray
    edge_js_usadel_A_m2: np.ndarray
    edge_js_gl_A_m2: np.ndarray
    edge_js_mismatch_A_m2: np.ndarray
    node_js_usadel_x_A_m2: np.ndarray
    node_js_usadel_y_A_m2: np.ndarray
    node_div_js_usadel_A_m3: np.ndarray
    node_div_js_gl_A_m3: np.ndarray
    node_mismatch_divergence_A_m3: np.ndarray


class UsadelCatalogWithSupercurrentTable:
    """Thin adapter exposing numeric PRE sidecar arrays on a loaded catalogue."""

    def __init__(self, base: Any, arrays: Mapping[str, np.ndarray]):
        self.base = base
        self.arrays = dict(arrays)

    @property
    def files(self) -> list[str]:
        names: list[str] = []
        try:
            names.extend(list(self.base.files))  # type: ignore[attr-defined]
        except Exception:
            pass
        names.extend(self.arrays.keys())
        return list(dict.fromkeys(names))

    def __getitem__(self, key: str) -> Any:
        if key in self.arrays:
            return self.arrays[key]
        return self.base[key]

    def __getattr__(self, name: str) -> Any:
        if name in self.arrays:
            return self.arrays[name]
        return getattr(self.base, name)


def attach_usadel_supercurrent_table_from_npz(catalog: Any, npz_path: str | bytes | Path) -> Any:
    """Attach strict PRE stiffness/current resources to a loaded catalogue."""

    arrays = load_usadel_supercurrent_table_arrays_npz(npz_path)
    if not arrays:
        return catalog
    return UsadelCatalogWithSupercurrentTable(catalog, arrays)


def load_usadel_supercurrent_table_arrays_npz(npz_path: str | bytes | Path) -> dict[str, np.ndarray]:
    """Load the canonical PRE Usadel stiffness/current arrays from an NPZ."""

    wanted = (
        "js_A_m2",
        "js_stiffness_A_per_m_J2",
        "q_axis_m_inv",
        "delta_axis_J",
        "delta2_axis_J2",
        "Te_axis_K",
        "js_table_layout",
        "js_table_backend",
        "js_table_n_matsubara",
    )
    out: dict[str, np.ndarray] = {}
    with np.load(npz_path, allow_pickle=False) as data:
        keys = set(data.files)
        for key in wanted:
            if key not in keys:
                continue
            arr = np.asarray(data[key])
            if arr.dtype.kind in "OUS" and key not in {"js_table_layout", "js_table_backend"}:
                continue
            out[key] = arr
    return out


def validate_strict_usadel_supercurrent_table_npz(npz_path: str | bytes | Path) -> dict[str, Any]:
    """Validate the v2 PRE stiffness/current constitutive contract.

    Returns a compact summary if valid; raises ``RuntimeError`` otherwise.
    """

    arrays = load_usadel_supercurrent_table_arrays_npz(npz_path)
    stiffness = arrays.get("js_stiffness_A_per_m_J2")
    current = arrays.get("js_A_m2")
    q_axis = arrays.get("q_axis_m_inv")
    delta_axis = arrays.get("delta_axis_J")
    delta2_axis = arrays.get("delta2_axis_J2")
    Te_axis = arrays.get("Te_axis_K")
    missing = [
        name
        for name, arr in (
            ("js_stiffness_A_per_m_J2", stiffness),
            ("js_A_m2", current),
            ("Te_axis_K", Te_axis),
            ("delta_axis_J", delta_axis),
            ("delta2_axis_J2", delta2_axis),
            ("q_axis_m_inv", q_axis),
        )
        if arr is None
    ]
    if missing:
        raise RuntimeError(
            "PRE Usadel stiffness table does not satisfy the v2 contract. Missing: "
            + ", ".join(missing)
            + ". Re-run 01_prerun_template.py with the Matsubara stiffness builder."
        )

    stiffness = np.asarray(stiffness, dtype=float)
    current = np.asarray(current, dtype=float)
    Te_axis = _clean_axis(Te_axis, name="Te_axis_K", positive=True)
    delta_axis = _clean_axis(delta_axis, name="delta_axis_J", nonnegative=True)
    delta2_axis = _clean_axis(delta2_axis, name="delta2_axis_J2", nonnegative=True)
    q_axis = _clean_axis(q_axis, name="q_axis_m_inv", nonnegative=True)
    expected = (Te_axis.size, delta_axis.size, q_axis.size)
    if delta2_axis.size != delta_axis.size or not np.allclose(
        delta2_axis,
        delta_axis * delta_axis,
        rtol=32.0 * np.finfo(float).eps,
        atol=0.0,
    ):
        raise RuntimeError("delta2_axis_J2 must equal delta_axis_J**2 point by point.")
    if stiffness.ndim != 3 or stiffness.shape != expected:
        raise RuntimeError(
            "PRE Usadel stiffness table must have layout "
            "js_stiffness_A_per_m_J2[Te,delta2,q] "
            f"with shape {expected}; got ndim={stiffness.ndim}, shape={stiffness.shape}."
        )
    if current.ndim != 3 or current.shape != expected:
        raise RuntimeError(f"Diagnostic js_A_m2 must have shape {expected}; got {current.shape}.")
    if np.any(~np.isfinite(stiffness)) or np.any(stiffness <= 0.0):
        raise RuntimeError("PRE Usadel stiffness table must be finite and strictly positive.")
    if np.any(~np.isfinite(current)):
        raise RuntimeError("PRE Usadel diagnostic current table contains non-finite values.")
    if Te_axis.size < 1 or delta_axis.size < 2 or q_axis.size < 2:
        raise RuntimeError(
            "PRE Usadel stiffness axes are too small for local interpolation: "
            f"n_Te={Te_axis.size}, n_delta={delta_axis.size}, n_q={q_axis.size}."
        )
    if delta_axis[0] != 0.0 or delta2_axis[0] != 0.0 or q_axis[0] != 0.0:
        raise RuntimeError("Stiffness axes must include exact delta=0, delta2=0, and q=0 nodes.")

    return {
        "valid": True,
        "layout": "Te,delta2,q",
        "shape": list(stiffness.shape),
        "n_Te": int(Te_axis.size),
        "n_delta": int(delta_axis.size),
        "n_q": int(q_axis.size),
        "Te_min_K": float(np.min(Te_axis)),
        "Te_max_K": float(np.max(Te_axis)),
        "delta_min_J": float(np.min(delta_axis)),
        "delta_max_J": float(np.max(delta_axis)),
        "q_min_m_inv": float(np.min(q_axis)),
        "q_max_m_inv": float(np.max(q_axis)),
        "stiffness_min_A_per_m_J2": float(np.min(stiffness)),
        "stiffness_max_A_per_m_J2": float(np.max(stiffness)),
        "backend": str(np.asarray(arrays.get("js_table_backend", np.array("unknown"))).reshape(()).item())
        if "js_table_backend" in arrays and np.asarray(arrays["js_table_backend"]).shape == ()
        else "unknown",
    }


def compute_usadel_supercurrent_diagnostic(
    *,
    usadel_catalog: Any | None,
    psi_dimensionless: np.ndarray,
    material: GTDGLMaterial,
    Te_K: np.ndarray,
    ops: FVOperators,
    blocked_edge_mask: np.ndarray | None = None,
    edge_link_variable: np.ndarray | None = None,
) -> UsadelSupercurrentDiagnostics:
    """Evaluate the v2 stiffness closure on regular gauge-invariant FV edges."""

    psi = np.asarray(psi_dimensionless, dtype=np.complex128)
    Te = np.asarray(Te_K, dtype=float)
    if psi.shape != (ops.n_nodes,) or Te.shape != (ops.n_nodes,):
        raise ValueError("psi_dimensionless and Te_K must match the FV node count.")
    if np.any(~np.isfinite(psi)) or np.any(~np.isfinite(Te)):
        raise FloatingPointError("Usadel edge-current inputs must be finite.")

    delta_node_J = psi * float(material.delta0_J)
    edge_q = gauge_invariant_edge_q_m_inv(
        psi_dimensionless=psi,
        ops=ops,
        edge_link_variable=edge_link_variable,
    )
    edge_pair_flow = gauge_invariant_edge_pair_flow_J2_m_inv(
        delta_node_J=delta_node_J,
        ops=ops,
        edge_link_variable=edge_link_variable,
    )
    edge_delta2 = np.abs(delta_node_J[ops.edge_i]) * np.abs(delta_node_J[ops.edge_j])
    edge_delta = np.sqrt(edge_delta2)
    edge_Te = edge_average(Te, ops)

    if usadel_catalog is None:
        return _unavailable(
            reason="Usadel catalogue was not supplied.",
            edge_q=edge_q,
            edge_delta=edge_delta,
            edge_delta2=edge_delta2,
            edge_Te=edge_Te,
            edge_pair_flow=edge_pair_flow,
            ops=ops,
        )

    table = _find_first_array(usadel_catalog, ("js_stiffness_A_per_m_J2",))
    q_axis = _find_first_array(usadel_catalog, ("q_axis_m_inv",))
    delta2_axis = _find_first_array(usadel_catalog, ("delta2_axis_J2",))
    Te_axis = _find_first_array(usadel_catalog, ("Te_axis_K",))

    try:
        edge_stiffness = interpolate_strict_usadel_stiffness_table(
            table=np.asarray(table, dtype=float) if table is not None else None,
            Te_axis_K=Te_axis,
            delta2_axis_J2=delta2_axis,
            q_axis_m_inv=q_axis,
            q_edge_m_inv=edge_q,
            delta2_edge_J2=edge_delta2,
            Te_edge_K=edge_Te,
        )
    except Exception as exc:
        return _unavailable(
            reason=f"Strict 3D Usadel stiffness table is unavailable/invalid: {exc}",
            edge_q=edge_q,
            edge_delta=edge_delta,
            edge_delta2=edge_delta2,
            edge_Te=edge_Te,
            edge_pair_flow=edge_pair_flow,
            ops=ops,
        )

    edge_js = edge_stiffness * edge_pair_flow
    gl_stiffness = (
        np.pi
        * float(material.sigma_n_S_m)
        / (4.0 * 1.602176634e-19 * 1.380649e-23 * float(material.Tc_K))
    )
    edge_js_gl = gl_stiffness * edge_pair_flow
    edge_js = _apply_blocked_edges(edge_js, blocked_edge_mask)
    edge_js_gl = _apply_blocked_edges(edge_js_gl, blocked_edge_mask)
    edge_mismatch = edge_js - edge_js_gl
    if (
        np.any(~np.isfinite(edge_stiffness))
        or np.any(~np.isfinite(edge_js))
        or np.any(~np.isfinite(edge_js_gl))
        or np.any(~np.isfinite(edge_mismatch))
    ):
        raise FloatingPointError("Regular Usadel edge-current closure produced non-finite values.")
    return _finish_available(
        backend="stiffness:Te,delta2,q:gauge_invariant_edge_pair_flow",
        edge_q=edge_q,
        edge_delta=edge_delta,
        edge_delta2=edge_delta2,
        edge_Te=edge_Te,
        edge_pair_flow=edge_pair_flow,
        edge_stiffness=edge_stiffness,
        edge_js=edge_js,
        edge_js_gl=edge_js_gl,
        edge_mismatch=edge_mismatch,
        ops=ops,
    )


def interpolate_strict_usadel_stiffness_table(
    *,
    table: np.ndarray | None,
    Te_axis_K: np.ndarray | None,
    delta2_axis_J2: np.ndarray | None,
    q_axis_m_inv: np.ndarray | None,
    q_edge_m_inv: np.ndarray,
    delta2_edge_J2: np.ndarray,
    Te_edge_K: np.ndarray,
) -> np.ndarray:
    """Vectorized trilinear interpolation of ``kappa[Te, delta2, |q|]``."""

    if table is None:
        raise ValueError("js_stiffness_A_per_m_J2 table not found")
    Te_axis = _clean_axis(Te_axis_K, name="Te_axis_K", positive=True)
    delta2_axis = _clean_axis(delta2_axis_J2, name="delta2_axis_J2", nonnegative=True)
    q_axis = _clean_axis(q_axis_m_inv, name="q_axis_m_inv", nonnegative=True)
    table = np.asarray(table, dtype=float)
    expected = (Te_axis.size, delta2_axis.size, q_axis.size)
    if table.ndim != 3 or table.shape != expected:
        raise ValueError(
            "expected js_stiffness_A_per_m_J2[Te,delta2,q] "
            f"shape {expected}, got {table.shape}"
        )
    if np.any(~np.isfinite(table)) or np.any(table <= 0.0):
        raise ValueError("stiffness table must be finite and strictly positive")

    q = np.abs(np.asarray(q_edge_m_inv, dtype=float).reshape(-1))
    delta2 = np.asarray(delta2_edge_J2, dtype=float).reshape(q.shape)
    Te = np.asarray(Te_edge_K, dtype=float).reshape(q.shape)
    if np.any(~np.isfinite(q)) or np.any(~np.isfinite(delta2)) or np.any(~np.isfinite(Te)):
        raise FloatingPointError("stiffness interpolation coordinates must be finite")
    if np.any(delta2 < 0.0):
        raise ValueError("delta2_edge_J2 must be nonnegative")

    t0, t1, wt = _bracket(Te_axis, Te)
    d0, d1, wd = _bracket(delta2_axis, delta2)
    q0, q1, wq = _bracket(q_axis, q)
    out = _trilinear(table, t0, t1, wt, d0, d1, wd, q0, q1, wq)
    if np.any(~np.isfinite(out)) or np.any(out <= 0.0):
        raise FloatingPointError("stiffness interpolation produced invalid values")
    return out


def gauge_invariant_edge_pair_flow_J2_m_inv(
    *,
    delta_node_J: np.ndarray,
    ops: FVOperators,
    edge_link_variable: np.ndarray | None = None,
) -> np.ndarray:
    """Return ``Im(conj(Delta_i) U_ij Delta_j)/ell_ij`` on oriented edges."""

    delta = np.asarray(delta_node_J, dtype=np.complex128).reshape(-1)
    if delta.shape != (ops.n_nodes,):
        raise ValueError(f"delta_node_J must have shape ({ops.n_nodes},), got {delta.shape}.")
    if np.any(~np.isfinite(delta)):
        raise FloatingPointError("delta_node_J must be finite")
    link = _edge_link_variable(edge_link_variable, ops)
    pair = np.conjugate(delta[ops.edge_i]) * link * delta[ops.edge_j]
    flow = np.imag(pair) / np.asarray(ops.edge_length_m, dtype=float)
    exact_zero = (np.abs(delta[ops.edge_i]) == 0.0) | (np.abs(delta[ops.edge_j]) == 0.0)
    if np.any(exact_zero):
        flow = np.asarray(flow, dtype=float)
        flow[exact_zero] = 0.0
    if np.any(~np.isfinite(flow)):
        raise FloatingPointError("gauge-invariant edge pair flow produced non-finite values")
    return np.asarray(flow, dtype=float)


def gauge_invariant_edge_q_m_inv(
    *,
    psi_dimensionless: np.ndarray,
    ops: FVOperators,
    edge_link_variable: np.ndarray | None = None,
) -> np.ndarray:
    """Return the wrapped gauge-invariant phase gradient used as table coordinate."""

    psi = np.asarray(psi_dimensionless, dtype=np.complex128).reshape(-1)
    if psi.shape != (ops.n_nodes,):
        raise ValueError(f"psi_dimensionless must have shape ({ops.n_nodes},), got {psi.shape}.")
    if np.any(~np.isfinite(psi)):
        raise FloatingPointError("psi_dimensionless must be finite")
    link = _edge_link_variable(edge_link_variable, ops)
    pair = np.conjugate(psi[ops.edge_i]) * link * psi[ops.edge_j]
    amplitude_product = np.abs(psi[ops.edge_i]) * np.abs(psi[ops.edge_j])
    q = np.zeros(ops.n_edges, dtype=float)
    resolved = amplitude_product > np.finfo(float).tiny
    q[resolved] = np.angle(pair[resolved]) / np.asarray(ops.edge_length_m, dtype=float)[resolved]
    if np.any(~np.isfinite(q)):
        raise FloatingPointError("gauge-invariant edge q coordinate produced non-finite values")
    return q


def _edge_link_variable(values: np.ndarray | None, ops: FVOperators) -> np.ndarray:
    if values is None:
        return np.ones(ops.n_edges, dtype=np.complex128)
    link = np.asarray(values, dtype=np.complex128).reshape(-1)
    if link.shape != (ops.n_edges,):
        raise ValueError(f"edge_link_variable must have shape ({ops.n_edges},), got {link.shape}.")
    if np.any(~np.isfinite(link)):
        raise FloatingPointError("edge_link_variable must be finite")
    if not np.allclose(np.abs(link), 1.0, rtol=0.0, atol=64.0 * np.finfo(float).eps):
        raise ValueError("edge_link_variable entries must have unit modulus")
    return link


def interpolate_strict_usadel_current_table(
    *,
    table: np.ndarray | None,
    Te_axis_K: np.ndarray | None,
    delta_axis_J: np.ndarray | None,
    q_axis_m_inv: np.ndarray | None,
    q_edge_m_inv: np.ndarray,
    delta_edge_J: np.ndarray,
    Te_edge_K: np.ndarray,
) -> np.ndarray:
    """Diagnostic interpolation of the superseded ``js[Te, delta, q]`` law.

    D1 intentionally uses this routine as the reference failure mechanism.  No
    SS or photon runtime path calls it.
    """

    if table is None:
        raise ValueError("js_A_m2 table not found")
    Te_axis = _clean_axis(Te_axis_K, name="Te_axis_K", positive=True)
    delta_axis = _clean_axis(delta_axis_J, name="delta_axis_J", nonnegative=True)
    q_axis = _clean_axis(q_axis_m_inv, name="q_axis_m_inv", nonnegative=True)
    table = np.asarray(table, dtype=float)
    expected = (Te_axis.size, delta_axis.size, q_axis.size)
    if table.ndim != 3 or table.shape != expected:
        raise ValueError(f"expected js_A_m2[Te,delta,q] shape {expected}, got {table.shape}")

    q = np.asarray(q_edge_m_inv, dtype=float).reshape(-1)
    delta = np.asarray(delta_edge_J, dtype=float).reshape(q.shape)
    Te = np.asarray(Te_edge_K, dtype=float).reshape(q.shape)
    sign = np.sign(q)
    q_abs = np.abs(q)

    t0, t1, wt = _bracket(Te_axis, Te)
    d0, d1, wd = _bracket(delta_axis, delta)
    q0, q1, wq = _bracket(q_axis, q_abs)

    out = _trilinear(table, t0, t1, wt, d0, d1, wd, q0, q1, wq) * sign
    if np.any(~np.isfinite(out)):
        raise FloatingPointError("diagnostic current interpolation produced non-finite values")
    return out


def _trilinear(
    table: np.ndarray,
    t0: np.ndarray,
    t1: np.ndarray,
    wt: np.ndarray,
    d0: np.ndarray,
    d1: np.ndarray,
    wd: np.ndarray,
    q0: np.ndarray,
    q1: np.ndarray,
    wq: np.ndarray,
) -> np.ndarray:
    c000 = table[t0, d0, q0]
    c001 = table[t0, d0, q1]
    c010 = table[t0, d1, q0]
    c011 = table[t0, d1, q1]
    c100 = table[t1, d0, q0]
    c101 = table[t1, d0, q1]
    c110 = table[t1, d1, q0]
    c111 = table[t1, d1, q1]
    c00 = c000 * (1.0 - wq) + c001 * wq
    c01 = c010 * (1.0 - wq) + c011 * wq
    c10 = c100 * (1.0 - wq) + c101 * wq
    c11 = c110 * (1.0 - wq) + c111 * wq
    c0 = c00 * (1.0 - wd) + c01 * wd
    c1 = c10 * (1.0 - wd) + c11 * wd
    return c0 * (1.0 - wt) + c1 * wt


def _bracket(axis: np.ndarray, points: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    axis = np.asarray(axis, dtype=float)
    x = np.asarray(points, dtype=float)
    if axis.size == 1:
        z = np.zeros_like(x, dtype=np.int64)
        w = np.zeros_like(x, dtype=float)
        return z, z, w
    xc = np.clip(x, axis[0], axis[-1])
    hi = np.searchsorted(axis, xc, side="right")
    hi = np.clip(hi, 1, axis.size - 1).astype(np.int64)
    lo = (hi - 1).astype(np.int64)
    denom = np.maximum(axis[hi] - axis[lo], 1.0e-300)
    w = np.clip((xc - axis[lo]) / denom, 0.0, 1.0)
    return lo, hi, w


def _apply_blocked_edges(edge_js: np.ndarray, blocked_edge_mask: np.ndarray | None) -> np.ndarray:
    out = np.asarray(edge_js, dtype=float).reshape(-1)
    if blocked_edge_mask is None:
        return out
    mask = np.asarray(blocked_edge_mask, dtype=bool).reshape(-1)
    if mask.size != out.size:
        raise ValueError(f"blocked_edge_mask has length {mask.size}, expected {out.size}.")
    if np.any(mask):
        out = out.copy()
        out[mask] = 0.0
    return out


def _finish_available(
    *,
    backend: str,
    edge_q: np.ndarray,
    edge_delta: np.ndarray,
    edge_delta2: np.ndarray,
    edge_Te: np.ndarray,
    edge_pair_flow: np.ndarray,
    edge_stiffness: np.ndarray,
    edge_js: np.ndarray,
    edge_js_gl: np.ndarray,
    edge_mismatch: np.ndarray,
    ops: FVOperators,
) -> UsadelSupercurrentDiagnostics:
    edge_js = np.asarray(edge_js, dtype=float).reshape(edge_q.shape)
    edge_js_gl = np.asarray(edge_js_gl, dtype=float).reshape(edge_q.shape)
    edge_mismatch = np.asarray(edge_mismatch, dtype=float).reshape(edge_q.shape)
    node_x, node_y = edge_scalar_to_node_vector_least_squares(edge_js, ops)
    div = divergence_from_edge_scalar(edge_js, ops)
    div_gl = divergence_from_edge_scalar(edge_js_gl, ops)
    mismatch = divergence_from_edge_scalar(edge_mismatch, ops)
    return UsadelSupercurrentDiagnostics(
        available=True,
        backend=backend,
        reason="ok",
        edge_q_m_inv=np.asarray(edge_q, dtype=float),
        edge_delta_J=np.asarray(edge_delta, dtype=float),
        edge_delta2_J2=np.asarray(edge_delta2, dtype=float),
        edge_Te_K=np.asarray(edge_Te, dtype=float),
        edge_pair_flow_J2_m_inv=np.asarray(edge_pair_flow, dtype=float),
        edge_stiffness_A_per_m_J2=np.asarray(edge_stiffness, dtype=float),
        edge_js_usadel_A_m2=edge_js,
        edge_js_gl_A_m2=edge_js_gl,
        edge_js_mismatch_A_m2=edge_mismatch,
        node_js_usadel_x_A_m2=node_x,
        node_js_usadel_y_A_m2=node_y,
        node_div_js_usadel_A_m3=div,
        node_div_js_gl_A_m3=div_gl,
        node_mismatch_divergence_A_m3=mismatch,
    )


def _unavailable(
    *,
    reason: str,
    edge_q: np.ndarray,
    edge_delta: np.ndarray,
    edge_delta2: np.ndarray,
    edge_Te: np.ndarray,
    edge_pair_flow: np.ndarray,
    ops: FVOperators,
) -> UsadelSupercurrentDiagnostics:
    edge_nan = np.full(ops.n_edges, np.nan, dtype=float)
    node_nan = np.full(ops.n_nodes, np.nan, dtype=float)
    return UsadelSupercurrentDiagnostics(
        available=False,
        backend="unavailable",
        reason=str(reason),
        edge_q_m_inv=np.asarray(edge_q, dtype=float),
        edge_delta_J=np.asarray(edge_delta, dtype=float),
        edge_delta2_J2=np.asarray(edge_delta2, dtype=float),
        edge_Te_K=np.asarray(edge_Te, dtype=float),
        edge_pair_flow_J2_m_inv=np.asarray(edge_pair_flow, dtype=float),
        edge_stiffness_A_per_m_J2=edge_nan.copy(),
        edge_js_usadel_A_m2=edge_nan,
        edge_js_gl_A_m2=edge_nan.copy(),
        edge_js_mismatch_A_m2=edge_nan.copy(),
        node_js_usadel_x_A_m2=node_nan.copy(),
        node_js_usadel_y_A_m2=node_nan.copy(),
        node_div_js_usadel_A_m3=node_nan.copy(),
        node_div_js_gl_A_m3=node_nan.copy(),
        node_mismatch_divergence_A_m3=node_nan.copy(),
    )


def _find_first_array(catalog: Any, names: tuple[str, ...]) -> np.ndarray | None:
    for name in names:
        value = _get(catalog, name)
        if value is None:
            continue
        try:
            arr = np.asarray(value)
        except Exception:
            continue
        if arr.dtype.kind in "OUS":
            continue
        return arr
    return None


def _get(obj: Any, name: str) -> Any:
    if obj is None:
        return None
    try:
        if hasattr(obj, name):
            return getattr(obj, name)
    except Exception:
        pass
    try:
        return obj[name]
    except Exception:
        return None


def _clean_axis(values: np.ndarray | None, *, name: str, positive: bool = False, nonnegative: bool = False) -> np.ndarray:
    if values is None:
        raise ValueError(f"{name} not found")
    arr = np.asarray(values, dtype=float).reshape(-1)
    if arr.size == 0:
        raise ValueError(f"{name} is empty")
    if np.any(~np.isfinite(arr)):
        raise ValueError(f"{name} contains non-finite values")
    if arr.size > 1 and np.any(np.diff(arr) <= 0.0):
        raise ValueError(f"{name} must be strictly increasing")
    if positive and np.any(arr <= 0.0):
        raise ValueError(f"{name} must be positive")
    if nonnegative and np.any(arr < 0.0):
        raise ValueError(f"{name} must be nonnegative")
    return arr
