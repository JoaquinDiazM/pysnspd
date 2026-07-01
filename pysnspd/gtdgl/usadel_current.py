"""Strict Usadel supercurrent interpolation for the flat gTDGL backend.

The production SS path requires a PRE table with the canonical layout

    js_A_m2[Te, |Delta|, q].

Older one-dimensional ``j_s(q)`` tables are intentionally rejected: they are not
valid near metallic contacts where |Delta| is independently suppressed by the
mesoscopic field.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from pysnspd.gtdgl.material import GTDGLMaterial
from pysnspd.gtdgl.operators import (
    FVOperators,
    divergence_from_edge_scalar,
    edge_average,
    edge_phase_gradient_from_psi,
    edge_scalar_to_node_vector_least_squares,
)


@dataclass(frozen=True)
class UsadelSupercurrentDiagnostics:
    """Edge/node diagnostic fields for Usadel supercurrent interpolation."""

    available: bool
    backend: str
    reason: str
    edge_q_m_inv: np.ndarray
    edge_delta_J: np.ndarray
    edge_Te_K: np.ndarray
    edge_js_usadel_A_m2: np.ndarray
    node_js_usadel_x_A_m2: np.ndarray
    node_js_usadel_y_A_m2: np.ndarray
    node_div_js_usadel_A_m3: np.ndarray


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
    """Attach the strict 3D PRE supercurrent table to a loaded catalogue."""

    arrays = load_usadel_supercurrent_table_arrays_npz(npz_path)
    if not arrays:
        return catalog
    return UsadelCatalogWithSupercurrentTable(catalog, arrays)


def load_usadel_supercurrent_table_arrays_npz(npz_path: str | bytes | Path) -> dict[str, np.ndarray]:
    """Load only numeric PRE Usadel supercurrent-table arrays from an NPZ."""

    wanted = (
        "js_A_m2",
        "j_s_A_m2",
        "js_T_delta_q_A_m2",
        "supercurrent_density_A_m2",
        "q_axis_m_inv",
        "q_values_m_inv",
        "q_grid_m_inv",
        "delta_axis_J",
        "delta_values_J",
        "Delta_axis_J",
        "Delta_values_J",
        "Te_axis_K",
        "T_axis_K",
        "temperature_axis_K",
        "Te_values_K",
        "T_values_K",
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

    # Canonical aliases used by the interpolator.
    if "js_A_m2" not in out:
        for key in ("j_s_A_m2", "js_T_delta_q_A_m2", "supercurrent_density_A_m2"):
            if key in out:
                out["js_A_m2"] = out[key]
                break
    if "q_axis_m_inv" not in out:
        for key in ("q_values_m_inv", "q_grid_m_inv"):
            if key in out:
                out["q_axis_m_inv"] = out[key]
                break
    if "delta_axis_J" not in out:
        for key in ("delta_values_J", "Delta_axis_J", "Delta_values_J"):
            if key in out:
                out["delta_axis_J"] = out[key]
                break
    if "Te_axis_K" not in out:
        for key in ("T_axis_K", "temperature_axis_K", "Te_values_K", "T_values_K"):
            if key in out:
                out["Te_axis_K"] = out[key]
                break
    return out


def validate_strict_usadel_supercurrent_table_npz(npz_path: str | bytes | Path) -> dict[str, Any]:
    """Validate that a PRE NPZ exposes ``js_A_m2[Te, delta, q]``.

    Returns a compact summary if valid; raises ``RuntimeError`` otherwise.
    """

    arrays = load_usadel_supercurrent_table_arrays_npz(npz_path)
    table = arrays.get("js_A_m2")
    q_axis = arrays.get("q_axis_m_inv")
    delta_axis = arrays.get("delta_axis_J")
    Te_axis = arrays.get("Te_axis_K")
    missing = [
        name
        for name, arr in (
            ("js_A_m2", table),
            ("Te_axis_K", Te_axis),
            ("delta_axis_J", delta_axis),
            ("q_axis_m_inv", q_axis),
        )
        if arr is None
    ]
    if missing:
        raise RuntimeError(
            "PRE Usadel current table is not strict 3D. Missing: "
            + ", ".join(missing)
            + ". Re-run 01_prerun_template.py with the 3D Matsubara table builder."
        )

    table = np.asarray(table, dtype=float)
    Te_axis = _clean_axis(Te_axis, name="Te_axis_K", positive=True)
    delta_axis = _clean_axis(delta_axis, name="delta_axis_J", nonnegative=True)
    q_axis = _clean_axis(q_axis, name="q_axis_m_inv", nonnegative=True)
    expected = (Te_axis.size, delta_axis.size, q_axis.size)
    if table.ndim != 3 or table.shape != expected:
        raise RuntimeError(
            "PRE Usadel current table must have layout js_A_m2[Te, delta, q] "
            f"with shape {expected}; got ndim={table.ndim}, shape={table.shape}."
        )
    if np.any(~np.isfinite(table)):
        raise RuntimeError("PRE Usadel current table contains non-finite values.")
    if Te_axis.size < 1 or delta_axis.size < 2 or q_axis.size < 2:
        raise RuntimeError(
            "PRE Usadel current table axes are too small for local interpolation: "
            f"n_Te={Te_axis.size}, n_delta={delta_axis.size}, n_q={q_axis.size}."
        )

    return {
        "valid": True,
        "layout": "Te,delta,q",
        "shape": list(table.shape),
        "n_Te": int(Te_axis.size),
        "n_delta": int(delta_axis.size),
        "n_q": int(q_axis.size),
        "Te_min_K": float(np.min(Te_axis)),
        "Te_max_K": float(np.max(Te_axis)),
        "delta_min_J": float(np.min(delta_axis)),
        "delta_max_J": float(np.max(delta_axis)),
        "q_min_m_inv": float(np.min(q_axis)),
        "q_max_m_inv": float(np.max(q_axis)),
        "backend": str(np.asarray(arrays.get("js_table_backend", np.array("unknown"))).reshape(()).item())
        if "js_table_backend" in arrays and np.asarray(arrays["js_table_backend"]).shape == ()
        else "matsubara_usadel_supercurrent_table_3d_v1",
    }


def compute_usadel_supercurrent_diagnostic(
    *,
    usadel_catalog: Any | None,
    psi_dimensionless: np.ndarray,
    material: GTDGLMaterial,
    Te_K: np.ndarray,
    ops: FVOperators,
    blocked_edge_mask: np.ndarray | None = None,
) -> UsadelSupercurrentDiagnostics:
    """Evaluate the strict 3D Usadel current table on FV edges."""

    psi = np.asarray(psi_dimensionless, dtype=np.complex128)
    Te = np.asarray(Te_K, dtype=float)
    edge_q = edge_phase_gradient_from_psi(psi, ops)
    edge_delta = edge_average(np.abs(psi) * float(material.delta0_J), ops)
    edge_Te = edge_average(Te, ops)

    if usadel_catalog is None:
        return _unavailable(
            reason="Usadel catalogue was not supplied.",
            edge_q=edge_q,
            edge_delta=edge_delta,
            edge_Te=edge_Te,
            ops=ops,
        )

    table = _find_first_array(usadel_catalog, ("js_A_m2", "j_s_A_m2", "js_T_delta_q_A_m2"))
    q_axis = _find_first_array(usadel_catalog, ("q_axis_m_inv", "q_values_m_inv", "q_grid_m_inv"))
    delta_axis = _find_first_array(usadel_catalog, ("delta_axis_J", "delta_values_J", "Delta_axis_J", "Delta_values_J"))
    Te_axis = _find_first_array(usadel_catalog, ("Te_axis_K", "T_axis_K", "temperature_axis_K", "Te_values_K", "T_values_K"))

    try:
        edge_js = interpolate_strict_usadel_current_table(
            table=np.asarray(table, dtype=float) if table is not None else None,
            Te_axis_K=Te_axis,
            delta_axis_J=delta_axis,
            q_axis_m_inv=q_axis,
            q_edge_m_inv=edge_q,
            delta_edge_J=edge_delta,
            Te_edge_K=edge_Te,
        )
    except Exception as exc:
        return _unavailable(
            reason=f"Strict 3D Usadel current table is unavailable/invalid: {exc}",
            edge_q=edge_q,
            edge_delta=edge_delta,
            edge_Te=edge_Te,
            ops=ops,
        )

    edge_js = _apply_blocked_edges(edge_js, blocked_edge_mask)
    return _finish_available(
        backend="table:Te,delta,q",
        edge_q=edge_q,
        edge_delta=edge_delta,
        edge_Te=edge_Te,
        edge_js=edge_js,
        ops=ops,
    )


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
    """Vectorized trilinear interpolation of ``js[Te, delta, q]``."""

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
    out = (c0 * (1.0 - wt) + c1 * wt) * sign
    out[~np.isfinite(out)] = 0.0
    return out


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
    edge_Te: np.ndarray,
    edge_js: np.ndarray,
    ops: FVOperators,
) -> UsadelSupercurrentDiagnostics:
    edge_js = np.asarray(edge_js, dtype=float).reshape(edge_q.shape)
    node_x, node_y = edge_scalar_to_node_vector_least_squares(edge_js, ops)
    div = divergence_from_edge_scalar(edge_js, ops)
    return UsadelSupercurrentDiagnostics(
        available=True,
        backend=backend,
        reason="ok",
        edge_q_m_inv=np.asarray(edge_q, dtype=float),
        edge_delta_J=np.asarray(edge_delta, dtype=float),
        edge_Te_K=np.asarray(edge_Te, dtype=float),
        edge_js_usadel_A_m2=edge_js,
        node_js_usadel_x_A_m2=node_x,
        node_js_usadel_y_A_m2=node_y,
        node_div_js_usadel_A_m3=div,
    )


def _unavailable(
    *,
    reason: str,
    edge_q: np.ndarray,
    edge_delta: np.ndarray,
    edge_Te: np.ndarray,
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
        edge_Te_K=np.asarray(edge_Te, dtype=float),
        edge_js_usadel_A_m2=edge_nan,
        node_js_usadel_x_A_m2=node_nan.copy(),
        node_js_usadel_y_A_m2=node_nan.copy(),
        node_div_js_usadel_A_m3=node_nan.copy(),
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
    if np.any(np.diff(arr) < 0.0):
        arr = np.sort(arr)
    if positive and np.any(arr <= 0.0):
        raise ValueError(f"{name} must be positive")
    if nonnegative and np.any(arr < 0.0):
        raise ValueError(f"{name} must be nonnegative")
    return arr
