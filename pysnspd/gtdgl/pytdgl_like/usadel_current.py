"""Usadel supercurrent diagnostics for the pyTDGL-like OE7 backend.

This module is deliberately diagnostic-only.  The pyTDGL-like solver keeps using
its native GL supercurrent in the evolution.  Here we reconstruct the same edge
phase gradients and amplitudes, then evaluate a supercurrent-density table or
callable provided by the PRE-run Usadel catalogue.

All quantities are SI at this boundary:

* q_edge              [m^-1]
* |Delta|_edge        [J]
* Te_edge             [K]
* j_s^Usadel,edge     [A m^-2]
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Mapping

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
    """Thin adapter that preserves PRE supercurrent-table arrays.

    ``load_usadel_catalog_npz(...)`` intentionally constructs the historical
    Usadel DOS catalogue object and may ignore newer sidecar arrays stored in
    the same ``.npz`` file.  This wrapper delegates every old attribute/key to
    that base catalogue and exposes only the numeric supercurrent-table arrays
    needed by the pyTDGL-like SS diagnostics.
    """

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


def attach_usadel_supercurrent_table_from_npz(catalog: Any, npz_path: str | bytes | "Path") -> Any:
    """Attach PRE-tabulated ``j_s^Usadel`` arrays to a loaded catalogue.

    The function is strict about physics and permissive only about catalogue
    plumbing: it never manufactures a current table.  It only exposes arrays
    that are already present in the PRE ``usadel_dos_catalog.npz``.  Known keys
    are loaded individually so unrelated object metadata in the same NPZ cannot
    trigger ``allow_pickle=False`` failures.
    """
    arrays = load_usadel_supercurrent_table_arrays_npz(npz_path)
    if not arrays:
        return catalog
    return UsadelCatalogWithSupercurrentTable(catalog, arrays)


def load_usadel_supercurrent_table_arrays_npz(npz_path: str | bytes | "Path") -> dict[str, np.ndarray]:
    """Load only numeric PRE Usadel supercurrent-table arrays from an NPZ."""
    wanted = (
        "js_A_m2",
        "j_s_A_m2",
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
        "js_table_T_K",
        "js_table_n_matsubara",
    )
    out: dict[str, np.ndarray] = {}
    with np.load(npz_path, allow_pickle=False) as data:
        keys = set(data.files)
        for key in wanted:
            if key not in keys:
                continue
            arr = np.asarray(data[key])
            if arr.dtype.kind in "OUS":
                # This should not happen for the table, but keep the SS loader
                # numeric-only and let diagnostics fail explicitly if the table
                # itself is malformed.
                continue
            out[key] = arr
    # Canonical aliases used by the interpolator.
    if "js_A_m2" not in out:
        for key in ("j_s_A_m2", "supercurrent_density_A_m2"):
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
    return out


def compute_usadel_supercurrent_diagnostic(
    *,
    usadel_catalog: Any | None,
    psi_dimensionless: np.ndarray,
    material: GTDGLMaterial,
    Te_K: np.ndarray,
    ops: FVOperators,
    blocked_edge_mask: np.ndarray | None = None,
) -> UsadelSupercurrentDiagnostics:
    """Evaluate a diagnostic Usadel supercurrent on FV edges.

    The function recycles the existing OE7 FV operators for

    * edge phase gradient ``q``;
    * edge-averaged amplitude and temperature;
    * edge-to-node vector reconstruction;
    * FV divergence.

    If ``blocked_edge_mask`` is supplied, the Usadel current is set to zero on
    those edges before node-vector reconstruction and divergence.  This is used
    for metallic normal terminals, where the order parameter is clamped to zero
    and no supercurrent should be injected through the contact.

    The evolution is untouched.  If the catalogue does not expose a usable
    supercurrent table/callable, the returned diagnostic is marked unavailable
    and its current arrays are filled with NaNs with the correct shapes.
    """

    psi = np.asarray(psi_dimensionless, dtype=np.complex128)
    Te = np.asarray(Te_K, dtype=float)
    edge_q = edge_phase_gradient_from_psi(psi, ops)
    edge_delta = edge_average(np.abs(psi) * float(material.delta0_J), ops)
    edge_Te = edge_average(Te, ops)

    unavailable = _unavailable(
        reason="Usadel catalogue was not supplied.",
        edge_q=edge_q,
        edge_delta=edge_delta,
        edge_Te=edge_Te,
        ops=ops,
    )
    if usadel_catalog is None:
        return unavailable

    callable_backend = _find_current_callable(usadel_catalog)
    if callable_backend is not None:
        name, func = callable_backend
        try:
            edge_js = _call_current_function(func, edge_q, edge_delta, edge_Te)
        except Exception as exc:  # pragma: no cover - defensive catalogue adapter
            return _unavailable(
                reason=f"Usadel current callable {name!r} failed: {exc}",
                edge_q=edge_q,
                edge_delta=edge_delta,
                edge_Te=edge_Te,
                ops=ops,
            )
        return _finish_available(
            backend=f"callable:{name}",
            edge_q=edge_q,
            edge_delta=edge_delta,
            edge_Te=edge_Te,
            edge_js=_apply_blocked_edges(edge_js, blocked_edge_mask),
            ops=ops,
        )

    table = _find_first_array(
        usadel_catalog,
        (
            "js_A_m2",
            "j_s_A_m2",
            "supercurrent_A_m2",
            "supercurrent_density_A_m2",
            "current_density_A_m2",
            "j_super_A_m2",
            "j_s_grid_A_m2",
            "js_grid_A_m2",
        ),
    )
    if table is None:
        return _unavailable(
            reason="No Usadel supercurrent-density table was found in the catalogue.",
            edge_q=edge_q,
            edge_delta=edge_delta,
            edge_Te=edge_Te,
            ops=ops,
        )

    axes = _catalog_axes(usadel_catalog)
    try:
        edge_js, backend = _interpolate_current_table(
            table=np.asarray(table, dtype=float),
            axes=axes,
            q=edge_q,
            delta=edge_delta,
            Te=edge_Te,
        )
    except Exception as exc:
        return _unavailable(
            reason=f"Could not interpolate Usadel current table: {exc}",
            edge_q=edge_q,
            edge_delta=edge_delta,
            edge_Te=edge_Te,
            ops=ops,
        )

    return _finish_available(
        backend=backend,
        edge_q=edge_q,
        edge_delta=edge_delta,
        edge_Te=edge_Te,
        edge_js=_apply_blocked_edges(edge_js, blocked_edge_mask),
        ops=ops,
    )


def _apply_blocked_edges(edge_js: np.ndarray, blocked_edge_mask: np.ndarray | None) -> np.ndarray:
    """Set selected edge currents to zero without modifying the input array."""
    out = np.asarray(edge_js, dtype=float)
    if blocked_edge_mask is None:
        return out
    mask = np.asarray(blocked_edge_mask, dtype=bool).reshape(-1)
    if mask.size != out.size:
        raise ValueError(
            f"blocked_edge_mask has length {mask.size}, expected {out.size}."
        )
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
        edge_q_m_inv=edge_q,
        edge_delta_J=edge_delta,
        edge_Te_K=edge_Te,
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


def _call_current_function(func: Callable[..., Any], q: np.ndarray, delta: np.ndarray, Te: np.ndarray) -> np.ndarray:
    """Call a catalogue current interpolator with tolerant keyword patterns."""
    call_patterns = (
        {"q_m_inv": q, "delta_J": delta, "Te_K": Te},
        {"q": q, "delta": delta, "Te": Te},
        {"q_m_inv": q, "Delta_J": delta, "T_K": Te},
        {"q": q, "Delta": delta, "T": Te},
    )
    last_exc: Exception | None = None
    for kwargs in call_patterns:
        try:
            out = func(**kwargs)
            return np.asarray(out, dtype=float)
        except TypeError as exc:
            last_exc = exc
    try:
        return np.asarray(func(q, delta, Te), dtype=float)
    except TypeError:
        return np.asarray(func(q, delta), dtype=float)
    except Exception as exc:
        raise exc from last_exc


def _find_current_callable(catalog: Any) -> tuple[str, Callable[..., Any]] | None:
    for name in (
        "interpolate_supercurrent_density_A_m2",
        "supercurrent_density_A_m2",
        "interpolate_js_A_m2",
        "js_A_m2",
        "j_s_A_m2",
    ):
        value = _get(catalog, name)
        if callable(value):
            return name, value
    return None


def _catalog_axes(catalog: Any) -> dict[str, np.ndarray]:
    return {
        "q": _find_first_array(
            catalog,
            (
                "q_axis_m_inv",
                "q_values_m_inv",
                "q_grid_m_inv",
                "q_m_inv",
                "q_axis",
                "q_values",
            ),
        ),
        "delta": _find_first_array(
            catalog,
            (
                "delta_axis_J",
                "delta_values_J",
                "Delta_axis_J",
                "Delta_values_J",
                "delta_J",
                "Delta_J",
            ),
        ),
        "Te": _find_first_array(
            catalog,
            (
                "Te_axis_K",
                "T_axis_K",
                "temperature_axis_K",
                "Te_values_K",
                "T_values_K",
            ),
        ),
    }


def _interpolate_current_table(
    *,
    table: np.ndarray,
    axes: dict[str, np.ndarray | None],
    q: np.ndarray,
    delta: np.ndarray,
    Te: np.ndarray,
) -> tuple[np.ndarray, str]:
    table = np.asarray(table, dtype=float)
    if table.ndim == 0:
        raise ValueError("current table is scalar")

    q_axis = _clean_axis(axes.get("q"))
    delta_axis = _clean_axis(axes.get("delta"))
    Te_axis = _clean_axis(axes.get("Te"))
    if q_axis is None:
        raise ValueError("q axis not found")

    if table.ndim == 1:
        if table.shape[0] != q_axis.size:
            raise ValueError(f"1D current table length {table.shape[0]} does not match q axis {q_axis.size}")
        return _interp_signed_q(q_axis, table, q), "table:q"

    # PRE-run Usadel diagnostics write js_A_m2[delta, q].  Handle this
    # canonical layout explicitly so we do not rely on ambiguous dimension
    # matching when two axes happen to have the same length.
    if table.ndim == 2 and delta_axis is not None:
        if table.shape == (delta_axis.size, q_axis.size):
            return _interp_nd_axiswise(
                table,
                [("delta", delta_axis, delta), ("q", q_axis, q)],
            ), "table:delta,q"
        if table.shape == (q_axis.size, delta_axis.size):
            return _interp_nd_axiswise(
                table,
                [("q", q_axis, q), ("delta", delta_axis, delta)],
            ), "table:q,delta"

    # Identify each table dimension by matching unique axis lengths.  This is
    # intentionally conservative: if a dimension cannot be matched, fail loudly
    # and record an unavailable diagnostic instead of pretending to know the
    # catalogue layout.
    axis_candidates: list[tuple[str, np.ndarray | None, np.ndarray]] = [
        ("Te", Te_axis, Te),
        ("delta", delta_axis, delta),
        ("q", q_axis, q),
    ]
    dim_specs: list[tuple[str, np.ndarray, np.ndarray]] = []
    used: set[str] = set()
    for dim_len in table.shape:
        matches = [(name, axis, points) for name, axis, points in axis_candidates if axis is not None and axis.size == dim_len and name not in used]
        if not matches:
            raise ValueError(f"could not match current-table dimension length {dim_len} to q/delta/Te axes")
        # Prefer exact q/delta/Te matching by declaration order above, while
        # avoiding reusing the same axis twice.
        name, axis, points = matches[-1] if len(matches) > 1 and any(m[0] == "q" for m in matches) else matches[0]
        dim_specs.append((name, np.asarray(axis, dtype=float), np.asarray(points, dtype=float)))
        used.add(name)

    values = _interp_nd_axiswise(table, dim_specs)
    backend = "table:" + ",".join(name for name, _, _ in dim_specs)
    return values, backend


def _interp_nd_axiswise(table: np.ndarray, dim_specs: list[tuple[str, np.ndarray, np.ndarray]]) -> np.ndarray:
    """Vectorized multilinear interpolation on q/delta/Te table axes.

    The first Usadel-Poisson implementation used a Python loop over edges and
    ``np.apply_along_axis`` for every interpolation point.  On a realistic FV
    mesh this made the SS run scale like

        n_steps * n_edges * Python-call-overhead.

    This implementation keeps the same clipping and odd-in-q convention, but
    evaluates all edge points at once.  It supports the layouts used in the PRE
    catalogue, especially ``js_A_m2[delta, q]``.
    """
    table = np.asarray(table, dtype=float)
    if table.ndim != len(dim_specs):
        raise ValueError(
            f"table has {table.ndim} dimensions but {len(dim_specs)} interpolation axes were supplied"
        )
    if not dim_specs:
        raise ValueError("at least one interpolation axis is required")

    n = int(np.asarray(dim_specs[0][2]).size)
    axes: list[np.ndarray] = []
    lo: list[np.ndarray] = []
    hi: list[np.ndarray] = []
    weights: list[np.ndarray] = []
    q_sign = np.ones(n, dtype=float)

    arr = table
    for axis_index, (name, axis_raw, points_raw) in enumerate(dim_specs):
        axis = np.asarray(axis_raw, dtype=float).reshape(-1)
        points = np.asarray(points_raw, dtype=float).reshape(-1)
        if points.size != n:
            raise ValueError(f"axis {name!r} points have length {points.size}, expected {n}")
        if axis.size < 1:
            raise ValueError(f"axis {name!r} is empty")
        if arr.shape[axis_index] != axis.size:
            raise ValueError(
                f"table dimension {axis_index} has length {arr.shape[axis_index]}, "
                f"but axis {name!r} has length {axis.size}"
            )

        order = np.argsort(axis)
        axis = axis[order]
        arr = np.take(arr, order, axis=axis_index)

        if name == "q" and axis[0] >= 0.0:
            q_sign = np.sign(points)
            points = np.abs(points)

        i0, i1, w = _axis_brackets_vectorized(axis, points)
        axes.append(axis)
        lo.append(i0)
        hi.append(i1)
        weights.append(w)

    out = np.zeros(n, dtype=float)
    ndim = len(dim_specs)
    for corner in range(1 << ndim):
        idx = []
        weight = np.ones(n, dtype=float)
        for dim in range(ndim):
            if (corner >> dim) & 1:
                idx.append(hi[dim])
                weight *= weights[dim]
            else:
                idx.append(lo[dim])
                weight *= 1.0 - weights[dim]
        out += weight * arr[tuple(idx)]

    return q_sign * out


def _axis_brackets_vectorized(axis: np.ndarray, points: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return lower/upper indices and interpolation weights for all points."""
    axis = np.asarray(axis, dtype=float).reshape(-1)
    points = np.asarray(points, dtype=float).reshape(-1)
    if axis.size == 1:
        z = np.zeros(points.size, dtype=np.int64)
        return z, z, np.zeros(points.size, dtype=float)
    x = np.clip(points, axis[0], axis[-1])
    hi = np.searchsorted(axis, x, side="right")
    hi = np.clip(hi, 1, axis.size - 1).astype(np.int64, copy=False)
    lo = hi - 1
    denom = axis[hi] - axis[lo]
    denom = np.where(np.abs(denom) > 0.0, denom, 1.0)
    w = (x - axis[lo]) / denom
    return lo, hi, np.clip(w, 0.0, 1.0)


def _interp_clipped(axis: np.ndarray, values: np.ndarray, x: float) -> float:
    axis = np.asarray(axis, dtype=float).reshape(-1)
    values = np.asarray(values, dtype=float).reshape(-1)
    order = np.argsort(axis)
    axis = axis[order]
    values = values[order]
    return float(np.interp(float(np.clip(x, axis[0], axis[-1])), axis, values))


def _interp_signed_q(q_axis: np.ndarray, values: np.ndarray, q: np.ndarray) -> np.ndarray:
    q_axis = np.asarray(q_axis, dtype=float).reshape(-1)
    values = np.asarray(values, dtype=float).reshape(-1)
    q = np.asarray(q, dtype=float)
    order = np.argsort(q_axis)
    q_axis = q_axis[order]
    values = values[order]
    if q_axis.size == 1:
        if q_axis[0] >= 0.0:
            return np.sign(q) * abs(float(values[0]))
        return np.full_like(q, float(values[0]), dtype=float)
    if q_axis[0] >= 0.0:
        q_abs = np.abs(q)
        interp = np.interp(np.clip(q_abs, q_axis[0], q_axis[-1]), q_axis, values)
        return np.sign(q) * np.abs(interp)
    return np.interp(np.clip(q, q_axis[0], q_axis[-1]), q_axis, values)


def _clean_axis(axis: np.ndarray | None) -> np.ndarray | None:
    if axis is None:
        return None
    arr = np.asarray(axis, dtype=float).reshape(-1)
    arr = arr[np.isfinite(arr)]
    if arr.size < 1:
        return None
    return arr


def _find_first_array(catalog: Any, names: tuple[str, ...]) -> np.ndarray | None:
    for name in names:
        value = _get(catalog, name)
        if value is not None and not callable(value):
            arr = np.asarray(value)
            if arr.size:
                return arr
    return None


def _get(obj: Any, name: str) -> Any:
    if obj is None:
        return None
    if isinstance(obj, Mapping) and name in obj:
        return obj[name]
    if hasattr(obj, name):
        return getattr(obj, name)
    # numpy.lib.npyio.NpzFile-style access
    try:
        keys = obj.files  # type: ignore[attr-defined]
    except Exception:
        keys = None
    if keys is not None and name in keys:
        return obj[name]
    return None
