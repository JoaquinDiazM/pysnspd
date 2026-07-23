"""Additional memory-quality stationary SS plots.

This module is plotting-only. It complements the existing
``02_plot_ss_run.py`` pipeline without replacing any current figure.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

import numpy as np

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.ticker as mticker
import matplotlib.tri as mtri

from pysnspd.plotting.ss_power_helpers import _snapshot_diffusion_power_density


def _prepare_output(output_path: str | Path) -> Path:
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    return output


def _load_npz_if_exists(path: Path) -> dict[str, np.ndarray]:
    if not path.exists():
        return {}
    with np.load(path, allow_pickle=True) as data:
        return {key: np.asarray(data[key]) for key in data.files}


def _triangulation(mesh: Any, dataset: Mapping[str, Any]) -> mtri.Triangulation:
    x = np.asarray(dataset.get("x_nm", np.asarray(mesh.nodes)[:, 0] * 1.0e9), dtype=float)
    y = np.asarray(dataset.get("y_nm", np.asarray(mesh.nodes)[:, 1] * 1.0e9), dtype=float)
    if "triangles" in dataset:
        triangles = np.asarray(dataset["triangles"], dtype=np.int64)
    elif hasattr(mesh, "triangles"):
        triangles = np.asarray(mesh.triangles, dtype=np.int64)
    elif hasattr(mesh, "elements"):
        triangles = np.asarray(mesh.elements, dtype=np.int64)
    else:
        raise ValueError("cannot build triangulation: missing triangles/elements")
    return mtri.Triangulation(x, y, triangles)


def _javg(dataset: Mapping[str, Any], snapshots: Mapping[str, np.ndarray] | None = None) -> float:
    for source in (dataset, snapshots or {}):
        if "javg_A_m2" in source:
            arr = np.asarray(source["javg_A_m2"], dtype=float).reshape(-1)
            if arr.size and np.isfinite(arr[-1]) and abs(arr[-1]) > 0.0:
                return abs(float(arr[-1]))
    return 1.0


def _snapshot_array(
    data: Mapping[str, np.ndarray] | None,
    keys: str | tuple[str, ...],
    *,
    shape_like: np.ndarray | None = None,
) -> np.ndarray:
    if isinstance(keys, str):
        keys = (keys,)
    if data:
        for key in keys:
            if key in data:
                arr = np.asarray(data[key], dtype=float)
                if arr.ndim == 1:
                    arr = arr[None, :]
                return arr
    if shape_like is not None and np.asarray(shape_like).size:
        return np.zeros_like(np.asarray(shape_like, dtype=float))
    return np.empty((0, 0), dtype=float)


def _snapshot_times_ps(data: Mapping[str, np.ndarray], *, preferred: tuple[str, ...], n: int) -> np.ndarray:
    if n <= 0:
        return np.array([], dtype=float)
    for key in preferred:
        if key in data:
            arr = np.asarray(data[key], dtype=float).reshape(-1)
            if arr.size:
                if arr.size != n:
                    arr = np.resize(arr, n)
                return arr / 1.0e-12
    return np.arange(int(n), dtype=float)


def _dataset_1d_array(dataset: Mapping[str, Any] | None, keys: tuple[str, ...]) -> np.ndarray:
    if dataset:
        for key in keys:
            if key in dataset:
                arr = np.asarray(dataset[key], dtype=float).reshape(-1)
                if arr.size:
                    return arr
    return np.empty(0, dtype=float)


def _last_snapshot_row(arr: np.ndarray) -> np.ndarray:
    a = np.asarray(arr, dtype=float)
    if a.size == 0:
        return np.empty(0, dtype=float)
    if a.ndim == 1:
        return a.reshape(-1)
    return a[-1].reshape(-1)


def _final_current_family_arrays(
    dataset: Mapping[str, Any],
    snapshots: Mapping[str, np.ndarray] | None,
    *,
    family: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if family == "jn":
        dataset_mag_keys = ("jn_mag_A_m2", "normal_current_density_A_m2")
        dataset_x_keys = ("jn_x_A_m2", "normal_current_density_x_A_m2")
        dataset_y_keys = ("jn_y_A_m2", "normal_current_density_y_A_m2")
        snapshot_mag_keys = ("jn_snapshot_mag_A_m2", "normal_current_density_snapshot_A_m2")
        snapshot_x_keys = ("jn_snapshot_x_A_m2", "normal_current_density_snapshot_x_A_m2")
        snapshot_y_keys = ("jn_snapshot_y_A_m2", "normal_current_density_snapshot_y_A_m2")
    elif family == "js_us":
        dataset_mag_keys = ("js_us_mag_A_m2", "js_usadel_mag_A_m2", "supercurrent_Usadel_density_A_m2")
        dataset_x_keys = (
            "js_us_x_A_m2",
            "js_usadel_x_A_m2",
            "supercurrent_Usadel_x_A_m2",
            "supercurrent_Usadel_density_x_A_m2",
        )
        dataset_y_keys = (
            "js_us_y_A_m2",
            "js_usadel_y_A_m2",
            "supercurrent_Usadel_y_A_m2",
            "supercurrent_Usadel_density_y_A_m2",
        )
        snapshot_mag_keys = ("supercurrent_Usadel_density_snapshot_A_m2", "js_usadel_snapshot_mag_A_m2")
        snapshot_x_keys = ("supercurrent_Usadel_density_snapshot_x_A_m2", "js_usadel_snapshot_x_A_m2")
        snapshot_y_keys = ("supercurrent_Usadel_density_snapshot_y_A_m2", "js_usadel_snapshot_y_A_m2")
    else:
        raise ValueError(f"unknown current family: {family}")

    mag = _dataset_1d_array(dataset, dataset_mag_keys)
    jx = _dataset_1d_array(dataset, dataset_x_keys)
    jy = _dataset_1d_array(dataset, dataset_y_keys)

    if snapshots:
        if mag.size == 0:
            mag = _last_snapshot_row(_snapshot_array(snapshots, snapshot_mag_keys))
        if jx.size == 0:
            jx = _last_snapshot_row(_snapshot_array(snapshots, snapshot_x_keys))
        if jy.size == 0:
            jy = _last_snapshot_row(_snapshot_array(snapshots, snapshot_y_keys))

    if mag.size == 0 and jx.size and jy.size:
        mag = np.sqrt(jx * jx + jy * jy)

    return mag, jx, jy


def _center_window(
    x_nm: np.ndarray,
    y_nm: np.ndarray,
    *,
    center_width_nm: float,
) -> tuple[tuple[float, float], tuple[float, float], np.ndarray]:
    xmid = 0.5 * (float(np.nanmin(x_nm)) + float(np.nanmax(x_nm)))
    half_width = 0.5 * float(center_width_nm)
    xlim = (xmid - half_width, xmid + half_width)
    ylim = (float(np.nanmin(y_nm)), float(np.nanmax(y_nm)))
    crop_mask = (x_nm >= xlim[0]) & (x_nm <= xlim[1])
    if not np.any(crop_mask):
        raise ValueError("central-strip plot has no nodes inside the requested x window")
    return xlim, ylim, crop_mask


def _center_mask_from_x(x_nm: np.ndarray, *, center_width_nm: float) -> np.ndarray:
    if x_nm.size == 0:
        return np.array([], dtype=bool)
    xmid = 0.5 * (float(np.nanmin(x_nm)) + float(np.nanmax(x_nm)))
    half_width = 0.5 * float(center_width_nm)
    return (x_nm >= xmid - half_width) & (x_nm <= xmid + half_width)


def _has_current_family_data(
    dataset: Mapping[str, Any],
    snapshots: Mapping[str, np.ndarray] | None,
    *,
    family: str,
) -> bool:
    mag, jx, jy = _final_current_family_arrays(dataset, snapshots, family=family)
    return bool((mag.size or (jx.size and jy.size)) and jx.size and jy.size)


def _has_thermal_scalar_data(
    power: Mapping[str, np.ndarray] | None,
    snapshots: Mapping[str, np.ndarray] | None,
    dataset: Mapping[str, Any],
) -> bool:
    return bool(
        _snapshot_array(power, "P_total_snapshot_W_m3").size
        or _snapshot_array(snapshots, "Te_snapshot_K").size
        or _snapshot_array(snapshots, "Tph_snapshot_K").size
        or _snapshot_array(power, "u_e_snapshot_J_m3").size
        or _snapshot_array(power, "u_ph_snapshot_J_m3").size
        or _has_nonempty(dataset, "thermal_mean_Te_K_history", "thermal_mean_P_J_W_m3_history", "thermal_mean_u_e_J_m3_history")
    )


def _has_field_map_data(dataset: Mapping[str, Any]) -> bool:
    return bool(np.asarray(dataset.get("delta_over_delta0", []), dtype=float).size)


def _snapshot_pair(
    t_ps: np.ndarray,
    values: np.ndarray,
    center_mask: np.ndarray,
    *,
    max_mode: str,
    abs_mean: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        return np.array([], dtype=float), np.array([], dtype=float), np.array([], dtype=float)
    if arr.ndim == 1:
        arr = arr[None, :]
    if center_mask.size != arr.shape[1]:
        center_mask = np.resize(center_mask.astype(bool), arr.shape[1])
    if not np.any(center_mask):
        center_mask = np.ones(arr.shape[1], dtype=bool)

    sub = arr[:, center_mask]
    with np.errstate(invalid="ignore"):
        if abs_mean:
            mean_vals = np.nanmean(np.abs(sub), axis=1)
        else:
            mean_vals = np.nanmean(sub, axis=1)

        if max_mode == "max_abs":
            max_vals = np.nanmax(np.abs(sub), axis=1)
        elif max_mode == "max":
            max_vals = np.nanmax(sub, axis=1)
        else:
            raise ValueError(f"unknown max_mode: {max_mode}")

    t = np.asarray(t_ps, dtype=float).reshape(-1)
    if t.size != arr.shape[0]:
        t = np.resize(t, arr.shape[0])
    return t, mean_vals, max_vals


def _snapshot_mean_series(
    t_ps: np.ndarray,
    values: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Spatial mean over the full snapshot map.

    This mirrors the energy extraction used by
    ``plot_ss_snapshot_runtime_metrics`` in its stored-energy panel.
    """

    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        return np.array([], dtype=float), np.array([], dtype=float)
    if arr.ndim == 1:
        arr = arr[None, :]

    with np.errstate(invalid="ignore"):
        mean_vals = np.nanmean(arr, axis=1)

    t = np.asarray(t_ps, dtype=float).reshape(-1)
    if t.size != arr.shape[0]:
        t = np.resize(t, arr.shape[0])
    return t, mean_vals


def _history_or_snapshot_pair(
    dataset: Mapping[str, Any],
    t_hist: np.ndarray,
    *,
    mean_key: str,
    max_key: str,
    snapshot_t: np.ndarray,
    snapshot_values: np.ndarray,
    center_mask: np.ndarray,
    max_mode: str = "max",
    abs_snapshot: bool = False,
    smooth_history: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    t, mean_vals, max_vals = _history_pair_if_compatible(dataset, t_hist, mean_key, max_key)
    if t.size:
        if smooth_history:
            t_mean, mean_vals = _filter_history_series(t, mean_vals)
            t_max, max_vals = _filter_history_series(t, max_vals)
            if t_mean.size != t_max.size:
                max_vals = np.resize(max_vals, t_mean.size)
            return t_mean, mean_vals, max_vals
        return t, mean_vals, max_vals

    return _snapshot_pair(
        snapshot_t,
        snapshot_values,
        center_mask,
        max_mode=max_mode,
        abs_mean=abs_snapshot,
    )


def _history_pair_if_compatible(
    dataset: Mapping[str, Any],
    t_hist: np.ndarray,
    mean_key: str,
    max_key: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    t = np.asarray(t_hist, dtype=float).reshape(-1)
    mean_vals = np.asarray(dataset.get(mean_key, []), dtype=float).reshape(-1)
    max_vals = np.asarray(dataset.get(max_key, []), dtype=float).reshape(-1)

    if t.size == 0 or max_vals.size == 0:
        return np.array([], dtype=float), np.array([], dtype=float), np.array([], dtype=float)

    if mean_vals.size == 0:
        mean_vals = max_vals.copy()

    n = min(t.size, mean_vals.size, max_vals.size)
    if n <= 1:
        return np.array([], dtype=float), np.array([], dtype=float), np.array([], dtype=float)

    return t[:n], mean_vals[:n], max_vals[:n]


def _overlay_current_arrows(
    ax,
    x_nm: np.ndarray,
    y_nm: np.ndarray,
    jx: np.ndarray,
    jy: np.ndarray,
    mask: np.ndarray,
    *,
    arrow_length_nm: float,
) -> None:
    jx = np.asarray(jx, dtype=float)
    jy = np.asarray(jy, dtype=float)
    mag = np.sqrt(jx * jx + jy * jy)
    keep = np.asarray(mask, dtype=bool) & np.isfinite(jx) & np.isfinite(jy) & np.isfinite(mag)
    keep &= mag > 1.0e-12 * max(float(np.nanmax(mag)) if np.any(np.isfinite(mag)) else 0.0, 1.0)
    idx = np.flatnonzero(keep)
    if idx.size == 0:
        return
    step = max(1, int(np.ceil(np.sqrt(idx.size / 140.0))))
    idx = idx[::step]
    mag_sel = mag[idx]
    ux = arrow_length_nm * jx[idx] / np.maximum(mag_sel, 1.0e-300)
    uy = arrow_length_nm * jy[idx] / np.maximum(mag_sel, 1.0e-300)
    ax.quiver(
        x_nm[idx],
        y_nm[idx],
        ux,
        uy,
        angles="xy",
        scale_units="xy",
        scale=1.0,
        color="black",
        width=0.0042,
        headwidth=3.8,
        headlength=4.8,
        headaxislength=4.1,
        pivot="middle",
        alpha=0.88,
        zorder=4,
    )


def _plot_series_pair(
    ax,
    t_ps: np.ndarray,
    mean_vals: np.ndarray,
    max_vals: np.ndarray,
    *,
    color: str,
    take_abs: bool = False,
):
    t = np.asarray(t_ps, dtype=float).reshape(-1)
    mean_arr = np.asarray(mean_vals, dtype=float).reshape(-1)
    max_arr = np.asarray(max_vals, dtype=float).reshape(-1)

    if mean_arr.size == 0 and max_arr.size == 0:
        h1, = ax.plot([], [], color=color, linestyle="-")
        h2, = ax.plot([], [], color=color, linestyle="--")
        return h1, h2

    if mean_arr.size == 0:
        mean_arr = max_arr.copy()
    if max_arr.size == 0:
        max_arr = mean_arr.copy()

    n = min(t.size if t.size else max(mean_arr.size, max_arr.size), mean_arr.size, max_arr.size)
    if n <= 1:
        h1, = ax.plot(t[:n], mean_arr[:n], color=color, linestyle="-")
        h2, = ax.plot(t[:n], max_arr[:n], color=color, linestyle="--")
        return h1, h2

    t = t[:n]
    mean_arr = mean_arr[:n]
    max_arr = max_arr[:n]

    if take_abs:
        mean_arr = np.abs(mean_arr)
        max_arr = np.abs(max_arr)

    t_mean, y_mean = _smooth_series(t, mean_arr)
    t_max, y_max = _smooth_series(t, max_arr)

    h_mean, = ax.plot(t_mean, y_mean, color=color, linestyle="-")
    h_max, = ax.plot(t_max, y_max, color=color, linestyle="--")
    return h_mean, h_max


def _plot_single_series(
    ax,
    t_ps: np.ndarray,
    values: np.ndarray,
    *,
    color: str,
):
    t = np.asarray(t_ps, dtype=float).reshape(-1)
    y = np.asarray(values, dtype=float).reshape(-1)
    n = min(t.size, y.size)
    if n <= 1:
        h, = ax.plot(t[:n], y[:n], color=color, linestyle="-")
        return h
    t_smooth, y_smooth = _smooth_series(t[:n], y[:n])
    h, = ax.plot(t_smooth, y_smooth, color=color, linestyle="-")
    return h


def _smooth_series(t: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    t = np.asarray(t, dtype=float).reshape(-1)
    y = np.asarray(y, dtype=float).reshape(-1)
    n = min(t.size, y.size)
    if n <= 1:
        return t[:n], y[:n]
    t = t[:n]
    y = y[:n]

    mask = np.isfinite(t) & np.isfinite(y)
    if np.count_nonzero(mask) < 2:
        return t[mask], y[mask]
    t = t[mask]
    y = y[mask]

    order = np.argsort(t)
    t = t[order]
    y = y[order]
    t, unique_idx = np.unique(t, return_index=True)
    y = y[unique_idx]

    if t.size < 3:
        return t, y
    if t.size < 12:
        return _monotone_cubic_smooth(t, y, n_dense=320)
    return _filter_history_series(t, y)


def _monotone_cubic_smooth(t: np.ndarray, y: np.ndarray, *, n_dense: int = 300) -> tuple[np.ndarray, np.ndarray]:
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    n = t.size
    if n < 2:
        return t, y

    h = np.diff(t)
    delta = np.diff(y) / np.maximum(h, 1.0e-300)
    m = np.zeros(n, dtype=float)
    m[0] = delta[0]
    m[-1] = delta[-1]
    for k in range(1, n - 1):
        if delta[k - 1] * delta[k] <= 0.0:
            m[k] = 0.0
        else:
            w1 = 2.0 * h[k] + h[k - 1]
            w2 = h[k] + 2.0 * h[k - 1]
            m[k] = (w1 + w2) / (w1 / delta[k - 1] + w2 / delta[k])

    td = np.linspace(float(t[0]), float(t[-1]), int(n_dense))
    yd = np.empty_like(td)

    j = 0
    for i, tv in enumerate(td):
        while j < n - 2 and tv > t[j + 1]:
            j += 1
        hj = t[j + 1] - t[j]
        if hj <= 0.0:
            yd[i] = y[j]
            continue
        s = (tv - t[j]) / hj
        h00 = 2.0 * s**3 - 3.0 * s**2 + 1.0
        h10 = s**3 - 2.0 * s**2 + s
        h01 = -2.0 * s**3 + 3.0 * s**2
        h11 = s**3 - s**2
        yd[i] = h00 * y[j] + h10 * hj * m[j] + h01 * y[j + 1] + h11 * hj * m[j + 1]
    return td, yd


def _filter_history_series(t: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    t = np.asarray(t, dtype=float).reshape(-1)
    y = np.asarray(y, dtype=float).reshape(-1)
    n = min(t.size, y.size)
    if n <= 1:
        return t[:n], y[:n]
    t = t[:n]
    y = y[:n]

    mask = np.isfinite(t) & np.isfinite(y)
    if np.count_nonzero(mask) < 2:
        return t[mask], y[mask]
    t = t[mask]
    y = y[mask]

    order = np.argsort(t)
    t = t[order]
    y = y[order]

    n = t.size
    if n < 7:
        return _monotone_cubic_smooth(t, y, n_dense=280)

    win = max(5, int(n // 180))
    if win % 2 == 0:
        win += 1
    kernel = np.ones(win, dtype=float) / float(win)
    y_pad = np.pad(y, (win // 2, win // 2), mode="edge")
    y_smooth = np.convolve(y_pad, kernel, mode="valid")

    target = 900
    if n > target:
        idx = np.linspace(0, n - 1, target).astype(int)
        return t[idx], y_smooth[idx]
    return t, y_smooth


def _shade_prethermal(ax, end_ps: float) -> None:
    ax.axvspan(0.0, float(end_ps), color="0.85", alpha=0.65, zorder=0)


def _clean_twin_axis(ax, ax_r) -> None:
    ax.spines["right"].set_visible(False)
    ax_r.spines["left"].set_visible(False)
    ax_r.patch.set_alpha(0.0)


def _apply_axis_limits(ax, arrays: list[np.ndarray], *, frac: float = 0.05) -> None:
    values = []
    for arr in arrays:
        a = np.asarray(arr, dtype=float).reshape(-1)
        a = a[np.isfinite(a)]
        if a.size:
            values.append(a)
    if not values:
        return
    vals = np.concatenate(values)
    vmin = float(np.nanmin(vals))
    vmax = float(np.nanmax(vals))
    if np.isclose(vmin, vmax):
        pad = max(0.01, 0.02 * max(abs(vmin), 1.0))
    else:
        pad = float(frac) * (vmax - vmin)
    ax.set_ylim(vmin - pad, vmax + pad)


def _center_linear_axis(ax, arrays: list[np.ndarray], *, frac: float = 0.10) -> None:
    values = []
    for arr in arrays:
        a = np.asarray(arr, dtype=float).reshape(-1)
        a = a[np.isfinite(a)]
        if a.size:
            values.append(a)
    if not values:
        return

    vals = np.concatenate(values)
    vmin = float(np.nanmin(vals))
    vmax = float(np.nanmax(vals))

    if np.isclose(vmin, vmax):
        center = 0.5 * (vmin + vmax)
        pad = max(1.0e-12, float(frac) * max(abs(center), 1.0))
        ax.set_ylim(center - pad, center + pad)
        return

    pad = float(frac) * (vmax - vmin)
    ax.set_ylim(vmin - pad, vmax + pad)


def _linthresh_from_axes_lines(ax) -> float:
    vals = []
    for line in ax.lines:
        y = np.asarray(line.get_ydata(), dtype=float)
        y = y[np.isfinite(y)]
        y = np.abs(y[y != 0.0])
        if y.size:
            vals.append(y)
    if not vals:
        return 1.0
    concat = np.concatenate(vals)
    return max(float(np.nanpercentile(concat, 12.0)), 1.0e-30)


def _autoscale_symlog_axis(ax) -> None:
    vals = []
    for line in ax.lines:
        y = np.asarray(line.get_ydata(), dtype=float)
        y = y[np.isfinite(y)]
        if y.size:
            vals.append(y)
    if not vals:
        return
    finite = np.concatenate(vals)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return
    ymin = float(np.nanmin(finite))
    ymax = float(np.nanmax(finite))
    if np.isclose(ymin, ymax):
        return
    pad = 0.10 * (ymax - ymin)
    ax.set_ylim(ymin - pad, ymax + pad)


def _has_nonempty(dataset: Mapping[str, Any], *keys: str) -> bool:
    for key in keys:
        arr = np.asarray(dataset.get(key, []), dtype=float).reshape(-1)
        if arr.size:
            return True
    return False
