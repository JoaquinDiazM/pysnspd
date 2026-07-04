"""Additional memory-quality stationary SS plots.

This module is plotting-only. It complements the existing
``02_plot_ss_run.py`` pipeline without replacing any current figure.

Added figures
-------------
1. ``ss_final_center_current_maps.png``
   Three final-state central-strip colormaps (100 nm centrales in x):
   ``j_s^Us``, ``j_n`` and ``j_s^Us + j_n``.  All share one colorbar and are
   normalized to ``j_avg``.  Current arrows are overlaid.

2. ``ss_snapshot_thermal_scalars.png``
   Three horizontal time-series panels: temperatures, powers and energies.
   Spatial means are solid and spatial maxima are dashed.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

import numpy as np

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import matplotlib.tri as mtri


def make_ss_memory_figures(
    *,
    mesh: Any,
    dataset: Mapping[str, Any],
    raw_ss: str | Path,
    output_dir: str | Path,
    dpi: int = 480,
    center_width_nm: float = 100.0,
) -> dict[str, Path]:
    """Create additional memory-quality SS figures.

    Missing snapshot diagnostic files are treated as non-errors so that the
    normal SS plotting pipeline remains backward compatible with older runs.
    """

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    raw = Path(raw_ss)
    snapshots = _load_npz_if_exists(raw / "stationary_snapshots.npz")
    power = _load_npz_if_exists(raw / "snapshot_power_energy_diagnostics.npz")

    saved: dict[str, Path] = {}

    if _has_current_family_data(dataset, snapshots, family="js_us") and _has_current_family_data(
        dataset,
        snapshots,
        family="jn",
    ):
        saved["final_center_current_maps"] = plot_ss_final_center_current_maps(
            mesh,
            dataset,
            snapshots,
            out / "ss_final_center_current_maps.png",
            center_width_nm=center_width_nm,
            dpi=dpi,
        )

    if _has_thermal_scalar_snapshot_data(power, snapshots):
        saved["snapshot_thermal_scalars"] = plot_ss_snapshot_thermal_scalars(
            power,
            snapshots,
            dataset,
            out / "ss_snapshot_thermal_scalars.png",
            dpi=dpi,
        )

    return saved


def plot_ss_final_center_current_maps(
    mesh: Any,
    dataset: Mapping[str, Any],
    snapshots: Mapping[str, np.ndarray] | None,
    output_path: str | Path,
    *,
    center_width_nm: float = 100.0,
    dpi: int = 480,
) -> Path:
    output = _prepare_output(output_path)
    tri = _triangulation(mesh, dataset)
    x_nm = np.asarray(tri.x, dtype=float)
    y_nm = np.asarray(tri.y, dtype=float)
    jscale = _javg(dataset, snapshots)

    js_mag, js_x, js_y = _final_current_family_arrays(dataset, snapshots, family="js_us")
    jn_mag, jn_x, jn_y = _final_current_family_arrays(dataset, snapshots, family="jn")

    n_nodes = x_nm.size
    js_mag = np.resize(js_mag, n_nodes)
    js_x = np.resize(js_x, n_nodes)
    js_y = np.resize(js_y, n_nodes)
    jn_mag = np.resize(jn_mag, n_nodes)
    jn_x = np.resize(jn_x, n_nodes)
    jn_y = np.resize(jn_y, n_nodes)

    jtot_x = js_x + jn_x
    jtot_y = js_y + jn_y
    jtot_mag = np.sqrt(jtot_x * jtot_x + jtot_y * jtot_y)

    xmid = 0.5 * (float(np.nanmin(x_nm)) + float(np.nanmax(x_nm)))
    half_width = 0.5 * float(center_width_nm)
    xlim = (xmid - half_width, xmid + half_width)
    ylim = (float(np.nanmin(y_nm)), float(np.nanmax(y_nm)))
    crop_mask = (x_nm >= xlim[0]) & (x_nm <= xlim[1])
    if not np.any(crop_mask):
        raise ValueError("central-strip current map has no nodes inside the requested x window")

    families = [
        (r"$j_s^{Us}$", js_mag / jscale, js_x / jscale, js_y / jscale),
        (r"$j_n$", jn_mag / jscale, jn_x / jscale, jn_y / jscale),
        (r"$j_s^{Us}+j_n$", jtot_mag / jscale, jtot_x / jscale, jtot_y / jscale),
    ]

    finite_for_scale: list[np.ndarray] = []
    for _, mag, _, _ in families:
        vis = np.asarray(mag[crop_mask], dtype=float)
        vis = vis[np.isfinite(vis)]
        if vis.size:
            finite_for_scale.append(vis)
    if not finite_for_scale:
        raise ValueError("central-strip current map does not contain finite current values")
    vmax = max(1.0, float(np.nanmax(np.concatenate(finite_for_scale))))

    fig, axes = plt.subplots(1, 3, figsize=(13.2, 3.45), constrained_layout=False)
    fig.subplots_adjust(left=0.060, right=0.920, bottom=0.150, top=0.860, wspace=0.18)
    fig.suptitle(f"SS final-state center-strip currents: {dataset.get('run_name', '')}", y=0.955)

    mappable = None
    arrow_length_nm = 0.10 * max(ylim[1] - ylim[0], 1.0)
    for idx, (ax, (title, mag, jx, jy)) in enumerate(zip(axes, families)):
        z = np.asarray(mag, dtype=float)
        mappable = ax.tripcolor(tri, z, shading="gouraud", vmin=0.0, vmax=vmax)
        _overlay_current_arrows(
            ax,
            x_nm,
            y_nm,
            jx,
            jy,
            crop_mask,
            arrow_length_nm=arrow_length_nm,
        )
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(title)
        ax.set_xlabel("x [nm]")
        if idx == 0:
            ax.set_ylabel("y [nm]")
        ax.grid(False)

    if mappable is not None:
        cbar = fig.colorbar(mappable, ax=axes[-1], fraction=0.075, pad=0.035)
        cbar.set_label(r"$|j|/j_{\rm avg}$")

    fig.savefig(output, dpi=dpi, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    return output


def plot_ss_snapshot_thermal_scalars(
    power: Mapping[str, np.ndarray],
    snapshots: Mapping[str, np.ndarray] | None,
    dataset: Mapping[str, Any],
    output_path: str | Path,
    *,
    dpi: int = 480,
) -> Path:
    output = _prepare_output(output_path)
    p_ep = _snapshot_array(power, ("P_total_snapshot_W_m3",))
    if p_ep.size == 0:
        raise ValueError("snapshot_power_energy_diagnostics.npz lacks P_total_snapshot_W_m3")

    t_ps = _snapshot_times_ps(power, preferred=("snapshot_t_s",), n=p_ep.shape[0])
    joule = _snapshot_array(power, ("joule_snapshot_W_m3",), shape_like=p_ep)
    p_esc = _snapshot_array(power, ("P_esc_snapshot_W_m3",), shape_like=p_ep)
    p_diff = _snapshot_array(
        power,
        ("P_diff_snapshot_W_m3", "P_diffusion_snapshot_W_m3", "diffusion_power_snapshot_W_m3"),
        shape_like=p_ep,
    )
    te = _snapshot_array(snapshots, ("Te_snapshot_K",))
    tph = _snapshot_array(snapshots, ("Tph_snapshot_K",), shape_like=te)
    u_e = _snapshot_array(power, ("u_e_snapshot_J_m3",), shape_like=p_ep)
    u_ph = _snapshot_array(power, ("u_ph_snapshot_J_m3",), shape_like=p_ep)

    fig, axes = plt.subplots(1, 3, figsize=(15.0, 3.8), constrained_layout=False)
    fig.subplots_adjust(left=0.060, right=0.985, bottom=0.165, top=0.865, wspace=0.28)
    fig.suptitle(f"SS thermal scalar diagnostics: {dataset.get('run_name', '')}", y=0.955)

    ax = axes[0]
    _plot_mean_max_series(ax, t_ps, te, label=r"$T_e$")
    _plot_mean_max_series(ax, t_ps, tph, label=r"$T_{ph}$")
    ax.set_xlabel("t [ps]")
    ax.set_ylabel("temperature [K]")
    ax.set_title("temperatures")
    ax.grid(False)
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1]
    _plot_mean_max_series(ax, t_ps, joule, label=r"$P_J$", use_abs=False)
    _plot_mean_max_series(ax, t_ps, p_ep, label=r"$P_{ep}$", use_abs=False)
    if p_esc.size:
        _plot_mean_max_series(ax, t_ps, p_esc, label=r"$P_{esc}$", use_abs=True)
    if np.any(np.isfinite(p_diff)):
        _plot_mean_max_series(ax, t_ps, p_diff, label=r"$P_{diff}$", use_abs=True)
    ax.set_xlabel("t [ps]")
    ax.set_ylabel(r"power density [W m$^{-3}$]")
    ax.set_title("power densities")
    ax.set_yscale("symlog", linthresh=_linthresh_from_fields(joule, p_ep, p_esc, p_diff, use_abs=True))
    ax.grid(False)
    ax.legend(frameon=False, fontsize=8, ncol=2)

    ax = axes[2]
    _plot_mean_max_series(ax, t_ps, u_e, label=r"$u_e$", use_abs=False)
    _plot_mean_max_series(ax, t_ps, u_ph, label=r"$u_{ph}$", use_abs=False)
    ax.set_xlabel("t [ps]")
    ax.set_ylabel(r"energy density [J m$^{-3}$]")
    ax.set_title("stored energies")
    if _all_positive(u_e, u_ph):
        ax.set_yscale("log")
    ax.grid(False)
    ax.legend(frameon=False, fontsize=8)

    fig.savefig(output, dpi=dpi, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    return output


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
    triangles = np.asarray(dataset.get("triangles", getattr(mesh, "triangles", mesh.elements)), dtype=np.int64)
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
    keys: tuple[str, ...],
    *,
    fallback: np.ndarray | None = None,
    shape_like: np.ndarray | None = None,
) -> np.ndarray:
    if data:
        for key in keys:
            if key in data:
                arr = np.asarray(data[key], dtype=float)
                if arr.ndim == 1:
                    arr = arr[None, :]
                return arr
    if fallback is not None:
        arr = np.asarray(fallback, dtype=float)
        if arr.ndim == 1:
            arr = arr[None, :]
        return arr
    if shape_like is not None and np.asarray(shape_like).size:
        return np.zeros_like(np.asarray(shape_like, dtype=float))
    return np.empty((0, 0), dtype=float)


def _resize_snapshot_field(arr: np.ndarray, shape: tuple[int, int]) -> np.ndarray:
    a = np.asarray(arr, dtype=float)
    if a.size == 0:
        return np.zeros(shape, dtype=float)
    if a.shape == shape:
        return a
    return np.resize(a, shape)


def _snapshot_times_ps(data: Mapping[str, np.ndarray], *, preferred: tuple[str, ...], n: int) -> np.ndarray:
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


def _has_current_family_data(
    dataset: Mapping[str, Any],
    snapshots: Mapping[str, np.ndarray] | None,
    *,
    family: str,
) -> bool:
    mag, jx, jy = _final_current_family_arrays(dataset, snapshots, family=family)
    return bool((mag.size or (jx.size and jy.size)) and jx.size and jy.size)


def _has_thermal_scalar_snapshot_data(
    power: Mapping[str, np.ndarray] | None,
    snapshots: Mapping[str, np.ndarray] | None,
) -> bool:
    return bool(
        _snapshot_array(power, ("P_total_snapshot_W_m3",)).size
        and (
            _snapshot_array(snapshots, ("Te_snapshot_K",)).size
            or _snapshot_array(snapshots, ("Tph_snapshot_K",)).size
            or _snapshot_array(power, ("u_e_snapshot_J_m3",)).size
            or _snapshot_array(power, ("u_ph_snapshot_J_m3",)).size
        )
    )


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


def _plot_mean_max_series(
    ax,
    t_ps: np.ndarray,
    fields: np.ndarray,
    *,
    label: str,
    use_abs: bool = False,
) -> None:
    data = np.asarray(fields, dtype=float)
    if data.size == 0:
        return
    if data.ndim == 1:
        data = data[None, :]
    if data.shape[0] != t_ps.size:
        data = _resize_snapshot_field(data, (t_ps.size, data.shape[-1]))
    if use_abs:
        data = np.abs(data)
    finite = np.isfinite(data)
    if not np.any(finite):
        return
    with np.errstate(invalid="ignore"):
        mean_vals = np.nanmean(np.where(finite, data, np.nan), axis=1)
        max_vals = np.nanmax(np.where(finite, data, np.nan), axis=1)
    line, = ax.plot(t_ps, mean_vals, label=rf"mean {label}")
    ax.plot(t_ps, max_vals, linestyle="--", color=line.get_color(), label=rf"max {label}")


def _linthresh_from_fields(*fields: np.ndarray, use_abs: bool = False) -> float:
    vals: list[np.ndarray] = []
    for field in fields:
        arr = np.asarray(field, dtype=float)
        if arr.size == 0:
            continue
        if use_abs:
            arr = np.abs(arr)
        arr = arr[np.isfinite(arr)]
        arr = np.abs(arr[arr != 0.0])
        if arr.size:
            vals.append(arr)
    if not vals:
        return 1.0
    concat = np.concatenate(vals)
    return max(float(np.nanpercentile(concat, 10.0)), 1.0e-30)


def _all_positive(*fields: np.ndarray) -> bool:
    for field in fields:
        arr = np.asarray(field, dtype=float)
        if arr.size == 0:
            continue
        finite = arr[np.isfinite(arr)]
        if finite.size == 0:
            continue
        if np.any(finite <= 0.0):
            return False
    return True
