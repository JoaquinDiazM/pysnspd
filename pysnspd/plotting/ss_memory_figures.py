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

    if _has_thermal_scalar_data(power, snapshots, dataset):
        saved["snapshot_thermal_scalars"] = plot_ss_snapshot_thermal_scalars(
            power,
            snapshots,
            dataset,
            out / "ss_snapshot_thermal_scalars.png",
            center_width_nm=center_width_nm,
            dpi=dpi,
        )

    if _has_field_map_data(dataset):
        saved["final_center_scalar_maps"] = plot_ss_final_center_scalar_maps(
            mesh,
            dataset,
            out / "ss_final_center_scalar_maps.png",
            center_width_nm=center_width_nm,
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

    xlim, ylim, crop_mask = _center_window(x_nm, y_nm, center_width_nm=center_width_nm)

    families = [
        (js_mag / jscale, js_x / jscale, js_y / jscale),
        (jn_mag / jscale, jn_x / jscale, jn_y / jscale),
        (jtot_mag / jscale, jtot_x / jscale, jtot_y / jscale),
    ]

    finite_for_scale: list[np.ndarray] = []
    for mag, _, _ in families:
        vis = np.asarray(mag[crop_mask], dtype=float)
        vis = vis[np.isfinite(vis)]
        if vis.size:
            finite_for_scale.append(vis)
    if not finite_for_scale:
        raise ValueError("central-strip current map does not contain finite current values")
    vmax = max(1.0, float(np.nanmax(np.concatenate(finite_for_scale))))

    fig = plt.figure(figsize=(10.8, 3.1), constrained_layout=False)
    gs = fig.add_gridspec(1, 4, width_ratios=[1.0, 1.0, 1.0, 0.06], wspace=0.06)
    axes = [fig.add_subplot(gs[0, k]) for k in range(3)]
    cax = fig.add_subplot(gs[0, 3])

    mappable = None
    arrow_length_nm = 0.10 * max(ylim[1] - ylim[0], 1.0)
    for idx, (ax, (mag, jx, jy)) in enumerate(zip(axes, families)):
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
        ax.set_xlabel("x [nm]")
        if idx == 0:
            ax.set_ylabel("y [nm]")
        ax.grid(False)

    if mappable is not None:
        cbar = fig.colorbar(mappable, cax=cax)
        cbar.set_label(r"$|j|/j_{\rm avg}$")

    fig.savefig(output, dpi=dpi, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    return output


def plot_ss_snapshot_thermal_scalars(
    power: Mapping[str, np.ndarray],
    snapshots: Mapping[str, np.ndarray] | None,
    dataset: Mapping[str, Any],
    output_path: str | Path,
    *,
    center_width_nm: float = 100.0,
    dpi: int = 480,
) -> Path:
    output = _prepare_output(output_path)

    x_nm = np.asarray(dataset.get("x_nm", []), dtype=float)
    center_mask = _center_mask_from_x(x_nm, center_width_nm=center_width_nm)
    if not np.any(center_mask):
        center_mask = np.ones_like(x_nm, dtype=bool)

    t_hist = np.asarray(dataset.get("t_ps", []), dtype=float).reshape(-1)
    t_snap = _snapshot_times_ps(power, preferred=("snapshot_t_s",), n=_infer_snapshot_count(power, snapshots))

    # Prefer central snapshot statistics so the curves represent the center strip.
    te_mean, te_max, te_t = _snapshot_series_pair(snapshots, "Te_snapshot_K", t_snap, center_mask)
    tph_mean, tph_max, tph_t = _snapshot_series_pair(snapshots, "Tph_snapshot_K", t_snap, center_mask)
    pJ_mean, pJ_max, pJ_t = _snapshot_series_pair(power, "joule_snapshot_W_m3", t_snap, center_mask)
    pep_mean, pep_max, pep_t = _snapshot_series_pair(power, "P_total_snapshot_W_m3", t_snap, center_mask)
    pesc_mean, pesc_max, pesc_t = _snapshot_series_pair(power, "P_esc_snapshot_W_m3", t_snap, center_mask)
    pdiff_mean, pdiff_max, pdiff_t = _snapshot_series_pair(
        power,
        ("P_diff_snapshot_W_m3", "P_diffusion_snapshot_W_m3", "diffusion_power_snapshot_W_m3"),
        t_snap,
        center_mask,
    )
    ue_mean, ue_max, ue_t = _snapshot_series_pair(power, "u_e_snapshot_J_m3", t_snap, center_mask)
    uph_mean, uph_max, uph_t = _snapshot_series_pair(power, "u_ph_snapshot_J_m3", t_snap, center_mask)

    # Fallback to histories if snapshots are missing.
    if te_mean.size == 0:
        te_t, te_mean, te_max = _history_series_pair(dataset, "thermal_mean_Te_K_history", "thermal_max_Te_K_history", t_hist)
    if tph_mean.size == 0:
        tph_t, tph_mean, tph_max = _history_series_pair(dataset, "thermal_mean_Tph_K_history", "thermal_max_Tph_K_history", t_hist)
    if pJ_mean.size == 0:
        pJ_t, pJ_mean, pJ_max = _history_series_pair(dataset, "thermal_mean_P_J_W_m3_history", "thermal_max_P_J_W_m3_history", t_hist)
    if pep_mean.size == 0:
        pep_t, pep_mean, pep_max = _history_series_pair(dataset, "thermal_mean_P_ep_W_m3_history", "thermal_max_P_ep_W_m3_history", t_hist)
    if pesc_mean.size == 0:
        pesc_t, pesc_mean, pesc_max = _history_series_pair(dataset, "thermal_mean_P_esc_W_m3_history", "thermal_max_P_esc_W_m3_history", t_hist)
    if pdiff_mean.size == 0:
        pdiff_t, pdiff_mean, pdiff_max = _history_series_pair(dataset, "thermal_mean_P_diff_W_m3_history", "thermal_max_P_diff_W_m3_history", t_hist)
    if ue_mean.size == 0:
        ue_t, ue_mean, ue_max = _history_series_pair(dataset, "thermal_mean_u_e_J_m3_history", "thermal_max_u_e_J_m3_history", t_hist)
    if uph_mean.size == 0:
        uph_t, uph_mean, uph_max = _history_series_pair(dataset, "thermal_mean_u_ph_J_m3_history", "thermal_max_u_ph_J_m3_history", t_hist)

    fig, axes = plt.subplots(1, 3, figsize=(13.8, 3.55), constrained_layout=False)
    fig.subplots_adjust(left=0.070, right=0.935, bottom=0.185, top=0.965, wspace=0.36)

    # temperatures
    ax = axes[0]
    _plot_smoothed_series_pair(ax, te_t, te_mean, te_max, label=r"$T_e$", color="tab:blue")
    _plot_smoothed_series_pair(ax, tph_t, tph_mean, tph_max, label=r"$T_{ph}$", color="tab:orange")
    ax.set_xlabel("t [ps]")
    ax.set_ylabel("temperature [K]")
    ax.grid(False)
    ax.legend(frameon=False, fontsize=8)

    # powers
    ax = axes[1]
    _plot_smoothed_series_pair(ax, pJ_t, pJ_mean, pJ_max, label=r"$P_J$", color="tab:blue", take_abs=False)
    _plot_smoothed_series_pair(ax, pep_t, pep_mean, pep_max, label=r"$P_{ep}$", color="tab:orange", take_abs=False)
    _plot_smoothed_series_pair(ax, pesc_t, pesc_mean, pesc_max, label=r"$P_{esc}$", color="tab:green", take_abs=True)
    _plot_smoothed_series_pair(ax, pdiff_t, pdiff_mean, pdiff_max, label=r"$P_{diff}$", color="tab:red", take_abs=True)
    ax.set_xlabel("t [ps]")
    ax.set_ylabel(r"power density [W m$^{-3}$]")
    ax.set_yscale("symlog", linthresh=_linthresh_from_axes_lines(ax))
    ax.grid(False)
    ax.legend(frameon=False, fontsize=8, ncol=2)

    # energies with twin axes and different colors
    ax = axes[2]
    ax_r = ax.twinx()
    _plot_smoothed_series_pair(ax, ue_t, ue_mean, ue_max, label=r"$u_e$", color="tab:blue", take_abs=False)
    _plot_smoothed_series_pair(ax_r, uph_t, uph_mean, uph_max, label=r"$u_{ph}$", color="tab:red", take_abs=False)

    ax.set_xlabel("t [ps]")
    ax.set_ylabel(r"$u_e$ [J m$^{-3}$]", color="tab:blue")
    ax_r.set_ylabel(r"$u_{ph}$ [J m$^{-3}$]", color="tab:red", labelpad=10)

    ax.tick_params(axis="y", colors="tab:blue")
    ax_r.tick_params(axis="y", colors="tab:red")

    if _all_positive_from_axis(ax):
        ax.set_yscale("log")
    elif _has_any_data(ax):
        ax.set_yscale("symlog", linthresh=_linthresh_from_axes_lines(ax))
    if _all_positive_from_axis(ax_r):
        ax_r.set_yscale("log")
    elif _has_any_data(ax_r):
        ax_r.set_yscale("symlog", linthresh=_linthresh_from_axes_lines(ax_r))

    ax.grid(False)
    ax_r.grid(False)

    ax.spines["right"].set_visible(False)
    ax_r.spines["left"].set_visible(False)
    ax_r.spines["bottom"].set_visible(False)
    ax_r.spines["top"].set_visible(False)
    ax_r.patch.set_alpha(0.0)

    handles_l, labels_l = ax.get_legend_handles_labels()
    handles_r, labels_r = ax_r.get_legend_handles_labels()
    ax.legend(handles_l + handles_r, labels_l + labels_r, frameon=False, fontsize=8, ncol=2, loc="best")

    fig.savefig(output, dpi=dpi, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    return output


def plot_ss_final_center_scalar_maps(
    mesh: Any,
    dataset: Mapping[str, Any],
    output_path: str | Path,
    *,
    center_width_nm: float = 100.0,
    dpi: int = 480,
) -> Path:
    output = _prepare_output(output_path)
    tri = _triangulation(mesh, dataset)
    x_nm = np.asarray(tri.x, dtype=float)
    y_nm = np.asarray(tri.y, dtype=float)
    xlim, ylim, _ = _center_window(x_nm, y_nm, center_width_nm=center_width_nm)

    fields = [
        (np.asarray(dataset.get("delta_over_delta0", np.zeros_like(x_nm)), dtype=float), r"$|\Delta|/\Delta_0$"),
        (np.asarray(dataset.get("phi_mV", np.zeros_like(x_nm)), dtype=float), r"$\phi$ [mV]"),
        (np.asarray(dataset.get("q_mag_m_inv", np.zeros_like(x_nm)), dtype=float), r"$|q|$ [m$^{-1}$]"),
        (np.asarray(dataset.get("Te_K", np.zeros_like(x_nm)), dtype=float), r"$T_e$ [K]"),
        (np.asarray(dataset.get("Tph_K", np.zeros_like(x_nm)), dtype=float), r"$T_{ph}$ [K]"),
    ]

    fig = plt.figure(figsize=(15.8, 3.1), constrained_layout=False)
    gs = fig.add_gridspec(1, 10, width_ratios=[1.0, 0.05] * 5, wspace=0.22)
    axes = [fig.add_subplot(gs[0, 2 * k]) for k in range(5)]
    caxes = [fig.add_subplot(gs[0, 2 * k + 1]) for k in range(5)]

    for idx, (ax, cax, (values, cbar_label)) in enumerate(zip(axes, caxes, fields)):
        finite = values[np.isfinite(values)]
        if finite.size == 0:
            values = np.zeros_like(x_nm)
            finite = values
        vmin = float(np.nanmin(finite))
        vmax = float(np.nanmax(finite))
        if np.isclose(vmin, vmax):
            vmax = vmin + 1.0

        mappable = ax.tripcolor(tri, values, shading="gouraud", vmin=vmin, vmax=vmax)
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("x [nm]")
        if idx == 0:
            ax.set_ylabel("y [nm]")
        ax.grid(False)

        cbar = fig.colorbar(mappable, cax=cax)
        cbar.set_label(cbar_label)

    fig.savefig(output, dpi=dpi, bbox_inches="tight", pad_inches=0.05)
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


def _center_window(x_nm: np.ndarray, y_nm: np.ndarray, *, center_width_nm: float) -> tuple[tuple[float, float], tuple[float, float], np.ndarray]:
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


def _snapshot_series_pair(
    data: Mapping[str, np.ndarray] | None,
    keys: str | tuple[str, ...],
    t_ps: np.ndarray,
    center_mask: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    arr = _snapshot_array(data, keys)
    if arr.size == 0:
        return np.array([], dtype=float), np.array([], dtype=float), np.array([], dtype=float)
    if arr.ndim == 1:
        arr = arr[None, :]
    if center_mask.size != arr.shape[1]:
        center_mask = np.resize(center_mask.astype(bool), arr.shape[1])
    if not np.any(center_mask):
        center_mask = np.ones(arr.shape[1], dtype=bool)
    sub = np.asarray(arr[:, center_mask], dtype=float)
    mean_vals = np.nanmean(sub, axis=1)
    max_vals = np.nanmax(sub, axis=1)
    t = np.asarray(t_ps, dtype=float).reshape(-1)
    if t.size != sub.shape[0]:
        t = np.resize(t, sub.shape[0])
    return mean_vals, max_vals, t


def _history_series_pair(
    dataset: Mapping[str, Any],
    mean_key: str,
    max_key: str,
    t_ps: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mean_vals = np.asarray(dataset.get(mean_key, []), dtype=float).reshape(-1)
    max_vals = np.asarray(dataset.get(max_key, []), dtype=float).reshape(-1)
    t = np.asarray(t_ps, dtype=float).reshape(-1)
    n = max(mean_vals.size, max_vals.size)
    if n == 0:
        return np.array([], dtype=float), np.array([], dtype=float), np.array([], dtype=float)
    if t.size != n:
        t = np.resize(t, n)
    if mean_vals.size != n:
        mean_vals = np.resize(mean_vals, n)
    if max_vals.size != n:
        max_vals = np.resize(max_vals, n)
    return t, mean_vals, max_vals


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


def _plot_smoothed_series_pair(
    ax,
    t_ps: np.ndarray,
    mean_vals: np.ndarray,
    max_vals: np.ndarray,
    *,
    label: str,
    color: str,
    take_abs: bool = False,
) -> None:
    t = np.asarray(t_ps, dtype=float).reshape(-1)
    mean_arr = np.asarray(mean_vals, dtype=float).reshape(-1)
    max_arr = np.asarray(max_vals, dtype=float).reshape(-1)
    if mean_arr.size == 0 and max_arr.size == 0:
        return
    n = max(mean_arr.size, max_arr.size)
    if n == 0:
        return
    if t.size != n:
        t = np.resize(t, n)
    if mean_arr.size != n:
        mean_arr = np.resize(mean_arr, n)
    if max_arr.size != n:
        max_arr = np.resize(max_arr, n)
    if take_abs:
        mean_arr = np.abs(mean_arr)
        max_arr = np.abs(max_arr)

    t_smooth_m, y_smooth_m = _smooth_series(t, mean_arr)
    t_smooth_M, y_smooth_M = _smooth_series(t, max_arr)
    if y_smooth_m.size:
        ax.plot(t_smooth_m, y_smooth_m, color=color, label=rf"mean {label}")
    if y_smooth_M.size:
        ax.plot(t_smooth_M, y_smooth_M, linestyle="--", color=color, label=rf"max {label}")


def _smooth_series(t: np.ndarray, y: np.ndarray, *, n_dense: int = 300) -> tuple[np.ndarray, np.ndarray]:
    t = np.asarray(t, dtype=float).reshape(-1)
    y = np.asarray(y, dtype=float).reshape(-1)
    mask = np.isfinite(t) & np.isfinite(y)
    if np.count_nonzero(mask) < 2:
        return t[mask], y[mask]
    t = t[mask]
    y = y[mask]
    order = np.argsort(t)
    t = t[order]
    y = y[order]
    unique_t, unique_idx = np.unique(t, return_index=True)
    t = unique_t
    y = y[unique_idx]
    if t.size < 2:
        return t, y
    t_dense = np.linspace(float(t[0]), float(t[-1]), int(n_dense))
    y_dense = np.interp(t_dense, t, y)
    return t_dense, y_dense


def _infer_snapshot_count(power: Mapping[str, np.ndarray], snapshots: Mapping[str, np.ndarray] | None) -> int:
    for data, key in (
        (power, "P_total_snapshot_W_m3"),
        (snapshots or {}, "Te_snapshot_K"),
        (snapshots or {}, "Tph_snapshot_K"),
        (power, "u_e_snapshot_J_m3"),
        (power, "u_ph_snapshot_J_m3"),
    ):
        if key in data:
            arr = np.asarray(data[key], dtype=float)
            if arr.ndim == 1:
                return 1
            if arr.ndim >= 2:
                return int(arr.shape[0])
    return 0


def _has_nonempty(dataset: Mapping[str, Any], *keys: str) -> bool:
    for key in keys:
        arr = np.asarray(dataset.get(key, []), dtype=float).reshape(-1)
        if arr.size:
            return True
    return False


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
    return max(float(np.nanpercentile(concat, 10.0)), 1.0e-30)


def _all_positive_from_axis(ax) -> bool:
    found = False
    for line in ax.lines:
        y = np.asarray(line.get_ydata(), dtype=float)
        y = y[np.isfinite(y)]
        if y.size:
            found = True
        if y.size and np.any(y <= 0.0):
            return False
    return found


def _has_any_data(ax) -> bool:
    for line in ax.lines:
        y = np.asarray(line.get_ydata(), dtype=float)
        if np.any(np.isfinite(y)):
            return True
    return False
