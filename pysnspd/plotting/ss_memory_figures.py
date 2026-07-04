"""Additional memory-quality stationary SS plots.

This module is plotting-only. It complements the existing
``02_plot_ss_run.py`` pipeline without replacing any current figure.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
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

    fig = plt.figure(figsize=(10.8, 3.0), constrained_layout=False)
    gs = fig.add_gridspec(1, 4, width_ratios=[1.0, 1.0, 1.0, 0.07], wspace=0.08)
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

    fig.savefig(output, dpi=dpi, bbox_inches="tight", pad_inches=0.04)
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

    # Central-strip snapshot statistics
    Te_mean, Te_max, Te_t = _snapshot_series_pair(snapshots, "Te_snapshot_K", t_snap, center_mask)
    Tph_mean, Tph_max, Tph_t = _snapshot_series_pair(snapshots, "Tph_snapshot_K", t_snap, center_mask)
    PJ_mean, PJ_max, PJ_t = _snapshot_series_pair(power, "joule_snapshot_W_m3", t_snap, center_mask)
    Pep_mean, Pep_max, Pep_t = _snapshot_series_pair(power, "P_total_snapshot_W_m3", t_snap, center_mask)
    Pesc_mean, Pesc_max, Pesc_t = _snapshot_series_pair(power, "P_esc_snapshot_W_m3", t_snap, center_mask)
    ue_mean, ue_max, ue_t = _snapshot_series_pair(power, "u_e_snapshot_J_m3", t_snap, center_mask)
    uph_mean, uph_max, uph_t = _snapshot_series_pair(power, "u_ph_snapshot_J_m3", t_snap, center_mask)

    # Diffusion: keep the smoother previous mean/max histories, but filter
    # high-frequency oscillations as requested.
    Pdiff_t, Pdiff_mean, Pdiff_max = _history_series_pair(
        dataset,
        "thermal_mean_P_diff_W_m3_history",
        "thermal_max_P_diff_W_m3_history",
        t_hist,
    )
    Pdiff_t, Pdiff_mean = _filter_history_series(Pdiff_t, Pdiff_mean)
    _, Pdiff_max = _filter_history_series(Pdiff_t, Pdiff_max)

    # Fallbacks if snapshots/histories are missing.
    if Te_mean.size == 0:
        Te_t, Te_mean, Te_max = _history_series_pair(dataset, "thermal_mean_Te_K_history", "thermal_max_Te_K_history", t_hist)
    if Tph_mean.size == 0:
        Tph_t, Tph_mean, Tph_max = _history_series_pair(dataset, "thermal_mean_Tph_K_history", "thermal_max_Tph_K_history", t_hist)
    if PJ_mean.size == 0:
        PJ_t, PJ_mean, PJ_max = _history_series_pair(dataset, "thermal_mean_P_J_W_m3_history", "thermal_max_P_J_W_m3_history", t_hist)
    if Pep_mean.size == 0:
        Pep_t, Pep_mean, Pep_max = _history_series_pair(dataset, "thermal_mean_P_ep_W_m3_history", "thermal_max_P_ep_W_m3_history", t_hist)
    if Pesc_mean.size == 0:
        Pesc_t, Pesc_mean, Pesc_max = _history_series_pair(dataset, "thermal_mean_P_esc_W_m3_history", "thermal_max_P_esc_W_m3_history", t_hist)
    if ue_mean.size == 0:
        ue_t, ue_mean, ue_max = _history_series_pair(dataset, "thermal_mean_u_e_J_m3_history", "thermal_max_u_e_J_m3_history", t_hist)
    if uph_mean.size == 0:
        uph_t, uph_mean, uph_max = _history_series_pair(dataset, "thermal_mean_u_ph_J_m3_history", "thermal_max_u_ph_J_m3_history", t_hist)

    fig, axes = plt.subplots(1, 3, figsize=(15.2, 4.1), constrained_layout=False)
    fig.subplots_adjust(left=0.065, right=0.935, bottom=0.18, top=0.86, wspace=0.55)

    # ------------------------------------------------------------------
    # 1) Temperatures: twin y-axis
    # ------------------------------------------------------------------
    ax1 = axes[0]
    ax1r = ax1.twinx()

    h1 = _plot_series_pair(
        ax1,
        t_ps=Te_t,
        mean_vals=Te_mean,
        max_vals=Te_max,
        color="tab:blue",
        label_mean=r"$T_e$",
        label_max=r"$T_e$",
    )
    h2 = _plot_series_pair(
        ax1r,
        t_ps=Tph_t,
        mean_vals=Tph_mean,
        max_vals=Tph_max,
        color="tab:orange",
        label_mean=r"$T_{ph}$",
        label_max=r"$T_{ph}$",
    )

    ax1.set_xlabel("t [ps]")
    ax1.set_ylabel(r"$T_e$ [K]", color="tab:blue")
    ax1r.set_ylabel(r"$T_{ph}$ [K]", color="tab:orange", labelpad=10)
    ax1.tick_params(axis="y", colors="tab:blue")
    ax1r.tick_params(axis="y", colors="tab:orange")
    ax1.grid(False)
    ax1r.grid(False)
    _clean_twin_axis(ax1, ax1r)
    _apply_temperature_axis_limits(ax1, [Te_mean, Te_max])
    _apply_temperature_axis_limits(ax1r, [Tph_mean, Tph_max])
    _add_dual_legend(
        ax1,
        variable_handles=[h1[0], h2[0]],
        variable_labels=[r"$T_e$", r"$T_{ph}$"],
        style_loc="upper right",
        variable_loc="upper center",
        variable_ncol=2,
    )

    # ------------------------------------------------------------------
    # 2) Powers: twin y-axis
    # left  -> P_J and P_diff
    # right -> P_ep and P_esc
    # ------------------------------------------------------------------
    ax2 = axes[1]
    ax2r = ax2.twinx()

    hPJ = _plot_series_pair(
        ax2,
        t_ps=PJ_t,
        mean_vals=PJ_mean,
        max_vals=PJ_max,
        color="tab:blue",
        label_mean=r"$P_J$",
        label_max=r"$P_J$",
        take_abs=False,
    )
    hPd = _plot_series_pair(
        ax2,
        t_ps=Pdiff_t,
        mean_vals=Pdiff_mean,
        max_vals=Pdiff_max,
        color="tab:red",
        label_mean=r"$P_{diff}$",
        label_max=r"$P_{diff}$",
        take_abs=True,
        force_history_smoothing=True,
    )
    hPep = _plot_series_pair(
        ax2r,
        t_ps=Pep_t,
        mean_vals=Pep_mean,
        max_vals=Pep_max,
        color="tab:orange",
        label_mean=r"$P_{ep}$",
        label_max=r"$P_{ep}$",
        take_abs=True,
    )
    hPesc = _plot_series_pair(
        ax2r,
        t_ps=Pesc_t,
        mean_vals=Pesc_mean,
        max_vals=Pesc_max,
        color="tab:green",
        label_mean=r"$P_{esc}$",
        label_max=r"$P_{esc}$",
        take_abs=True,
    )

    ax2.set_xlabel("t [ps]")
    ax2.set_ylabel(r"$P_J,\ P_{diff}$ [W m$^{-3}$]", color="black")
    ax2r.set_ylabel(r"$P_{ep},\ P_{esc}$ [W m$^{-3}$]", color="black", labelpad=10)
    ax2.set_yscale("symlog", linthresh=_linthresh_from_axes_lines(ax2))
    ax2r.set_yscale("symlog", linthresh=_linthresh_from_axes_lines(ax2r))
    ax2.grid(False)
    ax2r.grid(False)
    _clean_twin_axis(ax2, ax2r)
    _add_dual_legend(
        ax2,
        variable_handles=[hPJ[0], hPd[0], hPep[0], hPesc[0]],
        variable_labels=[r"$P_J$", r"$P_{diff}$", r"$P_{ep}$", r"$P_{esc}$"],
        style_loc="upper right",
        variable_loc="upper center",
        variable_ncol=4,
    )

    # ------------------------------------------------------------------
    # 3) Energies: twin y-axis
    # ------------------------------------------------------------------
    ax3 = axes[2]
    ax3r = ax3.twinx()

    hue = _plot_series_pair(
        ax3,
        t_ps=ue_t,
        mean_vals=ue_mean,
        max_vals=ue_max,
        color="tab:blue",
        label_mean=r"$u_e$",
        label_max=r"$u_e$",
    )
    huph = _plot_series_pair(
        ax3r,
        t_ps=uph_t,
        mean_vals=uph_mean,
        max_vals=uph_max,
        color="tab:red",
        label_mean=r"$u_{ph}$",
        label_max=r"$u_{ph}$",
    )

    ax3.set_xlabel("t [ps]")
    ax3.set_ylabel(r"$u_e$ [J m$^{-3}$]", color="tab:blue")
    ax3r.set_ylabel(r"$u_{ph}$ [J m$^{-3}$]", color="tab:red", labelpad=10)
    ax3.tick_params(axis="y", colors="tab:blue")
    ax3r.tick_params(axis="y", colors="tab:red")
    if _all_positive_arrays([ue_mean, ue_max]):
        ax3.set_yscale("log")
    elif _has_nonzero_arrays([ue_mean, ue_max]):
        ax3.set_yscale("symlog", linthresh=_linthresh_from_arrays([ue_mean, ue_max]))
    if _all_positive_arrays([uph_mean, uph_max]):
        ax3r.set_yscale("log")
    elif _has_nonzero_arrays([uph_mean, uph_max]):
        ax3r.set_yscale("symlog", linthresh=_linthresh_from_arrays([uph_mean, uph_max]))
    ax3.grid(False)
    ax3r.grid(False)
    _clean_twin_axis(ax3, ax3r)
    _add_dual_legend(
        ax3,
        variable_handles=[hue[0], huph[0]],
        variable_labels=[r"$u_e$", r"$u_{ph}$"],
        style_loc="upper right",
        variable_loc="upper center",
        variable_ncol=2,
    )

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
    xlim, ylim, crop_mask = _center_window(x_nm, y_nm, center_width_nm=center_width_nm)

    fields = [
        (np.asarray(dataset.get("delta_over_delta0", np.zeros_like(x_nm)), dtype=float), r"$|\Delta|/\Delta_0$"),
        (np.asarray(dataset.get("phi_mV", np.zeros_like(x_nm)), dtype=float), r"$\phi$ [mV]"),
        (np.asarray(dataset.get("q_mag_m_inv", np.zeros_like(x_nm)), dtype=float), r"$|q|$ [m$^{-1}$]"),
        (np.asarray(dataset.get("Te_K", np.zeros_like(x_nm)), dtype=float), r"$T_e$ [K]"),
        (np.asarray(dataset.get("Tph_K", np.zeros_like(x_nm)), dtype=float), r"$T_{ph}$ [K]"),
    ]

    # Horizontal colorbars above each subplot avoid the previous side-overlap.
    fig = plt.figure(figsize=(17.6, 3.65), constrained_layout=False)
    gs = fig.add_gridspec(
        2,
        5,
        height_ratios=[0.10, 1.00],
        wspace=0.18,
        hspace=0.10,
    )
    caxes = [fig.add_subplot(gs[0, k]) for k in range(5)]
    axes = [fig.add_subplot(gs[1, k]) for k in range(5)]

    for idx, (ax, cax, (values, cbar_label)) in enumerate(zip(axes, caxes, fields)):
        values = np.asarray(values, dtype=float)
        vis = values[crop_mask]
        vis = vis[np.isfinite(vis)]
        if vis.size == 0:
            vis = np.array([0.0, 1.0], dtype=float)
        vmin = float(np.nanmin(vis))
        vmax = float(np.nanmax(vis))
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

        cbar = fig.colorbar(mappable, cax=cax, orientation="horizontal")
        cbar.set_label(cbar_label, labelpad=2.0)
        cbar.ax.xaxis.set_ticks_position("top")
        cbar.ax.xaxis.set_label_position("top")

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


def _plot_series_pair(
    ax,
    *,
    t_ps: np.ndarray,
    mean_vals: np.ndarray,
    max_vals: np.ndarray,
    color: str,
    label_mean: str,
    label_max: str,
    take_abs: bool = False,
    force_history_smoothing: bool = False,
):
    t = np.asarray(t_ps, dtype=float).reshape(-1)
    mean_arr = np.asarray(mean_vals, dtype=float).reshape(-1)
    max_arr = np.asarray(max_vals, dtype=float).reshape(-1)
    if mean_arr.size == 0 and max_arr.size == 0:
        empty1, = ax.plot([], [], color=color, linestyle="-", label=label_mean)
        empty2, = ax.plot([], [], color=color, linestyle="--", label=label_max)
        return empty1, empty2

    n = max(mean_arr.size, max_arr.size)
    if t.size != n:
        t = np.resize(t, n)
    if mean_arr.size != n:
        mean_arr = np.resize(mean_arr, n)
    if max_arr.size != n:
        max_arr = np.resize(max_arr, n)
    if take_abs:
        mean_arr = np.abs(mean_arr)
        max_arr = np.abs(max_arr)

    if force_history_smoothing:
        t_mean, y_mean = _filter_history_series(t, mean_arr)
        t_max, y_max = _filter_history_series(t, max_arr)
    else:
        t_mean, y_mean = _smart_smooth_series(t, mean_arr)
        t_max, y_max = _smart_smooth_series(t, max_arr)

    h_mean, = ax.plot(t_mean, y_mean, color=color, linestyle="-", label=label_mean)
    h_max, = ax.plot(t_max, y_max, color=color, linestyle="--", label=label_max)
    return h_mean, h_max


def _smart_smooth_series(t: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
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
    t, unique_idx = np.unique(t, return_index=True)
    y = y[unique_idx]
    if t.size < 3:
        return t, y
    if t.size < 7:
        return _monotone_cubic_smooth(t, y, n_dense=320)
    return _filter_history_series(t, y, prefer_dense=False)


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


def _filter_history_series(
    t: np.ndarray,
    y: np.ndarray,
    *,
    prefer_dense: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
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

    n = t.size
    if n < 5:
        return _monotone_cubic_smooth(t, y, n_dense=280)

    # moving-average filter with scale-aware window
    win = max(5, int(n // 180))
    if win % 2 == 0:
        win += 1
    kernel = np.ones(win, dtype=float) / float(win)
    y_pad = np.pad(y, (win // 2, win // 2), mode="edge")
    y_smooth = np.convolve(y_pad, kernel, mode="valid")

    if prefer_dense and n < 1200:
        return _monotone_cubic_smooth(t, y_smooth, n_dense=min(700, max(300, 3 * n)))

    # decimate very long histories to keep PNGs light and readable
    target = 900
    if n > target:
        idx = np.linspace(0, n - 1, target).astype(int)
        t = t[idx]
        y_smooth = y_smooth[idx]
    return t, y_smooth


def _clean_twin_axis(ax, ax_r) -> None:
    ax.spines["right"].set_visible(False)
    ax_r.spines["left"].set_visible(False)
    ax_r.patch.set_alpha(0.0)


def _apply_temperature_axis_limits(ax, arrays: Sequence[np.ndarray]) -> None:
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
        pad = 0.06 * (vmax - vmin)
    ax.set_ylim(vmin - pad, vmax + pad)


def _add_dual_legend(
    ax,
    *,
    variable_handles,
    variable_labels,
    style_loc: str,
    variable_loc: str,
    variable_ncol: int,
) -> None:
    style_handles = [
        Line2D([0], [0], color="black", linestyle="-", label="mean"),
        Line2D([0], [0], color="black", linestyle="--", label="max"),
    ]
    leg_vars = ax.legend(
        variable_handles,
        variable_labels,
        loc=variable_loc,
        bbox_to_anchor=(0.50, 1.14),
        ncol=variable_ncol,
        frameon=False,
        columnspacing=1.2,
        handlelength=2.6,
        borderaxespad=0.0,
    )
    ax.add_artist(leg_vars)
    ax.legend(
        style_handles,
        ["mean", "max"],
        loc=style_loc,
        frameon=False,
        ncol=2,
        columnspacing=1.0,
        handlelength=2.4,
        borderaxespad=0.3,
    )


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


def _linthresh_from_arrays(arrays: Sequence[np.ndarray]) -> float:
    vals = []
    for arr in arrays:
        y = np.asarray(arr, dtype=float).reshape(-1)
        y = y[np.isfinite(y)]
        y = np.abs(y[y != 0.0])
        if y.size:
            vals.append(y)
    if not vals:
        return 1.0
    concat = np.concatenate(vals)
    return max(float(np.nanpercentile(concat, 12.0)), 1.0e-30)


def _all_positive_arrays(arrays: Sequence[np.ndarray]) -> bool:
    found = False
    for arr in arrays:
        y = np.asarray(arr, dtype=float).reshape(-1)
        y = y[np.isfinite(y)]
        if y.size:
            found = True
        if y.size and np.any(y <= 0.0):
            return False
    return found


def _has_nonzero_arrays(arrays: Sequence[np.ndarray]) -> bool:
    for arr in arrays:
        y = np.asarray(arr, dtype=float).reshape(-1)
        y = y[np.isfinite(y)]
        if y.size and np.any(np.abs(y) > 0.0):
            return True
    return False


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
