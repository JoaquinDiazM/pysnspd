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

from pysnspd.plotting.ss_power_figures import _snapshot_diffusion_power_density


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
        mappable = ax.tripcolor(tri, np.asarray(mag, dtype=float), shading="gouraud", vmin=0.0, vmax=vmax)
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
    """Plot compact memory-ready thermal scalar diagnostics.

    The power panel uses the same runtime-history power envelopes used by
    ``plot_ss_snapshot_runtime_metrics`` whenever those histories are present.
    Snapshot maps are kept only as fallbacks.
    """

    output = _prepare_output(output_path)

    p_ep = _snapshot_array(power, ("P_total_snapshot_W_m3",))
    if p_ep.size == 0:
        raise ValueError("snapshot_power_energy_diagnostics.npz lacks P_total_snapshot_W_m3")

    joule = _snapshot_array(power, ("joule_snapshot_W_m3",), shape_like=p_ep)
    p_esc = _snapshot_array(power, ("P_esc_snapshot_W_m3",), shape_like=p_ep)
    u_e = _snapshot_array(power, ("u_e_snapshot_J_m3",), shape_like=p_ep)
    u_ph = _snapshot_array(power, ("u_ph_snapshot_J_m3",), shape_like=p_ep)

    p_diff = _snapshot_diffusion_power_density(
        None,
        snapshots=snapshots,
        power=power,
        dataset=dataset,
        shape_like=p_ep,
    )

    t_ps_snap = _snapshot_times_ps(power, preferred=("snapshot_t_s",), n=p_ep.shape[0])
    t_ps_hist = np.asarray(dataset.get("t_ps", []), dtype=float).reshape(-1)

    x_nm = np.asarray(dataset.get("x_nm", []), dtype=float)
    center_mask = _center_mask_from_x(x_nm, center_width_nm=center_width_nm)
    if center_mask.size != p_ep.shape[1]:
        center_mask = np.resize(center_mask.astype(bool), p_ep.shape[1])
    if not np.any(center_mask):
        center_mask = np.ones(p_ep.shape[1], dtype=bool)

    te_snap = np.empty((0, 0), dtype=float)
    tph_snap = np.empty((0, 0), dtype=float)
    if snapshots:
        te_snap = _snapshot_array(snapshots, ("Te_snapshot_K",), shape_like=p_ep)
        tph_snap = _snapshot_array(snapshots, ("Tph_snapshot_K",), shape_like=p_ep)

    # Temperatures: history if present, snapshot fallback.
    te_t, te_mean, te_max = _history_or_snapshot_pair(
        dataset,
        t_ps_hist,
        mean_key="thermal_mean_Te_K_history",
        max_key="thermal_max_Te_K_history",
        snapshot_t=t_ps_snap,
        snapshot_values=te_snap,
        center_mask=center_mask,
    )
    tph_t, tph_mean, tph_max = _history_or_snapshot_pair(
        dataset,
        t_ps_hist,
        mean_key="thermal_mean_Tph_K_history",
        max_key="thermal_max_Tph_K_history",
        snapshot_t=t_ps_snap,
        snapshot_values=tph_snap,
        center_mask=center_mask,
    )

    # Powers: runtime histories first, exactly to match the runtime-metrics panel.
    pj_t, pj_mean, pj_max = _history_or_snapshot_pair(
        dataset,
        t_ps_hist,
        mean_key="thermal_mean_P_J_W_m3_history",
        max_key="thermal_max_P_J_W_m3_history",
        snapshot_t=t_ps_snap,
        snapshot_values=joule,
        center_mask=center_mask,
        max_mode="max",
        abs_snapshot=False,
    )
    pep_t, pep_mean, pep_max = _history_or_snapshot_pair(
        dataset,
        t_ps_hist,
        mean_key="thermal_mean_P_ep_W_m3_history",
        max_key="thermal_max_P_ep_W_m3_history",
        snapshot_t=t_ps_snap,
        snapshot_values=p_ep,
        center_mask=center_mask,
        max_mode="max_abs",
        abs_snapshot=True,
    )
    pesc_t, pesc_mean, pesc_max = _history_or_snapshot_pair(
        dataset,
        t_ps_hist,
        mean_key="thermal_mean_P_esc_W_m3_history",
        max_key="thermal_max_P_esc_W_m3_history",
        snapshot_t=t_ps_snap,
        snapshot_values=p_esc,
        center_mask=center_mask,
        max_mode="max_abs",
        abs_snapshot=True,
    )
    pdiff_t, pdiff_mean, pdiff_max = _history_or_snapshot_pair(
        dataset,
        t_ps_hist,
        mean_key="thermal_mean_P_diff_W_m3_history",
        max_key="thermal_max_P_diff_W_m3_history",
        snapshot_t=t_ps_snap,
        snapshot_values=p_diff,
        center_mask=center_mask,
        max_mode="max_abs",
        abs_snapshot=True,
        smooth_history=True,
    )

    # Energies: copy the data extraction used by
    # plot_ss_snapshot_runtime_metrics / stored energy diagnostics.
    # That panel plots spatial means over the full snapshot map, not center-strip
    # maxima.  We keep the memory-plot visual design, but separate u_e and u_ph
    # into two y axes because their operating ranges are very different.
    ue_t, ue_mean = _snapshot_mean_series(t_ps_snap, u_e)
    uph_t, uph_mean = _snapshot_mean_series(t_ps_snap, u_ph)

    # Colors
    te_color = "tab:blue"
    tph_color = "tab:orange"
    pj_color = "#d95f02"      # warm orange
    pdiff_color = "#b2182b"   # warm red
    pep_color = "#1f78b4"     # cool blue
    pesc_color = "#1b9e77"    # cool teal
    ue_color = "tab:blue"
    uph_color = "tab:red"

    prethermal_t_ps = 2.0

    fig, axes = plt.subplots(1, 3, figsize=(15.6, 4.10), constrained_layout=False)
    fig.subplots_adjust(left=0.070, right=0.955, bottom=0.18, top=0.86, wspace=0.56)

    # ------------------------------------------------------------------
    # 1) Temperatures
    # ------------------------------------------------------------------
    ax = axes[0]
    ax_r = ax.twinx()
    _shade_prethermal(ax, prethermal_t_ps)

    h_te = _plot_series_pair(ax, te_t, te_mean, te_max, color=te_color, take_abs=False)
    h_tph = _plot_series_pair(ax_r, tph_t, tph_mean, tph_max, color=tph_color, take_abs=False)

    ax.set_xlabel("t [ps]")
    ax.set_ylabel(r"$T_e$ [K]", color=te_color)
    ax_r.set_ylabel(r"$T_{ph}$ [K]", color=tph_color, labelpad=10)
    ax.tick_params(axis="y", colors=te_color)
    ax_r.tick_params(axis="y", colors=tph_color)
    ax.grid(False)
    ax_r.grid(False)
    _clean_twin_axis(ax, ax_r)
    _apply_axis_limits(ax, [te_mean, te_max], frac=0.07)
    _apply_axis_limits(ax_r, [tph_mean, tph_max], frac=0.07)
    _add_legends(
        ax,
        variable_handles=[h_te[0], h_tph[0]],
        variable_labels=[r"$T_e$", r"$T_{ph}$"],
        variable_ncol=2,
    )

    # ------------------------------------------------------------------
    # 2) Powers
    # ------------------------------------------------------------------
    ax = axes[1]
    ax_r = ax.twinx()
    _shade_prethermal(ax, prethermal_t_ps)

    h_pj = _plot_series_pair(ax, pj_t, pj_mean, pj_max, color=pj_color, take_abs=True)
    h_pd = _plot_series_pair(ax, pdiff_t, pdiff_mean, pdiff_max, color=pdiff_color, take_abs=True)
    h_pep = _plot_series_pair(ax_r, pep_t, pep_mean, pep_max, color=pep_color, take_abs=True)
    h_pesc = _plot_series_pair(ax_r, pesc_t, pesc_mean, pesc_max, color=pesc_color, take_abs=True)

    ax.set_xlabel("t [ps]")
    ax.set_ylabel(r"$P_J,\ P_{diff}$ [W m$^{-3}$]")
    ax_r.set_ylabel(r"$P_{ep},\ P_{esc}$ [W m$^{-3}$]", labelpad=10)
    ax.tick_params(axis="y", colors="black")
    ax_r.tick_params(axis="y", colors="black")

    ax.set_yscale("symlog", linthresh=_linthresh_from_axes_lines(ax))
    ax_r.set_yscale("symlog", linthresh=_linthresh_from_axes_lines(ax_r))
    _autoscale_symlog_axis(ax)
    _autoscale_symlog_axis(ax_r)

    ax.grid(False)
    ax_r.grid(False)
    _clean_twin_axis(ax, ax_r)
    _add_legends(
        ax,
        variable_handles=[h_pj[0], h_pd[0], h_pep[0], h_pesc[0]],
        variable_labels=[r"$P_J$", r"$P_{diff}$", r"$P_{ep}$", r"$P_{esc}$"],
        variable_ncol=4,
    )

    # ------------------------------------------------------------------
    # 3) Energies
    # ------------------------------------------------------------------
    ax = axes[2]
    ax_r = ax.twinx()
    _shade_prethermal(ax, prethermal_t_ps)

    h_ue = _plot_single_series(ax, ue_t, ue_mean, color=ue_color)
    h_uph = _plot_single_series(ax_r, uph_t, uph_mean, color=uph_color)

    ax.set_xlabel("t [ps]")
    ax.set_ylabel(r"$u_e$ [J m$^{-3}$]", color=ue_color)
    ax_r.set_ylabel(r"$u_{ph}$ [J m$^{-3}$]", color=uph_color, labelpad=10)
    ax.tick_params(axis="y", colors=ue_color)
    ax_r.tick_params(axis="y", colors=uph_color)

    # u_e can be negative in the current normalization; center the axis around
    # the negative operating value rather than forcing log/symlog.  u_ph is
    # shown independently on the right axis because it lives close to zero and
    # varies on a much smaller absolute scale.
    _center_linear_axis(ax, [ue_mean], frac=0.10)
    _center_linear_axis(ax_r, [uph_mean], frac=0.25)

    fmt_right = mticker.ScalarFormatter(useMathText=True)
    fmt_right.set_scientific(True)
    fmt_right.set_powerlimits((-3, 3))
    fmt_right.set_useOffset(True)
    ax_r.yaxis.set_major_formatter(fmt_right)
    ax_r.yaxis.set_major_locator(mticker.MaxNLocator(nbins=4))
    ax_r.yaxis.get_offset_text().set_color(uph_color)

    ax.grid(False)
    ax_r.grid(False)
    _clean_twin_axis(ax, ax_r)
    _add_variable_legend(
        ax,
        variable_handles=[h_ue, h_uph],
        variable_labels=[r"$u_e$", r"$u_{ph}$"],
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

    def _node_array_from_dataset(keys: tuple[str, ...], *, default: float = 0.0) -> np.ndarray:
        for key in keys:
            if key in dataset:
                arr = np.asarray(dataset[key], dtype=float).reshape(-1)
                if arr.size:
                    if arr.size != x_nm.size:
                        arr = np.resize(arr, x_nm.size)
                    return arr
        return np.full(x_nm.size, float(default), dtype=float)

    def _center_limits(
        values: np.ndarray,
        *,
        force_vmin: float | None = None,
        force_vmax: float | None = None,
        positive_floor: bool = False,
    ) -> tuple[float, float]:
        vals = np.asarray(values, dtype=float)
        if vals.size != x_nm.size:
            vals = np.resize(vals, x_nm.size)

        vis = vals[crop_mask]
        vis = vis[np.isfinite(vis)]

        if vis.size == 0:
            vmin = 0.0
            vmax = 1.0
        else:
            vmin = float(np.nanmin(vis))
            vmax = float(np.nanmax(vis))

        if positive_floor:
            vmin = 0.0
        if force_vmin is not None:
            vmin = float(force_vmin)
        if force_vmax is not None:
            vmax = float(force_vmax)

        if not np.isfinite(vmin) or not np.isfinite(vmax):
            vmin, vmax = 0.0, 1.0

        if np.isclose(vmin, vmax):
            pad = max(1.0e-12, 0.02 * max(abs(vmin), 1.0))
            if force_vmin is None:
                vmin -= pad
            if force_vmax is None:
                vmax += pad

        return vmin, vmax

    def _phase_gradient_q_abs_m_inv() -> np.ndarray:
        psi = np.asarray(dataset.get("psi_J", []), dtype=np.complex128).reshape(-1)
        if psi.size != x_nm.size:
            return np.zeros(x_nm.size, dtype=float)

        theta = np.angle(psi)
        x_m = x_nm * 1.0e-9
        y_m = y_nm * 1.0e-9

        triangles = np.asarray(tri.triangles, dtype=np.int64)
        if triangles.size == 0:
            return np.zeros(x_nm.size, dtype=float)

        edges = np.vstack((triangles[:, [0, 1]], triangles[:, [1, 2]], triangles[:, [2, 0]]))
        edges = np.sort(edges, axis=1)
        edges = np.unique(edges, axis=0)

        Axx = np.zeros(x_nm.size, dtype=float)
        Axy = np.zeros(x_nm.size, dtype=float)
        Ayy = np.zeros(x_nm.size, dtype=float)
        bx = np.zeros(x_nm.size, dtype=float)
        by = np.zeros(x_nm.size, dtype=float)

        for i, j in edges:
            dx = float(x_m[j] - x_m[i])
            dy = float(y_m[j] - y_m[i])
            if not np.isfinite(dx) or not np.isfinite(dy):
                continue

            dtheta = float(np.angle(np.exp(1j * (theta[j] - theta[i]))))

            Axx[i] += dx * dx
            Axy[i] += dx * dy
            Ayy[i] += dy * dy
            bx[i] += dtheta * dx
            by[i] += dtheta * dy

            Axx[j] += dx * dx
            Axy[j] += dx * dy
            Ayy[j] += dy * dy
            bx[j] += dtheta * dx
            by[j] += dtheta * dy

        det = Axx * Ayy - Axy * Axy
        good = np.isfinite(det) & (np.abs(det) > 1.0e-300)

        qx = np.zeros(x_nm.size, dtype=float)
        qy = np.zeros(x_nm.size, dtype=float)
        qx[good] = (Ayy[good] * bx[good] - Axy[good] * by[good]) / det[good]
        qy[good] = (-Axy[good] * bx[good] + Axx[good] * by[good]) / det[good]

        q_abs = np.sqrt(qx * qx + qy * qy)
        q_abs[~np.isfinite(q_abs)] = 0.0
        return q_abs

    def _q_abs_m_inv() -> np.ndarray:
        q_saved = _node_array_from_dataset(
            (
                "q_mag_m_inv",
                "q_abs_m_inv",
                "node_q_mag_m_inv",
                "node_q_abs_m_inv",
                "node_superfluid_momentum_mag_m_inv",
            )
        )
        q_center = q_saved[crop_mask]
        if np.any(np.isfinite(q_center)) and float(np.nanmax(np.abs(q_center))) > 0.0:
            return q_saved
        return _phase_gradient_q_abs_m_inv()

    delta_norm = _node_array_from_dataset(("delta_over_delta0",))

    phi_uV = 1.0e6 * _node_array_from_dataset(("phi_V",))
    if not np.any(np.isfinite(phi_uV)) or np.nanmax(np.abs(phi_uV)) == 0.0:
        phi_uV = 1.0e3 * _node_array_from_dataset(("phi_mV",))

    q_abs = _q_abs_m_inv()
    Te_K = _node_array_from_dataset(("Te_K",))
    Tph_K = _node_array_from_dataset(("Tph_K",))

    fields = [
        {
            "values": delta_norm,
            "label": r"$|\Delta|/\Delta_0$",
            "vmin": 0.0,
            "vmax": 1.0,
            "positive_floor": True,
        },
        {
            "values": phi_uV,
            "label": r"$\phi$ [$\mu$V]",
            "vmin": None,
            "vmax": None,
            "positive_floor": False,
        },
        {
            "values": q_abs,
            "label": r"$|q|$ [m$^{-1}$]",
            "vmin": 0.0,
            "vmax": None,
            "positive_floor": True,
        },
        {
            "values": Te_K,
            "label": r"$T_e$ [K]",
            "vmin": None,
            "vmax": None,
            "positive_floor": False,
        },
        {
            "values": Tph_K,
            "label": r"$T_{ph}$ [K]",
            "vmin": None,
            "vmax": None,
            "positive_floor": False,
        },
    ]

    fig = plt.figure(figsize=(18.4, 3.75), constrained_layout=False)
    gs = fig.add_gridspec(2, 5, height_ratios=[0.10, 1.00], wspace=0.24, hspace=0.12)
    caxes = [fig.add_subplot(gs[0, k]) for k in range(5)]
    axes = [fig.add_subplot(gs[1, k]) for k in range(5)]

    for idx, (ax, cax, spec) in enumerate(zip(axes, caxes, fields)):
        values = np.asarray(spec["values"], dtype=float).reshape(-1)
        if values.size != x_nm.size:
            values = np.resize(values, x_nm.size)

        vmin, vmax = _center_limits(
            values,
            force_vmin=spec["vmin"],
            force_vmax=spec["vmax"],
            positive_floor=bool(spec["positive_floor"]),
        )

        mappable = ax.tripcolor(tri, values, shading="gouraud", vmin=vmin, vmax=vmax)
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("x [nm]")
        if idx == 0:
            ax.set_ylabel("y [nm]")
        ax.grid(False)

        cbar = fig.colorbar(mappable, cax=cax, orientation="horizontal")
        cbar.set_label(str(spec["label"]), labelpad=2.0)
        cbar.ax.xaxis.set_ticks_position("top")
        cbar.ax.xaxis.set_label_position("top")
        cbar.ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=3))

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


def _add_legends(
    ax,
    *,
    variable_handles,
    variable_labels,
    variable_ncol: int,
) -> None:
    style_handles = [
        Line2D([0], [0], color="black", linestyle="-", label="mean"),
        Line2D([0], [0], color="black", linestyle="--", label="max"),
    ]
    leg_vars = ax.legend(
        variable_handles,
        variable_labels,
        loc="upper center",
        bbox_to_anchor=(0.50, 1.14),
        ncol=variable_ncol,
        frameon=False,
        columnspacing=1.1,
        handlelength=2.5,
        borderaxespad=0.0,
    )
    ax.add_artist(leg_vars)
    ax.legend(
        style_handles,
        ["mean", "max"],
        loc="upper right",
        frameon=False,
        ncol=2,
        columnspacing=1.0,
        handlelength=2.4,
        borderaxespad=0.3,
    )


def _add_variable_legend(
    ax,
    *,
    variable_handles,
    variable_labels,
    variable_ncol: int,
) -> None:
    ax.legend(
        variable_handles,
        variable_labels,
        loc="upper center",
        bbox_to_anchor=(0.50, 1.14),
        ncol=variable_ncol,
        frameon=False,
        columnspacing=1.2,
        handlelength=2.6,
        borderaxespad=0.0,
    )


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


def _linthresh_from_arrays(arrays: list[np.ndarray]) -> float:
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


def _all_positive_arrays(arrays: list[np.ndarray]) -> bool:
    found = False
    for arr in arrays:
        y = np.asarray(arr, dtype=float).reshape(-1)
        y = y[np.isfinite(y)]
        if y.size:
            found = True
        if y.size and np.any(y <= 0.0):
            return False
    return found


def _has_nonzero_arrays(arrays: list[np.ndarray]) -> bool:
    for arr in arrays:
        y = np.asarray(arr, dtype=float).reshape(-1)
        y = y[np.isfinite(y)]
        if y.size and np.any(np.abs(y) > 0.0):
            return True
    return False


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
