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

    panel_labels = [
        r"$j_s^{\mathrm{Us}}$",
        r"$j_n$",
        r"$j_{\mathrm{tot}}$",
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
        mappable = ax.tripcolor(
            tri,
            np.asarray(mag, dtype=float),
            shading="gouraud",
            vmin=0.0,
            vmax=vmax,
        )
        _overlay_current_arrows(
            ax,
            x_nm,
            y_nm,
            jx,
            jy,
            crop_mask,
            arrow_length_nm=arrow_length_nm,
        )

        ax.text(
            0.03,
            0.97,
            panel_labels[idx],
            transform=ax.transAxes,
            ha="left",
            va="top",
            color="white",
            fontsize=14,
            fontweight="bold",
            bbox=dict(
                boxstyle="round,pad=0.22",
                facecolor="#ff0000",
                edgecolor="#ff0000",
                linewidth=0.0,
            ),
            zorder=10,
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

    - Temperaturas: mean + max.
    - Potencias: solo máximos, con línea sólida.
    - Energías: medias espaciales tipo ss_snapshot_runtime_metrics.
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

    # ------------------------------------------------------------------
    # Temperaturas: history si existe, fallback snapshot.
    # ------------------------------------------------------------------
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

    # ------------------------------------------------------------------
    # Potencias: history primero, fallback snapshot.
    # Queremos conservar SOLO los máximos.
    # ------------------------------------------------------------------
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

    # ------------------------------------------------------------------
    # Energías: igual que stored energy diagnostics de
    # ss_snapshot_runtime_metrics: media espacial del mapa completo.
    # ------------------------------------------------------------------
    ue_t, ue_mean = _snapshot_mean_series(t_ps_snap, u_e)
    uph_t, uph_mean = _snapshot_mean_series(t_ps_snap, u_ph)

    # Colors
    te_color = "tab:blue"
    tph_color = "tab:orange"
    pj_color = "#ffa200"
    pdiff_color = "#ff0000"
    hot_color = "#ff5100"
    pep_color = "#7700ff"
    pesc_color = "#00c3ff"
    cold_color = "#0011ff"
    ue_color = "tab:blue"
    uph_color = "tab:red"

    prethermal_t_ps = 2.0

    fig, axes = plt.subplots(1, 3, figsize=(15.6, 4.10), constrained_layout=False)
    fig.subplots_adjust(left=0.070, right=0.955, bottom=0.18, top=0.86, wspace=0.66)

    # ------------------------------------------------------------------
    # 1) Temperatures
    # ------------------------------------------------------------------
    ax = axes[0]
    ax_r = ax.twinx()
    _shade_prethermal(ax, prethermal_t_ps)

    h_te = _plot_series_pair(ax, te_t, te_mean, te_max, color=te_color, take_abs=False)
    h_tph = _plot_series_pair(ax_r, tph_t, tph_mean, tph_max, color=tph_color, take_abs=False)

    ax.set_xlabel("t [ps]", fontsize=14)
    ax.set_ylabel(r"$T_e$ [K]", color=te_color, fontsize=14)
    ax_r.set_ylabel(r"$T_{ph}$ [K]", color=tph_color, fontsize=14)
    ax.tick_params(axis="x", labelsize=14)
    ax.tick_params(axis="y", colors=te_color, labelsize=14)
    ax_r.tick_params(axis="y", colors=tph_color, labelsize=14)
    ax.grid(False)
    ax_r.grid(False)
    _clean_twin_axis(ax, ax_r)
    _apply_axis_limits(ax, [te_mean, te_max], frac=0.07)
    _apply_axis_limits(ax_r, [tph_mean, tph_max], frac=0.07)

    ax.legend(
        [
            Line2D([0], [0], color="black", linestyle="-"),
            Line2D([0], [0], color="black", linestyle="--"),
        ],
        ["mean", "max"],
        loc="lower right",
        frameon=False,
        ncol=2,
        columnspacing=1.0,
        handlelength=2.4,
        borderaxespad=0.4,
        fontsize=14,
    )

    # ------------------------------------------------------------------
    # 2) Powers
    # SOLO máximos, todos con línea sólida.
    # ------------------------------------------------------------------
    ax = axes[1]
    ax_r = ax.twinx()
    _shade_prethermal(ax, prethermal_t_ps)

    h_pj = _plot_single_series(ax, pj_t, np.abs(pj_max), color=pj_color)
    h_pd = _plot_single_series(ax, pdiff_t, np.abs(pdiff_max), color=pdiff_color)
    h_pep = _plot_single_series(ax_r, pep_t, np.abs(pep_max), color=pep_color)
    h_pesc = _plot_single_series(ax_r, pesc_t, np.abs(pesc_max), color=pesc_color)

    ax.set_xlabel("t [ps]", fontsize=14)
    ax.set_ylabel(r"$P_J,\ P_{diff}$ [W m$^{-3}$]", fontsize=14, color=hot_color)
    ax_r.set_ylabel(r"$P_{ep},\ P_{esc}$ [W m$^{-3}$]", fontsize=14, color=cold_color)
    ax.tick_params(axis="x", labelsize=14)
    ax.tick_params(axis="y", colors=hot_color, labelsize=14)
    ax_r.tick_params(axis="y", colors=cold_color, labelsize=14)

    ax.set_yscale("symlog", linthresh=_linthresh_from_axes_lines(ax))
    ax_r.set_yscale("symlog", linthresh=_linthresh_from_axes_lines(ax_r))

    _autoscale_symlog_axis(ax)
    _autoscale_symlog_axis(ax_r)

    # Forzar que el eje izquierdo llegue al menos hasta 1e12.
    y0, y1 = ax.get_ylim()
    yr0, yr1 = ax_r.get_ylim()
    ax.set_ylim(min(y0, 0.0), max(y1, 1.0e12))
    ax_r.set_ylim(min(yr0, 0.0), max(yr1,1.0e-7))

    ax.grid(False)
    ax_r.grid(False)
    _clean_twin_axis(ax, ax_r)

    ax.legend(
        [h_pj, h_pd, h_pep, h_pesc],
        [r"$P_J$", r"$P_{diff}$", r"$P_{ep}$", r"$P_{esc}$"],
        loc="lower right",
        frameon=False,
        ncol=2,
        columnspacing=1.1,
        handlelength=2.5,
        borderaxespad=0.3,
        fontsize=14,
    )

    # ------------------------------------------------------------------
    # 3) Energies
    # ------------------------------------------------------------------
    ax = axes[2]
    ax_r = ax.twinx()
    _shade_prethermal(ax, prethermal_t_ps)

    h_uph = _plot_single_series(ax, uph_t, uph_mean, color=uph_color)
    h_ue = _plot_single_series(ax_r, ue_t, ue_mean, color=ue_color)

    ax.set_xlabel("t [ps]", fontsize=14)
    ax.set_ylabel(r"$u_{ph}$ [J m$^{-3}$]", color=uph_color, fontsize=14)
    ax_r.set_ylabel(r"$u_e$ [J m$^{-3}$]", color=ue_color, fontsize=14)

    ax.tick_params(axis="x", labelsize=14)
    ax.tick_params(axis="y", colors=uph_color, labelsize=14)
    ax_r.tick_params(axis="y", colors=ue_color, labelsize=14)

    _center_linear_axis(ax, [uph_mean], frac=0.25)
    _center_linear_axis(ax_r, [ue_mean], frac=0.10)

    fmt_right = mticker.ScalarFormatter(useMathText=True)
    fmt_right.set_scientific(True)
    fmt_right.set_powerlimits((-3, 3))
    fmt_right.set_useOffset(True)
    ax_r.yaxis.set_major_formatter(fmt_right)
    ax_r.yaxis.set_major_locator(mticker.MaxNLocator(nbins=4))
    ax_r.yaxis.get_offset_text().set_color(ue_color)

    ax.grid(False)
    ax_r.grid(False)
    _clean_twin_axis(ax, ax_r)

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

from pysnspd.plotting.ss_memory_helpers import (
    _apply_axis_limits,
    _autoscale_symlog_axis,
    _center_linear_axis,
    _center_mask_from_x,
    _center_window,
    _clean_twin_axis,
    _final_current_family_arrays,
    _has_current_family_data,
    _has_field_map_data,
    _has_thermal_scalar_data,
    _history_or_snapshot_pair,
    _javg,
    _linthresh_from_axes_lines,
    _load_npz_if_exists,
    _overlay_current_arrows,
    _plot_series_pair,
    _plot_single_series,
    _prepare_output,
    _shade_prethermal,
    _snapshot_array,
    _snapshot_mean_series,
    _snapshot_times_ps,
    _triangulation,
)
