"""Power, energy, and snapshot figures for stationary SS runs.

The SS solver writes two optional post-processing files:

``stationary_snapshots.npz``
    Mesoscopic fields sampled only at requested physical snapshot times.

``snapshot_power_energy_diagnostics.npz``
    Runtime lookup of PRE power/energy/transport catalogues at the same
    snapshot times.

This module deliberately treats both files as diagnostics.  It never changes
solver state and it gracefully returns no figures for older runs that do not
have the new files yet.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

import numpy as np

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.colors import LogNorm, SymLogNorm

from pysnspd.analysis.snapshots import compute_snapshot_joule_power_density
from pysnspd.plotting.style import THESIS_DPI, THESIS_WIDTH_IN, apply_thesis_style

apply_thesis_style()

MEV_J = 1.602176634e-22


def make_ss_snapshot_power_figures(
    *,
    mesh: Any,
    dataset: Mapping[str, Any],
    raw_ss: str | Path,
    output_dir: str | Path,
    dpi: int = THESIS_DPI,
) -> dict[str, Path]:
    """Create additional snapshot/power figures for a completed SS run.

    Missing diagnostic files are not an error because historical SS runs only
    contain final-state fields.  The returned dictionary is directly merged
    into the plotting-pipeline manifest.
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    raw = Path(raw_ss)
    snapshots = _load_npz_if_exists(raw / "stationary_snapshots.npz")
    power = _load_npz_if_exists(raw / "snapshot_power_energy_diagnostics.npz")
    power = _with_recomputed_joule_power(power, snapshots=snapshots, dataset=dataset)

    saved: dict[str, Path] = {}
    if snapshots:
        saved["snapshot_state_atlas"] = plot_ss_snapshot_state_atlas(
            mesh,
            snapshots,
            dataset,
            out / "ss_snapshot_state_atlas.png",
            dpi=dpi,
        )
        saved["snapshot_final_profile_comparison"] = plot_ss_snapshot_profile_comparison(
            mesh,
            snapshots,
            power,
            dataset,
            out / "ss_snapshot_profile_comparison.png",
            dpi=dpi,
        )

    if power:
        saved["snapshot_power_energy_maps"] = plot_ss_snapshot_power_energy_maps(
            mesh,
            power,
            dataset,
            out / "ss_snapshot_power_energy_maps.png",
            dpi=dpi,
        )
        saved["snapshot_power_balance_maps"] = plot_ss_snapshot_power_balance_maps(
            mesh,
            power,
            dataset,
            out / "ss_snapshot_power_balance_maps.png",
            snapshots=snapshots,
            dpi=dpi,
        )
        saved["snapshot_runtime_metrics"] = plot_ss_snapshot_runtime_metrics(
            power,
            snapshots,
            dataset,
            out / "ss_snapshot_runtime_metrics.png",
            dpi=dpi,
        )
        saved["snapshot_catalog_indices"] = plot_ss_snapshot_catalog_indices(
            mesh,
            power,
            dataset,
            out / "ss_snapshot_catalog_indices.png",
            dpi=dpi,
        )

    return saved


# -----------------------------------------------------------------------------
# Public figure functions
# -----------------------------------------------------------------------------


def plot_ss_snapshot_state_atlas(
    mesh: Any,
    snapshots: Mapping[str, np.ndarray],
    dataset: Mapping[str, Any],
    output_path: str | Path,
    *,
    dpi: int = THESIS_DPI,
) -> Path:
    """Plot the dynamic SS state at all stored snapshots.

    Default rows are |Delta|/Delta0, phi, and |j|/javg.  When runtime thermal
    coupling was enabled and the solver saved Te/Tph snapshots, those are added
    as extra rows so the atlas directly shows whether the thermal window is
    actually evolving.
    """
    output = _prepare_output(output_path)
    tri = _triangulation(mesh, dataset)

    delta_mev = _snapshot_array(
        snapshots,
        ("delta_snapshot_meV",),
        fallback=_delta_mev_from_psi_snapshots(snapshots),
    )
    phi_v = _snapshot_array(snapshots, ("phi_snapshot_V",))
    jmag = _snapshot_current_mag(snapshots, family="jtot")
    Te = _snapshot_array(snapshots, ("Te_snapshot_K",), shape_like=delta_mev)
    Tph = _snapshot_array(snapshots, ("Tph_snapshot_K",), shape_like=delta_mev)

    if delta_mev.size == 0:
        raise ValueError("stationary_snapshots.npz does not contain Delta snapshots")

    t_ps = _snapshot_times_ps(
        snapshots, preferred=("snapshot_t_s", "delta_snapshot_t_s"), n=delta_mev.shape[0]
    )
    indices = _representative_snapshot_indices(delta_mev.shape[0], max_panels=9)
    t_sel = t_ps[indices]

    delta0 = _delta0_mev(dataset, snapshots)
    delta_norm = delta_mev / delta0 if delta0 > 0.0 else delta_mev
    phi_mv = 1.0e3 * _resize_snapshot_field(phi_v, delta_mev.shape)
    jscale = _javg(dataset, snapshots)
    j_norm = _resize_snapshot_field(jmag, delta_mev.shape) / jscale

    fields: list[tuple[np.ndarray, str, bool, float | None, float | None, bool, str]] = [
        (delta_norm[indices], r"$|\Delta|/\Delta_0$", False, 0.0, None, False, r"order parameter"),
        (phi_mv[indices], r"$\phi$ [mV]", True, None, None, True, r"electrostatic potential"),
        (j_norm[indices], r"$|j|/j_{avg}$", False, 0.0, None, False, r"total current"),
    ]
    if Te.size:
        fields.append(
            (_resize_snapshot_field(Te, delta_mev.shape)[indices], r"$T_e$ [K]", False, None, None, False, r"electron temperature")
        )
    if Tph.size:
        fields.append(
            (_resize_snapshot_field(Tph, delta_mev.shape)[indices], r"$T_{ph}$ [K]", False, None, None, False, r"phonon temperature")
        )

    fig, axes = _snapshot_grid_figure(
        nrows=len(fields),
        ncols=len(indices),
        title=f"SS snapshots: {dataset.get('run_name', '')}",
        left=0.050,
        right=0.985,
        bottom=0.060,
        top=0.885,
        wspace=0.10,
        hspace=0.34,
    )

    for row, (z, label, symmetric, vmin, vmax, wrap_first, row_title) in enumerate(fields):
        scale_values = z[1:] if wrap_first and z.shape[0] > 1 else z
        norm, vmin_eff, vmax_eff = _node_color_limits(
            scale_values, symmetric=symmetric, vmin=vmin, vmax=vmax
        )
        z_plot = np.array(z, copy=True)
        if wrap_first and z_plot.shape[0] > 0:
            z_plot[0] = _wrap_values_to_range(z_plot[0], vmin_eff, vmax_eff)

        mappable = None
        for col, _idx in enumerate(indices):
            ax = axes[row, col]
            mappable = ax.tripcolor(
                tri, z_plot[col], shading="gouraud", vmin=vmin_eff, vmax=vmax_eff, norm=norm
            )
            if row == 0:
                _annotate_snapshot_time(ax, t_sel[col])
            if col == 0:
                ax.set_ylabel(row_title)
            _format_map_axis(ax, show_xlabel_top=(row == 0), show_ylabel=(col == 0))
        if mappable is not None:
            _add_row_colorbar(fig, axes[row, :], mappable, label)

    fig.savefig(output, dpi=dpi, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    return output


def plot_ss_snapshot_power_energy_maps(
    mesh: Any,
    power: Mapping[str, np.ndarray],
    dataset: Mapping[str, Any],
    output_path: str | Path,
    *,
    dpi: int = THESIS_DPI,
) -> Path:
    """Plot catalogue-derived power, Joule, escape, and transport maps."""
    output = _prepare_output(output_path)
    tri = _triangulation(mesh, dataset)

    p_ep = _snapshot_array(power, ("P_total_snapshot_W_m3",))
    joule = _snapshot_array(power, ("joule_snapshot_W_m3",), shape_like=p_ep)
    p_esc = _snapshot_array(power, ("P_esc_snapshot_W_m3",), shape_like=p_ep)
    kappa = _snapshot_array(power, ("kappa_s_snapshot_W_m_K",), shape_like=p_ep)

    if p_ep.size == 0:
        raise ValueError("snapshot_power_energy_diagnostics.npz lacks P_total_snapshot_W_m3")

    t_ps = _snapshot_times_ps(power, preferred=("snapshot_t_s",), n=p_ep.shape[0])
    indices = _representative_snapshot_indices(p_ep.shape[0], max_panels=9)
    t_sel = t_ps[indices]

    fields = [
        (p_ep[indices], r"$P_{ep}=P_S+P_R$ [W m$^{-3}$]", "signed", r"electron $\rightarrow$ phonon power"),
        (joule[indices], r"$P_J=|j_n|^2/\sigma_n$ [W m$^{-3}$]", "positive_log", r"Joule diagnostic"),
        (p_esc[indices], r"$P_{esc}$ [W m$^{-3}$]", "signed", r"phonon escape power"),
        (kappa[indices], r"$\kappa_s$ [W m$^{-1}$ K$^{-1}$]", "positive_log", r"thermal conductivity"),
    ]

    fig, axes = _snapshot_grid_figure(
        nrows=len(fields),
        ncols=len(indices),
        title=f"SS runtime power/transport maps: {dataset.get('run_name', '')}",
        left=0.050,
        right=0.985,
        bottom=0.060,
        top=0.860,
        wspace=0.10,
        hspace=0.36,
    )

    for row, (z, label, mode, row_title) in enumerate(fields):
        norm, vmin_eff, vmax_eff = _norm_for_mode(z, mode)
        mappable = None
        for col, _idx in enumerate(indices):
            ax = axes[row, col]
            z_panel = _plot_values_for_mode(z[col], mode=mode, norm=norm)
            mappable = ax.tripcolor(
                tri, z_panel, shading="gouraud", vmin=vmin_eff, vmax=vmax_eff, norm=norm
            )
            if row == 0:
                _annotate_snapshot_time(ax, t_sel[col])
            if col == 0:
                ax.set_ylabel(row_title)
            _format_map_axis(ax, show_xlabel_top=(row == 0), show_ylabel=(col == 0))
        if mappable is not None:
            _add_row_colorbar(fig, axes[row, :], mappable, label)

    fig.savefig(output, dpi=dpi, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    return output


def plot_ss_snapshot_power_balance_maps(
    mesh: Any,
    power: Mapping[str, np.ndarray],
    dataset: Mapping[str, Any],
    output_path: str | Path,
    *,
    snapshots: Mapping[str, np.ndarray] | None = None,
    dpi: int = THESIS_DPI,
) -> Path:
    """Plot diagnostic local energy-balance tendencies at snapshots.

    The electronic tendency is now the no-photon thermal-balance diagnostic

        P_J + P_diff - P_ep,

    where ``P_diff`` is the finite-volume graph approximation to
    ``div(kappa_s grad T_e)`` when a saved diffusion map is not already present
    in ``snapshot_power_energy_diagnostics.npz``.  The phonon tendency remains
    ``P_ep - P_esc``.
    """
    output = _prepare_output(output_path)
    tri = _triangulation(mesh, dataset)

    p_ep = _snapshot_array(power, ("P_total_snapshot_W_m3",))
    joule = _snapshot_array(power, ("joule_snapshot_W_m3",), shape_like=p_ep)
    p_esc = _snapshot_array(power, ("P_esc_snapshot_W_m3",), shape_like=p_ep)
    p_diff = _snapshot_diffusion_power_density(
        mesh,
        snapshots=snapshots,
        power=power,
        dataset=dataset,
        shape_like=p_ep,
    )

    if p_ep.size == 0:
        raise ValueError("snapshot_power_energy_diagnostics.npz lacks P_total_snapshot_W_m3")

    electron_balance = joule + p_diff - p_ep
    phonon_balance = p_ep - p_esc
    joule_plus_diff = joule + p_diff

    t_ps = _snapshot_times_ps(power, preferred=("snapshot_t_s",), n=p_ep.shape[0])
    indices = _representative_snapshot_indices(p_ep.shape[0], max_panels=9)
    t_sel = t_ps[indices]

    fields = [
        (joule[indices], r"$P_J$ [W m$^{-3}$]", "positive_log", r"Joule heating"),
        (p_diff[indices], r"$P_{diff}=\nabla\cdot(\kappa_s\nabla T_e)$ [W m$^{-3}$]", "signed", r"electron diffusion"),
        (joule_plus_diff[indices], r"$P_J+P_{diff}$ [W m$^{-3}$]", "signed", r"Joule + diffusion"),
        (electron_balance[indices], r"$P_J+P_{diff}-P_{ep}$ [W m$^{-3}$]", "signed", r"electronic tendency"),
        (phonon_balance[indices], r"$P_{ep}-P_{esc}$ [W m$^{-3}$]", "signed", r"phonon tendency"),
    ]

    fig, axes = _snapshot_grid_figure(
        nrows=len(fields),
        ncols=len(indices),
        title=f"SS diagnostic power-balance maps: {dataset.get('run_name', '')}",
        left=0.050,
        right=0.985,
        bottom=0.050,
        top=0.840,
        wspace=0.10,
        hspace=0.39,
    )

    for row, (z, label, mode, row_title) in enumerate(fields):
        norm, vmin_eff, vmax_eff = _norm_for_mode(z, mode)
        mappable = None
        for col, _idx in enumerate(indices):
            ax = axes[row, col]
            z_panel = _plot_values_for_mode(z[col], mode=mode, norm=norm)
            mappable = ax.tripcolor(
                tri, z_panel, shading="gouraud", vmin=vmin_eff, vmax=vmax_eff, norm=norm
            )
            if row == 0:
                _annotate_snapshot_time(ax, t_sel[col])
            if col == 0:
                ax.set_ylabel(row_title)
            _format_map_axis(ax, show_xlabel_top=(row == 0), show_ylabel=(col == 0))
        if mappable is not None:
            _add_row_colorbar(fig, axes[row, :], mappable, label)

    fig.savefig(output, dpi=dpi, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    return output


def plot_ss_snapshot_runtime_metrics(
    power: Mapping[str, np.ndarray],
    snapshots: Mapping[str, np.ndarray] | None,
    dataset: Mapping[str, Any],
    output_path: str | Path,
    *,
    dpi: int = THESIS_DPI,
) -> Path:
    """Plot compact scalar metrics extracted from snapshot maps and history."""
    output = _prepare_output(output_path)

    p_ep = _snapshot_array(power, ("P_total_snapshot_W_m3",))
    joule = _snapshot_array(power, ("joule_snapshot_W_m3",), shape_like=p_ep)
    p_esc = _snapshot_array(power, ("P_esc_snapshot_W_m3",), shape_like=p_ep)
    q_abs = _snapshot_array(power, ("q_abs_snapshot_m_inv",), shape_like=p_ep)
    u_e = _snapshot_array(power, ("u_e_snapshot_J_m3",), shape_like=p_ep)
    u_ph = _snapshot_array(power, ("u_ph_snapshot_J_m3",), shape_like=p_ep)
    c_e = _snapshot_array(power, ("C_e_snapshot_J_m3_K",), shape_like=p_ep)
    c_ph = _snapshot_array(power, ("C_ph_snapshot_J_m3_K",), shape_like=p_ep)

    if p_ep.size == 0:
        raise ValueError("snapshot_power_energy_diagnostics.npz lacks P_total_snapshot_W_m3")

    p_diff = _snapshot_diffusion_power_density(
        None,
        snapshots=snapshots,
        power=power,
        dataset=dataset,
        shape_like=p_ep,
    )

    t_ps_snap = _snapshot_times_ps(power, preferred=("snapshot_t_s",), n=p_ep.shape[0])
    t_ps_hist = np.asarray(dataset.get("t_ps", []), dtype=float)
    use_hist = t_ps_hist.size > 0

    delta_norm = np.empty((0, 0), dtype=float)
    te_snap = np.empty((0, 0), dtype=float)
    tph_snap = np.empty((0, 0), dtype=float)
    if snapshots:
        delta_mev = _snapshot_array(
            snapshots,
            ("delta_snapshot_meV",),
            fallback=_delta_mev_from_psi_snapshots(snapshots),
        )
        if delta_mev.size:
            delta0 = _delta0_mev(dataset, snapshots)
            delta_norm = delta_mev / delta0 if delta0 > 0.0 else delta_mev
        te_snap = _snapshot_array(snapshots, ("Te_snapshot_K",), shape_like=p_ep)
        tph_snap = _snapshot_array(snapshots, ("Tph_snapshot_K",), shape_like=p_ep)

    thermal_hist_present = any(
        np.asarray(dataset.get(key, []), dtype=float).size > 0
        for key in (
            "thermal_max_Te_K_history",
            "thermal_mean_Te_K_history",
            "thermal_max_rate_K_per_ps_history",
            "thermal_max_P_J_W_m3_history",
        )
    )

    nrows = 3 if thermal_hist_present or te_snap.size or tph_snap.size else 2
    fig, axes = plt.subplots(
        nrows,
        2,
        figsize=(THESIS_WIDTH_IN, 3.2 * nrows),
        constrained_layout=False,
    )
    fig.subplots_adjust(left=0.090, right=0.970, bottom=0.080, top=0.925, wspace=0.30, hspace=0.36)
    fig.suptitle(f"SS snapshot/runtime diagnostics: {dataset.get('run_name', '')}", y=0.975)

    ax = axes[0, 0]
    _plot_snapshot_metric(ax, t_ps_snap, p_ep, reducer="max_abs", label=r"max $|P_{ep}|$")
    _plot_snapshot_metric(ax, t_ps_snap, joule, reducer="max", label=r"max $P_J$")
    _plot_snapshot_metric(ax, t_ps_snap, joule + p_diff - p_ep, reducer="max_abs", label=r"max $|P_J+P_{diff}-P_{ep}|$")
    if p_esc.size:
        _plot_snapshot_metric(ax, t_ps_snap, p_esc, reducer="max_abs", label=r"max $|P_{esc}|$")
    if np.any(np.isfinite(p_diff)):
        _plot_snapshot_metric(ax, t_ps_snap, p_diff, reducer="max_abs", label=r"max $|P_{diff}|$")
    ax.set_yscale("symlog", linthresh=1.0e8)
    ax.set_ylabel(r"power density [W m$^{-3}$]")
    ax.set_title("power scales")
    ax.grid(False)
    ax.legend(frameon=False)

    ax = axes[0, 1]
    _plot_snapshot_metric(ax, t_ps_snap, q_abs / 1.0e7, reducer="p99", label=r"p99 $|q|$")
    _plot_snapshot_metric(ax, t_ps_snap, q_abs / 1.0e7, reducer="max", label=r"max $|q|$")
    if delta_norm.size:
        _plot_snapshot_metric(ax, t_ps_snap, delta_norm, reducer="min", label=r"min $|\Delta|/\Delta_0$")
        _plot_snapshot_metric(ax, t_ps_snap, delta_norm, reducer="mean", label=r"mean $|\Delta|/\Delta_0$")
    ax.set_ylabel(r"$q$ [$10^7$ m$^{-1}$] or normalized gap")
    ax.set_title("order parameter / momentum")
    ax.grid(False)
    ax.legend(frameon=False)

    row_offset = 1
    if nrows == 3:
        ax = axes[1, 0]
        if thermal_hist_present and use_hist:
            _plot_series_if_any(ax, t_ps_hist, dataset.get("thermal_mean_Te_K_history"), label=r"mean $T_e$")
            _plot_series_if_any(ax, t_ps_hist, dataset.get("thermal_max_Te_K_history"), label=r"max $T_e$")
            _plot_series_if_any(ax, t_ps_hist, dataset.get("thermal_mean_Tph_K_history"), label=r"mean $T_{ph}$")
            _plot_series_if_any(ax, t_ps_hist, dataset.get("thermal_max_Tph_K_history"), label=r"max $T_{ph}$")
            ax.set_xlabel("t [ps]")
        else:
            _plot_snapshot_metric(ax, t_ps_snap, te_snap, reducer="mean", label=r"mean $T_e$")
            _plot_snapshot_metric(ax, t_ps_snap, te_snap, reducer="max", label=r"max $T_e$")
            _plot_snapshot_metric(ax, t_ps_snap, tph_snap, reducer="mean", label=r"mean $T_{ph}$")
            _plot_snapshot_metric(ax, t_ps_snap, tph_snap, reducer="max", label=r"max $T_{ph}$")
            ax.set_xlabel("t [ps]")
        ax.set_ylabel("temperature [K]")
        ax.set_title("thermal-window temperatures")
        ax.grid(False)
        ax.legend(frameon=False)

        ax = axes[1, 1]
        if thermal_hist_present and use_hist:
            _plot_series_if_any(ax, t_ps_hist, dataset.get("thermal_max_rate_K_per_ps_history"), label=r"max $|dT/dt|$")
            _plot_series_if_any(ax, t_ps_hist, dataset.get("thermal_substeps_history"), label="thermal substeps")
            ax2 = ax.twinx()
            _plot_series_if_any(ax2, t_ps_hist, dataset.get("thermal_max_abs_dTe_K_history"), label=r"max $|\Delta T_e|$")
            _plot_series_if_any(ax2, t_ps_hist, dataset.get("thermal_max_abs_dTph_K_history"), label=r"max $|\Delta T_{ph}|$")
            ax.set_xlabel("t [ps]")
            ax.set_ylabel(r"rate [K ps$^{-1}$] / substeps")
            ax2.set_ylabel(r"per-step increment [K]")
            _legend_if_labels(ax, frameon=False, loc="upper left")
            _legend_if_labels(ax2, frameon=False, loc="upper right")
        else:
            _plot_snapshot_metric(ax, t_ps_snap, te_snap, reducer="max", label=r"max $T_e$")
            _plot_snapshot_metric(ax, t_ps_snap, tph_snap, reducer="max", label=r"max $T_{ph}$")
            ax.set_xlabel("t [ps]")
            ax.set_ylabel("temperature [K]")
            ax.legend(frameon=False)
        ax.set_title("thermal step diagnostics")
        ax.grid(False)
        row_offset = 2

    ax = axes[row_offset, 0]
    _plot_snapshot_metric(ax, t_ps_snap, u_e, reducer="mean", label=r"mean $u_e$")
    _plot_snapshot_metric(ax, t_ps_snap, u_ph, reducer="mean", label=r"mean $u_{ph}$")
    ax.set_ylabel(r"energy density [J m$^{-3}$]")
    ax.set_xlabel("t [ps]")
    ax.set_title("stored energy diagnostics")
    ax.grid(False)
    ax.legend(frameon=False)

    ax = axes[row_offset, 1]
    if thermal_hist_present and use_hist:
        _plot_series_if_any(ax, t_ps_hist, dataset.get("thermal_max_P_J_W_m3_history"), label=r"max $P_J$")
        _plot_series_if_any(ax, t_ps_hist, dataset.get("thermal_max_P_ep_W_m3_history"), label=r"max $P_{ep}$")
        _plot_series_if_any(ax, t_ps_hist, dataset.get("thermal_max_P_esc_W_m3_history"), label=r"max $P_{esc}$")
        _plot_series_if_any(ax, t_ps_hist, dataset.get("thermal_max_P_diff_W_m3_history"), label=r"max $P_{diff}$")
        ax.set_xlabel("t [ps]")
        ax.set_ylabel(r"power density [W m$^{-3}$]")
        ax.set_yscale("symlog", linthresh=1.0e8)
        ax.set_title("runtime thermal power envelopes")
        ax.grid(False)
        ax.legend(frameon=False)
    else:
        _plot_snapshot_metric(ax, t_ps_snap, c_e, reducer="mean", label=r"mean $C_e$")
        _plot_snapshot_metric(ax, t_ps_snap, c_ph, reducer="mean", label=r"mean $C_{ph}$")
        ax.set_yscale("log")
        ax.set_ylabel(r"heat capacity [J m$^{-3}$ K$^{-1}$]")
        ax.set_xlabel("t [ps]")
        ax.set_title("heat-capacity diagnostics")
        ax.grid(False)
        ax.legend(frameon=False)

    fig.savefig(output, dpi=dpi, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    return output


def plot_ss_snapshot_catalog_indices(
    mesh: Any,
    power: Mapping[str, np.ndarray],
    dataset: Mapping[str, Any],
    output_path: str | Path,
    *,
    dpi: int = THESIS_DPI,
) -> Path:
    """Visualize which PRE-table bins are being queried at the final snapshot."""
    output = _prepare_output(output_path)
    tri = _triangulation(mesh, dataset)

    keys = [
        ("power_table_iTe", r"$i_{T_e}$"),
        ("power_table_iTph", r"$i_{T_{ph}}$"),
        ("power_table_iDelta", r"$i_{\Delta}$"),
        ("power_table_iQ", r"$i_q$"),
    ]
    fields = []
    for key, label in keys:
        arr = _snapshot_array(power, (key,))
        if arr.size:
            fields.append((arr[-1], label, key))
    if not fields:
        raise ValueError("snapshot_power_energy_diagnostics.npz lacks catalogue index maps")

    fig, axes = plt.subplots(
        1,
        len(fields),
        figsize=(THESIS_WIDTH_IN, 3.2),
        constrained_layout=False,
    )
    axes = np.asarray(axes, dtype=object).reshape(1, -1)
    fig.subplots_adjust(left=0.055, right=0.965, bottom=0.145, top=0.845, wspace=0.32)
    fig.suptitle(f"PRE catalogue indices at final SS snapshot: {dataset.get('run_name', '')}", y=0.960)

    flat_axes = axes.ravel()
    for ax, (z, label, key) in zip(flat_axes, fields):
        z = np.asarray(z, dtype=float)
        finite = z[np.isfinite(z)]
        vmin = float(np.nanmin(finite)) if finite.size else 0.0
        vmax = float(np.nanmax(finite)) if finite.size else 1.0
        if vmax <= vmin:
            vmax = vmin + 1.0
        mappable = ax.tripcolor(tri, z, shading="gouraud", vmin=vmin, vmax=vmax)
        ax.set_title(label)
        _format_map_axis(ax, show_xlabel=True, show_ylabel=(ax is flat_axes[0]))
        cb = fig.colorbar(mappable, ax=ax, shrink=0.82)
        cb.set_label(key)

    fig.savefig(output, dpi=dpi, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    return output


def plot_ss_snapshot_profile_comparison(
    mesh: Any,
    snapshots: Mapping[str, np.ndarray],
    power: Mapping[str, np.ndarray] | None,
    dataset: Mapping[str, Any],
    output_path: str | Path,
    *,
    dpi: int = THESIS_DPI,
) -> Path:
    """Plot x-binned profiles at first/middle/final snapshot.

    When Te/Tph snapshots exist, they are appended as extra profile rows so the
    same figure can be used to judge both mesoscopic and thermal relaxation.
    """
    output = _prepare_output(output_path)
    x_nm = _mesh_x_nm(mesh, dataset)

    delta_mev = _snapshot_array(
        snapshots,
        ("delta_snapshot_meV",),
        fallback=_delta_mev_from_psi_snapshots(snapshots),
    )
    if delta_mev.size == 0:
        raise ValueError("stationary_snapshots.npz does not contain Delta snapshots")

    delta0 = _delta0_mev(dataset, snapshots)
    delta_norm = delta_mev / delta0 if delta0 > 0.0 else delta_mev
    phi = 1.0e3 * _resize_snapshot_field(_snapshot_array(snapshots, ("phi_snapshot_V",)), delta_mev.shape)
    j_norm = _resize_snapshot_field(_snapshot_current_mag(snapshots, family="jtot"), delta_mev.shape) / _javg(dataset, snapshots)
    te = _resize_snapshot_field(_snapshot_array(snapshots, ("Te_snapshot_K",)), delta_mev.shape)
    tph = _resize_snapshot_field(_snapshot_array(snapshots, ("Tph_snapshot_K",)), delta_mev.shape)

    q_abs = np.empty((0, 0), dtype=float)
    p_balance = np.empty((0, 0), dtype=float)
    if power:
        p_ep = _snapshot_array(power, ("P_total_snapshot_W_m3",))
        joule = _snapshot_array(power, ("joule_snapshot_W_m3",), shape_like=p_ep)
        q_abs = _snapshot_array(power, ("q_abs_snapshot_m_inv",), shape_like=p_ep) / 1.0e7
        p_diff = _snapshot_diffusion_power_density(mesh, snapshots=snapshots, power=power, dataset=dataset, shape_like=p_ep)
        p_balance = joule + p_diff - p_ep

    t_ps = _snapshot_times_ps(
        snapshots, preferred=("snapshot_t_s", "delta_snapshot_t_s"), n=delta_mev.shape[0]
    )
    indices = _representative_snapshot_indices(delta_mev.shape[0], max_panels=3)

    rows: list[tuple[str, np.ndarray]] = [
        (r"$|\Delta|/\Delta_0$", delta_norm),
        (r"$|j|/j_{avg}$", j_norm),
        (r"$\phi$ [mV]", phi),
    ]
    if np.any(np.isfinite(te)) and np.nanmax(np.abs(te)) > 0.0:
        rows.append((r"$T_e$ [K]", te))
    if np.any(np.isfinite(tph)) and np.nanmax(np.abs(tph)) > 0.0:
        rows.append((r"$T_{ph}$ [K]", tph))
    if q_abs.size:
        rows.append((r"$|q|$ [$10^7$ m$^{-1}$]", q_abs))
    elif p_balance.size:
        rows.append((r"$P_J+P_{diff}-P_{ep}$", p_balance))

    fig, axes = plt.subplots(
        len(rows),
        1,
        figsize=(THESIS_WIDTH_IN, 2.0 * len(rows) + 0.8),
        sharex=True,
        constrained_layout=False,
    )
    axes = np.atleast_1d(axes)
    fig.subplots_adjust(left=0.105, right=0.970, bottom=0.070, top=0.945, hspace=0.30)
    fig.suptitle(f"SS x-profile evolution: {dataset.get('run_name', '')}", y=0.985)

    for idx in indices:
        label = f"t={t_ps[idx]:.3g} ps"
        for ax, (_ylabel, values) in zip(axes, rows):
            _plot_binned_profile(ax, x_nm, values[idx], label=label)
    for ax, (ylabel, _values) in zip(axes, rows):
        ax.set_ylabel(ylabel)
        ax.grid(False)
        ax.legend(frameon=False, loc="best")
    axes[-1].set_xlabel("x [nm]")

    fig.savefig(output, dpi=dpi, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    return output


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

from pysnspd.plotting.ss_power_helpers import (
    _add_row_colorbar,
    _annotate_snapshot_time,
    _delta0_mev,
    _delta_mev_from_psi_snapshots,
    _format_map_axis,
    _javg,
    _legend_if_labels,
    _load_npz_if_exists,
    _mesh_x_nm,
    _node_color_limits,
    _norm_for_mode,
    _plot_binned_profile,
    _plot_series_if_any,
    _plot_snapshot_metric,
    _plot_values_for_mode,
    _prepare_output,
    _representative_snapshot_indices,
    _resize_snapshot_field,
    _snapshot_array,
    _snapshot_current_mag,
    _snapshot_diffusion_power_density,
    _snapshot_grid_figure,
    _snapshot_times_ps,
    _triangulation,
    _with_recomputed_joule_power,
    _wrap_values_to_range,
)
