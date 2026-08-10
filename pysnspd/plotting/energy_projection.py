"""Plots for the temporary D3 Simon-energy projection diagnostic."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping, Sequence

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import Normalize, SymLogNorm
from matplotlib.ticker import MaxNLocator
import numpy as np

from pysnspd.plotting.style import THESIS_DPI, THESIS_WIDTH_IN, apply_thesis_style

apply_thesis_style()

EV_J = 1.602176634e-19


def write_energy_projection_figures(
    dataset: Mapping[str, Any],
    output_dir: str | Path,
    *,
    dpi: int = THESIS_DPI,
) -> dict[str, Path]:
    """Create the three requested D3 PDF products."""

    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    saved = {
        "colormaps": plot_energy_projection_colormaps(
            dataset,
            output / "D3_energy_projection_colormaps.pdf",
            dpi=dpi,
        ),
        "temporal": plot_energy_projection_temporal(
            dataset,
            output / "D3_energy_projection_temporal.pdf",
            dpi=dpi,
        ),
        "profiles": plot_energy_projection_profiles(
            dataset,
            output / "D3_energy_projection_longitudinal_profiles.pdf",
            dpi=dpi,
        ),
    }
    return saved


def plot_energy_projection_colormaps(
    dataset: Mapping[str, Any],
    output_path: str | Path,
    *,
    dpi: int = THESIS_DPI,
) -> Path:
    """Write a multi-page atlas of state, energy, power and closure maps."""

    output = _pdf_path(output_path)
    groups = (
        (
            "Persisted state used by the projection",
            (
                ("delta_over_delta0", r"$|\Delta|/\Delta_0$", "viridis", False),
                ("q_xi", r"$|q|\xi$", "magma", False),
                ("Te_K", r"$T_e$ [K]", "inferno", False),
                ("Tph_K", r"$T_{ph}$ [K]", "cividis", False),
            ),
        ),
        (
            "Electronic ledger and phonon energy",
            (
                ("u_cond_J_m3", r"$u_{cond}$ [J m$^{-3}$]", "coolwarm", True),
                ("u_qp_J_m3", r"$u_{qp}$ [J m$^{-3}$]", "magma", False),
                ("u_e_J_m3", r"$u_e$ [J m$^{-3}$]", "coolwarm", True),
                ("u_ph_J_m3", r"$u_{ph}$ [J m$^{-3}$]", "cividis", False),
            ),
        ),
        (
            "Motion of the superconducting spectrum at fixed electron temperature",
            (
                ("P_delta_W_m3", r"$P_\Delta$ [W m$^{-3}$]", "coolwarm", True),
                ("P_delta_cond_W_m3", r"$P_{\Delta,cond}$ [W m$^{-3}$]", "coolwarm", True),
                ("P_delta_qp_W_m3", r"$P_{\Delta,qp}$ [W m$^{-3}$]", "coolwarm", True),
                ("P_q_W_m3", r"$P_q$ (spectral only) [W m$^{-3}$]", "coolwarm", True),
            ),
        ),
        (
            "Omitted storage versus powers retained by the current temperature update",
            (
                ("P_spec_W_m3", r"$P_{spec}$ [W m$^{-3}$]", "coolwarm", True),
                ("P_J_W_m3", r"$P_J$ [W m$^{-3}$]", "magma", False),
                ("minus_P_ep_W_m3", r"$-P_{e-ph}$ [W m$^{-3}$]", "coolwarm", True),
                ("P_diff_W_m3", r"$P_{diff}$ [W m$^{-3}$]", "PuOr_r", True),
            ),
        ),
        (
            "Energy closure and finite-step sensitivity",
            (
                ("Q_ret_W_m3", r"$Q_{ret}=P_J+P_{diff}-P_{e-ph}$", "coolwarm", True),
                ("du_e_dt_W_m3", r"$\Delta u_e/\Delta t$", "coolwarm", True),
                ("residual_W_m3", r"$R_e=\Delta u_e/\Delta t-Q_{ret}$", "coolwarm", True),
                ("P_path_W_m3", r"Path ambiguity in $P_\Delta$", "coolwarm", True),
            ),
        ),
        (
            "Temperature-coordinate update and catalogue support",
            (
                ("P_thermal_W_m3", r"$P_T$ from exact $u_e$ difference", "coolwarm", True),
                ("P_CdT_W_m3", r"$\bar C_e\Delta T_e/\Delta t$", "coolwarm", True),
                ("P_T_nonlinear_W_m3", r"$P_T-\bar C_e\Delta T_e/\Delta t$", "coolwarm", True),
                ("catalog_clipped", "Catalogue clipping mask", "gray_r", False),
            ),
        ),
    )
    with PdfPages(output) as pdf:
        for title, fields in groups:
            figure = _map_page(dataset, title=title, fields=fields)
            pdf.savefig(figure, dpi=dpi)
            plt.close(figure)
    return output


def plot_energy_projection_temporal(
    dataset: Mapping[str, Any],
    output_path: str | Path,
    *,
    dpi: int = THESIS_DPI,
) -> Path:
    """Plot four central-window temporal diagnostics."""

    output = _pdf_path(output_path)
    time_ps = _array(dataset, "time_ps")
    if time_ps.size < 1:
        raise ValueError("D3 temporal diagnostics require at least one interval.")
    fig, axes = plt.subplots(2, 2, figsize=(THESIS_WIDTH_IN, 5.75), sharex=True)
    fig.subplots_adjust(left=0.105, right=0.975, bottom=0.095, top=0.925, wspace=0.28, hspace=0.24)

    axis = axes[0, 0]
    for key, label, color, style in (
        ("p99_abs_P_spec_W_m3", r"$P_{spec}$", "black", "-"),
        ("p99_abs_P_delta_W_m3", r"$P_\Delta$", "tab:blue", "-"),
        ("p99_abs_P_delta_cond_W_m3", r"$P_{\Delta,cond}$", "tab:cyan", "--"),
        ("p99_abs_P_q_W_m3", r"$P_q$", "tab:orange", "-"),
        ("p99_abs_Q_ret_W_m3", r"$Q_{ret}$", "tab:green", "-"),
        ("p99_abs_residual_W_m3", r"$R_e$", "tab:red", ":"),
    ):
        axis.plot(time_ps, np.maximum(_array(dataset, key), 1.0e-30), color=color, linestyle=style, linewidth=0.9, label=label)
    axis.set_yscale("log")
    axis.set_ylabel(r"Central p99 $|P|$ [W m$^{-3}$]")
    axis.set_title("Local energetic scales")
    axis.legend(frameon=False, fontsize=6.6, ncol=2)

    axis = axes[0, 1]
    signed_power_series = []
    for key, label, color in (
        ("integrated_P_spec_W", r"$P_{spec}$", "black"),
        ("integrated_P_delta_W", r"$P_\Delta$", "tab:blue"),
        ("integrated_P_q_W", r"$P_q$", "tab:orange"),
        ("integrated_Q_ret_W", r"$Q_{ret}$", "tab:green"),
        ("integrated_residual_W", r"$R_e$", "tab:red"),
    ):
        values_nW = 1.0e9 * _array(dataset, key)
        signed_power_series.append(values_nW)
        axis.plot(time_ps, values_nW, color=color, linewidth=0.9, label=label)
    signed_limit = max(
        float(np.nanmax(np.abs(values))) for values in signed_power_series
    )
    signed_limit = max(signed_limit, 1.0e-12)
    axis.set_yscale("symlog", linthresh=max(1.0e-4 * signed_limit, 1.0e-12))
    signed_ticks = [-signed_limit, -0.1 * signed_limit, 0.0, 0.1 * signed_limit, signed_limit]
    axis.set_yticks(signed_ticks)
    axis.set_yticklabels([f"{value:.2g}" for value in signed_ticks])
    axis.set_ylabel("Central integrated power [nW]")
    axis.set_title("Signed central power")
    axis.legend(frameon=False, fontsize=6.6, ncol=2)

    axis = axes[1, 0]
    dt_s = 1.0e-12 * _array(dataset, "dt_ps")
    for key, label, color, style in (
        ("integrated_P_spec_W", r"$E_{spec}$", "black", "-"),
        ("integrated_P_delta_W", r"$E_\Delta$", "tab:blue", "-"),
        ("integrated_P_delta_cond_W", r"$E_{\Delta,cond}$", "tab:cyan", "--"),
        ("integrated_P_q_W", r"$E_q$", "tab:orange", "-"),
        ("integrated_P_J_W", r"$E_J$", "tab:red", ":"),
        ("integrated_minus_P_ep_W", r"$-E_{e-ph}$", "tab:purple", ":"),
    ):
        energy_eV = np.cumsum(_array(dataset, key) * dt_s) / EV_J
        axis.plot(time_ps, energy_eV, color=color, linestyle=style, linewidth=0.9, label=label)
    axis.set_ylabel("Cumulative central energy [eV]")
    axis.set_xlabel(r"$t$ [ps]")
    axis.set_title("Integrated energetic impact")
    axis.legend(frameon=False, fontsize=6.4, ncol=2)

    axis = axes[1, 1]
    Tc_K = _scalar(dataset, "Tc_K")
    voltage = _array(dataset, "V_tdgl_V")
    voltage_scale = max(float(np.nanmax(np.abs(voltage))), 1.0e-300)
    axis.plot(time_ps, _array(dataset, "min_delta_over_delta0"), color="tab:green", label=r"$\min |\Delta|/\Delta_0$")
    axis.plot(time_ps, _array(dataset, "max_Te_K") / Tc_K, color="tab:red", label=r"$\max T_e/T_c$")
    axis.plot(time_ps, _array(dataset, "max_q_xi"), color="tab:orange", label=r"$\max |q|\xi$")
    axis.plot(time_ps, voltage / voltage_scale, color="tab:blue", linestyle="--", label=r"$V_{TDGL}/\max|V_{TDGL}|$")
    axis.plot(time_ps, _array(dataset, "catalog_clipped_fraction"), color="0.35", linestyle=":", label="any catalogue clipping")
    axis.plot(time_ps, _array(dataset, "catalog_q_clipped_fraction"), color="tab:purple", linestyle="-.", label=r"$q$ catalogue clipping")
    axis.set_ylabel("Dimensionless diagnostic")
    axis.set_xlabel(r"$t$ [ps]")
    axis.set_title("Detector variables and catalogue coverage")
    axis.legend(frameon=False, fontsize=6.2, ncol=2)

    for axis in axes.flat:
        _event_lines(axis, dataset)
        axis.grid(True)
        axis.set_xlim(float(time_ps[0]), float(time_ps[-1]))
    fig.suptitle(_run_scope_title(dataset), y=0.985, fontsize=9.5)
    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    return output


def plot_energy_projection_profiles(
    dataset: Mapping[str, Any],
    output_path: str | Path,
    *,
    dpi: int = THESIS_DPI,
) -> Path:
    """Plot cross-width integrated longitudinal profiles at three times."""

    output = _pdf_path(output_path)
    groups = (
        (
            "State profiles in the central diagnostic window",
            (
                ("delta_over_delta0", r"$|\Delta|/\Delta_0$", False),
                ("q_xi", r"$|q|\xi$", False),
                ("Te_K", r"$T_e$ [K]", False),
                ("Tph_K", r"$T_{ph}$ [K]", False),
            ),
        ),
        (
            "Energy changes relative to the initial persisted state",
            (
                ("u_cond_J_m3", r"$\Delta u_{cond}$ [J m$^{-3}$]", True),
                ("u_qp_J_m3", r"$\Delta u_{qp}$ [J m$^{-3}$]", True),
                ("u_e_J_m3", r"$\Delta u_e$ [J m$^{-3}$]", True),
                ("u_ph_J_m3", r"$\Delta u_{ph}$ [J m$^{-3}$]", True),
            ),
        ),
        (
            "Normalized energetic profiles (one scale per panel)",
            (
                ("P_delta_W_m3", r"$P_\Delta/\max|P_\Delta|$", False),
                ("P_q_W_m3", r"$P_q/\max|P_q|$", False),
                ("P_spec_W_m3", r"$P_{spec}/\max|P_{spec}|$", False),
                ("residual_W_m3", r"$R_e/\max|R_e|$", False),
            ),
        ),
    )
    with PdfPages(output) as pdf:
        for page_index, (title, fields) in enumerate(groups):
            fig = _profile_page(
                dataset,
                title=title,
                fields=fields,
                subtract_initial=page_index == 1,
                normalize=page_index == 2,
            )
            pdf.savefig(fig, dpi=dpi)
            plt.close(fig)
    return output


def _map_page(
    dataset: Mapping[str, Any],
    *,
    title: str,
    fields: Sequence[tuple[str, str, str, bool]],
) -> plt.Figure:
    x = _array(dataset, "nodes_x_nm")
    y = _array(dataset, "nodes_y_nm")
    triangles = np.asarray(dataset["triangles"], dtype=np.int64)
    tri = mtri.Triangulation(x, y, triangles)
    times = _array(dataset, "selected_times_ps")
    n_rows = times.size
    fig, axes = plt.subplots(
        n_rows,
        len(fields),
        figsize=(THESIS_WIDTH_IN, 2.0 + 1.65 * n_rows),
        squeeze=False,
    )
    fig.subplots_adjust(left=0.055, right=0.985, bottom=0.075, top=0.875, wspace=0.18, hspace=0.13)
    central = np.asarray(dataset["central_mask"], dtype=bool)
    window_nm = _scalar(dataset, "window_nm")
    center_x = 0.5 * (float(np.nanmin(x)) + float(np.nanmax(x)))
    for col, (key, label, cmap, diverging) in enumerate(fields):
        values = np.asarray(dataset[f"selected_{key}"], dtype=float)
        norm = _map_norm(values[:, central], diverging=diverging, binary=(key == "catalog_clipped"))
        mappable = None
        for row in range(n_rows):
            axis = axes[row, col]
            mappable = axis.tripcolor(tri, values[row], shading="gouraud", cmap=cmap, norm=norm, rasterized=True)
            axis.axvline(center_x - 0.5 * window_nm, color="white", linewidth=0.45, alpha=0.75)
            axis.axvline(center_x + 0.5 * window_nm, color="white", linewidth=0.45, alpha=0.75)
            axis.set_aspect("equal", adjustable="box")
            axis.set_xlim(float(np.nanmin(x)), float(np.nanmax(x)))
            axis.set_ylim(float(np.nanmin(y)), float(np.nanmax(y)))
            axis.grid(False)
            if row < n_rows - 1:
                axis.tick_params(labelbottom=False)
            if col:
                axis.tick_params(labelleft=False)
            if col == len(fields) - 1:
                axis.text(1.025, 0.5, rf"$t={times[row]:.3g}$ ps", transform=axis.transAxes, rotation=-90, va="center", fontsize=7.0)
        position = axes[0, col].get_position()
        color_axis = fig.add_axes([position.x0, position.y1 + 0.018, position.width, 0.014])
        colorbar = fig.colorbar(mappable, cax=color_axis, orientation="horizontal")
        if isinstance(norm, SymLogNorm):
            limit = max(abs(float(norm.vmin)), abs(float(norm.vmax)))
            colorbar.set_ticks([-limit, -float(norm.linthresh), 0.0, float(norm.linthresh), limit])
        else:
            colorbar.locator = MaxNLocator(nbins=4)
            colorbar.update_ticks()
        color_axis.xaxis.set_ticks_position("top")
        color_axis.xaxis.set_label_position("top")
        colorbar.set_label(label, fontsize=7.5, labelpad=1.5)
        color_axis.tick_params(labelsize=6.0, pad=0.5)
    fig.supxlabel(r"$x$ [nm]", y=0.018)
    fig.supylabel(r"$y$ [nm]", x=0.012)
    fig.suptitle(f"{title}\n{_run_scope_title(dataset)}", y=0.985, fontsize=9.5)
    return fig


def _profile_page(
    dataset: Mapping[str, Any],
    *,
    title: str,
    fields: Sequence[tuple[str, str, bool]],
    subtract_initial: bool,
    normalize: bool,
) -> plt.Figure:
    x = _array(dataset, "nodes_x_nm")
    weights = _array(dataset, "node_area_m2")
    mask = np.asarray(dataset["central_mask"], dtype=bool)
    times = _array(dataset, "selected_times_ps")
    colors = ("tab:blue", "tab:orange", "tab:green", "tab:red")
    fig, axes = plt.subplots(2, 2, figsize=(THESIS_WIDTH_IN, 5.25), sharex=True)
    fig.subplots_adjust(left=0.10, right=0.975, bottom=0.095, top=0.89, wspace=0.28, hspace=0.22)
    for axis, (key, label, _) in zip(axes.flat, fields):
        values = np.asarray(dataset[f"selected_{key}"], dtype=float).copy()
        if subtract_initial:
            values -= np.asarray(dataset[f"initial_{key}"], dtype=float)[None, :]
        scale = 1.0
        if normalize:
            finite = np.abs(values[:, mask])
            scale = max(float(np.nanmax(finite)), 1.0e-300)
        for row, time_ps in enumerate(times):
            xp, yp = _binned_profile(x[mask], values[row, mask] / scale, weights[mask], n_bins=72)
            axis.plot(xp, yp, color=colors[row % len(colors)], linewidth=1.0, label=rf"$t={time_ps:.3g}$ ps")
        axis.set_ylabel(label)
        axis.grid(True)
        axis.legend(frameon=False, fontsize=6.8)
    for axis in axes[1]:
        axis.set_xlabel(r"$x$ [nm]")
    fig.suptitle(f"{title}\nCross-width area-weighted mean; {_run_scope_title(dataset)}", y=0.985, fontsize=9.3)
    return fig


def _binned_profile(
    x: np.ndarray,
    values: np.ndarray,
    weights: np.ndarray,
    *,
    n_bins: int,
) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray(x, dtype=float)
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    edges = np.linspace(float(np.nanmin(x)), float(np.nanmax(x)), max(2, int(n_bins)) + 1)
    bins = np.clip(np.digitize(x, edges) - 1, 0, edges.size - 2)
    valid = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    weight_sum = np.bincount(bins[valid], weights=weights[valid], minlength=edges.size - 1)
    value_sum = np.bincount(bins[valid], weights=weights[valid] * values[valid], minlength=edges.size - 1)
    profile = np.divide(value_sum, weight_sum, out=np.full_like(value_sum, np.nan), where=weight_sum > 0.0)
    centers = 0.5 * (edges[:-1] + edges[1:])
    keep = np.isfinite(profile)
    return centers[keep], profile[keep]


def _map_norm(values: np.ndarray, *, diverging: bool, binary: bool) -> Normalize:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if binary:
        return Normalize(vmin=0.0, vmax=1.0)
    if finite.size == 0:
        return Normalize(vmin=-1.0 if diverging else 0.0, vmax=1.0)
    if diverging:
        limit = max(float(np.nanpercentile(np.abs(finite), 99.0)), 1.0e-30)
        return SymLogNorm(linthresh=max(1.0e-3 * limit, 1.0e-30), vmin=-limit, vmax=limit, base=10)
    low = float(np.nanpercentile(finite, 1.0))
    high = float(np.nanpercentile(finite, 99.0))
    if high <= low:
        high = low + max(abs(low), 1.0)
    return Normalize(vmin=low, vmax=high)


def _event_lines(axis: plt.Axes, dataset: Mapping[str, Any]) -> None:
    photon = _scalar(dataset, "photon_time_ps", default=np.nan)
    vmax = _scalar(dataset, "vmax_time_ps", default=np.nan)
    if np.isfinite(photon):
        axis.axvline(photon, color="0.25", linestyle="--", linewidth=0.7)
    if np.isfinite(vmax):
        axis.axvline(vmax, color="tab:blue", linestyle=":", linewidth=0.7)


def _run_scope_title(dataset: Mapping[str, Any]) -> str:
    window = _scalar(dataset, "window_nm")
    processed = int(_scalar(dataset, "processed_snapshot_count"))
    total = int(_scalar(dataset, "total_snapshot_count"))
    stride = int(_scalar(dataset, "snapshot_stride"))
    prefix = "SMOKE/TRUNCATED; " if bool(np.asarray(dataset["truncated"]).reshape(-1)[0]) else ""
    return f"{prefix}central window={window:g} nm; snapshots={processed}/{total}; stride={stride}"


def _array(dataset: Mapping[str, Any], key: str) -> np.ndarray:
    return np.asarray(dataset[key], dtype=float).reshape(-1)


def _scalar(dataset: Mapping[str, Any], key: str, *, default: float = np.nan) -> float:
    values = np.asarray(dataset.get(key, []), dtype=float).reshape(-1)
    return float(values[0]) if values.size else float(default)


def _pdf_path(path: str | Path) -> Path:
    output = Path(path)
    if output.suffix.lower() != ".pdf":
        output = output.with_suffix(".pdf")
    output.parent.mkdir(parents=True, exist_ok=True)
    return output


__all__ = [
    "plot_energy_projection_colormaps",
    "plot_energy_projection_profiles",
    "plot_energy_projection_temporal",
    "write_energy_projection_figures",
]
