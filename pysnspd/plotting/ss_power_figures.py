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

MEV_J = 1.602176634e-22


def make_ss_snapshot_power_figures(
    *,
    mesh: Any,
    dataset: Mapping[str, Any],
    raw_ss: str | Path,
    output_dir: str | Path,
    dpi: int = 480,
) -> dict[str, Path]:
    """Create additional snapshot/power figures for a completed SS run.

    Missing diagnostic files are not an error because historical SS runs only
    contain final-state fields.  The returned dictionary is directly merged into
    the plotting-pipeline manifest.
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    raw = Path(raw_ss)

    snapshots = _load_npz_if_exists(raw / "stationary_snapshots.npz")
    power = _load_npz_if_exists(raw / "snapshot_power_energy_diagnostics.npz")

    saved: dict[str, Path] = {}
    if snapshots:
        saved["snapshot_state_atlas"] = plot_ss_snapshot_state_atlas(
            mesh, snapshots, dataset, out / "ss_snapshot_state_atlas.png", dpi=dpi
        )
        saved["snapshot_final_profile_comparison"] = plot_ss_snapshot_profile_comparison(
            mesh, snapshots, power, dataset, out / "ss_snapshot_profile_comparison.png", dpi=dpi
        )

    if power:
        saved["snapshot_power_energy_maps"] = plot_ss_snapshot_power_energy_maps(
            mesh, power, dataset, out / "ss_snapshot_power_energy_maps.png", dpi=dpi
        )
        saved["snapshot_power_balance_maps"] = plot_ss_snapshot_power_balance_maps(
            mesh, power, dataset, out / "ss_snapshot_power_balance_maps.png", dpi=dpi
        )
        saved["snapshot_runtime_metrics"] = plot_ss_snapshot_runtime_metrics(
            power, snapshots, dataset, out / "ss_snapshot_runtime_metrics.png", dpi=dpi
        )
        saved["snapshot_catalog_indices"] = plot_ss_snapshot_catalog_indices(
            mesh, power, dataset, out / "ss_snapshot_catalog_indices.png", dpi=dpi
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
    dpi: int = 480,
) -> Path:
    """Plot the dynamic SS state at all stored snapshots.

    Rows are |Delta|/Delta0, phi, and |j|/javg.  The intent is to see whether a
    60 ps SS relaxation really stopped changing, not just whether the final
    state looks reasonable.
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

    if delta_mev.size == 0:
        raise ValueError("stationary_snapshots.npz does not contain Delta snapshots")

    t_ps = _snapshot_times_ps(snapshots, preferred=("snapshot_t_s", "delta_snapshot_t_s"), n=delta_mev.shape[0])
    indices = _representative_snapshot_indices(delta_mev.shape[0], max_panels=9)
    t_sel = t_ps[indices]

    delta0 = _delta0_mev(dataset, snapshots)
    delta_norm = delta_mev / delta0 if delta0 > 0.0 else delta_mev
    phi_mv = 1.0e3 * _resize_snapshot_field(phi_v, delta_mev.shape)
    jscale = _javg(dataset, snapshots)
    j_norm = _resize_snapshot_field(jmag, delta_mev.shape) / jscale

    fields = [
        (delta_norm[indices], r"$|\Delta|/\Delta_0$", False, 0.0, None),
        (phi_mv[indices], r"$\phi$ [mV]", True, None, None),
        (j_norm[indices], r"$|j|/j_{avg}$", False, 0.0, None),
    ]
    row_titles = [r"order parameter", r"electrostatic potential", r"total current"]

    fig, axes = _snapshot_grid_figure(nrows=3, ncols=len(indices), title=f"SS snapshots: {dataset.get('run_name', '')}")
    for row, (z, label, symmetric, vmin, vmax) in enumerate(fields):
        norm, vmin_eff, vmax_eff = _node_color_limits(z, symmetric=symmetric, vmin=vmin, vmax=vmax)
        mappable = None
        for col, _idx in enumerate(indices):
            ax = axes[row, col]
            mappable = ax.tripcolor(tri, z[col], shading="gouraud", vmin=vmin_eff, vmax=vmax_eff, norm=norm)
            if row == 0:
                ax.set_title(f"t={t_sel[col]:.3g} ps")
            if col == 0:
                ax.set_ylabel(row_titles[row])
            _format_map_axis(ax, show_xlabel=(row == 2), show_ylabel=(col == 0))
        if mappable is not None:
            cb = fig.colorbar(mappable, ax=list(axes[row, :]), shrink=0.72, pad=0.012)
            cb.set_label(label)

    fig.savefig(output, dpi=dpi, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    return output


def plot_ss_snapshot_power_energy_maps(
    mesh: Any,
    power: Mapping[str, np.ndarray],
    dataset: Mapping[str, Any],
    output_path: str | Path,
    *,
    dpi: int = 480,
) -> Path:
    """Plot catalogue-derived power, Joule, and thermal transport maps."""
    output = _prepare_output(output_path)
    tri = _triangulation(mesh, dataset)

    p_ep = _snapshot_array(power, ("P_total_snapshot_W_m3",))
    joule = _snapshot_array(power, ("joule_snapshot_W_m3",), shape_like=p_ep)
    kappa = _snapshot_array(power, ("kappa_s_snapshot_W_m_K",), shape_like=p_ep)
    if p_ep.size == 0:
        raise ValueError("snapshot_power_energy_diagnostics.npz lacks P_total_snapshot_W_m3")

    t_ps = _snapshot_times_ps(power, preferred=("snapshot_t_s",), n=p_ep.shape[0])
    indices = _representative_snapshot_indices(p_ep.shape[0], max_panels=9)
    t_sel = t_ps[indices]

    fields = [
        (p_ep[indices], r"$P_{ep}=P_S+P_R$ [W m$^{-3}$]", "signed", r"electron $\rightarrow$ phonon power"),
        (joule[indices], r"$P_J$ [W m$^{-3}$]", "positive_log", r"Joule diagnostic"),
        (kappa[indices], r"$\kappa_s$ [W m$^{-1}$ K$^{-1}$]", "positive_log", r"thermal conductivity"),
    ]

    fig, axes = _snapshot_grid_figure(nrows=3, ncols=len(indices), title=f"SS runtime power/transport maps: {dataset.get('run_name', '')}")
    for row, (z, label, mode, row_title) in enumerate(fields):
        norm, vmin_eff, vmax_eff = _norm_for_mode(z, mode)
        mappable = None
        for col, _idx in enumerate(indices):
            ax = axes[row, col]
            mappable = ax.tripcolor(tri, z[col], shading="gouraud", vmin=vmin_eff, vmax=vmax_eff, norm=norm)
            if row == 0:
                ax.set_title(f"t={t_sel[col]:.3g} ps")
            if col == 0:
                ax.set_ylabel(row_title)
            _format_map_axis(ax, show_xlabel=(row == 2), show_ylabel=(col == 0))
        if mappable is not None:
            cb = fig.colorbar(mappable, ax=list(axes[row, :]), shrink=0.72, pad=0.012)
            cb.set_label(label)

    fig.savefig(output, dpi=dpi, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    return output


def plot_ss_snapshot_power_balance_maps(
    mesh: Any,
    power: Mapping[str, np.ndarray],
    dataset: Mapping[str, Any],
    output_path: str | Path,
    *,
    dpi: int = 480,
) -> Path:
    """Plot diagnostic local energy-balance tendencies at snapshots.

    ``P_J - P_ep`` is the local electronic heating tendency before temperature
    dynamics are actually coupled.  ``P_ep - P_esc`` is the analogous phonon
    tendency.  These are diagnostics only; they are not fed back into the SS
    solver yet.
    """
    output = _prepare_output(output_path)
    tri = _triangulation(mesh, dataset)

    p_ep = _snapshot_array(power, ("P_total_snapshot_W_m3",))
    joule = _snapshot_array(power, ("joule_snapshot_W_m3",), shape_like=p_ep)
    p_esc = _snapshot_array(power, ("P_esc_snapshot_W_m3",), shape_like=p_ep)
    if p_ep.size == 0:
        raise ValueError("snapshot_power_energy_diagnostics.npz lacks P_total_snapshot_W_m3")

    electron_balance = joule - p_ep
    phonon_balance = p_ep - p_esc

    t_ps = _snapshot_times_ps(power, preferred=("snapshot_t_s",), n=p_ep.shape[0])
    indices = _representative_snapshot_indices(p_ep.shape[0], max_panels=9)
    t_sel = t_ps[indices]

    fields = [
        (electron_balance[indices], r"$P_J-P_{ep}$ [W m$^{-3}$]", "electronic tendency"),
        (phonon_balance[indices], r"$P_{ep}-P_{esc}$ [W m$^{-3}$]", "phonon tendency"),
    ]

    fig, axes = _snapshot_grid_figure(nrows=2, ncols=len(indices), title=f"SS diagnostic power-balance maps: {dataset.get('run_name', '')}")
    for row, (z, label, row_title) in enumerate(fields):
        norm, vmin_eff, vmax_eff = _norm_for_mode(z, "signed")
        mappable = None
        for col, _idx in enumerate(indices):
            ax = axes[row, col]
            mappable = ax.tripcolor(tri, z[col], shading="gouraud", vmin=vmin_eff, vmax=vmax_eff, norm=norm)
            if row == 0:
                ax.set_title(f"t={t_sel[col]:.3g} ps")
            if col == 0:
                ax.set_ylabel(row_title)
            _format_map_axis(ax, show_xlabel=(row == 1), show_ylabel=(col == 0))
        if mappable is not None:
            cb = fig.colorbar(mappable, ax=list(axes[row, :]), shrink=0.72, pad=0.012)
            cb.set_label(label)

    fig.savefig(output, dpi=dpi, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    return output


def plot_ss_snapshot_runtime_metrics(
    power: Mapping[str, np.ndarray],
    snapshots: Mapping[str, np.ndarray] | None,
    dataset: Mapping[str, Any],
    output_path: str | Path,
    *,
    dpi: int = 480,
) -> Path:
    """Plot compact scalar metrics extracted from snapshot maps."""
    output = _prepare_output(output_path)

    p_ep = _snapshot_array(power, ("P_total_snapshot_W_m3",))
    joule = _snapshot_array(power, ("joule_snapshot_W_m3",), shape_like=p_ep)
    q_abs = _snapshot_array(power, ("q_abs_snapshot_m_inv",), shape_like=p_ep)
    u_e = _snapshot_array(power, ("u_e_snapshot_J_m3",), shape_like=p_ep)
    u_ph = _snapshot_array(power, ("u_ph_snapshot_J_m3",), shape_like=p_ep)
    c_e = _snapshot_array(power, ("C_e_snapshot_J_m3_K",), shape_like=p_ep)
    c_ph = _snapshot_array(power, ("C_ph_snapshot_J_m3_K",), shape_like=p_ep)

    if p_ep.size == 0:
        raise ValueError("snapshot_power_energy_diagnostics.npz lacks P_total_snapshot_W_m3")

    t_ps = _snapshot_times_ps(power, preferred=("snapshot_t_s",), n=p_ep.shape[0])
    delta_norm = np.empty((0, 0), dtype=float)
    if snapshots:
        delta_mev = _snapshot_array(snapshots, ("delta_snapshot_meV",), fallback=_delta_mev_from_psi_snapshots(snapshots))
        if delta_mev.size:
            delta0 = _delta0_mev(dataset, snapshots)
            delta_norm = delta_mev / delta0 if delta0 > 0.0 else delta_mev

    fig, axes = plt.subplots(2, 2, figsize=(10.0, 7.0), constrained_layout=False)
    fig.subplots_adjust(left=0.095, right=0.965, bottom=0.095, top=0.905, wspace=0.32, hspace=0.36)
    fig.suptitle(f"SS snapshot runtime diagnostics: {dataset.get('run_name', '')}", y=0.975)

    ax = axes[0, 0]
    _plot_snapshot_metric(ax, t_ps, p_ep, reducer="max_abs", label=r"max $|P_{ep}|$")
    _plot_snapshot_metric(ax, t_ps, joule, reducer="max", label=r"max $P_J$")
    _plot_snapshot_metric(ax, t_ps, joule - p_ep, reducer="max_abs", label=r"max $|P_J-P_{ep}|$")
    ax.set_yscale("symlog", linthresh=1.0e8)
    ax.set_ylabel(r"power density [W m$^{-3}$]")
    ax.set_title("power scales")
    ax.grid(False)
    ax.legend(frameon=False)

    ax = axes[0, 1]
    _plot_snapshot_metric(ax, t_ps, q_abs / 1.0e7, reducer="p99", label=r"p99 $|q|$")
    _plot_snapshot_metric(ax, t_ps, q_abs / 1.0e7, reducer="max", label=r"max $|q|$")
    if delta_norm.size:
        _plot_snapshot_metric(ax, t_ps, delta_norm, reducer="min", label=r"min $|\Delta|/\Delta_0$")
        _plot_snapshot_metric(ax, t_ps, delta_norm, reducer="mean", label=r"mean $|\Delta|/\Delta_0$")
    ax.set_ylabel(r"$q$ [$10^7$ m$^{-1}$] or normalized gap")
    ax.set_title("order parameter / momentum")
    ax.grid(False)
    ax.legend(frameon=False)

    ax = axes[1, 0]
    _plot_snapshot_metric(ax, t_ps, u_e, reducer="mean", label=r"mean $u_e$")
    _plot_snapshot_metric(ax, t_ps, u_ph, reducer="mean", label=r"mean $u_{ph}$")
    ax.set_ylabel(r"energy density [J m$^{-3}$]")
    ax.set_xlabel("t [ps]")
    ax.set_title("stored energy diagnostics")
    ax.grid(False)
    ax.legend(frameon=False)

    ax = axes[1, 1]
    _plot_snapshot_metric(ax, t_ps, c_e, reducer="mean", label=r"mean $C_e$")
    _plot_snapshot_metric(ax, t_ps, c_ph, reducer="mean", label=r"mean $C_{ph}$")
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
    dpi: int = 480,
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

    fig, axes = plt.subplots(1, len(fields), figsize=(3.1 * len(fields), 3.2), constrained_layout=False)
    axes = np.asarray(axes, dtype=object).reshape(1, -1)
    fig.subplots_adjust(left=0.055, right=0.965, bottom=0.145, top=0.845, wspace=0.32)
    fig.suptitle(f"PRE catalogue indices at final SS snapshot: {dataset.get('run_name', '')}", y=0.960)
    for ax, (z, label, key) in zip(axes.ravel(), fields):
        z = np.asarray(z, dtype=float)
        finite = z[np.isfinite(z)]
        vmin = float(np.nanmin(finite)) if finite.size else 0.0
        vmax = float(np.nanmax(finite)) if finite.size else 1.0
        if vmax <= vmin:
            vmax = vmin + 1.0
        mappable = ax.tripcolor(tri, z, shading="gouraud", vmin=vmin, vmax=vmax)
        ax.set_title(label)
        _format_map_axis(ax, show_xlabel=True, show_ylabel=(ax is axes.ravel()[0]))
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
    dpi: int = 480,
) -> Path:
    """Plot x-binned profiles at first/middle/final snapshot."""
    output = _prepare_output(output_path)
    x_nm = _mesh_x_nm(mesh, dataset)
    delta_mev = _snapshot_array(snapshots, ("delta_snapshot_meV",), fallback=_delta_mev_from_psi_snapshots(snapshots))
    if delta_mev.size == 0:
        raise ValueError("stationary_snapshots.npz does not contain Delta snapshots")
    delta0 = _delta0_mev(dataset, snapshots)
    delta_norm = delta_mev / delta0 if delta0 > 0.0 else delta_mev
    phi = 1.0e3 * _resize_snapshot_field(_snapshot_array(snapshots, ("phi_snapshot_V",)), delta_mev.shape)
    j_norm = _resize_snapshot_field(_snapshot_current_mag(snapshots, family="jtot"), delta_mev.shape) / _javg(dataset, snapshots)

    q_abs = np.empty((0, 0), dtype=float)
    p_balance = np.empty((0, 0), dtype=float)
    if power:
        p_ep = _snapshot_array(power, ("P_total_snapshot_W_m3",))
        joule = _snapshot_array(power, ("joule_snapshot_W_m3",), shape_like=p_ep)
        q_abs = _snapshot_array(power, ("q_abs_snapshot_m_inv",), shape_like=p_ep) / 1.0e7
        p_balance = joule - p_ep

    t_ps = _snapshot_times_ps(snapshots, preferred=("snapshot_t_s", "delta_snapshot_t_s"), n=delta_mev.shape[0])
    indices = _representative_snapshot_indices(delta_mev.shape[0], max_panels=3)

    fig, axes = plt.subplots(4, 1, figsize=(9.2, 8.4), sharex=True, constrained_layout=False)
    fig.subplots_adjust(left=0.105, right=0.970, bottom=0.080, top=0.930, hspace=0.30)
    fig.suptitle(f"SS x-profile evolution: {dataset.get('run_name', '')}", y=0.985)

    for idx in indices:
        label = f"t={t_ps[idx]:.3g} ps"
        _plot_binned_profile(axes[0], x_nm, delta_norm[idx], label=label)
        _plot_binned_profile(axes[1], x_nm, j_norm[idx], label=label)
        _plot_binned_profile(axes[2], x_nm, phi[idx], label=label)
        if q_abs.size:
            _plot_binned_profile(axes[3], x_nm, q_abs[idx], label=label)
        elif p_balance.size:
            _plot_binned_profile(axes[3], x_nm, p_balance[idx], label=label)

    axes[0].set_ylabel(r"$|\Delta|/\Delta_0$")
    axes[1].set_ylabel(r"$|j|/j_{avg}$")
    axes[2].set_ylabel(r"$\phi$ [mV]")
    axes[3].set_ylabel(r"$|q|$ [$10^7$ m$^{-1}$]" if q_abs.size else r"$P_J-P_{ep}$")
    axes[3].set_xlabel("x [nm]")
    for ax in axes:
        ax.grid(False)
        ax.legend(frameon=False, loc="best")
    fig.savefig(output, dpi=dpi, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    return output


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------


def _load_npz_if_exists(path: Path) -> dict[str, np.ndarray]:
    if not path.exists():
        return {}
    with np.load(path, allow_pickle=True) as data:
        return {key: np.asarray(data[key]) for key in data.files}


def _prepare_output(output_path: str | Path) -> Path:
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    return output


def _mesh_x_nm(mesh: Any, dataset: Mapping[str, Any]) -> np.ndarray:
    if "x_nm" in dataset:
        return np.asarray(dataset["x_nm"], dtype=float)
    nodes = np.asarray(getattr(mesh, "nodes", getattr(mesh, "sites", [])), dtype=float)
    return nodes[:, 0] * 1.0e9


def _mesh_y_nm(mesh: Any, dataset: Mapping[str, Any]) -> np.ndarray:
    if "y_nm" in dataset:
        return np.asarray(dataset["y_nm"], dtype=float)
    nodes = np.asarray(getattr(mesh, "nodes", getattr(mesh, "sites", [])), dtype=float)
    return nodes[:, 1] * 1.0e9


def _triangulation(mesh: Any, dataset: Mapping[str, Any]) -> mtri.Triangulation:
    x = _mesh_x_nm(mesh, dataset)
    y = _mesh_y_nm(mesh, dataset)
    triangles = np.asarray(dataset.get("triangles", getattr(mesh, "triangles", getattr(mesh, "elements", []))), dtype=np.int64)
    return mtri.Triangulation(x, y, triangles)


def _delta0_mev(dataset: Mapping[str, Any], snapshots: Mapping[str, np.ndarray] | None = None) -> float:
    for source in (dataset, snapshots or {}):
        if "delta0_meV" in source:
            arr = np.asarray(source["delta0_meV"], dtype=float).reshape(-1)
            if arr.size and np.isfinite(arr[-1]) and arr[-1] > 0.0:
                return float(arr[-1])
    return 1.0


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
    if shape_like is not None and shape_like.size:
        return np.zeros_like(np.asarray(shape_like, dtype=float))
    return np.empty((0, 0), dtype=float)


def _resize_snapshot_field(arr: np.ndarray, shape: tuple[int, int]) -> np.ndarray:
    a = np.asarray(arr, dtype=float)
    if a.size == 0:
        return np.zeros(shape, dtype=float)
    if a.shape == shape:
        return a
    return np.resize(a, shape)


def _delta_mev_from_psi_snapshots(snapshots: Mapping[str, np.ndarray] | None) -> np.ndarray:
    if not snapshots:
        return np.empty((0, 0), dtype=float)
    real = snapshots.get("psi_snapshot_real_J")
    imag = snapshots.get("psi_snapshot_imag_J")
    if real is None:
        return np.empty((0, 0), dtype=float)
    r = np.asarray(real, dtype=float)
    i = np.asarray(imag if imag is not None else np.zeros_like(r), dtype=float)
    return np.sqrt(r * r + i * i) / MEV_J


def _snapshot_current_mag(snapshots: Mapping[str, np.ndarray], *, family: str) -> np.ndarray:
    if family == "jtot":
        mag_keys = ("jtot_snapshot_mag_A_m2", "jmag_snapshot_A_m2", "current_density_snapshot_A_m2")
        x_keys = ("jtot_snapshot_x_A_m2", "current_density_snapshot_x_A_m2", "node_jtot_x_snapshot_A_m2", "jx_snapshot_A_m2")
        y_keys = ("jtot_snapshot_y_A_m2", "current_density_snapshot_y_A_m2", "node_jtot_y_snapshot_A_m2", "jy_snapshot_A_m2")
    elif family == "jn":
        mag_keys = ("jn_snapshot_mag_A_m2", "normal_current_density_snapshot_A_m2")
        x_keys = ("jn_snapshot_x_A_m2", "normal_current_density_snapshot_x_A_m2")
        y_keys = ("jn_snapshot_y_A_m2", "normal_current_density_snapshot_y_A_m2")
    else:
        raise ValueError(f"unknown current family: {family}")

    mag = _snapshot_array(snapshots, mag_keys)
    if mag.size:
        return mag
    x = _snapshot_array(snapshots, x_keys)
    y = _snapshot_array(snapshots, y_keys, shape_like=x)
    if x.size:
        return np.sqrt(x * x + y * y)
    return np.empty((0, 0), dtype=float)


def _snapshot_times_ps(data: Mapping[str, np.ndarray], *, preferred: tuple[str, ...], n: int) -> np.ndarray:
    for key in preferred:
        if key in data:
            arr = np.asarray(data[key], dtype=float).reshape(-1)
            if arr.size:
                if arr.size != n:
                    arr = np.resize(arr, n)
                return arr / 1.0e-12
    return np.arange(int(n), dtype=float)


def _representative_snapshot_indices(n: int, *, max_panels: int) -> np.ndarray:
    n = int(n)
    if n <= 0:
        return np.array([], dtype=int)
    if n <= max_panels:
        return np.arange(n, dtype=int)
    return np.unique(np.linspace(0, n - 1, int(max_panels)).round().astype(int))


def _snapshot_grid_figure(*, nrows: int, ncols: int, title: str):
    width = max(7.5, 2.12 * max(ncols, 1) + 1.2)
    height = max(3.0, 2.12 * max(nrows, 1) + 0.9)
    fig, axes = plt.subplots(nrows, ncols, figsize=(width, height), squeeze=False, constrained_layout=False)
    fig.subplots_adjust(left=0.055, right=0.910, bottom=0.080, top=0.900, wspace=0.10, hspace=0.24)
    fig.suptitle(title, y=0.975)
    return fig, axes


def _format_map_axis(ax, *, show_xlabel: bool, show_ylabel: bool) -> None:
    ax.set_aspect("equal", adjustable="box")
    ax.grid(False)
    if show_xlabel:
        ax.set_xlabel("x [nm]")
    else:
        ax.set_xticklabels([])
    if show_ylabel:
        ax.set_ylabel("y [nm]")
    else:
        ax.set_yticklabels([])


def _node_color_limits(
    values: np.ndarray,
    *,
    symmetric: bool,
    vmin: float | None,
    vmax: float | None,
):
    z = np.asarray(values, dtype=float)
    finite = z[np.isfinite(z)]
    if finite.size == 0:
        finite = np.array([0.0])
    norm = None
    if symmetric:
        vm = float(np.nanpercentile(np.abs(finite), 99.5))
        vm = max(vm, 1.0e-30)
        return norm, -vm, vm
    if vmax is None:
        vmax = float(np.nanpercentile(finite, 99.5))
    if vmin is None:
        vmin = float(np.nanpercentile(finite, 0.5))
    if not np.isfinite(vmax) or not np.isfinite(vmin) or vmax <= vmin:
        vmax = float(vmin) + 1.0
    return norm, vmin, vmax


def _norm_for_mode(values: np.ndarray, mode: str):
    z = np.asarray(values, dtype=float)
    finite = z[np.isfinite(z)]
    if finite.size == 0:
        finite = np.array([0.0])
    if mode == "signed":
        vmax = float(np.nanpercentile(np.abs(finite), 99.2))
        vmax = max(vmax, 1.0e-30)
        linthresh = max(vmax * 1.0e-6, 1.0e8)
        return SymLogNorm(linthresh=linthresh, vmin=-vmax, vmax=vmax), None, None
    if mode == "positive_log":
        pos = finite[finite > 0.0]
        if pos.size == 0:
            return None, 0.0, 1.0
        vmin = float(np.nanpercentile(pos, 1.0))
        vmax = float(np.nanpercentile(pos, 99.2))
        if not np.isfinite(vmin) or vmin <= 0.0:
            vmin = max(float(np.nanmin(pos)), 1.0e-300)
        if not np.isfinite(vmax) or vmax <= vmin:
            vmax = vmin * 10.0
        return LogNorm(vmin=vmin, vmax=vmax), None, None
    return _node_color_limits(z, symmetric=False, vmin=None, vmax=None)


def _plot_snapshot_metric(ax, t_ps: np.ndarray, values: np.ndarray, *, reducer: str, label: str) -> None:
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 2 or arr.size == 0:
        return
    if reducer == "max":
        y = np.nanmax(arr, axis=1)
    elif reducer == "max_abs":
        y = np.nanmax(np.abs(arr), axis=1)
    elif reducer == "min":
        y = np.nanmin(arr, axis=1)
    elif reducer == "mean":
        y = np.nanmean(arr, axis=1)
    elif reducer == "p99":
        y = np.nanpercentile(arr, 99.0, axis=1)
    else:
        raise ValueError(f"unknown reducer: {reducer}")
    n = min(np.asarray(t_ps).size, y.size)
    if n:
        ax.plot(np.asarray(t_ps)[:n], y[:n], marker="o", linewidth=1.2, label=label)


def _plot_binned_profile(ax, x_nm: np.ndarray, z: np.ndarray, *, label: str, n_bins: int = 80) -> None:
    x = np.asarray(x_nm, dtype=float)
    y = np.asarray(z, dtype=float)
    if x.size == 0 or y.size == 0:
        return
    if y.size != x.size:
        y = np.resize(y, x.size)
    bins = np.linspace(float(np.nanmin(x)), float(np.nanmax(x)), int(n_bins) + 1)
    centers = 0.5 * (bins[:-1] + bins[1:])
    which = np.clip(np.digitize(x, bins) - 1, 0, centers.size - 1)
    prof = np.full(centers.size, np.nan, dtype=float)
    for k in range(centers.size):
        mask = which == k
        if np.any(mask):
            prof[k] = float(np.nanmean(y[mask]))
    ax.plot(centers, prof, marker="o", markersize=2.2, linewidth=1.2, label=label)
