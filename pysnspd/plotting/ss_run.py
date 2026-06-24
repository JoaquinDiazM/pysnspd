"""Diagnostic plots for OE7 stationary gTDGL/Poisson runs.

The snapshot plotting keeps pySNSPD's one-figure-per-field design, but uses the
notebook diagnostic conventions when the needed scales are present in history:

* |Delta| is plotted as |Delta|/Delta0;
* current-density snapshots are plotted as |j|/j_avg, |j_s|/j_avg, |j_n|/j_avg;
* potential is plotted in mV;
* arrows are sampled on a regular spatial grid, not by mesh-node ordering.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg", force=True)

import matplotlib.pyplot as plt
import matplotlib.tri as mtri

from pysnspd.gtdgl.operators import (
    strip_transport_current_profile_from_node_vectors,
    unwrap_phase_graph,
)

MEV_J = 1.602176634e-22


# =============================================================================
# Snapshot diagnostics
# =============================================================================


def plot_ss_available_snapshots(
    mesh,
    history: dict,
    output_dir: str | Path,
    *,
    dpi: int = 480,
    ncols: int = 3,
) -> dict[str, Path]:
    """Plot all snapshot fields currently present in the relaxation history."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    saved: dict[str, Path] = {}
    maybe_plots = [
        ("phi", plot_ss_phi_snapshots, output_dir / "ss_phi_snapshots.png"),
        ("delta", plot_ss_delta_snapshots, output_dir / "ss_delta_snapshots.png"),
        ("phase", plot_ss_phase_snapshots, output_dir / "ss_phase_snapshots.png"),
        ("current_density", plot_ss_current_density_snapshots, output_dir / "ss_current_density_snapshots.png"),
        ("supercurrent_density", plot_ss_supercurrent_density_snapshots, output_dir / "ss_supercurrent_density_snapshots.png"),
        ("normal_current_density", plot_ss_normal_current_density_snapshots, output_dir / "ss_normal_current_density_snapshots.png"),
        ("divergence", plot_ss_divergence_snapshots, output_dir / "ss_divergence_snapshots.png"),
        ("pairbreaking", plot_ss_pairbreaking_snapshots, output_dir / "ss_pairbreaking_snapshots.png"),
    ]

    for name, func, path in maybe_plots:
        try:
            saved[name] = func(mesh, history, path, dpi=dpi, ncols=ncols)
        except (KeyError, ValueError):
            continue

    # Extra edge-level zoom diagnostics. These do not replace the mesoscopic
    # snapshots; they live in a dedicated subfolder and inspect the literal
    # edge projections j_s,e and j_n,e between neighboring nodes.
    try:
        edge_saved = plot_ss_edge_current_zoom_snapshots(
            mesh,
            history,
            output_dir / "edge_current_diagnostics",
            dpi=dpi,
        )
        saved.update({f"edge_zoom_{key}": path for key, path in edge_saved.items()})
    except (KeyError, ValueError):
        pass

    return saved


def plot_ss_phi_snapshots(
    mesh,
    history: dict,
    output_path: str | Path,
    *,
    dpi: int = 480,
    ncols: int = 3,
) -> Path:
    """Plot electrostatic-potential snapshots in mV."""
    phi = 1.0e3 * _require_history_array(history, ("phi_snapshot_V",))
    t_s = _snapshot_times(history, ("phi_snapshot_t_s",), phi.shape[0])
    return _plot_snapshot_grid(
        mesh,
        phi,
        t_s,
        output_path,
        title="OE7 SS: electrostatic-potential snapshots",
        label="φ [mV]",
        symmetric=True,
        dpi=dpi,
        ncols=ncols,
    )


def plot_ss_delta_snapshots(
    mesh,
    history: dict,
    output_path: str | Path,
    *,
    dpi: int = 480,
    ncols: int = 3,
) -> Path:
    """Plot |Delta|/Delta0 snapshots."""
    if "delta_snapshot_meV" in history:
        delta_meV = np.asarray(history["delta_snapshot_meV"], dtype=float)
        t_s = _snapshot_times(history, ("delta_snapshot_t_s", "phi_snapshot_t_s"), delta_meV.shape[0])
    else:
        psi = _psi_snapshots_from_history(history)
        delta_meV = np.abs(psi) / MEV_J
        t_s = _snapshot_times(history, ("psi_snapshot_t_s", "delta_snapshot_t_s", "phi_snapshot_t_s"), delta_meV.shape[0])
    delta0_meV = _history_scalar(history, "delta0_meV")
    if delta0_meV is not None and delta0_meV > 0.0:
        z = delta_meV / delta0_meV
        label = r"$|\Delta|/\Delta_0$"
        title = rf"OE7 SS: $|\Delta|/\Delta_0$ snapshots ($\Delta_0={delta0_meV:.4g}$ meV)"
        vmax = max(1.05, float(np.nanmax(z)))
    else:
        z = delta_meV
        label = "|Δ| [meV]"
        title = "OE7 SS: |Δ| snapshots"
        vmax = None
    return _plot_snapshot_grid(
        mesh,
        z,
        t_s,
        output_path,
        title=title,
        label=label,
        vmin=0.0,
        vmax=vmax,
        dpi=dpi,
        ncols=ncols,
    )


def plot_ss_phase_snapshots(
    mesh,
    history: dict,
    output_path: str | Path,
    *,
    dpi: int = 480,
    ncols: int = 3,
) -> Path:
    """Plot phase residual snapshots after removing best linear ramp."""
    psi = _psi_snapshots_from_history(history)
    t_s = _snapshot_times(history, ("psi_snapshot_t_s", "phase_snapshot_t_s", "phi_snapshot_t_s"), psi.shape[0])
    edges = _mesh_edges_from_triangles(mesh)
    theta = []
    for k in range(psi.shape[0]):
        th = _unwrap_phase_safe(psi[k], edges, seed_index=_center_node_index(mesh))
        th = _remove_best_linear_ramp(mesh, th)
        th -= float(np.nanmean(th))
        theta.append(th)
    theta = np.vstack(theta)
    return _plot_snapshot_grid(
        mesh,
        theta,
        t_s,
        output_path,
        title="OE7 SS: phase residual snapshots",
        label="phase residual [rad]",
        symmetric=True,
        robust_percentile=99.5,
        min_vmax=1.0e-12,
        dpi=dpi,
        ncols=ncols,
    )


def plot_ss_current_density_snapshots(
    mesh,
    history: dict,
    output_path: str | Path,
    *,
    dpi: int = 480,
    ncols: int = 3,
) -> Path:
    """Plot |j_tot|/j_avg snapshots with vector arrows."""
    return _plot_current_family_snapshots(
        mesh,
        history,
        output_path,
        mag_names=("jtot_snapshot_mag_A_m2", "jmag_snapshot_A_m2", "current_density_snapshot_A_m2"),
        x_names=("jtot_snapshot_x_A_m2", "current_density_snapshot_x_A_m2", "node_jtot_x_snapshot_A_m2", "jx_snapshot_A_m2"),
        y_names=("jtot_snapshot_y_A_m2", "current_density_snapshot_y_A_m2", "node_jtot_y_snapshot_A_m2", "jy_snapshot_A_m2"),
        t_names=("jtot_snapshot_t_s", "current_snapshot_t_s", "phi_snapshot_t_s"),
        title="OE7 SS: total current-density snapshots",
        abs_label=r"|j| [A m$^{-2}$]",
        norm_label=r"$|j|/j_{\rm avg}$",
        dpi=dpi,
        ncols=ncols,
    )


def plot_ss_supercurrent_density_snapshots(
    mesh,
    history: dict,
    output_path: str | Path,
    *,
    dpi: int = 480,
    ncols: int = 3,
) -> Path:
    """Plot |j_s|/j_avg snapshots with vector arrows."""
    return _plot_current_family_snapshots(
        mesh,
        history,
        output_path,
        mag_names=("js_us_snapshot_mag_A_m2", "supercurrent_density_snapshot_A_m2"),
        x_names=("js_us_snapshot_x_A_m2", "supercurrent_density_snapshot_x_A_m2"),
        y_names=("js_us_snapshot_y_A_m2", "supercurrent_density_snapshot_y_A_m2"),
        t_names=("js_us_snapshot_t_s", "supercurrent_snapshot_t_s", "current_snapshot_t_s", "phi_snapshot_t_s"),
        title="OE7 SS: supercurrent-density snapshots",
        abs_label=r"|j_s| [A m$^{-2}$]",
        norm_label=r"$|j_s|/j_{\rm avg}$",
        dpi=dpi,
        ncols=ncols,
    )


def plot_ss_normal_current_density_snapshots(
    mesh,
    history: dict,
    output_path: str | Path,
    *,
    dpi: int = 480,
    ncols: int = 3,
) -> Path:
    """Plot |j_n|/j_avg snapshots with vector arrows."""
    return _plot_current_family_snapshots(
        mesh,
        history,
        output_path,
        mag_names=("jn_snapshot_mag_A_m2", "normal_current_density_snapshot_A_m2"),
        x_names=("jn_snapshot_x_A_m2", "normal_current_density_snapshot_x_A_m2"),
        y_names=("jn_snapshot_y_A_m2", "normal_current_density_snapshot_y_A_m2"),
        t_names=("jn_snapshot_t_s", "normal_current_snapshot_t_s", "current_snapshot_t_s", "phi_snapshot_t_s"),
        title="OE7 SS: normal-current-density snapshots",
        abs_label=r"|j_n| [A m$^{-2}$]",
        norm_label=r"$|j_n|/j_{\rm avg}$",
        dpi=dpi,
        ncols=ncols,
    )


def plot_ss_edge_current_zoom_snapshots(
    mesh,
    history: dict,
    output_dir: str | Path,
    *,
    dpi: int = 480,
) -> dict[str, Path]:
    """Plot node-to-node edge currents in three diagnostic zoom windows.

    One PNG is emitted per stored snapshot inside ``output_dir``. Each figure has
    three columns: a top insulating-edge window, a left terminal-edge window,
    and a bottom-right corner window. The upper row shows the signed
    supercurrent edge projection ``j_s,e/j_avg``; the lower row shows the signed
    normal-current edge projection ``j_n,e/j_avg``.

    These plots intentionally use the exact solver edge ordering saved in the
    relaxation history; they do not rebuild edges from triangles.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    nodes_m = np.asarray(mesh.nodes, dtype=float)
    x_nm = nodes_m[:, 0] * 1.0e9
    y_nm = nodes_m[:, 1] * 1.0e9

    edge_i = _require_history_array(history, ("edge_i",)).astype(np.int64)
    edge_j = _require_history_array(history, ("edge_j",)).astype(np.int64)
    js = _require_history_array(history, ("edge_js_us_snapshot_A_m2",))
    jn = _require_history_array(history, ("edge_jn_snapshot_A_m2",))
    t_s = _snapshot_times(history, ("edge_snapshot_t_s", "snapshot_t_s", "phi_snapshot_t_s"), js.shape[0])

    if js.shape != jn.shape:
        raise ValueError(f"edge current snapshots must have matching shapes, got {js.shape} and {jn.shape}.")
    if js.ndim != 2:
        raise ValueError(f"edge current snapshots must be 2D, got shape {js.shape}.")
    if js.shape[1] != edge_i.size:
        raise ValueError(
            f"edge current snapshot has {js.shape[1]} edges, but edge_i has {edge_i.size}."
        )

    javg = _history_scalar(history, "javg_A_m2")
    scale = abs(javg) if javg is not None and abs(javg) > 0.0 else 1.0
    scale_label = r"/ $j_{\rm avg}$" if javg is not None and abs(javg) > 0.0 else r" [A m$^{-2}$]"

    regions = _edge_zoom_regions(x_nm, y_nm, edge_i, edge_j)
    paths: dict[str, Path] = {}

    for k in range(js.shape[0]):
        out = output_dir / f"edge_currents_snapshot_{k:03d}.png"
        _plot_one_edge_current_zoom_snapshot(
            x_nm=x_nm,
            y_nm=y_nm,
            edge_i=edge_i,
            edge_j=edge_j,
            js=js[k] / scale,
            jn=jn[k] / scale,
            t_s=float(t_s[k]),
            regions=regions,
            output_path=out,
            scale_label=scale_label,
            dpi=dpi,
        )
        paths[f"snapshot_{k:03d}"] = out

    return paths


def plot_ss_divergence_snapshots(
    mesh,
    history: dict,
    output_path: str | Path,
    *,
    dpi: int = 480,
    ncols: int = 3,
) -> Path:
    """Plot div(j) snapshots."""
    div = _require_history_array(
        history,
        ("div_jtot_snapshot_A_m3", "node_div_jtot_snapshot_A_m3", "divergence_snapshot_A_m3"),
    )
    t_s = _snapshot_times(history, ("divergence_snapshot_t_s", "div_jtot_snapshot_t_s", "phi_snapshot_t_s"), div.shape[0])
    scale = _history_scalar(history, "javg_A_m2")
    if scale is not None and scale > 0.0:
        h = float(getattr(mesh, "target_spacing_m", 1.0e-9))
        z = div / max(scale / max(h, 1.0e-300), 1.0e-300)
        label = r"div(j)/(j_avg/h)"
    else:
        z = div
        label = r"div(j) [A m$^{-3}$]"
    return _plot_snapshot_grid(
        mesh,
        z,
        t_s,
        output_path,
        title="OE7 SS: finite-volume div(j) snapshots",
        label=label,
        symmetric=True,
        robust_percentile=99.5,
        dpi=dpi,
        ncols=ncols,
    )


def plot_ss_pairbreaking_snapshots(
    mesh,
    history: dict,
    output_path: str | Path,
    *,
    dpi: int = 480,
    ncols: int = 3,
) -> Path:
    """Plot pair-breaking-ratio snapshots."""
    chi = _require_history_array(
        history,
        ("pairbreaking_snapshot", "pairbreaking_ratio_snapshot", "node_pairbreaking_snapshot"),
    )
    t_s = _snapshot_times(history, ("pairbreaking_snapshot_t_s", "pairbreaking_ratio_snapshot_t_s", "phi_snapshot_t_s"), chi.shape[0])
    return _plot_snapshot_grid(
        mesh,
        chi,
        t_s,
        output_path,
        title="OE7 SS: pair-breaking ratio snapshots",
        label=r"$\chi_{\rm pb}$",
        vmin=0.0,
        robust_percentile=99.5,
        min_vmax=1.0,
        dpi=dpi,
        ncols=ncols,
    )


# =============================================================================
# Backward-compatible final-field diagnostics
# =============================================================================


def plot_ss_state_delta(mesh, state, output_path: str | Path, *, dpi: int = 480) -> Path:
    delta_meV = np.abs(state.psi_J) / MEV_J
    vmax = max(float(np.nanmax(delta_meV)), 1.0e-30)
    return _plot_node_scalar(mesh, delta_meV, output_path, title="OE7 SS: relaxed Δ", label="Δ [meV]", vmin=0.0, vmax=vmax, dpi=dpi)


def plot_ss_state_phase(mesh, state, output_path: str | Path, *, dpi: int = 480) -> Path:
    theta = _unwrap_phase_safe(np.asarray(state.psi_J, dtype=np.complex128), _mesh_edges_from_triangles(mesh), seed_index=_center_node_index(mesh))
    theta = _remove_best_linear_ramp(mesh, theta)
    return _plot_node_scalar(mesh, theta, output_path, title="OE7 SS: phase residual", label="phase residual [rad]", dpi=dpi)


def plot_ss_state_phi(mesh, state, output_path: str | Path, *, dpi: int = 480) -> Path:
    return _plot_node_scalar(mesh, 1.0e3 * state.phi_V, output_path, title="OE7 SS: electrostatic potential φ", label="φ [mV]", dpi=dpi)


def plot_ss_state_divergence(mesh, state, output_path: str | Path, *, dpi: int = 480) -> Path:
    div = np.asarray(state.currents.node_div_jtot_A_m3, dtype=float)
    vmax = max(float(np.nanmax(np.abs(div))) if div.size else 1.0, 1.0e-30)
    return _plot_node_scalar(mesh, div, output_path, title="OE7 SS: finite-volume div(j)", label=r"div(j) [A m$^{-3}$]", vmin=-vmax, vmax=vmax, dpi=dpi)


def plot_ss_state_current_density(mesh, state, output_path: str | Path, *, ops=None, dpi: int = 480) -> Path:
    del ops
    jx = np.asarray(state.currents.node_jtot_x_A_m2, dtype=float)
    jy = np.asarray(state.currents.node_jtot_y_A_m2, dtype=float)
    jmag = np.sqrt(jx * jx + jy * jy)
    return _plot_vector_snapshot_grid(mesh, jmag[None, :], jx[None, :], jy[None, :], np.array([0.0]), output_path, title="OE7 SS: total current density", label=r"|j| [A m$^{-2}$]", vmin=0.0, dpi=dpi, ncols=1)


def plot_ss_pairbreaking_ratio(mesh, state, output_path: str | Path, *, dpi: int = 480) -> Path:
    chi = np.asarray(state.currents.node_pairbreaking_ratio, dtype=float)
    finite = chi[np.isfinite(chi)]
    vmax = float(np.nanpercentile(finite, 99.5)) if finite.size else 1.0
    vmax = max(vmax, 1.0)
    return _plot_node_scalar(mesh, chi, output_path, title="OE7 SS: pairbreaking diagnostic", label=r"$\chi_{\rm pb}$", vmin=0.0, vmax=vmax, dpi=dpi)


# =============================================================================
# Scalar diagnostics
# =============================================================================


def plot_ss_boundary_currents(summary: dict, output_path: str | Path, *, dpi: int = 480) -> Path:
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    currents = summary["boundary_currents_A"]
    labels = ["left", "right", "bottom", "top"]
    values = [currents[f"{label}_A"] for label in labels]
    fig, ax = plt.subplots(figsize=(6.0, 3.6))
    ax.bar(labels, values)
    ax.axhline(0.0, linewidth=0.8)
    ax.set_title("OE7 SS: integrated boundary currents")
    ax.set_ylabel("current [A]")
    ax.grid(False)
    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output


def plot_ss_transport_current_profile(
    *,
    mesh,
    ops,
    state,
    output_path: str | Path,
    target_current_A: float,
    thickness_m: float,
    dpi: int = 480,
    n_bins: int = 41,
) -> Path:
    del ops
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    x_m, current_A = strip_transport_current_profile_from_node_vectors(mesh=mesh, jx_A_m2=np.asarray(state.currents.node_jtot_x_A_m2, dtype=float), thickness_m=float(thickness_m), n_bins=int(n_bins))
    fig, ax = plt.subplots(figsize=(8.0, 3.2))
    ax.plot(x_m * 1.0e9, current_A, marker="o", label="node-LS profile")
    ax.axhline(float(target_current_A), linestyle="--", label=r"$I_{\rm target}$")
    ax.set_title("OE7 SS: transport-current profile")
    ax.set_xlabel("x [nm]")
    ax.set_ylabel("I(x) [A]")
    ax.grid(False)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output


def plot_ss_relaxation_history(history: dict, output_path: str | Path, *, dpi: int = 480) -> Path:
    """Plot compact relaxation diagnostics in two panels."""
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    t_ps = np.asarray(history.get("t_s", []), dtype=float) / 1.0e-12
    fig, axes = plt.subplots(2, 1, figsize=(8.8, 6.2), sharex=True, constrained_layout=False)
    fig.subplots_adjust(left=0.105, right=0.970, bottom=0.090, top=0.925, hspace=0.42)
    ax_log, ax_lin = axes
    if t_ps.size:
        _plot_log_history_curve(ax_log, t_ps, history, "eta_R", r"$\max\Delta |\Delta|^2/\Delta_0^2$")
        _plot_log_history_curve(ax_log, t_ps, history, "current_residual", r"rms div(j)")
        _plot_log_history_curve(ax_log, t_ps, history, "pairbreaking_max", r"$\max\chi_{\rm pb}$")
        _plot_normalized_history_curve(ax_lin, t_ps, history, "terminal_voltage_V", r"$|V_{\rm TDGL}|$", unit="V", absolute=True)
        delta0_meV = _history_scalar(history, "delta0_meV")
        delta_label = r"$\min|\Delta|/\Delta_0$" if delta0_meV is None else rf"$\min|\Delta|/\Delta_0$; $\Delta_0={delta0_meV:.4g}$ meV"
        _plot_direct_history_curve(ax_lin, t_ps, history, "delta_min_over_delta0", delta_label)
        _plot_normalized_by_history_scale(ax_lin, t_ps, history, value_key="normal_current_max_A_m2", scale_key="total_current_max_A_m2", label=r"$\max |j_n|$", scale_label=r"$\max |j|$", unit=r"A m$^{-2}$", absolute=True)
    ax_log.set_title("OE7 SS: stiff relaxation diagnostics")
    ax_log.set_ylabel("diagnostic value")
    ax_log.grid(False)
    ax_log.legend(frameon=False)
    ax_lin.set_title("OE7 SS: normalized physical monitors")
    ax_lin.set_xlabel("t [ps]")
    ax_lin.set_ylabel("normalized value")
    ax_lin.set_ylim(0.0, 1.05)
    ax_lin.grid(False)
    ax_lin.legend(frameon=False)
    fig.savefig(output, dpi=dpi, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    return output


# =============================================================================
# Internal helpers
# =============================================================================


def _edge_zoom_regions(
    x_nm: np.ndarray,
    y_nm: np.ndarray,
    edge_i: np.ndarray,
    edge_j: np.ndarray,
) -> list[dict[str, object]]:
    """Return the three standard local edge-current diagnostic regions."""
    xmin = float(np.nanmin(x_nm)); xmax = float(np.nanmax(x_nm))
    ymin = float(np.nanmin(y_nm)); ymax = float(np.nanmax(y_nm))
    lx = max(xmax - xmin, 1.0)
    wy = max(ymax - ymin, 1.0)

    specs = [
        {
            "name": "top edge",
            "slug": "top_edge",
            "center": (0.40 * (xmin + xmax), ymax),
            "half_width": 0.045 * lx,
            "half_height": 0.095 * wy,
            "kind": "top",
        },
        {
            "name": "left edge",
            "slug": "left_edge",
            "center": (xmin, 0.60 * ymax),
            "half_width": 0.045 * lx,
            "half_height": 0.095 * wy,
            "kind": "left",
        },
        {
            "name": "bottom-right corner",
            "slug": "bottom_right_corner",
            "center": (xmax, ymin),
            "half_width": 0.055 * lx,
            "half_height": 0.115 * wy,
            "kind": "bottom_right",
        },
    ]

    regions: list[dict[str, object]] = []
    for spec in specs:
        cx, cy = spec["center"]
        hw = float(spec["half_width"])
        hh = float(spec["half_height"])
        xlo = cx - hw; xhi = cx + hw
        ylo = cy - hh; yhi = cy + hh

        xi = x_nm[edge_i]; xj = x_nm[edge_j]
        yi = y_nm[edge_i]; yj = y_nm[edge_j]
        mx = 0.5 * (xi + xj); my = 0.5 * (yi + yj)
        in_box = (mx >= xlo) & (mx <= xhi) & (my >= ylo) & (my <= yhi)
        in_box |= (
            ((xi >= xlo) & (xi <= xhi) & (yi >= ylo) & (yi <= yhi))
            | ((xj >= xlo) & (xj <= xhi) & (yj >= ylo) & (yj <= yhi))
        )
        idx = np.flatnonzero(in_box)
        if idx.size == 0:
            # Fallback to closest edge midpoint, so the diagnostic never silently
            # disappears on very small test meshes.
            d2 = (mx - cx) ** 2 + (my - cy) ** 2
            idx = np.asarray([int(np.nanargmin(d2))], dtype=int)

        highlight = _edge_zoom_highlight_edge(
            kind=str(spec["kind"]),
            x_nm=x_nm,
            y_nm=y_nm,
            edge_i=edge_i,
            edge_j=edge_j,
            edge_idx=idx,
            center_x=float(cx),
            center_y=float(cy),
        )
        regions.append(
            {
                "name": spec["name"],
                "slug": spec["slug"],
                "edge_idx": idx,
                "highlight_edge": highlight,
                "xlim": (xlo, xhi),
                "ylim": (ylo, yhi),
            }
        )
    return regions


def _edge_zoom_highlight_edge(
    *,
    kind: str,
    x_nm: np.ndarray,
    y_nm: np.ndarray,
    edge_i: np.ndarray,
    edge_j: np.ndarray,
    edge_idx: np.ndarray,
    center_x: float,
    center_y: float,
) -> int | None:
    """Pick the representative boundary-normal edge in a zoom region."""
    if edge_idx.size == 0:
        return None
    xi = x_nm[edge_i[edge_idx]]; xj = x_nm[edge_j[edge_idx]]
    yi = y_nm[edge_i[edge_idx]]; yj = y_nm[edge_j[edge_idx]]
    dx = xj - xi; dy = yj - yi
    ell = np.maximum(np.sqrt(dx * dx + dy * dy), 1.0e-300)
    mx = 0.5 * (xi + xj); my = 0.5 * (yi + yj)

    xmin = float(np.nanmin(x_nm)); xmax = float(np.nanmax(x_nm))
    ymin = float(np.nanmin(y_nm)); ymax = float(np.nanmax(y_nm))
    tol_x = max(1.0e-6, 0.02 * max(xmax - xmin, 1.0))
    tol_y = max(1.0e-6, 0.02 * max(ymax - ymin, 1.0))

    if kind == "top":
        touches = (np.abs(yi - ymax) < tol_y) | (np.abs(yj - ymax) < tol_y)
        score = touches.astype(float) * 10.0 + np.abs(dy) / ell - 0.02 * np.abs(mx - center_x)
    elif kind == "left":
        touches = (np.abs(xi - xmin) < tol_x) | (np.abs(xj - xmin) < tol_x)
        score = touches.astype(float) * 10.0 + np.abs(dx) / ell - 0.02 * np.abs(my - center_y)
    else:
        touches_right = (np.abs(xi - xmax) < tol_x) | (np.abs(xj - xmax) < tol_x)
        touches_bottom = (np.abs(yi - ymin) < tol_y) | (np.abs(yj - ymin) < tol_y)
        score = (touches_right | touches_bottom).astype(float) * 10.0
        score += 0.5 * (np.abs(dx) + np.abs(dy)) / ell
        score -= 0.02 * ((mx - center_x) ** 2 + (my - center_y) ** 2) ** 0.5
    return int(edge_idx[int(np.nanargmax(score))])


def _plot_one_edge_current_zoom_snapshot(
    *,
    x_nm: np.ndarray,
    y_nm: np.ndarray,
    edge_i: np.ndarray,
    edge_j: np.ndarray,
    js: np.ndarray,
    jn: np.ndarray,
    t_s: float,
    regions: list[dict[str, object]],
    output_path: str | Path,
    scale_label: str,
    dpi: int,
) -> Path:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(
        2,
        3,
        figsize=(14.3, 6.4),
        squeeze=False,
        constrained_layout=False,
    )
    fig.subplots_adjust(left=0.055, right=0.925, bottom=0.090, top=0.875, wspace=0.25, hspace=0.30)
    fig.suptitle(f"OE7 SS: edge-resolved currents, t = {t_s / 1.0e-12:.4g} ps", y=0.970)

    js_lim = max(1.0, float(np.nanpercentile(np.abs(js[np.isfinite(js)]), 99.0)) if np.any(np.isfinite(js)) else 1.0)
    jn_lim = max(0.05, float(np.nanpercentile(np.abs(jn[np.isfinite(jn)]), 99.0)) if np.any(np.isfinite(jn)) else 0.05)

    first_js = None
    first_jn = None
    for col, region in enumerate(regions):
        edge_idx = np.asarray(region["edge_idx"], dtype=int)
        highlight = region["highlight_edge"]
        for row, (field, clim, title, label) in enumerate(
            [
                (js, js_lim, r"signed $j_{s,e}$" + scale_label, r"$j_{s,e}$" + scale_label),
                (jn, jn_lim, r"signed $j_{n,e}$" + scale_label, r"$j_{n,e}$" + scale_label),
            ]
        ):
            ax = axes[row, col]
            mappable = _draw_edge_zoom_panel(
                ax=ax,
                x_nm=x_nm,
                y_nm=y_nm,
                edge_i=edge_i,
                edge_j=edge_j,
                edge_idx=edge_idx,
                values=field,
                clim=clim,
                region=region,
                highlight_edge=highlight,
            )
            if row == 0:
                ax.set_title(str(region["name"]), pad=9)
            if col == 0:
                ax.set_ylabel("y [nm]\n" + title)
            else:
                ax.set_ylabel("y [nm]")
            ax.set_xlabel("x [nm]")
            if row == 0 and first_js is None:
                first_js = mappable
            if row == 1 and first_jn is None:
                first_jn = mappable

    if first_js is not None:
        cax = fig.add_axes([0.942, 0.545, 0.018, 0.300])
        cb = fig.colorbar(first_js, cax=cax)
        cb.set_label(r"$j_{s,e}$" + scale_label, labelpad=8)
    if first_jn is not None:
        cax = fig.add_axes([0.942, 0.145, 0.018, 0.300])
        cb = fig.colorbar(first_jn, cax=cax)
        cb.set_label(r"$j_{n,e}$" + scale_label, labelpad=8)

    fig.savefig(output_path, dpi=dpi, bbox_inches="tight", pad_inches=0.10)
    plt.close(fig)
    return output_path


def _draw_edge_zoom_panel(
    *,
    ax,
    x_nm: np.ndarray,
    y_nm: np.ndarray,
    edge_i: np.ndarray,
    edge_j: np.ndarray,
    edge_idx: np.ndarray,
    values: np.ndarray,
    clim: float,
    region: dict[str, object],
    highlight_edge: int | None,
):
    import matplotlib as mpl
    from matplotlib.collections import LineCollection

    if edge_idx.size == 0:
        ax.axis("off")
        return mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=-clim, vmax=clim), cmap="coolwarm")

    segs = np.stack(
        [
            np.column_stack([x_nm[edge_i[edge_idx]], y_nm[edge_i[edge_idx]]]),
            np.column_stack([x_nm[edge_j[edge_idx]], y_nm[edge_j[edge_idx]]]),
        ],
        axis=1,
    )
    vals = np.asarray(values[edge_idx], dtype=float)
    norm = mpl.colors.Normalize(vmin=-float(clim), vmax=float(clim))
    lc = LineCollection(segs, array=vals, cmap="coolwarm", norm=norm, linewidths=2.1, alpha=0.95, zorder=2)
    ax.add_collection(lc)

    # Thin gray mesh skeleton below the current-carrying colored segments.
    skel = LineCollection(segs, colors="0.78", linewidths=0.55, alpha=0.65, zorder=1)
    ax.add_collection(skel)

    xi = x_nm[edge_i[edge_idx]]; xj = x_nm[edge_j[edge_idx]]
    yi = y_nm[edge_i[edge_idx]]; yj = y_nm[edge_j[edge_idx]]
    dx = xj - xi; dy = yj - yi
    ell = np.maximum(np.sqrt(dx * dx + dy * dy), 1.0e-300)
    mx = 0.5 * (xi + xj); my = 0.5 * (yi + yj)
    signed = vals
    local_len = 0.20 * np.nanmedian(ell) if ell.size else 1.0
    local_len = max(float(local_len), 0.25)
    qx = local_len * np.sign(signed) * dx / ell
    qy = local_len * np.sign(signed) * dy / ell
    keep = np.isfinite(signed) & (np.abs(signed) > 1.0e-10 * max(float(clim), 1.0e-300))
    if np.any(keep):
        ax.quiver(
            mx[keep],
            my[keep],
            qx[keep],
            qy[keep],
            np.clip(signed[keep], -clim, clim),
            cmap="coolwarm",
            norm=norm,
            angles="xy",
            scale_units="xy",
            scale=1.0,
            width=0.0042,
            headwidth=3.8,
            headlength=4.8,
            headaxislength=4.1,
            pivot="mid",
            zorder=3,
        )

    node_idx = np.unique(np.concatenate([edge_i[edge_idx], edge_j[edge_idx]]))
    ax.scatter(x_nm[node_idx], y_nm[node_idx], s=13, c="0.18", edgecolors="white", linewidths=0.35, zorder=4)

    if highlight_edge is not None:
        hi = int(highlight_edge)
        ax.plot(
            [x_nm[edge_i[hi]], x_nm[edge_j[hi]]],
            [y_nm[edge_i[hi]], y_nm[edge_j[hi]]],
            color="black",
            linewidth=3.0,
            alpha=0.90,
            zorder=5,
        )
        val = float(values[hi])
        hx = 0.5 * (x_nm[edge_i[hi]] + x_nm[edge_j[hi]])
        hy = 0.5 * (y_nm[edge_i[hi]] + y_nm[edge_j[hi]])
        ax.text(
            hx,
            hy,
            f"  {val:+.3g}",
            fontsize=8,
            va="center",
            ha="left",
            bbox={"boxstyle": "round,pad=0.18", "facecolor": "white", "alpha": 0.78, "edgecolor": "0.70"},
            zorder=6,
        )

    xlim = region["xlim"]; ylim = region["ylim"]
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(False)
    ax.tick_params(labelsize=8)
    return lc


def _plot_current_family_snapshots(
    mesh,
    history: dict,
    output_path: str | Path,
    *,
    mag_names: tuple[str, ...],
    x_names: tuple[str, ...],
    y_names: tuple[str, ...],
    t_names: tuple[str, ...],
    title: str,
    abs_label: str,
    norm_label: str,
    dpi: int,
    ncols: int,
) -> Path:
    jmag = _optional_history_array(history, mag_names)
    jx = _optional_history_array(history, x_names)
    jy = _optional_history_array(history, y_names)
    if jmag is None:
        if jx is None or jy is None:
            raise KeyError("history does not contain current snapshot magnitudes or vector components.")
        jmag = np.sqrt(jx * jx + jy * jy)
    t_s = _snapshot_times(history, t_names, jmag.shape[0])
    scale = _history_scalar(history, "javg_A_m2")
    if scale is not None and abs(scale) > 0.0:
        z = jmag / abs(scale)
        ux = None if jx is None else jx / abs(scale)
        uy = None if jy is None else jy / abs(scale)
        label = norm_label
    else:
        z = jmag
        ux = jx
        uy = jy
        label = abs_label
    if ux is not None and uy is not None:
        return _plot_vector_snapshot_grid(mesh, z, ux, uy, t_s, output_path, title=title, label=label, vmin=0.0, dpi=dpi, ncols=ncols)
    return _plot_snapshot_grid(mesh, z, t_s, output_path, title=title, label=label, vmin=0.0, dpi=dpi, ncols=ncols)


def _plot_snapshot_grid(
    mesh,
    values,
    t_s,
    output_path: str | Path,
    *,
    title: str,
    label: str,
    symmetric: bool = False,
    vmin: float | None = None,
    vmax: float | None = None,
    robust_percentile: float | None = None,
    min_vmax: float | None = None,
    dpi: int = 480,
    ncols: int = 3,
) -> Path:
    arr = np.asarray(values, dtype=float)
    if arr.ndim == 1:
        arr = arr[None, :]
    if arr.ndim != 2:
        raise ValueError(f"Snapshot array must be 2D, got shape {arr.shape}.")
    nodes = np.asarray(mesh.nodes, dtype=float)
    if arr.shape[1] != nodes.shape[0]:
        raise ValueError(f"Snapshot array has {arr.shape[1]} nodes, but mesh has {nodes.shape[0]} nodes.")
    t_s = np.asarray(t_s, dtype=float).reshape(-1)
    if t_s.size != arr.shape[0]:
        raise ValueError(f"Snapshot time axis has {t_s.size} entries, but data has {arr.shape[0]} snapshots.")
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    vmin, vmax = _resolve_limits(arr, symmetric=symmetric, vmin=vmin, vmax=vmax, robust_percentile=robust_percentile, min_vmax=min_vmax)
    x_nm = nodes[:, 0] * 1.0e9
    y_nm = nodes[:, 1] * 1.0e9
    tri = mtri.Triangulation(x_nm, y_nm, np.asarray(mesh.triangles, dtype=np.int64))
    n_snap = int(arr.shape[0])
    ncols = max(1, int(ncols))
    nrows = int(np.ceil(n_snap / ncols))
    fig, axes = _snapshot_figure(nrows=nrows, ncols=ncols)
    last_im = None
    for k in range(nrows * ncols):
        ax = axes.flat[k]
        if k >= n_snap:
            ax.axis("off")
            continue
        last_im = ax.tripcolor(tri, arr[k], shading="gouraud", vmin=vmin, vmax=vmax)
        _format_snapshot_axis(ax, t_s[k])
    if last_im is not None:
        cax = fig.add_axes([0.910, 0.165, 0.020, 0.675])
        cbar = fig.colorbar(last_im, cax=cax)
        cbar.set_label(label, labelpad=10)
    fig.suptitle(title, y=0.975)
    fig.savefig(output, dpi=dpi, bbox_inches="tight", pad_inches=0.12)
    plt.close(fig)
    return output


def _plot_vector_snapshot_grid(
    mesh,
    values,
    vx,
    vy,
    t_s,
    output_path: str | Path,
    *,
    title: str,
    label: str,
    symmetric: bool = False,
    vmin: float | None = None,
    vmax: float | None = None,
    robust_percentile: float | None = None,
    min_vmax: float | None = None,
    dpi: int = 480,
    ncols: int = 3,
    arrow_nx: int = 18,
    arrow_ny: int = 7,
) -> Path:
    arr = np.asarray(values, dtype=float)
    ux_raw = np.asarray(vx, dtype=float)
    uy_raw = np.asarray(vy, dtype=float)
    if arr.ndim == 1:
        arr = arr[None, :]
    if ux_raw.ndim == 1:
        ux_raw = ux_raw[None, :]
    if uy_raw.ndim == 1:
        uy_raw = uy_raw[None, :]
    if arr.shape != ux_raw.shape or arr.shape != uy_raw.shape:
        raise ValueError(f"Vector snapshot arrays must have matching shapes, got {arr.shape}, {ux_raw.shape}, {uy_raw.shape}.")
    nodes = np.asarray(mesh.nodes, dtype=float)
    if arr.shape[1] != nodes.shape[0]:
        raise ValueError(f"Snapshot array has {arr.shape[1]} nodes, but mesh has {nodes.shape[0]} nodes.")
    t_s = np.asarray(t_s, dtype=float).reshape(-1)
    if t_s.size != arr.shape[0]:
        raise ValueError(f"Snapshot time axis has {t_s.size} entries, but data has {arr.shape[0]} snapshots.")
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    vmin, vmax = _resolve_limits(arr, symmetric=symmetric, vmin=vmin, vmax=vmax, robust_percentile=robust_percentile, min_vmax=min_vmax)
    x_nm = nodes[:, 0] * 1.0e9
    y_nm = nodes[:, 1] * 1.0e9
    tri = mtri.Triangulation(x_nm, y_nm, np.asarray(mesh.triangles, dtype=np.int64))
    n_snap = int(arr.shape[0])
    ncols = max(1, int(ncols))
    nrows = int(np.ceil(n_snap / ncols))
    fig, axes = _snapshot_figure(nrows=nrows, ncols=ncols)
    x_span_nm = max(float(np.nanmax(x_nm) - np.nanmin(x_nm)), 1.0)
    y_span_nm = max(float(np.nanmax(y_nm) - np.nanmin(y_nm)), 1.0)
    arrow_len_nm = 0.050 * min(x_span_nm, y_span_nm)
    last_im = None
    for k in range(nrows * ncols):
        ax = axes.flat[k]
        if k >= n_snap:
            ax.axis("off")
            continue
        last_im = ax.tripcolor(tri, arr[k], shading="gouraud", vmin=vmin, vmax=vmax)
        jx = ux_raw[k]
        jy = uy_raw[k]
        jmag = np.sqrt(jx * jx + jy * jy)
        idx = _regular_arrow_indices(x_nm, y_nm, jmag, nx=arrow_nx, ny=arrow_ny, min_relative_weight=1.0e-10)
        if idx.size:
            qx = arrow_len_nm * jx[idx] / np.maximum(jmag[idx], 1.0e-300)
            qy = arrow_len_nm * jy[idx] / np.maximum(jmag[idx], 1.0e-300)
            ax.quiver(x_nm[idx], y_nm[idx], qx, qy, angles="xy", scale_units="xy", scale=1.0, width=0.0020, headwidth=3.5, headlength=4.5, headaxislength=3.8, pivot="mid")
        _format_snapshot_axis(ax, t_s[k])
    if last_im is not None:
        cax = fig.add_axes([0.910, 0.165, 0.020, 0.675])
        cbar = fig.colorbar(last_im, cax=cax)
        cbar.set_label(label, labelpad=10)
    fig.suptitle(title, y=0.975)
    fig.savefig(output, dpi=dpi, bbox_inches="tight", pad_inches=0.12)
    plt.close(fig)
    return output


def _snapshot_figure(*, nrows: int, ncols: int):
    fig_width = 4.35 * ncols + 1.10
    fig_height = 2.75 * nrows + 0.55
    fig, axes = plt.subplots(nrows, ncols, figsize=(fig_width, fig_height), squeeze=False, constrained_layout=False)
    fig.subplots_adjust(left=0.070, right=0.875, bottom=0.085, top=0.895, wspace=0.42, hspace=0.30)
    return fig, axes


def _format_snapshot_axis(ax, t_s: float) -> None:
    ax.set_title(f"t = {t_s / 1.0e-12:.4g} ps", pad=8)
    ax.set_xlabel("x [nm]", labelpad=4)
    ax.set_ylabel("y [nm]", labelpad=6)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(False)


def _regular_arrow_indices(x_nm: np.ndarray, y_nm: np.ndarray, weight: np.ndarray, *, nx: int = 18, ny: int = 7, min_relative_weight: float = 1.0e-10) -> np.ndarray:
    x_nm = np.asarray(x_nm, dtype=float)
    y_nm = np.asarray(y_nm, dtype=float)
    weight = np.asarray(weight, dtype=float)
    finite = np.isfinite(x_nm) & np.isfinite(y_nm) & np.isfinite(weight)
    if not np.any(finite):
        return np.array([], dtype=int)
    wmax = float(np.nanmax(np.abs(weight[finite])))
    if not np.isfinite(wmax) or wmax <= 0.0:
        return np.array([], dtype=int)
    finite &= np.abs(weight) > min_relative_weight * wmax
    if not np.any(finite):
        return np.array([], dtype=int)
    xmin = float(np.nanmin(x_nm[finite])); xmax = float(np.nanmax(x_nm[finite]))
    ymin = float(np.nanmin(y_nm[finite])); ymax = float(np.nanmax(y_nm[finite]))
    if xmax <= xmin or ymax <= ymin:
        return np.flatnonzero(finite)
    gx = np.linspace(xmin, xmax, int(nx))
    gy = np.linspace(ymin, ymax, int(ny))
    valid_idx = np.flatnonzero(finite)
    selected: list[int] = []
    for yy in gy:
        for xx in gx:
            dx = x_nm[valid_idx] - xx
            dy = y_nm[valid_idx] - yy
            d2 = (dx / max(xmax - xmin, 1.0e-30)) ** 2 + (dy / max(ymax - ymin, 1.0e-30)) ** 2
            selected.append(int(valid_idx[int(np.argmin(d2))]))
    return np.unique(np.asarray(selected, dtype=int))


def _resolve_limits(arr, *, symmetric, vmin, vmax, robust_percentile, min_vmax):
    finite = np.asarray(arr, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        finite = np.array([0.0], dtype=float)
    if vmin is not None and vmax is not None:
        return vmin, vmax
    if symmetric:
        scale = float(np.nanmax(np.abs(finite))) if robust_percentile is None else float(np.nanpercentile(np.abs(finite), robust_percentile))
        scale = max(scale, 1.0e-30)
        if min_vmax is not None:
            scale = max(scale, float(min_vmax))
        return (-scale if vmin is None else vmin, scale if vmax is None else vmax)
    local_vmin = float(np.nanmin(finite))
    local_vmax = float(np.nanmax(finite)) if robust_percentile is None else float(np.nanpercentile(finite, robust_percentile))
    if min_vmax is not None:
        local_vmax = max(local_vmax, float(min_vmax))
    if local_vmax <= local_vmin:
        pad = max(abs(local_vmax), 1.0) * 1.0e-12
        local_vmin -= pad
        local_vmax += pad
    return (local_vmin if vmin is None else vmin, local_vmax if vmax is None else vmax)


def _plot_node_scalar(mesh, values, output_path: str | Path, *, title: str, label: str, vmin=None, vmax=None, dpi: int = 480) -> Path:
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    nodes = np.asarray(mesh.nodes, dtype=float)
    x_nm = nodes[:, 0] * 1.0e9
    y_nm = nodes[:, 1] * 1.0e9
    tri = mtri.Triangulation(x_nm, y_nm, np.asarray(mesh.triangles, dtype=np.int64))
    z = np.asarray(values, dtype=float).reshape(-1)
    fig, ax = plt.subplots(figsize=(8.0, 3.2))
    im = ax.tripcolor(tri, z, shading="gouraud", vmin=vmin, vmax=vmax)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(label)
    ax.set_title(title)
    ax.set_xlabel("x [nm]")
    ax.set_ylabel("y [nm]")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(False)
    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output


def _plot_log_history_curve(ax, t_ps: np.ndarray, history: dict, key: str, label: str) -> None:
    if key not in history:
        return
    y = np.asarray(history[key], dtype=float).reshape(-1)
    t, y = _match_time_and_series(t_ps, y)
    if t.size == 0:
        return
    y_plot = np.maximum(np.abs(y), 1.0e-300)
    ax.semilogy(t, y_plot, label=f"{label}, final={y[-1]:.3e}")


def _plot_normalized_history_curve(ax, t_ps: np.ndarray, history: dict, key: str, label: str, *, unit: str, absolute: bool) -> None:
    if key not in history:
        return
    y = np.asarray(history[key], dtype=float).reshape(-1)
    t, y = _match_time_and_series(t_ps, y)
    if t.size == 0:
        return
    y_real = np.abs(y) if absolute else y
    scale = _normalization_scale(y_real)
    ax.plot(t, np.clip(y_real / scale, 0.0, 1.0), label=f"{label} / {_format_scale(scale, unit)}")


def _plot_direct_history_curve(ax, t_ps: np.ndarray, history: dict, key: str, label: str) -> None:
    if key not in history:
        return
    y = np.asarray(history[key], dtype=float).reshape(-1)
    t, y = _match_time_and_series(t_ps, y)
    if t.size == 0:
        return
    ax.plot(t, np.clip(y, 0.0, 1.0), label=f"{label}, final={y[-1]:.4g}")


def _plot_normalized_by_history_scale(ax, t_ps: np.ndarray, history: dict, *, value_key: str, scale_key: str, label: str, scale_label: str, unit: str, absolute: bool) -> None:
    if value_key not in history or scale_key not in history:
        return
    y = np.asarray(history[value_key], dtype=float).reshape(-1)
    scale_series = np.asarray(history[scale_key], dtype=float).reshape(-1)
    t, y = _match_time_and_series(t_ps, y)
    _, scale_series = _match_time_and_series(t_ps, scale_series)
    n = min(t.size, y.size, scale_series.size)
    if n <= 0:
        return
    y_real = np.abs(y[:n]) if absolute else y[:n]
    scale = _normalization_scale(np.abs(scale_series[:n]) if absolute else scale_series[:n])
    ax.plot(t[:n], np.clip(y_real / scale, 0.0, 1.0), label=f"{label} / {scale_label}={_format_scale(scale, unit)}")


def _match_time_and_series(t: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n = min(int(t.size), int(y.size))
    if n <= 0:
        return np.array([], dtype=float), np.array([], dtype=float)
    return t[:n], y[:n]


def _normalization_scale(y: np.ndarray) -> float:
    finite = np.asarray(y, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return 1.0
    scale = float(np.nanmax(np.abs(finite)))
    if not np.isfinite(scale) or scale <= 0.0:
        return 1.0
    return scale


def _format_scale(value: float, unit: str) -> str:
    value = float(value)
    if unit == "V":
        av = abs(value)
        if av >= 1.0:
            return f"{value:.3g} V"
        if av >= 1.0e-3:
            return f"{value * 1.0e3:.3g} mV"
        if av >= 1.0e-6:
            return f"{value * 1.0e6:.3g} µV"
        if av >= 1.0e-9:
            return f"{value * 1.0e9:.3g} nV"
        return f"{value:.3e} V"
    if unit:
        return f"{value:.3g} {unit}"
    return f"{value:.3g}"


def _optional_history_array(history: dict, names: tuple[str, ...]) -> np.ndarray | None:
    for name in names:
        if name in history:
            return np.asarray(history[name])
    return None


def _require_history_array(history: dict, names: tuple[str, ...]) -> np.ndarray:
    out = _optional_history_array(history, names)
    if out is None:
        raise KeyError(f"history does not contain any of these keys: {names}")
    return out


def _psi_snapshots_from_history(history: dict) -> np.ndarray:
    if "psi_snapshot_J" in history:
        psi = np.asarray(history["psi_snapshot_J"], dtype=np.complex128)
    else:
        real = _require_history_array(history, ("psi_snapshot_real_J", "psi_real_snapshot_J", "snapshot_psi_real_J"))
        imag = _require_history_array(history, ("psi_snapshot_imag_J", "psi_imag_snapshot_J", "snapshot_psi_imag_J"))
        psi = np.asarray(real, dtype=float) + 1j * np.asarray(imag, dtype=float)
    if psi.ndim != 2:
        raise ValueError(f"psi snapshots must be 2D, got shape {psi.shape}.")
    return psi


def _snapshot_times(history: dict, names: tuple[str, ...], n_snap: int) -> np.ndarray:
    for key in names:
        if key in history:
            t_s = np.asarray(history[key], dtype=float).reshape(-1)
            if t_s.size == n_snap:
                return t_s
    raise KeyError(f"No compatible snapshot time axis found. Tried {names}; expected {n_snap} entries.")


def _unwrap_phase_safe(psi: np.ndarray, edges: np.ndarray, *, seed_index: int) -> np.ndarray:
    try:
        return unwrap_phase_graph(np.asarray(psi, dtype=np.complex128), np.asarray(edges, dtype=np.int64), seed_index=int(seed_index), subtract_mean=False)
    except TypeError:
        return unwrap_phase_graph(np.asarray(psi, dtype=np.complex128), np.asarray(edges, dtype=np.int64))


def _center_node_index(mesh) -> int:
    nodes = np.asarray(mesh.nodes, dtype=float)
    center = np.array([0.5 * (float(np.min(nodes[:, 0])) + float(np.max(nodes[:, 0]))), 0.5 * (float(np.min(nodes[:, 1])) + float(np.max(nodes[:, 1])))], dtype=float)
    dist2 = np.sum((nodes[:, :2] - center[None, :]) ** 2, axis=1)
    return int(np.argmin(dist2))


def _mesh_edges_from_triangles(mesh) -> np.ndarray:
    tri = np.asarray(mesh.triangles, dtype=np.int64)
    edges = np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]])
    edges = np.sort(edges, axis=1)
    return np.unique(edges, axis=0)


def _remove_best_linear_ramp(mesh, phase: np.ndarray) -> np.ndarray:
    nodes = np.asarray(mesh.nodes, dtype=float)
    x = nodes[:, 0]
    y = nodes[:, 1]
    A = np.column_stack([x, y, np.ones_like(x)])
    good = np.isfinite(phase)
    if np.count_nonzero(good) < 3:
        return np.asarray(phase, dtype=float)
    coeff, *_ = np.linalg.lstsq(A[good], np.asarray(phase, dtype=float)[good], rcond=None)
    return np.asarray(phase, dtype=float) - A @ coeff


def _history_scalar(history: dict, key: str) -> float | None:
    if key not in history:
        return None
    arr = np.asarray(history[key], dtype=float).reshape(-1)
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return None
    return float(finite[-1])
