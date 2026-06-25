"""Diagnostic plots for OE7 stationary gTDGL/Poisson runs.

This module is intentionally plotting-only.  The functions accept the light
data containers used in the pySNSPD pipeline and save PNG diagnostics without
opening interactive windows.

OE7 note:
    The edge-current zoom diagnostic now uses four columns:
    top insulating edge, left terminal edge, bottom-right corner, and an
    interior strip.  The fourth column is meant to separate true interior
    current-distribution problems from boundary/corner artifacts.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable

import numpy as np

import matplotlib
matplotlib.use("Agg", force=True)

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.collections import LineCollection

from pysnspd.gtdgl.operators import unwrap_phase_graph
from pysnspd.gtdgl.operators import (
    boundary_currents_from_edge_scalar_least_squares,
    boundary_currents_from_node_vectors,
    edge_scalar_to_node_vector_least_squares,
    strip_transport_current_profile_from_node_vectors,
)

MEV_J = 1.602176634e-22


# ---------------------------------------------------------------------------
# Final-state scalar/vector diagnostics
# ---------------------------------------------------------------------------

def plot_ss_state_delta(mesh, state, output_path: str | Path, *, dpi: int = 480) -> Path:
    """Plot the relaxed order-parameter amplitude in meV.

    The color scale starts at zero to avoid visually amplifying tiny numerical
    variations around the nearly uniform superconducting gap.
    """
    delta_meV = np.abs(state.psi_J) / MEV_J
    vmax = max(_safe_nanmax(delta_meV), 1.0e-30)

    return _plot_node_scalar(
        mesh,
        delta_meV,
        output_path,
        title="OE7 SS: relaxed Δ",
        label="Δ [meV]",
        vmin=0.0,
        vmax=vmax,
        dpi=dpi,
    )


def plot_ss_state_phase(mesh, state, output_path: str | Path, *, dpi: int = 480) -> Path:
    """Plot an unwrapped phase diagnostic."""
    theta = _unwrap_phase_by_x(mesh, np.angle(state.psi_J))
    return _plot_node_scalar(
        mesh,
        theta,
        output_path,
        title="OE7 SS: unwrapped phase θ",
        label="θ [rad]",
        dpi=dpi,
    )


def plot_ss_state_phi(mesh, state, output_path: str | Path, *, dpi: int = 480) -> Path:
    """Plot electrostatic potential."""
    return _plot_node_scalar(
        mesh,
        state.phi_V,
        output_path,
        title="OE7 SS: electrostatic potential φ",
        label="φ [V]",
        dpi=dpi,
    )


def plot_ss_state_divergence(mesh, state, output_path: str | Path, *, dpi: int = 480) -> Path:
    """Plot finite-volume current divergence."""
    div = np.asarray(state.currents.node_div_jtot_A_m3, dtype=float)
    vmax = max(_safe_nanmax(np.abs(div)), 1.0e-30)

    return _plot_node_scalar(
        mesh,
        div,
        output_path,
        title="OE7 SS: finite-volume div(j)",
        label="div(j) [A m$^{-3}$]",
        vmin=-vmax,
        vmax=vmax,
        dpi=dpi,
    )


def plot_ss_state_current_density(
    mesh,
    state,
    output_path: str | Path,
    *,
    ops=None,
    dpi: int = 480,
) -> Path:
    """Plot total current-density magnitude and vectors.

    If FV operators are provided, reconstruct the node vector field from edge
    current projections using local least squares.  This matches the diagnostic
    philosophy of the older notebook better than the simple node average.
    """
    nodes = np.asarray(mesh.nodes, dtype=float)
    x_nm = nodes[:, 0] * 1.0e9
    y_nm = nodes[:, 1] * 1.0e9
    triangles = np.asarray(mesh.triangles, dtype=np.int64)

    if ops is not None:
        jx, jy = edge_scalar_to_node_vector_least_squares(
            state.currents.edge_jtot_A_m2,
            ops,
        )
    else:
        jx = np.asarray(state.currents.node_jtot_x_A_m2, dtype=float)
        jy = np.asarray(state.currents.node_jtot_y_A_m2, dtype=float)

    mag = np.sqrt(jx**2 + jy**2)
    vmax = max(_safe_nanmax(mag), 1.0e-30)

    output = _prepare_output(output_path)
    fig, ax = plt.subplots(figsize=(7.0, 3.2), constrained_layout=True)
    tri = mtri.Triangulation(x_nm, y_nm, triangles)

    im = ax.tripcolor(tri, mag, shading="gouraud", vmin=0.0, vmax=vmax)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(r"$|\vec{j}|$ [A m$^{-2}$]")

    step = max(1, mag.size // 120)
    scale = _safe_nanmax(mag)
    if np.isfinite(scale) and scale > 0.0:
        ax.quiver(
            x_nm[::step],
            y_nm[::step],
            jx[::step] / scale,
            jy[::step] / scale,
            angles="xy",
            scale_units="xy",
            scale=0.040,
            width=0.0025,
        )

    ax.set_title("OE7 SS: total current density")
    ax.set_xlabel("x [nm]")
    ax.set_ylabel("y [nm]")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(False)
    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    return output


def _plot_node_scalar(
    mesh,
    values,
    output_path: str | Path,
    *,
    title: str,
    label: str,
    vmin=None,
    vmax=None,
    dpi: int = 480,
) -> Path:
    """Plot a nodal scalar field on the Delaunay triangulation."""
    nodes = np.asarray(mesh.nodes, dtype=float)
    x_nm = nodes[:, 0] * 1.0e9
    y_nm = nodes[:, 1] * 1.0e9
    triangles = np.asarray(mesh.triangles, dtype=np.int64)
    z = np.asarray(values, dtype=float).reshape(-1)

    output = _prepare_output(output_path)

    fig, ax = plt.subplots(figsize=(7.0, 3.2), constrained_layout=True)
    tri = mtri.Triangulation(x_nm, y_nm, triangles)
    im = ax.tripcolor(tri, z, shading="gouraud", vmin=vmin, vmax=vmax)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(label)

    ax.set_title(title)
    ax.set_xlabel("x [nm]")
    ax.set_ylabel("y [nm]")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(False)

    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    return output


# ---------------------------------------------------------------------------
# Relaxation-history diagnostics
# ---------------------------------------------------------------------------

def plot_ss_boundary_currents(
    summary: dict,
    output_path: str | Path,
    *,
    dpi: int = 480,
) -> Path:
    """Plot integrated terminal and transverse boundary currents."""
    boundary = dict(summary["boundary_currents_A"])
    labels = ["left", "right", "bottom", "top"]
    values = [boundary.get(f"{name}_A", 0.0) for name in labels]

    output = _prepare_output(output_path)

    fig, ax = plt.subplots(figsize=(6.2, 3.2), constrained_layout=True)
    ax.bar(labels, values)
    ax.axhline(0.0, linewidth=0.8)
    ax.set_title("OE7 SS: integrated boundary currents")
    ax.set_ylabel("current [A]")
    ax.grid(False)

    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    return output


def plot_ss_relaxation_history(
    history: dict,
    output_path: str | Path,
    *,
    dpi: int = 480,
) -> Path:
    """Plot stationary relaxation residual history."""
    t_ps = _history_array(history, "t_s") / 1.0e-12
    eta = _history_array(history, "eta_R")
    residual = _history_array(history, "current_residual")
    voltage = np.abs(_history_array(history, "terminal_voltage_V"))

    output = _prepare_output(output_path)

    fig, ax = plt.subplots(figsize=(6.4, 3.4), constrained_layout=True)
    if t_ps.size:
        ax.semilogy(t_ps, np.maximum(eta, 1.0e-300), label=r"$\eta_R$")
        ax.semilogy(t_ps, np.maximum(residual, 1.0e-300), label=r"$\epsilon_{\nabla\cdot j}$")
        ax.semilogy(t_ps, np.maximum(voltage, 1.0e-300), label=r"$|V_{\rm TDGL}|$ [V]")

    ax.set_title("OE7 SS: relaxation diagnostics")
    ax.set_xlabel("t [ps]")
    ax.set_ylabel("diagnostic value")
    ax.grid(False)
    ax.legend(frameon=False)

    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    return output


def plot_ss_relaxation_diagnostics(
    history: dict,
    output_path: str | Path,
    *,
    delta0_J: float | None = None,
    voltage_scale_V: float | None = None,
    j_scale_A_m2: float | None = None,
    dpi: int = 480,
) -> Path:
    """Plot the two-panel stiff/physical monitor diagnostic used in OE7.

    This is deliberately tolerant to missing keys, because not every SS run stores
    the same history arrays while the OE7 backend is still evolving.
    """
    t_ps = _history_array(history, "t_s") / 1.0e-12
    output = _prepare_output(output_path)

    amp_rel = _history_first_available(
        history,
        [
            "max_delta_abs2_rel",
            "max_delta2_rel",
            "max_abs_delta2_rel",
            "maxΔ|Δ|2_over_delta0_2",
        ],
    )
    div = _history_first_available(
        history,
        [
            "current_residual",
            "rms_div_j",
            "rms_divj_rel",
            "divergence_residual",
        ],
    )
    xpb = _history_first_available(
        history,
        [
            "max_pairbreaking_ratio",
            "pairbreaking_ratio_max",
            "max_x_pb",
            "xpb_max",
        ],
    )

    voltage = np.abs(_history_first_available(
        history,
        ["terminal_voltage_V", "V_TDGL_V", "Vtdgl_V", "voltage_V"],
    ))
    delta_min = _history_first_available(
        history,
        [
            "min_delta_over_delta0",
            "delta_min_over_delta0",
            "min_abs_delta_over_delta0",
        ],
    )
    normal_fraction = _history_first_available(
        history,
        [
            "normal_current_fraction_max",
            "max_normal_current_fraction",
            "jn_over_j_max",
            "max_jn_over_j",
        ],
    )

    fig, axes = plt.subplots(2, 1, figsize=(12.5, 7.8), constrained_layout=True)

    ax = axes[0]
    if t_ps.size:
        _plot_semilogy_if_present(ax, t_ps, amp_rel, r"max $\Delta|\Delta|^2/\Delta_0^2$")
        _plot_semilogy_if_present(ax, t_ps, div, r"rms div(j)")
        _plot_semilogy_if_present(ax, t_ps, xpb, r"max $x_{\rm pb}$")
    ax.set_title("OE7 SS: stiff relaxation diagnostics")
    ax.set_ylabel("diagnostic value")
    ax.grid(False)
    ax.legend(frameon=False, loc="best")

    ax = axes[1]
    if t_ps.size:
        if voltage.size:
            scale = voltage_scale_V if voltage_scale_V is not None else max(_safe_nanmax(voltage), 1.0e-30)
            ax.plot(t_ps, voltage / scale, label=rf"$|V_{{\rm TDGL}}|$ / {scale * 1.0e3:.3g} mV")
        _plot_if_present(ax, t_ps, delta_min, r"min $|\Delta|/\Delta_0$")
        _plot_if_present(ax, t_ps, normal_fraction, r"max $|j_n|$/max $|j|$")
    ax.set_title("OE7 SS: normalized physical monitors")
    ax.set_xlabel("t [ps]")
    ax.set_ylabel("normalized value")
    ax.set_ylim(bottom=0.0)
    ax.grid(False)
    ax.legend(frameon=False, loc="best")

    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    return output


# ---------------------------------------------------------------------------
# Snapshot mosaic diagnostics
# ---------------------------------------------------------------------------

def plot_ss_normal_current_density_snapshots(
    *,
    mesh,
    history: dict,
    output_path: str | Path,
    state=None,
    ops=None,
    max_snapshots: int = 3,
    dpi: int = 480,
) -> Path:
    """Plot normal-current-density snapshots as field magnitude plus arrows.

    The function accepts several possible history key names and falls back to the
    final state if only the final current field exists.
    """
    jx_hist = _history_first_available(
        history,
        [
            "node_jn_x_A_m2",
            "jn_x_A_m2",
            "snapshot_node_jn_x_A_m2",
            "node_jn_x_snapshots_A_m2",
        ],
    )
    jy_hist = _history_first_available(
        history,
        [
            "node_jn_y_A_m2",
            "jn_y_A_m2",
            "snapshot_node_jn_y_A_m2",
            "node_jn_y_snapshots_A_m2",
        ],
    )

    if (not jx_hist.size or not jy_hist.size) and state is not None:
        currents = getattr(state, "currents", state)
        if hasattr(currents, "node_jn_x_A_m2") and hasattr(currents, "node_jn_y_A_m2"):
            jx_hist = np.asarray(getattr(currents, "node_jn_x_A_m2"), dtype=float)[None, :]
            jy_hist = np.asarray(getattr(currents, "node_jn_y_A_m2"), dtype=float)[None, :]

    if not jx_hist.size or not jy_hist.size:
        raise ValueError("No normal-current node-vector snapshots found in history/state.")

    return _plot_node_vector_snapshots(
        mesh=mesh,
        vx_history=jx_hist,
        vy_history=jy_hist,
        t_s=_history_array(history, "snapshot_t_s", fallback=_history_array(history, "t_s")),
        output_path=output_path,
        title="OE7 SS: normal-current-density snapshots",
        cbar_label=r"$|j_n|/j_{\rm avg}$",
        max_snapshots=max_snapshots,
        dpi=dpi,
    )


def plot_ss_total_current_density_snapshots(
    *,
    mesh,
    history: dict,
    output_path: str | Path,
    state=None,
    ops=None,
    max_snapshots: int = 3,
    dpi: int = 480,
) -> Path:
    """Plot total-current-density snapshots as field magnitude plus arrows."""
    jx_hist = _history_first_available(
        history,
        [
            "node_jtot_x_A_m2",
            "jtot_x_A_m2",
            "snapshot_node_jtot_x_A_m2",
            "node_jtot_x_snapshots_A_m2",
        ],
    )
    jy_hist = _history_first_available(
        history,
        [
            "node_jtot_y_A_m2",
            "jtot_y_A_m2",
            "snapshot_node_jtot_y_A_m2",
            "node_jtot_y_snapshots_A_m2",
        ],
    )

    if (not jx_hist.size or not jy_hist.size) and state is not None:
        currents = getattr(state, "currents", state)
        if hasattr(currents, "node_jtot_x_A_m2") and hasattr(currents, "node_jtot_y_A_m2"):
            jx_hist = np.asarray(getattr(currents, "node_jtot_x_A_m2"), dtype=float)[None, :]
            jy_hist = np.asarray(getattr(currents, "node_jtot_y_A_m2"), dtype=float)[None, :]

    if not jx_hist.size or not jy_hist.size:
        raise ValueError("No total-current node-vector snapshots found in history/state.")

    return _plot_node_vector_snapshots(
        mesh=mesh,
        vx_history=jx_hist,
        vy_history=jy_hist,
        t_s=_history_array(history, "snapshot_t_s", fallback=_history_array(history, "t_s")),
        output_path=output_path,
        title="OE7 SS: total-current-density snapshots",
        cbar_label=r"$|j|/j_{\rm avg}$",
        max_snapshots=max_snapshots,
        dpi=dpi,
    )


def _plot_node_vector_snapshots(
    *,
    mesh,
    vx_history,
    vy_history,
    t_s,
    output_path: str | Path,
    title: str,
    cbar_label: str,
    max_snapshots: int,
    dpi: int,
) -> Path:
    nodes = np.asarray(mesh.nodes, dtype=float)
    triangles = np.asarray(mesh.triangles, dtype=np.int64)
    x_nm = nodes[:, 0] * 1.0e9
    y_nm = nodes[:, 1] * 1.0e9

    vx = _as_snapshot_matrix(vx_history, nodes.shape[0])
    vy = _as_snapshot_matrix(vy_history, nodes.shape[0])
    idx = _snapshot_indices(vx.shape[0], max_snapshots)

    t = np.asarray(t_s, dtype=float).reshape(-1)
    if t.size != vx.shape[0]:
        t = np.arange(vx.shape[0], dtype=float)

    mag_all = np.sqrt(vx**2 + vy**2)
    javg = max(float(np.nanmean(mag_all[np.isfinite(mag_all)])) if np.any(np.isfinite(mag_all)) else 0.0, 1.0e-30)

    fig, axes = plt.subplots(
        1,
        len(idx),
        figsize=(5.0 * len(idx) + 1.2, 3.6),
        squeeze=False,
        constrained_layout=True,
    )
    tri = mtri.Triangulation(x_nm, y_nm, triangles)

    last_im = None
    for ax, k in zip(axes[0], idx):
        mag = np.sqrt(vx[k] ** 2 + vy[k] ** 2) / javg
        vmax = max(_safe_nanmax(mag), 1.0e-30)
        last_im = ax.tripcolor(tri, mag, shading="gouraud", vmin=0.0, vmax=vmax)

        step = max(1, nodes.shape[0] // 130)
        scale = _safe_nanmax(np.sqrt(vx[k] ** 2 + vy[k] ** 2))
        if np.isfinite(scale) and scale > 0.0:
            ax.quiver(
                x_nm[::step],
                y_nm[::step],
                vx[k, ::step] / scale,
                vy[k, ::step] / scale,
                angles="xy",
                scale_units="xy",
                scale=0.040,
                width=0.0024,
            )

        ax.set_title(f"t = {_format_ps(t[k])} ps")
        ax.set_xlabel("x [nm]")
        ax.set_ylabel("y [nm]")
        ax.set_aspect("equal", adjustable="box")
        ax.grid(False)

    if last_im is not None:
        cbar = fig.colorbar(last_im, ax=axes.ravel().tolist(), shrink=0.86)
        cbar.set_label(cbar_label)

    fig.suptitle(title)
    output = _prepare_output(output_path)
    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    return output


# ---------------------------------------------------------------------------
# Boundary/current-reconstruction diagnostics
# ---------------------------------------------------------------------------

def plot_ss_boundary_current_reconstruction_comparison(
    *,
    mesh,
    edge_data,
    ops,
    state,
    output_path: str | Path,
    target_current_A: float | None = None,
    thickness_m: float,
    dpi: int = 480,
) -> Path:
    """Compare terminal currents from different diagnostic reconstructions."""
    node_avg = boundary_currents_from_node_vectors(
        mesh=mesh,
        edge_data=edge_data,
        jx_A_m2=state.currents.node_jtot_x_A_m2,
        jy_A_m2=state.currents.node_jtot_y_A_m2,
        thickness_m=thickness_m,
    )

    ls = boundary_currents_from_edge_scalar_least_squares(
        mesh=mesh,
        edge_data=edge_data,
        ops=ops,
        edge_current_i_to_j=state.currents.edge_jtot_A_m2,
        thickness_m=thickness_m,
    )

    labels = [
        "left\nnode avg",
        "left\nLS",
        "right\nnode avg",
        "right\nLS",
    ]
    values = [
        node_avg.get("left_A", 0.0),
        ls.get("left_A", 0.0),
        node_avg.get("right_A", 0.0),
        ls.get("right_A", 0.0),
    ]

    output = _prepare_output(output_path)

    fig, ax = plt.subplots(figsize=(6.8, 3.4), constrained_layout=True)
    ax.bar(labels, values)
    ax.axhline(0.0, linewidth=0.8)

    if target_current_A is not None:
        target = float(target_current_A)
        ax.axhline(+target, linestyle="--", linewidth=0.9, label=r"$+I_{\rm target}$")
        ax.axhline(-target, linestyle="--", linewidth=0.9, label=r"$-I_{\rm target}$")
        ax.legend(frameon=False)

    ax.set_title("OE7 SS: boundary-current reconstruction")
    ax.set_ylabel("current [A]")
    ax.grid(False)

    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    return output


def plot_ss_transport_current_profile(
    *,
    mesh,
    ops,
    state,
    output_path: str | Path,
    target_current_A: float | None = None,
    thickness_m: float,
    n_bins: int = 41,
    dpi: int = 480,
) -> Path:
    """Plot longitudinal transport-current profile from LS reconstructed jx."""
    jx_ls, _ = edge_scalar_to_node_vector_least_squares(
        state.currents.edge_jtot_A_m2,
        ops,
    )

    x_m, I_A = strip_transport_current_profile_from_node_vectors(
        mesh=mesh,
        jx_A_m2=jx_ls,
        thickness_m=thickness_m,
        n_bins=n_bins,
    )

    output = _prepare_output(output_path)

    fig, ax = plt.subplots(figsize=(6.6, 3.4), constrained_layout=True)
    ax.plot(x_m * 1.0e9, I_A, marker="o", markersize=2.5, linewidth=1.0)

    if target_current_A is not None:
        ax.axhline(float(target_current_A), linestyle="--", linewidth=0.9, label=r"$I_{\rm target}$")
        ax.legend(frameon=False)

    ax.set_title("OE7 SS: LS transport-current profile")
    ax.set_xlabel("x [nm]")
    ax.set_ylabel("I(x) [A]")
    ax.grid(False)

    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    return output


# ---------------------------------------------------------------------------
# Edge-current zoom diagnostics
# ---------------------------------------------------------------------------

def plot_ss_edge_current_diagnostics(
    mesh,
    edge_data=None,
    history: dict | None = None,
    output_dir: str | Path | None = None,
    *,
    state=None,
    ops=None,
    output_path: str | Path | None = None,
    edge_current_i_to_j=None,
    max_snapshots: int = 3,
    dpi: int = 480,
) -> list[Path] | Path:
    """Plot node-to-node edge currents in four diagnostic zoom windows.

    This is the main diagnostic for the directory

        diagnostics/edge_current_diagnostics

    Each figure has four columns: a top insulating-edge window, a left
    terminal-edge window, a bottom-right corner window, and an interior strip
    window.  The upper row shows the signed current projection on the stored
    edge orientation.  The lower row shows the same region in absolute value,
    normalized by the global mean absolute edge current for that snapshot.

    Parameters are permissive on purpose.  During OE7, the pipeline has stored
    edge snapshots under a few different key names.  This routine first tries
    the explicit `edge_current_i_to_j`, then history arrays, then final state.
    """
    edges = _edge_array_from_any(mesh, edge_data=edge_data, ops=ops)
    edge_values, t_s = _extract_edge_current_snapshots(
        mesh=mesh,
        edges=edges,
        history=history,
        state=state,
        explicit=edge_current_i_to_j,
    )

    edge_values = _as_snapshot_matrix(edge_values, edges.shape[0])
    idx = _snapshot_indices(edge_values.shape[0], max_snapshots)

    if output_path is not None:
        k = int(idx[-1])
        return plot_ss_edge_current_zoom_windows(
            mesh=mesh,
            edge_data=edge_data,
            ops=ops,
            edge_current_i_to_j=edge_values[k],
            output_path=output_path,
            t_s=float(t_s[k]) if t_s.size == edge_values.shape[0] else None,
            dpi=dpi,
        )

    if output_dir is None:
        raise ValueError("Either output_dir or output_path must be provided.")

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    outputs: list[Path] = []
    for n, k in enumerate(idx):
        t_val = float(t_s[k]) if t_s.size == edge_values.shape[0] else None
        suffix = f"{n:02d}"
        if t_val is not None:
            suffix += f"_t_{t_val / 1.0e-12:.6g}ps"
        out = out_dir / f"edge_current_zoom_{suffix}.png"
        outputs.append(
            plot_ss_edge_current_zoom_windows(
                mesh=mesh,
                edge_data=edge_data,
                ops=ops,
                edge_current_i_to_j=edge_values[k],
                output_path=out,
                t_s=t_val,
                dpi=dpi,
            )
        )

    return outputs


def plot_ss_edge_current_zoom_windows(
    *,
    mesh,
    edge_data=None,
    ops=None,
    edge_current_i_to_j,
    output_path: str | Path,
    t_s: float | None = None,
    dpi: int = 480,
) -> Path:
    """Plot a single edge-current snapshot in four zoom windows."""
    nodes = np.asarray(mesh.nodes, dtype=float)
    nodes_nm = nodes * 1.0e9
    edges = _edge_array_from_any(mesh, edge_data=edge_data, ops=ops)
    values = np.asarray(edge_current_i_to_j, dtype=float).reshape(-1)

    if values.size != edges.shape[0]:
        raise ValueError(
            f"edge_current_i_to_j has length {values.size}, but {edges.shape[0]} edges were found."
        )

    regions = _edge_current_zoom_regions(mesh)
    n_regions = max(1, len(regions))
    fig_width = 4.45 * n_regions + 0.95

    fig, axes = plt.subplots(
        2,
        n_regions,
        figsize=(fig_width, 6.4),
        squeeze=False,
        constrained_layout=False,
    )
    fig.subplots_adjust(
        left=0.047,
        right=0.930,
        bottom=0.090,
        top=0.875,
        wspace=0.23,
        hspace=0.30,
    )

    abs_values = np.abs(values)
    signed_vmax = max(_safe_nanmax(abs_values), 1.0e-30)
    mean_abs = max(float(np.nanmean(abs_values[np.isfinite(abs_values)])) if np.any(np.isfinite(abs_values)) else 0.0, 1.0e-30)
    normalized_abs = abs_values / mean_abs
    abs_vmax = max(_safe_percentile(normalized_abs, 99.0), 1.0e-30)

    signed_mappable = None
    abs_mappable = None

    for col, region in enumerate(regions):
        mask = _edge_mask_for_region(nodes_nm, edges, region)
        highlight_edge = _select_representative_edge(
            nodes_nm=nodes_nm,
            edges=edges,
            mask=mask,
            region=region,
        )

        ax = axes[0, col]
        signed_mappable = _draw_edge_zoom(
            ax=ax,
            nodes_nm=nodes_nm,
            edges=edges,
            values=values,
            mask=mask,
            title=f"{region['name']}\nsigned edge current",
            label="signed",
            cmap="coolwarm",
            vmin=-signed_vmax,
            vmax=signed_vmax,
            highlight_edge=highlight_edge,
        )

        ax = axes[1, col]
        abs_mappable = _draw_edge_zoom(
            ax=ax,
            nodes_nm=nodes_nm,
            edges=edges,
            values=normalized_abs,
            mask=mask,
            title=f"{region['name']}\n|edge current| / mean",
            label="absolute",
            cmap="viridis",
            vmin=0.0,
            vmax=abs_vmax,
            highlight_edge=highlight_edge,
        )

    if signed_mappable is not None:
        cbar = fig.colorbar(signed_mappable, ax=axes[0, :].ravel().tolist(), shrink=0.82, pad=0.012)
        cbar.set_label(r"$j_{ij}$ [A m$^{-2}$]")

    if abs_mappable is not None:
        cbar = fig.colorbar(abs_mappable, ax=axes[1, :].ravel().tolist(), shrink=0.82, pad=0.012)
        cbar.set_label(r"$|j_{ij}|/\langle |j_{ij}| \rangle$")

    if t_s is None:
        fig.suptitle("OE7 SS: edge-current diagnostic zooms")
    else:
        fig.suptitle(f"OE7 SS: edge-current diagnostic zooms, t = {_format_ps(t_s)} ps")

    output = _prepare_output(output_path)
    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    return output


# Backward-compatible aliases for pipeline/import experiments during OE7.
plot_ss_edge_current_zoom_diagnostics = plot_ss_edge_current_diagnostics
plot_ss_edge_current_snapshots = plot_ss_edge_current_diagnostics


def _edge_current_zoom_regions(mesh) -> list[dict[str, Any]]:
    nodes_nm = np.asarray(mesh.nodes, dtype=float) * 1.0e9
    xmin = float(np.nanmin(nodes_nm[:, 0]))
    xmax = float(np.nanmax(nodes_nm[:, 0]))
    ymin = float(np.nanmin(nodes_nm[:, 1]))
    ymax = float(np.nanmax(nodes_nm[:, 1]))

    lx = max(xmax - xmin, 1.0e-30)
    wy = max(ymax - ymin, 1.0e-30)

    return [
        {
            "name": "top insulating edge",
            "slug": "top_edge",
            "center": (0.50 * (xmin + xmax), ymax),
            "half_width": 0.105 * lx,
            "half_height": 0.115 * wy,
            "kind": "top",
        },
        {
            "name": "left terminal edge",
            "slug": "left_terminal",
            "center": (xmin, 0.50 * (ymin + ymax)),
            "half_width": 0.055 * lx,
            "half_height": 0.165 * wy,
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
        {
            "name": "interior strip",
            "slug": "interior_strip",
            "center": (0.50 * (xmin + xmax), 0.50 * (ymin + ymax)),
            "half_width": 0.075 * lx,
            "half_height": 0.120 * wy,
            "kind": "interior",
        },
    ]


def _edge_mask_for_region(nodes_nm: np.ndarray, edges: np.ndarray, region: dict[str, Any]) -> np.ndarray:
    center_x, center_y = region["center"]
    half_width = float(region["half_width"])
    half_height = float(region["half_height"])

    p0 = nodes_nm[edges[:, 0]]
    p1 = nodes_nm[edges[:, 1]]
    mx = 0.5 * (p0[:, 0] + p1[:, 0])
    my = 0.5 * (p0[:, 1] + p1[:, 1])

    mask = (
        (np.abs(mx - center_x) <= half_width)
        & (np.abs(my - center_y) <= half_height)
    )

    if np.any(mask):
        return mask

    # Robust fallback for very coarse meshes: choose nearest local edges.
    dist = ((mx - center_x) ** 2 + (my - center_y) ** 2) ** 0.5
    n_keep = min(max(8, edges.shape[0] // 120), edges.shape[0])
    idx = np.argsort(dist)[:n_keep]
    mask = np.zeros(edges.shape[0], dtype=bool)
    mask[idx] = True
    return mask


def _select_representative_edge(
    *,
    nodes_nm: np.ndarray,
    edges: np.ndarray,
    mask: np.ndarray,
    region: dict[str, Any],
) -> int | None:
    edge_idx = np.flatnonzero(mask)
    if edge_idx.size == 0:
        return None

    local_edges = edges[edge_idx]
    p0 = nodes_nm[local_edges[:, 0]]
    p1 = nodes_nm[local_edges[:, 1]]

    xi = p0[:, 0]
    yi = p0[:, 1]
    xj = p1[:, 0]
    yj = p1[:, 1]
    mx = 0.5 * (xi + xj)
    my = 0.5 * (yi + yj)
    dx = xj - xi
    dy = yj - yi

    ell = np.sqrt(dx**2 + dy**2)
    ell = np.maximum(ell, 1.0e-300)

    xmin = float(np.nanmin(nodes_nm[:, 0]))
    xmax = float(np.nanmax(nodes_nm[:, 0]))
    ymin = float(np.nanmin(nodes_nm[:, 1]))
    ymax = float(np.nanmax(nodes_nm[:, 1]))
    lx = max(xmax - xmin, 1.0e-30)
    wy = max(ymax - ymin, 1.0e-30)
    tol_x = 0.020 * lx
    tol_y = 0.020 * wy

    center_x, center_y = region["center"]
    kind = region.get("kind", "")

    if kind == "top":
        touches_top = (np.abs(yi - ymax) < tol_y) | (np.abs(yj - ymax) < tol_y)
        score = touches_top.astype(float) * 10.0
        score += 0.5 * np.abs(dx) / ell
        score -= 0.02 * np.abs(mx - center_x)
    elif kind == "left":
        touches_left = (np.abs(xi - xmin) < tol_x) | (np.abs(xj - xmin) < tol_x)
        score = touches_left.astype(float) * 10.0
        score += 0.5 * np.abs(dy) / ell
        score -= 0.02 * np.abs(my - center_y)
    elif kind == "bottom_right":
        touches_right = (np.abs(xi - xmax) < tol_x) | (np.abs(xj - xmax) < tol_x)
        touches_bottom = (np.abs(yi - ymin) < tol_y) | (np.abs(yj - ymin) < tol_y)
        score = (touches_right | touches_bottom).astype(float) * 10.0
        score += 0.5 * (np.abs(dx) + np.abs(dy)) / ell
        score -= 0.02 * ((mx - center_x) ** 2 + (my - center_y) ** 2) ** 0.5
    else:
        # Interior diagnostic: pick a representative edge close to the middle
        # of the strip, with a mild preference for transport-aligned edges.
        # This avoids highlighting boundary artifacts in the fourth column.
        d = ((mx - center_x) ** 2 + (my - center_y) ** 2) ** 0.5
        d_scale = max(_safe_nanmax(d), 1.0e-300)
        score = 0.35 * np.abs(dx) / ell - d / d_scale

    return int(edge_idx[int(np.nanargmax(score))])


def _draw_edge_zoom(
    *,
    ax,
    nodes_nm: np.ndarray,
    edges: np.ndarray,
    values: np.ndarray,
    mask: np.ndarray,
    title: str,
    label: str,
    cmap: str,
    vmin: float,
    vmax: float,
    highlight_edge: int | None,
):
    segments = nodes_nm[edges[mask]]
    vals = np.asarray(values, dtype=float)[mask]

    collection = LineCollection(
        segments,
        array=vals,
        cmap=cmap,
        norm=plt.Normalize(vmin=vmin, vmax=vmax),
        linewidths=1.7,
        alpha=0.95,
    )
    ax.add_collection(collection)

    if highlight_edge is not None:
        hseg = nodes_nm[edges[[highlight_edge]]]
        hcollection = LineCollection(
            hseg,
            colors="black",
            linewidths=3.0,
            alpha=0.95,
        )
        ax.add_collection(hcollection)

    if segments.size:
        x = segments[:, :, 0].reshape(-1)
        y = segments[:, :, 1].reshape(-1)
        dx = max(float(np.nanmax(x) - np.nanmin(x)), 1.0e-30)
        dy = max(float(np.nanmax(y) - np.nanmin(y)), 1.0e-30)
        ax.set_xlim(float(np.nanmin(x) - 0.12 * dx), float(np.nanmax(x) + 0.12 * dx))
        ax.set_ylim(float(np.nanmin(y) - 0.18 * dy), float(np.nanmax(y) + 0.18 * dy))

    ax.set_title(title)
    ax.set_xlabel("x [nm]")
    ax.set_ylabel("y [nm]")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(False)

    return collection


# ---------------------------------------------------------------------------
# Phase/mesh helpers
# ---------------------------------------------------------------------------

def _unwrap_phase_by_x(mesh, theta_wrapped: np.ndarray) -> np.ndarray:
    psi = np.exp(1j * np.asarray(theta_wrapped, dtype=float))
    edges = _edges_from_triangles(mesh.triangles)
    return unwrap_phase_graph(psi, edges)


def _edges_from_triangles(triangles: np.ndarray) -> np.ndarray:
    triangles = np.asarray(triangles, dtype=np.int64)
    pairs = np.vstack(
        [
            triangles[:, [0, 1]],
            triangles[:, [1, 2]],
            triangles[:, [2, 0]],
        ]
    )
    pairs.sort(axis=1)
    return np.unique(pairs, axis=0)


def _edge_array_from_any(mesh, *, edge_data=None, ops=None) -> np.ndarray:
    if edge_data is not None and hasattr(edge_data, "edges"):
        return np.asarray(edge_data.edges, dtype=np.int64)

    if ops is not None:
        for name in ("edges", "edge_nodes", "edge_node_indices"):
            if hasattr(ops, name):
                return np.asarray(getattr(ops, name), dtype=np.int64)

    if hasattr(mesh, "edges"):
        return np.asarray(mesh.edges, dtype=np.int64)

    return _edges_from_triangles(np.asarray(mesh.triangles, dtype=np.int64))


# ---------------------------------------------------------------------------
# Data extraction helpers
# ---------------------------------------------------------------------------

def _extract_edge_current_snapshots(
    *,
    mesh,
    edges: np.ndarray,
    history: dict | None,
    state,
    explicit,
) -> tuple[np.ndarray, np.ndarray]:
    n_edges = edges.shape[0]

    if explicit is not None:
        arr = np.asarray(explicit, dtype=float)
        t = _history_array(history or {}, "snapshot_t_s", fallback=_history_array(history or {}, "t_s"))
        if t.size == 0:
            t = np.arange(arr.shape[0] if arr.ndim == 2 else 1, dtype=float)
        return arr, t

    if history is not None:
        candidates = [
            "edge_jtot_A_m2",
            "edge_current_A_m2",
            "snapshot_edge_jtot_A_m2",
            "edge_jtot_snapshots_A_m2",
            "edge_jtot_A_m2_snapshots",
            "edge_current_i_to_j_A_m2",
            "edge_current_i_to_j",
        ]
        for key in candidates:
            if key not in history:
                continue
            arr = np.asarray(history[key], dtype=float)
            if arr.size == 0:
                continue
            if arr.ndim == 1 and arr.size == n_edges:
                t = _history_array(history, "snapshot_t_s", fallback=_history_array(history, "t_s"))
                return arr[None, :], t[:1] if t.size else np.array([0.0])
            if arr.ndim == 2 and n_edges in arr.shape:
                if arr.shape[1] != n_edges and arr.shape[0] == n_edges:
                    arr = arr.T
                t = _history_array(history, "snapshot_t_s", fallback=_history_array(history, "t_s"))
                if t.size != arr.shape[0]:
                    t = np.arange(arr.shape[0], dtype=float)
                return arr, t

    if state is not None:
        currents = getattr(state, "currents", state)
        for key in ("edge_jtot_A_m2", "edge_current_A_m2", "edge_current_i_to_j"):
            if hasattr(currents, key):
                arr = np.asarray(getattr(currents, key), dtype=float)
                if arr.size == n_edges:
                    return arr[None, :], np.array([0.0])

    raise ValueError("No edge current snapshots found in explicit/history/state inputs.")


def _history_array(history: dict | None, key: str, fallback=None) -> np.ndarray:
    if history is not None and key in history:
        return np.asarray(history[key], dtype=float)
    if fallback is None:
        return np.asarray([], dtype=float)
    return np.asarray(fallback, dtype=float)


def _history_first_available(history: dict | None, keys: Iterable[str]) -> np.ndarray:
    if history is None:
        return np.asarray([], dtype=float)

    for key in keys:
        if key in history:
            arr = np.asarray(history[key], dtype=float)
            if arr.size:
                return arr

    return np.asarray([], dtype=float)


def _as_snapshot_matrix(arr, n_values: int) -> np.ndarray:
    out = np.asarray(arr, dtype=float)

    if out.ndim == 1:
        if out.size != n_values:
            raise ValueError(f"Expected {n_values} values, got {out.size}.")
        return out[None, :]

    if out.ndim != 2:
        raise ValueError(f"Expected 1D or 2D snapshot array, got shape {out.shape}.")

    if out.shape[1] == n_values:
        return out

    if out.shape[0] == n_values:
        return out.T

    raise ValueError(f"Could not orient snapshot matrix with shape {out.shape} and n_values={n_values}.")


def _snapshot_indices(n_snapshots: int, max_snapshots: int) -> np.ndarray:
    n = int(max(1, n_snapshots))
    m = int(max(1, min(max_snapshots, n)))
    if m == 1:
        return np.asarray([n - 1], dtype=int)
    return np.unique(np.linspace(0, n - 1, m).round().astype(int))


# ---------------------------------------------------------------------------
# Small plotting helpers
# ---------------------------------------------------------------------------

def _plot_semilogy_if_present(ax, x: np.ndarray, y: np.ndarray, label: str) -> None:
    y = np.asarray(y, dtype=float)
    if not y.size:
        return
    if y.size != x.size:
        y = _resample_like(y, x.size)
    final = y[-1] if y.size else np.nan
    ax.semilogy(x, np.maximum(np.abs(y), 1.0e-300), label=f"{label}, final={final:.3e}")


def _plot_if_present(ax, x: np.ndarray, y: np.ndarray, label: str) -> None:
    y = np.asarray(y, dtype=float)
    if not y.size:
        return
    if y.size != x.size:
        y = _resample_like(y, x.size)
    final = y[-1] if y.size else np.nan
    ax.plot(x, y, label=f"{label}, final={final:.4g}")


def _resample_like(y: np.ndarray, n: int) -> np.ndarray:
    y = np.asarray(y, dtype=float).reshape(-1)
    if y.size == n:
        return y
    if y.size == 1:
        return np.full(n, y[0], dtype=float)
    xp = np.linspace(0.0, 1.0, y.size)
    x = np.linspace(0.0, 1.0, n)
    return np.interp(x, xp, y)


def _format_ps(t_s: float) -> str:
    return f"{float(t_s) / 1.0e-12:.6g}"


def _prepare_output(output_path: str | Path) -> Path:
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    return output


def _safe_nanmax(x) -> float:
    arr = np.asarray(x, dtype=float)
    if arr.size == 0:
        return 0.0
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return 0.0
    return float(np.max(finite))


def _safe_percentile(x, q: float) -> float:
    arr = np.asarray(x, dtype=float)
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return 0.0
    return float(np.percentile(finite, q))
