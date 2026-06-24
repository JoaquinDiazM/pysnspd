"""Diagnostic plots for OE7 stationary gTDGL/Poisson runs."""
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
    """Plot all snapshot fields currently present in the relaxation history.

    Important:
    This function does not modify or infer solver data. If the solver currently
    stores only phi snapshots, only phi snapshots are emitted. Extra fields are
    plotted automatically if future solver diagnostics add them to history.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    saved: dict[str, Path] = {}

    maybe_plots = [
        ("phi", plot_ss_phi_snapshots, output_dir / "ss_phi_snapshots.png"),
        ("delta", plot_ss_delta_snapshots, output_dir / "ss_delta_snapshots.png"),
        ("phase", plot_ss_phase_snapshots, output_dir / "ss_phase_snapshots.png"),
        ("current_density", plot_ss_current_density_snapshots, output_dir / "ss_current_density_snapshots.png"),
        ("divergence", plot_ss_divergence_snapshots, output_dir / "ss_divergence_snapshots.png"),
        ("pairbreaking", plot_ss_pairbreaking_snapshots, output_dir / "ss_pairbreaking_snapshots.png"),
    ]

    for name, func, path in maybe_plots:
        try:
            saved[name] = func(mesh, history, path, dpi=dpi, ncols=ncols)
        except KeyError:
            continue
        except ValueError:
            continue

    return saved


def plot_ss_phi_snapshots(
    mesh,
    history: dict,
    output_path: str | Path,
    *,
    dpi: int = 480,
    ncols: int = 3,
) -> Path:
    """Plot electrostatic-potential snapshots during OE7 relaxation."""
    phi = _require_history_array(history, ("phi_snapshot_V",))
    t_s = _snapshot_times(history, ("phi_snapshot_t_s",), phi.shape[0])

    return _plot_snapshot_grid(
        mesh,
        phi,
        t_s,
        output_path,
        title="OE7 SS: electrostatic potential snapshots",
        label="φ [V]",
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
    """Plot |Delta| snapshots in meV if present in history."""
    if "delta_snapshot_meV" in history:
        delta_meV = np.asarray(history["delta_snapshot_meV"], dtype=float)
        t_s = _snapshot_times(history, ("delta_snapshot_t_s", "phi_snapshot_t_s"), delta_meV.shape[0])
    else:
        psi = _psi_snapshots_from_history(history)
        delta_meV = np.abs(psi) / MEV_J
        t_s = _snapshot_times(history, ("psi_snapshot_t_s", "delta_snapshot_t_s", "phi_snapshot_t_s"), delta_meV.shape[0])

    return _plot_snapshot_grid(
        mesh,
        delta_meV,
        t_s,
        output_path,
        title="OE7 SS: |Δ| snapshots",
        label="|Δ| [meV]",
        vmin=0.0,
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
    """Plot unwrapped phase snapshots if psi snapshots are present."""
    psi = _psi_snapshots_from_history(history)
    t_s = _snapshot_times(history, ("psi_snapshot_t_s", "phase_snapshot_t_s", "phi_snapshot_t_s"), psi.shape[0])

    edges = _mesh_edges_from_triangles(mesh)
    seed_index = _center_node_index(mesh)

    theta = np.vstack(
        [
            _unwrap_phase_safe(psi[k], edges, seed_index=seed_index)
            for k in range(psi.shape[0])
        ]
    )

    return _plot_snapshot_grid(
        mesh,
        theta,
        t_s,
        output_path,
        title="OE7 SS: graph-unwrapped phase snapshots",
        label="θ [rad]",
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
    """Plot |j| snapshots if current-density snapshots are present."""
    jmag = _optional_history_array(
        history,
        ("jtot_snapshot_mag_A_m2", "jmag_snapshot_A_m2", "current_density_snapshot_A_m2"),
    )

    if jmag is None:
        jx = _require_history_array(
            history,
            ("jtot_snapshot_x_A_m2", "node_jtot_x_snapshot_A_m2", "jx_snapshot_A_m2"),
        )
        jy = _require_history_array(
            history,
            ("jtot_snapshot_y_A_m2", "node_jtot_y_snapshot_A_m2", "jy_snapshot_A_m2"),
        )
        jmag = np.sqrt(jx * jx + jy * jy)

    t_s = _snapshot_times(
        history,
        ("jtot_snapshot_t_s", "current_snapshot_t_s", "phi_snapshot_t_s"),
        jmag.shape[0],
    )

    return _plot_snapshot_grid(
        mesh,
        jmag,
        t_s,
        output_path,
        title="OE7 SS: total current-density snapshots",
        label=r"|j| [A m$^{-2}$]",
        vmin=0.0,
        dpi=dpi,
        ncols=ncols,
    )


def plot_ss_divergence_snapshots(
    mesh,
    history: dict,
    output_path: str | Path,
    *,
    dpi: int = 480,
    ncols: int = 3,
) -> Path:
    """Plot div(j) snapshots if present in history."""
    div = _require_history_array(
        history,
        ("div_jtot_snapshot_A_m3", "node_div_jtot_snapshot_A_m3", "divergence_snapshot_A_m3"),
    )
    t_s = _snapshot_times(
        history,
        ("divergence_snapshot_t_s", "div_jtot_snapshot_t_s", "phi_snapshot_t_s"),
        div.shape[0],
    )

    return _plot_snapshot_grid(
        mesh,
        div,
        t_s,
        output_path,
        title="OE7 SS: finite-volume div(j) snapshots",
        label=r"div(j) [A m$^{-3}$]",
        symmetric=True,
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
    """Plot pair-breaking-ratio snapshots if present in history."""
    chi = _require_history_array(
        history,
        ("pairbreaking_snapshot", "pairbreaking_ratio_snapshot", "node_pairbreaking_snapshot"),
    )
    t_s = _snapshot_times(
        history,
        ("pairbreaking_snapshot_t_s", "pairbreaking_ratio_snapshot_t_s", "phi_snapshot_t_s"),
        chi.shape[0],
    )

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
# Kept intentionally so old imports/tests/scripts do not break.
# The updated pipeline below no longer calls these final-field colormaps.


def plot_ss_state_delta(mesh, state, output_path: str | Path, *, dpi: int = 480) -> Path:
    """Plot the relaxed order-parameter amplitude in meV."""
    delta_meV = np.abs(state.psi_J) / MEV_J
    vmax = max(float(np.nanmax(delta_meV)), 1.0e-30)

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
    """Plot graph-unwrapped phase."""
    theta = _unwrap_phase_safe(
        np.asarray(state.psi_J, dtype=np.complex128),
        np.asarray(_mesh_edges_from_triangles(mesh), dtype=np.int64),
        seed_index=_center_node_index(mesh),
    )

    return _plot_node_scalar(
        mesh,
        theta,
        output_path,
        title="OE7 SS: unwrapped phase θ",
        label="θ [rad]",
        dpi=dpi,
    )


def plot_ss_state_phi(mesh, state, output_path: str | Path, *, dpi: int = 480) -> Path:
    """Plot final electrostatic potential."""
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
    vmax = float(np.nanmax(np.abs(div))) if div.size else 1.0
    vmax = max(vmax, 1.0e-30)

    return _plot_node_scalar(
        mesh,
        div,
        output_path,
        title="OE7 SS: finite-volume div(j)",
        label=r"div(j) [A m$^{-3}$]",
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
    """Plot total current-density magnitude and sparse vectors."""
    del ops

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    nodes = np.asarray(mesh.nodes, dtype=float)
    x_nm = nodes[:, 0] * 1.0e9
    y_nm = nodes[:, 1] * 1.0e9
    triangles = np.asarray(mesh.triangles, dtype=np.int64)
    tri = mtri.Triangulation(x_nm, y_nm, triangles)

    jx = np.asarray(state.currents.node_jtot_x_A_m2, dtype=float)
    jy = np.asarray(state.currents.node_jtot_y_A_m2, dtype=float)
    jmag = np.sqrt(jx * jx + jy * jy)
    vmax = max(float(np.nanmax(jmag)), 1.0e-30)

    fig, ax = plt.subplots(figsize=(8.0, 3.2))
    im = ax.tripcolor(tri, jmag, shading="gouraud", vmin=0.0, vmax=vmax)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(r"|j| [A m$^{-2}$]")

    step = max(1, nodes.shape[0] // 150)
    ax.quiver(
        x_nm[::step],
        y_nm[::step],
        jx[::step],
        jy[::step],
        angles="xy",
        scale_units="xy",
        scale=None,
        width=0.002,
    )

    ax.set_title("OE7 SS: total current density")
    ax.set_xlabel("x [nm]")
    ax.set_ylabel("y [nm]")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(False)

    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output


def plot_ss_pairbreaking_ratio(
    mesh,
    state,
    output_path: str | Path,
    *,
    dpi: int = 480,
) -> Path:
    """Plot chi_pb = xi^2 Q^2 / (1 - T/Tc)."""
    chi = np.asarray(state.currents.node_pairbreaking_ratio, dtype=float)
    finite = chi[np.isfinite(chi)]
    vmax = float(np.nanpercentile(finite, 99.5)) if finite.size else 1.0
    vmax = max(vmax, 1.0)

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    nodes = np.asarray(mesh.nodes, dtype=float)
    x_nm = nodes[:, 0] * 1.0e9
    y_nm = nodes[:, 1] * 1.0e9
    tri = mtri.Triangulation(
        x_nm,
        y_nm,
        np.asarray(mesh.triangles, dtype=np.int64),
    )

    fig, ax = plt.subplots(figsize=(8.0, 3.2))
    im = ax.tripcolor(tri, chi, shading="gouraud", vmin=0.0, vmax=vmax)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(r"$\chi_{\rm pb}=\xi^2Q^2/(1-T/T_c)$")

    if finite.size and np.nanmin(chi) <= 1.0 <= np.nanmax(chi):
        ax.tricontour(tri, chi, levels=[1.0], linewidths=1.0)

    ax.set_title("OE7 SS: pairbreaking diagnostic")
    ax.set_xlabel("x [nm]")
    ax.set_ylabel("y [nm]")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(False)

    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output


# =============================================================================
# Scalar / non-field diagnostics
# =============================================================================


def plot_ss_boundary_currents(
    summary: dict,
    output_path: str | Path,
    *,
    dpi: int = 480,
) -> Path:
    """Plot integrated boundary currents from the final state summary."""
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
    """Plot longitudinal transport-current profile."""
    del ops

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    x_m, current_A = strip_transport_current_profile_from_node_vectors(
        mesh=mesh,
        jx_A_m2=np.asarray(state.currents.node_jtot_x_A_m2, dtype=float),
        thickness_m=float(thickness_m),
        n_bins=int(n_bins),
    )

    fig, ax = plt.subplots(figsize=(8.0, 3.2))
    ax.plot(x_m * 1.0e9, current_A, marker="o", label="node-avg profile")
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


def plot_ss_relaxation_history(
    history: dict,
    output_path: str | Path,
    *,
    dpi: int = 480,
) -> Path:
    """Plot compact relaxation diagnostics in two panels.

    Panel 1:
        log-scale stiff diagnostics: eta_R, current_residual, pairbreaking_max.

    Panel 2:
        normalized linear diagnostics in [0, 1]:
        terminal voltage, min Delta/Delta0, and normal-current RMS fraction.
    """
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    t_ps = np.asarray(history.get("t_s", []), dtype=float) / 1.0e-12

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(8.4, 6.0),
        sharex=True,
    )
    ax_log, ax_lin = axes

    if t_ps.size:
        _plot_log_history_curve(
            ax_log,
            t_ps,
            history,
            "eta_R",
            r"$\eta_R$",
        )
        _plot_log_history_curve(
            ax_log,
            t_ps,
            history,
            "current_residual",
            r"$\epsilon_{\nabla\cdot j}$",
        )
        _plot_log_history_curve(
            ax_log,
            t_ps,
            history,
            "pairbreaking_max",
            r"$\max\chi_{\rm pb}$",
        )

        _plot_normalized_history_curve(
            ax_lin,
            t_ps,
            history,
            "terminal_voltage_V",
            r"$|V_{\rm TDGL}|$",
            unit="V",
            absolute=True,
        )
        _plot_normalized_history_curve(
            ax_lin,
            t_ps,
            history,
            "delta_min_over_delta0",
            r"$\min|\Delta|/\Delta_0$",
            unit="",
            absolute=False,
        )
        _plot_normalized_history_curve(
            ax_lin,
            t_ps,
            history,
            "normal_current_fraction_rms",
            r"$\mathrm{rms}(|j_n|)/\mathrm{rms}(|j|)$",
            unit="",
            absolute=False,
        )

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

    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output


# =============================================================================
# Internal helpers
# =============================================================================


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
    n_nodes = int(nodes.shape[0])
    if arr.shape[1] != n_nodes:
        raise ValueError(
            f"Snapshot array has {arr.shape[1]} nodes, but mesh has {n_nodes} nodes."
        )

    t_s = np.asarray(t_s, dtype=float).reshape(-1)
    if t_s.size != arr.shape[0]:
        raise ValueError(
            f"Snapshot time axis has {t_s.size} entries, but data has {arr.shape[0]} snapshots."
        )

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        finite = np.array([0.0], dtype=float)

    if vmin is None or vmax is None:
        if symmetric:
            if robust_percentile is None:
                scale = float(np.nanmax(np.abs(finite)))
            else:
                scale = float(np.nanpercentile(np.abs(finite), robust_percentile))
            scale = max(scale, 1.0e-30)
            if min_vmax is not None:
                scale = max(scale, float(min_vmax))
            if vmin is None:
                vmin = -scale
            if vmax is None:
                vmax = scale
        else:
            local_vmin = float(np.nanmin(finite))
            if robust_percentile is None:
                local_vmax = float(np.nanmax(finite))
            else:
                local_vmax = float(np.nanpercentile(finite, robust_percentile))
            if min_vmax is not None:
                local_vmax = max(local_vmax, float(min_vmax))
            if local_vmax <= local_vmin:
                pad = max(abs(local_vmax), 1.0) * 1.0e-12
                local_vmin -= pad
                local_vmax += pad
            if vmin is None:
                vmin = local_vmin
            if vmax is None:
                vmax = local_vmax

    x_nm = nodes[:, 0] * 1.0e9
    y_nm = nodes[:, 1] * 1.0e9
    tri = mtri.Triangulation(
        x_nm,
        y_nm,
        np.asarray(mesh.triangles, dtype=np.int64),
    )

    n_snap = int(arr.shape[0])
    ncols = max(1, int(ncols))
    nrows = int(np.ceil(n_snap / ncols))

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(3.7 * ncols, 2.7 * nrows),
        squeeze=False,
    )

    last_im = None
    for k in range(nrows * ncols):
        ax = axes.flat[k]

        if k >= n_snap:
            ax.axis("off")
            continue

        last_im = ax.tripcolor(
            tri,
            arr[k],
            shading="gouraud",
            vmin=vmin,
            vmax=vmax,
        )
        ax.set_title(f"t = {t_s[k] / 1.0e-12:.4g} ps")
        ax.set_xlabel("x [nm]")
        ax.set_ylabel("y [nm]")
        ax.set_aspect("equal", adjustable="box")
        ax.grid(False)

    if last_im is not None:
        cbar = fig.colorbar(last_im, ax=axes.ravel().tolist())
        cbar.set_label(label)

    fig.suptitle(title)
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
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
    """Common triangular node-scalar plot."""
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    nodes = np.asarray(mesh.nodes, dtype=float)
    x_nm = nodes[:, 0] * 1.0e9
    y_nm = nodes[:, 1] * 1.0e9
    tri = mtri.Triangulation(
        x_nm,
        y_nm,
        np.asarray(mesh.triangles, dtype=np.int64),
    )

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


def _plot_normalized_history_curve(
    ax,
    t_ps: np.ndarray,
    history: dict,
    key: str,
    label: str,
    *,
    unit: str,
    absolute: bool,
) -> None:
    if key not in history:
        return

    y = np.asarray(history[key], dtype=float).reshape(-1)
    t, y = _match_time_and_series(t_ps, y)

    if t.size == 0:
        return

    if absolute:
        y_real = np.abs(y)
    else:
        y_real = y

    scale = _normalization_scale(y_real)
    y_norm = y_real / scale
    y_norm = np.clip(y_norm, 0.0, 1.0)

    ax.plot(t, y_norm, label=f"{label} / {_format_scale(scale, unit)}")


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
        real = _require_history_array(
            history,
            ("psi_snapshot_real_J", "psi_real_snapshot_J", "snapshot_psi_real_J"),
        )
        imag = _require_history_array(
            history,
            ("psi_snapshot_imag_J", "psi_imag_snapshot_J", "snapshot_psi_imag_J"),
        )
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

    raise KeyError(
        f"No compatible snapshot time axis found. Tried {names}; expected {n_snap} entries."
    )


def _unwrap_phase_safe(
    psi: np.ndarray,
    edges: np.ndarray,
    *,
    seed_index: int,
) -> np.ndarray:
    try:
        return unwrap_phase_graph(
            np.asarray(psi, dtype=np.complex128),
            np.asarray(edges, dtype=np.int64),
            seed_index=int(seed_index),
            subtract_mean=False,
        )
    except TypeError:
        return unwrap_phase_graph(
            np.asarray(psi, dtype=np.complex128),
            np.asarray(edges, dtype=np.int64),
        )


def _center_node_index(mesh) -> int:
    """Return node closest to the geometric center."""
    nodes = np.asarray(mesh.nodes, dtype=float)
    center = np.array(
        [
            0.5 * (float(np.min(nodes[:, 0])) + float(np.max(nodes[:, 0]))),
            0.5 * (float(np.min(nodes[:, 1])) + float(np.max(nodes[:, 1]))),
        ],
        dtype=float,
    )
    dist2 = np.sum((nodes[:, :2] - center[None, :]) ** 2, axis=1)
    return int(np.argmin(dist2))


def _mesh_edges_from_triangles(mesh) -> np.ndarray:
    """Build unique undirected edges from mesh triangles for phase unwrapping."""
    tri = np.asarray(mesh.triangles, dtype=np.int64)
    edges = np.vstack(
        [
            tri[:, [0, 1]],
            tri[:, [1, 2]],
            tri[:, [2, 0]],
        ]
    )
    edges = np.sort(edges, axis=1)
    edges = np.unique(edges, axis=0)
    return edges