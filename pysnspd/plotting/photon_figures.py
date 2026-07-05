"""Plotting helpers for pipeline 03 photon/circuit transient runs."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.tri as mtri

EV_J = 1.602176634e-19


def make_photon_run_figures(
    *,
    history: Mapping[str, np.ndarray],
    summary: Mapping[str, Any] | None,
    output_dir: str | Path,
    dpi: int = 480,
    mesh: Any | None = None,
    snapshots: Mapping[str, np.ndarray] | None = None,
    scalar_times_ps: Sequence[float] | None = None,
    center_width_nm: float = 100.0,
) -> dict[str, Path]:
    """Create photon/circuit transient figures."""

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    saved: dict[str, Path] = {}
    saved["photon_circuit_response"] = plot_photon_circuit_response(
        history=history,
        summary=summary or {},
        output_path=out / "photon_circuit_response.png",
        dpi=dpi,
    )

    if scalar_times_ps is not None and len(list(scalar_times_ps)) > 0:
        if mesh is None:
            raise ValueError("mesh is required to plot photon center scalar snapshots.")
        if snapshots is None:
            raise ValueError("transient snapshots are required to plot photon center scalar snapshots.")
        saved["photon_center_scalar_snapshots"] = plot_photon_center_scalar_snapshot_rows(
            mesh=mesh,
            snapshots=snapshots,
            summary=summary or {},
            requested_times_ps=list(scalar_times_ps),
            output_path=out / "photon_center_scalar_snapshots.png",
            center_width_nm=center_width_nm,
            dpi=dpi,
        )

    return saved


def plot_photon_circuit_response(
    *,
    history: Mapping[str, np.ndarray],
    summary: Mapping[str, Any] | None,
    output_path: str | Path,
    dpi: int = 480,
) -> Path:
    """Plot coupled-circuit response for a pipeline 03 transient."""

    output = _prepare_output(output_path)

    t_ps = _history_time_ps(history)
    I_s_A = _history_array(history, "I_s_A", required=True)
    I_rf_A = _history_array(history, "I_rf_A", required=False)
    V_tdgl_V = _history_array(history, "V_tdgl_center_V", required=True)
    V_out_V = _history_array(history, "V_out_V", required=False)

    if I_rf_A.size == 0:
        I_b_A = _history_array(history, "I_b_A", required=True)
        I_rf_A = I_b_A - I_s_A

    if V_out_V.size == 0:
        R_load = _summary_float(summary, ("circuit", "params", "R_load_ohm"), default=50.0)
        V_out_V = R_load * I_rf_A

    n = min(t_ps.size, I_s_A.size, I_rf_A.size, V_tdgl_V.size, V_out_V.size)
    if n <= 1:
        raise ValueError("transient_history.npz does not contain enough circuit samples to plot.")

    t_ps = t_ps[:n]
    I_s_A = I_s_A[:n]
    I_rf_A = I_rf_A[:n]
    V_tdgl_V = V_tdgl_V[:n]
    V_out_V = V_out_V[:n]

    dI_s_nA = (I_s_A - float(I_s_A[0])) * 1.0e9
    I_rf_nA = I_rf_A * 1.0e9
    V_tdgl_uV = V_tdgl_V * 1.0e6
    V_out_uV = V_out_V * 1.0e6

    dIs_color = "#ffa200"
    Irf_color = "#ff0000"
    Vtdgl_color = "#7700ff"
    Vout_color = "#00c3ff"
    hot_axis_color = "#ff5100"
    cold_axis_color = "#0011ff"

    photon_time_ps = _summary_float(summary, ("photon", "time_s"), default=np.nan)
    if np.isfinite(photon_time_ps):
        photon_time_ps *= 1.0e12
    else:
        photon_time_ps = _event_time_from_history(history)

    fig, ax = plt.subplots(1, 1, figsize=(8.4, 4.35), constrained_layout=False)
    fig.subplots_adjust(left=0.120, right=0.870, bottom=0.170, top=0.900)
    ax_r = ax.twinx()

    if np.isfinite(photon_time_ps):
        t0 = float(np.nanmin(t_ps))
        if photon_time_ps > t0:
            ax.axvspan(t0, photon_time_ps, color="0.88", alpha=0.62, zorder=0)
        ax.axvline(photon_time_ps, color="0.20", linestyle=":", linewidth=1.2, alpha=0.85, zorder=2)

    h_dIs, = ax.plot(t_ps, dI_s_nA, color=dIs_color, linewidth=1.9, label=r"$\Delta I_s$")
    h_Irf, = ax.plot(t_ps, I_rf_nA, color=Irf_color, linewidth=1.9, label=r"$I_{\rm RF}$")
    h_Vtdgl, = ax_r.plot(
        t_ps,
        V_tdgl_uV,
        color=Vtdgl_color,
        linewidth=1.9,
        label=r"$V_{\rm TDGL}^{center}$",
    )
    h_Vout, = ax_r.plot(t_ps, V_out_uV, color=Vout_color, linewidth=1.9, label=r"$V_{\rm out}$")

    ax.set_xlabel("t [ps]", fontsize=14)
    ax.set_ylabel(r"$\Delta I_s,\ I_{\rm RF}$ [nA]", color=hot_axis_color, fontsize=14)
    ax_r.set_ylabel(r"$V_{\rm TDGL}^{center},\ V_{\rm out}$ [$\mu$V]", color=cold_axis_color, fontsize=14)

    ax.tick_params(axis="x", labelsize=14)
    ax.tick_params(axis="y", colors=hot_axis_color, labelsize=14)
    ax_r.tick_params(axis="y", colors=cold_axis_color, labelsize=14)

    ax.grid(False)
    ax_r.grid(False)
    _clean_twin_axis(ax, ax_r)

    _center_linear_axis(ax, [dI_s_nA, I_rf_nA], frac=0.14)
    _center_linear_axis(ax_r, [V_tdgl_uV, V_out_uV], frac=0.14)

    ax.legend(
        [h_dIs, h_Irf, h_Vtdgl, h_Vout],
        [r"$\Delta I_s$", r"$I_{\rm RF}$", r"$V_{\rm TDGL}^{center}$", r"$V_{\rm out}$"],
        loc="lower right",
        frameon=False,
        ncol=2,
        columnspacing=1.1,
        handlelength=2.5,
        borderaxespad=0.3,
        fontsize=14,
    )

    fig.savefig(output, dpi=dpi, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    return output


def plot_photon_center_scalar_snapshot_rows(
    *,
    mesh: Any,
    snapshots: Mapping[str, np.ndarray],
    summary: Mapping[str, Any] | None,
    requested_times_ps: Sequence[float],
    output_path: str | Path,
    center_width_nm: float = 100.0,
    dpi: int = 480,
) -> Path:
    """Plot center-strip scalar maps for the snapshots nearest requested times."""

    output = _prepare_output(output_path)
    tri = _triangulation(mesh)
    x_nm = np.asarray(tri.x, dtype=float)
    y_nm = np.asarray(tri.y, dtype=float)
    n_nodes = x_nm.size
    xlim, ylim, crop_mask = _center_window(x_nm, y_nm, center_width_nm=center_width_nm)

    t_snap_ps = _snapshot_times_ps(snapshots)
    if t_snap_ps.size == 0:
        raise ValueError("transient_snapshots.npz lacks snapshot_t_ps/snapshot_t_s.")
    requested = np.asarray(list(requested_times_ps), dtype=float).reshape(-1)
    requested = requested[np.isfinite(requested)]
    if requested.size == 0:
        raise ValueError("No finite scalar snapshot times were requested.")

    indices = [int(np.nanargmin(np.abs(t_snap_ps - t_req))) for t_req in requested]
    actual_t = np.asarray([float(t_snap_ps[idx]) for idx in indices], dtype=float)

    fields_per_row: list[list[dict[str, Any]]] = []
    for idx in indices:
        row = _scalar_fields_for_snapshot(
            tri=tri,
            snapshots=snapshots,
            snapshot_index=idx,
            n_nodes=n_nodes,
            x_nm=x_nm,
            y_nm=y_nm,
            summary=summary or {},
        )
        fields_per_row.append(row)

    column_specs = [
        {"label": r"$|\Delta|/\Delta_0$", "force_vmin": 0.0, "force_vmax": 1.0, "positive_floor": True},
        {"label": r"$\phi$ [$\mu$V]", "force_vmin": None, "force_vmax": None, "positive_floor": False},
        {"label": r"$|q|$ [m$^{-1}$]", "force_vmin": 0.0, "force_vmax": None, "positive_floor": True},
        {"label": r"$T_e$ [K]", "force_vmin": None, "force_vmax": None, "positive_floor": False},
        {"label": r"$T_{ph}$ [K]", "force_vmin": None, "force_vmax": None, "positive_floor": False},
    ]

    limits: list[tuple[float, float]] = []
    for col, spec in enumerate(column_specs):
        values = [row[col]["values"] for row in fields_per_row]
        limits.append(
            _center_limits_for_many(
                values,
                crop_mask=crop_mask,
                force_vmin=spec["force_vmin"],
                force_vmax=spec["force_vmax"],
                positive_floor=bool(spec["positive_floor"]),
                n_nodes=n_nodes,
            )
        )

    n_rows = len(indices)
    fig_height = max(3.15, 0.50 + 2.70 * n_rows)
    fig = plt.figure(figsize=(18.4, fig_height), constrained_layout=False)
    gs = fig.add_gridspec(
        n_rows + 1,
        5,
        height_ratios=[0.11] + [1.0] * n_rows,
        wspace=0.24,
        hspace=0.18,
    )
    caxes = [fig.add_subplot(gs[0, k]) for k in range(5)]
    axes = [[fig.add_subplot(gs[r + 1, k]) for k in range(5)] for r in range(n_rows)]

    for r in range(n_rows):
        for col in range(5):
            ax = axes[r][col]
            values = np.asarray(fields_per_row[r][col]["values"], dtype=float).reshape(-1)
            if values.size != n_nodes:
                values = np.resize(values, n_nodes)

            vmin, vmax = limits[col]
            mappable = ax.tripcolor(
                tri,
                values,
                shading="gouraud",
                vmin=vmin,
                vmax=vmax,
            )
            if r == 0:
                cbar = fig.colorbar(mappable, cax=caxes[col], orientation="horizontal")
                cbar.set_label(str(column_specs[col]["label"]), labelpad=2.0)
                cbar.ax.xaxis.set_ticks_position("top")
                cbar.ax.xaxis.set_label_position("top")
                cbar.ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=3))

            ax.set_xlim(*xlim)
            ax.set_ylim(*ylim)
            ax.set_aspect("equal", adjustable="box")
            ax.grid(False)

            if r == n_rows - 1:
                ax.set_xlabel("x [nm]", fontsize=12)
            else:
                ax.set_xlabel("")
                ax.tick_params(axis="x", labelbottom=False)

            if col == 0:
                ax.set_ylabel("y [nm]", fontsize=12)
                ax.text(
                    -0.26,
                    0.50,
                    rf"$t={actual_t[r]:.3g}$ ps",
                    transform=ax.transAxes,
                    ha="right",
                    va="center",
                    rotation=90,
                    fontsize=13,
                    color="black",
                    clip_on=False,
                )
            else:
                ax.set_ylabel("")
                ax.tick_params(axis="y", labelleft=False)

    fig.savefig(output, dpi=dpi, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    return output


def _scalar_fields_for_snapshot(
    *,
    tri: mtri.Triangulation,
    snapshots: Mapping[str, np.ndarray],
    snapshot_index: int,
    n_nodes: int,
    x_nm: np.ndarray,
    y_nm: np.ndarray,
    summary: Mapping[str, Any],
) -> list[dict[str, np.ndarray]]:
    psi = _snapshot_complex_psi(snapshots, snapshot_index, n_nodes=n_nodes)
    delta0_J = _resolve_delta0_J(summary, snapshots, psi, n_nodes=n_nodes)
    delta_norm = np.abs(psi) / max(delta0_J, 1.0e-300)
    delta_norm[~np.isfinite(delta_norm)] = 0.0

    phi_uV = 1.0e6 * _snapshot_node_array(snapshots, ("phi_snapshot_V", "phi_V"), snapshot_index, n_nodes=n_nodes)
    q_abs = _phase_gradient_q_abs_m_inv(tri, psi, x_nm=x_nm, y_nm=y_nm)
    Te_K = _snapshot_node_array(snapshots, ("Te_snapshot_K", "Te_K"), snapshot_index, n_nodes=n_nodes)
    Tph_K = _snapshot_node_array(snapshots, ("Tph_snapshot_K", "Tph_K"), snapshot_index, n_nodes=n_nodes)

    return [
        {"values": delta_norm},
        {"values": phi_uV},
        {"values": q_abs},
        {"values": Te_K},
        {"values": Tph_K},
    ]


def _snapshot_complex_psi(
    snapshots: Mapping[str, np.ndarray],
    snapshot_index: int,
    *,
    n_nodes: int,
) -> np.ndarray:
    if "psi_real_snapshot_J" in snapshots and "psi_imag_snapshot_J" in snapshots:
        real = _snapshot_node_array(snapshots, ("psi_real_snapshot_J",), snapshot_index, n_nodes=n_nodes)
        imag = _snapshot_node_array(snapshots, ("psi_imag_snapshot_J",), snapshot_index, n_nodes=n_nodes)
        return real + 1j * imag

    if "delta_snapshot_meV" in snapshots:
        delta_meV = _snapshot_node_array(snapshots, ("delta_snapshot_meV",), snapshot_index, n_nodes=n_nodes)
        return delta_meV * 1.0e-3 * EV_J + 0.0j

    return np.zeros(n_nodes, dtype=np.complex128)


def _snapshot_node_array(
    snapshots: Mapping[str, np.ndarray],
    keys: tuple[str, ...],
    snapshot_index: int,
    *,
    n_nodes: int,
    default: float = 0.0,
) -> np.ndarray:
    for key in keys:
        if key not in snapshots:
            continue
        arr = np.asarray(snapshots[key], dtype=float)
        if arr.ndim == 1:
            out = arr.reshape(-1)
        elif arr.ndim >= 2:
            idx = min(max(int(snapshot_index), 0), arr.shape[0] - 1)
            out = np.asarray(arr[idx], dtype=float).reshape(-1)
        else:
            continue
        if out.size:
            if out.size != n_nodes:
                out = np.resize(out, n_nodes)
            return out.astype(float, copy=False)
    return np.full(n_nodes, float(default), dtype=float)


def _resolve_delta0_J(
    summary: Mapping[str, Any],
    snapshots: Mapping[str, np.ndarray],
    psi_this: np.ndarray,
    *,
    n_nodes: int,
) -> float:
    for path in (
        ("material", "delta0_J"),
        ("gtdgl_material", "delta0_J"),
        ("initial_state", "delta0_J"),
        ("solver", "delta0_J"),
    ):
        value = _summary_float(summary, path, default=np.nan)
        if np.isfinite(value) and value > 0.0:
            return float(value)

    if "delta0_J" in snapshots:
        arr = np.asarray(snapshots["delta0_J"], dtype=float).reshape(-1)
        arr = arr[np.isfinite(arr) & (arr > 0.0)]
        if arr.size:
            return float(arr[0])

    if "psi_real_snapshot_J" in snapshots and "psi_imag_snapshot_J" in snapshots:
        real = np.asarray(snapshots["psi_real_snapshot_J"], dtype=float)
        imag = np.asarray(snapshots["psi_imag_snapshot_J"], dtype=float)
        if real.ndim >= 2 and imag.ndim >= 2 and real.shape[0] and imag.shape[0]:
            psi0 = np.asarray(real[0], dtype=float).reshape(-1) + 1j * np.asarray(imag[0], dtype=float).reshape(-1)
            if psi0.size != n_nodes:
                psi0 = np.resize(psi0, n_nodes)
            vals = np.abs(psi0)
        else:
            vals = np.abs(psi_this)
    else:
        vals = np.abs(psi_this)

    vals = np.asarray(vals, dtype=float).reshape(-1)
    vals = vals[np.isfinite(vals) & (vals > 0.0)]
    if vals.size:
        return float(np.nanpercentile(vals, 99.5))
    return 1.0


def _phase_gradient_q_abs_m_inv(
    tri: mtri.Triangulation,
    psi: np.ndarray,
    *,
    x_nm: np.ndarray,
    y_nm: np.ndarray,
) -> np.ndarray:
    psi = np.asarray(psi, dtype=np.complex128).reshape(-1)
    n_nodes = x_nm.size
    if psi.size != n_nodes:
        return np.zeros(n_nodes, dtype=float)

    theta = np.angle(psi)
    x_m = np.asarray(x_nm, dtype=float) * 1.0e-9
    y_m = np.asarray(y_nm, dtype=float) * 1.0e-9

    triangles = np.asarray(tri.triangles, dtype=np.int64)
    if triangles.size == 0:
        return np.zeros(n_nodes, dtype=float)

    edges = np.vstack((triangles[:, [0, 1]], triangles[:, [1, 2]], triangles[:, [2, 0]]))
    edges = np.sort(edges, axis=1)
    edges = np.unique(edges, axis=0)

    Axx = np.zeros(n_nodes, dtype=float)
    Axy = np.zeros(n_nodes, dtype=float)
    Ayy = np.zeros(n_nodes, dtype=float)
    bx = np.zeros(n_nodes, dtype=float)
    by = np.zeros(n_nodes, dtype=float)

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

    qx = np.zeros(n_nodes, dtype=float)
    qy = np.zeros(n_nodes, dtype=float)
    qx[good] = (Ayy[good] * bx[good] - Axy[good] * by[good]) / det[good]
    qy[good] = (-Axy[good] * bx[good] + Axx[good] * by[good]) / det[good]

    q_abs = np.sqrt(qx * qx + qy * qy)
    q_abs[~np.isfinite(q_abs)] = 0.0
    return q_abs


def _triangulation(mesh: Any) -> mtri.Triangulation:
    nodes = np.asarray(mesh.nodes, dtype=float)
    if nodes.ndim != 2 or nodes.shape[1] < 2:
        raise ValueError("mesh.nodes must have shape (n_nodes, >=2).")
    if hasattr(mesh, "triangles"):
        triangles = np.asarray(mesh.triangles, dtype=np.int64)
    elif hasattr(mesh, "elements"):
        triangles = np.asarray(mesh.elements, dtype=np.int64)
    else:
        raise AttributeError("mesh must expose triangles or elements.")
    x_nm = nodes[:, 0] * 1.0e9
    y_nm = nodes[:, 1] * 1.0e9
    return mtri.Triangulation(x_nm, y_nm, triangles)


def _center_window(
    x_nm: np.ndarray,
    y_nm: np.ndarray,
    *,
    center_width_nm: float,
) -> tuple[tuple[float, float], tuple[float, float], np.ndarray]:
    x = np.asarray(x_nm, dtype=float)
    y = np.asarray(y_nm, dtype=float)
    x_mid = 0.5 * (float(np.nanmin(x)) + float(np.nanmax(x)))
    half = 0.5 * float(center_width_nm)
    xlim = (x_mid - half, x_mid + half)
    ylim = (float(np.nanmin(y)), float(np.nanmax(y)))
    mask = (x >= xlim[0]) & (x <= xlim[1]) & np.isfinite(x) & np.isfinite(y)
    if not np.any(mask):
        mask = np.ones_like(x, dtype=bool)
    return xlim, ylim, mask


def _center_limits_for_many(
    arrays: Sequence[np.ndarray],
    *,
    crop_mask: np.ndarray,
    force_vmin: float | None,
    force_vmax: float | None,
    positive_floor: bool,
    n_nodes: int,
) -> tuple[float, float]:
    finite: list[np.ndarray] = []
    for arr in arrays:
        vals = np.asarray(arr, dtype=float).reshape(-1)
        if vals.size != n_nodes:
            vals = np.resize(vals, n_nodes)
        vis = vals[crop_mask]
        vis = vis[np.isfinite(vis)]
        if vis.size:
            finite.append(vis)

    if finite:
        all_vals = np.concatenate(finite)
        vmin = float(np.nanmin(all_vals))
        vmax = float(np.nanmax(all_vals))
    else:
        vmin, vmax = 0.0, 1.0

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


def _prepare_output(output_path: str | Path) -> Path:
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    return output


def load_npz_dict(path: str | Path) -> dict[str, np.ndarray]:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Missing NPZ file: {p}")
    with np.load(p, allow_pickle=True) as data:
        return {key: np.asarray(data[key]) for key in data.files}


def _history_array(history: Mapping[str, np.ndarray], key: str, *, required: bool) -> np.ndarray:
    if key not in history:
        if required:
            raise KeyError(f"transient_history.npz lacks required key: {key}")
        return np.empty(0, dtype=float)
    arr = np.asarray(history[key], dtype=float).reshape(-1)
    if arr.size == 0 and required:
        raise ValueError(f"transient_history.npz key {key!r} is empty.")
    return arr


def _history_time_ps(history: Mapping[str, np.ndarray]) -> np.ndarray:
    if "t_ps" in history:
        t = np.asarray(history["t_ps"], dtype=float).reshape(-1)
        if t.size:
            return t
    if "t_s" in history:
        t = np.asarray(history["t_s"], dtype=float).reshape(-1)
        if t.size:
            return t / 1.0e-12
    raise KeyError("transient_history.npz lacks t_ps/t_s.")


def _snapshot_times_ps(snapshots: Mapping[str, np.ndarray]) -> np.ndarray:
    if "snapshot_t_ps" in snapshots:
        t = np.asarray(snapshots["snapshot_t_ps"], dtype=float).reshape(-1)
        if t.size:
            return t
    if "snapshot_t_s" in snapshots:
        t = np.asarray(snapshots["snapshot_t_s"], dtype=float).reshape(-1)
        if t.size:
            return t / 1.0e-12
    raise KeyError("transient_snapshots.npz lacks snapshot_t_ps/snapshot_t_s.")


def _summary_float(summary: Mapping[str, Any] | None, path: tuple[str, ...], *, default: float) -> float:
    cur: Any = summary or {}
    for key in path:
        if isinstance(cur, Mapping) and key in cur:
            cur = cur[key]
        else:
            return float(default)
    try:
        return float(cur)
    except Exception:
        return float(default)


def _event_time_from_history(history: Mapping[str, np.ndarray]) -> float:
    if "photon_applied" not in history:
        return float("nan")
    applied = np.asarray(history["photon_applied"]).reshape(-1).astype(bool)
    if not np.any(applied):
        return float("nan")
    t = _history_time_ps(history)
    n = min(t.size, applied.size)
    idx = np.flatnonzero(applied[:n])
    if idx.size == 0:
        return float("nan")
    return float(t[idx[0]])


def _clean_twin_axis(ax, ax_r) -> None:
    ax.spines["right"].set_visible(False)
    ax_r.spines["left"].set_visible(False)
    ax_r.patch.set_alpha(0.0)


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
