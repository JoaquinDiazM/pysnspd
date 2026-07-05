"""Plotting helpers for pipeline 03 photon/circuit transient runs."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

import numpy as np

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt


def make_photon_run_figures(
    *,
    history: Mapping[str, np.ndarray],
    summary: Mapping[str, Any] | None,
    output_dir: str | Path,
    dpi: int = 480,
) -> dict[str, Path]:
    """Create the first set of photon/circuit transient figures."""

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    saved: dict[str, Path] = {}
    saved["photon_circuit_response"] = plot_photon_circuit_response(
        history=history,
        summary=summary or {},
        output_path=out / "photon_circuit_response.png",
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
    """Plot coupled-circuit response for a pipeline 03 transient.

    Left y-axis:
        Delta I_s = I_s(t) - I_s(0), and I_RF = I_b - I_s, in nA.

    Right y-axis:
        V_TDGL^center and V_out, in microvolts.
    """

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


def _history_array(
    history: Mapping[str, np.ndarray],
    key: str,
    *,
    required: bool,
) -> np.ndarray:
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


def _summary_float(
    summary: Mapping[str, Any] | None,
    path: tuple[str, ...],
    *,
    default: float,
) -> float:
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
