"""Presentation-quality Usadel DOS curve plots for PRE/E-type pipelines.

The PRE catalogue stores the real-axis dirty-limit Usadel DOS over independent
``(|Delta|, Gamma_q, E)`` axes. This module turns those catalogue slices into
line plots that are more useful for thesis figures than the raw 2-D colormaps.

Two slices are provided:

1. ``Delta = Delta_eq(q)``: the equilibrium branch used by the existing
   ``usadel_equilibrium_dos_map`` diagnostic.
2. ``Delta = Delta_0``: a fixed-gap section through the same catalogue, useful
   to isolate the spectral depairing caused by ``Gamma_q`` from the additional
   self-consistent depression of the order parameter.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Sequence

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np

from pysnspd.plotting.style import THESIS_DOUBLE_FIGSIZE, apply_thesis_style

MEV_J = 1.602176634e-22
_DEFAULT_CURRENT_FRACTIONS: tuple[float, ...] = (0.0, 0.25, 0.50, 0.65, 0.80, 0.95)


def plot_usadel_dos_curves_equilibrium_gap(
    usadel_catalog: Any,
    output_path: str | Path,
    *,
    current_fractions: Sequence[float] = _DEFAULT_CURRENT_FRACTIONS,
    dpi: int = 480,
    energy_max_meV: float | None = 4.0,
    energy_window: bool = True,
) -> Path:
    """Plot DOS curves ``rho(E; Delta_eq(q), q)`` for selected currents.

    The selected currents are interpreted as fractions of the model critical
    current obtained from the PRE Matsubara calibration branch. Only the stable
    pre-maximum branch of ``I_s(q)`` is used to convert ``I_s/I_c`` to ``q``.
    """
    return _plot_dos_curves(
        usadel_catalog,
        output_path,
        mode="equilibrium_gap",
        current_fractions=current_fractions,
        dpi=dpi,
        energy_max_meV=energy_max_meV,
        energy_window=energy_window,
    )


def plot_usadel_dos_curves_fixed_delta0(
    usadel_catalog: Any,
    output_path: str | Path,
    *,
    current_fractions: Sequence[float] = _DEFAULT_CURRENT_FRACTIONS,
    dpi: int = 480,
    energy_max_meV: float | None = 4.0,
    energy_window: bool = True,
) -> Path:
    """Plot DOS curves ``rho(E; Delta_BCS(0), q)`` for selected currents.

    This is not a self-consistent finite-current branch. It is a fixed-gap
    catalogue section that keeps ``Delta`` at ``Delta_BCS(0)`` while changing the
    depairing parameter through ``q``.
    """
    return _plot_dos_curves(
        usadel_catalog,
        output_path,
        mode="fixed_delta0",
        current_fractions=current_fractions,
        dpi=dpi,
        energy_max_meV=energy_max_meV,
        energy_window=energy_window,
    )


def _plot_dos_curves(
    usadel_catalog: Any,
    output_path: str | Path,
    *,
    mode: str,
    current_fractions: Sequence[float],
    dpi: int,
    energy_max_meV: float | None,
    energy_window: bool,
) -> Path:
    apply_thesis_style()
    output = _prepare_output(output_path)
    fractions = _normalize_current_fractions(current_fractions)
    energy_meV_full = _energy_axis_meV(usadel_catalog)
    positive_energy = np.isfinite(energy_meV_full) & (energy_meV_full >= 0.0)
    if not np.any(positive_energy):
        raise ValueError("energy_values_J has no finite non-negative energy points.")

    curves = _extract_dos_curves_for_current_fractions(
        usadel_catalog,
        current_fractions=fractions,
        mode=mode,
    )

    energy_meV = energy_meV_full[positive_energy]
    for curve in curves:
        curve["rho"] = np.asarray(curve["rho"], dtype=float)[positive_energy]

    metadata = getattr(usadel_catalog, "metadata", {})
    T_bias_K = _metadata_float(metadata, "T_bias_K")
    Tc_K = _metadata_float(metadata, "Tc_K")
    T_ratio = T_bias_K / Tc_K if np.isfinite(T_bias_K) and np.isfinite(Tc_K) and Tc_K > 0.0 else np.nan

    if mode == "equilibrium_gap":
        state_label = r"$\Delta=\Delta_{\mathrm{eq}}(q)$"
    elif mode == "fixed_delta0":
        state_label = r"$\Delta=\Delta_{\mathrm{BCS}}(0)$"
    else:  # pragma: no cover - internal defensive branch
        raise ValueError(f"Unknown DOS-curve mode: {mode!r}")

    fig, ax = plt.subplots(figsize=THESIS_DOUBLE_FIGSIZE)
    label_fs = 12
    tick_fs = 10
    legend_fs = 9.0
    legend_title_fs = 9.0

    colors = _curve_colors(len(curves))
    handles: list[Any] = []
    labels: list[str] = []

    for curve, color in zip(curves, colors):
        line, = ax.plot(
            energy_meV,
            curve["rho"],
            linewidth=1.25,
            color=color,
            label=rf"{_fraction_label(curve['fraction'])}",
        )
        handles.append(line)
        labels.append(rf"{_fraction_label(curve['fraction'])}")

    rho_stack = np.vstack([np.asarray(curve["rho"], dtype=float) for curve in curves])
    xlim = _energy_xlim(
        energy_meV,
        rho_stack,
        energy_max_meV=energy_max_meV,
        energy_window=energy_window,
    )
    ax.set_xlim(*xlim)

    visible = (energy_meV >= xlim[0]) & (energy_meV <= xlim[1])
    y_values = rho_stack[:, visible] if np.any(visible) else rho_stack
    y_values = y_values[np.isfinite(y_values)]
    y_top = 1.05 * float(np.nanmax(y_values)) if y_values.size else 1.0
    ax.set_ylim(bottom=0.0, top=max(1.05, y_top))

    ax.set_xlabel(r"Energy $E$ [meV]", fontsize=label_fs)
    ax.set_ylabel(r"$\rho(E)$", fontsize=label_fs)
    ax.tick_params(axis="both", which="major", labelsize=tick_fs, direction="in", length=7.0, width=1.15)
    ax.tick_params(axis="both", which="minor", direction="in", length=4.0, width=0.9)
    ax.minorticks_on()
    ax.grid(True, which="major", linewidth=0.50, alpha=0.20)
    ax.grid(True, which="minor", linewidth=0.30, alpha=0.08)

    if np.isfinite(T_ratio):
        temp_line = rf"$T/T_c={T_ratio:.3f}$"
    elif np.isfinite(T_bias_K):
        temp_line = rf"$T={T_bias_K:.2f}$ [K]"
    else:
        temp_line = ""

    legend_title = state_label + "\n" + r"$I_s/I_c$"
    if temp_line:
        legend_title = state_label + "\n" + temp_line + "\n" + r"$I_s/I_c$"

    legend = ax.legend(
        handles,
        labels,
        title=legend_title,
        fontsize=legend_fs,
        title_fontsize=legend_title_fs,
        loc="upper right",
        frameon=True,
        fancybox=False,
        framealpha=1.0,
        facecolor="white",
        edgecolor="0.35",
        handlelength=3.0,
        borderpad=0.75,
        labelspacing=0.52,
    )
    legend.get_frame().set_linewidth(1.05)
    legend.set_zorder(10)
    legend.get_title().set_multialignment("center")
    legend._legend_box.align = "center"

    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output


def _extract_dos_curves_for_current_fractions(
    usadel_catalog: Any,
    *,
    current_fractions: Sequence[float],
    mode: str,
) -> list[dict[str, Any]]:
    q_targets = _q_targets_from_current_fractions(usadel_catalog, current_fractions)
    q_axis = np.asarray(usadel_catalog.q_values_m_inv, dtype=float)
    delta_axis = np.asarray(usadel_catalog.delta_values_J, dtype=float)
    rho = np.asarray(usadel_catalog.rho_delta_gamma_E, dtype=float)

    if rho.ndim != 3:
        raise ValueError(f"rho_delta_gamma_E must be 3-D, got shape {rho.shape}.")
    if rho.shape[0] != delta_axis.size:
        raise ValueError("rho_delta_gamma_E first axis does not match delta_values_J.")
    if rho.shape[1] != q_axis.size:
        raise ValueError("rho_delta_gamma_E second axis does not match q_values_m_inv.")

    out: list[dict[str, Any]] = []
    delta0_J = _resolve_delta0_J(usadel_catalog)
    for fraction, q_target in zip(current_fractions, q_targets):
        iq = _nearest_index(q_axis, q_target)
        q_actual = float(q_axis[iq])

        if mode == "equilibrium_gap":
            delta_target = _delta_eq_at_q(usadel_catalog, q_actual)
        elif mode == "fixed_delta0":
            delta_target = delta0_J
        else:  # pragma: no cover - internal defensive branch
            raise ValueError(f"Unknown DOS-curve mode: {mode!r}")

        idelta = _nearest_index(delta_axis, delta_target)
        out.append(
            {
                "fraction": float(fraction),
                "q_target_m_inv": float(q_target),
                "q_actual_m_inv": q_actual,
                "delta_target_J": float(delta_target),
                "delta_actual_J": float(delta_axis[idelta]),
                "iq": int(iq),
                "idelta": int(idelta),
                "rho": np.asarray(rho[idelta, iq, :], dtype=float),
            }
        )
    return out


def _q_targets_from_current_fractions(
    usadel_catalog: Any,
    fractions: Sequence[float],
) -> np.ndarray:
    q = np.asarray(usadel_catalog.calibration_q_values_m_inv, dtype=float)
    current = np.asarray(usadel_catalog.calibration_current_values_A, dtype=float)
    finite = np.isfinite(q) & np.isfinite(current)
    q = q[finite]
    current = current[finite]
    if q.size < 2:
        raise ValueError("Usadel calibration table needs at least two finite points.")

    order_q = np.argsort(q)
    q = q[order_q]
    current = current[order_q]

    i_peak = int(np.nanargmax(current))
    if i_peak < 1:
        raise ValueError("Could not identify a usable pre-critical current branch.")

    q_branch = q[: i_peak + 1]
    current_branch = current[: i_peak + 1]
    Ic = float(np.nanmax(current_branch))
    if not np.isfinite(Ic) or Ic <= 0.0:
        raise ValueError(f"Invalid model critical current from calibration branch: {Ic!r}")

    current_mono = np.maximum.accumulate(current_branch)
    unique_current, unique_idx = np.unique(current_mono, return_index=True)
    q_unique = q_branch[unique_idx]
    if unique_current.size < 2:
        raise ValueError("Current branch is not usable for I/Ic interpolation.")

    targets = np.clip(np.asarray(fractions, dtype=float), 0.0, 1.0) * Ic
    return np.interp(targets, unique_current, q_unique)


def _delta_eq_at_q(usadel_catalog: Any, q_m_inv: float) -> float:
    q_cal = np.asarray(usadel_catalog.calibration_q_values_m_inv, dtype=float)
    delta_cal = np.asarray(usadel_catalog.calibration_delta_eq_values_J, dtype=float)
    finite = np.isfinite(q_cal) & np.isfinite(delta_cal)
    q_cal = q_cal[finite]
    delta_cal = delta_cal[finite]
    if q_cal.size == 0:
        raise ValueError("calibration_delta_eq_values_J has no finite points.")

    order = np.argsort(q_cal)
    q_cal = q_cal[order]
    delta_cal = delta_cal[order]
    return float(np.interp(float(q_m_inv), q_cal, delta_cal, left=delta_cal[0], right=delta_cal[-1]))


def _resolve_delta0_J(usadel_catalog: Any) -> float:
    metadata = getattr(usadel_catalog, "metadata", {})
    for value in _metadata_values_recursive(metadata, ("delta0_J", "Delta0_J", "delta_ref_J", "Delta_ref_J")):
        try:
            out = float(value)
        except Exception:
            continue
        if np.isfinite(out) and out > 0.0:
            return out

    delta_axis = np.asarray(usadel_catalog.delta_values_J, dtype=float)
    finite = delta_axis[np.isfinite(delta_axis)]
    if finite.size == 0:
        raise ValueError("delta_values_J has no finite values, cannot resolve Delta0.")
    return float(np.nanmax(finite))


def _metadata_values_recursive(metadata: Any, keys: tuple[str, ...]) -> list[Any]:
    if not isinstance(metadata, dict):
        return []
    found: list[Any] = []
    for key, value in metadata.items():
        if key in keys:
            found.append(value)
        if isinstance(value, dict):
            found.extend(_metadata_values_recursive(value, keys))
    return found


def _energy_axis_meV(usadel_catalog: Any) -> np.ndarray:
    energy_J = np.asarray(usadel_catalog.energy_values_J, dtype=float)
    energy_meV = energy_J / MEV_J
    if energy_meV.ndim != 1 or energy_meV.size == 0:
        raise ValueError("energy_values_J must be a non-empty 1-D axis.")
    return energy_meV


def _energy_xlim(
    energy_meV: np.ndarray,
    curves: np.ndarray,
    *,
    energy_max_meV: float | None,
    energy_window: bool,
) -> tuple[float, float]:
    x_min = max(0.0, float(np.nanmin(energy_meV)))
    x_max_available = float(np.nanmax(energy_meV))

    # E-type DOS figures are thesis figures, not exploratory colormaps.
    # Keep the default positive-energy domain fixed to 4 meV unless the
    # pipeline explicitly requests another value.
    if energy_max_meV is None:
        energy_max_meV = 4.0

    x_max = min(float(energy_max_meV), x_max_available)
    if np.isfinite(x_max) and x_max > x_min:
        return (x_min, x_max)

    if not energy_window:
        return (x_min, x_max_available)

    return _auto_energy_xlim_from_tail(
        energy_meV,
        curves,
        rel_threshold=0.010,
        tail_fraction=0.15,
        pad_fraction=0.12,
        min_visible_fraction=0.35,
    )


def _auto_energy_xlim_from_tail(
    E_meV: np.ndarray,
    values_qE: np.ndarray,
    *,
    rel_threshold: float = 0.015,
    tail_fraction: float = 0.15,
    pad_fraction: float = 0.10,
    min_visible_fraction: float = 0.20,
) -> tuple[float, float]:
    """Choose a compact energy window from deviations from the high-energy tail."""
    E = np.asarray(E_meV, dtype=float)
    values = np.asarray(values_qE, dtype=float)

    if E.ndim != 1 or E.size < 3 or values.size == 0:
        return (float(np.nanmin(E)), float(np.nanmax(E)))
    if values.ndim == 1:
        values = values[None, :]
    if values.shape[-1] != E.size:
        return (float(np.nanmin(E)), float(np.nanmax(E)))

    finite = np.isfinite(E)
    if not np.any(finite):
        return (0.0, 1.0)

    n_tail = max(3, int(round(tail_fraction * E.size)))
    n_tail = min(n_tail, E.size)
    tail = values[..., -n_tail:]
    baseline = np.nanmedian(tail, axis=-1, keepdims=True)
    scale = np.maximum(np.nanmax(np.abs(values - baseline), axis=-1, keepdims=True), 1.0)
    deviation = np.nanmax(np.abs(values - baseline) / scale, axis=0)

    active = np.isfinite(deviation) & (deviation > float(rel_threshold))
    if not np.any(active):
        return (max(0.0, float(np.nanmin(E))), float(np.nanmax(E)))

    last = int(np.max(np.flatnonzero(active)))
    i_max = min(E.size - 1, last + max(1, int(round(pad_fraction * E.size))))
    min_i_max = max(1, int(round(min_visible_fraction * (E.size - 1))))
    i_max = max(i_max, min_i_max)

    x0 = max(0.0, float(np.nanmin(E)))
    x1 = float(E[i_max])
    if not np.isfinite(x1) or x1 <= x0:
        x1 = float(np.nanmax(E))
    return (x0, x1)


def _curve_colors(n: int) -> list[Any]:
    if n <= 1:
        return ["tab:blue"]
    cmap = plt.get_cmap("plasma")
    return [cmap(v) for v in np.linspace(0.10, 0.82, n)]


def _fraction_label(value: float) -> str:
    v = float(value)
    if abs(v) < 5.0e-13:
        return "0"
    return f"{v:.2f}"


def _normalize_current_fractions(fractions: Sequence[float]) -> tuple[float, ...]:
    out = tuple(float(v) for v in fractions)
    if len(out) == 0:
        raise ValueError("At least one current fraction must be provided.")
    for value in out:
        if not np.isfinite(value):
            raise ValueError(f"Current fractions must be finite, got {value!r}.")
        if value < 0.0 or value > 1.0:
            raise ValueError(f"Current fractions must lie in [0, 1], got {value!r}.")
    return out


def _nearest_index(axis: np.ndarray, value: float) -> int:
    arr = np.asarray(axis, dtype=float)
    if arr.size == 0:
        raise ValueError("Cannot find nearest index on an empty axis.")
    return int(np.nanargmin(np.abs(arr - float(value))))


def _metadata_float(metadata: Any, key: str) -> float:
    if not isinstance(metadata, dict) or key not in metadata:
        return float("nan")
    try:
        return float(metadata[key])
    except Exception:
        return float("nan")


def _prepare_output(output_path: str | Path) -> Path:
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    return output


__all__ = [
    "plot_usadel_dos_curves_equilibrium_gap",
    "plot_usadel_dos_curves_fixed_delta0",
]
