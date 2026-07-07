"""Extra plotting utilities for Usadel equilibrium-gap diagnostics.

The official PRE-run does *not* store a two-dimensional Delta_eq(T, q)
array.  It stores the DOS catalogue plus calibration metadata at T_bias and
then appends the strict local-current table js_A_m2[Te, delta, q].  Therefore
this module reconstructs Delta_eq(T, q) for plotting by reusing the same
Matsubara self-consistency solver used by the calibration layer.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from pysnspd.usadel.calibration import matsubara_energy_axis_J, solve_gap_for_gamma_J
from pysnspd.usadel.parameters import E_CHARGE_C, HBAR_J_S
from scipy.constants import Boltzmann, e

_MEV_J = E_CHARGE_C * 1.0e-3


@dataclass(frozen=True)
class UsadelGapCatalog:
    """Minimal representation needed for E1 pre-run gap plots.

    Parameters
    ----------
    temperature_K:
        One-dimensional electron-temperature grid in K.
    q_m_inv:
        One-dimensional superfluid-momentum grid in m^-1.  For E1 this is the
        selected six-curve grid by default, not necessarily the full DOS axis.
    gap_meV:
        Equilibrium gap array with shape ``(n_temperature, n_q)`` in meV.
    source_key:
        Human-readable description of how the field was obtained.
    q_critical_m_inv:
        Critical q used for normalization in the legend.
    metadata:
        Compact metadata copied/derived from the PRE-run Usadel catalogue.
    """

    temperature_K: np.ndarray
    q_m_inv: np.ndarray
    gap_meV: np.ndarray
    source_key: str
    q_critical_m_inv: float
    metadata: dict[str, Any]


def _load_npz_arrays(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"Usadel catalog not found: {path}")
    with np.load(path, allow_pickle=True) as data:
        return {key: data[key] for key in data.files}


def _metadata_from_arrays(arrays: Mapping[str, Any]) -> dict[str, Any]:
    raw = arrays.get("metadata")
    if raw is None:
        return {}
    try:
        item = raw.item()
    except Exception:
        return {}
    if isinstance(item, Mapping):
        return dict(item)
    return {}


def _as_float(value: Any, *, name: str) -> float:
    try:
        out = float(np.asarray(value).squeeze())
    except Exception as exc:  # pragma: no cover - defensive error detail
        raise ValueError(f"Could not read scalar {name!r} from Usadel catalog.") from exc
    if not np.isfinite(out):
        raise ValueError(f"Scalar {name!r} is not finite: {out!r}")
    return out


def _as_int(value: Any, *, name: str) -> int:
    out = int(round(_as_float(value, name=name)))
    if out <= 0:
        raise ValueError(f"Integer {name!r} must be positive, got {out}.")
    return out


def _catalog_scalar(
    arrays: Mapping[str, Any],
    metadata: Mapping[str, Any],
    keys: tuple[str, ...],
    *,
    required: bool = True,
    default: float | None = None,
) -> float:
    for key in keys:
        if key in metadata:
            return _as_float(metadata[key], name=key)
        if key in arrays:
            return _as_float(arrays[key], name=key)
    calibration = metadata.get("calibration")
    if isinstance(calibration, Mapping):
        for key in keys:
            if key in calibration:
                return _as_float(calibration[key], name=f"calibration.{key}")
    if default is not None:
        return float(default)
    if required:
        joined = ", ".join(keys)
        raise KeyError(f"Could not find any of these scalar fields in the Usadel catalog: {joined}")
    return float("nan")


def _catalog_int(
    arrays: Mapping[str, Any],
    metadata: Mapping[str, Any],
    keys: tuple[str, ...],
    *,
    default: int,
) -> int:
    for key in keys:
        if key in metadata:
            return _as_int(metadata[key], name=key)
        if key in arrays:
            return _as_int(arrays[key], name=key)
    return int(default)


def _array_1d(arrays: Mapping[str, Any], keys: tuple[str, ...]) -> np.ndarray | None:
    for key in keys:
        if key in arrays:
            out = np.asarray(arrays[key], dtype=float).squeeze()
            if out.ndim == 1 and out.size > 0:
                return out
    return None


def _resolve_q_critical(
    arrays: Mapping[str, Any],
    metadata: Mapping[str, Any],
    q_critical_m_inv: float | None,
) -> float:
    if q_critical_m_inv is not None:
        qcrit = float(q_critical_m_inv)
        if not np.isfinite(qcrit) or qcrit <= 0.0:
            raise ValueError(f"--q-critical-m-inv must be finite and positive, got {qcrit!r}")
        return qcrit

    try:
        return _catalog_scalar(
            arrays,
            metadata,
            ("q_critical_m_inv", "qcrit_m_inv", "q_c_m_inv"),
            required=True,
        )
    except KeyError:
        pass

    q_cal = _array_1d(arrays, ("calibration_q_values_m_inv",))
    current = _array_1d(arrays, ("calibration_current_values_A", "calibration_current_density_values_A_m2"))
    if q_cal is not None and current is not None and q_cal.size == current.size:
        idx = int(np.nanargmax(np.asarray(current, dtype=float)))
        qcrit = float(q_cal[idx])
        if np.isfinite(qcrit) and qcrit > 0.0:
            return qcrit

    q_axis = _array_1d(arrays, ("q_values_m_inv", "q_axis_m_inv", "calibration_q_values_m_inv"))
    if q_axis is not None:
        qcrit = float(np.nanmax(q_axis))
        if np.isfinite(qcrit) and qcrit > 0.0:
            return qcrit

    raise KeyError("Could not determine q_c from metadata, calibration arrays, or q axis.")


def _temperature_axis(
    *,
    Tc_K: float,
    T_bias_K: float,
    n_temperature: int,
    T_min_K: float | None,
    T_max_K: float | None,
) -> np.ndarray:
    n = int(n_temperature)
    if n < 3:
        raise ValueError("n_temperature must be at least 3.")
    if not np.isfinite(Tc_K) or Tc_K <= 0.0:
        raise ValueError(f"Tc_K must be finite and positive, got {Tc_K!r}")

    lo = float(T_min_K) if T_min_K is not None else float(T_bias_K)
    if not np.isfinite(lo) or lo <= 0.0:
        lo = max(1.0e-3, 0.01 * Tc_K)

    hi = float(T_max_K) if T_max_K is not None else float(Tc_K)
    if not np.isfinite(hi) or hi <= 0.0:
        hi = float(Tc_K)
    hi = min(hi, float(Tc_K))

    if lo >= hi:
        raise ValueError(f"Invalid temperature range: T_min_K={lo:.6g}, T_max_K={hi:.6g}")
    return np.linspace(lo, hi, n)


def build_gap_eq_temperature_catalog(
    npz_path: str | Path,
    *,
    n_curves: int = 6,
    n_temperature: int = 160,
    T_min_K: float | None = None,
    T_max_K: float | None = None,
    q_critical_m_inv: float | None = None,
    n_matsubara: int | None = None,
) -> UsadelGapCatalog:
    """Reconstruct ``Delta_eq(T, q)`` from an official PRE-run Usadel NPZ.

    The official file written by ``save_usadel_catalog_npz`` stores
    ``metadata``, ``q_values_m_inv`` and calibration arrays, but not a full
    temperature-dependent equilibrium-gap table.  This function computes the
    requested curves using ``solve_gap_for_gamma_J`` with
    ``Gamma_q = hbar D q^2 / 2``.
    """

    path = Path(npz_path)
    arrays = _load_npz_arrays(path)
    metadata = _metadata_from_arrays(arrays)

    Tc_K = _catalog_scalar(arrays, metadata, ("Tc_K", "Tc"), required=True)
    T_bias_K = _catalog_scalar(arrays, metadata, ("T_bias_K", "T_K"), required=False, default=0.01 * Tc_K)
    D_m2_s = _catalog_scalar(arrays, metadata, ("D_m2_s", "D"), required=True)
    n_m = int(n_matsubara) if n_matsubara is not None else _catalog_int(
        arrays,
        metadata,
        ("n_matsubara_configured", "n_matsubara", "js_table_n_matsubara"),
        default=500,
    )
    if n_m <= 0:
        raise ValueError(f"n_matsubara must be positive, got {n_m}.")

    qcrit = _resolve_q_critical(arrays, metadata, q_critical_m_inv)
    q_values = np.linspace(0.0, qcrit, int(n_curves))
    if q_values.size < 2:
        raise ValueError("n_curves must be at least 2.")

    temperature = _temperature_axis(
        Tc_K=Tc_K,
        T_bias_K=T_bias_K,
        n_temperature=int(n_temperature),
        T_min_K=T_min_K,
        T_max_K=T_max_K,
    )

    gap_J = np.zeros((temperature.size, q_values.size), dtype=float)
    gamma_values = 0.5 * HBAR_J_S * D_m2_s * q_values * q_values

    for iT, T_K in enumerate(temperature):
        eps_n_J = matsubara_energy_axis_J(T_K=float(T_K), n_matsubara=int(n_m))
        for iq, gamma_J in enumerate(gamma_values):
            gap_J[iT, iq] = solve_gap_for_gamma_J(
                gamma_J=float(gamma_J),
                T_K=float(T_K),
                Tc_K=float(Tc_K),
                eps_n_J=eps_n_J,
            )

    compact_metadata: dict[str, Any] = {
        "source_npz": str(path),
        "Tc_K": float(Tc_K),
        "T_bias_K": float(T_bias_K),
        "D_m2_s": float(D_m2_s),
        "n_matsubara": int(n_m),
        "n_temperature": int(temperature.size),
        "n_curves": int(q_values.size),
    }

    return UsadelGapCatalog(
        temperature_K=temperature,
        q_m_inv=q_values,
        gap_meV=gap_J / _MEV_J,
        source_key="computed_from_matsubara_self_consistency",
        q_critical_m_inv=float(qcrit),
        metadata=compact_metadata,
    )


def load_usadel_gap_catalog(
    npz_path: str | Path,
    *,
    n_curves: int = 6,
    n_temperature: int = 160,
    T_min_K: float | None = None,
    T_max_K: float | None = None,
    q_critical_m_inv: float | None = None,
    n_matsubara: int | None = None,
) -> UsadelGapCatalog:
    """Compatibility wrapper for the E1 pipeline.

    The name says "load" because E1 reads an existing PRE-run, but the
    temperature-dependent equilibrium gap is reconstructed on demand.
    """

    return build_gap_eq_temperature_catalog(
        npz_path,
        n_curves=n_curves,
        n_temperature=n_temperature,
        T_min_K=T_min_K,
        T_max_K=T_max_K,
        q_critical_m_inv=q_critical_m_inv,
        n_matsubara=n_matsubara,
    )


def interpolate_gap_curves(
    catalog: UsadelGapCatalog,
    *,
    n_curves: int | None = None,
    q_critical_m_inv: float | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return q values, temperature grid and curves in ``(n_curves, n_T)`` layout."""

    temperature = np.asarray(catalog.temperature_K, dtype=float)
    q = np.asarray(catalog.q_m_inv, dtype=float)
    gap = np.asarray(catalog.gap_meV, dtype=float)
    if gap.shape != (temperature.size, q.size):
        raise ValueError(
            f"gap_meV shape {gap.shape} is not compatible with temperature/q axes "
            f"({temperature.size}, {q.size})."
        )

    qcrit = float(q_critical_m_inv) if q_critical_m_inv is not None else float(catalog.q_critical_m_inv)
    if n_curves is None or int(n_curves) == q.size:
        return q, temperature, gap.T

    q_targets = np.linspace(float(np.nanmin(q)), min(qcrit, float(np.nanmax(q))), int(n_curves))
    curves = np.empty((q_targets.size, temperature.size), dtype=float)
    q_order = np.argsort(q)
    q_sorted = q[q_order]
    gap_sorted = gap[:, q_order]
    for iT in range(temperature.size):
        curves[:, iT] = np.interp(q_targets, q_sorted, gap_sorted[iT, :])
    return q_targets, temperature, curves


def plot_gap_eq_vs_temperature(
    catalog: UsadelGapCatalog,
    output_path: str | Path,
    *,
    n_curves: int | None = None,
    q_critical_m_inv: float | None = None,
    dpi: int = 480,
    title: str | None = None,
) -> Path:
    """Plot ``Delta_eq`` versus temperature for multiple q values.

    E1 memory-ready style:
    - show 4 q-curves by default, even if the reconstructed catalog has 6,
    - extend the temperature axis beyond Tc so the Tc marker is visible,
    - add vertical Tc and horizontal BCS Delta(0) reference lines,
    - use larger labels/ticks/legend,
    - arrange legend as 2 columns x 3 rows.
    """

    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    # Requested: remove two of the original six q-curves.
    # If the caller asks for fewer than 4, respect that; otherwise cap at 4.
    n_plot_curves = 4 if n_curves is None else min(int(n_curves), 4)

    qcrit = float(q_critical_m_inv) if q_critical_m_inv is not None else float(catalog.q_critical_m_inv)

    q_targets, temperature, curves = interpolate_gap_curves(
        catalog,
        n_curves=n_plot_curves,
        q_critical_m_inv=qcrit,
    )

    try:
        Tc_K = float(catalog.metadata["Tc_K"])
    except Exception as exc:
        raise KeyError(
            "plot_gap_eq_vs_temperature requires catalog.metadata['Tc_K'] "
            "to draw the Tc and Delta_BCS reference lines."
        ) from exc

    if not np.isfinite(Tc_K) or Tc_K <= 0.0:
        raise ValueError(f"Invalid Tc_K in catalog metadata: {Tc_K!r}")

    delta_bcs_0_meV = 1.764 * Boltzmann * Tc_K / e * 1.0e3

    fig, ax = plt.subplots(figsize=(8.4, 5.8), constrained_layout=True)

    # Larger typography for thesis-ready PDF figures.
    label_fs = 17
    tick_fs = 14
    legend_fs = 12

    handles = []
    labels = []

    for q_value, curve in zip(q_targets, curves):
        ratio = q_value / qcrit if qcrit != 0.0 else np.nan
        (line,) = ax.plot(
            temperature,
            curve,
            linewidth=2.3,
            label=rf"$q/q_c={ratio:.2f}$",
        )
        handles.append(line)
        labels.append(rf"$q/q_c={ratio:.2f}$")

    tc_line = ax.axvline(
        Tc_K,
        linestyle="--",
        linewidth=2.1,
        label=rf"$T_c={Tc_K:.2f}\,\mathrm{{K}}$",
    )

    bcs_line = ax.axhline(
        delta_bcs_0_meV,
        linestyle="-.",
        linewidth=2.1,
        label=r"$\Delta_{\mathrm{BCS}}(0)=1.764\,k_B T_c$",
    )

    handles.extend([tc_line, bcs_line])
    labels.extend(
        [
            rf"$T_c={Tc_K:.2f}\,\mathrm{{K}}$",
            r"$\Delta_{\mathrm{BCS}}(0)=1.764\,k_B T_c$",
        ]
    )

    # Wider than Tc so the vertical dashed line is inside the plot, not on
    # the right border. The curves still stop at the computed temperature grid.
    T_min = float(np.nanmin(temperature))
    T_max_data = float(np.nanmax(temperature))
    T_max_plot = max(T_max_data, 1.08 * Tc_K)
    ax.set_xlim(T_min, T_max_plot)

    y_max_curves = float(np.nanmax(curves)) if curves.size else 0.0
    y_max_plot = max(y_max_curves, delta_bcs_0_meV) * 1.08
    ax.set_ylim(bottom=0.0, top=y_max_plot)

    ax.set_xlabel(r"Electron temperature $T_e$ (K)", fontsize=label_fs)
    ax.set_ylabel(r"Equilibrium gap $\Delta_{\mathrm{eq}}$ (meV)", fontsize=label_fs)
    ax.tick_params(axis="both", which="major", labelsize=tick_fs)

    ax.grid(True, alpha=0.25)

    # 2 columns x 3 rows:
    # row 1: q-curve 1, q-curve 2
    # row 2: q-curve 3, q-curve 4
    # row 3: Tc, Delta_BCS(0)
    ax.legend(
        handles,
        labels,
        frameon=False,
        fontsize=legend_fs,
        ncol=2,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.02),
        columnspacing=1.5,
        handlelength=2.8,
        handletextpad=0.7,
    )

    if title:
        ax.set_title(title, fontsize=label_fs)

    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output

__all__ = [
    "UsadelGapCatalog",
    "build_gap_eq_temperature_catalog",
    "interpolate_gap_curves",
    "load_usadel_gap_catalog",
    "plot_gap_eq_vs_temperature",
]
