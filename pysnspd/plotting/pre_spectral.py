"""PRE-run diagnostic plots for mesh, edges and Usadel calibration.

These plots are intentionally deterministic and inexpensive compared with the
catalogue construction itself. The normal PRE-run uses this module for compact
raw diagnostics. Presentation-quality mesh figures are handled by the dedicated
PRE plotting pipeline through ``pysnspd.plotting.mesh.plot_mesh_pytdgl_style``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
import numpy as np

from pysnspd.usadel.calibration import matsubara_energy_axis_J, solve_gap_for_gamma_J
from pysnspd.usadel.parameters import HBAR_J_S
from pysnspd.mesh.delaunay import MeshData, triangle_areas
from pysnspd.mesh.edges import EdgeData
from pysnspd.plotting.style import THESIS_DOUBLE_FIGSIZE, apply_thesis_style

MEV_J = 1.602176634e-22

def _auto_energy_xlim_from_tail(
    E_meV: np.ndarray,
    values_qE: np.ndarray,
    *,
    rel_threshold: float = 0.015,
    tail_fraction: float = 0.15,
    pad_fraction: float = 0.10,
    min_visible_fraction: float = 0.20,
) -> tuple[float, float]:
    """
    Choose an automatic energy window by detecting where the data still differs
    from its high-energy asymptotic tail.

    This is intended for Usadel DOS/anomalous maps where a long high-energy
    domain becomes visually uninformative.
    """
    E = np.asarray(E_meV, dtype=float)
    values = np.asarray(values_qE, dtype=float)

    if E.ndim != 1 or values.ndim < 2 or values.shape[-1] != E.size:
        return float(np.nanmin(E)), float(np.nanmax(E))

    finite_E = np.isfinite(E)
    if np.count_nonzero(finite_E) < 4:
        return float(np.nanmin(E)), float(np.nanmax(E))

    # Assume the last part of the energy grid is the asymptotic tail.
    nE = E.size
    n_tail = max(3, int(np.ceil(tail_fraction * nE)))
    n_tail = min(n_tail, nE)

    tail = values[..., -n_tail:]
    tail_level = np.nanmedian(tail, axis=-1, keepdims=True)

    # Activity as a function of energy: max deviation from the high-energy tail.
    activity = np.nanmax(np.abs(values - tail_level), axis=0)

    if not np.any(np.isfinite(activity)):
        return float(np.nanmin(E)), float(np.nanmax(E))

    max_activity = float(np.nanmax(activity))
    if not np.isfinite(max_activity) or max_activity <= 0.0:
        return float(np.nanmin(E)), float(np.nanmax(E))

    threshold = rel_threshold * max_activity
    active = np.flatnonzero(activity >= threshold)

    if active.size == 0:
        return float(np.nanmin(E)), float(np.nanmax(E))

    i0 = int(active[0])
    i1 = int(active[-1])

    # For positive energy grids, keep the plot anchored at the minimum energy.
    # For symmetric grids, keep both lower and upper active bounds.
    pad_points = max(2, int(np.ceil(pad_fraction * max(1, i1 - i0 + 1))))

    if np.nanmin(E) >= 0.0:
        i0_plot = 0
    else:
        i0_plot = max(0, i0 - pad_points)

    i1_plot = min(nE - 1, i1 + pad_points)

    # Avoid over-cropping if the threshold is too aggressive.
    min_visible_points = max(4, int(np.ceil(min_visible_fraction * nE)))
    if (i1_plot - i0_plot + 1) < min_visible_points:
        missing = min_visible_points - (i1_plot - i0_plot + 1)
        i1_plot = min(nE - 1, i1_plot + missing)

    return float(E[i0_plot]), float(E[i1_plot])


def _energy_visible_mask(E_meV: np.ndarray, xlim: tuple[float, float]) -> np.ndarray:
    """Mask energy columns inside the plotted x-limits."""
    E = np.asarray(E_meV, dtype=float)
    lo, hi = xlim
    mask = np.isfinite(E) & (E >= lo) & (E <= hi)
    if np.count_nonzero(mask) < 3:
        return np.isfinite(E)
    return mask

def plot_usadel_equilibrium_dos_map(
    usadel_catalog: Any,
    output_path: str | Path,
    *,
    dpi: int = 480,
    energy_window: bool = True,
    energy_tail_rel_threshold: float = 0.1
) -> Path:
    """Plot rho(E,q) along the equilibrium gap branch Delta_eq(q)."""
    output = _prepare_output(output_path)

    rho_eq = _catalog_field_on_equilibrium_gap(
        usadel_catalog,
        field_name="rho_delta_gamma_E",
    )

    # Use the q-axis that matches the equilibrium branch map. In normal PRE-runs
    # this is q_values_m_inv, but smoke/legacy objects may be closer to the
    # calibration axis.
    q_axis_m_inv = np.asarray(getattr(usadel_catalog, "q_values_m_inv"), dtype=float)
    if q_axis_m_inv.size != rho_eq.shape[0] and hasattr(
        usadel_catalog,
        "calibration_q_values_m_inv",
    ):
        q_axis_m_inv = np.asarray(
            usadel_catalog.calibration_q_values_m_inv,
            dtype=float,
        )

    q_1e7 = q_axis_m_inv / 1.0e7
    E_meV = _joule_to_mev(np.asarray(usadel_catalog.energy_values_J, dtype=float))
    metadata = getattr(usadel_catalog, "metadata", {})

    T_bias_K = _metadata_float(metadata, "T_bias_K")
    Tc_K = _metadata_float(metadata, "Tc_K")

    if np.isfinite(T_bias_K) and np.isfinite(Tc_K) and Tc_K > 0.0:
        title = (
            rf"Usadel DOS along $\Delta_{{\rm eq}}(q)$ "
            rf"at $T={T_bias_K:.2f}$ [K] "
            rf"$(T/T_c={T_bias_K / Tc_K:.3f})$"
        )
    elif np.isfinite(T_bias_K):
        title = rf"Usadel DOS along $\Delta_{{\rm eq}}(q)$ at $T={T_bias_K:.2f}$ [K]"
    else:
        title = r"Usadel DOS along $\Delta_{\rm eq}(q)$"

    if energy_window:
        xlim = _auto_energy_xlim_from_tail(
            E_meV,
            rho_eq,
            rel_threshold=0.05,
            tail_fraction=0.15,
            pad_fraction=0.04,
            min_visible_fraction=0.08,
        )
    else:
        xlim = (float(np.nanmin(E_meV)), float(np.nanmax(E_meV)))

    visible_E = _energy_visible_mask(E_meV, xlim)

    fig, ax = plt.subplots(figsize=(7.1, 4.35))

    # Normalize using only the visible energy range, so the long rho -> 1 tail
    # does not dominate the visual scale.
    norm = _positive_power_norm(rho_eq[:, visible_E], gamma=0.35)

    im = ax.pcolormesh(
        E_meV,
        q_1e7,
        rho_eq,
        shading="auto",
        norm=norm,
    )

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(r"$\rho(E,q)$")

    ax.set_title(title)
    ax.set_xlabel(r"energy $E$ [meV]")
    ax.set_ylabel(r"$q$ [$10^7$ m$^{-1}$]")
    ax.set_xlim(*xlim)
    ax.grid(False)

    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output

def plot_usadel_equilibrium_anomalous_map(
    usadel_catalog: Any,
    output_path: str | Path,
    *,
    dpi: int = 480,
    energy_window: bool = True,
    energy_tail_rel_threshold: float = 0.1,
) -> Path:
    """Plot |anomalous(E,q)| along the equilibrium gap branch Delta_eq(q)."""
    output = _prepare_output(output_path)

    anomalous_eq = np.abs(
        _catalog_field_on_equilibrium_gap(
            usadel_catalog,
            field_name="anomalous_delta_gamma_E",
        )
    )

    q_1e7 = np.asarray(usadel_catalog.q_values_m_inv, dtype=float) / 1.0e7
    E_meV = _joule_to_mev(np.asarray(usadel_catalog.energy_values_J, dtype=float))
    metadata = getattr(usadel_catalog, "metadata", {})

    if energy_window:
        xlim = _auto_energy_xlim_from_tail(
            E_meV,
            anomalous_eq,
            rel_threshold=0.05,
            tail_fraction=0.15,
            pad_fraction=0.04,
            min_visible_fraction=0.08,
        )
    else:
        xlim = (float(np.nanmin(E_meV)), float(np.nanmax(E_meV)))

    visible_E = _energy_visible_mask(E_meV, xlim)

    fig, ax = plt.subplots(figsize=(7.1, 4.35))

    # Normalize using only the visible energy range. This prevents the hidden
    # asymptotic tail from affecting the color scale.
    norm = _positive_power_norm(anomalous_eq[:, visible_E], gamma=0.35)

    im = ax.pcolormesh(
        E_meV,
        q_1e7,
        anomalous_eq,
        shading="auto",
        norm=norm,
    )

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(r"$|F(E,q)|$ (power scale, $\gamma=0.35$)")

    ax.set_title(
        r"Usadel anomalous amplitude along $\Delta_{\rm eq}(q)$"
        + _temperature_suffix(metadata)
    )
    ax.set_xlabel(r"energy $E$ [meV]")
    ax.set_ylabel(r"$q$ [$10^7$ m$^{-1}$]")
    ax.set_xlim(*xlim)
    ax.grid(False)

    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output


def plot_usadel_zero_energy_dos_map(
    usadel_catalog: Any,
    output_path: str | Path,
    *,
    dpi: int = 480,
) -> Path:
    """Plot the zero-energy DOS over the independent (Delta,q) catalogue axes.

    Use an interpolated regular-image rendering instead of raw cell colors so the
    map is less blocky and easier to read for coarse catalogue grids.
    """
    output = _prepare_output(output_path)
    rho = np.asarray(usadel_catalog.rho_delta_gamma_E, dtype=float)
    energy = np.asarray(usadel_catalog.energy_values_J, dtype=float)
    idx0 = int(np.nanargmin(np.abs(energy)))
    q_1e7 = np.asarray(usadel_catalog.q_values_m_inv, dtype=float) / 1.0e7
    delta_meV = _joule_to_mev(np.asarray(usadel_catalog.delta_values_J, dtype=float))
    zero_map = rho[:, :, idx0]

    fig, ax = plt.subplots(figsize=(7.1, 4.35))
    im = ax.imshow(
        zero_map,
        origin="lower",
        aspect="auto",
        interpolation="bicubic",
        extent=_imshow_extent(q_1e7, delta_meV),
        vmin=0.0,
        vmax=float(np.nanmax(zero_map)) if np.isfinite(np.nanmax(zero_map)) else 1.0,
        cmap="viridis",
    )
    # Add faint contours to recover the underlying catalogue structure without the
    # visually harsh pixelation of a plain nearest-cell rendering.
    if np.isfinite(zero_map).any():
        levels = np.linspace(float(np.nanmin(zero_map)), float(np.nanmax(zero_map)), 7)
        if np.nanmax(levels) > np.nanmin(levels):
            ax.contour(q_1e7, delta_meV, zero_map, levels=levels[1:-1], colors="white", linewidths=0.35, alpha=0.35)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(r"$\rho(E\approx0;|\Delta|,q)$")
    ax.set_title("Usadel zero-energy DOS over catalogue grid")
    ax.set_xlabel(r"$q$ [$10^7$ m$^{-1}$]")
    ax.set_ylabel(r"$|\Delta|$ [meV]")
    ax.grid(False)
    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output


def plot_usadel_equilibrium_gap_Tq_map(
    usadel_catalog: Any,
    output_path: str | Path,
    *,
    dpi: int = 480,
    usadel_npz_path: str | Path | None = None,
) -> Path:
    """Plot the equilibrium gap over the available temperature--q grid.

    The PRE catalogue may carry a strict 3D Matsubara supercurrent table. When a
    corresponding equilibrium-gap table is available, this helper visualizes the
    equilibrium gap magnitude over that grid. For legacy/lean smoke objects, it
    falls back to a single-temperature row built from the calibration branch.

    The full table can extend far above Tc or beyond the depairing momentum,
    where |Delta_eq| is exactly zero. For visualization, the axes are cropped to
    the active superconducting region.
    """
    output = _prepare_output(output_path)
    T_vals_K, q_vals_m_inv, delta_eq_values_J = _extract_equilibrium_gap_Tq_data(
        usadel_catalog,
        usadel_npz_path=usadel_npz_path,
    )

    T_vals_K = np.asarray(T_vals_K, dtype=float)
    q_1e7 = np.asarray(q_vals_m_inv, dtype=float) / 1.0e7
    delta_eq_meV = _joule_to_mev(np.asarray(delta_eq_values_J, dtype=float))

    fig, ax = plt.subplots(figsize=(7.1, 4.35))
    vmax = float(np.nanmax(delta_eq_meV)) if np.isfinite(delta_eq_meV).any() else 1.0
    vmax = max(vmax, 1.0e-12)

    im = ax.imshow(
        delta_eq_meV,
        origin="lower",
        aspect="auto",
        interpolation="bilinear",
        extent=_imshow_extent(q_1e7, T_vals_K),
        vmin=0.0,
        vmax=vmax,
        cmap="viridis",
    )

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(r"$|\Delta_{eq}(T,q)|$ [meV]")

    ax.set_title(r"Equilibrium gap over the $(T,q)$ calibration grid")
    ax.set_xlabel(r"$q$ [$10^7$ m$^{-1}$]")
    ax.set_ylabel(r"temperature $T$ [K]")

    # Crop to the active superconducting region. This removes the large
    # normal-state area where |Delta_eq| = 0 and the plot otherwise wastes space.
    active_threshold = 1.0e-3 * vmax
    active = np.isfinite(delta_eq_meV) & (delta_eq_meV > active_threshold)

    if np.any(active) and delta_eq_meV.ndim == 2:
        active_T_mask = np.any(active, axis=1)
        active_q_mask = np.any(active, axis=0)

        active_T = T_vals_K[active_T_mask]
        active_q = q_1e7[active_q_mask]

        if active_T.size and active_q.size:
            dT = (
                float(np.nanmedian(np.diff(T_vals_K)))
                if T_vals_K.size > 1
                else 0.5
            )
            dq = (
                float(np.nanmedian(np.diff(q_1e7)))
                if q_1e7.size > 1
                else 0.5
            )

            T_min = max(
                float(np.nanmin(T_vals_K)),
                float(np.nanmin(active_T)) - dT,
            )
            T_max = min(
                float(np.nanmax(T_vals_K)),
                float(np.nanmax(active_T)) + dT,
            )
            q_min = max(
                float(np.nanmin(q_1e7)),
                float(np.nanmin(active_q)) - dq,
            )
            q_max = min(
                float(np.nanmax(q_1e7)),
                float(np.nanmax(active_q)) + dq,
            )

            if T_max > T_min:
                ax.set_ylim(T_min, T_max)
            if q_max > q_min:
                ax.set_xlim(q_min, q_max)

    ax.grid(False)
    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output

def _calibration_delta_eq_mev(usadel_catalog: Any) -> np.ndarray:
    """Return calibration equilibrium gap in meV, or NaNs for legacy smoke objects."""
    if hasattr(usadel_catalog, "calibration_delta_eq_values_J"):
        return _joule_to_mev(np.asarray(usadel_catalog.calibration_delta_eq_values_J, dtype=float))
    q = np.asarray(getattr(usadel_catalog, "calibration_q_values_m_inv"), dtype=float)
    return np.full(q.shape, np.nan, dtype=float)


def _positive_power_norm(values: np.ndarray, *, gamma: float = 0.35) -> PowerNorm:
    """Nonlinear positive normalization that keeps zero finite and compresses singular peaks."""
    arr = np.asarray(values, dtype=float)
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return PowerNorm(gamma=gamma, vmin=0.0, vmax=1.0)
    vmax = float(np.nanpercentile(finite, 99.7))
    vmax = max(vmax, float(np.nanmax(finite)), 1.0e-12)
    return PowerNorm(gamma=gamma, vmin=0.0, vmax=vmax)


def _catalog_field_on_equilibrium_gap(usadel_catalog: Any, *, field_name: str) -> np.ndarray:
    """Return field[q,E] by sampling the catalogue at nearest Delta_eq(q)."""
    field = np.asarray(getattr(usadel_catalog, field_name), dtype=float)
    if field.ndim != 3:
        raise ValueError(f"{field_name} must have shape (n_delta, n_q, n_energy).")

    delta_axis = np.asarray(usadel_catalog.delta_values_J, dtype=float)
    q_axis = np.asarray(usadel_catalog.q_values_m_inv, dtype=float)
    q_cal = np.asarray(usadel_catalog.calibration_q_values_m_inv, dtype=float)
    if hasattr(usadel_catalog, "calibration_delta_eq_values_J"):
        delta_cal = np.asarray(usadel_catalog.calibration_delta_eq_values_J, dtype=float)
    else:
        delta_cal = np.asarray([], dtype=float)

    if q_cal.size and delta_cal.size and q_cal.size == delta_cal.size:
        order = np.argsort(q_cal)
        delta_eq = np.interp(q_axis, q_cal[order], delta_cal[order], left=delta_cal[order][0], right=delta_cal[order][-1])
    else:
        delta_eq = np.full(q_axis.shape, float(np.nanmax(delta_axis)), dtype=float)

    sampled = np.empty((q_axis.size, field.shape[2]), dtype=float)
    for iq, delta in enumerate(delta_eq):
        idelta = int(np.nanargmin(np.abs(delta_axis - delta)))
        sampled[iq, :] = field[idelta, iq, :]
    return sampled


def _extract_equilibrium_gap_Tq_data(
    usadel_catalog: Any,
    *,
    usadel_npz_path: str | Path | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Extract or reconstruct the self-consistent equilibrium gap on the table axes.

    ``load_usadel_catalog_npz`` intentionally returns the core OE3 spectral
    catalogue dataclass. The strict 3D current table appended later to the same
    NPZ (``Te_axis_K``, ``delta_axis_J``, ``q_axis_m_inv``, ``js_A_m2`` and
    aliases) is not exposed as dataclass attributes. The PRE plotting pipeline
    therefore passes the raw NPZ path so this plot can read those table axes
    directly, without changing the PRE-run file format or the catalogue loader.
    """
    temperature_candidates = [
        "js_table_temperature_values_K",
        "js_table_temperatures_K",
        "Te_axis_K",
        "temperature_values_K",
        "Te_values_K",
        "temperatures_K",
        "supercurrent_table_temperature_values_K",
        "strict_3d_temperature_values_K",
        "js_temperature_values_K",
    ]
    q_candidates = [
        "js_table_q_values_m_inv",
        "q_axis_m_inv",
        "q_values_m_inv",
        "calibration_q_values_m_inv",
        "supercurrent_table_q_values_m_inv",
        "strict_3d_q_values_m_inv",
    ]
    table_delta_axis_candidates = [
        "delta_axis_J",
        "js_table_delta_axis_J",
        "delta_values_J",
        "supercurrent_table_delta_axis_J",
        "strict_3d_delta_axis_J",
    ]
    delta_table_candidates = [
        "js_table_delta_eq_values_J",
        "js_table_equilibrium_gap_values_J",
        "delta_eq_Tq_values_J",
        "equilibrium_gap_Tq_values_J",
        "supercurrent_table_delta_eq_values_J",
        "strict_3d_delta_eq_values_J",
        "js_delta_eq_values_J",
        "js_table_delta_eq_J",
    ]

    npz_payload = _load_npz_payload_for_gap_map(usadel_npz_path)
    metadata = dict(getattr(usadel_catalog, "metadata", {}) or {})
    if isinstance(npz_payload.get("metadata"), dict):
        metadata.update(npz_payload["metadata"])

    # First try a future explicit Delta_eq(T,q) table, either attached to the
    # object or stored directly in the NPZ.
    T_vals = _find_first_attr(usadel_catalog, temperature_candidates)
    q_vals = _find_first_attr(usadel_catalog, q_candidates)
    delta_map = _find_first_attr(usadel_catalog, delta_table_candidates)
    if T_vals is None:
        T_vals = _find_first_mapping_value(npz_payload, temperature_candidates)
    if q_vals is None:
        q_vals = _find_first_mapping_value(npz_payload, q_candidates)
    if delta_map is None:
        delta_map = _find_first_mapping_value(npz_payload, delta_table_candidates)

    if T_vals is not None and q_vals is not None and delta_map is not None:
        T_vals_arr = np.asarray(T_vals, dtype=float)
        q_vals_arr = np.asarray(q_vals, dtype=float)
        delta_map_arr = np.asarray(delta_map, dtype=float)
        if delta_map_arr.ndim == 2 and delta_map_arr.shape == (T_vals_arr.size, q_vals_arr.size):
            return T_vals_arr, q_vals_arr, delta_map_arr
        if delta_map_arr.ndim == 2 and delta_map_arr.shape == (q_vals_arr.size, T_vals_arr.size):
            return T_vals_arr, q_vals_arr, delta_map_arr.T

    # Current PRE-runs store j_s(T,Delta,q), not Delta_eq(T,q). Reconstruct the
    # equilibrium branch on the real table axes and snap it to the stored Delta
    # axis. Reading the axes from the raw NPZ is the important step here: the
    # frozen UsadelCatalog dataclass does not carry them.
    delta_axis = _find_first_attr(usadel_catalog, table_delta_axis_candidates)
    if delta_axis is None:
        delta_axis = _find_first_mapping_value(npz_payload, table_delta_axis_candidates)

    if T_vals is not None and q_vals is not None and delta_axis is not None:
        T_vals_arr = np.asarray(T_vals, dtype=float)
        q_vals_arr = np.asarray(q_vals, dtype=float)
        delta_axis_arr = np.asarray(delta_axis, dtype=float)
        inferred = _infer_self_consistent_gap_from_table_axes(
            usadel_catalog=usadel_catalog,
            T_vals_K=T_vals_arr,
            q_vals_m_inv=q_vals_arr,
            delta_axis_J=delta_axis_arr,
            metadata=metadata,
        )
        if inferred is not None:
            return T_vals_arr, q_vals_arr, inferred

    # Legacy fallback: use only the calibration branch at the configured bias temperature.
    q_cal = np.asarray(getattr(usadel_catalog, "calibration_q_values_m_inv"), dtype=float)
    delta_cal = np.asarray(_calibration_delta_eq_mev(usadel_catalog), dtype=float) * MEV_J
    T_bias = _metadata_float(metadata, "T_bias_K")
    if not np.isfinite(T_bias):
        T_bias = 0.0
    return np.asarray([T_bias], dtype=float), q_cal, delta_cal[None, :]


def _infer_self_consistent_gap_from_table_axes(
    *,
    usadel_catalog: Any,
    T_vals_K: np.ndarray,
    q_vals_m_inv: np.ndarray,
    delta_axis_J: np.ndarray,
    metadata: Any,
) -> np.ndarray | None:
    """Compute Delta_eq(T,q) on existing table axes and snap to the Delta grid."""
    D_m2_s = _metadata_float(metadata, "D_m2_s")
    Tc_K = _metadata_float(metadata, "Tc_K")
    n_matsubara = _integer_attr_or_metadata(
        usadel_catalog,
        metadata,
        attr_candidates=("js_table_n_matsubara", "n_matsubara", "n_matsubara_configured"),
        metadata_candidates=("js_table_n_matsubara", "n_matsubara", "n_matsubara_configured"),
        default=500,
    )
    if not (np.isfinite(D_m2_s) and D_m2_s > 0.0 and np.isfinite(Tc_K) and Tc_K > 0.0):
        return None

    T_axis = np.asarray(T_vals_K, dtype=float).reshape(-1)
    q_axis = np.asarray(q_vals_m_inv, dtype=float).reshape(-1)
    delta_axis = np.asarray(delta_axis_J, dtype=float).reshape(-1)
    delta_axis = delta_axis[np.isfinite(delta_axis) & (delta_axis >= 0.0)]
    if T_axis.size == 0 or q_axis.size == 0 or delta_axis.size == 0:
        return None
    delta_axis = np.unique(delta_axis)
    delta_axis.sort()

    out = np.zeros((T_axis.size, q_axis.size), dtype=float)
    gamma_axis_J = 0.5 * HBAR_J_S * float(D_m2_s) * q_axis * q_axis
    for iT, T in enumerate(T_axis):
        if not np.isfinite(T) or T <= 0.0 or T >= Tc_K:
            continue
        eps_n_J = matsubara_energy_axis_J(T_K=float(T), n_matsubara=int(n_matsubara))
        for iq, gamma in enumerate(gamma_axis_J):
            if not np.isfinite(gamma) or gamma < 0.0:
                continue
            delta_cont = solve_gap_for_gamma_J(
                gamma_J=float(gamma),
                T_K=float(T),
                Tc_K=float(Tc_K),
                eps_n_J=eps_n_J,
            )
            if delta_cont <= 0.0:
                out[iT, iq] = 0.0
                continue
            # The map represents the equilibrium gap available to interpolation on the
            # stored table, so snap the continuous Matsubara solution to the nearest
            # PRE-run Delta-axis value.
            idx = int(np.nanargmin(np.abs(delta_axis - delta_cont)))
            out[iT, iq] = float(delta_axis[idx])
    return out


def _integer_attr_or_metadata(
    obj: Any,
    metadata: Any,
    *,
    attr_candidates: tuple[str, ...],
    metadata_candidates: tuple[str, ...],
    default: int,
) -> int:
    for name in attr_candidates:
        if hasattr(obj, name):
            try:
                value = np.asarray(getattr(obj, name)).reshape(-1)[0]
                return int(value)
            except Exception:
                pass
    if isinstance(metadata, dict):
        for name in metadata_candidates:
            if name in metadata:
                try:
                    return int(metadata[name])
                except Exception:
                    pass
        supercurrent = metadata.get("supercurrent_table")
        if isinstance(supercurrent, dict):
            for name in metadata_candidates:
                if name in supercurrent:
                    try:
                        return int(supercurrent[name])
                    except Exception:
                        pass
    return int(default)


def _load_npz_payload_for_gap_map(path: str | Path | None) -> dict[str, Any]:
    """Read only small table axes/metadata needed by the gap map from the raw NPZ."""
    if path is None:
        return {}
    source = Path(path)
    if not source.exists():
        return {}
    wanted = {
        "Te_axis_K",
        "delta_axis_J",
        "q_axis_m_inv",
        "js_A_m2",
        "j_s_A_m2",
        "js_T_delta_q_A_m2",
        "metadata",
    }
    out: dict[str, Any] = {}
    with np.load(source, allow_pickle=True) as data:
        for key in data.files:
            if key in wanted or key.endswith("_axis_K") or key.endswith("_axis_J") or key.endswith("_axis_m_inv") or "delta_eq" in key:
                value = data[key]
                if key == "metadata":
                    try:
                        out[key] = value.item()
                    except Exception:
                        pass
                else:
                    out[key] = np.asarray(value)
    return out


def _find_first_mapping_value(mapping: dict[str, Any], candidates: list[str]) -> Any | None:
    for name in candidates:
        if name in mapping:
            value = mapping[name]
            if value is not None:
                return value
    return None


def _find_first_attr(obj: Any, candidates: list[str]) -> Any | None:
    for name in candidates:
        if hasattr(obj, name):
            value = getattr(obj, name)
            if value is not None:
                return value
    return None


def _imshow_extent(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float, float]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size == 1:
        dx = max(abs(float(x[0])) * 0.05, 0.5)
        xmin, xmax = float(x[0] - dx), float(x[0] + dx)
    else:
        dx = np.diff(x)
        xmin = float(x[0] - 0.5 * dx[0])
        xmax = float(x[-1] + 0.5 * dx[-1])
    if y.size == 1:
        dy = max(abs(float(y[0])) * 0.05, 0.5)
        ymin, ymax = float(y[0] - dy), float(y[0] + dy)
    else:
        dy = np.diff(y)
        ymin = float(y[0] - 0.5 * dy[0])
        ymax = float(y[-1] + 0.5 * dy[-1])
    return xmin, xmax, ymin, ymax


def _safe_histogram_bins(values: np.ndarray, *, max_bins: int = 60) -> np.ndarray:
    """Return finite histogram edges, including for constant large-valued data."""
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        raise ValueError("Cannot plot histogram: no finite values were provided.")
    vmin = float(np.min(finite))
    vmax = float(np.max(finite))
    scale = max(abs(vmin), abs(vmax), 1.0)
    span = vmax - vmin
    if (not np.isfinite(span)) or span <= 16.0 * np.finfo(float).eps * scale:
        center = float(np.mean(finite))
        pad = max(1.0e-6 * max(abs(center), 1.0), 1.0e-12)
        return np.array([center - pad, center + pad], dtype=float)
    n_bins = int(min(max_bins, max(1, finite.size)))
    return np.linspace(vmin, vmax, n_bins + 1, dtype=float)


def _metadata_float(metadata: Any, key: str) -> float:
    if not isinstance(metadata, dict) or key not in metadata:
        return float("nan")
    try:
        return float(metadata[key])
    except Exception:
        return float("nan")


def _temperature_label(T_bias_K: float, Tc_K: float) -> str:
    if np.isfinite(T_bias_K) and np.isfinite(Tc_K) and Tc_K > 0.0:
        #return rf"$T={T_bias_K:.2f}$ K ($T/T_c={T_bias_K / Tc_K:.3f}$)"
        return rf"$T={T_bias_K:.2f}$ [K]"
    if np.isfinite(T_bias_K):
        return rf"$T={T_bias_K:.2f}$ [K]"
    return "configured bias temperature"


def _temperature_suffix(metadata: Any) -> str:
    T_bias_K = _metadata_float(metadata, "T_bias_K")
    Tc_K = _metadata_float(metadata, "Tc_K")
    if np.isfinite(T_bias_K) and np.isfinite(Tc_K) and Tc_K > 0.0:
        return rf" at $T={T_bias_K:.2f}$ [K] ($T/T_c={T_bias_K / Tc_K:.3f}$)"
    if np.isfinite(T_bias_K):
        return rf" at $T={T_bias_K:.2f}$ [K]"
    return ""


def _joule_to_mev(values_J: np.ndarray | float) -> np.ndarray | float:
    return np.asarray(values_J, dtype=float) / MEV_J


def _prepare_output(path: str | Path) -> Path:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    return output
