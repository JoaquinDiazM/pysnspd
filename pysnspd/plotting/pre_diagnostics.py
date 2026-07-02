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
import numpy as np

from pysnspd.mesh.delaunay import MeshData, triangle_areas
from pysnspd.mesh.edges import EdgeData

MEV_J = 1.602176634e-22


def write_pre_diagnostic_plots(
    *,
    mesh: MeshData,
    edge_data: EdgeData,
    usadel_catalog: Any,
    output_dir: str | Path,
    dpi: int = 480,
) -> dict[str, str]:
    """Write standard PRE diagnostic plots and return their paths.

    The mesh presentation figure is intentionally not created here. It belongs
    to ``plot_pipelines/01_plot_prerun.py``, which can be rerun without touching
    the raw PRE stage. This keeps the normal PRE-run focused on catalogue
    generation and lightweight sanity diagnostics.
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    paths = {
        "mesh_triangle_area_hist_png": plot_triangle_area_histogram(
            mesh,
            out / "mesh_triangle_area_hist.png",
            dpi=dpi,
        ),
        "mesh_edge_length_hist_png": plot_edge_length_histogram(
            edge_data,
            out / "mesh_edge_length_hist.png",
            dpi=dpi,
        ),
        "usadel_supercurrent_curve_png": plot_usadel_supercurrent_curve(
            usadel_catalog,
            out / "usadel_supercurrent_curve.png",
            dpi=dpi,
        ),
        "usadel_equilibrium_dos_map_png": plot_usadel_equilibrium_dos_map(
            usadel_catalog,
            out / "usadel_equilibrium_dos_map.png",
            dpi=dpi,
        ),
        "usadel_zero_energy_dos_map_png": plot_usadel_zero_energy_dos_map(
            usadel_catalog,
            out / "usadel_zero_energy_dos_map.png",
            dpi=dpi,
        ),
        "usadel_equilibrium_anomalous_map_png": plot_usadel_equilibrium_anomalous_map(
            usadel_catalog,
            out / "usadel_equilibrium_anomalous_map.png",
            dpi=dpi,
        ),
    }
    return {key: str(value) for key, value in paths.items()}


def plot_mesh_boundary_tags(
    mesh: MeshData,
    edge_data: EdgeData,
    output_path: str | Path,
    *,
    dpi: int = 480,
) -> Path:
    """Plot the triangulation and overlay boundary-edge tags.

    Kept as a callable helper for ad-hoc debugging, but no longer part of the
    default PRE-run diagnostic set. The official PRE plotting pipeline now uses
    the pyTDGL-style mesh figure instead.
    """
    output = _prepare_output(output_path)
    nodes_nm = np.asarray(mesh.nodes, dtype=float) * 1.0e9

    fig, ax = plt.subplots(figsize=(8.0, 4.2))
    ax.set_aspect("equal")
    ax.triplot(
        nodes_nm[:, 0],
        nodes_nm[:, 1],
        np.asarray(mesh.triangles, dtype=np.int64),
        linewidth=0.25,
        alpha=0.35,
    )

    boundary = np.asarray(edge_data.is_boundary, dtype=bool)
    tags = np.asarray(edge_data.tags).astype(str)
    edges = np.asarray(edge_data.edges, dtype=np.int64)

    ordered_tags = ["left", "right", "bottom", "top", "boundary_unknown"]
    for tag in ordered_tags:
        mask = boundary & (tags == tag)
        if not np.any(mask):
            continue
        _plot_edge_segments(
            ax,
            nodes_nm,
            edges[mask],
            label=f"{tag} ({int(np.count_nonzero(mask))})",
            linewidth=1.0,
        )

    unknown_boundary = boundary & ~np.isin(tags, ordered_tags)
    if np.any(unknown_boundary):
        _plot_edge_segments(
            ax,
            nodes_nm,
            edges[unknown_boundary],
            label=f"other boundary ({int(np.count_nonzero(unknown_boundary))})",
            linewidth=1.0,
        )

    ax.set_title("PRE mesh: Delaunay triangulation and boundary tags")
    ax.set_xlabel("x [nm]")
    ax.set_ylabel("y [nm]")
    ax.legend(loc="best", fontsize=8)
    ax.grid(False)
    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output


def plot_triangle_area_histogram(
    mesh: MeshData,
    output_path: str | Path,
    *,
    dpi: int = 480,
) -> Path:
    """Plot the distribution of triangle areas in nm^2."""
    output = _prepare_output(output_path)
    areas_nm2 = triangle_areas(mesh.nodes, mesh.triangles) * 1.0e18

    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    ax.hist(areas_nm2, bins=_safe_histogram_bins(areas_nm2, max_bins=60))
    ax.axvline(float(np.mean(areas_nm2)), linestyle="--", linewidth=1.0, label="mean")
    ax.set_title("PRE mesh: triangle area distribution")
    ax.set_xlabel(r"triangle area [nm$^2$]")
    ax.set_ylabel("count")
    ax.legend(loc="best", fontsize=8)
    ax.grid(True, linewidth=0.35, alpha=0.25)
    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output


def plot_edge_length_histogram(
    edge_data: EdgeData,
    output_path: str | Path,
    *,
    dpi: int = 480,
) -> Path:
    """Plot the interior/boundary edge-length distributions in nm."""
    output = _prepare_output(output_path)
    lengths_nm = np.asarray(edge_data.lengths, dtype=float) * 1.0e9
    boundary = np.asarray(edge_data.is_boundary, dtype=bool)

    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    if np.any(~boundary):
        interior_lengths = lengths_nm[~boundary]
        ax.hist(
            interior_lengths,
            bins=_safe_histogram_bins(interior_lengths, max_bins=60),
            alpha=0.65,
            label="interior",
        )
    if np.any(boundary):
        boundary_lengths = lengths_nm[boundary]
        ax.hist(
            boundary_lengths,
            bins=_safe_histogram_bins(boundary_lengths, max_bins=60),
            alpha=0.65,
            label="boundary",
        )
    ax.axvline(float(np.mean(lengths_nm)), linestyle="--", linewidth=1.0, label="mean")
    ax.set_title("PRE mesh: edge-length distribution")
    ax.set_xlabel("edge length [nm]")
    ax.set_ylabel("count")
    ax.legend(loc="best", fontsize=8)
    ax.grid(True, linewidth=0.35, alpha=0.25)
    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output


def plot_usadel_supercurrent_curve(
    usadel_catalog: Any,
    output_path: str | Path,
    *,
    dpi: int = 480,
) -> Path:
    """Plot Usadel calibration current and equilibrium gap versus superfluid momentum."""
    output = _prepare_output(output_path)

    q = np.asarray(usadel_catalog.calibration_q_values_m_inv, dtype=float)
    current_uA = 1.0e6 * np.asarray(usadel_catalog.calibration_current_values_A, dtype=float)
    delta_eq_meV = _joule_to_mev(np.asarray(usadel_catalog.calibration_delta_eq_values_J, dtype=float))
    q_1e7 = q / 1.0e7

    finite_i = np.isfinite(q_1e7) & np.isfinite(current_uA)
    finite_d = np.isfinite(q_1e7) & np.isfinite(delta_eq_meV)
    if not np.any(finite_i):
        raise ValueError("Usadel calibration current table has no finite points.")

    metadata = getattr(usadel_catalog, "metadata", {})
    T_bias_K = _metadata_float(metadata, "T_bias_K")
    Tc_K = _metadata_float(metadata, "Tc_K")
    temperature_label = _temperature_label(T_bias_K, Tc_K)

    fig, ax_i = plt.subplots(figsize=(7.1, 4.35))
    line_i, = ax_i.plot(
        q_1e7[finite_i],
        current_uA[finite_i],
        marker=".",
        markersize=1.9,
        linewidth=1.05,
        label=rf"$I_s(q)$ at {temperature_label}",
    )

    q_plot = q_1e7[finite_i]
    i_plot = current_uA[finite_i]
    i_max = int(np.nanargmax(i_plot))
    peak_line, = ax_i.plot(
        q_plot[i_max],
        i_plot[i_max],
        marker="o",
        markersize=4.0,
        linestyle="None",
        label=rf"model $I_c$ = {i_plot[i_max]:.2f} $\mu$A",
    )

    target_line = None
    target = metadata.get("Ic_target_A") if isinstance(metadata, dict) else None
    if target is not None:
        target_uA = 1.0e6 * float(target)
        target_line = ax_i.axhline(
            target_uA,
            linestyle="--",
            linewidth=1.0,
            label=rf"target $I_c$ = {target_uA:.2f} $\mu$A",
        )

    ax_d = ax_i.twinx()
    gap_line = None
    if np.any(finite_d):
        gap_line, = ax_d.plot(
            q_1e7[finite_d],
            delta_eq_meV[finite_d],
            marker=".",
            markersize=1.7,
            linewidth=1.0,
            label=rf"$|\Delta_{{eq}}(q)|$ at {temperature_label}",
        )

    ax_i.set_title("Usadel/Matsubara calibration: current and equilibrium gap")
    ax_i.set_xlabel(r"superfluid momentum $q$ [$10^7$ m$^{-1}$]")
    ax_i.set_ylabel(r"$I_s$ [$\mu$A]")
    ax_d.set_ylabel(r"$|\Delta_{eq}|$ [meV]")
    ax_i.grid(True, linewidth=0.35, alpha=0.28)

    handles = [line_i, peak_line]
    if target_line is not None:
        handles.append(target_line)
    if gap_line is not None:
        handles.append(gap_line)
    labels = [h.get_label() for h in handles]

    extra = _usadel_metadata_summary(metadata)
    legend_title = "Usadel calibration"
    if extra:
        legend_title += "\n" + extra
    ax_i.legend(handles, labels, loc="best", fontsize=7.6, title=legend_title, title_fontsize=7.5)

    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output


def plot_usadel_equilibrium_dos_map(
    usadel_catalog: Any,
    output_path: str | Path,
    *,
    dpi: int = 480,
) -> Path:
    """Plot rho(E,q) along the equilibrium gap branch Delta_eq(q)."""
    output = _prepare_output(output_path)
    rho_eq = _catalog_field_on_equilibrium_gap(usadel_catalog, field_name="rho_delta_gamma_E")
    q_1e7 = np.asarray(usadel_catalog.q_values_m_inv, dtype=float) / 1.0e7
    E_meV = _joule_to_mev(np.asarray(usadel_catalog.energy_values_J, dtype=float))
    metadata = getattr(usadel_catalog, "metadata", {})

    fig, ax = plt.subplots(figsize=(7.1, 4.35))
    im = ax.pcolormesh(E_meV, q_1e7, rho_eq, shading="auto", vmin=0.0)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(r"$\rho(E,q)$")
    ax.set_title(r"Usadel DOS along $\Delta_{eq}(q)$" + _temperature_suffix(metadata))
    ax.set_xlabel(r"energy $E$ [meV]")
    ax.set_ylabel(r"$q$ [$10^7$ m$^{-1}$]")
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
    """Plot the zero-energy DOS over the independent (Delta,q) catalogue axes."""
    output = _prepare_output(output_path)
    rho = np.asarray(usadel_catalog.rho_delta_gamma_E, dtype=float)
    energy = np.asarray(usadel_catalog.energy_values_J, dtype=float)
    idx0 = int(np.nanargmin(np.abs(energy)))
    q_1e7 = np.asarray(usadel_catalog.q_values_m_inv, dtype=float) / 1.0e7
    delta_meV = _joule_to_mev(np.asarray(usadel_catalog.delta_values_J, dtype=float))
    zero_map = rho[:, :, idx0]

    fig, ax = plt.subplots(figsize=(7.1, 4.35))
    im = ax.pcolormesh(q_1e7, delta_meV, zero_map, shading="auto", vmin=0.0)
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


def plot_usadel_equilibrium_anomalous_map(
    usadel_catalog: Any,
    output_path: str | Path,
    *,
    dpi: int = 480,
) -> Path:
    """Plot |anomalous(E,q)| along the equilibrium gap branch Delta_eq(q)."""
    output = _prepare_output(output_path)
    anomalous_eq = np.abs(
        _catalog_field_on_equilibrium_gap(usadel_catalog, field_name="anomalous_delta_gamma_E")
    )
    q_1e7 = np.asarray(usadel_catalog.q_values_m_inv, dtype=float) / 1.0e7
    E_meV = _joule_to_mev(np.asarray(usadel_catalog.energy_values_J, dtype=float))
    metadata = getattr(usadel_catalog, "metadata", {})

    fig, ax = plt.subplots(figsize=(7.1, 4.35))
    im = ax.pcolormesh(E_meV, q_1e7, anomalous_eq, shading="auto", vmin=0.0)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(r"$|F(E,q)|$")
    ax.set_title(r"Usadel anomalous amplitude along $\Delta_{eq}(q)$" + _temperature_suffix(metadata))
    ax.set_xlabel(r"energy $E$ [meV]")
    ax.set_ylabel(r"$q$ [$10^7$ m$^{-1}$]")
    ax.grid(False)
    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output


def _catalog_field_on_equilibrium_gap(usadel_catalog: Any, *, field_name: str) -> np.ndarray:
    """Return field[q,E] by sampling the catalogue at nearest Delta_eq(q)."""
    field = np.asarray(getattr(usadel_catalog, field_name), dtype=float)
    if field.ndim != 3:
        raise ValueError(f"{field_name} must have shape (n_delta, n_q, n_energy).")

    delta_axis = np.asarray(usadel_catalog.delta_values_J, dtype=float)
    q_axis = np.asarray(usadel_catalog.q_values_m_inv, dtype=float)
    q_cal = np.asarray(usadel_catalog.calibration_q_values_m_inv, dtype=float)
    delta_cal = np.asarray(usadel_catalog.calibration_delta_eq_values_J, dtype=float)

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


def _plot_edge_segments(
    ax: Any,
    nodes_nm: np.ndarray,
    edges: np.ndarray,
    *,
    label: str,
    linewidth: float,
) -> None:
    """Plot many edge segments using the current Matplotlib property cycle."""
    first = True
    for i, j in np.asarray(edges, dtype=np.int64):
        p = nodes_nm[[int(i), int(j)]]
        ax.plot(
            p[:, 0],
            p[:, 1],
            linewidth=linewidth,
            label=label if first else None,
        )
        first = False


def _metadata_float(metadata: Any, key: str) -> float:
    if not isinstance(metadata, dict) or key not in metadata:
        return float("nan")
    try:
        return float(metadata[key])
    except Exception:
        return float("nan")


def _temperature_label(T_bias_K: float, Tc_K: float) -> str:
    if np.isfinite(T_bias_K) and np.isfinite(Tc_K) and Tc_K > 0.0:
        return rf"$T={T_bias_K:.2f}$ K ($T/T_c={T_bias_K / Tc_K:.3f}$)"
    if np.isfinite(T_bias_K):
        return rf"$T={T_bias_K:.2f}$ K"
    return "configured bias temperature"


def _temperature_suffix(metadata: Any) -> str:
    T_bias_K = _metadata_float(metadata, "T_bias_K")
    Tc_K = _metadata_float(metadata, "Tc_K")
    if np.isfinite(T_bias_K) and np.isfinite(Tc_K) and Tc_K > 0.0:
        return rf" at $T={T_bias_K:.2f}$ K ($T/T_c={T_bias_K / Tc_K:.3f}$)"
    if np.isfinite(T_bias_K):
        return rf" at $T={T_bias_K:.2f}$ K"
    return ""


def _usadel_metadata_summary(metadata: Any) -> str:
    if not isinstance(metadata, dict):
        return ""
    parts: list[str] = []
    D = _metadata_float(metadata, "D_m2_s")
    sigma = _metadata_float(metadata, "sigma_n_S_m")
    delta0 = _metadata_float(metadata, "delta0_meV")
    if np.isfinite(D):
        parts.append(rf"$D={D:.3g}$ m$^2$/s")
    if np.isfinite(sigma):
        parts.append(rf"$\sigma_n={sigma:.3g}$ S/m")
    if np.isfinite(delta0):
        parts.append(rf"$\Delta_0={delta0:.3f}$ meV")
    return ", ".join(parts)


def _joule_to_mev(values_J: np.ndarray | float) -> np.ndarray | float:
    return np.asarray(values_J, dtype=float) / MEV_J


def _prepare_output(path: str | Path) -> Path:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    return output
