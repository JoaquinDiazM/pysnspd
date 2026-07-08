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

MEV_J = 1.602176634e-22


def write_pre_diagnostic_plots(
    *,
    mesh: MeshData,
    edge_data: EdgeData,
    usadel_catalog: Any,
    output_dir: str | Path,
    dpi: int = 480,
    usadel_npz_path: str | Path | None = None,
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
        "usadel_equilibrium_gap_Tq_map_png": plot_usadel_equilibrium_gap_Tq_map(
            usadel_catalog,
            out / "usadel_equilibrium_gap_Tq_map.png",
            dpi=dpi,
            usadel_npz_path=usadel_npz_path,
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
    delta_eq_meV = _calibration_delta_eq_mev(usadel_catalog)

    q_1e7 = q / 1.0e7
    finite_i = np.isfinite(q_1e7) & np.isfinite(current_uA)
    finite_d = np.isfinite(q_1e7) & np.isfinite(delta_eq_meV)
    if not np.any(finite_i):
        raise ValueError("Usadel calibration current table has no finite points.")

    metadata = getattr(usadel_catalog, "metadata", {})
    T_bias_K = _metadata_float(metadata, "T_bias_K")
    Tc_K = _metadata_float(metadata, "Tc_K")
    temperature_label = _temperature_label(T_bias_K, Tc_K)

    current_color = "tab:blue"
    gap_color = "tab:purple"

    fig, ax_i = plt.subplots(figsize=(9.4, 5.9))
    label_fs = 18
    tick_fs = 16
    legend_fs = 14.0
    legend_title_fs = 14.0

    q_plot = q_1e7[finite_i]
    i_plot = current_uA[finite_i]
    i_max = int(np.nanargmax(i_plot))
    ic_uA = float(i_plot[i_max])

    line_i, = ax_i.plot(
        q_plot,
        i_plot,
        linewidth=2.2,
        color=current_color,
        zorder=3,
        label=rf"$I_s(q)$ at {temperature_label}",
    )

    target_line = None
    target = metadata.get("Ic_target_A") if isinstance(metadata, dict) else None
    if target is not None:
        target_uA = 1.0e6 * float(target)
        target_line = ax_i.axhline(
            target_uA,
            linestyle="--",
            linewidth=1.6,
            color=current_color,
            alpha=0.95,
            zorder=2,
            label=rf"target $I_c$ = {target_uA:.2f} $\mu$A",
        )

    ax_d = ax_i.twinx()
    gap_line = None
    if np.any(finite_d):
        gap_line, = ax_d.plot(
            q_1e7[finite_d],
            delta_eq_meV[finite_d],
            linewidth=2.0,
            color=gap_color,
            zorder=3,
            label=rf"$|\Delta_{{eq}}(q)|$ at {temperature_label}",
        )

    ax_i.set_xlabel(r"superfluid momentum $q$ [$10^7$ m$^{-1}$]", fontsize=label_fs)
    ax_i.set_ylabel(r"$I_s$ [$\mu$A]", fontsize=label_fs, color=current_color)
    ax_d.set_ylabel(r"$|\Delta_{eq}|$ [meV]", fontsize=label_fs, color=gap_color)

    ax_i.tick_params(axis="x", which="both", direction="in", labelsize=tick_fs)
    ax_i.tick_params(axis="y", which="both", direction="in", labelsize=tick_fs, colors=current_color)
    ax_d.tick_params(axis="y", which="both", direction="in", labelsize=tick_fs, colors=gap_color)

    ax_i.minorticks_on()
    ax_d.minorticks_on()

    ax_i.spines["left"].set_color(current_color)
    ax_i.spines["left"].set_linewidth(1.2)
    ax_d.spines["right"].set_color(gap_color)
    ax_d.spines["right"].set_linewidth(1.2)

    ax_i.yaxis.label.set_color(current_color)
    ax_d.yaxis.label.set_color(gap_color)

    ax_i.grid(True, linewidth=0.4, alpha=0.25, zorder=1)

    handles = [line_i]
    if target_line is not None:
        handles.append(target_line)
    if gap_line is not None:
        handles.append(gap_line)

    labels = [h.get_label() for h in handles]

    def _find_diffusivity_value(meta: Any) -> float | None:
        if not isinstance(meta, dict):
            return None

        candidate_keys = (
            "D_m2_s",
            "D",
            "diffusivity_m2_s",
            "diffusion_constant_m2_s",
        )

        for key in candidate_keys:
            if key in meta:
                try:
                    value = float(meta[key])
                    if np.isfinite(value):
                        return value
                except Exception:
                    pass

        for block_key in ("calibration", "material", "supercurrent_table"):
            block = meta.get(block_key)
            if isinstance(block, dict):
                for key in candidate_keys:
                    if key in block:
                        try:
                            value = float(block[key])
                            if np.isfinite(value):
                                return value
                        except Exception:
                            pass

        return None

    legend_lines = []

    D_value = _find_diffusivity_value(metadata)
    if D_value is not None:
        D_cm2_s = 1.0e4 * D_value
        legend_lines.append(rf"Calibration $\Rightarrow$ $D={D_cm2_s:.3f}$ cm$^2$/s")

    # Esta es la línea que agregaba el resumen:
    # D = ..., sigma_n = ..., Delta_0 = ...
    # Se deja comentada para no mostrar esa segunda línea en la leyenda.
    #
    # extra = _usadel_metadata_summary(metadata)
    # if extra:
    #     for line in extra.splitlines():
    #         clean = line.strip()
    #         if not clean:
    #             continue
    #         clean = clean.replace("(", "").replace(")", "")
    #         if clean:
    #             legend_lines.append(clean)

    legend_title = "\n".join(legend_lines) if legend_lines else None

    legend = ax_i.legend(
        handles,
        labels,
        loc="lower center",
        bbox_to_anchor=(0.50, 0.055),
        fontsize=legend_fs,
        title=legend_title,
        title_fontsize=legend_title_fs,
        frameon=True,
        fancybox=False,
        framealpha=1.0,
        facecolor="white",
        edgecolor="black",
        handlelength=2.8,
        borderpad=0.60,
        labelspacing=0.45,
    )

    legend.get_frame().set_linewidth(1.0)
    legend.set_zorder(10)

    # Centra horizontalmente la primera línea, que en Matplotlib corresponde
    # al título de la leyenda.
    legend.get_title().set_ha("center")
    legend.get_title().set_multialignment("center")
    legend._legend_box.align = "center"

    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output

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
            rf"at $T={T_bias_K:.2f}$ K "
            rf"$(T/T_c={T_bias_K / Tc_K:.3f})$"
        )
    elif np.isfinite(T_bias_K):
        title = rf"Usadel DOS along $\Delta_{{\rm eq}}(q)$ at $T={T_bias_K:.2f}$ K"
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
        #return rf"$T={T_bias_K:.2f}$ K ($T/T_c={T_bias_K / Tc_K:.3f}$)"
        return rf"$T={T_bias_K:.2f}$ K"
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
        parts.append(rf"$D={1.0e4 * D:.3g}$ cm$^2$/s")
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
