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
from pysnspd.plotting.style import THESIS_DOUBLE_FIGSIZE, THESIS_DPI, apply_thesis_style

apply_thesis_style()

MEV_J = 1.602176634e-22


def write_pre_diagnostic_plots(
    *,
    mesh: MeshData,
    edge_data: EdgeData,
    usadel_catalog: Any,
    output_dir: str | Path,
    dpi: int = THESIS_DPI,
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


def plot_triangle_area_histogram(
    mesh: MeshData,
    output_path: str | Path,
    *,
    dpi: int = THESIS_DPI,
) -> Path:
    """Plot the distribution of triangle areas in nm^2."""
    output = _prepare_output(output_path)
    areas_nm2 = triangle_areas(mesh.nodes, mesh.triangles) * 1.0e18

    fig, ax = plt.subplots(figsize=THESIS_DOUBLE_FIGSIZE)
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
    dpi: int = THESIS_DPI,
) -> Path:
    """Plot the interior/boundary edge-length distributions in nm."""
    output = _prepare_output(output_path)
    lengths_nm = np.asarray(edge_data.lengths, dtype=float) * 1.0e9
    boundary = np.asarray(edge_data.is_boundary, dtype=bool)

    fig, ax = plt.subplots(figsize=THESIS_DOUBLE_FIGSIZE)
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
    dpi: int = THESIS_DPI,
) -> Path:
    """Plot Usadel calibration current and equilibrium gap versus superfluid momentum."""
    apply_thesis_style()
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

    fig, ax_i = plt.subplots(figsize=THESIS_DOUBLE_FIGSIZE)
    label_fs = 12
    tick_fs = 10
    legend_fs = 9.0
    legend_title_fs = 9.0

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
            label=rf"Target $I_c={target_uA:.2f}$ [$\mu$A]",
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

    ax_i.set_xlabel(r"Superfluid momentum $q$ [$10^7$ m$^{-1}$]", fontsize=label_fs)
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
        legend_lines.append(rf"Calibrated $D={D_cm2_s:.3f}$ [cm$^2$ s$^{{-1}}$]")

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
        bbox_to_anchor=(0.50, 0.035),
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

from pysnspd.plotting.pre_spectral import (
    _calibration_delta_eq_mev,
    _metadata_float,
    _prepare_output,
    _safe_histogram_bins,
    _temperature_label,
    plot_usadel_equilibrium_anomalous_map,
    plot_usadel_equilibrium_dos_map,
    plot_usadel_equilibrium_gap_Tq_map,
    plot_usadel_zero_energy_dos_map,
)
