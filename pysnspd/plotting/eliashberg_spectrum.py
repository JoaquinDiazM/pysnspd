"""Presentation-quality Eliashberg/PhDOS plots for E-type pipelines.

The material spectrum is loaded by ``pysnspd.kinetic.eliashberg`` from the
Simon/MIT three-column file and is plotted here without imposing the high-energy
cutoffs used by phase-space or power-table integrations.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np


def plot_eliashberg_spectrum(
    spectrum: Any,
    output_path: str | Path,
    *,
    dpi: int = 480,
) -> Path:
    """Plot the full Eliashberg spectral function and phonon DOS.

    Parameters
    ----------
    spectrum:
        ``EliashbergSpectrum`` loaded by ``load_simon_eliashberg_dat``.
    output_path:
        Destination PDF path. The complete loaded frequency/energy axis is used;
        no high-frequency cut is applied here.
    dpi:
        Rasterization DPI for PDF backends.
    """
    output = _prepare_output(output_path)

    omega_meV = np.asarray(spectrum.omega_meV, dtype=float)
    alpha2F = np.asarray(spectrum.alpha2F, dtype=float)
    phdos = np.asarray(spectrum.phdos_states_per_THz, dtype=float)

    finite = np.isfinite(omega_meV) & np.isfinite(alpha2F) & np.isfinite(phdos)
    if not np.any(finite):
        raise ValueError("Eliashberg spectrum has no finite points to plot.")

    omega_meV = omega_meV[finite]
    alpha2F = alpha2F[finite]
    phdos = phdos[finite]

    order = np.argsort(omega_meV)
    omega_meV = omega_meV[order]
    alpha2F = alpha2F[order]
    phdos = phdos[order]

    alpha_color = "tab:blue"
    phdos_color = "tab:purple"

    fig, ax_a = plt.subplots(figsize=(9.4, 5.9))
    ax_p = ax_a.twinx()

    label_fs = 20
    tick_fs = 16
    legend_fs = 13.5
    legend_title_fs = 13.0

    line_a, = ax_a.plot(
        omega_meV,
        alpha2F,
        linewidth=2.2,
        color=alpha_color,
        label=r"$\alpha^2F(\Omega)$",
        zorder=3,
    )
    line_p, = ax_p.plot(
        omega_meV,
        phdos,
        linewidth=2.0,
        color=phdos_color,
        label=r"PhDOS $F(\Omega)$",
        zorder=3,
    )

    ax_a.set_xlabel(r"phonon energy $\Omega$ [meV]", fontsize=label_fs)
    ax_a.set_ylabel(r"$\alpha^2F(\Omega)$", fontsize=label_fs, color=alpha_color)
    ax_p.set_ylabel(r"PhDOS [states/THz]", fontsize=label_fs, color=phdos_color)

    ax_a.tick_params(axis="x", which="both", direction="in", labelsize=tick_fs)
    ax_a.tick_params(axis="y", which="both", direction="in", labelsize=tick_fs, colors=alpha_color)
    ax_p.tick_params(axis="y", which="both", direction="in", labelsize=tick_fs, colors=phdos_color)
    ax_a.minorticks_on()
    ax_p.minorticks_on()

    ax_a.spines["left"].set_color(alpha_color)
    ax_a.spines["left"].set_linewidth(1.2)
    ax_p.spines["right"].set_color(phdos_color)
    ax_p.spines["right"].set_linewidth(1.2)
    ax_a.yaxis.label.set_color(alpha_color)
    ax_p.yaxis.label.set_color(phdos_color)

    ax_a.grid(True, which="major", linewidth=0.45, alpha=0.22, zorder=1)
    ax_a.grid(True, which="minor", linewidth=0.25, alpha=0.08, zorder=1)

    x_min = max(0.0, float(np.nanmin(omega_meV)))
    x_max = float(np.nanmax(omega_meV))
    ax_a.set_xlim(x_min, x_max)

    high_energy_cut_meV = 35.0
    high_energy_band = None

    if x_max > high_energy_cut_meV:
        high_energy_band = ax_a.axvspan(
            high_energy_cut_meV,
            x_max,
            facecolor="0.88",
            edgecolor="none",
            alpha=0.55,
            zorder=0,
            label=r"high-energy band, $\Omega \geq 35$ meV",
        )

        ax_a.axvline(
            high_energy_cut_meV,
            color="0.35",
            linestyle="--",
            linewidth=1.1,
            alpha=0.75,
            zorder=2,
        )

    alpha_top = 1.08 * float(np.nanmax(alpha2F)) if alpha2F.size else 1.0
    phdos_top = 1.08 * float(np.nanmax(phdos)) if phdos.size else 1.0
    ax_a.set_ylim(bottom=0.0, top=max(alpha_top, 1.0e-12))
    ax_p.set_ylim(bottom=0.0, top=max(phdos_top, 1.0e-12))

    metadata = getattr(spectrum, "metadata", {})
    title_lines: list[str] = []
    lambda_ep = float(getattr(spectrum, "lambda_ep", np.nan))
    if np.isfinite(lambda_ep):
        title_lines.append(rf"$\lambda_{{ep}}={lambda_ep:.3f}$")
    title_lines.append("Full Simon/MIT NbN spectrum")
    legend_title = "\n".join(title_lines)

    legend_handles = [line_a, line_p]
    legend_labels = [line_a.get_label(), line_p.get_label()]

    if high_energy_band is not None:
        legend_handles.append(high_energy_band)
        legend_labels.append(r"high-energy band")

    legend = ax_a.legend(
        legend_handles,
        legend_labels,
        loc="upper right",
        fontsize=legend_fs,
        title=legend_title,
        title_fontsize=legend_title_fs,
        frameon=True,
        fancybox=False,
        framealpha=1.0,
        facecolor="white",
        edgecolor="black",
        handlelength=2.7,
        borderpad=0.58,
        labelspacing=0.45,
    )
    legend.get_frame().set_linewidth(1.0)
    legend.set_zorder(10)
    legend.get_title().set_multialignment("center")
    legend._legend_box.align = "center"

    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output


def _prepare_output(output_path: str | Path) -> Path:
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    return output


__all__ = ["plot_eliashberg_spectrum"]
