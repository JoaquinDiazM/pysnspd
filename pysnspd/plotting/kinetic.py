"""Diagnostic plots for kinetic material functions and projected powers."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping

import numpy as np
import matplotlib.pyplot as plt


def plot_eliashberg_spectrum(spectrum, output_path: str | Path, *, dpi: int = 480) -> Path:
    """Plot normalized alpha^2F and PhDOS from a Simon/MIT material file."""
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    omega = spectrum.omega_meV
    alpha = np.asarray(spectrum.alpha2F, dtype=float)
    phdos = np.asarray(spectrum.phdos_states_per_THz, dtype=float)

    alpha_norm = alpha / np.max(alpha) if np.max(alpha) > 0.0 else alpha
    phdos_norm = phdos / np.max(phdos) if np.max(phdos) > 0.0 else phdos

    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    ax.plot(omega, alpha_norm, linewidth=1.3, label=r"$\alpha^2F(\Omega)$ normalized")
    ax.plot(omega, phdos_norm, linewidth=1.3, label="PhDOS normalized")

    ax.set_title("Kinetic phonon material functions")
    ax.set_xlabel(r"phonon energy $\Omega$ [meV]")
    ax.set_ylabel("normalized value")
    ax.grid(True, linewidth=0.25, alpha=0.35)
    ax.legend(frameon=True)

    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output


def plot_power_curve(
    power_curve: Mapping[str, np.ndarray],
    output_path: str | Path,
    *,
    dpi: int = 480,
) -> Path:
    """Plot projected electron-phonon powers versus Te."""
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    Te = np.asarray(power_curve["Te_values_K"], dtype=float)

    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    ax.plot(Te, power_curve["P_S_W_m3"], linewidth=1.3, label=r"$P_{ep}^{S}$")
    ax.plot(Te, power_curve["P_R_W_m3"], linewidth=1.3, label=r"$P_{ep}^{R}$")
    ax.plot(Te, power_curve["P_total_W_m3"], linewidth=1.5, label=r"$P_{ep}^{S}+P_{ep}^{R}$")
    ax.plot(
        Te,
        power_curve["P_Debye_Vodolazov_W_m3"],
        linewidth=1.2,
        linestyle="--",
        label=r"Vodolazov/Allmaras Debye $T^5$",
    )

    ax.axhline(0.0, linewidth=0.8)
    ax.set_title("Projected electron-phonon power density")
    ax.set_xlabel(r"$T_e$ [K]")
    ax.set_ylabel(r"power density [W m$^{-3}$]")
    ax.grid(True, linewidth=0.25, alpha=0.35)
    ax.legend(frameon=True)

    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output


def plot_spectral_support(
    support: Mapping[str, np.ndarray],
    output_path: str | Path,
    *,
    dpi: int = 480,
) -> Path:
    """Plot cumulative spectral support of OE5 integrands."""
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    omega = np.asarray(support["omega_meV"], dtype=float)

    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    ax.plot(
        omega,
        support["cumulative_alpha_omega"],
        linewidth=1.3,
        label=r"cumulative $|\Omega\alpha^2F|$",
    )
    ax.plot(
        omega,
        support["cumulative_scattering"],
        linewidth=1.3,
        label=r"cumulative $|P^S$ integrand|",
    )
    ax.plot(
        omega,
        support["cumulative_recombination"],
        linewidth=1.3,
        label=r"cumulative $|P^R$ integrand|",
    )

    ax.set_ylim(-0.02, 1.02)
    ax.set_title("Cumulative spectral support")
    ax.set_xlabel(r"phonon energy $\Omega$ [meV]")
    ax.set_ylabel("cumulative fraction")
    ax.grid(True, linewidth=0.25, alpha=0.35)
    ax.legend(frameon=True)

    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output