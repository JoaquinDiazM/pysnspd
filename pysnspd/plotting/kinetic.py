"""Diagnostic plots for kinetic material functions and projected powers."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping

import numpy as np
import matplotlib.pyplot as plt


MEV_J = 1.602176634e-22


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


def plot_usadel_self_consistent_trajectory(
    trajectory,
    output_path: str | Path,
    *,
    dpi: int = 480,
) -> Path:
    """Plot the Usadel self-consistent Delta(q) trajectory used in OE5."""
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    frac = np.asarray(trajectory.current_fraction, dtype=float)
    delta_meV = np.asarray(trajectory.delta_eq_values_J, dtype=float) / MEV_J
    q = np.asarray(trajectory.q_values_m_inv, dtype=float)
    q_norm = q / np.max(q) if np.max(q) > 0.0 else q

    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    ax.plot(frac, delta_meV, linewidth=1.5, label=r"$\Delta_{\rm eq}(q)$ [meV]")
    ax.plot(frac, q_norm, linewidth=1.2, linestyle="--", label=r"$q/q_{\max}$")

    ax.set_title("Usadel self-consistent stable branch")
    ax.set_xlabel(r"normalized current $I/I_c$")
    ax.set_ylabel(r"$\Delta_{\rm eq}$ [meV] / normalized $q$")
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
    tau_label: str = "",
    title_suffix: str = "",
    dpi: int = 480,
) -> Path:
    """Plot projected electron-phonon powers versus Te at one Usadel state."""
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    Te = np.asarray(power_curve["Te_values_K"], dtype=float)

    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    ax.plot(Te, power_curve["P_S_W_m3"], linewidth=1.3, label=r"$P_{ep}^{S}$")
    ax.plot(Te, power_curve["P_R_W_m3"], linewidth=1.3, label=r"$P_{ep}^{R}$")
    ax.plot(
        Te,
        power_curve["P_total_W_m3"],
        linewidth=1.5,
        label=r"$P_{ep}^{S}+P_{ep}^{R}$",
    )

    debye_label = r"Vodolazov/Allmaras Debye $T^5$"
    if tau_label:
        debye_label += f" ({tau_label})"

    ax.plot(
        Te,
        power_curve["P_Debye_Vodolazov_W_m3"],
        linewidth=1.2,
        linestyle="--",
        label=debye_label,
    )

    ax.axhline(0.0, linewidth=0.8)
    title = "Projected electron-phonon power density"
    if title_suffix:
        title += f"\n{title_suffix}"
    ax.set_title(title)
    ax.set_xlabel(r"$T_e$ [K]")
    ax.set_ylabel(r"power density [W m$^{-3}$]")
    ax.grid(True, linewidth=0.25, alpha=0.35)
    ax.legend(frameon=True)

    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output


def plot_power_vs_usadel_current(
    q_scan: Mapping[str, np.ndarray],
    output_path: str | Path,
    *,
    dpi: int = 480,
) -> Path:
    """Plot projected powers along the self-consistent Usadel current branch."""
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    x = np.asarray(q_scan["current_fraction"], dtype=float)
    Te_values = np.asarray(q_scan["Te_values_K"], dtype=float)
    P_total = np.asarray(q_scan["P_total_W_m3"], dtype=float)
    P_R = np.asarray(q_scan["P_R_W_m3"], dtype=float)

    fig, ax = plt.subplots(figsize=(7.8, 5.0))

    for i, Te in enumerate(Te_values):
        ax.plot(
            x,
            P_total[i, :],
            linewidth=1.4,
            label=rf"total, $T_e={Te:.2f}$ K",
        )
        ax.plot(
            x,
            P_R[i, :],
            linewidth=1.0,
            linestyle="--",
            label=rf"$R$, $T_e={Te:.2f}$ K",
        )

    ax.axhline(0.0, linewidth=0.8)
    ax.set_title("Projected powers along Usadel self-consistent branch")
    ax.set_xlabel(r"normalized current $I/I_c$")
    ax.set_ylabel(r"power density [W m$^{-3}$]")
    ax.grid(True, linewidth=0.25, alpha=0.35)
    ax.legend(frameon=True, fontsize=8)

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