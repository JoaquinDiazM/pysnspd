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


def plot_thermal_usadel_gap_grid(
    thermal_grid,
    output_path: str | Path,
    *,
    target_fraction: float | None = None,
    dpi: int = 480,
) -> Path:
    """Plot Delta_eq(Te,q) for representative current fractions."""
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    Te = np.asarray(thermal_grid.Te_values_K, dtype=float)
    frac = np.asarray(thermal_grid.reference_current_fraction, dtype=float)
    delta = np.asarray(thermal_grid.delta_eq_Tq_J, dtype=float) / MEV_J

    targets = [0.0, 0.5, 0.8, 0.95]
    if target_fraction is not None:
        targets.append(float(target_fraction))

    indices = []
    for f in targets:
        idx = int(np.argmin(np.abs(frac - f)))
        if idx not in indices:
            indices.append(idx)

    fig, ax = plt.subplots(figsize=(7.5, 4.8))

    for idx in indices:
        ax.plot(
            Te,
            delta[:, idx],
            linewidth=1.3,
            label=rf"$I/I_c^{{bias}}\approx {frac[idx]:.3f}$",
        )

    ax.axvline(
        float(thermal_grid.metadata["Tc_K"]),
        linewidth=1.0,
        linestyle=":",
        label=rf"$T_c={thermal_grid.metadata['Tc_K']:.2f}$ K",
    )

    ax.set_title(r"Thermal Usadel self-consistency: $\Delta_{\rm eq}(T_e,q)$")
    ax.set_xlabel(r"$T_e$ [K]")
    ax.set_ylabel(r"$\Delta_{\rm eq}$ [meV]")
    ax.grid(True, linewidth=0.25, alpha=0.35)
    ax.legend(frameon=True, fontsize=8)

    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output


def plot_power_curve_thermal_usadel(
    power_curve: Mapping[str, np.ndarray],
    output_path: str | Path,
    *,
    tau_label: str = "",
    title_suffix: str = "",
    dpi: int = 480,
) -> Path:
    """Plot projected powers versus Te using Delta_eq(Te,q)."""
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    Te = np.asarray(power_curve["Te_values_K"], dtype=float)

    fig, ax = plt.subplots(figsize=(7.6, 4.9))

    ax.plot(Te, power_curve["P_S_W_m3"], linewidth=1.3, label=r"$P_{ep}^{S}$")
    ax.plot(Te, power_curve["P_R_W_m3"], linewidth=1.3, label=r"$P_{ep}^{R}$")
    ax.plot(
        Te,
        power_curve["P_total_W_m3"],
        linewidth=1.5,
        label=r"$P_{ep}^{S}+P_{ep}^{R}$",
    )
    ax.plot(
        Te,
        power_curve["P_normal_Eliashberg_W_m3"],
        linewidth=1.2,
        linestyle="-.",
        label=r"normal Eliashberg $\Delta=0$",
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
    ax.legend(frameon=True, fontsize=8)

    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output


def plot_power_ratios_thermal_usadel(
    power_curve: Mapping[str, np.ndarray],
    output_path: str | Path,
    *,
    dpi: int = 480,
) -> Path:
    """Plot ratios against the Vodolazov/Allmaras Debye reference."""
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    Te = np.asarray(power_curve["Te_values_K"], dtype=float)
    debye = np.asarray(power_curve["P_Debye_Vodolazov_W_m3"], dtype=float)

    def ratio(y):
        y = np.asarray(y, dtype=float)
        out = np.full_like(y, np.nan)
        mask = np.abs(debye) > 0.0
        out[mask] = y[mask] / debye[mask]
        return out

    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    ax.plot(Te, ratio(power_curve["P_S_W_m3"]), linewidth=1.3, label=r"$P_S/P_D$")
    ax.plot(
        Te,
        ratio(power_curve["P_total_W_m3"]),
        linewidth=1.3,
        label=r"$(P_S+P_R)/P_D$",
    )
    ax.plot(
        Te,
        ratio(power_curve["P_normal_Eliashberg_W_m3"]),
        linewidth=1.2,
        linestyle="-.",
        label=r"$P_{\Delta=0}^{\rm Eliashberg}/P_D$",
    )

    ax.axhline(1.0, linewidth=0.8, linestyle=":")
    ax.set_title("Power-density ratios relative to Debye reference")
    ax.set_xlabel(r"$T_e$ [K]")
    ax.set_ylabel("ratio")
    ax.grid(True, linewidth=0.25, alpha=0.35)
    ax.legend(frameon=True, fontsize=8)

    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output


def plot_power_scan_thermal_usadel(
    scan: Mapping[str, np.ndarray],
    output_path: str | Path,
    *,
    dpi: int = 480,
) -> Path:
    """Plot powers along the thermal Usadel q branch.

    A semilog scale is used because low-temperature curves otherwise collapse
    visually under the largest Te curve.
    """
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    x = np.asarray(scan["reference_current_fraction"], dtype=float)
    Te_values = np.asarray(scan["Te_values_K"], dtype=float)
    P_total = np.asarray(scan["P_total_W_m3"], dtype=float)
    P_R = np.asarray(scan["P_R_W_m3"], dtype=float)

    floor = 1.0e-100

    fig, ax = plt.subplots(figsize=(7.8, 5.0))

    for i, Te in enumerate(Te_values):
        total = np.where(np.abs(P_total[i, :]) > floor, np.abs(P_total[i, :]), np.nan)
        recomb = np.where(np.abs(P_R[i, :]) > floor, np.abs(P_R[i, :]), np.nan)

        ax.semilogy(
            x,
            total,
            linewidth=1.4,
            label=rf"total, $T_e={Te:.2f}$ K",
        )
        ax.semilogy(
            x,
            recomb,
            linewidth=1.0,
            linestyle="--",
            label=rf"$R$, $T_e={Te:.2f}$ K",
        )

    ax.set_title("Projected powers along thermal Usadel branch")
    ax.set_xlabel(r"reference normalized current $I(q,T_{bias})/I_c(T_{bias})$")
    ax.set_ylabel(r"$|P|$ [W m$^{-3}$]")
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


def plot_low_energy_recombination_scattering_band(
    result,
    output_path: str | Path,
    *,
    omega_max_meV: float | None = None,
    dpi: int = 480,
) -> Path:
    """Compare cumulative S/R powers in the low-energy gap-scale band.

    The goal is not to compare the full integrated powers, but to check whether
    the superconducting recombination channel is comparable to scattering in
    the low-energy region where Omega is of order a few Delta.

    This plot is only meaningful when Delta > 0.
    """
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    omega_J = np.asarray(result.omega_J, dtype=float)
    omega_meV = omega_J / MEV_J
    delta_meV = float(result.delta_J / MEV_J)

    if omega_max_meV is None:
        omega_max_meV = max(16.0, 4.0 * delta_meV + 6.0, 8.0 * delta_meV)

    mask = omega_meV <= float(omega_max_meV)
    if np.sum(mask) < 3:
        mask = np.ones_like(omega_meV, dtype=bool)

    om = omega_J[mask]
    x = omega_meV[mask]

    s_integrand = np.asarray(result.integrand_S_J2[mask], dtype=float)
    r_integrand = np.asarray(result.integrand_R_J2[mask], dtype=float)

    # Local cumulative trapezoids without scipy.
    d_om = np.diff(om)
    cS = np.concatenate(
        [[0.0], np.cumsum(0.5 * (s_integrand[1:] + s_integrand[:-1]) * d_om)]
    )
    cR = np.concatenate(
        [[0.0], np.cumsum(0.5 * (r_integrand[1:] + r_integrand[:-1]) * d_om)]
    )

    pref_S = 8.0 * np.pi * float(result.N0_J_m3) / 1.054571817e-34
    pref_R = 4.0 * np.pi * float(result.N0_J_m3) / 1.054571817e-34

    PS_partial = pref_S * cS
    PR_partial = pref_R * cR

    fig, ax = plt.subplots(figsize=(7.6, 4.8))

    ax.plot(x, PS_partial, linewidth=1.4, label=r"partial $P_{ep}^{S}(\Omega_c)$")
    ax.plot(x, PR_partial, linewidth=1.4, label=r"partial $P_{ep}^{R}(\Omega_c)$")

    if delta_meV > 0.0:
        ax.axvline(
            2.0 * delta_meV,
            linewidth=1.0,
            linestyle=":",
            label=rf"$2\Delta={2.0 * delta_meV:.2f}$ meV",
        )

    ax.axhline(0.0, linewidth=0.8)
    ax.set_title(
        "Low-energy gap-scale power band\n"
        rf"$T_e={result.Te_K:.2f}$ K, $\Delta={delta_meV:.3f}$ meV"
    )
    ax.set_xlabel(r"upper phonon energy $\Omega_c$ [meV]")
    ax.set_ylabel(r"partial power density [W m$^{-3}$]")
    ax.grid(True, linewidth=0.25, alpha=0.35)
    ax.legend(frameon=True, fontsize=8)

    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output