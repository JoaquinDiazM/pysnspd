"""Diagnostic plots for PRE-run projected power and energy catalogues.

The power table is the runtime-oriented reduction of the QP--phonon phase-space
catalogue. These figures are meant to answer three sanity questions before OE6:

1. Do the projected powers have the expected antisymmetric sign with Te and Tph?
2. Are the scattering/recombination channels finite and comparable on reasonable scales?
3. Do the electronic/phononic energy tables and transport coefficients behave
   smoothly enough for runtime interpolation?
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
import numpy as np
from scipy.constants import Boltzmann, hbar
from scipy.special import zeta

from pysnspd.plotting.style import (
    THESIS_DOUBLE_FIGSIZE,
    THESIS_DPI,
    THESIS_WIDTH_IN,
    apply_thesis_style,
)

apply_thesis_style()

MEV_J = 1.602176634e-22


@dataclass(frozen=True)
class DebyeReferenceParameters:
    """Explicit parameters for the Debye/Vodolazov comparison curves.

    ``N0_J_m3`` follows the single-spin convention used by the Simon and
    Vodolazov equations quoted in the thesis.  The reference is deliberately
    separate from the production spectral catalogue.
    """

    N0_J_m3: float
    ion_density_m3: float
    Tc_K: float
    omega_D_J: float
    lambda_ep: float
    tau0_s: float
    normal_dos_spin_convention: str
    lambda_provenance: dict[str, Any]
    runtime_catalog_spectrum: dict[str, Any]
    sources: dict[str, str]

    def manifest_dict(self) -> dict[str, Any]:
        return {
            "enabled": True,
            "model": "normal-state Debye/Vodolazov limiting reference",
            "N0_J_m3": float(self.N0_J_m3),
            "ion_density_m3": float(self.ion_density_m3),
            "Tc_K": float(self.Tc_K),
            "omega_D_J": float(self.omega_D_J),
            "omega_D_meV": float(self.omega_D_J / MEV_J),
            "lambda_ep": float(self.lambda_ep),
            "tau0_s": float(self.tau0_s),
            "normal_dos_spin_convention": self.normal_dos_spin_convention,
            "lambda_provenance": dict(self.lambda_provenance),
            "runtime_catalog_spectrum": dict(self.runtime_catalog_spectrum),
            "sources": dict(self.sources),
            "energy_formula": "u_ph_D=(3*pi^4/5)*Ni*kB*T*(kB*T/OmegaD)^3",
            "heat_capacity_formula": "C_ph_D=(12*pi^4/5)*Ni*kB*(kB*T/OmegaD)^3",
            "power_formula": "P_D=96*zeta(5)*N0*kB^2*(Te^5-Tph^5)/(tau0*Tc^3)",
        }


def resolve_debye_reference_parameters(
    catalog: "PowerTablePlotCatalog",
    *,
    ion_density_m3: float,
    Tc_K: float,
    omega_D_J: float,
    normal_dos_spin_convention: str,
    lambda_ep: float,
    lambda_provenance: dict[str, Any],
) -> DebyeReferenceParameters:
    """Recover every Debye-reference constant from declared stored sources.

    The routine raises instead of silently substituting a literature value for
    a missing production constant.  ``omega_D_J`` is the one intentionally
    external comparison parameter and must be supplied by the caller.
    """

    convention = str(normal_dos_spin_convention).strip().lower().replace("-", "_")
    if convention != "single_spin":
        raise ValueError(
            "The Vodolazov reference requires the single-spin N(0) convention; "
            f"received {normal_dos_spin_convention!r}."
        )
    N0_J_m3 = _metadata_float_recursive(catalog.metadata, "N0_J_m3")
    if not np.isfinite(N0_J_m3) or N0_J_m3 <= 0.0:
        raise ValueError("power-table metadata do not contain a positive N0_J_m3")
    if not np.isfinite(ion_density_m3) or ion_density_m3 <= 0.0:
        raise ValueError("ion_density_m3 must be recovered explicitly from the material configuration")
    if not np.isfinite(Tc_K) or Tc_K <= 0.0:
        raise ValueError("Tc_K must be recovered explicitly from the material configuration")
    if not np.isfinite(omega_D_J) or omega_D_J <= 0.0:
        raise ValueError("omega_D_J must be finite and positive")

    if not np.isfinite(lambda_ep) or lambda_ep <= 0.0:
        raise ValueError("lambda_ep from the complete Eliashberg source must be positive")
    required_provenance = {
        "source_path",
        "sha256",
        "n_points",
        "frequency_min_THz",
        "frequency_max_THz",
        "definition",
    }
    missing_provenance = sorted(required_provenance.difference(lambda_provenance))
    if missing_provenance:
        raise ValueError(
            "lambda_provenance is missing required fields: "
            + ", ".join(missing_provenance)
        )

    runtime_catalog_spectrum = _runtime_catalog_spectrum_diagnostic(catalog)

    tau0_s = float(
        1.0
        / (
            np.pi
            * lambda_ep
            * (Boltzmann * float(Tc_K) / float(omega_D_J)) ** 2
            * (Boltzmann * float(Tc_K) / hbar)
        )
    )
    return DebyeReferenceParameters(
        N0_J_m3=float(N0_J_m3),
        ion_density_m3=float(ion_density_m3),
        Tc_K=float(Tc_K),
        omega_D_J=float(omega_D_J),
        lambda_ep=lambda_ep,
        tau0_s=tau0_s,
        normal_dos_spin_convention="single_spin",
        lambda_provenance=dict(lambda_provenance),
        runtime_catalog_spectrum=runtime_catalog_spectrum,
        sources={
            "N0_J_m3": "power_table_catalog.npz metadata",
            "normal_dos_spin_convention": "declared single-spin catalogue contract",
            "ion_density_m3": "validated project material configuration",
            "Tc_K": "validated project material configuration",
            "lambda_ep": "complete Simon Eliashberg DAT loaded by load_simon_eliashberg_dat",
            "omega_D_J": "explicit Vodolazov comparison cutoff supplied to E1",
            "tau0_s": "Annex B.5 mapping from lambda_ep, Omega_D, and Tc",
        },
    )


def _runtime_catalog_spectrum_diagnostic(
    catalog: "PowerTablePlotCatalog",
) -> dict[str, Any]:
    """Describe the truncated runtime quadrature grid without using it for lambda."""

    omega = np.asarray(catalog.omega_values_J, dtype=float).reshape(-1)
    alpha2F = np.asarray(catalog.alpha2F, dtype=float).reshape(-1)
    count = min(omega.size, alpha2F.size)
    valid = (
        np.isfinite(omega[:count])
        & np.isfinite(alpha2F[:count])
        & (omega[:count] > 0.0)
    )
    diagnostic_lambda = float("nan")
    if np.count_nonzero(valid) >= 2:
        order = np.argsort(omega[:count][valid])
        omega_valid = omega[:count][valid][order]
        alpha_valid = alpha2F[:count][valid][order]
        diagnostic_lambda = float(
            2.0 * np.trapezoid(alpha_valid / omega_valid, omega_valid)
        )
    finite_omega = omega[np.isfinite(omega)]
    return {
        "role": "truncated runtime quadrature grid; never authoritative for lambda_ep",
        "n_points": int(count),
        "energy_min_J": (
            float(np.nanmin(finite_omega)) if finite_omega.size else None
        ),
        "energy_max_J": (
            float(np.nanmax(finite_omega)) if finite_omega.size else None
        ),
        "lambda_if_reintegrated_diagnostic_only": (
            diagnostic_lambda if np.isfinite(diagnostic_lambda) else None
        ),
    }


@dataclass(frozen=True)
class PowerTablePlotCatalog:
    """Small plotting-facing view of ``power_table_catalog.npz``."""

    Te_values_K: np.ndarray
    Tph_values_K: np.ndarray
    delta_values_J: np.ndarray
    q_values_m_inv: np.ndarray
    P_S_W_m3: np.ndarray
    P_R_W_m3: np.ndarray
    P_total_W_m3: np.ndarray
    u_e_J_m3: np.ndarray
    C_e_J_m3_K: np.ndarray
    kappa_s_W_m_K: np.ndarray
    u_ph_J_m3: np.ndarray
    C_ph_J_m3_K: np.ndarray
    u_ph_weighted_J: np.ndarray
    C_ph_weighted_J_K: np.ndarray
    omega_values_J: np.ndarray
    alpha2F: np.ndarray
    phdos_states_per_THz: np.ndarray
    metadata: dict[str, Any]


def write_power_table_diagnostic_plots(
    *,
    power_table_npz: str | Path,
    output_dir: str | Path,
    dpi: int = THESIS_DPI,
) -> dict[str, str]:
    """Write diagnostic plots for a PRE-run projected power table."""
    cat = load_power_table_plot_catalog(power_table_npz)
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    paths = {
        "power_channels_Te_Tph_maps_png": plot_power_channels_Te_Tph_maps(
            cat,
            out / "power_channels_Te_Tph_maps.png",
            dpi=dpi,
        ),
        "power_total_Delta_q_maps_png": plot_power_total_Delta_q_maps(
            cat,
            out / "power_total_Delta_q_maps.png",
            dpi=dpi,
        ),
        "power_total_Te_curves_png": plot_power_total_Te_curves(
            cat,
            out / "power_total_Te_curves.png",
            dpi=dpi,
        ),
        "energy_heat_capacity_curves_png": plot_energy_heat_capacity_curves(
            cat,
            out / "energy_heat_capacity_curves.png",
            dpi=dpi,
        ),
        "electronic_thermal_conductivity_curves_png": plot_electronic_thermal_conductivity_curves(
            cat,
            out / "electronic_thermal_conductivity_curves.png",
            dpi=dpi,
        ),
        "power_equal_temperature_residual_png": plot_equal_temperature_residual(
            cat,
            out / "power_equal_temperature_residual.png",
            dpi=dpi,
        ),
    }
    return {key: str(value) for key, value in paths.items()}


def load_power_table_plot_catalog(path: str | Path) -> PowerTablePlotCatalog:
    """Load a plotting-facing view of ``power_table_catalog.npz``."""
    with np.load(Path(path), allow_pickle=True) as data:
        metadata = _metadata_from_npz(data)
        Te_values = np.asarray(data["Te_values_K"], dtype=float)
        Tph_values = np.asarray(data["Tph_values_K"], dtype=float)
        delta_values = np.asarray(data["delta_values_J"], dtype=float)
        return PowerTablePlotCatalog(
            Te_values_K=Te_values,
            Tph_values_K=Tph_values,
            delta_values_J=delta_values,
            q_values_m_inv=np.asarray(data["q_values_m_inv"], dtype=float),
            P_S_W_m3=np.asarray(data["P_S_W_m3"], dtype=float),
            P_R_W_m3=np.asarray(data["P_R_W_m3"], dtype=float),
            P_total_W_m3=np.asarray(data["P_total_W_m3"], dtype=float),
            u_e_J_m3=np.asarray(data["u_e_J_m3"], dtype=float),
            C_e_J_m3_K=np.asarray(data["C_e_J_m3_K"], dtype=float),
            kappa_s_W_m_K=np.asarray(
                data.get("kappa_s_W_m_K", np.zeros((Te_values.size, delta_values.size), dtype=float)),
                dtype=float,
            ),
            u_ph_J_m3=np.asarray(
                data.get("u_ph_J_m3", data.get("u_ph_weighted_J", np.array([], dtype=float))),
                dtype=float,
            ),
            C_ph_J_m3_K=np.asarray(
                data.get("C_ph_J_m3_K", data.get("C_ph_weighted_J_K", np.array([], dtype=float))),
                dtype=float,
            ),
            u_ph_weighted_J=np.asarray(data.get("u_ph_weighted_J", np.array([], dtype=float)), dtype=float),
            C_ph_weighted_J_K=np.asarray(data.get("C_ph_weighted_J_K", np.array([], dtype=float)), dtype=float),
            omega_values_J=np.asarray(data.get("omega_values_J", np.array([], dtype=float)), dtype=float),
            alpha2F=np.asarray(data.get("alpha2F", np.array([], dtype=float)), dtype=float),
            phdos_states_per_THz=np.asarray(data.get("phdos_states_per_THz", np.array([], dtype=float)), dtype=float),
            metadata=metadata,
        )


def plot_power_channels_Te_Tph_maps(
    catalog: PowerTablePlotCatalog,
    output_path: str | Path,
    *,
    dpi: int = THESIS_DPI,
) -> Path:
    """Plot scattering, recombination and total powers over the (Te,Tph) plane.

    The slice is taken at the largest tabulated gap and q=0. This is the cleanest
    superconducting reference state for checking signs and relative channel size.
    """
    apply_thesis_style()
    output = _prepare_output(output_path)
    i_delta = int(np.nanargmax(catalog.delta_values_J))
    i_q = _nearest_index(catalog.q_values_m_inv, 0.0)

    channels = [
        (catalog.P_S_W_m3[:, :, i_delta, i_q], r"$P_{e\text{-}ph}^{S}$ scattering"),
        (catalog.P_R_W_m3[:, :, i_delta, i_q], r"$P_{e\text{-}ph}^{R}$ recombination"),
        (catalog.P_total_W_m3[:, :, i_delta, i_q], r"$P_{e\text{-}ph}$ total"),
    ]
    vmax = _robust_symmetric_vmax([arr for arr, _ in channels])
    norm = _symmetric_log_norm(vmax)
    extent = _imshow_extent(catalog.Tph_values_K, catalog.Te_values_K)

    fig, axes = plt.subplots(1, 3, figsize=THESIS_DOUBLE_FIGSIZE, constrained_layout=True)
    for ax, (arr, title) in zip(axes, channels):
        im = ax.imshow(
            arr,
            origin="lower",
            aspect="auto",
            interpolation="bilinear",
            extent=extent,
            cmap="coolwarm",
            norm=norm,
        )
        ax.plot(catalog.Tph_values_K, catalog.Tph_values_K, color="black", linewidth=0.75, alpha=0.75)
        ax.set_title(title)
        ax.set_xlabel(r"$T_{ph}$ [K]")
        ax.set_ylabel(r"$T_e$ [K]")
        ax.grid(False)
    cbar = fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.92)
    cbar.set_label(r"Power density [W m$^{-3}$]; positive: electrons $\rightarrow$ phonons")
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output


def plot_power_total_Delta_q_maps(
    catalog: PowerTablePlotCatalog,
    output_path: str | Path,
    *,
    dpi: int = THESIS_DPI,
) -> Path:
    """Plot total projected power over (Delta,q) for four thermal states."""
    output = _prepare_output(output_path)
    Tb = float(np.nanmin(catalog.Tph_values_K))
    Tc_like = _metadata_float(catalog.metadata, "Tc_K")
    if not np.isfinite(Tc_like):
        Tc_like = 8.65
    Te_targets = [Tb, min(Tc_like, float(np.nanmax(catalog.Te_values_K))), 2.0 * Tc_like, float(np.nanmax(catalog.Te_values_K))]
    Te_targets = [float(np.clip(v, np.nanmin(catalog.Te_values_K), np.nanmax(catalog.Te_values_K))) for v in Te_targets]
    Tph_targets = [Tb, Tb, Tb, Tb]

    slices: list[tuple[np.ndarray, str]] = []
    for Te_target, Tph_target in zip(Te_targets, Tph_targets):
        iT = _nearest_index(catalog.Te_values_K, Te_target)
        iP = _nearest_index(catalog.Tph_values_K, Tph_target)
        arr = catalog.P_total_W_m3[iT, iP, :, :]
        label = rf"$T_e={catalog.Te_values_K[iT]:.2f}$ [K], $T_{{ph}}={catalog.Tph_values_K[iP]:.2f}$ [K]"
        slices.append((arr, label))

    vmax = _robust_symmetric_vmax([arr for arr, _ in slices])
    norm = _symmetric_log_norm(vmax)
    q_1e7 = catalog.q_values_m_inv / 1.0e7
    delta_meV = _joule_to_mev(catalog.delta_values_J)
    extent = _imshow_extent(q_1e7, delta_meV)

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(THESIS_WIDTH_IN, 6.0),
        constrained_layout=True,
    )
    for ax, (arr, label) in zip(axes.ravel(), slices):
        im = ax.imshow(
            arr,
            origin="lower",
            aspect="auto",
            interpolation="bilinear",
            extent=extent,
            cmap="coolwarm",
            norm=norm,
        )
        ax.set_title(label)
        ax.set_xlabel(r"$q$ [$10^7$ m$^{-1}$]")
        ax.set_ylabel(r"$|\Delta|$ [meV]")
        ax.grid(False)
    cbar = fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.92)
    cbar.set_label(r"$P_S+P_R$ [W m$^{-3}$]")
    fig.suptitle(r"Projected total electron--phonon power over the $(|\Delta|,q)$ state grid")
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output


def plot_power_total_Te_curves(
    catalog: PowerTablePlotCatalog,
    output_path: str | Path,
    *,
    debye_reference: DebyeReferenceParameters | None = None,
    dpi: int = THESIS_DPI,
) -> Path:
    """Plot production power slices and one optional Debye limiting reference."""
    apply_thesis_style()
    output = _prepare_output(output_path)
    iTph = _nearest_index(catalog.Tph_values_K, float(np.nanmin(catalog.Tph_values_K)))
    states = _representative_state_indices(catalog)

    fig, ax = plt.subplots(figsize=(THESIS_WIDTH_IN, 3.15))
    positive_values: list[np.ndarray] = []
    for label, i_delta, i_q in states:
        y = catalog.P_total_W_m3[:, iTph, i_delta, i_q]
        mask = np.isfinite(y) & (y > 0.0)
        if np.any(mask):
            positive_values.append(y[mask])
            ax.plot(
                catalog.Te_values_K[mask],
                y[mask],
                marker=".",
                markersize=2.2,
                linewidth=1.0,
                label=label,
            )
    T_bath_K = float(catalog.Tph_values_K[iTph])
    if debye_reference is not None:
        reference_power = _debye_power_density(
            catalog.Te_values_K,
            T_bath_K,
            debye_reference,
        )
        mask = np.isfinite(reference_power) & (reference_power > 0.0)
        if np.any(mask):
            positive_values.append(reference_power[mask])
            ax.plot(
                catalog.Te_values_K[mask],
                reference_power[mask],
                color="black",
                linestyle=":",
                linewidth=1.35,
                label="Debye/Vodolazov reference",
            )
    Tc_K = _critical_temperature_K(catalog)
    ax.axvline(
        T_bath_K,
        color="0.20",
        linewidth=0.9,
        linestyle=":",
        label=rf"$T_b={T_bath_K:.2f}$ [K]",
    )
    if np.isfinite(Tc_K):
        ax.axvline(
            Tc_K,
            color="0.35",
            linewidth=0.9,
            linestyle="--",
            label=rf"$T_c={Tc_K:.2f}$ [K]",
        )
    if positive_values:
        positive = np.concatenate(positive_values)
        ax.set_yscale("symlog", linthresh=1.0e11)
        ax.set_ylim(0.0, 1.15 * float(np.nanmax(positive)))
    ax.set_xlabel(r"$T_e$ [K]")
    ax.set_ylabel(r"$P_S+P_R$ [W m$^{-3}$]")
    T_min_K = float(np.nanmin(catalog.Te_values_K))
    T_max_K = float(np.nanmax(catalog.Te_values_K))
    T_margin_K = 0.02 * max(T_max_K - T_min_K, 1.0)
    ax.set_xlim(max(0.0, T_min_K - T_margin_K), T_max_K)
    ax.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, 1.02),
        borderaxespad=0.0,
        ncol=3,
    )
    ax.grid(True, linewidth=0.35, alpha=0.28)
    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output


def plot_energy_heat_capacity_curves(
    catalog: PowerTablePlotCatalog,
    output_path: str | Path,
    *,
    debye_reference: DebyeReferenceParameters | None = None,
    dpi: int = THESIS_DPI,
) -> Path:
    r"""Plot electronic/phononic energy and heat capacity curves.

    Electronic curves are shown for representative $(|\Delta|, q)$ states, while
    the phonon subsystem is overplotted as a single volume-normalized reference
    curve versus $T_{ph}$.  The horizontal axis is therefore labelled as a
    generic temperature variable rather than purely $T_e$.
    """
    apply_thesis_style()
    output = _prepare_output(output_path)
    states = _representative_state_indices(catalog)

    fig, (ax_u, ax_c) = plt.subplots(1, 2, figsize=THESIS_DOUBLE_FIGSIZE, constrained_layout=True)
    for label, i_delta, i_q in states:
        ax_u.plot(
            catalog.Te_values_K,
            catalog.u_e_J_m3[:, i_delta, i_q],
            marker=".",
            markersize=2.0,
            linewidth=1.0,
            label=label,
        )
        ax_c.plot(
            catalog.Te_values_K,
            catalog.C_e_J_m3_K[:, i_delta, i_q],
            marker=".",
            markersize=2.0,
            linewidth=1.0,
            label=label,
        )

    if catalog.u_ph_J_m3.size:
        ax_u.plot(
            catalog.Tph_values_K,
            catalog.u_ph_J_m3,
            color="black",
            linestyle="--",
            linewidth=1.35,
            label=r"Phonons: $u_{ph}(T_{ph})$",
        )
    if catalog.C_ph_J_m3_K.size:
        ax_c.plot(
            catalog.Tph_values_K,
            catalog.C_ph_J_m3_K,
            color="black",
            linestyle="--",
            linewidth=1.35,
            label=r"Phonons: $C_{ph}(T_{ph})$",
        )
    if debye_reference is not None:
        u_debye, C_debye = _debye_phonon_storage(
            catalog.Tph_values_K,
            debye_reference,
        )
        ax_u.plot(
            catalog.Tph_values_K,
            u_debye,
            color="black",
            linestyle=":",
            linewidth=1.35,
            label=r"Debye reference: $u_{ph}^{D}$",
        )
        ax_c.plot(
            catalog.Tph_values_K,
            C_debye,
            color="black",
            linestyle=":",
            linewidth=1.35,
            label=r"Debye reference: $C_{ph}^{D}$",
        )

    ax_u.set_title(r"Energy densities")
    ax_u.set_xlabel(r"Temperature [$T_e$ or $T_{ph}$] [K]")
    ax_u.set_ylabel(r"Energy density [J m$^{-3}$]")
    ax_u.grid(True, linewidth=0.35, alpha=0.28)
    # ax_u.legend(loc="best", fontsize=7.0)
    ax_u.legend(loc="best")

    ax_c.set_title(r"Heat capacities")
    ax_c.set_xlabel(r"Temperature [$T_e$ or $T_{ph}$] [K]")
    ax_c.set_ylabel(r"Heat capacity [J m$^{-3}$ K$^{-1}$]")
    finite_c = np.abs(catalog.C_e_J_m3_K[np.isfinite(catalog.C_e_J_m3_K)])
    if catalog.C_ph_J_m3_K.size:
        finite_c = np.concatenate([finite_c, np.abs(catalog.C_ph_J_m3_K[np.isfinite(catalog.C_ph_J_m3_K)])]) if finite_c.size else np.abs(catalog.C_ph_J_m3_K[np.isfinite(catalog.C_ph_J_m3_K)])
    max_c = float(np.nanmax(finite_c)) if finite_c.size else 1.0
    ax_c.set_yscale("symlog", linthresh=max(1.0e-3, 1.0e-4 * max_c))
    ax_c.grid(True, linewidth=0.35, alpha=0.28)
    # ax_c.legend(loc="best", fontsize=7.0)
    ax_c.legend(loc="best")

    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output


def plot_electronic_thermal_conductivity_curves(
    catalog: PowerTablePlotCatalog,
    output_path: str | Path,
    *,
    dpi: int = THESIS_DPI,
) -> Path:
    """Plot the Bardeen/Allmaras superconducting electronic thermal conductivity."""
    apply_thesis_style()
    output = _prepare_output(output_path)

    fig, ax = plt.subplots(figsize=(THESIS_WIDTH_IN, 3.0))
    delta_states = _representative_delta_indices(catalog)
    Tc_K = _critical_temperature_K(catalog)
    T_max_K = float(np.nanmax(catalog.Te_values_K))
    for prefix, i_delta in delta_states:
        y = catalog.kappa_s_W_m_K[:, i_delta]
        kappa_at_tc = _interpolate_finite(catalog.Te_values_K, y, Tc_K)
        kappa_at_tmax = _interpolate_finite(catalog.Te_values_K, y, T_max_K)
        # label = (
        #     _delta_label(catalog, prefix, i_delta)
        #     + "\n"
        #     + rf"$\kappa_s(T_c)={_format_compact_value(kappa_at_tc)}$, "
        #     + rf"$\kappa_s(T_{{\max}})={_format_compact_value(kappa_at_tmax)}$ [W m$^{{-1}}$ K$^{{-1}}$]"
        # )
        label = prefix
        ax.plot(
            catalog.Te_values_K,
            y,
            marker=".",
            markersize=2.0,
            linewidth=1.0,
            label=label,
        )

    ax.set_xlabel(r"$T_e$ [K]")
    ax.set_ylabel(r"$\kappa_s$ [W m$^{-1}$ K$^{-1}$]")
    positive = catalog.kappa_s_W_m_K[np.isfinite(catalog.kappa_s_W_m_K) & (catalog.kappa_s_W_m_K > 0.0)]
    if positive.size and float(np.nanmax(positive) / max(np.nanmin(positive), 1.0e-300)) > 50.0:
        ax.set_yscale("log")
    ax.grid(True, linewidth=0.35, alpha=0.28)
    # ax.legend(loc="best", fontsize=6.8, labelspacing=0.55)
    ax.legend(loc="best", labelspacing=0.55)
    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output


def plot_equal_temperature_residual(
    catalog: PowerTablePlotCatalog,
    output_path: str | Path,
    *,
    dpi: int = THESIS_DPI,
) -> Path:
    """Plot detailed-balance residual max |P_total(Te=Tph)| over the state grid."""
    output = _prepare_output(output_path)
    common_T, iTe, iTph = _matched_temperature_indices(catalog.Te_values_K, catalog.Tph_values_K)
    residual = np.empty(common_T.size, dtype=float)
    for k, (it, ip) in enumerate(zip(iTe, iTph)):
        residual[k] = float(np.nanmax(np.abs(catalog.P_total_W_m3[it, ip, :, :])))

    fig, ax = plt.subplots(figsize=THESIS_DOUBLE_FIGSIZE)
    ax.plot(common_T, residual, marker=".", markersize=3.0, linewidth=1.0)
    ax.set_title(r"Detailed-balance check: max $|P_S+P_R|$ at $T_e=T_{ph}$")
    ax.set_xlabel("temperature [K]")
    ax.set_ylabel(r"max state residual [W m$^{-3}$]")
    if np.nanmax(residual) > 0.0:
        ax.set_yscale("symlog", linthresh=max(1.0, 1.0e-6 * float(np.nanmax(residual))))
    ax.grid(True, linewidth=0.35, alpha=0.28)
    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output


def _representative_state_indices(catalog: PowerTablePlotCatalog) -> list[tuple[str, int, int]]:
    i_delta0 = _nearest_index(catalog.delta_values_J, 0.0)
    i_delta_max = int(np.nanargmax(catalog.delta_values_J))
    i_q0 = _nearest_index(catalog.q_values_m_inv, 0.0)
    return [
        ("Normal-like", i_delta0, i_q0),
        ("SC, q=0", i_delta_max, i_q0),
    ]


def _representative_delta_indices(catalog: PowerTablePlotCatalog) -> list[tuple[str, int]]:
    i_delta0 = _nearest_index(catalog.delta_values_J, 0.0)
    i_delta_half = _nearest_index(catalog.delta_values_J, 0.5 * float(np.nanmax(catalog.delta_values_J)))
    i_delta_max = int(np.nanargmax(catalog.delta_values_J))
    return [
        ("Normal-like", i_delta0),
        ("Intermediate gap", i_delta_half),
        ("Maximum gap", i_delta_max),
    ]


def _matched_temperature_indices(Te: np.ndarray, Tph: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    Te = np.asarray(Te, dtype=float)
    Tph = np.asarray(Tph, dtype=float)
    n = min(Te.size, Tph.size)
    if n == 0:
        return np.array([], dtype=float), np.array([], dtype=np.int64), np.array([], dtype=np.int64)
    if Te.size == Tph.size and np.allclose(Te, Tph, rtol=1.0e-12, atol=1.0e-12):
        idx = np.arange(Te.size, dtype=np.int64)
        return Te.copy(), idx, idx.copy()
    iTe = np.arange(Te.size, dtype=np.int64)
    iTph = np.asarray([_nearest_index(Tph, value) for value in Te], dtype=np.int64)
    common = 0.5 * (Te + Tph[iTph])
    return common, iTe, iTph


def _metadata_from_npz(data: Any) -> dict[str, Any]:
    if "metadata" not in data.files:
        return {}
    raw = data["metadata"]
    try:
        value = raw.item()
    except Exception:
        value = raw
    return value if isinstance(value, dict) else {}


def _nearest_index(values: np.ndarray, target: float) -> int:
    arr = np.asarray(values, dtype=float)
    return int(np.nanargmin(np.abs(arr - float(target))))


def _robust_symmetric_vmax(arrays: list[np.ndarray]) -> float:
    finite_parts = []
    for arr in arrays:
        a = np.asarray(arr, dtype=float)
        finite = a[np.isfinite(a)]
        if finite.size:
            finite_parts.append(np.abs(finite))
    if not finite_parts:
        return 1.0
    all_abs = np.concatenate(finite_parts)
    vmax = float(np.nanpercentile(all_abs, 99.5))
    return max(vmax, float(np.nanmax(all_abs)), 1.0)


def _symmetric_log_norm(vmax: float) -> SymLogNorm:
    vmax = max(float(vmax), 1.0)
    return SymLogNorm(linthresh=max(1.0e-6 * vmax, 1.0), vmin=-vmax, vmax=vmax, base=10.0)


def _imshow_extent(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float, float]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size <= 1:
        dx = max(abs(float(x[0])) * 0.05, 0.5) if x.size else 0.5
        xmin, xmax = float(x[0] - dx), float(x[0] + dx)
    else:
        dx = np.diff(x)
        xmin = float(x[0] - 0.5 * dx[0])
        xmax = float(x[-1] + 0.5 * dx[-1])
    if y.size <= 1:
        dy = max(abs(float(y[0])) * 0.05, 0.5) if y.size else 0.5
        ymin, ymax = float(y[0] - dy), float(y[0] + dy)
    else:
        dy = np.diff(y)
        ymin = float(y[0] - 0.5 * dy[0])
        ymax = float(y[-1] + 0.5 * dy[-1])
    return xmin, xmax, ymin, ymax


def _joule_to_mev(values_J: np.ndarray | float) -> np.ndarray | float:
    return np.asarray(values_J, dtype=float) / MEV_J


def _metadata_float(metadata: dict[str, Any], key: str) -> float:
    if not isinstance(metadata, dict) or key not in metadata:
        return float("nan")
    try:
        return float(metadata[key])
    except Exception:
        return float("nan")


def _metadata_float_recursive(metadata: dict[str, Any], key: str) -> float:
    """Find a scalar in nested catalogue metadata without guessing aliases."""

    if not isinstance(metadata, dict):
        return float("nan")
    if key in metadata:
        try:
            return float(metadata[key])
        except (TypeError, ValueError):
            return float("nan")
    for value in metadata.values():
        if isinstance(value, dict):
            found = _metadata_float_recursive(value, key)
            if np.isfinite(found):
                return found
    return float("nan")


def _debye_phonon_storage(
    temperature_K: np.ndarray,
    reference: DebyeReferenceParameters,
) -> tuple[np.ndarray, np.ndarray]:
    temperature = np.asarray(temperature_K, dtype=float)
    reduced = Boltzmann * temperature / float(reference.omega_D_J)
    energy = (
        3.0
        * np.pi**4
        / 5.0
        * float(reference.ion_density_m3)
        * Boltzmann
        * temperature
        * reduced**3
    )
    heat_capacity = (
        12.0
        * np.pi**4
        / 5.0
        * float(reference.ion_density_m3)
        * Boltzmann
        * reduced**3
    )
    return np.asarray(energy, dtype=float), np.asarray(heat_capacity, dtype=float)


def _debye_power_density(
    electron_temperature_K: np.ndarray,
    phonon_temperature_K: float,
    reference: DebyeReferenceParameters,
) -> np.ndarray:
    electron_temperature = np.asarray(electron_temperature_K, dtype=float)
    coefficient = (
        96.0
        * float(zeta(5.0, 1.0))
        * float(reference.N0_J_m3)
        * Boltzmann**2
        / (float(reference.tau0_s) * float(reference.Tc_K) ** 3)
    )
    return coefficient * (electron_temperature**5 - float(phonon_temperature_K) ** 5)


def _critical_temperature_K(catalog: PowerTablePlotCatalog) -> float:
    Tc_K = _metadata_float(catalog.metadata, "Tc_K")
    if np.isfinite(Tc_K) and Tc_K > 0.0:
        return Tc_K
    delta0_J = _metadata_float(catalog.metadata, "delta0_J")
    if np.isfinite(delta0_J) and delta0_J > 0.0:
        return float(delta0_J / (1.764 * Boltzmann))
    return float("nan")


def _interpolate_finite(x: np.ndarray, y: np.ndarray, target: float) -> float:
    x_values = np.asarray(x, dtype=float)
    y_values = np.asarray(y, dtype=float)
    mask = np.isfinite(x_values) & np.isfinite(y_values)
    if not np.any(mask) or not np.isfinite(target):
        return float("nan")
    order = np.argsort(x_values[mask])
    return float(np.interp(float(target), x_values[mask][order], y_values[mask][order]))


def _prepare_output(path: str | Path) -> Path:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    return output


__all__ = [
    "DebyeReferenceParameters",
    "PowerTablePlotCatalog",
    "load_power_table_plot_catalog",
    "resolve_debye_reference_parameters",
    "write_power_table_diagnostic_plots",
    "plot_power_channels_Te_Tph_maps",
    "plot_power_total_Delta_q_maps",
    "plot_power_total_Te_curves",
    "plot_energy_heat_capacity_curves",
    "plot_electronic_thermal_conductivity_curves",
    "plot_equal_temperature_residual",
]
