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

MEV_J = 1.602176634e-22


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
    dpi: int = 480,
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
    dpi: int = 480,
) -> Path:
    """Plot scattering, recombination and total powers over the (Te,Tph) plane.

    The slice is taken at the largest tabulated gap and q=0. This is the cleanest
    superconducting reference state for checking signs and relative channel size.
    """
    output = _prepare_output(output_path)
    i_delta = int(np.nanargmax(catalog.delta_values_J))
    i_q = _nearest_index(catalog.q_values_m_inv, 0.0)
    delta_meV = _joule_to_mev(catalog.delta_values_J[i_delta])
    q_1e7 = catalog.q_values_m_inv[i_q] / 1.0e7

    channels = [
        (catalog.P_S_W_m3[:, :, i_delta, i_q], r"$P_S$ scattering"),
        (catalog.P_R_W_m3[:, :, i_delta, i_q], r"$P_R$ recombination"),
        (catalog.P_total_W_m3[:, :, i_delta, i_q], r"$P_S+P_R$ total"),
    ]
    vmax = _robust_symmetric_vmax([arr for arr, _ in channels])
    norm = _symmetric_log_norm(vmax)
    extent = _imshow_extent(catalog.Tph_values_K, catalog.Te_values_K)

    fig, axes = plt.subplots(1, 3, figsize=(12.3, 3.75), constrained_layout=True)
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
    cbar.set_label(r"power density [W m$^{-3}$], positive: electrons $\rightarrow$ phonons")
    fig.suptitle(rf"Projected powers at $|\Delta|={delta_meV:.3f}$ meV, $q={q_1e7:.2f}\times10^7$ m$^{{-1}}$")
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output



def plot_power_total_Delta_q_maps(
    catalog: PowerTablePlotCatalog,
    output_path: str | Path,
    *,
    dpi: int = 480,
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
        label = rf"$T_e={catalog.Te_values_K[iT]:.2f}$ K, $T_{{ph}}={catalog.Tph_values_K[iP]:.2f}$ K"
        slices.append((arr, label))

    vmax = _robust_symmetric_vmax([arr for arr, _ in slices])
    norm = _symmetric_log_norm(vmax)
    q_1e7 = catalog.q_values_m_inv / 1.0e7
    delta_meV = _joule_to_mev(catalog.delta_values_J)
    extent = _imshow_extent(q_1e7, delta_meV)

    fig, axes = plt.subplots(2, 2, figsize=(9.2, 7.0), constrained_layout=True)
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
    dpi: int = 480,
) -> Path:
    """Plot total power versus Te at bath phonon temperature for representative states."""
    output = _prepare_output(output_path)
    iTph = _nearest_index(catalog.Tph_values_K, float(np.nanmin(catalog.Tph_values_K)))
    states = _representative_state_indices(catalog)

    fig, ax = plt.subplots(figsize=(7.4, 4.45))
    max_abs = 0.0
    for label, i_delta, i_q in states:
        y = catalog.P_total_W_m3[:, iTph, i_delta, i_q]
        max_abs = max(max_abs, float(np.nanmax(np.abs(y)))) if np.isfinite(y).any() else max_abs
        ax.plot(catalog.Te_values_K, y, marker=".", markersize=2.2, linewidth=1.0, label=label)
    ax.axhline(0.0, color="black", linewidth=0.75, alpha=0.75)
    ax.axvline(catalog.Tph_values_K[iTph], color="black", linewidth=0.75, alpha=0.55, linestyle="--")
    ax.set_yscale("symlog", linthresh=max(1.0e-6 * max_abs, 1.0))
    ax.set_title(rf"Projected total power vs $T_e$ at $T_{{ph}}={catalog.Tph_values_K[iTph]:.2f}$ K")
    ax.set_xlabel(r"$T_e$ [K]")
    ax.set_ylabel(r"$P_S+P_R$ [W m$^{-3}$]")
    ax.legend(loc="best", fontsize=7.2)
    ax.grid(True, linewidth=0.35, alpha=0.28)
    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output



def plot_energy_heat_capacity_curves(
    catalog: PowerTablePlotCatalog,
    output_path: str | Path,
    *,
    dpi: int = 480,
) -> Path:
    """Plot electronic/phononic energy and heat capacity curves.

    Electronic curves are shown for representative $(|\Delta|, q)$ states, while
    the phonon subsystem is overplotted as a single volume-normalized reference
    curve versus $T_{ph}$.  The horizontal axis is therefore labelled as a
    generic temperature variable rather than purely $T_e$.
    """
    output = _prepare_output(output_path)
    states = _representative_state_indices(catalog)

    fig, (ax_u, ax_c) = plt.subplots(1, 2, figsize=(11.4, 4.35), constrained_layout=True)
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
            label=r"phonons: $u_{ph}(T_{ph})$",
        )
    if catalog.C_ph_J_m3_K.size:
        ax_c.plot(
            catalog.Tph_values_K,
            catalog.C_ph_J_m3_K,
            color="black",
            linestyle="--",
            linewidth=1.35,
            label=r"phonons: $C_{ph}(T_{ph})$",
        )

    ax_u.set_title(r"Energy densities")
    ax_u.set_xlabel(r"temperature variable [$T_e$ or $T_{ph}$] [K]")
    ax_u.set_ylabel(r"energy density [J m$^{-3}$]")
    ax_u.grid(True, linewidth=0.35, alpha=0.28)
    ax_u.legend(loc="best", fontsize=7.0)

    ax_c.set_title(r"Heat capacities")
    ax_c.set_xlabel(r"temperature variable [$T_e$ or $T_{ph}$] [K]")
    ax_c.set_ylabel(r"heat capacity [J m$^{-3}$ K$^{-1}$]")
    finite_c = np.abs(catalog.C_e_J_m3_K[np.isfinite(catalog.C_e_J_m3_K)])
    if catalog.C_ph_J_m3_K.size:
        finite_c = np.concatenate([finite_c, np.abs(catalog.C_ph_J_m3_K[np.isfinite(catalog.C_ph_J_m3_K)])]) if finite_c.size else np.abs(catalog.C_ph_J_m3_K[np.isfinite(catalog.C_ph_J_m3_K)])
    max_c = float(np.nanmax(finite_c)) if finite_c.size else 1.0
    ax_c.set_yscale("symlog", linthresh=max(1.0e-3, 1.0e-4 * max_c))
    ax_c.grid(True, linewidth=0.35, alpha=0.28)
    ax_c.legend(loc="best", fontsize=7.0)

    fig.suptitle("Runtime energy and heat-capacity tables extracted from the PRE power catalogue")
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output



def plot_electronic_thermal_conductivity_curves(
    catalog: PowerTablePlotCatalog,
    output_path: str | Path,
    *,
    dpi: int = 480,
) -> Path:
    """Plot the Bardeen/Allmaras superconducting electronic thermal conductivity."""
    output = _prepare_output(output_path)

    fig, ax = plt.subplots(figsize=(7.2, 4.35))
    delta_states = _representative_delta_indices(catalog)
    for prefix, i_delta in delta_states:
        y = catalog.kappa_s_W_m_K[:, i_delta]
        ax.plot(
            catalog.Te_values_K,
            y,
            marker=".",
            markersize=2.0,
            linewidth=1.0,
            label=_delta_label(catalog, prefix, i_delta),
        )

    ax.set_title(r"Superconducting electronic thermal conductivity")
    ax.set_xlabel(r"$T_e$ [K]")
    ax.set_ylabel(r"$\kappa_s$ [W m$^{-1}$ K$^{-1}$]")
    positive = catalog.kappa_s_W_m_K[np.isfinite(catalog.kappa_s_W_m_K) & (catalog.kappa_s_W_m_K > 0.0)]
    if positive.size and float(np.nanmax(positive) / max(np.nanmin(positive), 1.0e-300)) > 50.0:
        ax.set_yscale("log")
    ax.grid(True, linewidth=0.35, alpha=0.28)
    ax.legend(loc="best", fontsize=7.2)
    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output



def plot_equal_temperature_residual(
    catalog: PowerTablePlotCatalog,
    output_path: str | Path,
    *,
    dpi: int = 480,
) -> Path:
    """Plot detailed-balance residual max |P_total(Te=Tph)| over the state grid."""
    output = _prepare_output(output_path)
    common_T, iTe, iTph = _matched_temperature_indices(catalog.Te_values_K, catalog.Tph_values_K)
    residual = np.empty(common_T.size, dtype=float)
    for k, (it, ip) in enumerate(zip(iTe, iTph)):
        residual[k] = float(np.nanmax(np.abs(catalog.P_total_W_m3[it, ip, :, :])))

    fig, ax = plt.subplots(figsize=(6.7, 4.15))
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
    i_delta_half = _nearest_index(catalog.delta_values_J, 0.5 * float(np.nanmax(catalog.delta_values_J)))
    i_q0 = _nearest_index(catalog.q_values_m_inv, 0.0)
    i_q_mid = _nearest_index(catalog.q_values_m_inv, 0.5 * float(np.nanmax(catalog.q_values_m_inv)))
    i_q_high = _nearest_index(catalog.q_values_m_inv, 0.85 * float(np.nanmax(catalog.q_values_m_inv)))
    return [
        (_state_label(catalog, "normal-like", i_delta0, i_q0), i_delta0, i_q0),
        (_state_label(catalog, "SC q=0", i_delta_max, i_q0), i_delta_max, i_q0),
        (_state_label(catalog, "SC mid-q", i_delta_max, i_q_mid), i_delta_max, i_q_mid),
        (_state_label(catalog, "reduced gap high-q", i_delta_half, i_q_high), i_delta_half, i_q_high),
    ]



def _representative_delta_indices(catalog: PowerTablePlotCatalog) -> list[tuple[str, int]]:
    i_delta0 = _nearest_index(catalog.delta_values_J, 0.0)
    i_delta_half = _nearest_index(catalog.delta_values_J, 0.5 * float(np.nanmax(catalog.delta_values_J)))
    i_delta_max = int(np.nanargmax(catalog.delta_values_J))
    return [
        ("normal-like", i_delta0),
        ("intermediate gap", i_delta_half),
        ("max gap", i_delta_max),
    ]



def _delta_label(catalog: PowerTablePlotCatalog, prefix: str, i_delta: int) -> str:
    delta_meV = float(_joule_to_mev(catalog.delta_values_J[i_delta]))
    return rf"{prefix}: $|\Delta|={delta_meV:.2f}$ meV"



def _state_label(catalog: PowerTablePlotCatalog, prefix: str, i_delta: int, i_q: int) -> str:
    delta_meV = float(_joule_to_mev(catalog.delta_values_J[i_delta]))
    q_1e7 = float(catalog.q_values_m_inv[i_q] / 1.0e7)
    return rf"{prefix}: $|\Delta|={delta_meV:.2f}$ meV, $q={q_1e7:.2f}$"



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



def _prepare_output(path: str | Path) -> Path:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    return output


__all__ = [
    "PowerTablePlotCatalog",
    "load_power_table_plot_catalog",
    "write_power_table_diagnostic_plots",
    "plot_power_channels_Te_Tph_maps",
    "plot_power_total_Delta_q_maps",
    "plot_power_total_Te_curves",
    "plot_energy_heat_capacity_curves",
    "plot_electronic_thermal_conductivity_curves",
    "plot_equal_temperature_residual",
]
