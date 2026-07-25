#!/usr/bin/env python3
"""Temporary D1 diagnostic for low-amplitude Usadel current regularization.

The pipeline reads an existing PRE catalogue and never modifies it.  It compares
the production interpolation of ``j_s`` with a candidate stiffness
interpolation on ``|Delta|^2`` and with direct Matsubara evaluations.  The
second figure passes those closures through a smooth one-dimensional notch to
expose the phase-source factor that can amplify constitutive errors.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any, Iterable

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import yaml

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from pysnspd.analysis.usadel_low_amplitude_diagnostic import (
    ConstitutiveCurves,
    DiagnosticCurrentCatalog,
    NotchDiagnostic,
    build_notch_diagnostic,
    compare_constitutive_curves,
    load_diagnostic_current_catalog,
)
from pysnspd.config import load_config, validate_config
from pysnspd.io.manager import create_run_layout
from pysnspd.plotting.style import THESIS_DPI, THESIS_DOUBLE_FIGSIZE, apply_thesis_style


DEFAULT_PRE_RUN = (
    "pre_oe6_v3_ultra_L360nm_mesh4p0nm_smooth50_"
    "js81T101D121Q_phase200T31D41Q2400W_power200Tph_01"
)
Q_COLORS = ("#0072B2", "#D55E00", "#009E73")
METHOD_STYLES = {
    "current": ("--", 1.4, "interpolación actual de $j_s$"),
    "stiffness": ("-", 1.6, "candidato $\\kappa(|\\Delta|^2)$"),
    "direct": (":", 1.7, "Matsubara directo"),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Diagnose the R->0 ambiguity in the strict PRE Usadel current table "
            "without changing PRE generation or solver code."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--pre-run-name", default=DEFAULT_PRE_RUN)
    parser.add_argument(
        "--catalog-npz",
        type=Path,
        default=None,
        help="Direct PRE catalogue path; otherwise use raw/<pre-run>/pre/usadel_dos_catalog.npz.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Defaults to plots/<pre-run>/figures/D1_usadel_low_amplitude.",
    )
    parser.add_argument(
        "--temperatures-K",
        nargs="+",
        type=float,
        default=None,
        help="Temperatures to compare; automatic values span the physical PRE range.",
    )
    parser.add_argument(
        "--q-fractions-of-qc",
        nargs="+",
        type=float,
        default=(0.2, 0.5, 0.8),
        help="Phase-gradient values as fractions of the PRE critical q.",
    )
    parser.add_argument("--n-amplitude", type=int, default=180)
    parser.add_argument("--notch-points", type=int, default=801)
    parser.add_argument("--notch-half-width-nm", type=float, default=36.0)
    parser.add_argument("--notch-sigma-nm", type=float, default=8.0)
    parser.add_argument("--dpi", type=int, default=THESIS_DPI)
    parser.add_argument(
        "--formats",
        nargs="+",
        choices=("pdf", "png"),
        default=("pdf", "png"),
        help="Figure formats written for each diagnostic.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    cfg = validate_config(load_config(args.config))
    layout = create_run_layout(cfg, args.pre_run_name)
    raw_pre = Path(layout["raw_pre"])
    catalog_path = (
        args.catalog_npz.expanduser().resolve()
        if args.catalog_npz is not None
        else raw_pre / "usadel_dos_catalog.npz"
    )
    if not catalog_path.exists():
        raise FileNotFoundError(f"Usadel PRE catalogue not found: {catalog_path}")
    figures_dir = (
        args.output_dir.expanduser().resolve()
        if args.output_dir is not None
        else Path(layout["plots_figures"]) / "D1_usadel_low_amplitude"
    )
    figures_dir.mkdir(parents=True, exist_ok=True)

    catalog = load_diagnostic_current_catalog(catalog_path)
    temperatures = _resolve_temperatures(catalog, args.temperatures_K)
    q_values, q_critical = _resolve_q_values(catalog, args.q_fractions_of_qc)
    first = catalog.first_positive_delta_J
    amplitude = np.geomspace(1.0e-4 * first, 1.25 * first, max(80, int(args.n_amplitude)))
    curves = [
        compare_constitutive_curves(catalog, Te_K=T, q_m_inv=q, delta_J=amplitude)
        for T in temperatures
        for q in q_values
    ]

    apply_thesis_style()
    saved: dict[str, list[str]] = {}
    saved["constitutive"] = _save_formats(
        _plot_constitutive_grid(catalog, temperatures, q_values, curves, q_critical),
        figures_dir / "D1_low_amplitude_constitutive_comparison",
        args.formats,
        int(args.dpi),
    )

    notch_T = float(temperatures[len(temperatures) // 2])
    notch_q = float(q_values[len(q_values) // 2])
    notch = build_notch_diagnostic(
        catalog,
        Te_K=notch_T,
        q_m_inv=notch_q,
        n_points=max(201, int(args.notch_points)),
        half_width_nm=float(args.notch_half_width_nm),
        notch_sigma_nm=float(args.notch_sigma_nm),
    )
    saved["notch"] = _save_formats(
        _plot_notch(catalog, notch),
        figures_dir / "D1_notch_phase_source_comparison",
        args.formats,
        int(args.dpi),
    )

    summary = _build_summary(
        args=args,
        catalog_path=catalog_path,
        figures_dir=figures_dir,
        catalog=catalog,
        q_critical=q_critical,
        curves=curves,
        notch=notch,
        saved=saved,
    )
    summary_path = figures_dir / "D1_usadel_low_amplitude_summary.yaml"
    with summary_path.open("w", encoding="utf-8") as stream:
        yaml.safe_dump(summary, stream, sort_keys=False, allow_unicode=True)

    print("D1 low-amplitude Usadel diagnostic")
    print(f" PRE: {args.pre_run_name}")
    print(f" catalogue: {catalog_path}")
    print(f" table shape: {catalog.js_A_m2.shape}")
    print(f" first positive |Delta|/Delta0: {first / catalog.delta0_J:.6g}")
    print(f" q_c: {q_critical * 1.0e-6:.6g} um^-1")
    print()
    for metric in summary["aggregate_metrics"].items():
        print(f" {metric[0]}: {metric[1]:.6g}")
    print()
    for group, paths in saved.items():
        print(f" {group}:")
        for path in paths:
            print(f"   {path}")
    print(f" summary: {summary_path}")
    print("Status: diagnostic only; production code unchanged")
    return 0


def _plot_constitutive_grid(
    catalog: DiagnosticCurrentCatalog,
    temperatures: np.ndarray,
    q_values: np.ndarray,
    curves: list[ConstitutiveCurves],
    q_critical: float,
) -> plt.Figure:
    ncols = len(temperatures)
    fig, axes = plt.subplots(
        2,
        ncols,
        figsize=(THESIS_DOUBLE_FIGSIZE[0], THESIS_DOUBLE_FIGSIZE[1]),
        squeeze=False,
        sharex=True,
    )
    curve_map = {(round(item.Te_K, 12), round(item.q_m_inv, 3)): item for item in curves}
    first = catalog.first_positive_delta_J

    for col, temperature in enumerate(temperatures):
        ax_current = axes[0, col]
        ax_exponent = axes[1, col]
        common_scale = max(
            float(
                np.interp(
                    first,
                    curve_map[(round(float(temperature), 12), round(float(q_value), 3))].delta_J,
                    curve_map[(round(float(temperature), 12), round(float(q_value), 3))].direct_A_m2,
                )
            )
            for q_value in q_values
        )
        common_scale = max(abs(common_scale), np.finfo(float).tiny)
        for q_index, q_value in enumerate(q_values):
            item = curve_map[(round(float(temperature), 12), round(float(q_value), 3))]
            x = item.delta_J / first
            for key, values in (
                ("current", item.current_A_m2),
                ("stiffness", item.stiffness_A_m2),
                ("direct", item.direct_A_m2),
            ):
                linestyle, linewidth, _ = METHOD_STYLES[key]
                ax_current.loglog(
                    x,
                    np.abs(values) / common_scale,
                    color=Q_COLORS[q_index],
                    linestyle=linestyle,
                    linewidth=linewidth,
                )
            for key, values in (
                ("current", item.current_exponent),
                ("stiffness", item.stiffness_exponent),
                ("direct", item.direct_exponent),
            ):
                linestyle, linewidth, _ = METHOD_STYLES[key]
                ax_exponent.semilogx(
                    x,
                    values,
                    color=Q_COLORS[q_index],
                    linestyle=linestyle,
                    linewidth=linewidth,
                )

        ax_current.axvline(1.0, color="0.55", linewidth=0.7)
        ax_exponent.axvline(1.0, color="0.55", linewidth=0.7)
        ax_exponent.axhline(2.0, color="0.35", linewidth=0.7, linestyle=":")
        ax_current.set_title(f"$T_e={temperature:.3g}$ K")
        ax_current.set_ylim(bottom=1.0e-10)
        ax_exponent.set_ylim(0.5, 2.35)
        ax_exponent.set_xlabel(r"$|\Delta|/|\Delta|_1$")
        if col == 0:
            ax_current.set_ylabel(r"$|j_s|/j_s^{\mathrm{dir}}(|\Delta|_1)$")
            ax_exponent.set_ylabel(r"$p=d\ln |j_s|/d\ln |\Delta|$")

    method_handles = [
        Line2D(
            [0],
            [0],
            color="0.15",
            linestyle=style,
            linewidth=width,
            label=label,
        )
        for style, width, label in METHOD_STYLES.values()
    ]
    q_handles = [
        Line2D(
            [0],
            [0],
            color=Q_COLORS[index],
            linewidth=1.7,
            label=f"$q={q * 1.0e-6:.3g}$ $\\mu$m$^{{-1}}$ ({q / q_critical:.2g}$q_c$)",
        )
        for index, q in enumerate(q_values)
    ]
    fig.legend(
        handles=method_handles + q_handles,
        loc="lower center",
        ncol=3,
        bbox_to_anchor=(0.5, 0.005),
        frameon=False,
    )
    fig.suptitle(
        "Diagnóstico D1: ley constitutiva cerca de amplitud nula\n"
        f"$|\\Delta|_1/\\Delta_0={first / catalog.delta0_J:.4g}$; "
        "la línea vertical marca el primer nodo positivo de la PRE",
        y=0.995,
    )
    fig.tight_layout(rect=(0.0, 0.12, 1.0, 0.93))
    return fig


def _plot_notch(catalog: DiagnosticCurrentCatalog, notch: NotchDiagnostic) -> plt.Figure:
    fig, axes = plt.subplots(3, 1, figsize=THESIS_DOUBLE_FIGSIZE, sharex=True)
    x_nm = notch.x_m * 1.0e9
    axes[0].plot(x_nm, notch.delta_J / catalog.first_positive_delta_J, color="#222222")
    axes[0].axhline(1.0, color="0.55", linestyle=":", linewidth=0.8)
    axes[0].set_ylabel(r"$|\Delta|/|\Delta|_1$")

    direct_current_scale = max(float(np.max(np.abs(notch.direct_A_m2))), np.finfo(float).tiny)
    for key, values in (
        ("current", notch.current_A_m2 - notch.gl_A_m2),
        ("stiffness", notch.stiffness_A_m2 - notch.gl_A_m2),
        ("direct", notch.direct_A_m2 - notch.gl_A_m2),
    ):
        style, width, label = METHOD_STYLES[key]
        axes[1].plot(
            x_nm,
            values / direct_current_scale,
            linestyle=style,
            linewidth=width,
            label=label,
        )
    axes[1].set_ylabel(r"$(j_s^{Us}-j_s^{GL})/\max|j_s^{dir}|$")

    direct_source_scale = max(
        float(np.max(np.abs(notch.source_direct_A_m3_J_inv))),
        np.finfo(float).tiny,
    )
    source_series = (
        (
            "current",
            notch.source_current_A_m3_J_inv,
            notch.metrics["current_source_peak_over_direct_peak"],
        ),
        (
            "stiffness",
            notch.source_stiffness_A_m3_J_inv,
            notch.metrics["stiffness_source_peak_over_direct_peak"],
        ),
        ("direct", notch.source_direct_A_m3_J_inv, 1.0),
    )
    for key, values, peak_ratio in source_series:
        style, width, label = METHOD_STYLES[key]
        axes[2].plot(
            x_nm,
            values / direct_source_scale,
            linestyle=style,
            linewidth=width,
            label=f"{label}; pico={peak_ratio:.3g}$\\times$ directo",
        )
    axes[2].set_yscale("symlog", linthresh=0.05)
    axes[2].set_yticks((-100.0, -1.0, 0.0, 1.0, 100.0))
    axes[2].set_ylabel(r"$[\partial_x(j_s^{Us}-j_s^{GL})/|\Delta|]/S_{dir}$")
    axes[2].set_xlabel("$x$ [nm]")
    axes[2].legend(loc="upper right", fontsize=6.5)
    axes[0].set_title(
        "Muesca sintética: propagación del error hacia la fuente de fase\n"
        f"$T_e={notch.metrics['Te_K']:.3g}$ K, "
        f"$q={notch.metrics['q_um_inv']:.3g}$ $\\mu$m$^{{-1}}$, "
        f"$|\\Delta|_{{min}}/\\Delta_0={notch.metrics['notch_min_delta_over_delta0']:.3g}$"
    )
    fig.tight_layout()
    return fig


def _resolve_temperatures(
    catalog: DiagnosticCurrentCatalog,
    requested: Iterable[float] | None,
) -> np.ndarray:
    if requested is None:
        T_bias = float(catalog.metadata.get("T_bias_K", catalog.Te_axis_K[0]))
        raw = np.array((T_bias, 0.5 * catalog.Tc_K, 0.85 * catalog.Tc_K), dtype=float)
    else:
        raw = np.asarray(tuple(requested), dtype=float)
    clipped = np.clip(raw, catalog.Te_axis_K[0], catalog.Te_axis_K[-1])
    values = np.unique(np.round(clipped, decimals=12))
    if values.size < 1 or values.size > 3:
        raise ValueError("D1 supports one to three distinct diagnostic temperatures.")
    return values


def _resolve_q_values(
    catalog: DiagnosticCurrentCatalog,
    fractions: Iterable[float],
) -> tuple[np.ndarray, float]:
    calibration = catalog.metadata.get("calibration", {})
    q_critical = float(calibration.get("q_critical_m_inv", np.nan)) if isinstance(calibration, dict) else np.nan
    if not np.isfinite(q_critical) or q_critical <= 0.0:
        calibration_q = np.asarray(
            catalog.metadata.get("calibration_q_values_m_inv", []),
            dtype=float,
        )
        q_critical = float(np.max(calibration_q)) if calibration_q.size else 0.65 * catalog.q_axis_m_inv[-1]
    q_critical = min(q_critical, float(catalog.q_axis_m_inv[-1]))
    raw = np.asarray(tuple(fractions), dtype=float)
    if raw.size < 1 or raw.size > len(Q_COLORS) or np.any(raw <= 0.0):
        raise ValueError("Supply one to three positive q fractions.")
    q_min = float(catalog.q_axis_m_inv[1])
    values = np.clip(raw * q_critical, q_min, catalog.q_axis_m_inv[-1])
    return values, q_critical


def _save_formats(
    figure: plt.Figure,
    stem: Path,
    formats: Iterable[str],
    dpi: int,
) -> list[str]:
    paths: list[str] = []
    for extension in dict.fromkeys(str(value) for value in formats):
        path = stem.with_suffix(f".{extension}")
        figure.savefig(path, dpi=int(dpi), bbox_inches="tight")
        paths.append(str(path))
    plt.close(figure)
    return paths


def _build_summary(
    *,
    args: argparse.Namespace,
    catalog_path: Path,
    figures_dir: Path,
    catalog: DiagnosticCurrentCatalog,
    q_critical: float,
    curves: list[ConstitutiveCurves],
    notch: NotchDiagnostic,
    saved: dict[str, list[str]],
) -> dict[str, Any]:
    current_exponents = np.array(
        [item.metrics["current_low_amplitude_exponent"] for item in curves],
        dtype=float,
    )
    stiffness_exponents = np.array(
        [item.metrics["stiffness_low_amplitude_exponent"] for item in curves],
        dtype=float,
    )
    direct_exponents = np.array(
        [item.metrics["direct_low_amplitude_exponent"] for item in curves],
        dtype=float,
    )
    current_errors = np.array(
        [item.metrics["current_max_relative_error_below_first_node"] for item in curves],
        dtype=float,
    )
    stiffness_errors = np.array(
        [item.metrics["stiffness_max_relative_error_below_first_node"] for item in curves],
        dtype=float,
    )
    return {
        "schema_version": 1,
        "pipeline": "plot_pipelines/D1_usadel_low_amplitude_diagnostic.py",
        "scope": "temporary read-only diagnostic; no PRE or solver implementation changed",
        "pre_run_name": args.pre_run_name,
        "catalog_npz": str(catalog_path),
        "figures_dir": str(figures_dir),
        "figures": saved,
        "table": {
            "shape": list(catalog.js_A_m2.shape),
            "n_matsubara": int(catalog.n_matsubara),
            "D_m2_s": float(catalog.D_m2_s),
            "sigma_n_S_m": float(catalog.sigma_n_S_m),
            "Tc_K": float(catalog.Tc_K),
            "delta0_J": float(catalog.delta0_J),
            "first_positive_delta_J": float(catalog.first_positive_delta_J),
            "first_positive_delta_over_delta0": float(
                catalog.first_positive_delta_J / catalog.delta0_J
            ),
            "q_critical_m_inv": float(q_critical),
        },
        "constitutive_cases": [dict(item.metrics) for item in curves],
        "notch_case": dict(notch.metrics),
        "aggregate_metrics": {
            "current_exponent_median": float(np.median(current_exponents)),
            "stiffness_exponent_median": float(np.median(stiffness_exponents)),
            "direct_exponent_median": float(np.median(direct_exponents)),
            "current_error_max_below_first_node": float(np.max(current_errors)),
            "stiffness_error_max_below_first_node": float(np.max(stiffness_errors)),
            "current_notch_source_peak_over_direct": float(
                notch.metrics["current_source_peak_over_direct_peak"]
            ),
            "stiffness_notch_source_peak_over_direct": float(
                notch.metrics["stiffness_source_peak_over_direct_peak"]
            ),
        },
        "interpretation_guardrail": (
            "The 1D notch is a mechanism test, not proof that every notch in a "
            "coupled 2D run has this single cause."
        ),
    }


if __name__ == "__main__":
    raise SystemExit(main())
