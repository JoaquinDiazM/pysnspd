#!/usr/bin/env python3
"""E1 PRE plotting pipeline.

This plotting-only pipeline reads an existing PRE-run Usadel catalogue and writes
E-type thesis figures in PDF format.

Default output is intentionally fast and does not reconstruct the expensive
Delta_eq(T) figure:

- usadel_supercurrent_curve.pdf
- E1_usadel_dos_curves_delta_eq.pdf
- E1_usadel_dos_curves_delta0.pdf

Use --with-gap-plot only when the additional Delta_eq(T) figure is needed.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any

import yaml

# Allow running as ``python plot_pipelines/E1_plot_prerun.py`` from repo root.
_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from pysnspd.config import load_config, validate_config
from pysnspd.io.manager import create_run_layout
from pysnspd.plotting.pre_diagnostics import plot_usadel_supercurrent_curve
from pysnspd.plotting.usadel_dos_curves import (
    plot_usadel_dos_curves_equilibrium_gap,
    plot_usadel_dos_curves_fixed_delta0,
)
from pysnspd.plotting.usadel_gap import load_usadel_gap_catalog, plot_gap_eq_vs_temperature
from pysnspd.usadel.catalog import load_usadel_catalog_npz


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "E1 PRE plotting pipeline: write the Usadel supercurrent PDF, two DOS-curve PDFs, "
            "and optionally the expensive Delta_eq(T) PDF."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--config",
        required=True,
        type=Path,
        help="pySNSPD YAML config used by the PRE-run.",
    )
    parser.add_argument(
        "--pre-run-name",
        required=True,
        help="Existing PRE-run name produced by 01_prerun_template.py.",
    )
    parser.add_argument(
        "--catalog-npz",
        type=Path,
        default=None,
        help=(
            "Optional direct path to the Usadel .npz catalog; defaults to "
            "raw/<pre-run-name>/pre/usadel_dos_catalog.npz."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Optional output directory; defaults to plots/<pre-run-name>/figures/E1_prerun.",
    )

    parser.add_argument(
        "--supercurrent-pdf-name",
        default="usadel_supercurrent_curve.pdf",
        help="Output PDF filename for the supercurrent calibration figure.",
    )
    parser.add_argument(
        "--dos-eq-pdf-name",
        default="E1_usadel_dos_curves_delta_eq.pdf",
        help="Output PDF filename for the DOS curves evaluated at Delta_eq(q).",
    )
    parser.add_argument(
        "--dos-delta0-pdf-name",
        default="E1_usadel_dos_curves_delta0.pdf",
        help="Output PDF filename for the DOS curves evaluated at fixed Delta_0.",
    )
    parser.add_argument(
        "--skip-dos-curves",
        action="store_true",
        help="Only write the supercurrent figure; skip the two DOS-curve PDFs.",
    )
    parser.add_argument(
        "--dos-current-fractions",
        nargs="+",
        type=float,
        default=(0.0, 0.50, 0.65, 0.80, 0.95),
        help="Current fractions I_s/I_c used for both DOS-curve PDFs.",
    )
    parser.add_argument(
        "--dos-energy-max-meV",
        type=float,
        default=None,
        help="Optional upper energy limit for the DOS-curve PDFs.",
    )
    parser.add_argument(
        "--dos-energy-window",
        dest="dos_energy_window",
        action="store_true",
        default=True,
        help="Use an automatic compact energy window for the DOS-curve PDFs.",
    )
    parser.add_argument(
        "--no-dos-energy-window",
        dest="dos_energy_window",
        action="store_false",
        help="Use the full catalogue energy axis for the DOS-curve PDFs.",
    )

    parser.add_argument(
        "--with-gap-plot",
        action="store_true",
        help="Also reconstruct and write the expensive Delta_eq(T) PDF.",
    )
    parser.add_argument(
        "--gap-pdf-name",
        default="E1_usadel_gap_eq_vs_temperature.pdf",
        help="Output PDF filename for the optional Delta_eq(T) figure.",
    )
    parser.add_argument(
        "--n-curves",
        type=int,
        default=4,
        help="Number of q curves for the optional Delta_eq(T) plot, including q=0 and q=q_c.",
    )
    parser.add_argument(
        "--n-temperature",
        type=int,
        default=64,
        help=(
            "Number of temperature samples for the optional Delta_eq(T) plot. "
            "Default is a fast smoke-style resolution."
        ),
    )
    parser.add_argument(
        "--T-min-K",
        type=float,
        default=None,
        help="Optional minimum temperature for the optional Delta_eq(T) plot. Defaults to PRE metadata T_bias_K.",
    )
    parser.add_argument(
        "--T-max-K",
        type=float,
        default=None,
        help=(
            "Optional maximum temperature for the optional Delta_eq(T) plot. Values at or below Tc "
            "are automatically extended to show Tc inside the axis."
        ),
    )
    parser.add_argument(
        "--q-critical-m-inv",
        type=float,
        default=None,
        help="Optional q_c override in m^-1 for the optional Delta_eq(T) plot.",
    )
    parser.add_argument(
        "--n-matsubara",
        type=int,
        default=96,
        help=(
            "Optional Matsubara cutoff for the optional Delta_eq(T) reconstruction. "
            "Default is a fast smoke-style value."
        ),
    )
    parser.add_argument(
        "--progress",
        dest="progress",
        action="store_true",
        default=False,
        help="Show a progress bar while solving the optional Delta_eq(T) Matsubara self-consistency points.",
    )
    parser.add_argument(
        "--no-progress",
        dest="progress",
        action="store_false",
        help="Disable the progress bar for the optional Delta_eq(T) solve.",
    )

    parser.add_argument(
        "--dpi",
        type=int,
        default=480,
        help="PDF rasterization DPI for any rasterized artists.",
    )
    parser.add_argument(
        "--title",
        default=None,
        help="Optional title for the optional Delta_eq(T) plot. Default keeps the figure title-free.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    cfg = validate_config(load_config(args.config))
    layout = create_run_layout(cfg, args.pre_run_name)

    raw_pre = Path(layout["raw_pre"])
    figures_dir = (
        Path(layout["plots_figures"]) / "E1_prerun"
        if args.output_dir is None
        else args.output_dir.expanduser().resolve()
    )
    figures_dir.mkdir(parents=True, exist_ok=True)

    catalog_path = (
        args.catalog_npz.expanduser().resolve()
        if args.catalog_npz is not None
        else raw_pre / "usadel_dos_catalog.npz"
    )
    if not catalog_path.exists():
        raise FileNotFoundError(f"Usadel catalog not found: {catalog_path}")

    saved: dict[str, Path] = {}
    usadel_catalog = load_usadel_catalog_npz(catalog_path)

    supercurrent_output = plot_usadel_supercurrent_curve(
        usadel_catalog,
        figures_dir / _ensure_pdf_name(args.supercurrent_pdf_name),
        dpi=int(args.dpi),
    )
    saved["usadel_supercurrent_curve_pdf"] = supercurrent_output

    if not args.skip_dos_curves:
        dos_eq_output = plot_usadel_dos_curves_equilibrium_gap(
            usadel_catalog,
            figures_dir / _ensure_pdf_name(args.dos_eq_pdf_name),
            current_fractions=tuple(float(v) for v in args.dos_current_fractions),
            dpi=int(args.dpi),
            energy_max_meV=args.dos_energy_max_meV,
            energy_window=bool(args.dos_energy_window),
        )
        saved["usadel_dos_curves_delta_eq_pdf"] = dos_eq_output

        dos_delta0_output = plot_usadel_dos_curves_fixed_delta0(
            usadel_catalog,
            figures_dir / _ensure_pdf_name(args.dos_delta0_pdf_name),
            current_fractions=tuple(float(v) for v in args.dos_current_fractions),
            dpi=int(args.dpi),
            energy_max_meV=args.dos_energy_max_meV,
            energy_window=bool(args.dos_energy_window),
        )
        saved["usadel_dos_curves_delta0_pdf"] = dos_delta0_output

    gap_source = "not_requested"
    q_critical_m_inv = None
    gap_metadata: dict[str, Any] = {}
    if args.with_gap_plot:
        gap_catalog = load_usadel_gap_catalog(
            catalog_path,
            n_curves=int(args.n_curves),
            n_temperature=int(args.n_temperature),
            T_min_K=args.T_min_K,
            T_max_K=args.T_max_K,
            q_critical_m_inv=args.q_critical_m_inv,
            n_matsubara=args.n_matsubara,
            progress=bool(args.progress),
        )
        gap_output = plot_gap_eq_vs_temperature(
            gap_catalog,
            figures_dir / _ensure_pdf_name(args.gap_pdf_name),
            dpi=int(args.dpi),
            title=args.title,
        )
        saved["usadel_gap_eq_vs_temperature_pdf"] = gap_output
        gap_source = gap_catalog.source_key
        q_critical_m_inv = float(gap_catalog.q_critical_m_inv)
        gap_metadata = dict(gap_catalog.metadata)

    manifest_path = _write_manifest(
        pre_run_name=args.pre_run_name,
        raw_pre=raw_pre,
        figures_dir=figures_dir,
        catalog_path=catalog_path,
        saved=saved,
        skip_dos_curves=bool(args.skip_dos_curves),
        dos_settings={
            "current_fractions": [float(v) for v in args.dos_current_fractions],
            "energy_max_meV": args.dos_energy_max_meV,
            "energy_window": bool(args.dos_energy_window),
        },
        with_gap_plot=bool(args.with_gap_plot),
        gap_source=gap_source,
        q_critical_m_inv=q_critical_m_inv,
        gap_metadata=gap_metadata,
        gap_settings={
            "n_curves": int(args.n_curves),
            "n_temperature": int(args.n_temperature),
            "n_matsubara": int(args.n_matsubara) if args.n_matsubara is not None else None,
            "T_min_K": args.T_min_K,
            "T_max_K": args.T_max_K,
            "progress": bool(args.progress),
        },
    )

    print("E1 pre-run Usadel plots")
    print(f" pre_run_name: {args.pre_run_name}")
    print(f" raw_pre: {raw_pre}")
    print(f" figures_dir: {figures_dir}")
    print(f" catalog_npz: {catalog_path}")
    print()
    print("Figures")
    for key, path in saved.items():
        print(f" {key}: {path}")
    print(f" manifest: {manifest_path}")
    print("Status: OK")
    return 0


def _ensure_pdf_name(name: str | Path) -> str:
    raw = str(name)
    return raw if raw.lower().endswith(".pdf") else f"{raw}.pdf"


def _write_manifest(
    *,
    pre_run_name: str,
    raw_pre: Path,
    figures_dir: Path,
    catalog_path: Path,
    saved: dict[str, Path],
    skip_dos_curves: bool,
    dos_settings: dict[str, Any],
    with_gap_plot: bool,
    gap_source: str,
    q_critical_m_inv: float | None,
    gap_metadata: dict[str, Any],
    gap_settings: dict[str, Any],
) -> Path:
    manifest: dict[str, Any] = {
        "schema_version": 3,
        "pipeline": "plot_pipelines/E1_plot_prerun.py",
        "purpose": "E-type PRE figures in PDF format: supercurrent curve, DOS curves, and optional Delta_eq(T,q).",
        "pre_run_name": pre_run_name,
        "raw_pre": str(raw_pre),
        "figures_dir": str(figures_dir),
        "catalog_npz": str(catalog_path),
        "figures": {key: str(path) for key, path in saved.items()},
        "skip_dos_curves": bool(skip_dos_curves),
        "dos_settings": dos_settings,
        "with_gap_plot": bool(with_gap_plot),
        "gap_source": gap_source,
        "q_critical_m_inv": None if q_critical_m_inv is None else float(q_critical_m_inv),
        "gap_settings": gap_settings,
        "gap_metadata": gap_metadata,
    }
    path = figures_dir / "E1_plot_prerun_manifest.yaml"
    with path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(manifest, f, sort_keys=False, allow_unicode=True, default_flow_style=False)
    return path


if __name__ == "__main__":
    raise SystemExit(main())
