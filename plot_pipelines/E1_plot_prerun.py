#!/usr/bin/env python3
"""E1 extra PRE plotting pipeline: Usadel equilibrium gap and supercurrent curve.

Input: an existing PRE-run produced by ``pipelines/01_prerun_template.py``.

Outputs, by default, inside ``plots/<pre-run-name>/figures/E1_prerun``:

- ``E1_usadel_gap_eq_vs_temperature.pdf``
- ``usadel_supercurrent_curve.png``

The Delta_eq(T, q) plot keeps the same lightweight/smoke-style defaults already
used by E1. The extra supercurrent curve is the same PRE diagnostic plot used by
``plot_pipelines/01_plot_prerun.py`` through ``plot_usadel_supercurrent_curve``.
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
from pysnspd.plotting.usadel_gap import load_usadel_gap_catalog, plot_gap_eq_vs_temperature
from pysnspd.usadel.catalog import load_usadel_catalog_npz


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "E1 extra plots: Delta_eq(T) and Usadel supercurrent curve "
            "reconstructed from an existing PRE-run Usadel catalog."
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
        "--pdf-name",
        default="E1_usadel_gap_eq_vs_temperature.pdf",
        help="Output PDF filename for the existing E1 Delta_eq(T, q) plot.",
    )
    parser.add_argument(
        "--supercurrent-name",
        default="usadel_supercurrent_curve.png",
        help="Output PNG filename for the added Usadel supercurrent curve plot.",
    )
    parser.add_argument(
        "--n-curves",
        type=int,
        default=4,
        help="Number of q curves in the Delta_eq(T, q) PDF, including q=0 and q=q_c.",
    )
    parser.add_argument(
        "--n-temperature",
        type=int,
        default=240,
        help="Number of temperature samples from T_min to slightly above Tc.",
    )
    parser.add_argument(
        "--T-min-K",
        type=float,
        default=None,
        help="Optional minimum temperature. Defaults to PRE metadata T_bias_K.",
    )
    parser.add_argument(
        "--T-max-K",
        type=float,
        default=None,
        help=(
            "Optional maximum temperature. Values at or below Tc are automatically "
            "extended to show Tc inside the axis."
        ),
    )
    parser.add_argument(
        "--q-critical-m-inv",
        type=float,
        default=None,
        help="Optional q_c override in m^-1. Defaults to PRE calibration metadata.",
    )
    parser.add_argument(
        "--n-matsubara",
        type=int,
        default=None,
        help="Optional Matsubara cutoff override. Defaults to PRE metadata n_matsubara_configured.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=480,
        help="Figure DPI.",
    )
    parser.add_argument(
        "--title",
        default=None,
        help="Optional title for the Delta_eq(T, q) PDF. Default keeps it memory-ready/title-free.",
    )
    parser.add_argument(
        "--progress",
        dest="progress",
        action="store_true",
        default=True,
        help="Show a progress bar while solving the Matsubara self-consistency points.",
    )
    parser.add_argument(
        "--no-progress",
        dest="progress",
        action="store_false",
        help="Disable the E1 solver progress bar.",
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

    pdf_name = args.pdf_name if args.pdf_name.lower().endswith(".pdf") else f"{args.pdf_name}.pdf"
    supercurrent_name = (
        args.supercurrent_name
        if args.supercurrent_name.lower().endswith(".png")
        else f"{args.supercurrent_name}.png"
    )

    gap_output_path = figures_dir / pdf_name
    supercurrent_output_path = figures_dir / supercurrent_name

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
        gap_output_path,
        dpi=int(args.dpi),
        title=args.title,
    )

    usadel_catalog = load_usadel_catalog_npz(catalog_path)
    supercurrent_output = plot_usadel_supercurrent_curve(
        usadel_catalog,
        supercurrent_output_path,
        dpi=int(args.dpi),
    )

    saved = {
        "usadel_gap_eq_vs_temperature_pdf": Path(gap_output),
        "usadel_supercurrent_curve_png": Path(supercurrent_output),
    }

    manifest_path = _write_manifest(
        pre_run_name=args.pre_run_name,
        raw_pre=raw_pre,
        figures_dir=figures_dir,
        catalog_path=catalog_path,
        saved=saved,
        metadata=gap_catalog.metadata,
        q_critical_m_inv=gap_catalog.q_critical_m_inv,
    )

    print("E1 pre-run Usadel plots")
    print(f" pre_run_name: {args.pre_run_name}")
    print(f" raw_pre: {raw_pre}")
    print(f" figures_dir: {figures_dir}")
    print(f" catalog_npz: {catalog_path}")
    print(f" gap_source: {gap_catalog.source_key}")
    print(f" q_c_m_inv: {gap_catalog.q_critical_m_inv:.8e}")
    print()
    print("Figures")
    for key, path in saved.items():
        print(f" {key}: {path}")
    print(f" manifest: {manifest_path}")
    print("Status: OK")
    return 0


def _write_manifest(
    *,
    pre_run_name: str,
    raw_pre: Path,
    figures_dir: Path,
    catalog_path: Path,
    saved: dict[str, Path],
    metadata: dict[str, Any],
    q_critical_m_inv: float,
) -> Path:
    manifest: dict[str, Any] = {
        "schema_version": 2,
        "pipeline": "plot_pipelines/E1_plot_prerun.py",
        "purpose": (
            "Extra PRE figures: Usadel equilibrium gap Delta_eq(T, q) and "
            "the standard Usadel supercurrent calibration curve."
        ),
        "pre_run_name": pre_run_name,
        "raw_pre": str(raw_pre),
        "figures_dir": str(figures_dir),
        "catalog_npz": str(catalog_path),
        "figures": {key: str(path) for key, path in saved.items()},
        "q_critical_m_inv": float(q_critical_m_inv),
        "gap_metadata": dict(metadata),
    }

    path = figures_dir / "E1_plot_prerun_manifest.yaml"
    with path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(manifest, f, sort_keys=False, allow_unicode=True, default_flow_style=False)
    return path


if __name__ == "__main__":
    raise SystemExit(main())
