#!/usr/bin/env python3
"""E1 extra PRE plotting pipeline: Usadel equilibrium gap versus temperature.

Input: an existing PRE-run produced by ``pipelines/01_prerun_template.py``.
Output: a memory-ready PDF with four ``Delta_eq(T)`` curves from q=0 to q_c.
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
from pysnspd.plotting.usadel_gap import load_usadel_gap_catalog, plot_gap_eq_vs_temperature


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="E1 extra plot: Delta_eq(T) reconstructed from a PRE-run Usadel catalog, saved as PDF.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--config", required=True, type=Path, help="pySNSPD YAML config used by the PRE-run.")
    parser.add_argument("--pre-run-name", required=True, help="Existing PRE-run name produced by 01_prerun_template.py.")
    parser.add_argument(
        "--catalog-npz",
        type=Path,
        default=None,
        help="Optional direct path to the Usadel .npz catalog; defaults to raw/<pre-run-name>/pre/usadel_dos_catalog.npz.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Optional output directory; defaults to plots/<pre-run-name>/figures/E1_prerun.",
    )
    parser.add_argument("--pdf-name", default="E1_usadel_gap_eq_vs_temperature.pdf", help="Output PDF filename.")
    parser.add_argument("--n-curves", type=int, default=4, help="Number of q curves, including q=0 and q=q_c.")
    parser.add_argument("--n-temperature", type=int, default=240, help="Number of temperature samples from T_min to slightly above Tc.")
    parser.add_argument("--T-min-K", type=float, default=None, help="Optional minimum temperature. Defaults to PRE metadata T_bias_K.")
    parser.add_argument(
        "--T-max-K",
        type=float,
        default=None,
        help="Optional maximum temperature. Values at or below Tc are automatically extended to show Tc inside the axis.",
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
    parser.add_argument("--dpi", type=int, default=480, help="PDF rasterization DPI for any rasterized artists.")
    parser.add_argument("--title", default=None, help="Optional figure title. Default keeps the memory-ready figure title-free.")
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
    output_path = figures_dir / pdf_name

    catalog = load_usadel_gap_catalog(
        catalog_path,
        n_curves=int(args.n_curves),
        n_temperature=int(args.n_temperature),
        T_min_K=args.T_min_K,
        T_max_K=args.T_max_K,
        q_critical_m_inv=args.q_critical_m_inv,
        n_matsubara=args.n_matsubara,
        progress=bool(args.progress),
    )
    output = plot_gap_eq_vs_temperature(
        catalog,
        output_path,
        dpi=int(args.dpi),
        title=args.title,
    )

    manifest_path = _write_manifest(
        pre_run_name=args.pre_run_name,
        raw_pre=raw_pre,
        figures_dir=figures_dir,
        catalog_path=catalog_path,
        output=output,
        metadata=catalog.metadata,
        q_critical_m_inv=catalog.q_critical_m_inv,
    )

    print("E1 pre-run Usadel gap plot")
    print(f" pre_run_name: {args.pre_run_name}")
    print(f" raw_pre:      {raw_pre}")
    print(f" figures_dir:  {figures_dir}")
    print(f" catalog_npz:  {catalog_path}")
    print(f" gap_source:   {catalog.source_key}")
    print(f" q_c_m_inv:    {catalog.q_critical_m_inv:.8e}")
    print(f" output_pdf:   {output}")
    print(f" manifest:     {manifest_path}")
    print("Status: OK")
    return 0


def _write_manifest(
    *,
    pre_run_name: str,
    raw_pre: Path,
    figures_dir: Path,
    catalog_path: Path,
    output: Path,
    metadata: dict[str, Any],
    q_critical_m_inv: float,
) -> Path:
    manifest: dict[str, Any] = {
        "schema_version": 1,
        "pipeline": "plot_pipelines/E1_plot_prerun.py",
        "purpose": "Extra PDF figure: Usadel equilibrium gap Delta_eq(T, q) reconstructed from the PRE-run catalog.",
        "pre_run_name": pre_run_name,
        "raw_pre": str(raw_pre),
        "figures_dir": str(figures_dir),
        "catalog_npz": str(catalog_path),
        "output_pdf": str(output),
        "q_critical_m_inv": float(q_critical_m_inv),
        "gap_metadata": dict(metadata),
    }
    path = figures_dir / "E1_plot_prerun_manifest.yaml"
    with path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(manifest, f, sort_keys=False, allow_unicode=True, default_flow_style=False)
    return path


if __name__ == "__main__":
    raise SystemExit(main())
