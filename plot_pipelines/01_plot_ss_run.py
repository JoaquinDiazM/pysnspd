"""Create presentation figures from an existing stationary SS run.

This is a post-processing pipeline.  It reads raw SS outputs from
``scratch/big_data/raw/<run>/ss`` and writes figures to
``scratch/big_data/plots/<run>/figures``.  It does not modify the raw run.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import yaml

from pysnspd.analysis.ss_run import build_ss_plot_dataset, load_ss_run
from pysnspd.plotting.ss_figures import make_ss_run_figures


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot figures for a completed stationary SS run.")
    parser.add_argument("--config", required=True, help="Path to YAML project config.")
    parser.add_argument("--run-name", required=True, help="Existing SS run name to plot.")
    parser.add_argument(
        "--pre-run-name",
        default=None,
        help="PRE run containing mesh/edge files. If omitted, read from ss_summary.yaml.",
    )
    parser.add_argument("--dpi", type=int, default=480)
    parser.add_argument(
        "--figures-subdir",
        default=None,
        help=(
            "Optional subdirectory inside plots/<run>/figures. "
            "By default figures are written directly to plots/<run>/figures."
        ),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    run = load_ss_run(
        config_path=args.config,
        run_name=args.run_name,
        pre_run_name=args.pre_run_name,
    )
    dataset = build_ss_plot_dataset(run)
    figures_dir = run.figures_dir
    if args.figures_subdir:
        figures_dir = figures_dir / str(args.figures_subdir)
    figures_dir.mkdir(parents=True, exist_ok=True)

    saved = make_ss_run_figures(mesh=run.mesh, dataset=dataset, output_dir=figures_dir, dpi=int(args.dpi))
    manifest_path = _write_plot_manifest(run=run, figures_dir=figures_dir, saved=saved, dataset=dataset)

    print("SS plotting pipeline")
    print(f"  run_name:      {run.run_name}")
    print(f"  pre_run_name:  {run.pre_run_name}")
    print(f"  raw_ss:        {run.raw_ss}")
    print(f"  figures_dir:   {figures_dir}")
    print()
    print("Figures")
    for key, path in saved.items():
        print(f"  {key}: {path}")
    print(f"  manifest: {manifest_path}")
    print("Status: OK")
    return 0


def _write_plot_manifest(*, run, figures_dir: Path, saved: dict[str, Path], dataset: dict[str, Any]) -> Path:
    manifest = {
        "schema_version": 1,
        "pipeline": "plot_pipelines/01_plot_ss_run.py",
        "purpose": "Presentation figures from an existing stationary SS run.",
        "run_name": run.run_name,
        "pre_run_name": run.pre_run_name,
        "raw_ss": str(run.raw_ss),
        "figures_dir": str(figures_dir),
        "figures": {key: str(path) for key, path in saved.items()},
        "summary_scalars": dataset.get("summary_scalars", {}),
        "npz_keys": dataset.get("npz_keys", {}),
    }
    out = figures_dir / "plot_manifest.yaml"
    with out.open("w", encoding="utf-8") as f:
        yaml.safe_dump(manifest, f, sort_keys=False, allow_unicode=True, default_flow_style=False)
    return out


if __name__ == "__main__":
    raise SystemExit(main())
