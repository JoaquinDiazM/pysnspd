"""Create presentation figures from an existing PRE-run.

This post-processing pipeline reads raw PRE outputs from

    scratch/big_data/raw/<run_name>/pre

and writes figures to

    scratch/big_data/plots/<run_name>/figures

It does not modify the raw PRE-run.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import yaml

from pysnspd.config import load_config, validate_config
from pysnspd.io.manager import create_run_layout
from pysnspd.mesh.delaunay import load_mesh_npz
from pysnspd.mesh.edges import load_edges_npz
from pysnspd.usadel.catalog import load_usadel_catalog_npz
from pysnspd.plotting.pre_diagnostics import write_pre_diagnostic_plots


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot figures for a completed PRE-run."
    )
    parser.add_argument("--config", required=True, help="Path to YAML project config.")
    parser.add_argument("--run-name", required=True, help="Existing PRE-run name to plot.")
    parser.add_argument("--dpi", type=int, default=480)
    parser.add_argument(
        "--figures-subdir",
        default=None,
        help=(
            "Optional subdirectory inside plots/<run_name>/figures. "
            "By default figures are written directly to plots/<run_name>/figures."
        ),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    cfg = validate_config(load_config(args.config))
    layout = create_run_layout(cfg, args.run_name)

    raw_pre = Path(layout["raw_pre"])
    figures_dir = Path(layout["plots_figures"])
    if args.figures_subdir:
        figures_dir = figures_dir / str(args.figures_subdir)
    figures_dir.mkdir(parents=True, exist_ok=True)

    mesh = load_mesh_npz(raw_pre / "mesh.npz")
    edge_data = load_edges_npz(raw_pre / "edges.npz")
    usadel_catalog = load_usadel_catalog_npz(raw_pre / "usadel_dos_catalog.npz")

    saved_raw = write_pre_diagnostic_plots(
        mesh=mesh,
        edge_data=edge_data,
        usadel_catalog=usadel_catalog,
        output_dir=figures_dir,
        dpi=int(args.dpi),
    )
    saved = {key: Path(value) for key, value in saved_raw.items()}

    manifest_path = _write_plot_manifest(
        run_name=args.run_name,
        raw_pre=raw_pre,
        figures_dir=figures_dir,
        saved=saved,
    )

    print("PRE plotting pipeline")
    print(f" run_name: {args.run_name}")
    print(f" raw_pre: {raw_pre}")
    print(f" figures_dir: {figures_dir}")
    print()
    print("Figures")
    for key, path in saved.items():
        print(f" {key}: {path}")
    print(f" manifest: {manifest_path}")
    print("Status: OK")
    return 0


def _write_plot_manifest(
    *,
    run_name: str,
    raw_pre: Path,
    figures_dir: Path,
    saved: dict[str, Path],
) -> Path:
    manifest: dict[str, Any] = {
        "schema_version": 1,
        "pipeline": "plot_pipelines/01_plot_prerun.py",
        "purpose": "Presentation figures from an existing PRE-run.",
        "run_name": run_name,
        "raw_pre": str(raw_pre),
        "figures_dir": str(figures_dir),
        "figures": {key: str(path) for key, path in saved.items()},
    }
    out = figures_dir / "plot_prerun_manifest.yaml"
    with out.open("w", encoding="utf-8") as f:
        yaml.safe_dump(
            manifest,
            f,
            sort_keys=False,
            allow_unicode=True,
            default_flow_style=False,
        )
    return out


if __name__ == "__main__":
    raise SystemExit(main())
