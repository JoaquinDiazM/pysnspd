"""Z1 multi-run current-sweep inventory pipeline.

Z-series plot pipelines compare multiple runs.  This first Z1 pass does not yet
make current-sweep physics plots.  Instead it verifies that the project can
access the raw database: every requested run folder, every stage folder, every
NPZ file, and every YAML/JSON summary or manifest.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import yaml

from pysnspd.io.manager import create_run_layout
from pysnspd.io.run_database import (
    discover_raw_run_records,
    load_database_config,
    summarize_inventory,
    write_database_inventory,
)
from pysnspd.plotting.current_sweep import plot_current_sweep_placeholder


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Index raw data for a multi-run current-sweep analysis."
    )
    parser.add_argument("--config", required=True, help="Path to YAML project config.")
    parser.add_argument(
        "--run-names",
        nargs="*",
        default=None,
        help="Explicit run names to inspect. If omitted, scan raw/ folders.",
    )
    parser.add_argument(
        "--run-prefix",
        action="append",
        default=[],
        help=(
            "Only include discovered run names starting with this prefix. "
            "Can be given multiple times. Ignored for explicit --run-names."
        ),
    )
    parser.add_argument(
        "--stage",
        action="append",
        choices=("all", "pre", "ss", "photon"),
        default=None,
        help="Stage to inspect. Can be repeated. Default: all.",
    )
    parser.add_argument(
        "--output-run-name",
        default="Z1_current_sweep_analysis",
        help="Run-like folder name under big_data_root/plots for Z1 outputs.",
    )
    parser.add_argument(
        "--figures-subdir",
        default="figures",
        help="Subdirectory under plots/<output-run-name> for inventory files.",
    )
    parser.add_argument(
        "--no-npz-keys",
        action="store_true",
        help="List NPZ files but do not open them to read key/shape/dtype metadata.",
    )
    parser.add_argument(
        "--no-yaml-data",
        action="store_true",
        help="List YAML/JSON files but do not parse their contents.",
    )
    parser.add_argument("--dpi", type=int, default=480)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    cfg = load_database_config(args.config)
    layout = create_run_layout(cfg, args.output_run_name)
    plots_run = Path(layout["plots_run"])
    output_dir = plots_run / str(args.figures_subdir)
    output_dir.mkdir(parents=True, exist_ok=True)

    records = discover_raw_run_records(
        cfg,
        run_names=args.run_names,
        run_prefixes=[] if args.run_names else args.run_prefix,
        stages=args.stage,
        include_npz_keys=not args.no_npz_keys,
        include_yaml_data=not args.no_yaml_data,
    )
    summary = summarize_inventory(records)
    inventory_paths = write_database_inventory(records, output_dir)

    placeholder_status = plot_current_sweep_placeholder(
        records,
        output_dir,
        dpi=int(args.dpi),
    )

    manifest_path = _write_z1_manifest(
        output_dir=output_dir,
        config_path=args.config,
        args=vars(args),
        summary=summary,
        inventory_paths=inventory_paths,
        placeholder_status=placeholder_status,
    )

    print("Z1 current sweep analysis")
    print(f" output_run_name: {args.output_run_name}")
    print(f" output_dir: {output_dir}")
    print(f" runs indexed: {summary['n_runs']}")
    print(f" stage dirs found: {summary['n_stage_dirs']}")
    print(f" npz files found: {summary['n_npz_files']}")
    print(f" summary/manifest files found: {summary['n_summary_files']}")
    print()
    print("Inventory")
    for key, path in inventory_paths.items():
        print(f" {key}: {path}")
    print(f" manifest: {manifest_path}")
    print(f" placeholder_status: {placeholder_status}")
    print("Status: OK")
    return 0


def _write_z1_manifest(
    *,
    output_dir: Path,
    config_path: str | Path,
    args: dict[str, Any],
    summary: dict[str, Any],
    inventory_paths: dict[str, Path],
    placeholder_status: int,
) -> Path:
    manifest = {
        "schema_version": 1,
        "pipeline": "plot_pipelines/Z1_current_sweep_analysis.py",
        "purpose": "Multi-run current-sweep raw database inventory.",
        "config_path": str(config_path),
        "args": args,
        "summary": summary,
        "inventory_files": {key: str(path) for key, path in inventory_paths.items()},
        "plotting_placeholder_status": int(placeholder_status),
    }
    out = output_dir / "Z1_current_sweep_manifest.yaml"
    with out.open("w", encoding="utf-8") as f:
        yaml.safe_dump(manifest, f, sort_keys=False, allow_unicode=True, default_flow_style=False)
    return out


if __name__ == "__main__":
    raise SystemExit(main())
