"""Z2 multi-run current-sweep analysis pipeline."""

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
from pysnspd.plotting.current_sweep import make_current_sweep_figures
from pysnspd.plotting.style import THESIS_DPI


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Index raw data and create current-sweep IV figures."
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
        default="Z2_current_sweep_analysis",
        help="Run-like folder name under big_data_root/plots for Z2 outputs.",
    )
    parser.add_argument(
        "--figures-subdir",
        default="figures",
        help="Subdirectory under plots/<output-run-name> for inventory files and figures.",
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
    parser.add_argument(
        "--voltage-probe-offset-nm",
        type=float,
        default=50.0,
        help=(
            "Probe the IV voltage from the x-profile as "
            "phi(x_center + offset) - phi(x_center - offset)."
        ),
    )
    parser.add_argument(
        "--voltage-probe-half-window-nm",
        type=float,
        default=None,
        help=(
            "Optional half-width of the averaging window around each probe x. "
            "Default: infer from x-profile bin spacing."
        ),
    )
    parser.add_argument(
        "--delta-inset-currents-uA",
        nargs=4,
        type=float,
        default=None,
        metavar=("I1", "I2", "I3", "I4"),
        help=(
            "Exactly four currents (in microampere) used to select the four "
            "|Delta| inset colormaps on the IV figure. Each requested current "
            "is matched to the nearest available SS run."
        ),
    )
    parser.add_argument(
        "--terminal-delta-inset-currents-uA",
        nargs=3,
        type=float,
        default=None,
        metavar=("I1", "I2", "I3"),
        help=(
            "Exactly three currents (in microampere) used for the vertically "
            "stacked full-strip |Delta| maps beside the terminal-voltage IV curve."
        ),
    )
    parser.add_argument(
        "--no-origin",
        action="store_true",
        help="Do not prepend the synthetic (I,V)=(0,0) point.",
    )
    parser.add_argument("--dpi", type=int, default=THESIS_DPI)
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

    figure_outputs = make_current_sweep_figures(
        config_path=args.config,
        records=records,
        output_dir=output_dir,
        dpi=int(args.dpi),
        voltage_probe_offset_nm=float(args.voltage_probe_offset_nm),
        voltage_probe_half_window_nm=args.voltage_probe_half_window_nm,
        include_origin=not args.no_origin,
        delta_inset_currents_uA=args.delta_inset_currents_uA,
        terminal_delta_inset_currents_uA=args.terminal_delta_inset_currents_uA,
    )

    manifest_path = _write_z2_manifest(
        output_dir=output_dir,
        config_path=args.config,
        args=vars(args),
        summary=summary,
        inventory_paths=inventory_paths,
        figure_outputs=figure_outputs,
    )

    iv_summary = figure_outputs.get("iv_summary", {}) if isinstance(figure_outputs, dict) else {}

    print("Z2 current sweep analysis")
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
    print()
    print("Figures / tables")
    for key in (
        "iv_curve",
        "terminal_iv_curve",
        "iv_points_csv",
        "iv_points_yaml",
        "iv_skipped_yaml",
        "iv_insets_yaml",
        "terminal_iv_insets_yaml",
    ):
        path = figure_outputs.get(key)
        if path:
            print(f" {key}: {path}")
    print(f" manifest: {manifest_path}")
    if iv_summary:
        print("IV summary")
        print(f" points: {iv_summary.get('n_points', 0)}")
        print(f" runs loaded: {iv_summary.get('n_runs_loaded', 0)}")
        print(f" runs skipped: {iv_summary.get('n_runs_skipped', 0)}")
        print(f" probe offset [nm]: {iv_summary.get('voltage_probe_offset_nm')}")
        print(f" probe half-window [nm]: {iv_summary.get('voltage_probe_half_window_nm')}")
        print(f" sign flipped: {iv_summary.get('voltage_sign_flipped')}")
        print(f" terminal sign flipped: {iv_summary.get('terminal_voltage_sign_flipped')}")
        print(f" terminal normal resistance [ohm]: {iv_summary.get('normal_resistance_terminal_ohm')}")
        print(f" delta inset requests [uA]: {iv_summary.get('delta_inset_currents_uA', [])}")
        print(f" delta inset resolved [uA]: {iv_summary.get('delta_inset_resolved_currents_uA', [])}")
        print(f" terminal delta requests [uA]: {iv_summary.get('terminal_delta_inset_currents_uA', [])}")
        print(f" terminal delta resolved [uA]: {iv_summary.get('terminal_delta_inset_resolved_currents_uA', [])}")
    print("Status: OK")
    return 0



def _write_z2_manifest(
    *,
    output_dir: Path,
    config_path: str | Path,
    args: dict[str, Any],
    summary: dict[str, Any],
    inventory_paths: dict[str, Path],
    figure_outputs: dict[str, Any],
) -> Path:
    manifest = {
        "schema_version": 4,
        "pipeline": "plot_pipelines/Z2_current_sweep_analysis.py",
        "purpose": (
            "Multi-run current-sweep inventory with central TDGL and terminal-voltage "
            "IV figures, full-length normal-state reference, and |Delta| snapshots."
        ),
        "config_path": str(config_path),
        "args": args,
        "summary": summary,
        "inventory_files": {key: str(path) for key, path in inventory_paths.items()},
        "figure_outputs": {
            key: (str(value) if isinstance(value, Path) else value)
            for key, value in figure_outputs.items()
        },
    }
    out = output_dir / "Z2_current_sweep_manifest.yaml"
    with out.open("w", encoding="utf-8") as f:
        yaml.safe_dump(manifest, f, sort_keys=False, allow_unicode=True, default_flow_style=False)
    return out


if __name__ == "__main__":
    raise SystemExit(main())
