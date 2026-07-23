"""Plot figures for a completed pipeline 03 photon/circuit run."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import yaml

from pysnspd.config import load_config, validate_config
from pysnspd.io.manager import create_run_layout
from pysnspd.mesh.delaunay import load_mesh_npz
from pysnspd.plotting.photon_figures import load_npz_dict, make_photon_run_figures
from pysnspd.plotting.style import THESIS_DPI


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot figures for a completed photon/circuit transient run.")
    parser.add_argument("--config", required=True, help="Path to YAML project config.")
    parser.add_argument("--run-name", required=True, help="Existing pipeline 03 run name to plot.")
    parser.add_argument(
        "--pre-run-name",
        default=None,
        help="PRE run name. Required when --scalar-times-ps is used because the mesh is loaded from PRE.",
    )
    parser.add_argument("--dpi", type=int, default=THESIS_DPI)
    parser.add_argument("--center-width-nm", type=float, default=100.0)
    parser.add_argument(
        "--scalar-times-ps",
        nargs="*",
        type=float,
        default=None,
        help=(
            "Requested times in ps for the multi-row center scalar map. "
            "The closest stored transient snapshot is used for each requested time."
        ),
    )
    parser.add_argument(
        "--figures-subdir",
        default=None,
        help="Optional subdirectory inside plots/<run-name>/figures.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    cfg = validate_config(load_config(args.config))
    layout = create_run_layout(cfg, args.run_name)

    raw_photon = Path(layout["raw_photon"])
    figures_dir = Path(layout["plots_figures"])
    if args.figures_subdir:
        figures_dir = figures_dir / str(args.figures_subdir)
    figures_dir.mkdir(parents=True, exist_ok=True)

    history_path = raw_photon / "transient_history.npz"
    snapshots_path = raw_photon / "transient_snapshots.npz"
    summary_path = raw_photon / "photon_summary.yaml"

    history = load_npz_dict(history_path)
    summary = _read_yaml(summary_path)

    mesh = None
    snapshots = None
    scalar_times_ps = args.scalar_times_ps

    if scalar_times_ps is not None and len(scalar_times_ps) > 0:
        if not args.pre_run_name:
            raise ValueError("--pre-run-name is required when --scalar-times-ps is used.")
        pre_layout = create_run_layout(cfg, args.pre_run_name)
        raw_pre = Path(pre_layout["raw_pre"])
        mesh_path = raw_pre / "mesh.npz"
        if not mesh_path.exists():
            raise FileNotFoundError(f"Missing PRE mesh: {mesh_path}")
        if not snapshots_path.exists():
            raise FileNotFoundError(f"Missing transient snapshots: {snapshots_path}")

        mesh = load_mesh_npz(mesh_path)
        snapshots = load_npz_dict(snapshots_path)

    saved = make_photon_run_figures(
        history=history,
        summary=summary,
        output_dir=figures_dir,
        dpi=int(args.dpi),
        mesh=mesh,
        snapshots=snapshots,
        scalar_times_ps=scalar_times_ps,
        center_width_nm=float(args.center_width_nm),
    )

    manifest_path = _write_plot_manifest(
        run_name=args.run_name,
        raw_photon=raw_photon,
        figures_dir=figures_dir,
        saved=saved,
        history=history,
        summary=summary,
        scalar_times_ps=scalar_times_ps,
        center_width_nm=float(args.center_width_nm),
        pre_run_name=args.pre_run_name,
    )

    print("Photon plotting pipeline")
    print(f" run_name: {args.run_name}")
    if args.pre_run_name:
        print(f" pre_run_name: {args.pre_run_name}")
    print(f" raw_photon: {raw_photon}")
    print(f" figures_dir: {figures_dir}")
    if scalar_times_ps:
        print(f" scalar_times_ps: {', '.join(f'{t:g}' for t in scalar_times_ps)}")
    print()
    print("Figures")
    for key, path in saved.items():
        print(f" {key}: {path}")
    print(f" manifest: {manifest_path}")
    print("Status: OK")
    return 0


def _read_yaml(path: str | Path) -> dict[str, Any]:
    p = Path(path)
    if not p.exists():
        return {}
    with p.open("r", encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
    return data if isinstance(data, dict) else {}


def _write_plot_manifest(
    *,
    run_name: str,
    raw_photon: Path,
    figures_dir: Path,
    saved: dict[str, Path],
    history: dict[str, Any],
    summary: dict[str, Any],
    scalar_times_ps: list[float] | None,
    center_width_nm: float,
    pre_run_name: str | None,
) -> Path:
    manifest = {
        "schema_version": 2,
        "pipeline": "plot_pipelines/03_plot_photon_run.py",
        "purpose": "Presentation figures from an existing pipeline 03 photon/circuit transient run.",
        "run_name": str(run_name),
        "pre_run_name": None if pre_run_name is None else str(pre_run_name),
        "raw_photon": str(raw_photon),
        "figures_dir": str(figures_dir),
        "center_width_nm": float(center_width_nm),
        "scalar_times_ps": None if scalar_times_ps is None else [float(t) for t in scalar_times_ps],
        "figures": {key: str(path) for key, path in saved.items()},
        "history_keys": sorted(str(key) for key in history.keys()),
        "summary_keys": sorted(str(key) for key in summary.keys()),
    }
    out = figures_dir / "plot_manifest.yaml"
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
