#!/usr/bin/env python3
"""Compare two completed photon runs that differ only in impact position."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import yaml

from pysnspd.config import load_config, validate_config
from pysnspd.io.manager import create_run_layout
from pysnspd.mesh.delaunay import load_mesh_npz
from pysnspd.plotting.photon_comparison import make_photon_position_figures
from pysnspd.plotting.style import THESIS_DPI
from pysnspd.analysis.timing import analyze_photon_timing
from pysnspd.analysis.timing_cli import (
    add_timing_analysis_arguments,
    timing_criteria_from_args,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create thesis-ready PDF comparisons from two completed photon runs."
    )
    parser.add_argument("--config", required=True, help="Absolute YAML project configuration.")
    parser.add_argument("--pre-run-name", required=True, help="PRE run containing the common mesh.")
    parser.add_argument("--center-run-name", required=True, help="Completed central-impact photon run.")
    parser.add_argument("--edge-run-name", required=True, help="Completed edge-impact photon run.")
    parser.add_argument(
        "--times-ps",
        nargs="+",
        type=float,
        default=(50.0, 55.0, 60.0, 100.0),
        help="Requested matched field-map times; the nearest stored snapshot is used.",
    )
    parser.add_argument(
        "--output-run-name",
        default="E3_photon_position_comparison",
        help="Run-like directory under big_data_root/plots.",
    )
    parser.add_argument("--figures-subdir", default="figures")
    parser.add_argument("--dpi", type=int, default=THESIS_DPI)
    add_timing_analysis_arguments(parser)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    cfg = validate_config(load_config(args.config))
    pre_layout = create_run_layout(cfg, args.pre_run_name)
    center_layout = create_run_layout(cfg, args.center_run_name)
    edge_layout = create_run_layout(cfg, args.edge_run_name)
    output_layout = create_run_layout(cfg, args.output_run_name)

    raw_pre = Path(pre_layout["raw_pre"])
    center_raw = Path(center_layout["raw_photon"])
    edge_raw = Path(edge_layout["raw_photon"])
    output_dir = Path(output_layout["plots_run"]) / str(args.figures_subdir)
    output_dir.mkdir(parents=True, exist_ok=True)

    mesh = load_mesh_npz(_require_file(raw_pre / "mesh.npz", "PRE mesh"))
    center_history = _load_npz(_require_file(center_raw / "transient_history.npz", "center history"))
    center_snapshots = _load_npz(_require_file(center_raw / "transient_snapshots.npz", "center snapshots"))
    center_summary = _read_yaml(_require_file(center_raw / "photon_summary.yaml", "center summary"))
    edge_history = _load_npz(_require_file(edge_raw / "transient_history.npz", "edge history"))
    edge_snapshots = _load_npz(_require_file(edge_raw / "transient_snapshots.npz", "edge snapshots"))
    edge_summary = _read_yaml(_require_file(edge_raw / "photon_summary.yaml", "edge summary"))
    delta0_meV = _read_delta0_meV(raw_pre)
    detection_criteria, recovery_criteria = timing_criteria_from_args(args)
    center_timing = analyze_photon_timing(
        center_history,
        snapshots=center_snapshots,
        detection=detection_criteria,
        recovery=recovery_criteria,
    )
    edge_timing = analyze_photon_timing(
        edge_history,
        snapshots=edge_snapshots,
        detection=detection_criteria,
        recovery=recovery_criteria,
    )

    saved = make_photon_position_figures(
        mesh=mesh,
        center_history=center_history,
        center_snapshots=center_snapshots,
        center_summary=center_summary,
        edge_history=edge_history,
        edge_snapshots=edge_snapshots,
        edge_summary=edge_summary,
        delta0_meV=delta0_meV,
        requested_times_ps=args.times_ps,
        output_dir=output_dir,
        dpi=int(args.dpi),
        center_timing=center_timing,
        edge_timing=edge_timing,
    )
    manifest_path = _write_manifest(
        args=args,
        raw_pre=raw_pre,
        center_raw=center_raw,
        edge_raw=edge_raw,
        output_dir=output_dir,
        saved=saved,
        delta0_meV=delta0_meV,
        center_history=center_history,
        edge_history=edge_history,
        center_snapshots=center_snapshots,
        edge_snapshots=edge_snapshots,
        center_timing=center_timing,
        edge_timing=edge_timing,
    )

    print("E3 photon-position comparison")
    print(f" center_run: {args.center_run_name}")
    print(f" edge_run:   {args.edge_run_name}")
    print(f" output_dir: {output_dir}")
    print(f" Delta_BCS(0): {delta0_meV:.9g} meV")
    print(
        " center timing: "
        f"t_lat={dict(center_timing.get('latency', {})).get('t_lat_ps', 'censored')} ps, "
        f"t_rec={dict(dict(center_timing.get('recovery', {})).get('selected', {})).get('t_rec_ps', 'censored')} ps"
    )
    print(
        " edge timing:   "
        f"t_lat={dict(edge_timing.get('latency', {})).get('t_lat_ps', 'censored')} ps, "
        f"t_rec={dict(dict(edge_timing.get('recovery', {})).get('selected', {})).get('t_rec_ps', 'censored')} ps"
    )
    print("Figures")
    for key, path in saved.items():
        print(f" {key}: {path}")
    print(f" manifest: {manifest_path}")
    print("Status: OK")
    return 0


def _write_manifest(
    *,
    args: argparse.Namespace,
    raw_pre: Path,
    center_raw: Path,
    edge_raw: Path,
    output_dir: Path,
    saved: Mapping[str, Path],
    delta0_meV: float,
    center_history: Mapping[str, Any],
    edge_history: Mapping[str, Any],
    center_snapshots: Mapping[str, Any],
    edge_snapshots: Mapping[str, Any],
    center_timing: Mapping[str, Any],
    edge_timing: Mapping[str, Any],
) -> Path:
    manifest = {
        "schema_version": 1,
        "pipeline": "plot_pipelines/E3_photon_position_comparison.py",
        "purpose": "Matched center/edge photon-impact field and circuit-response PDFs.",
        "pre_run_name": str(args.pre_run_name),
        "center_run_name": str(args.center_run_name),
        "edge_run_name": str(args.edge_run_name),
        "raw_pre": str(raw_pre),
        "center_raw_photon": str(center_raw),
        "edge_raw_photon": str(edge_raw),
        "output_dir": str(output_dir),
        "requested_times_ps": [float(value) for value in args.times_ps],
        "center_resolved_times_ps": _resolved_times(center_snapshots, args.times_ps),
        "edge_resolved_times_ps": _resolved_times(edge_snapshots, args.times_ps),
        "delta0_meV": float(delta0_meV),
        "center_photon_time_ps": _photon_time(center_history),
        "edge_photon_time_ps": _photon_time(edge_history),
        "center_timing": dict(center_timing),
        "edge_timing": dict(edge_timing),
        "figures": {key: str(path) for key, path in saved.items()},
    }
    path = output_dir / "E3_photon_position_manifest.yaml"
    with path.open("w", encoding="utf-8") as stream:
        yaml.safe_dump(manifest, stream, sort_keys=False, allow_unicode=True, default_flow_style=False)
    return path


def _read_delta0_meV(raw_pre: Path) -> float:
    summary = _read_yaml(_require_file(raw_pre / "usadel_dos_summary.yaml", "Usadel summary"))
    value = _find_numeric(summary, "delta0_meV")
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError("PRE metadata do not provide a positive delta0_meV.")
    return float(value)


def _find_numeric(value: Any, target: str) -> float:
    if isinstance(value, Mapping):
        if target in value:
            try:
                return float(value[target])
            except Exception:
                pass
        for nested in value.values():
            found = _find_numeric(nested, target)
            if np.isfinite(found):
                return found
    return np.nan


def _resolved_times(snapshots: Mapping[str, Any], requested: list[float]) -> list[float]:
    stored = np.asarray(snapshots.get("snapshot_t_ps", []), dtype=float)
    if stored.size == 0:
        return []
    return [float(stored[int(np.nanargmin(np.abs(stored - float(value))))]) for value in requested]


def _photon_time(history: Mapping[str, Any]) -> float | None:
    time = np.asarray(history.get("t_ps", []), dtype=float)
    applied = np.asarray(history.get("photon_applied", []), dtype=bool)
    if time.size == 0 or applied.size == 0:
        return None
    if applied.size != time.size:
        applied = np.resize(applied, time.size)
    indices = np.flatnonzero(applied)
    return float(time[indices[0]]) if indices.size else None


def _load_npz(path: Path) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=True) as data:
        return {key: np.asarray(data[key]) for key in data.files}


def _read_yaml(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as stream:
        value = yaml.safe_load(stream) or {}
    return value if isinstance(value, dict) else {}


def _require_file(path: Path, label: str) -> Path:
    if not path.exists():
        raise FileNotFoundError(f"Missing {label}: {path}")
    return path


if __name__ == "__main__":
    raise SystemExit(main())
