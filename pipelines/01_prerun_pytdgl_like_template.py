#!/usr/bin/env python3
"""pyTDGL-faithful PRE-run entry point in SI units.

This wrapper delegates the heavy PRE computation to the standard pySNSPD
``01_prerun_template.py`` after removing the extra pyTDGL-only CLI flags from
``sys.argv``.  After the standard PRE-run finishes, it writes the full
pyTDGL-like finite-volume sidecar: Mesh, EdgeMesh and Voronoi arrays.

The meshing order follows pyTDGL's ``Device.make_mesh`` pattern:

    generate_mesh(...) -> Mesh.from_triangulation(..., create_submesh=False)
    -> Mesh.smooth(..., create_submesh=False)
    -> Mesh.from_triangulation(..., create_submesh=True)

Coordinates remain in SI meters; no pyTDGL coherence-length normalization is
applied here.
"""
from __future__ import annotations

import argparse
import runpy
import sys
from pathlib import Path
from typing import Any, Mapping

import yaml

from pysnspd.config import load_config
from pysnspd.mesh.pytdgl_like import (
    build_pytdgl_like_mesh_summary,
    generate_rectangular_pytdgl_fvm_mesh_from_parameters,
    parameters_from_config,
    save_pytdgl_like_mesh_npz,
)
from pysnspd.plotting.pytdgl_mesh import plot_pytdgl_fvm_mesh


def _parse_wrapper_args(argv: list[str]) -> tuple[argparse.Namespace, list[str]]:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--config", required=True)
    parser.add_argument("--run-name", required=True)
    parser.add_argument("--jitter-fraction", type=float, default=0.0)
    parser.add_argument("--boundary-guard-layers", type=int, default=1)
    parser.add_argument("--pytdgl-min-angle", type=float, default=32.5)
    parser.add_argument("--pytdgl-smooth", type=int, default=None)
    parser.add_argument("--pytdgl-max-edge-length-m", type=float, default=None)
    return parser.parse_known_args(argv)


def _config_section(config: Mapping[str, Any], name: str) -> Mapping[str, Any]:
    section = config.get(name, {})
    if not isinstance(section, Mapping):
        raise TypeError(f"config[{name!r}] must be a mapping.")
    return section


def _run_dirs_from_config(config: Mapping[str, Any], run_name: str) -> tuple[Path, Path]:
    paths = _config_section(config, "paths")
    root_value = paths.get("big_data_root")
    if root_value is None:
        raise KeyError("Missing config paths.big_data_root.")
    root = Path(str(root_value)).expanduser()
    raw_pre = root / "raw" / run_name / "pre"
    plots_mesh = root / "plots" / run_name / "mesh"
    raw_pre.mkdir(parents=True, exist_ok=True)
    plots_mesh.mkdir(parents=True, exist_ok=True)
    return raw_pre, plots_mesh


def _run_standard_prerun(
    *,
    standard_pipeline: Path,
    args: argparse.Namespace,
    passthrough: list[str],
) -> None:
    """Run the standard PRE pipeline without exposing pyTDGL-only flags to it."""
    forwarded_argv = [
        str(standard_pipeline),
        "--config",
        str(args.config),
        "--run-name",
        str(args.run_name),
        "--jitter-fraction",
        str(args.jitter_fraction),
        "--boundary-guard-layers",
        str(args.boundary_guard_layers),
        *passthrough,
    ]
    old_argv = sys.argv[:]
    try:
        sys.argv = forwarded_argv
        try:
            runpy.run_path(str(standard_pipeline), run_name="__main__")
        except SystemExit as exc:
            if exc.code not in (None, 0):
                raise
    finally:
        sys.argv = old_argv


def main() -> int:
    args, passthrough = _parse_wrapper_args(sys.argv[1:])
    standard_pipeline = Path(__file__).with_name("01_prerun_template.py")

    _run_standard_prerun(
        standard_pipeline=standard_pipeline,
        args=args,
        passthrough=passthrough,
    )

    cfg = load_config(args.config)
    raw_pre, plots_mesh = _run_dirs_from_config(cfg, args.run_name)

    params = parameters_from_config(
        cfg,
        jitter_fraction=args.jitter_fraction,
        boundary_guard_layers=args.boundary_guard_layers,
        max_edge_length_m=args.pytdgl_max_edge_length_m,
        min_angle_deg=args.pytdgl_min_angle,
        smooth=args.pytdgl_smooth,
    )
    fvm_mesh = generate_rectangular_pytdgl_fvm_mesh_from_parameters(params)
    sidecar = save_pytdgl_like_mesh_npz(fvm_mesh, raw_pre / "pytdgl_fvm_mesh.npz")
    summary = build_pytdgl_like_mesh_summary(fvm_mesh)
    summary.update(
        {
            "max_edge_length_m": float(params.max_edge_length_m),
            "min_angle_deg": float(params.min_angle_deg),
            "smooth": int(params.smooth),
            "units_policy": "SI meters; no pyTDGL coherence-length normalization in PRE mesh.",
        }
    )
    summary_path = raw_pre / "pytdgl_fvm_mesh_summary.yaml"
    with summary_path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(summary, f, sort_keys=False)
    plot_path = plot_pytdgl_fvm_mesh(fvm_mesh, plots_mesh / "pytdgl_fvm_mesh.png")

    print("pyTDGL-like finite-volume sidecar")
    print(f"  pytdgl_fvm_mesh_npz    : {sidecar}")
    print(f"  pytdgl_fvm_mesh_summary: {summary_path}")
    print(f"  pytdgl_fvm_mesh_plot   : {plot_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
