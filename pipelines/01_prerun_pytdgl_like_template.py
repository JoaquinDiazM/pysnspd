#!/usr/bin/env python3
"""pyTDGL-faithful PRE-run entry point.

This wrapper runs the standard pySNSPD PRE-run after replacing the rectangular
mesh generator with a pyTDGL-style meshpy/triangle generator.  It then writes a
full finite-volume sidecar containing the pyTDGL-like Mesh/EdgeMesh/Voronoi
arrays so later solver work can compare directly against pyTDGL's data model.
"""
from __future__ import annotations

import argparse
import runpy
import sys
from pathlib import Path

import yaml

from pysnspd.config import load_config
from pysnspd.io.manager import RunManager
from pysnspd.mesh.pytdgl_like import (
    build_pytdgl_like_mesh_summary,
    generate_rectangular_pytdgl_fvm_mesh_from_parameters,
    parameters_from_config,
    save_pytdgl_like_mesh_npz,
)
from pysnspd.plotting.pytdgl_mesh import plot_pytdgl_fvm_mesh


def _extract_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--config", required=True)
    parser.add_argument("--run-name", required=True)
    parser.add_argument("--jitter-fraction", type=float, default=0.0)
    parser.add_argument("--boundary-guard-layers", type=int, default=1)
    parser.add_argument("--pytdgl-min-angle", type=float, default=32.5)
    parser.add_argument("--pytdgl-smooth", type=int, default=None)
    parser.add_argument("--pytdgl-max-edge-length-m", type=float, default=None)
    return parser.parse_known_args(argv)[0]


def main() -> int:
    args = _extract_args(sys.argv[1:])

    # First run the official PRE-run pipeline.  It will call
    # generate_rectangular_delaunay_mesh(), which this branch routes through the
    # pyTDGL-like meshpy generator.
    runpy.run_path(str(Path(__file__).with_name("01_prerun_template.py")), run_name="__main__")

    cfg = load_config(args.config)
    manager = RunManager.from_config(cfg, run_name=args.run_name)
    raw_pre = manager.raw_dir("pre")
    plots_mesh = manager.plots_dir("mesh")

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
