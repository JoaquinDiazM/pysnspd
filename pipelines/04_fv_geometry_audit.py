#!/usr/bin/env python3
"""OE7 finite-volume geometry audit.

This diagnostic does not run the stationary relaxation.  It loads the PRE mesh,
builds the current FV operators, and writes one plot plus one YAML report that
compare the current operator control volumes against an independent
circumcentric/Voronoi reference.
"""
from __future__ import annotations

import argparse
from pathlib import Path

from pysnspd.config import load_config, validate_config
from pysnspd.io.manager import create_run_layout
from pysnspd.mesh.delaunay import load_mesh_npz
from pysnspd.mesh.edges import load_edges_npz
from pysnspd.gtdgl.operators import build_fv_operators
from pysnspd.plotting.fv_geometry import (
    compute_fv_geometry_audit,
    plot_fv_geometry_audit,
    write_fv_geometry_audit_yaml,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Audit OE7 finite-volume geometry against a circumcentric/Voronoi reference."
    )
    parser.add_argument("--config", required=True, help="Path to YAML config.")
    parser.add_argument("--run-name", required=True, help="Run name where diagnostics are written.")
    parser.add_argument(
        "--pre-run-name",
        default=None,
        help="Run name containing PRE outputs. If omitted, use --run-name.",
    )
    parser.add_argument("--dpi", type=int, default=480, help="DPI for diagnostic plots.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    cfg = validate_config(load_config(args.config))
    run_name = args.run_name
    pre_run_name = args.pre_run_name or run_name

    layout = create_run_layout(cfg, run_name)
    raw_ss = Path(layout["raw_ss"])
    plots_diag = Path(layout["plots_diagnostics"])

    big_root = Path(cfg["project"]["big_data_root"]).expanduser().resolve()
    raw_pre = big_root / "raw" / pre_run_name / "pre"
    mesh_path = raw_pre / "mesh.npz"
    edges_path = raw_pre / "edges.npz"
    _require_existing(mesh_path)
    _require_existing(edges_path)

    mesh = load_mesh_npz(mesh_path)
    edge_data = load_edges_npz(edges_path)
    fv_ops = build_fv_operators(mesh, edge_data)

    audit = compute_fv_geometry_audit(mesh, edge_data=edge_data, ops=fv_ops)
    plot_path = plot_fv_geometry_audit(
        mesh,
        edge_data=edge_data,
        ops=fv_ops,
        output_path=plots_diag / "fv_geometry_audit.png",
        dpi=args.dpi,
    )
    report_path = write_fv_geometry_audit_yaml(audit, raw_ss / "fv_geometry_audit.yaml")

    print("FV geometry audit: done")
    print(f"Plot:   {plot_path}")
    print(f"Report: {report_path}")
    print(
        "negative_reference_area_fraction: "
        f"{audit['summary']['negative_reference_area_fraction']:.6g}"
    )
    ratio = audit["summary"]["reference_to_current_area_ratio_signed"]
    print(
        "reference/current area ratio percentiles: "
        f"p05={ratio['p05']:.6g}, p50={ratio['p50']:.6g}, p95={ratio['p95']:.6g}"
    )
    scale = audit["summary"]["laplace_scale_current_over_median"]
    print(
        "current Laplacian scale / median percentiles: "
        f"p05={scale['p05']:.6g}, p50={scale['p50']:.6g}, p95={scale['p95']:.6g}"
    )
    return 0


def _require_existing(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required input does not exist: {path}")


if __name__ == "__main__":
    raise SystemExit(main())
