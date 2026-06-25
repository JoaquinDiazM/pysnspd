"""Audit PRE-run mesh quality for OE7 gTDGL/Poisson stability."""
from __future__ import annotations

import argparse
from pathlib import Path

from pysnspd.config import load_config, validate_config
from pysnspd.mesh.delaunay import load_mesh_npz
from pysnspd.mesh.edges import load_edges_npz
from pysnspd.gtdgl.operators import build_fv_operators
from pysnspd.mesh.quality import (
    assert_mesh_quality,
    build_mesh_quality_report,
    save_mesh_quality_report,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit FV mesh quality for OE7.")
    parser.add_argument("--config", required=True, help="Path to YAML project config.")
    parser.add_argument("--run-name", required=True, help="Run name containing raw/<run>/pre/mesh.npz.")
    parser.add_argument("--pre-run-name", default=None, help="Alias for --run-name if supplied.")
    parser.add_argument("--fail-on-bad-mesh", action="store_true", help="Raise if quality status is fail.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    cfg = validate_config(load_config(args.config))
    run_name = args.pre_run_name or args.run_name
    big_root = Path(cfg["project"]["big_data_root"]).expanduser().resolve()
    raw_pre = big_root / "raw" / run_name / "pre"
    mesh_path = raw_pre / "mesh.npz"
    edges_path = raw_pre / "edges.npz"

    if not mesh_path.is_file():
        raise FileNotFoundError(mesh_path)
    if not edges_path.is_file():
        raise FileNotFoundError(edges_path)

    mesh = load_mesh_npz(mesh_path)
    edge_data = load_edges_npz(edges_path)
    ops = build_fv_operators(mesh, edge_data)
    report = build_mesh_quality_report(mesh, ops)
    out = save_mesh_quality_report(report, raw_pre / "mesh_quality_report.yaml")

    print(f"Mesh-quality status: {report.status}")
    print(f"Report: {out}")
    for msg in report.warnings:
        print(msg)
    for msg in report.recommendations:
        print("RECOMMENDATION:", msg)

    if args.fail_on_bad_mesh:
        assert_mesh_quality(report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
