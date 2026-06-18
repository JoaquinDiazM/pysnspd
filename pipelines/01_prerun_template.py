"""
PRE-run template for pySNSPD.

OE2 implementation:
- Load project configuration.
- Create the run folder layout.
- Generate a protected rectangular mesh.
- Extract edges and boundary tags.
- Save mesh/edge arrays and summaries.
- Save diagnostic mesh plots.

No Usadel, phase-space catalog or gTDGL calculation is performed here yet.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import yaml

from pysnspd.config import load_config, validate_config
from pysnspd.io.manager import create_run_layout, write_manifest
from pysnspd.mesh.delaunay import (
    generate_rectangular_delaunay_mesh,
    mesh_summary,
    save_mesh_npz,
)
from pysnspd.mesh.edges import (
    assert_edge_data_consistent,
    build_edge_data,
    edge_summary,
    save_edges_npz,
)
from pysnspd.plotting.figures import (
    plot_boundary_tags,
    plot_mesh_geometry,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate OE2 mesh, edge data and diagnostic plots."
    )
    parser.add_argument(
        "--config",
        required=True,
        help="Path to YAML project configuration.",
    )
    parser.add_argument(
        "--run-name",
        default=None,
        help="Optional run name. If omitted, project.default_run_name is used.",
    )
    parser.add_argument(
        "--jitter-fraction",
        type=float,
        default=0.20,
        help="Interior mesh jitter as fraction of nominal spacing.",
    )
    parser.add_argument(
        "--boundary-guard-layers",
        type=int,
        default=1,
        help=(
            "Number of grid layers near each boundary kept unjittered. "
            "Use 1 as the default for stable boundary operators."
        ),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    cfg = load_config(args.config)
    cfg = validate_config(cfg)

    layout = create_run_layout(cfg, args.run_name)
    run_name = layout["run_name"]

    raw_pre = Path(layout["raw_pre"])
    plots_mesh = Path(layout["plots_mesh"])

    mesh = generate_rectangular_delaunay_mesh(
        cfg,
        jitter_fraction=args.jitter_fraction,
        boundary_guard_layers=args.boundary_guard_layers,
    )

    edge_data = build_edge_data(
        mesh.nodes,
        mesh.triangles,
        length_m=mesh.length_m,
        width_m=mesh.width_m,
    )
    assert_edge_data_consistent(edge_data)

    mesh_npz = save_mesh_npz(mesh, raw_pre / "mesh.npz")
    edges_npz = save_edges_npz(edge_data, raw_pre / "edges.npz")

    summary = {
        "run_name": run_name,
        "mesh": mesh_summary(mesh),
        "edges": edge_summary(edge_data),
    }

    summary_path = raw_pre / "mesh_summary.yaml"
    with summary_path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(
            summary,
            f,
            sort_keys=False,
            allow_unicode=True,
            default_flow_style=False,
        )

    mesh_plot = plot_mesh_geometry(
        mesh,
        edge_data,
        plots_mesh / "mesh_nodes_edges.png",
    )
    tags_plot = plot_boundary_tags(
        mesh,
        edge_data,
        plots_mesh / "mesh_boundary_tags.png",
    )

    manifest_path = write_manifest(
        cfg,
        run_name,
        stage="pre",
        extra={
            "pipeline": "01_prerun_template.py",
            "purpose": "OE2 mesh, edges, boundary tags and diagnostic plots",
            "outputs": {
                "mesh_npz": str(mesh_npz),
                "edges_npz": str(edges_npz),
                "mesh_summary": str(summary_path),
                "mesh_plot": str(mesh_plot),
                "boundary_tags_plot": str(tags_plot),
            },
            "summary": summary,
        },
    )

    print("PRE-run mesh generation")
    print(f"run_name              : {run_name}")
    print(f"raw_pre               : {raw_pre}")
    print(f"plots_mesh            : {plots_mesh}")
    print()
    print("Mesh summary")
    for key, value in summary["mesh"].items():
        print(f"  {key}: {value}")
    print()
    print("Edge summary")
    for key, value in summary["edges"].items():
        print(f"  {key}: {value}")
    print()
    print("Outputs")
    print(f"  mesh_npz            : {mesh_npz}")
    print(f"  edges_npz           : {edges_npz}")
    print(f"  mesh_summary        : {summary_path}")
    print(f"  mesh_plot           : {mesh_plot}")
    print(f"  boundary_tags_plot  : {tags_plot}")
    print(f"  pre_manifest        : {manifest_path}")
    print("Status: OK")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())