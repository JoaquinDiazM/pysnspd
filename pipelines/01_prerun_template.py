"""
PRE-run template for pySNSPD.

Current implementation:
- OE2: generate protected mesh, edges, boundary tags and mesh plots.
- OE3: generate first Usadel/DOS catalogue and diagnostic DOS plot.

No phase-space catalogue, gTDGL or photon dynamics is performed here yet.
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
    plot_phase_space_slices,
    plot_usadel_calibration_sweep,
    plot_usadel_dos_slices,
)
from pysnspd.usadel.catalog import (
    build_usadel_catalog_from_config,
    catalog_summary,
    save_usadel_catalog_npz,
)
from pysnspd.kinetic.phase_space import (
    build_phase_space_catalog_from_usadel_catalog,
    phase_space_summary,
    save_phase_space_catalog_npz,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate PRE-run mesh and first Usadel/DOS catalogue."
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
        help="Number of grid layers near each boundary kept unjittered.",
    )
    parser.add_argument(
        "--eta-fraction",
        type=float,
        default=1.0e-3,
        help="DOS numerical broadening as fraction of Delta0.",
    )
    parser.add_argument(
        "--gamma-max-fraction",
        type=float,
        default=0.80,
        help=(
            "Maximum depairing proxy Gamma_q as fraction of Delta0. "
            "The default 0.80 is chosen to cover current-crowding regions "
            "above the uniform critical current."
        ),
    )
    parser.add_argument(
        "--energy-max-factor",
        type=float,
        default=6.0,
        help="Maximum DOS energy as multiple of Delta0.",
    )
    parser.add_argument(
        "--phase-n-Te",
        type=int,
        default=6,
        help="Number of Te grid points for OE4 first-attempt phase-space catalogue.",
    )
    parser.add_argument(
        "--phase-n-delta",
        type=int,
        default=6,
        help="Number of Delta grid points selected from the Usadel catalogue.",
    )
    parser.add_argument(
        "--phase-n-q",
        type=int,
        default=6,
        help="Number of q grid points selected from the Usadel catalogue.",
    )
    parser.add_argument(
        "--phase-n-omega",
        type=int,
        default=160,
        help="Number of Omega grid points for OE4 first-attempt phase-space catalogue.",
    )
    parser.add_argument(
        "--phase-omega-max-meV",
        type=float,
        default=None,
        help=(
            "Maximum Omega axis for the OE4 phase-space catalogue in meV. "
            "If omitted, Omega_max defaults to the parent DOS E_max. "
            "For acoustic NbN tests, a typical value is 30 meV, while "
            "--energy-max-factor should be large enough to keep E+Omega "
            "inside the DOS catalogue."
        ),
    )
    parser.add_argument(
        "--phase-Te-min-K",
        type=float,
        default=None,
        help="Optional minimum Te for the phase-space catalogue.",
    )
    parser.add_argument(
        "--phase-Te-max-K",
        type=float,
        default=None,
        help="Optional maximum Te for the phase-space catalogue.",
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
    plots_diagnostics = Path(layout["plots_diagnostics"])

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

    mesh_edge_summary = {
        "run_name": run_name,
        "mesh": mesh_summary(mesh),
        "edges": edge_summary(edge_data),
    }

    mesh_summary_path = raw_pre / "mesh_summary.yaml"
    with mesh_summary_path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(
            mesh_edge_summary,
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

    usadel_catalog = build_usadel_catalog_from_config(
        cfg,
        eta_fraction=args.eta_fraction,
        gamma_max_fraction=args.gamma_max_fraction,
        energy_max_factor=args.energy_max_factor,
    )

    usadel_npz = save_usadel_catalog_npz(
        usadel_catalog,
        raw_pre / "usadel_dos_catalog.npz",
    )

    usadel_summary = catalog_summary(usadel_catalog)
    usadel_summary_path = raw_pre / "usadel_dos_summary.yaml"
    with usadel_summary_path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(
            {
                "run_name": run_name,
                "usadel": usadel_summary,
                "metadata": usadel_catalog.metadata,
            },
            f,
            sort_keys=False,
            allow_unicode=True,
            default_flow_style=False,
        )

    usadel_plot = plot_usadel_dos_slices(
        usadel_catalog,
        plots_diagnostics / "usadel_dos_slices.png",
    )
    calibration_plot = plot_usadel_calibration_sweep(
        usadel_catalog,
        plots_diagnostics / "usadel_calibration_sweep.png",
    )

    phase_catalog = build_phase_space_catalog_from_usadel_catalog(
        usadel_catalog,
        cfg,
        n_Te=args.phase_n_Te,
        n_delta=args.phase_n_delta,
        n_q=args.phase_n_q,
        n_omega=args.phase_n_omega,
        Te_min_K=args.phase_Te_min_K,
        Te_max_K=args.phase_Te_max_K,
        omega_max_meV=args.phase_omega_max_meV,
    )

    phase_npz = save_phase_space_catalog_npz(
        phase_catalog,
        raw_pre / "phase_space_catalog.npz",
    )

    phase_summary_data = phase_space_summary(phase_catalog)
    phase_summary_path = raw_pre / "phase_space_summary.yaml"

    with phase_summary_path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(
            {
                "run_name": run_name,
                "phase_space": phase_summary_data,
                "metadata": phase_catalog.metadata,
            },
            f,
            sort_keys=False,
            allow_unicode=True,
            default_flow_style=False,
        )

    phase_plot = plot_phase_space_slices(
        phase_catalog,
        plots_diagnostics / "phase_space_slices.png",
    )

    manifest_path = write_manifest(
        cfg,
        run_name,
        stage="pre",
        extra={
            "pipeline": "01_prerun_template.py",
            "purpose": "OE2 mesh plus OE3 Usadel/DOS catalogue plus OE4 phase-space catalogue",
            "outputs": {
                "mesh_npz": str(mesh_npz),
                "edges_npz": str(edges_npz),
                "mesh_summary": str(mesh_summary_path),
                "mesh_plot": str(mesh_plot),
                "boundary_tags_plot": str(tags_plot),
                "usadel_npz": str(usadel_npz),
                "usadel_summary": str(usadel_summary_path),
                "usadel_plot": str(usadel_plot),
                "calibration_plot": str(calibration_plot),
                "phase_space_npz": str(phase_npz),
                "phase_space_summary": str(phase_summary_path),
                "phase_space_plot": str(phase_plot),
            },
            "mesh_edge_summary": mesh_edge_summary,
            "usadel_summary": usadel_summary,
            "phase_space_summary": phase_summary_data,
        },
    )

    print("PRE-run generation")
    print(f"run_name              : {run_name}")
    print(f"raw_pre               : {raw_pre}")
    print(f"plots_mesh            : {plots_mesh}")
    print(f"plots_diagnostics     : {plots_diagnostics}")
    print()

    print("Mesh summary")
    for key, value in mesh_edge_summary["mesh"].items():
        print(f"  {key}: {value}")

    print()
    print("Edge summary")
    for key, value in mesh_edge_summary["edges"].items():
        print(f"  {key}: {value}")

    print()
    print("Usadel/DOS summary")
    for key, value in usadel_summary.items():
        print(f"  {key}: {value}")
    calibration_warnings = usadel_summary.get("calibration_warnings", [])
    print()
    if calibration_warnings:
        print("Calibration warnings")
        for warning in calibration_warnings:
            print(f"  WARNING: {warning}")
    else:
        print("Calibration warnings: none")

    print()
    print("Phase-space summary")
    for key, value in phase_summary_data.items():
        print(f"  {key}: {value}")

    print()
    print("Outputs")
    print(f"  mesh_npz            : {mesh_npz}")
    print(f"  edges_npz           : {edges_npz}")
    print(f"  mesh_summary        : {mesh_summary_path}")
    print(f"  mesh_plot           : {mesh_plot}")
    print(f"  boundary_tags_plot  : {tags_plot}")
    print(f"  usadel_npz          : {usadel_npz}")
    print(f"  usadel_summary      : {usadel_summary_path}")
    print(f"  usadel_plot         : {usadel_plot}")
    print(f"  calibration_plot    : {calibration_plot}")
    print(f"  pre_manifest        : {manifest_path}")
    print(f"  phase_space_npz    : {phase_npz}")
    print(f"  phase_space_summary: {phase_summary_path}")
    print(f"  phase_space_plot   : {phase_plot}")
    print("Status: OK")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())