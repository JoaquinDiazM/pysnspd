"""OE6 first attempt: analytic stationary seed for gTDGL relaxation.

This pipeline loads PRE-run outputs and writes an SS-stage seed:

    raw/<run_name>/ss/ss_seed.npz
    raw/<run_name>/ss/ss_seed_summary.yaml

It does not perform gTDGL time evolution. That belongs to OE7.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import yaml

from pysnspd.config import load_config, validate_config
from pysnspd.io.manager import create_run_layout, write_manifest
from pysnspd.mesh.delaunay import load_mesh_npz
from pysnspd.mesh.edges import load_edges_npz
from pysnspd.usadel.catalog import load_usadel_catalog_npz
from pysnspd.gtdgl.seed import (
    build_stationary_seed,
    save_stationary_seed_npz,
    seed_summary,
)
from pysnspd.plotting.ss_seed import (
    plot_ss_seed_boundary_currents,
    plot_ss_seed_current_density,
    plot_ss_seed_delta,
    plot_ss_seed_divergence,
    plot_ss_seed_phase,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="OE6 analytic stationary seed from PRE-run catalogues."
    )

    parser.add_argument("--config", required=True, help="Path to YAML config.")
    parser.add_argument(
        "--run-name",
        required=True,
        help="Run name where the SS seed will be written.",
    )
    parser.add_argument(
        "--pre-run-name",
        default=None,
        help=(
            "Run name containing PRE outputs. If omitted, use --run-name. "
            "Expected files are raw/<pre_run_name>/pre/mesh.npz, edges.npz "
            "and usadel_dos_catalog.npz."
        ),
    )
    parser.add_argument(
        "--I-bias-A",
        type=float,
        default=None,
        help="Override bias current. If omitted, use Usadel catalogue metadata.",
    )
    parser.add_argument(
        "--T-bias-K",
        type=float,
        default=None,
        help="Override bias temperature. If omitted, use Usadel catalogue metadata.",
    )
    parser.add_argument(
        "--phase-origin",
        choices=["center", "left"],
        default="center",
        help="Origin used in theta=q*(x-x0).",
    )

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
    usadel_path = raw_pre / "usadel_dos_catalog.npz"

    _require_existing(mesh_path)
    _require_existing(edges_path)
    _require_existing(usadel_path)

    mesh = load_mesh_npz(mesh_path)
    edge_data = load_edges_npz(edges_path)
    usadel_catalog = load_usadel_catalog_npz(usadel_path)

    seed = build_stationary_seed(
        mesh=mesh,
        edge_data=edge_data,
        usadel_catalog=usadel_catalog,
        I_bias_A=args.I_bias_A,
        T_bias_K=args.T_bias_K,
        phase_origin=args.phase_origin,
    )

    seed_npz = save_stationary_seed_npz(seed, raw_ss / "ss_seed.npz")
    summary = seed_summary(seed)

    summary_path = raw_ss / "ss_seed_summary.yaml"
    with summary_path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(
            {
                "run_name": run_name,
                "pre_run_name": pre_run_name,
                "seed": summary,
                "metadata": seed.metadata,
                "inputs": {
                    "mesh_npz": str(mesh_path),
                    "edges_npz": str(edges_path),
                    "usadel_npz": str(usadel_path),
                },
            },
            f,
            sort_keys=False,
            allow_unicode=True,
        )

    delta_plot = plot_ss_seed_delta(
        mesh,
        seed,
        plots_diag / "ss_seed_delta.png",
    )
    phase_plot = plot_ss_seed_phase(
        mesh,
        seed,
        plots_diag / "ss_seed_phase.png",
    )
    current_plot = plot_ss_seed_current_density(
        mesh,
        seed,
        plots_diag / "ss_seed_current_density.png",
    )
    div_plot = plot_ss_seed_divergence(
        mesh,
        seed,
        plots_diag / "ss_seed_divergence.png",
    )
    boundary_plot = plot_ss_seed_boundary_currents(
        seed,
        plots_diag / "ss_seed_boundary_currents.png",
    )

    manifest = write_manifest(
        cfg,
        run_name,
        stage="ss",
        extra={
            "pipeline": "02_oe6_ss_seed.py",
            "purpose": (
                "OE6 analytic stationary seed. No gTDGL relaxation is performed."
            ),
            "pre_run_name": pre_run_name,
            "inputs": {
                "mesh_npz": str(mesh_path),
                "edges_npz": str(edges_path),
                "usadel_npz": str(usadel_path),
            },
            "outputs": {
                "seed_npz": str(seed_npz),
                "seed_summary": str(summary_path),
                "delta_plot": str(delta_plot),
                "phase_plot": str(phase_plot),
                "current_plot": str(current_plot),
                "divergence_plot": str(div_plot),
                "boundary_currents_plot": str(boundary_plot),
            },
            "seed_summary": summary,
        },
    )

    print("OE6 analytic stationary seed")
    print(f"run_name                  : {run_name}")
    print(f"pre_run_name              : {pre_run_name}")
    print(f"raw_pre                   : {raw_pre}")
    print(f"raw_ss                    : {raw_ss}")
    print(f"plots_diagnostics         : {plots_diag}")
    print()
    print("Selected bias state")
    for key in [
        "I_bias_A",
        "Ic_A",
        "I_bias_over_Ic",
        "q_bias_m_inv",
        "q_critical_m_inv",
        "gamma_bias_meV",
        "delta_bias_meV",
        "current_density_bias_A_m2",
    ]:
        print(f"  {key}: {summary[key]}")
    print()
    print("Stationary diagnostics")
    for key in [
        "terminal_voltage_V",
        "integrated_left_current_A",
        "integrated_right_current_A",
        "net_boundary_current_A",
        "left_current_error_rel",
        "right_current_error_rel",
        "divergence_rms_A_m3",
        "phase_span_rad",
    ]:
        print(f"  {key}: {summary[key]}")
    print()
    print("Outputs")
    print(f"  seed_npz                 : {seed_npz}")
    print(f"  seed_summary             : {summary_path}")
    print(f"  delta_plot               : {delta_plot}")
    print(f"  phase_plot               : {phase_plot}")
    print(f"  current_plot             : {current_plot}")
    print(f"  divergence_plot          : {div_plot}")
    print(f"  boundary_currents_plot   : {boundary_plot}")
    print(f"  manifest                 : {manifest}")
    print("Status: OK")

    return 0


def _require_existing(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required OE6 input does not exist: {path}")
    if not path.is_file():
        raise FileNotFoundError(f"Required OE6 input is not a file: {path}")


if __name__ == "__main__":
    raise SystemExit(main())