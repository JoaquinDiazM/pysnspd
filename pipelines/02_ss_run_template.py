"""SS-run template: OE6 analytic seed plus OE7 stationary gTDGL/Poisson."""
from __future__ import annotations

import argparse
from pathlib import Path

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
from pysnspd.gtdgl.material import build_gtdgl_material
from pysnspd.gtdgl.operators import build_fv_operators
from pysnspd.gtdgl.relax import (
    relax_stationary_gtdgl,
    save_relaxation_history_npz,
    save_stationary_state_npz,
)

from pysnspd.plotting.ss_seed import (
    plot_ss_seed_boundary_currents,
    plot_ss_seed_current_density,
    plot_ss_seed_delta,
    plot_ss_seed_divergence,
    plot_ss_seed_phase,
)
from pysnspd.plotting.ss_run import (
    plot_ss_available_snapshots,
    plot_ss_boundary_currents,
    plot_ss_relaxation_history,
    plot_ss_transport_current_profile,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="SS-run template: analytic seed plus stationary gTDGL/Poisson relaxation."
    )

    parser.add_argument("--config", required=True, help="Path to YAML config.")
    parser.add_argument("--run-name", required=True, help="Run name where SS outputs are written.")
    parser.add_argument(
        "--pre-run-name",
        default=None,
        help="Run name containing PRE outputs. If omitted, use --run-name.",
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
        help="Origin used in theta=q*(x-x0) for the analytic seed.",
    )
    parser.add_argument(
        "--ss-steps",
        type=int,
        default=2000,
        help="Maximum number of stationary gTDGL/Poisson steps.",
    )
    parser.add_argument(
        "--ss-min-steps",
        type=int,
        default=10,
        help="Minimum accepted steps before allowing convergence.",
    )
    parser.add_argument(
        "--ss-dt-fs",
        type=float,
        default=0.25,
        help="Initial SS gTDGL time step in femtoseconds.",
    )
    parser.add_argument(
        "--ss-tau-scale",
        type=float,
        default=0.10,
        help=(
            "Scale tau_ee(Tc) and tau_ep(Tc) during SS relaxation. "
            "This accelerates approach to the same stationary branch."
        ),
    )
    parser.add_argument(
        "--ss-tolerance-eta",
        type=float,
        default=1.0e-9,
        help="Convergence tolerance for max relative |Delta|^2 change.",
    )
    parser.add_argument(
        "--ss-tolerance-current-residual",
        type=float,
        default=1.0e-6,
        help="Convergence tolerance for dimensionless div(j) residual.",
    )
    parser.add_argument(
        "--ss-no-adapt-dt",
        action="store_true",
        help="Disable simple adaptive time-step rejection/growth.",
    )
    parser.add_argument(
        "--ss-unlock-terminals",
        action="store_true",
        help="Do not impose stationary terminal boundary conditions.",
    )
    parser.add_argument(
        "--ss-delta-boundary-policy",
        choices=["current_inversion", "vacuum_only", "normal_terminal", "none"],
        default="current_inversion",
        help=(
            "Delta boundary policy for OE7 diagnostics. current_inversion is "
            "the current OE7 behavior; vacuum_only removes terminal Delta "
            "forcing; normal_terminal imposes psi=0 on left/right, closer "
            "to pyTDGL terminal diagnostics."
        ),
    )
    parser.add_argument(
        "--ss-poisson-terminal-policy",
        choices=["target_flux", "zero_flux"],
        default="target_flux",
        help=(
            "Poisson terminal-current policy. target_flux is the current OE7 "
            "behavior; zero_flux removes the external terminal RHS flux."
        ),
    )
    parser.add_argument(
        "--ss-n-phi-snapshots",
        type=int,
        default=6,
        help="Number of electrostatic-potential snapshots to save and plot.",
    )
    parser.add_argument(
        "--ss-n-snapshots",
        type=int,
        default=None,
        help=(
            "Alias for --ss-n-phi-snapshots. This controls the shared snapshot "
            "times used for phi, |Delta|, phase, current density, divergence, "
            "supercurrent, normal-current, divergence, and pair-breaking diagnostics."
        ),
    )
    parser.add_argument(
        "--ss-progress",
        action="store_true",
        help="Show a progress bar during OE7 stationary relaxation.",
    )
    parser.add_argument("--dpi", type=int, default=480, help="DPI for diagnostic plots.")

    return parser.parse_args()


def main() -> int:
    args = parse_args()

    cfg = validate_config(load_config(args.config))
    run_name = args.run_name
    pre_run_name = args.pre_run_name or run_name
    n_phi_snapshots = args.ss_n_snapshots or args.ss_n_phi_snapshots

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
    seed_sum = seed_summary(seed)

    material = build_gtdgl_material(
        cfg,
        usadel_catalog,
        tau_scale=args.ss_tau_scale,
    )
    fv_ops = build_fv_operators(mesh, edge_data)

    result = relax_stationary_gtdgl(
        mesh=mesh,
        edge_data=edge_data,
        seed=seed,
        material=material,
        ops=fv_ops,
        steps=args.ss_steps,
        min_steps=args.ss_min_steps,
        dt_s=args.ss_dt_fs * 1.0e-15,
        tolerance_eta=args.ss_tolerance_eta,
        tolerance_current_residual=args.ss_tolerance_current_residual,
        adapt_dt=not args.ss_no_adapt_dt,
        lock_terminals=not args.ss_unlock_terminals,
        delta_boundary_policy=args.ss_delta_boundary_policy,
        poisson_terminal_policy=args.ss_poisson_terminal_policy,
        progress=args.ss_progress,
        n_phi_snapshots=n_phi_snapshots,
    )

    state_npz = save_stationary_state_npz(
        result.state,
        raw_ss / "ss_state_relaxed.npz",
    )
    history_npz = save_relaxation_history_npz(
        result.history,
        raw_ss / "ss_relaxation_history.npz",
    )

    summary_path = raw_ss / "ss_run_summary.yaml"
    summary_doc = {
        "run_name": run_name,
        "pre_run_name": pre_run_name,
        "stage": "ss",
        "seed": seed_sum,
        "stationary": result.summary,
        "metadata": {
            "backend": "oe7_ss_run_template_snapshot_diagnostics_v2",
            "description": (
                "SS-run template with OE6 analytic seed followed by OE7 "
                "frozen-temperature stationary gTDGL/Poisson relaxation."
            ),
            "thermal_policy": "frozen_Te_Tph",
            "circuit_policy": "inactive",
            "tau_policy": "tau_ee and tau_ep multiplied by --ss-tau-scale during SS only",
            "boundary_policy": args.ss_delta_boundary_policy,
            "poisson_terminal_policy": args.ss_poisson_terminal_policy,
            "diagnostic_plot_policy": (
                "Seed diagnostics are still emitted. Relaxed final-field "
                "colormaps are not emitted. Field diagnostics are emitted as "
                "snapshot mosaics when the corresponding arrays exist in history."
            ),
            "inputs": {
                "mesh_npz": str(mesh_path),
                "edges_npz": str(edges_path),
                "usadel_npz": str(usadel_path),
            },
        },
    }

    with summary_path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(summary_doc, f, sort_keys=False, allow_unicode=True)

    seed_delta_plot = plot_ss_seed_delta(
        mesh,
        seed,
        plots_diag / "ss_seed_delta.png",
        dpi=args.dpi,
    )
    seed_phase_plot = plot_ss_seed_phase(
        mesh,
        seed,
        plots_diag / "ss_seed_phase.png",
        dpi=args.dpi,
    )
    seed_current_plot = plot_ss_seed_current_density(
        mesh,
        seed,
        plots_diag / "ss_seed_current_density.png",
        dpi=args.dpi,
    )
    seed_div_plot = plot_ss_seed_divergence(
        mesh,
        seed,
        plots_diag / "ss_seed_divergence.png",
        dpi=args.dpi,
    )
    seed_boundary_plot = plot_ss_seed_boundary_currents(
        seed,
        plots_diag / "ss_seed_boundary_currents.png",
        dpi=args.dpi,
    )

    snapshot_plots = plot_ss_available_snapshots(
        mesh,
        result.history,
        plots_diag,
        dpi=args.dpi,
        ncols=3,
    )

    relaxed_boundary_plot = plot_ss_boundary_currents(
        result.summary,
        plots_diag / "ss_relaxed_boundary_currents.png",
        dpi=args.dpi,
    )

    transport_profile_plot = plot_ss_transport_current_profile(
        mesh=mesh,
        ops=fv_ops,
        state=result.state,
        output_path=plots_diag / "ss_transport_current_profile.png",
        target_current_A=seed_sum["I_bias_A"],
        thickness_m=material.thickness_m,
        dpi=args.dpi,
    )

    history_plot = plot_ss_relaxation_history(
        result.history,
        plots_diag / "ss_relaxation_history.png",
        dpi=args.dpi,
    )

    manifest = write_manifest(
        cfg,
        run_name,
        stage="ss",
        extra={
            "pipeline": "02_ss_run_template.py",
            "purpose": "OE7 stationary gTDGL/Poisson relaxation from analytic seed.",
            "pre_run_name": pre_run_name,
            "inputs": {
                "mesh_npz": str(mesh_path),
                "edges_npz": str(edges_path),
                "usadel_npz": str(usadel_path),
            },
            "outputs": {
                "seed_npz": str(seed_npz),
                "state_npz": str(state_npz),
                "history_npz": str(history_npz),
                "summary_yaml": str(summary_path),
                "seed_delta_plot": str(seed_delta_plot),
                "seed_phase_plot": str(seed_phase_plot),
                "seed_current_plot": str(seed_current_plot),
                "seed_divergence_plot": str(seed_div_plot),
                "seed_boundary_currents_plot": str(seed_boundary_plot),
                "snapshot_plots": {
                    key: str(path) for key, path in snapshot_plots.items()
                },
                "relaxed_boundary_currents_plot": str(relaxed_boundary_plot),
                "transport_profile_plot": str(transport_profile_plot),
                "history_plot": str(history_plot),
            },
            "seed_summary": seed_sum,
            "stationary_summary": result.summary,
        },
    )

    print("SS-run template: OE6 seed + OE7 stationary gTDGL/Poisson")
    print(f"run_name             : {run_name}")
    print(f"pre_run_name         : {pre_run_name}")
    print(f"raw_pre              : {raw_pre}")
    print(f"raw_ss               : {raw_ss}")
    print(f"plots_diagnostics    : {plots_diag}")
    print()

    print("Seed")
    for key in (
        "I_bias_A",
        "Ic_A",
        "I_bias_over_Ic",
        "q_bias_m_inv",
        "q_critical_m_inv",
        "delta_bias_meV",
        "terminal_voltage_V",
    ):
        print(f"  {key}: {seed_sum[key]}")

    print()
    print("Stationary relaxation")
    for key in (
        "converged",
        "accepted_steps",
        "rejected_steps",
        "final_time_ps",
        "tau_scale",
        "tau_ee_Tc_effective_ps",
        "tau_ep_Tc_effective_ps",
        "terminal_voltage_V",
        "normal_ohmic_voltage_V",
        "terminal_voltage_over_normal",
        "normal_current_fraction_rms",
        "normal_current_fraction_max",
        "normal_current_max_A_m2",
        "total_current_max_A_m2",
        "current_residual",
        "eta_R_final",
        "divergence_rms_A_m3",
        "min_delta_over_delta0",
        "mean_delta_over_delta0",
        "max_pairbreaking_ratio",
        "p99_pairbreaking_ratio",
        "edge_Q_max_m_inv",
    ):
        print(f"  {key}: {result.summary[key]}")

    print()
    print("Boundary currents [A]")
    for key, value in result.summary["boundary_currents_A"].items():
        print(f"  {key}: {value}")

    print()
    print("Snapshot plots")
    if snapshot_plots:
        for key, path in snapshot_plots.items():
            print(f"  {key}: {path}")
    else:
        print("  none")

    print()
    print("Outputs")
    print(f"  seed_npz        : {seed_npz}")
    print(f"  state_npz       : {state_npz}")
    print(f"  history_npz     : {history_npz}")
    print(f"  summary_yaml    : {summary_path}")
    print(f"  manifest        : {manifest}")
    print("Status: OK")

    return 0


def _require_existing(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required SS-run input does not exist: {path}")
    if not path.is_file():
        raise FileNotFoundError(f"Required SS-run input is not a file: {path}")


if __name__ == "__main__":
    raise SystemExit(main())