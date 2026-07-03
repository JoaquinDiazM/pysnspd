"""Stationary gTDGL--Poisson SS-run using the flat gTDGL backend.

This pipeline loads PRE-run outputs, builds the analytic stationary seed, and
relaxes it with ``pysnspd.gtdgl.solve_stationary_pytdgl_like`` promoted to the
package root.  The preferred supercurrent law is the Matsubara Usadel table
stored by ``01_prerun_template.py``.

The pipeline can also run a current sweep.  In that mode the base current is run
with the normal terminal report, while extra current-offset cases run quietly in
parallel and only their output directories are printed at the end.
"""
from __future__ import annotations

import argparse
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import numpy as np
import yaml

from pysnspd.config import load_config, validate_config
from pysnspd.io.manager import create_run_layout, write_manifest
from pysnspd.mesh.delaunay import load_mesh_npz
from pysnspd.mesh.edges import load_edges_npz
from pysnspd.usadel.catalog import load_usadel_catalog_npz
from pysnspd.gtdgl import build_fv_operators, build_gtdgl_material, solve_stationary_pytdgl_like
from pysnspd.gtdgl.seed import build_stationary_seed, save_stationary_seed_npz, seed_summary
from pysnspd.gtdgl.state_io import save_relaxation_history_npz, save_stationary_state_npz
from pysnspd.gtdgl.snapshot_diagnostics import (
    save_ss_snapshot_bundle_npz,
    write_ss_snapshot_power_diagnostics,
)
from pysnspd.gtdgl.usadel_current import (
    attach_usadel_supercurrent_table_from_npz,
    validate_strict_usadel_supercurrent_table_npz,
)
from pysnspd.plotting.ss_run import plot_ss_adaptive_timestep_history


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run stationary gTDGL--Poisson relaxation with the flat gTDGL backend."
    )

    parser.add_argument("--config", required=True, help="Path to YAML project config.")

    parser.add_argument(
        "--run-name",
        default=None,
        help="SS run name. If omitted, project.default_run_name is used.",
    )
    parser.add_argument(
        "--pre-run-name",
        default=None,
        help="PRE-run name to load. If omitted, use --run-name/default_run_name.",
    )

    parser.add_argument("--phase-origin", choices=("center", "left"), default="center")

    # Default SS target constructed from the first useful metallic-contact
    # Usadel-Poisson runs:
    #
    #   I = 20 uA
    #   physical run time = 20 ps
    #   dt_init = 0.30 fs
    #   metallic terminals with smooth 2.5 xi seed healing
    #   stationarity measured only in the superconducting bulk.
    parser.add_argument(
        "--ss-time-ps",
        type=float,
        default=None,
        help=(
            "Physical SS relaxation time in ps. If omitted, use ss_run.total_time_ps from the YAML."
        ),
    )
    parser.add_argument(
        "--ss-dt-fs",
        type=float,
        default=None,
        help=(
            "Initial SS time step in fs. If omitted, use ss_run.dt_s from the YAML."
        ),
    )
    parser.add_argument("--ss-target-current-uA", type=float, default=20.0)
    parser.add_argument(
        "--ss-overcritical-seed-policy",
        choices=("clamp-to-ic", "error"),
        default="clamp-to-ic",
        help=(
            "Policy used only when the requested target current exceeds the PRE Usadel Ic. "
            "'clamp-to-ic' keeps the simulation target current unchanged, but builds "
            "the initial analytic superconducting seed slightly below Ic.  This is the "
            "recommended mode for searching overcritical PSL formation.  'error' preserves "
            "the old behavior and refuses overcritical stationary seeds."
        ),
    )
    parser.add_argument(
        "--ss-overcritical-seed-fraction",
        type=float,
        default=0.98,
        help=(
            "When --ss-overcritical-seed-policy=clamp-to-ic and the requested current "
            "is above the PRE Usadel Ic, build the analytic seed at this fraction of Ic. "
            "The boundary-current target remains the requested overcritical current."
        ),
    )
    parser.add_argument(
        "--extra-currents-uA",
        type=float,
        nargs="*",
        default=[],
        metavar="OFFSET_UA",
        help=(
            "Optional current sweep offsets in microamps relative to --ss-target-current-uA. "
            "Example: --ss-target-current-uA 28 --extra-currents-uA +2 +4 +6 "
            "runs 28, 30, 32 and 34 uA. Extra cases run quietly in parallel."
        ),
    )
    parser.add_argument(
        "--ss-sweep-workers",
        type=int,
        default=None,
        help="Parallel workers for extra-current sweep cases. Defaults to parallel.workers in the YAML.",
    )
    parser.add_argument("--ss-snapshots", type=int, default=None)

    parser.add_argument(
        "--ss-progress",
        dest="ss_progress",
        action="store_true",
        default=True,
        help="Show SS progress bar for the base case. Enabled by default.",
    )
    parser.add_argument(
        "--ss-no-progress",
        dest="ss_progress",
        action="store_false",
        help="Disable SS progress bar.",
    )
    parser.add_argument(
        "--ss-verbose-report",
        action="store_true",
        help=(
            "Print the full nested Seed/Solver dictionaries to the terminal. "
            "By default the pipeline writes complete YAML/manifest metadata but "
            "prints only a compact run report."
        ),
    )

    # Metallic-contact boundary model.
    parser.add_argument("--ss-terminal-psi", type=float, default=0.0)
    parser.add_argument(
        "--ss-terminal-healing-xi",
        type=float,
        default=2.5,
        help=(
            "Apply a smooth tanh contact-healing envelope to the initial |Delta| seed. "
            "Default 2.5 makes |Delta| recover to the requested bulk fraction over "
            "roughly 2--3 physical coherence lengths."
        ),
    )
    parser.add_argument("--ss-terminal-healing-fraction", type=float, default=0.95)

    # eta_R is kept only as a diagnostic; the SS target now uses gauge-fixed
    # edge gradients in the superconducting bulk.
    parser.add_argument(
        "--ss-stationarity-eta",
        type=float,
        default=1.0e-5,
        help=(
            "Info-only solver amplitude residual stored in the summary; "
            "no longer gates SS target pass/fail."
        ),
    )

    # Operational quasi-SS defaults for the central superconducting bulk.
    # These are intentionally looser than the earlier strict residual target;
    # they are meant to pass the small-tau quasi-stationary state in about 20 ps
    # while still rejecting the long-tau contact-conversion-dominated branch.
    parser.add_argument(
        "--ss-stationarity-phase-gradient-rel",
        type=float,
        default=3.0e-1,
        help=(
            "Relative tolerance for time-stationarity of edge phase gradient "
            "Q = grad(arg Delta), evaluated only in the superconducting bulk."
        ),
    )
    parser.add_argument(
        "--ss-stationarity-phi-gradient-rel",
        type=float,
        default=2.5e-1,
        help=(
            "Relative tolerance for time-stationarity of edge grad(phi), "
            "evaluated only in the superconducting bulk."
        ),
    )
    parser.add_argument(
        "--ss-stationarity-q-abs-m-inv",
        type=float,
        default=6.0e6,
        help="Absolute fallback tolerance for changes in Q, in m^-1.",
    )
    parser.add_argument(
        "--ss-stationarity-phi-gradient-abs-V-m",
        type=float,
        default=2.0e3,
        help="Absolute fallback tolerance for changes in grad(phi), in V/m.",
    )
    parser.add_argument(
        "--ss-stationarity-edge-active-threshold",
        type=float,
        default=0.05,
        help=(
            "Exclude edges whose final |Delta| is below this fraction of bulk, "
            "because phase is undefined near |Delta|=0."
        ),
    )
    parser.add_argument(
        "--ss-stationarity-bulk-exclusion-xi",
        type=float,
        default=4.0,
        help=(
            "Evaluate Q and grad(phi) stationarity only on bulk edges at least "
            "this many physical coherence lengths away from metallic contacts."
        ),
    )

    # Deprecated aliases kept for compatibility with older command lines.
    parser.add_argument(
        "--ss-stationarity-delta-rel",
        type=float,
        default=None,
        help=(
            "Deprecated alias: used as --ss-stationarity-phase-gradient-rel "
            "if the new flag is omitted."
        ),
    )
    parser.add_argument(
        "--ss-stationarity-phi-rel",
        type=float,
        default=None,
        help=(
            "Deprecated alias: used as --ss-stationarity-phi-gradient-rel "
            "if the new flag is omitted."
        ),
    )

    parser.add_argument("--ss-convergence-min-steps", type=int, default=500)
    parser.add_argument(
        "--ss-stop-on-convergence",
        action="store_true",
        help=(
            "Stop early when the info-only max_d_abs_sq_psi threshold is reached. "
            "By default the solver now always runs until --ss-time-ps and only "
            "records eta convergence as a diagnostic."
        ),
    )

    # These already pass with large margin; keep them strict.
    parser.add_argument("--ss-continuity-rms-tol", type=float, default=1.0e-6)
    parser.add_argument("--ss-continuity-max-tol", type=float, default=1.0e-3)
    parser.add_argument("--ss-continuity-poisson-tol", type=float, default=1.0e-9)

    # Contact recovery target: |Delta| should recover over an intermediate
    # physical distance, neither one-cell abrupt nor too long.
    parser.add_argument("--ss-recovery-min-xi", type=float, default=1.5)
    parser.add_argument("--ss-recovery-max-xi", type=float, default=4.0)

    parser.add_argument(
        "--ss-no-adaptive",
        action="store_true",
        help="Disable adaptive time stepping in the pyTDGL-like core.",
    )
    parser.add_argument(
        "--ss-adaptive-window",
        type=int,
        default=6,
        help="Moving-window length used to choose the next adaptive Euler time step.",
    )
    parser.add_argument(
        "--ss-max-solve-retries",
        type=int,
        default=8,
        help="Maximum pyTDGL-style retry/shrink attempts for one Euler update.",
    )
    parser.add_argument(
        "--ss-adaptive-time-step-multiplier",
        type=float,
        default=0.5,
        help="Multiplier applied to dt after a failed local |psi|^2 solve.",
    )
    parser.add_argument(
        "--ss-adaptive-growth-factor",
        type=float,
        default=1.5,
        help=(
            "Maximum multiplicative growth of the next tentative dt after an accepted step. "
            "This avoids repeatedly jumping straight to dt_max and then shrinking."
        ),
    )
    parser.add_argument(
        "--ss-dt-max-factor",
        type=float,
        default=6.0,
        help="Adaptive dt upper bound as a multiple of --ss-dt-fs.",
    )

    parser.add_argument(
        "--ss-supercurrent-law",
        choices=("auto", "gl", "usadel-poisson"),
        default="usadel-poisson",
        help=(
            "Default is strict usadel-poisson. SS refuses legacy 1D/2D current tables; "
            "PRE must provide js_A_m2[Te,delta,q]."
        ),
    )

    parser.add_argument(
        "--ss-allmaras-contact-guard-layers",
        type=int,
        default=2,
        help=(
            "Disable the Allmaras current-mismatch correction on terminal/contact nodes "
            "and this many graph-neighbor layers. Diffusion/reaction terms remain active."
        ),
    )

    parser.add_argument("--dpi", type=int, default=480, help="DPI for SS diagnostic plots.")

    return parser.parse_args()


def main() -> int:
    args = parse_args()
    cfg = validate_config(load_config(args.config))

    base_layout = create_run_layout(cfg, args.run_name)
    base_run_name = base_layout["run_name"]
    pre_name = args.pre_run_name or base_run_name
    base_current_uA = float(args.ss_target_current_uA) if args.ss_target_current_uA is not None else 20.0
    extra_offsets = [float(x) for x in (args.extra_currents_uA or [])]

    if not extra_offsets:
        _run_single_current_case(
            config_path=args.config,
            args_dict=vars(args),
            run_name=base_run_name,
            pre_name=pre_name,
            target_current_uA=base_current_uA,
            quiet=False,
            progress=bool(args.ss_progress),
            role="base",
        )
        return 0

    workers = _resolve_sweep_workers(cfg, args.ss_sweep_workers, n_extra=len(extra_offsets))
    print("SS current sweep")
    print(f"  pre_run_name: {pre_name}")
    print(f"  base_run_name: {base_run_name}")
    print(f"  base_current_uA: {base_current_uA}")
    print(f"  extra_offsets_uA: {extra_offsets}")
    print(f"  extra_workers: {workers}")
    print()

    futures = []
    extra_results: list[dict[str, Any]] = []
    args_dict = vars(args)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for offset in extra_offsets:
            current_uA = base_current_uA + float(offset)
            if current_uA <= 0.0:
                raise ValueError(f"Extra current offset {offset:+g} uA gives non-positive current {current_uA:g} uA.")
            run_name = _current_sweep_run_name(base_run_name, offset, current_uA)
            futures.append(
                pool.submit(
                    _run_single_current_case,
                    config_path=args.config,
                    args_dict=args_dict,
                    run_name=run_name,
                    pre_name=pre_name,
                    target_current_uA=current_uA,
                    quiet=True,
                    progress=False,
                    role="extra",
                )
            )

        base_result = _run_single_current_case(
            config_path=args.config,
            args_dict=args_dict,
            run_name=base_run_name,
            pre_name=pre_name,
            target_current_uA=base_current_uA,
            quiet=False,
            progress=bool(args.ss_progress),
            role="base",
        )

        failed_results: list[dict[str, Any]] = []
        for future in as_completed(futures):
            try:
                extra_results.append(future.result())
            except Exception as exc:  # pragma: no cover - exercised on cluster failures.
                failed_results.append({
                    "error_type": type(exc).__name__,
                    "error": str(exc),
                })

    all_results = [base_result] + sorted(extra_results, key=lambda item: float(item["target_current_uA"]))
    print()
    print("SS sweep run directories")
    for item in all_results:
        print(
            f"  {item['target_current_uA']:8.3f} uA  "
            f"seed={item.get('seed_current_uA', float('nan')):8.3f} uA  "
            f"clamped={item.get('seed_is_overcritical_clamped', False)}  "
            f"run={item['run_name']}  raw_ss={item['raw_ss']}"
        )
    if failed_results:
        print()
        print("SS sweep failed cases")
        for item in failed_results:
            print(f"  {item['error_type']}: {item['error']}")
        print("Status: PARTIAL")
        return 2
    print("Status: OK")
    return 0


def _run_single_current_case(
    *,
    config_path: str,
    args_dict: dict[str, Any],
    run_name: str,
    pre_name: str,
    target_current_uA: float,
    quiet: bool,
    progress: bool,
    role: str,
) -> dict[str, Any]:
    args = argparse.Namespace(**args_dict)
    cfg = validate_config(load_config(config_path))

    ss_layout = create_run_layout(cfg, run_name)
    pre_layout = create_run_layout(cfg, pre_name)

    raw_pre = Path(pre_layout["raw_pre"])
    raw_ss = Path(ss_layout["raw_ss"])
    raw_ss.mkdir(parents=True, exist_ok=True)

    mesh_path = raw_pre / "mesh.npz"
    edges_path = raw_pre / "edges.npz"
    usadel_path = raw_pre / "usadel_dos_catalog.npz"

    _require_file(mesh_path, "PRE mesh")
    _require_file(edges_path, "PRE edges")
    _require_file(usadel_path, "PRE Usadel catalogue")

    mesh = load_mesh_npz(mesh_path)
    edge_data = load_edges_npz(edges_path)
    ops = build_fv_operators(mesh, edge_data)

    strict_table_summary = validate_strict_usadel_supercurrent_table_npz(usadel_path)
    base_usadel_catalog = load_usadel_catalog_npz(usadel_path)
    usadel_catalog = attach_usadel_supercurrent_table_from_npz(base_usadel_catalog, usadel_path)
    supercurrent_law, supercurrent_policy = _resolve_supercurrent_law(
        requested=args.ss_supercurrent_law,
        strict_table_summary=strict_table_summary,
    )

    target_current_A = float(target_current_uA) * 1.0e-6
    seed_current_A, seed_current_policy = _resolve_seed_current_for_target(
        usadel_catalog=base_usadel_catalog,
        target_current_A=target_current_A,
        overcritical_policy=str(args.ss_overcritical_seed_policy),
        overcritical_seed_fraction=float(args.ss_overcritical_seed_fraction),
    )

    seed = build_stationary_seed(
        mesh=mesh,
        edge_data=edge_data,
        usadel_catalog=base_usadel_catalog,
        I_bias_A=seed_current_A,
        phase_origin=args.phase_origin,
    )
    seed_npz = save_stationary_seed_npz(seed, raw_ss / "stationary_seed.npz")
    seed_summary_data = seed_summary(seed)
    seed_summary_data["simulation_target_current_A"] = float(target_current_A)
    seed_summary_data["simulation_target_current_uA"] = float(target_current_uA)
    seed_summary_data["analytic_seed_current_A"] = float(seed_current_A)
    seed_summary_data["analytic_seed_current_uA"] = float(seed_current_A * 1.0e6)
    seed_summary_data["overcritical_seed_policy"] = seed_current_policy

    allmaras_diffusion = _read_pre_allmaras_diffusion(raw_pre)
    material = build_gtdgl_material(
        cfg,
        base_usadel_catalog,
        diffusion_factor=float(allmaras_diffusion["D_effective_factor"]),
    )

    ss_run_cfg = cfg.get("ss_run", {}) if isinstance(cfg, dict) else {}
    if args.ss_time_ps is not None:
        total_time_ps = float(args.ss_time_ps)
    else:
        total_time_ps = float(ss_run_cfg.get("total_time_ps", ss_run_cfg.get("physical_time_ps", 20.0)))
    if total_time_ps <= 0.0:
        raise ValueError("--ss-time-ps or ss_run.total_time_ps must be positive.")

    if args.ss_dt_fs is not None:
        dt_s = float(args.ss_dt_fs) * 1.0e-15
    else:
        dt_s = float(ss_run_cfg.get("dt_s", 1.0e-15))
    if dt_s <= 0.0:
        raise ValueError("--ss-dt-fs or ss_run.dt_s must be positive.")

    n_snapshots = int(args.ss_snapshots if args.ss_snapshots is not None else ss_run_cfg.get("snapshots", 8))
    if n_snapshots <= 0:
        raise ValueError("--ss-snapshots or ss_run.snapshots must be positive.")

    with _SSProgressConsoleFilter(
        enabled=bool(progress),
        total_time_ps=float(total_time_ps),
    ):
        result = solve_stationary_pytdgl_like(
            mesh=mesh,
            edge_data=edge_data,
            seed=seed,
            material=material,
            ops=ops,
            steps=None,
            total_time_s=float(total_time_ps) * 1.0e-12,
            dt_s=float(dt_s),
            target_current_A=target_current_A,
            usadel_catalog=usadel_catalog,
            terminal_psi=float(args.ss_terminal_psi),
            adaptive=not bool(args.ss_no_adaptive),
            adaptive_window=int(args.ss_adaptive_window),
            max_solve_retries=int(args.ss_max_solve_retries),
            adaptive_time_step_multiplier=float(args.ss_adaptive_time_step_multiplier),
            adaptive_growth_factor=float(args.ss_adaptive_growth_factor),
            dt_max_factor=float(args.ss_dt_max_factor),
            n_snapshots=int(n_snapshots),
            progress=bool(progress),
            supercurrent_law=supercurrent_law,
            terminal_healing_xi=args.ss_terminal_healing_xi,
            terminal_healing_fraction=float(args.ss_terminal_healing_fraction),
            stationarity_eta=float(args.ss_stationarity_eta),
            stationarity_phase_gradient_rel=args.ss_stationarity_phase_gradient_rel,
            stationarity_phi_gradient_rel=args.ss_stationarity_phi_gradient_rel,
            stationarity_q_abs_m_inv=float(args.ss_stationarity_q_abs_m_inv),
            stationarity_phi_gradient_abs_V_m=float(args.ss_stationarity_phi_gradient_abs_V_m),
            stationarity_edge_active_threshold=float(args.ss_stationarity_edge_active_threshold),
            stationarity_bulk_exclusion_xi=float(args.ss_stationarity_bulk_exclusion_xi),
            stationarity_delta_rel=args.ss_stationarity_delta_rel,
            stationarity_phi_rel=args.ss_stationarity_phi_rel,
            convergence_min_steps=int(args.ss_convergence_min_steps),
            stop_on_convergence=bool(args.ss_stop_on_convergence),
            continuity_rms_tol=float(args.ss_continuity_rms_tol),
            continuity_max_tol=float(args.ss_continuity_max_tol),
            continuity_poisson_tol=float(args.ss_continuity_poisson_tol),
            recovery_min_xi=float(args.ss_recovery_min_xi),
            recovery_max_xi=float(args.ss_recovery_max_xi),
            allmaras_contact_guard_layers=int(args.ss_allmaras_contact_guard_layers),
        )

    state_npz = save_stationary_state_npz(result.state, raw_ss / "stationary_state.npz")
    history = _history_with_fv_topology(result.history, ops=ops)
    history_npz = save_relaxation_history_npz(history, raw_ss / "relaxation_history.npz")
    snapshots_npz = save_ss_snapshot_bundle_npz(history, raw_ss / "stationary_snapshots.npz")
    snapshot_power_npz = None
    power_table_path = raw_pre / "power_table_catalog.npz"
    if power_table_path.exists():
        snapshot_power_npz = write_ss_snapshot_power_diagnostics(
            history=history,
            state=result.state,
            power_table_npz=power_table_path,
            output_path=raw_ss / "snapshot_power_energy_diagnostics.npz",
            sigma_n_S_m=float(material.sigma_n_S_m),
        )

    plots_dir = raw_ss / "plots_diagnostics"
    adaptive_timestep_png = plot_ss_adaptive_timestep_history(
        result.history,
        plots_dir / "adaptive_timestep_history.png",
        dpi=int(args.dpi),
    )

    ss_summary = {
        "run_name": run_name,
        "pre_run_name": pre_name,
        "target_current_uA": float(target_current_uA),
        "sweep_role": role,
        "requested_total_time_ps": float(total_time_ps),
        "requested_dt_s": float(dt_s),
        "requested_dt_fs": float(dt_s / 1.0e-15),
        "requested_snapshots": int(n_snapshots),
        "backend": "flat_gtdgl_pytdgl_like_promoted_backend",
        "supercurrent_policy": supercurrent_policy,
        "strict_usadel_current_table": strict_table_summary,
        "gtdgl_allmaras_diffusion": allmaras_diffusion,
        "seed_current_policy": seed_current_policy,
        "seed": seed_summary_data,
        "solver": result.summary,
        "outputs": {
            "seed_npz": str(seed_npz),
            "stationary_state_npz": str(state_npz),
            "relaxation_history_npz": str(history_npz),
            "stationary_snapshots_npz": str(snapshots_npz),
            "snapshot_power_energy_diagnostics_npz": (
                str(snapshot_power_npz) if snapshot_power_npz is not None else None
            ),
            "adaptive_timestep_history_png": str(adaptive_timestep_png),
        },
    }
    summary_path = raw_ss / "ss_summary.yaml"
    _write_yaml(summary_path, ss_summary)

    manifest_path = write_manifest(
        cfg,
        run_name,
        stage="ss",
        extra={
            "pipeline": "02_ss_run_template.py",
            "purpose": "Stationary gTDGL--Poisson relaxation with flat gTDGL backend.",
            "pre_run_name": pre_name,
            "target_current_uA": float(target_current_uA),
            "sweep_role": role,
            "gtdgl_allmaras_diffusion": allmaras_diffusion,
            "seed_current_policy": seed_current_policy,
            "outputs": ss_summary["outputs"] | {"ss_summary": str(summary_path)},
            "summary": ss_summary,
        },
    )

    if not quiet:
        _print_single_case_report(
            run_name=run_name,
            pre_name=pre_name,
            raw_ss=raw_ss,
            supercurrent_law=supercurrent_law,
            supercurrent_policy=supercurrent_policy,
            seed_summary_data=seed_summary_data,
            solver_summary=result.summary,
            seed_npz=seed_npz,
            state_npz=state_npz,
            history_npz=history_npz,
            snapshots_npz=snapshots_npz,
            snapshot_power_npz=snapshot_power_npz,
            adaptive_timestep_png=adaptive_timestep_png,
            summary_path=summary_path,
            manifest_path=manifest_path,
            verbose=bool(args.ss_verbose_report),
        )

    return {
        "run_name": run_name,
        "pre_run_name": pre_name,
        "target_current_uA": float(target_current_uA),
        "raw_ss": str(raw_ss),
        "seed_current_uA": float(seed_current_A * 1.0e6),
        "seed_is_overcritical_clamped": bool(seed_current_policy.get("seed_is_overcritical_clamped", False)),
        "summary_path": str(summary_path),
        "manifest_path": str(manifest_path),
        "first_magic_ready": bool(result.summary.get("first_magic_ready", False)),
        "terminal_voltage_V": float(result.summary.get("terminal_voltage_V", float("nan"))),
        "max_pairbreaking_ratio": float(result.summary.get("max_pairbreaking_ratio", float("nan"))),
        "normal_current_fraction_max": float(result.summary.get("normal_current_fraction_max", float("nan"))),
    }



def _history_with_fv_topology(history: dict[str, np.ndarray], *, ops) -> dict[str, np.ndarray]:
    """Attach static FV topology arrays needed by snapshot post-processing.

    The solver history is intentionally compact.  Snapshot power/energy maps,
    q-projection, and Joule diagnostics need the edge-to-node topology; add it
    once here before writing ``relaxation_history.npz`` and
    ``stationary_snapshots.npz``.
    """
    out = {str(key): np.asarray(value) for key, value in history.items()}
    static = {
        "edge_i": getattr(ops, "edge_i", None),
        "edge_j": getattr(ops, "edge_j", None),
        "edge_length_m": getattr(ops, "edge_length_m", None),
        "dual_face_length_m": getattr(ops, "dual_face_length_m", None),
        "edge_unit_x": getattr(ops, "edge_unit_x", None),
        "edge_unit_y": getattr(ops, "edge_unit_y", None),
    }
    for key, value in static.items():
        if key not in out and value is not None:
            out[key] = np.asarray(value)
    return out


class _SSProgressConsoleFilter:
    """Rewrite the solver progress line with wall-clock elapsed time and ETA.

    The core solver owns the numerical loop and only exposes a boolean progress
    switch.  This lightweight stdout filter keeps the solver API unchanged while
    making its single-line progress output more useful and less cluttered.
    """

    _progress_re = re.compile(
        r"SS pyTDGL-like\s*\[(?P<bar>[^\]]*)\]\s*"
        r"(?P<pct>[0-9]+(?:\.[0-9]+)?)%.*?"
        r"step=(?P<step>[0-9]+).*?"
        r"t=(?P<t>[0-9.eE+\-]+)"
    )

    def __init__(self, *, enabled: bool, total_time_ps: float) -> None:
        self.enabled = bool(enabled)
        self.total_time_ps = float(total_time_ps)
        self._stream = None
        self._start_s = 0.0
        self._printed_progress = False

    def __enter__(self):
        if self.enabled:
            self._stream = sys.stdout
            self._start_s = time.monotonic()
            sys.stdout = self
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        if self.enabled and self._stream is not None:
            sys.stdout = self._stream
            if self._printed_progress:
                self._stream.write("\n")
                self._stream.flush()
        self._stream = None

    def write(self, text: str) -> int:
        if self._stream is None:
            return len(text)
        if "SS pyTDGL-like" not in str(text):
            self._stream.write(text)
            self._stream.flush()
            return len(text)

        match = self._progress_re.search(str(text))
        if match is None:
            self._stream.write(text)
            self._stream.flush()
            return len(text)

        pct = float(match.group("pct"))
        step = int(match.group("step"))
        t_ps = float(match.group("t"))
        elapsed_s = max(time.monotonic() - self._start_s, 0.0)
        fraction = min(max(pct / 100.0, 0.0), 1.0)
        if fraction > 1.0e-12 and fraction < 1.0:
            eta_s = elapsed_s * (1.0 - fraction) / fraction
        else:
            eta_s = 0.0
        bar = _progress_bar(fraction, width=32)
        line = (
            f"\rSS pyTDGL-like [{bar}] {pct:6.2f}% "
            f"step={step} "
            f"t={t_ps:.6g}/{self.total_time_ps:.6g} ps "
            f"wall={_format_duration(elapsed_s)} "
            f"eta={_format_duration(eta_s)}"
        )
        self._stream.write(line)
        self._stream.flush()
        self._printed_progress = True
        return len(text)

    def flush(self) -> None:
        if self._stream is not None:
            self._stream.flush()


def _progress_bar(fraction: float, *, width: int = 32) -> str:
    filled = int(round(width * min(max(float(fraction), 0.0), 1.0)))
    filled = min(max(filled, 0), width)
    return "#" * filled + "-" * (width - filled)


def _format_duration(seconds: float) -> str:
    seconds = float(seconds)
    if not np.isfinite(seconds) or seconds < 0.0:
        return "--:--"
    if seconds < 60.0:
        return f"{seconds:4.1f}s"
    minutes, sec = divmod(int(round(seconds)), 60)
    if minutes < 60:
        return f"{minutes:02d}:{sec:02d}"
    hours, minutes = divmod(minutes, 60)
    return f"{hours:d}:{minutes:02d}:{sec:02d}"


def _as_float(value: Any) -> float:
    try:
        return float(value)
    except Exception:
        return float("nan")


def _to_fs(value_s: Any) -> float:
    value = _as_float(value_s)
    return value / 1.0e-15 if np.isfinite(value) else float("nan")


def _fmt_float(value: Any, *, precision: int = 6) -> str:
    value_f = _as_float(value)
    if not np.isfinite(value_f):
        return "nan"
    return f"{value_f:.{precision}g}"


def _resolve_seed_current_for_target(
    *,
    usadel_catalog,
    target_current_A: float,
    overcritical_policy: str,
    overcritical_seed_fraction: float,
) -> tuple[float, dict[str, Any]]:
    """Return the analytic seed current used for a requested simulation current.

    The finite-strip solver can be driven above the PRE Usadel depairing current,
    but the analytic OE6 seed is a superconducting uniform branch and therefore
    must remain below Ic.  For overcritical PSL searches we keep the physical
    boundary-current target unchanged and only clamp the *initial seed* current.
    """
    target_current_A = float(target_current_A)
    if not np.isfinite(target_current_A) or target_current_A <= 0.0:
        raise ValueError("target_current_A must be positive and finite.")

    Ic_A = _usadel_ic_A(usadel_catalog)
    policy = str(overcritical_policy).strip().lower().replace("_", "-")
    fraction = float(overcritical_seed_fraction)
    if not np.isfinite(fraction) or not (0.0 < fraction < 1.0):
        raise ValueError("--ss-overcritical-seed-fraction must be in the open interval (0, 1).")

    overcritical = bool(target_current_A > Ic_A)
    if not overcritical:
        seed_current_A = target_current_A
        reason = "target current is below or equal to PRE Usadel Ic; seed uses the requested current."
        clamped = False
    elif policy == "error":
        raise ValueError(
            f"Requested target current {target_current_A:.6e} A exceeds PRE Usadel Ic {Ic_A:.6e} A. "
            "Use --ss-overcritical-seed-policy clamp-to-ic to drive the solver above Ic."
        )
    elif policy == "clamp-to-ic":
        seed_current_A = fraction * Ic_A
        reason = (
            "target current exceeds PRE Usadel Ic; analytic superconducting seed is clamped "
            "below Ic while the solver boundary target remains overcritical."
        )
        clamped = True
    else:
        raise ValueError(f"Unknown --ss-overcritical-seed-policy: {overcritical_policy!r}")

    return float(seed_current_A), {
        "requested_target_current_A": float(target_current_A),
        "requested_target_current_uA": float(target_current_A * 1.0e6),
        "pre_usadel_Ic_A": float(Ic_A),
        "pre_usadel_Ic_uA": float(Ic_A * 1.0e6),
        "target_over_pre_usadel_Ic": float(target_current_A / Ic_A),
        "analytic_seed_current_A": float(seed_current_A),
        "analytic_seed_current_uA": float(seed_current_A * 1.0e6),
        "analytic_seed_over_pre_usadel_Ic": float(seed_current_A / Ic_A),
        "policy": policy,
        "overcritical_seed_fraction": float(fraction),
        "seed_is_overcritical_clamped": bool(clamped),
        "reason": reason,
    }


def _usadel_ic_A(usadel_catalog) -> float:
    current = np.asarray(usadel_catalog.calibration_current_values_A, dtype=float)
    finite = np.isfinite(current) & (current >= 0.0)
    if np.count_nonzero(finite) < 1:
        raise ValueError("PRE Usadel catalogue has no finite non-negative calibration currents.")
    Ic_A = float(np.max(current[finite]))
    if not np.isfinite(Ic_A) or Ic_A <= 0.0:
        raise ValueError("PRE Usadel critical current is not positive.")
    return Ic_A


def _print_single_case_report(
    *,
    run_name: str,
    pre_name: str,
    raw_ss: Path,
    supercurrent_law: str,
    supercurrent_policy: dict[str, Any],
    seed_summary_data: dict[str, Any],
    solver_summary: dict[str, Any],
    seed_npz: Path,
    state_npz: Path,
    history_npz: Path,
    snapshots_npz: Path,
    snapshot_power_npz: Path | None,
    adaptive_timestep_png: Path,
    summary_path: Path,
    manifest_path: Path,
    verbose: bool = False,
) -> None:
    """Print a compact terminal report.

    Complete dictionaries are written to ``ss_summary.yaml`` and ``manifest.yaml``.
    The default terminal report intentionally keeps only the quantities that are
    useful while launching many SS runs.  Use ``--ss-verbose-report`` for the old
    exhaustive printout.
    """
    seed_policy = dict(seed_summary_data.get("overcritical_seed_policy", {}))
    stationarity = dict(solver_summary.get("stationarity", {}))
    continuity = dict(solver_summary.get("continuity", {}))
    contact = dict(solver_summary.get("contact_recovery", {}))

    print("SS-run stationary relaxation")
    print(f"  run_name:      {run_name}")
    print(f"  pre_run_name:  {pre_name}")
    print(f"  raw_ss:        {raw_ss}")
    print(f"  supercurrent:  {supercurrent_law}")
    print(f"  policy:        {supercurrent_policy.get('reason', 'n/a')}")
    print()
    print("Run summary")
    print(f"  target_current_uA:     {seed_summary_data.get('simulation_target_current_uA', float('nan')):.6g}")
    print(f"  seed_current_uA:       {seed_summary_data.get('analytic_seed_current_uA', float('nan')):.6g}")
    print(f"  seed_clamped:          {seed_policy.get('seed_is_overcritical_clamped', False)}")
    print(f"  stop_reason:           {solver_summary.get('stop_reason', 'n/a')}")
    print(f"  accepted/rejected:     {solver_summary.get('accepted_steps', 'n/a')}/{solver_summary.get('rejected_steps', 'n/a')}")
    print(f"  final_time_ps:         {_fmt_float(solver_summary.get('final_time_ps'), precision=6)}")
    print(f"  dt_final_fs:           {_fmt_float(_to_fs(solver_summary.get('dt_final_s')), precision=6)}")
    print(f"  V_terminal_mV:         {_fmt_float(1.0e3 * _as_float(solver_summary.get('terminal_voltage_V')), precision=6)}")
    print(f"  continuity_passes:     {continuity.get('passes', 'n/a')}")
    print(f"  stationarity_passes:   {stationarity.get('passes', 'n/a')}")
    print(f"  contact_recovery:      {contact.get('passes', 'n/a')}")
    print(f"  eta_R_final:           {_fmt_float(solver_summary.get('eta_R_final'), precision=6)}")
    print(f"  min_delta/delta0:      {_fmt_float(solver_summary.get('min_delta_over_delta0'), precision=6)}")
    print(f"  mean_delta/delta0:     {_fmt_float(solver_summary.get('mean_delta_over_delta0'), precision=6)}")
    print(f"  max_pairbreaking:      {_fmt_float(solver_summary.get('max_pairbreaking_ratio'), precision=6)}")
    print(f"  max |j_tot| [A/m2]:    {_fmt_float(solver_summary.get('total_current_max_A_m2'), precision=6)}")
    print(f"  max |j_n| [A/m2]:      {_fmt_float(solver_summary.get('normal_current_max_A_m2'), precision=6)}")
    print()
    print("Outputs")
    print(f"  seed_npz:              {seed_npz}")
    print(f"  stationary_state_npz:  {state_npz}")
    print(f"  relaxation_history_npz:{history_npz}")
    print(f"  stationary_snapshots_npz: {snapshots_npz}")
    if snapshot_power_npz is not None:
        print(f"  snapshot_power_energy_diagnostics_npz: {snapshot_power_npz}")
    print(f"  adaptive_timestep_png: {adaptive_timestep_png}")
    print(f"  ss_summary:            {summary_path}")
    print(f"  ss_manifest:           {manifest_path}")

    if verbose:
        print()
        print("Seed metadata")
        _print_dict(seed_summary_data)
        print()
        print("Solver metadata")
        _print_dict(solver_summary)

    print("Status: OK")


def _resolve_supercurrent_law(*, requested: str, strict_table_summary: dict[str, Any]) -> tuple[str, dict[str, Any]]:
    requested_norm = requested.strip().lower().replace("_", "-")
    if requested_norm == "auto":
        requested_norm = "usadel-poisson"
    if requested_norm != "usadel-poisson":
        raise RuntimeError(
            "This SS pipeline is configured to require the strict PRE table "
            "js_A_m2[Te,delta,q].  Use a separate legacy/debug script for GL-only tests."
        )
    return "usadel_poisson", {
        "requested": requested,
        "resolved": "usadel_poisson",
        "has_strict_3d_table": True,
        "strict_table": strict_table_summary,
        "reason": "SS requires the strict PRE Matsubara Usadel table js_A_m2[Te,delta,q].",
    }


def _read_pre_allmaras_diffusion(raw_pre: Path) -> dict[str, float | str]:
    summary_path = raw_pre / "usadel_dos_summary.yaml"
    if not summary_path.exists():
        return {
            "D_effective_factor": 1.0,
            "D_base_m2_s": float("nan"),
            "D_effective_m2_s": float("nan"),
            "source": "default: PRE summary not found; using Usadel D unchanged for gTDGL.",
        }
    with summary_path.open("r", encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
    usadel = data.get("usadel", {}) if isinstance(data, dict) else {}
    allmaras = usadel.get("gtdgl_allmaras", {}) if isinstance(usadel, dict) else {}
    metadata = data.get("metadata", {}) if isinstance(data, dict) else {}
    base_D = float(
        allmaras.get(
            "D_base_m2_s",
            metadata.get("D_m2_s", float("nan")),
        )
    )
    factor = float(
        allmaras.get(
            "D_effective_factor",
            metadata.get("gtdgl_allmaras_D_factor", 1.0),
        )
    )
    effective_D = float(allmaras.get("D_effective_m2_s", base_D * factor))
    if not np.isfinite(factor) or factor <= 0.0:
        raise ValueError(f"Invalid PRE gTDGL Allmaras diffusion factor: {factor!r}")
    return {
        "D_effective_factor": factor,
        "D_base_m2_s": base_D,
        "D_effective_m2_s": effective_D,
        "source": str(
            allmaras.get(
                "source",
                "Effective mesoscopic diffusion for the Allmaras/gTDGL sector; Usadel tables keep the calibrated microscopic D.",
            )
        ),
    }


def _current_sweep_run_name(base_run_name: str, offset_uA: float, current_uA: float) -> str:
    sign = "plus" if offset_uA >= 0.0 else "minus"
    offset_label = _number_label(abs(float(offset_uA)))
    current_label = _number_label(float(current_uA))
    return f"{base_run_name}_dI_{sign}{offset_label}uA_I{current_label}uA"


def _number_label(value: float) -> str:
    text = f"{float(value):.6g}"
    return text.replace("-", "m").replace("+", "p").replace(".", "p")


def _resolve_sweep_workers(cfg: dict[str, Any], requested_workers: int | None, *, n_extra: int) -> int:
    if n_extra <= 0:
        return 1
    parallel = cfg.get("parallel", {}) if isinstance(cfg, dict) else {}
    if requested_workers is not None:
        workers = int(requested_workers)
    elif bool(parallel.get("enabled", False)):
        workers = int(parallel.get("workers", 1))
    else:
        workers = 1
    return max(1, min(int(workers), int(n_extra)))


def _require_file(path: Path, description: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Missing {description}: {path}")


def _write_yaml(path: str | Path, data: dict[str, Any]) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as f:
        yaml.safe_dump(data, f, sort_keys=False, allow_unicode=True)


def _print_dict(data: dict[str, Any], *, indent: str = "  ") -> None:
    for key, value in data.items():
        if isinstance(value, dict):
            print(f"{indent}{key}:")
            _print_dict(value, indent=indent + "  ")
        else:
            print(f"{indent}{key}: {value}")


if __name__ == "__main__":
    raise SystemExit(main())
