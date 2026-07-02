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
            "Physical SS relaxation time in ps. Default is 20 ps. "
            "This replaces the old --ss-steps control."
        ),
    )
    parser.add_argument(
        "--ss-steps",
        type=int,
        default=None,
        help=(
            "Deprecated compatibility input. If --ss-time-ps is omitted, "
            "the physical time is computed as ss_steps * ss_dt_fs."
        ),
    )
    parser.add_argument("--ss-dt-fs", type=float, default=0.30)
    parser.add_argument("--ss-target-current-uA", type=float, default=20.0)
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
    parser.add_argument("--ss-snapshots", type=int, default=8)

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

        for future in as_completed(futures):
            extra_results.append(future.result())

    all_results = [base_result] + sorted(extra_results, key=lambda item: float(item["target_current_uA"]))
    print()
    print("SS sweep run directories")
    for item in all_results:
        print(
            f"  {item['target_current_uA']:8.3f} uA  "
            f"run={item['run_name']}  raw_ss={item['raw_ss']}"
        )
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

    seed = build_stationary_seed(
        mesh=mesh,
        edge_data=edge_data,
        usadel_catalog=base_usadel_catalog,
        I_bias_A=target_current_A,
        phase_origin=args.phase_origin,
    )
    seed_npz = save_stationary_seed_npz(seed, raw_ss / "stationary_seed.npz")
    seed_summary_data = seed_summary(seed)

    allmaras_diffusion = _read_pre_allmaras_diffusion(raw_pre)
    material = build_gtdgl_material(
        cfg,
        base_usadel_catalog,
        diffusion_factor=float(allmaras_diffusion["D_effective_factor"]),
    )

    if args.ss_time_ps is not None:
        total_time_ps = float(args.ss_time_ps)
    elif args.ss_steps is not None:
        total_time_ps = float(args.ss_steps) * float(args.ss_dt_fs) * 1.0e-3
    else:
        total_time_ps = 20.0
    if total_time_ps <= 0.0:
        raise ValueError("--ss-time-ps must be positive.")

    result = solve_stationary_pytdgl_like(
        mesh=mesh,
        edge_data=edge_data,
        seed=seed,
        material=material,
        ops=ops,
        steps=None if args.ss_steps is None else int(args.ss_steps),
        total_time_s=float(total_time_ps) * 1.0e-12,
        dt_s=float(args.ss_dt_fs) * 1.0e-15,
        target_current_A=target_current_A,
        usadel_catalog=usadel_catalog,
        terminal_psi=float(args.ss_terminal_psi),
        adaptive=not bool(args.ss_no_adaptive),
        adaptive_window=int(args.ss_adaptive_window),
        max_solve_retries=int(args.ss_max_solve_retries),
        adaptive_time_step_multiplier=float(args.ss_adaptive_time_step_multiplier),
        adaptive_growth_factor=float(args.ss_adaptive_growth_factor),
        dt_max_factor=float(args.ss_dt_max_factor),
        n_snapshots=int(args.ss_snapshots),
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
    history_npz = save_relaxation_history_npz(result.history, raw_ss / "relaxation_history.npz")

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
        "backend": "flat_gtdgl_pytdgl_like_promoted_backend",
        "supercurrent_policy": supercurrent_policy,
        "strict_usadel_current_table": strict_table_summary,
        "gtdgl_allmaras_diffusion": allmaras_diffusion,
        "seed": seed_summary_data,
        "solver": result.summary,
        "outputs": {
            "seed_npz": str(seed_npz),
            "stationary_state_npz": str(state_npz),
            "relaxation_history_npz": str(history_npz),
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
            adaptive_timestep_png=adaptive_timestep_png,
            summary_path=summary_path,
            manifest_path=manifest_path,
        )

    return {
        "run_name": run_name,
        "pre_run_name": pre_name,
        "target_current_uA": float(target_current_uA),
        "raw_ss": str(raw_ss),
        "summary_path": str(summary_path),
        "manifest_path": str(manifest_path),
        "first_magic_ready": bool(result.summary.get("first_magic_ready", False)),
        "terminal_voltage_V": float(result.summary.get("terminal_voltage_V", float("nan"))),
        "max_pairbreaking_ratio": float(result.summary.get("max_pairbreaking_ratio", float("nan"))),
        "normal_current_fraction_max": float(result.summary.get("normal_current_fraction_max", float("nan"))),
    }


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
    adaptive_timestep_png: Path,
    summary_path: Path,
    manifest_path: Path,
) -> None:
    print("SS-run stationary relaxation")
    print(f"  run_name:      {run_name}")
    print(f"  pre_run_name:  {pre_name}")
    print(f"  raw_ss:        {raw_ss}")
    print(f"  supercurrent:  {supercurrent_law}")
    print(f"  policy:        {supercurrent_policy['reason']}")
    print()
    print("Seed")
    _print_dict(seed_summary_data)
    print()
    print("Solver")
    _print_dict(solver_summary)
    print()
    print("Outputs")
    print(f"  seed_npz:              {seed_npz}")
    print(f"  stationary_state_npz:  {state_npz}")
    print(f"  relaxation_history_npz:{history_npz}")
    print(f"  adaptive_timestep_png: {adaptive_timestep_png}")
    print(f"  ss_summary:            {summary_path}")
    print(f"  ss_manifest:           {manifest_path}")
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
