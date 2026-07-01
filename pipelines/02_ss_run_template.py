"""Stationary gTDGL--Poisson SS-run using the flat gTDGL backend.

This pipeline loads PRE-run outputs, builds the analytic stationary seed, and
relaxes it with ``pysnspd.gtdgl.solve_stationary_pytdgl_like`` promoted to the
package root.  The preferred supercurrent law is the Matsubara Usadel table
stored by ``01_prerun_template.py``.  If an older PRE-run is used, ``auto`` falls
back to the native GL current and records that choice in the manifest.
"""
from __future__ import annotations

import argparse
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
from pysnspd.gtdgl.usadel_current import attach_usadel_supercurrent_table_from_npz


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
    parser.add_argument("--ss-steps", type=int, default=2000)
    parser.add_argument("--ss-dt-fs", type=float, default=0.25)
    parser.add_argument("--ss-tau-scale", type=float, default=0.10)
    parser.add_argument("--ss-target-current-uA", type=float, default=None)
    parser.add_argument("--ss-snapshots", type=int, default=6)
    parser.add_argument("--ss-progress", action="store_true")
    parser.add_argument("--ss-terminal-psi", type=float, default=0.0)
    parser.add_argument(
        "--ss-terminal-healing-xi",
        type=float,
        default=None,
        help=(
            "Apply a smooth tanh contact-healing envelope to the initial |Delta| seed. "
            "Use 2.5 for the first metallic-contact SS objective."
        ),
    )
    parser.add_argument("--ss-terminal-healing-fraction", type=float, default=0.95)
    parser.add_argument("--ss-stationarity-eta", type=float, default=1.0e-5)
    parser.add_argument("--ss-stationarity-delta-rel", type=float, default=1.0e-4)
    parser.add_argument("--ss-stationarity-phi-rel", type=float, default=1.0e-4)
    parser.add_argument("--ss-convergence-min-steps", type=int, default=50)
    parser.add_argument("--ss-continuity-rms-tol", type=float, default=1.0e-6)
    parser.add_argument("--ss-continuity-max-tol", type=float, default=1.0e-3)
    parser.add_argument("--ss-continuity-poisson-tol", type=float, default=1.0e-9)
    parser.add_argument("--ss-recovery-min-xi", type=float, default=1.5)
    parser.add_argument("--ss-recovery-max-xi", type=float, default=4.0)
    parser.add_argument(
        "--ss-no-adaptive",
        action="store_true",
        help="Disable adaptive time stepping in the pyTDGL-like core.",
    )
    parser.add_argument(
        "--ss-supercurrent-law",
        choices=("auto", "gl", "usadel-poisson"),
        default="auto",
        help=(
            "auto uses the PRE Matsubara Usadel table when present, otherwise GL. "
            "Use usadel-poisson to require the table and fail loudly if absent."
        ),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    cfg = validate_config(load_config(args.config))

    ss_layout = create_run_layout(cfg, args.run_name)
    run_name = ss_layout["run_name"]
    pre_name = args.pre_run_name or run_name
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

    base_usadel_catalog = load_usadel_catalog_npz(usadel_path)
    usadel_catalog = attach_usadel_supercurrent_table_from_npz(base_usadel_catalog, usadel_path)
    supercurrent_law, supercurrent_policy = _resolve_supercurrent_law(
        requested=args.ss_supercurrent_law,
        usadel_npz=usadel_path,
    )

    target_current_A = None
    if args.ss_target_current_uA is not None:
        target_current_A = float(args.ss_target_current_uA) * 1.0e-6

    seed = build_stationary_seed(
        mesh=mesh,
        edge_data=edge_data,
        usadel_catalog=base_usadel_catalog,
        I_bias_A=target_current_A,
        phase_origin=args.phase_origin,
    )
    seed_npz = save_stationary_seed_npz(seed, raw_ss / "stationary_seed.npz")
    seed_summary_data = seed_summary(seed)

    material = build_gtdgl_material(
        cfg,
        base_usadel_catalog,
        tau_scale=float(args.ss_tau_scale),
    )
    if target_current_A is None:
        target_current_A = float(seed.metadata["target_current_A"])

    result = solve_stationary_pytdgl_like(
        mesh=mesh,
        edge_data=edge_data,
        seed=seed,
        material=material,
        ops=ops,
        steps=int(args.ss_steps),
        dt_s=float(args.ss_dt_fs) * 1.0e-15,
        target_current_A=target_current_A,
        usadel_catalog=usadel_catalog,
        terminal_psi=float(args.ss_terminal_psi),
        adaptive=not bool(args.ss_no_adaptive),
        n_snapshots=int(args.ss_snapshots),
        progress=bool(args.ss_progress),
        supercurrent_law=supercurrent_law,
        terminal_healing_xi=args.ss_terminal_healing_xi,
        terminal_healing_fraction=float(args.ss_terminal_healing_fraction),
        stationarity_eta=float(args.ss_stationarity_eta),
        stationarity_delta_rel=float(args.ss_stationarity_delta_rel),
        stationarity_phi_rel=float(args.ss_stationarity_phi_rel),
        convergence_min_steps=int(args.ss_convergence_min_steps),
        continuity_rms_tol=float(args.ss_continuity_rms_tol),
        continuity_max_tol=float(args.ss_continuity_max_tol),
        continuity_poisson_tol=float(args.ss_continuity_poisson_tol),
        recovery_min_xi=float(args.ss_recovery_min_xi),
        recovery_max_xi=float(args.ss_recovery_max_xi),
    )

    state_npz = save_stationary_state_npz(result.state, raw_ss / "stationary_state.npz")
    history_npz = save_relaxation_history_npz(result.history, raw_ss / "relaxation_history.npz")

    ss_summary = {
        "run_name": run_name,
        "pre_run_name": pre_name,
        "backend": "flat_gtdgl_pytdgl_like_promoted_backend",
        "supercurrent_policy": supercurrent_policy,
        "seed": seed_summary_data,
        "solver": result.summary,
        "outputs": {
            "seed_npz": str(seed_npz),
            "stationary_state_npz": str(state_npz),
            "relaxation_history_npz": str(history_npz),
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
            "outputs": ss_summary["outputs"] | {"ss_summary": str(summary_path)},
            "summary": ss_summary,
        },
    )

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
    _print_dict(result.summary)
    print()
    print("Outputs")
    print(f"  seed_npz:              {seed_npz}")
    print(f"  stationary_state_npz:  {state_npz}")
    print(f"  relaxation_history_npz:{history_npz}")
    print(f"  ss_summary:            {summary_path}")
    print(f"  ss_manifest:           {manifest_path}")
    print("Status: OK")
    return 0


def _resolve_supercurrent_law(*, requested: str, usadel_npz: Path) -> tuple[str, dict[str, Any]]:
    requested_norm = requested.strip().lower().replace("_", "-")
    has_table = _npz_has_usadel_current_table(usadel_npz)
    if requested_norm == "gl":
        return "gl", {
            "requested": requested,
            "resolved": "gl",
            "has_matsubara_table": has_table,
            "reason": "User explicitly requested the native GL current law.",
        }
    if requested_norm == "usadel-poisson":
        if not has_table:
            raise RuntimeError(
                "--ss-supercurrent-law usadel-poisson requires a PRE Usadel current table. "
                "Run the updated 01_prerun_template.py first, or use --ss-supercurrent-law gl."
            )
        return "usadel_poisson", {
            "requested": requested,
            "resolved": "usadel_poisson",
            "has_matsubara_table": True,
            "reason": "User explicitly required the PRE Matsubara Usadel supercurrent table.",
        }
    if requested_norm != "auto":
        raise ValueError(f"Unknown supercurrent law: {requested!r}")
    if has_table:
        return "usadel_poisson", {
            "requested": requested,
            "resolved": "usadel_poisson",
            "has_matsubara_table": True,
            "reason": "AUTO selected the PRE Matsubara Usadel supercurrent table.",
        }
    return "gl", {
        "requested": requested,
        "resolved": "gl",
        "has_matsubara_table": False,
        "reason": "AUTO fell back to GL because the PRE NPZ has no Usadel current table.",
    }


def _npz_has_usadel_current_table(path: Path) -> bool:
    current_keys = {
        "js_A_m2",
        "j_s_A_m2",
        "supercurrent_density_A_m2",
        "current_density_A_m2",
        "j_s_grid_A_m2",
        "js_grid_A_m2",
    }
    q_keys = {"q_axis_m_inv", "q_values_m_inv", "q_grid_m_inv"}
    with np.load(path, allow_pickle=False) as data:
        keys = set(data.files)
    return bool(keys & current_keys) and bool(keys & q_keys)


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
