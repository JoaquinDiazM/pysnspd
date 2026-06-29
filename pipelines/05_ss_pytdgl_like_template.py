"""pyTDGL-like stationary SS backend for OE7 comparisons.

This pipeline intentionally leaves the existing OE7 solver untouched.  It loads
an OE6/PRE seed exactly like ``02_ss_run_template.py`` and evolves it with the
new ``pysnspd.gtdgl.pytdgl_like`` backend, whose solver class/method names mirror
pyTDGL's CPU no-screening path.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import yaml

from pysnspd.config import load_config, validate_config
from pysnspd.io.manager import create_run_layout, write_manifest
from pysnspd.mesh.delaunay import load_mesh_npz
from pysnspd.mesh.edges import load_edges_npz
from pysnspd.usadel.catalog import load_usadel_catalog_npz
from pysnspd.gtdgl.seed import build_stationary_seed, save_stationary_seed_npz, seed_summary
from pysnspd.gtdgl.material import build_gtdgl_material
from pysnspd.gtdgl.operators import build_fv_operators
from pysnspd.gtdgl.state_io import save_relaxation_history_npz, save_stationary_state_npz
from pysnspd.gtdgl.pytdgl_like import solve_stationary_pytdgl_like

try:
    from pysnspd.plotting.ss_run import plot_ss_available_snapshots, plot_ss_relaxation_history
except Exception:  # pragma: no cover
    plot_ss_available_snapshots = None
    plot_ss_relaxation_history = None

try:
    from pysnspd.plotting.pytdgl_like import (
        plot_pytdgl_like_native_history,
        plot_pytdgl_like_native_edge_currents,
        plot_pytdgl_like_poisson_snapshots,
    )
except Exception:  # pragma: no cover
    plot_pytdgl_like_native_history = None
    plot_pytdgl_like_native_edge_currents = None
    plot_pytdgl_like_poisson_snapshots = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the pyTDGL-like OE7 stationary comparison backend."
    )
    parser.add_argument("--config", required=True, help="Path to YAML config.")
    parser.add_argument("--run-name", required=True, help="Run name where outputs are written.")
    parser.add_argument("--pre-run-name", default=None, help="Run name containing PRE outputs.")
    parser.add_argument("--I-bias-A", type=float, default=None, help="Override bias current.")
    parser.add_argument("--T-bias-K", type=float, default=None, help="Override bias temperature.")
    parser.add_argument(
        "--phase-origin",
        choices=["center", "left"],
        default="center",
        help="Origin used in theta=q*(x-x0) for the analytic seed.",
    )
    parser.add_argument("--ss-steps", type=int, default=2000, help="Number of pyTDGL-like steps.")
    parser.add_argument("--ss-dt-fs", type=float, default=0.01, help="Initial time step in fs.")
    parser.add_argument("--ss-tau-scale", type=float, default=1.0, help="Material tau scale.")
    parser.add_argument(
        "--ss-terminal-psi",
        default="0",
        help="pyTDGL terminal_psi value: 0, none, or a Python complex literal.",
    )
    parser.add_argument("--ss-no-adapt-dt", action="store_true", help="Disable pyTDGL adaptive dt.")
    parser.add_argument("--ss-n-snapshots", type=int, default=6, help="Number of final-state snapshots to emit.")
    parser.add_argument("--dpi", type=int, default=480, help="DPI for optional diagnostic plots.")
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
    for path in (mesh_path, edges_path, usadel_path):
        _require_existing(path)

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
    material = build_gtdgl_material(cfg, usadel_catalog, tau_scale=args.ss_tau_scale)
    ops = build_fv_operators(mesh, edge_data)

    terminal_psi = _parse_terminal_psi(args.ss_terminal_psi)
    result = solve_stationary_pytdgl_like(
        mesh=mesh,
        edge_data=edge_data,
        seed=seed,
        material=material,
        ops=ops,
        steps=args.ss_steps,
        dt_s=args.ss_dt_fs * 1.0e-15,
        target_current_A=seed_sum["I_bias_A"],
        terminal_psi=terminal_psi,
        adaptive=not args.ss_no_adapt_dt,
        n_snapshots=args.ss_n_snapshots,
    )

    state_npz = save_stationary_state_npz(result.state, raw_ss / "ss_state_pytdgl_like.npz")
    history_npz = save_relaxation_history_npz(result.history, raw_ss / "ss_pytdgl_like_history.npz")
    summary_path = raw_ss / "ss_pytdgl_like_summary.yaml"
    summary_doc = {
        "run_name": run_name,
        "pre_run_name": pre_run_name,
        "stage": "ss_pytdgl_like",
        "seed": seed_sum,
        "stationary": result.summary,
        "metadata": {
            "backend": "pytdgl_like_minimal_no_screening",
            "reference_commit": "fc18de6",
            "pytdgl_reference": "loganbvh/py-tdgl CPU no-screening solver structure",
            "terminal_psi": args.ss_terminal_psi,
            "inputs": {
                "mesh_npz": str(mesh_path),
                "edges_npz": str(edges_path),
                "usadel_npz": str(usadel_path),
            },
        },
    }
    with summary_path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(summary_doc, f, sort_keys=False, allow_unicode=True)

    snapshot_plots = {}
    history_plot = None
    if plot_ss_available_snapshots is not None:
        snapshot_plots = plot_ss_available_snapshots(mesh, result.history, plots_diag, dpi=args.dpi, ncols=3)
    if plot_ss_relaxation_history is not None:
        history_plot = plot_ss_relaxation_history(
            result.history,
            plots_diag / "ss_pytdgl_like_relaxation_history.png",
            dpi=args.dpi,
        )

    native_plots = {}
    if plot_pytdgl_like_native_history is not None:
        native_plots["native_history"] = plot_pytdgl_like_native_history(
            result.history,
            plots_diag / "pytdgl_like_native_history.png",
            dpi=args.dpi,
        )
    if plot_pytdgl_like_poisson_snapshots is not None:
        native_plots["poisson_snapshots"] = plot_pytdgl_like_poisson_snapshots(
            mesh,
            result.history,
            plots_diag / "pytdgl_like_poisson_terms_snapshots.png",
            dpi=args.dpi,
            ncols=3,
        )
    if plot_pytdgl_like_native_edge_currents is not None:
        native_plots["native_edge_currents"] = plot_pytdgl_like_native_edge_currents(
            mesh,
            result.history,
            plots_diag / "pytdgl_like_native_edge_currents.png",
            dpi=args.dpi,
        )

    manifest = write_manifest(
        cfg,
        run_name,
        stage="ss",
        extra={
            "pipeline": "05_ss_pytdgl_like_template.py",
            "purpose": "pyTDGL-like stationary solver comparison backend.",
            "pre_run_name": pre_run_name,
            "outputs": {
                "seed_npz": str(seed_npz),
                "state_npz": str(state_npz),
                "history_npz": str(history_npz),
                "summary_yaml": str(summary_path),
                "history_plot": None if history_plot is None else str(history_plot),
                "snapshot_plots": {k: str(v) for k, v in snapshot_plots.items()},
                "native_plots": {k: str(v) for k, v in native_plots.items()},
            },
            "seed_summary": seed_sum,
            "stationary_summary": result.summary,
        },
    )

    print("SS pyTDGL-like backend")
    print(f"run_name          : {run_name}")
    print(f"pre_run_name      : {pre_run_name}")
    print(f"raw_ss            : {raw_ss}")
    print(f"plots_diagnostics : {plots_diag}")
    print(f"state_npz         : {state_npz}")
    print(f"history_npz       : {history_npz}")
    print(f"summary_yaml      : {summary_path}")
    print(f"manifest          : {manifest}")
    print("Stationary summary")
    for key in (
        "accepted_steps",
        "final_time_ps",
        "terminal_voltage_V",
        "current_residual",
        "min_delta_over_delta0",
        "mean_delta_over_delta0",
        "max_pairbreaking_ratio",
        "normal_current_fraction_max",
        "native_poisson_residual_rel_final",
        "native_boundary_rhs_norm_final",
        "pytdgl_u",
        "pytdgl_gamma",
    ):
        print(f"  {key}: {result.summary.get(key)}")
    print("Status: OK")
    return 0


def _parse_terminal_psi(text: str):
    value = str(text).strip().lower()
    if value in {"none", "null"}:
        return None
    if value in {"0", "0.0"}:
        return 0.0
    return complex(text)


def _require_existing(path: Path) -> None:
    if not path.exists() or not path.is_file():
        raise FileNotFoundError(f"Required input does not exist: {path}")


if __name__ == "__main__":
    raise SystemExit(main())
