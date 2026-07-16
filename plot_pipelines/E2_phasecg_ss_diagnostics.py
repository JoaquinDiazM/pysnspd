"""Focused diagnostics for a completed normalized phase-CG SS run.

The pipeline only reads existing SS raw data.  It produces a snapshot atlas,
physical scalar histories, and numerical diagnostics for the corrected
Allmaras phase drive without modifying the simulation state.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import yaml

from pysnspd.analysis.ss_phasecg import build_phasecg_diagnostic_dataset
from pysnspd.analysis.ss_run import load_ss_run
from pysnspd.config import load_config, validate_config
from pysnspd.plotting.ss_phasecg_figures import make_phasecg_ss_figures


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot focused phase-CG diagnostics for a completed stationary SS run."
    )
    parser.add_argument("--config", required=True, help="Absolute path to the YAML project config.")
    parser.add_argument("--run-name", required=True, help="Existing stationary SS run name.")
    parser.add_argument(
        "--pre-run-name",
        default=None,
        help="PRE run containing the mesh and edge data; otherwise use ss_summary.yaml.",
    )
    parser.add_argument(
        "--center-width-nm",
        type=float,
        default=100.0,
        help="Longitudinal width of the central voltage/current probe in nm.",
    )
    parser.add_argument(
        "--wall-time-seconds",
        type=float,
        default=None,
        help=(
            "Measured total wall time. Per-step time is only an explicitly labelled "
            "attempt-count estimate because this run format stores no step timings."
        ),
    )
    parser.add_argument("--dpi", type=int, default=240)
    parser.add_argument(
        "--figures-subdir",
        default="E2_phasecg_diagnostics",
        help="Subdirectory below plots/<run>/figures.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    cfg = validate_config(load_config(args.config))
    run = load_ss_run(
        config_path=args.config,
        run_name=args.run_name,
        pre_run_name=args.pre_run_name,
    )
    dataset = build_phasecg_diagnostic_dataset(
        run,
        thickness_m=float(cfg["material"]["thickness_m"]),
        center_width_m=float(args.center_width_nm) * 1.0e-9,
        measured_wall_time_s=args.wall_time_seconds,
    )
    _attach_snapshot_power_diagnostics(dataset, run.raw_ss)

    figures_dir = run.figures_dir / str(args.figures_subdir)
    saved = make_phasecg_ss_figures(
        dataset=dataset,
        output_dir=figures_dir,
        dpi=int(args.dpi),
    )
    manifest_path = _write_manifest(
        run=run,
        figures_dir=figures_dir,
        saved=saved,
        dataset=dataset,
    )

    print("E2 phase-CG stationary-run diagnostics")
    print(f" run_name:       {run.run_name}")
    print(f" pre_run_name:   {run.pre_run_name}")
    print(f" raw_ss:         {run.raw_ss}")
    print(f" figures_dir:    {figures_dir}")
    print(f" center_probe:   {float(args.center_width_nm):.6g} nm")
    if args.wall_time_seconds is not None:
        print(
            " wall_time:      "
            f"{float(args.wall_time_seconds):.6g} s measured total; "
            "per-step curve is an attempt-count estimate"
        )
    print()
    print("Figures")
    for key, path in saved.items():
        print(f" {key}: {path}")
    print(f" manifest: {manifest_path}")
    print("Status: OK")
    return 0


def _write_manifest(
    *,
    run: Any,
    figures_dir: Path,
    saved: Mapping[str, Path],
    dataset: Mapping[str, Any],
) -> Path:
    phase_converged = np.asarray(
        dataset.get("allmaras_phase_convergence_converged", []),
        dtype=bool,
    )
    phase_residual = np.asarray(
        dataset.get("allmaras_phase_convergence_residual_rel", []),
        dtype=float,
    )
    phase_iterations = np.asarray(
        dataset.get("allmaras_phase_convergence_iterations", []),
        dtype=float,
    )
    div_max = np.asarray(dataset.get("div_j_normalized_max_snapshot", []), dtype=float)
    target_current = float(dataset.get("target_current_uA", np.nan))
    total_current = np.asarray(dataset.get("current_total_snapshot_uA", []), dtype=float)
    terminal_voltage = np.asarray(dataset.get("voltage_terminal_snapshot_mV", []), dtype=float)
    center_voltage = np.asarray(dataset.get("voltage_center_snapshot_mV", []), dtype=float)

    wall_total = dataset.get("measured_wall_time_s")
    wall_note = (
        "No measured total wall time was supplied; solve-attempt counts are plotted directly."
        if wall_total is None
        else (
            "The total wall time is measured externally. Per-step time is estimated in "
            "proportion to 1 + rejected attempts and integrates exactly to that measured total."
        )
    )
    manifest = {
        "schema_version": 2,
        "pipeline": "plot_pipelines/E2_phasecg_ss_diagnostics.py",
        "purpose": (
            "Focused validation of the normalized Allmaras phase-drive continuation, "
            "physical SS response, and adaptive-solver behavior."
        ),
        "run_name": run.run_name,
        "pre_run_name": run.pre_run_name,
        "raw_ss": str(run.raw_ss),
        "figures_dir": str(figures_dir),
        "figures": {key: str(path) for key, path in saved.items()},
        "normalizations": {
            "current_density": "j / j_avg",
            "current_divergence": "xi * div(j_total) / j_avg",
            "phase_gradient": "abs(q) * xi",
            "order_parameter": "abs(Delta) / Delta_BCS(0)",
        },
        "central_probe_width_nm": float(dataset.get("center_width_nm", np.nan)),
        "wall_time": {
            "measured_total_s": None if wall_total is None else float(wall_total),
            "per_step_is_estimated": wall_total is not None,
            "note": wall_note,
        },
        "summary": {
            "target_current_uA": target_current,
            "final_snapshot_total_current_uA": _last_finite(total_current),
            "final_snapshot_terminal_voltage_mV": _last_finite(terminal_voltage),
            "final_snapshot_center_voltage_mV": _last_finite(center_voltage),
            "max_normalized_bulk_current_divergence": _max_finite(div_max),
            "continuity_passes": bool(dataset.get("continuity_passes", False)),
            "stationarity_passes": bool(dataset.get("stationarity_passes", False)),
            "dynamic_stationarity_passes": bool(
                dataset.get("dynamic_stationarity_passes", False)
            ),
            "thermal_stationarity_passes": bool(
                dataset.get("thermal_stationarity_passes", False)
            ),
            "thermal_enabled": bool(dataset.get("thermal_enabled", False)),
            "phase_cg_converged_all_accepted_steps": bool(
                phase_converged.size and np.all(phase_converged)
            ),
            "phase_cg_max_relative_residual": _max_finite(phase_residual),
            "phase_cg_max_iterations": _max_finite(phase_iterations),
        },
    }
    figures_dir.mkdir(parents=True, exist_ok=True)
    out = figures_dir / "E2_phasecg_diagnostics_manifest.yaml"
    with out.open("w", encoding="utf-8") as stream:
        yaml.safe_dump(
            manifest,
            stream,
            sort_keys=False,
            allow_unicode=True,
            default_flow_style=False,
        )
    return out


def _attach_snapshot_power_diagnostics(dataset: dict[str, Any], raw_ss: Path) -> None:
    """Attach the stored thermal balance without recomputing catalogue lookups."""

    path = Path(raw_ss) / "snapshot_power_energy_diagnostics.npz"
    if not path.exists():
        return
    requested = (
        "snapshot_t_ps",
        "joule_snapshot_W_m3",
        "P_S_snapshot_W_m3",
        "P_R_snapshot_W_m3",
        "P_total_snapshot_W_m3",
        "P_esc_snapshot_W_m3",
        "u_e_snapshot_J_m3",
        "u_ph_snapshot_J_m3",
        "C_e_snapshot_J_m3_K",
        "C_ph_snapshot_J_m3_K",
    )
    with np.load(path, allow_pickle=True) as stored:
        for key in requested:
            if key in stored.files:
                dataset[key] = np.asarray(stored[key])


def _last_finite(values: np.ndarray) -> float | None:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    return float(finite[-1]) if finite.size else None


def _max_finite(values: np.ndarray) -> float | None:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    return float(np.nanmax(finite)) if finite.size else None


if __name__ == "__main__":
    raise SystemExit(main())
