"""Photon + circuit transient pipeline.

This pipeline starts from a validated SS run, initializes a lumped readout
circuit at the SS fixed point, optionally injects a phonon bubble at a requested
time, and continues the gTDGL/thermal evolution with a circuit-updated terminal
current.

The first implementation is intentionally conservative: it reuses the existing
``solve_stationary_pytdgl_like`` adapter in short chunks instead of modifying
the validated SS core.
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
from pysnspd.gtdgl.material import build_gtdgl_material
from pysnspd.mesh.operators import build_fv_operators
from pysnspd.excitation.photon import PhotonBubbleParams
from pysnspd.solver.transient import CoupledTransientConfig, run_coupled_transient
from pysnspd.analysis.timing import DetectionCriteria, RecoveryCriteria
from pysnspd.gtdgl.usadel_current import (
    attach_usadel_supercurrent_table_from_npz,
    validate_strict_usadel_supercurrent_table_npz,
)
from pysnspd.circuit.readout import CircuitParams


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run coupled gTDGL/thermal + circuit + photon transient from an SS state."
    )
    parser.add_argument("--config", required=True, help="Path to YAML project config.")
    parser.add_argument("--run-name", required=True, help="Output run name for pipeline 03.")
    parser.add_argument("--pre-run-name", required=True, help="PRE run name used by the SS run.")
    parser.add_argument("--ss-run-name", required=True, help="Existing SS run name to initialize from.")

    parser.add_argument("--total-time-ps", type=float, default=1500.0)
    parser.add_argument("--gtdgl-dt-fs", type=float, default=0.75)
    parser.add_argument(
        "--coupling-step-fs",
        type=float,
        default=100.0,
        help="Chunk size for gTDGL/circuit splitting. Keep this small for production.",
    )
    parser.add_argument(
        "--snapshots",
        type=int,
        default=None,
        help="Stored snapshots. Default: total_time_ps * 10 (about 10 snapshots/ps).",
    )
    parser.add_argument("--progress", dest="progress", action="store_true", default=True)
    parser.add_argument("--no-progress", dest="progress", action="store_false")

    parser.add_argument("--center-voltage-width-nm", type=float, default=100.0)
    parser.add_argument("--center-voltage-probe-band-nm", type=float, default=None)

    parser.add_argument(
        "--thermal-enable",
        dest="thermal_enable",
        action="store_true",
        default=True,
        help="Enable runtime Te/Tph evolution. Enabled by default.",
    )
    parser.add_argument(
        "--thermal-disable",
        dest="thermal_enable",
        action="store_false",
        help="Freeze Te/Tph during mesoscopic chunks; a photon event can still set the initial Tph bubble.",
    )
    parser.add_argument("--thermal-window-nm", type=float, default=100.0)
    parser.add_argument("--thermal-max-step-K", type=float, default=0.20)
    parser.add_argument("--thermal-max-substeps", type=int, default=32)

    parser.add_argument("--allmaras-direct-amplitude-fraction", type=float, default=2.0e-2)
    parser.add_argument("--allmaras-convergence-tol", type=float, default=3.0e-3)
    parser.add_argument("--allmaras-convergence-max-iterations", type=int, default=32)

    parser.add_argument("--terminal-psi", type=float, default=0.0)
    parser.add_argument(
        "--terminal-healing-xi",
        type=float,
        default=None,
        help="Default None avoids re-healing the already converged SS initial state.",
    )
    parser.add_argument("--terminal-healing-fraction", type=float, default=0.95)

    parser.add_argument("--circuit-Rload-ohm", type=float, default=50.0)
    parser.add_argument("--circuit-Rbias-ohm", type=float, default=1.0e4)
    parser.add_argument("--circuit-Lbias-nH", type=float, default=1000.0)
    parser.add_argument("--circuit-Lk-nH", type=float, default=10.0)
    parser.add_argument("--circuit-Ccouple-pF", type=float, default=100.0)
    parser.add_argument(
        "--circuit-Vbias-uV",
        type=float,
        default=None,
        help="Optional explicit Vbias. If omitted, initialize at exact SS fixed point.",
    )

    parser.add_argument("--photon-enable", dest="photon_enable", action="store_true", default=True)
    parser.add_argument("--photon-disable", dest="photon_enable", action="store_false")
    parser.add_argument("--photon-time-ps", type=float, default=50.0)
    parser.add_argument("--photon-energy-eV", type=float, default=0.8)
    parser.add_argument("--photon-x-nm", type=float, default=None, help="Default: mesh center.")
    parser.add_argument("--photon-y-nm", type=float, default=0.0)
    parser.add_argument("--photon-sigma-nm", type=float, default=10.0)

    parser.add_argument(
        "--early-stop-mode",
        choices=("none", "latency", "recovery"),
        default="recovery",
        help=(
            "none: run to total time; latency: stop after a confirmed V_out peak plus "
            "the safety tail; recovery: stop when the selected recovery criterion is held."
        ),
    )
    parser.add_argument(
        "--recovery-mode",
        choices=("electrical", "efficiency90", "state"),
        default="electrical",
        help="Recovery classifier used by --early-stop-mode recovery.",
    )
    parser.add_argument("--detection-threshold-uV", type=float, default=100.0)
    parser.add_argument(
        "--detection-polarity",
        choices=("positive", "negative", "auto"),
        default="positive",
    )
    parser.add_argument("--detection-baseline-window-ps", type=float, default=10.0)
    parser.add_argument("--detection-confirmation-ps", type=float, default=0.5)
    parser.add_argument("--detection-hysteresis-fraction", type=float, default=0.10)
    parser.add_argument("--peak-confirmation-ps", type=float, default=2.0)
    parser.add_argument("--post-peak-safety-ps", type=float, default=10.0)
    parser.add_argument("--recovery-hold-ps", type=float, default=10.0)
    parser.add_argument("--recovery-efficiency-fraction", type=float, default=0.90)
    parser.add_argument("--recovery-current-rel-tol", type=float, default=1.0e-2)
    parser.add_argument("--recovery-current-abs-uA", type=float, default=0.05)
    parser.add_argument("--recovery-voltage-rel-tol", type=float, default=1.0e-2)
    parser.add_argument("--recovery-voltage-abs-uV", type=float, default=10.0)
    parser.add_argument("--recovery-temperature-abs-K", type=float, default=0.05)
    parser.add_argument("--recovery-condensate-rel-tol", type=float, default=2.0e-2)
    parser.add_argument("--recovery-spatial-quantile", type=float, default=0.995)
    parser.add_argument("--recovery-spatial-max-guard-factor", type=float, default=4.0)
    parser.add_argument("--timing-evaluation-interval-ps", type=float, default=0.5)

    return parser.parse_args()


def main() -> int:
    args = parse_args()
    cfg = validate_config(load_config(args.config))
    n_snapshots = (
        int(args.snapshots)
        if args.snapshots is not None
        else max(2, int(round(float(args.total_time_ps) * 10.0)))
    )

    out_layout = create_run_layout(cfg, args.run_name)
    pre_layout = create_run_layout(cfg, args.pre_run_name)
    ss_layout = create_run_layout(cfg, args.ss_run_name)

    raw_pre = Path(pre_layout["raw_pre"])
    raw_ss = Path(ss_layout["raw_ss"])
    raw_photon = Path(out_layout["raw_ss"]).parent / "photon"
    raw_photon.mkdir(parents=True, exist_ok=True)

    mesh_path = raw_pre / "mesh.npz"
    edges_path = raw_pre / "edges.npz"
    usadel_path = raw_pre / "usadel_dos_catalog.npz"
    power_table_path = raw_pre / "power_table_catalog.npz"
    initial_state_path = raw_ss / "stationary_state.npz"
    ss_summary_path = raw_ss / "ss_summary.yaml"

    _require_file(mesh_path, "PRE mesh")
    _require_file(edges_path, "PRE edges")
    _require_file(usadel_path, "PRE Usadel catalogue")
    _require_file(initial_state_path, "SS stationary state")

    mesh = load_mesh_npz(mesh_path)
    edge_data = load_edges_npz(edges_path)
    ops = build_fv_operators(mesh, edge_data)

    strict_table_summary = validate_strict_usadel_supercurrent_table_npz(usadel_path)
    base_usadel_catalog = load_usadel_catalog_npz(usadel_path)
    usadel_catalog = attach_usadel_supercurrent_table_from_npz(base_usadel_catalog, usadel_path)

    allmaras_diffusion = _read_pre_allmaras_diffusion(raw_pre)
    material = build_gtdgl_material(
        cfg,
        base_usadel_catalog,
        diffusion_factor=float(allmaras_diffusion["D_effective_factor"]),
    )

    ss_summary = _read_yaml(ss_summary_path)
    initial_current_A = _initial_current_A(ss_summary)

    circuit_params = CircuitParams(
        R_load_ohm=float(args.circuit_Rload_ohm),
        R_bias_ohm=float(args.circuit_Rbias_ohm),
        L_bias_H=float(args.circuit_Lbias_nH) * 1.0e-9,
        L_k_H=float(args.circuit_Lk_nH) * 1.0e-9,
        C_couple_F=float(args.circuit_Ccouple_pF) * 1.0e-12,
        V_bias_V=(None if args.circuit_Vbias_uV is None else float(args.circuit_Vbias_uV) * 1.0e-6),
    )

    photon_params = PhotonBubbleParams(
        enabled=bool(args.photon_enable),
        energy_eV=float(args.photon_energy_eV),
        time_s=float(args.photon_time_ps) * 1.0e-12,
        x_m=(None if args.photon_x_nm is None else float(args.photon_x_nm) * 1.0e-9),
        y_m=float(args.photon_y_nm) * 1.0e-9,
        sigma_m=float(args.photon_sigma_nm) * 1.0e-9,
    )

    transient_config = CoupledTransientConfig(
        total_time_s=float(args.total_time_ps) * 1.0e-12,
        mesoscopic_dt_s=float(args.gtdgl_dt_fs) * 1.0e-15,
        chunk_time_s=float(args.coupling_step_fs) * 1.0e-15,
        n_snapshots=n_snapshots,
        center_voltage_width_m=float(args.center_voltage_width_nm) * 1.0e-9,
        center_voltage_probe_band_m=(
            None if args.center_voltage_probe_band_nm is None else float(args.center_voltage_probe_band_nm) * 1.0e-9
        ),
        thermal_enabled=bool(args.thermal_enable),
        thermal_window_m=float(args.thermal_window_nm) * 1.0e-9,
        thermal_max_step_K=float(args.thermal_max_step_K),
        thermal_max_substeps=int(args.thermal_max_substeps),
        terminal_psi=float(args.terminal_psi),
        terminal_healing_xi=args.terminal_healing_xi,
        terminal_healing_fraction=float(args.terminal_healing_fraction),
        supercurrent_law="usadel_poisson",
        allmaras_phase_direct_amplitude_fraction=float(args.allmaras_direct_amplitude_fraction),
        allmaras_phase_convergence_tol=float(args.allmaras_convergence_tol),
        allmaras_phase_convergence_max_iterations=int(args.allmaras_convergence_max_iterations),
        early_stop_mode=str(args.early_stop_mode),
        timing_evaluation_interval_s=float(args.timing_evaluation_interval_ps) * 1.0e-12,
        detection_criteria=DetectionCriteria(
            threshold_V=float(args.detection_threshold_uV) * 1.0e-6,
            polarity=str(args.detection_polarity),
            baseline_window_s=float(args.detection_baseline_window_ps) * 1.0e-12,
            confirmation_s=float(args.detection_confirmation_ps) * 1.0e-12,
            hysteresis_fraction=float(args.detection_hysteresis_fraction),
            peak_confirmation_s=float(args.peak_confirmation_ps) * 1.0e-12,
            post_peak_safety_s=float(args.post_peak_safety_ps) * 1.0e-12,
        ),
        recovery_criteria=RecoveryCriteria(
            mode=str(args.recovery_mode),
            hold_s=float(args.recovery_hold_ps) * 1.0e-12,
            efficiency_fraction=float(args.recovery_efficiency_fraction),
            current_relative_tolerance=float(args.recovery_current_rel_tol),
            current_absolute_tolerance_A=float(args.recovery_current_abs_uA) * 1.0e-6,
            voltage_relative_tolerance=float(args.recovery_voltage_rel_tol),
            voltage_absolute_tolerance_V=float(args.recovery_voltage_abs_uV) * 1.0e-6,
            temperature_absolute_tolerance_K=float(args.recovery_temperature_abs_K),
            condensate_relative_tolerance=float(args.recovery_condensate_rel_tol),
            spatial_quantile=float(args.recovery_spatial_quantile),
            spatial_max_guard_factor=float(args.recovery_spatial_max_guard_factor),
        ),
        progress=bool(args.progress),
    )

    summary = run_coupled_transient(
        mesh=mesh,
        edge_data=edge_data,
        ops=ops,
        material=material,
        initial_state_npz=initial_state_path,
        initial_current_A=float(initial_current_A),
        usadel_catalog=usadel_catalog,
        power_table_npz=(power_table_path if power_table_path.exists() else None),
        output_dir=raw_photon,
        config=transient_config,
        circuit_params=circuit_params,
        photon_params=photon_params,
    )

    manifest_path = write_manifest(
        cfg,
        args.run_name,
        stage="photon",
        extra={
            "pipeline": "03_photon_run_template.py",
            "purpose": "Coupled gTDGL/thermal transient with readout circuit and optional phonon bubble.",
            "pre_run_name": args.pre_run_name,
            "ss_run_name": args.ss_run_name,
            "raw_photon": str(raw_photon),
            "strict_usadel_current_table": strict_table_summary,
            "gtdgl_allmaras_diffusion": allmaras_diffusion,
            "summary": summary,
        },
    )

    print("Photon/circuit transient pipeline")
    print(f" run_name: {args.run_name}")
    print(f" pre_run_name: {args.pre_run_name}")
    print(f" ss_run_name: {args.ss_run_name}")
    print(f" raw_photon: {raw_photon}")
    print()
    print("Run summary")
    print(f" total_time_ps: {float(args.total_time_ps):.6g}")
    print(f" coupling_step_fs: {float(args.coupling_step_fs):.6g}")
    print(f" gtdgl_dt_fs: {float(args.gtdgl_dt_fs):.6g}")
    print(f" snapshots: {n_snapshots}")
    print(f" thermal_enabled: {bool(args.thermal_enable)}")
    print(f" initial_current_uA: {float(initial_current_A) * 1.0e6:.6g}")
    print(f" initial_V_tdgl_center_uV: {float(summary['initial_V_tdgl_center_V']) * 1.0e6:.6g}")
    print(f" photon_energy_eV: {float(args.photon_energy_eV):.6g}")
    print(f" photon_time_ps: {float(args.photon_time_ps):.6g}")
    print(f" early_stop_mode: {args.early_stop_mode}")
    print(f" recovery_mode: {args.recovery_mode}")
    print(f" stop_reason: {summary.get('stop_reason', 'n/a')}")
    timing = dict(summary.get("timing", {}))
    print(f" t_lat_ps: {dict(timing.get('latency', {})).get('t_lat_ps', 'n/a')}")
    print(
        " t_rec_ps: "
        f"{dict(dict(timing.get('recovery', {})).get('selected', {})).get('t_rec_ps', 'censored')}"
    )
    print()
    print("Outputs")
    for key, value in dict(summary.get("outputs", {})).items():
        print(f" {key}: {value}")
    print(f" photon_summary_yaml: {raw_photon / 'photon_summary.yaml'}")
    print(f" manifest: {manifest_path}")
    print("Status: OK")
    return 0


def _initial_current_A(ss_summary: dict[str, Any]) -> float:
    for path in (
        ("solver", "circuit_runtime", "final_state", "I_s_A"),
        ("solver", "target_current_A"),
        ("seed", "simulation_target_current_A"),
        ("target_current_A",),
    ):
        cur: Any = ss_summary
        ok = True
        for key in path:
            if isinstance(cur, dict) and key in cur:
                cur = cur[key]
            else:
                ok = False
                break
        if ok:
            value = float(cur)
            if np.isfinite(value):
                return value
    raise ValueError("Could not determine target current from ss_summary.yaml.")


def _read_pre_allmaras_diffusion(raw_pre: Path) -> dict[str, float | str]:
    summary_path = raw_pre / "usadel_dos_summary.yaml"
    if not summary_path.exists():
        return {
            "D_effective_factor": 1.0,
            "D_base_m2_s": float("nan"),
            "D_effective_m2_s": float("nan"),
            "source": "default: PRE summary not found; using Usadel D unchanged for gTDGL.",
        }
    data = _read_yaml(summary_path)
    usadel = data.get("usadel", {}) if isinstance(data, dict) else {}
    allmaras = usadel.get("gtdgl_allmaras", {}) if isinstance(usadel, dict) else {}
    metadata = data.get("metadata", {}) if isinstance(data, dict) else {}
    base_D = float(allmaras.get("D_base_m2_s", metadata.get("D_m2_s", float("nan"))))
    factor = float(allmaras.get("D_effective_factor", metadata.get("gtdgl_allmaras_D_factor", 1.0)))
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


def _read_yaml(path: str | Path) -> dict[str, Any]:
    p = Path(path)
    if not p.exists():
        return {}
    with p.open("r", encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
    return data if isinstance(data, dict) else {}


def _require_file(path: Path, description: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Missing {description}: {path}")


if __name__ == "__main__":
    raise SystemExit(main())
