"""PHOTON-run placeholder for the future coupled SNSPD transient.

This file intentionally does not claim to simulate detection yet.  It reserves
the official pipeline location and records the inputs that the coupled transient
will need: PRE catalogues, SS stationary state, photon parameters, thermal
closure and circuit parameters.
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import numpy as np
import yaml

from pysnspd.config import load_config, validate_config
from pysnspd.io.manager import create_run_layout, write_manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare a placeholder PHOTON-run manifest for the coupled transient."
    )
    parser.add_argument("--config", required=True, help="Path to YAML project config.")
    parser.add_argument(
        "--run-name",
        default=None,
        help="Photon run name. If omitted, project.default_run_name is used.",
    )
    parser.add_argument(
        "--pre-run-name",
        default=None,
        help="PRE-run name to load. If omitted, use --run-name/default_run_name.",
    )
    parser.add_argument(
        "--ss-run-name",
        default=None,
        help="SS-run name to load. If omitted, use --run-name/default_run_name.",
    )
    parser.add_argument("--photon-energy-eV", type=float, default=1.165, help="Photon energy in eV.")
    parser.add_argument("--absorption-x-fraction", type=float, default=0.50)
    parser.add_argument("--absorption-y-fraction", type=float, default=0.50)
    parser.add_argument("--eta-abs", type=float, default=1.0, help="Absorbed-energy efficiency placeholder.")
    parser.add_argument("--t-final-ps", type=float, default=200.0)
    parser.add_argument("--dt-fs", type=float, default=0.25)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    cfg = validate_config(load_config(args.config))

    photon_layout = create_run_layout(cfg, args.run_name)
    run_name = photon_layout["run_name"]
    pre_name = args.pre_run_name or run_name
    ss_name = args.ss_run_name or run_name

    pre_layout = create_run_layout(cfg, pre_name)
    ss_layout = create_run_layout(cfg, ss_name)

    raw_photon = Path(photon_layout["raw_photon"])
    raw_photon.mkdir(parents=True, exist_ok=True)

    expected_inputs = {
        "pre_mesh_npz": str(Path(pre_layout["raw_pre"]) / "mesh.npz"),
        "pre_edges_npz": str(Path(pre_layout["raw_pre"]) / "edges.npz"),
        "pre_usadel_npz": str(Path(pre_layout["raw_pre"]) / "usadel_dos_catalog.npz"),
        "pre_phase_space_npz": str(Path(pre_layout["raw_pre"]) / "phase_space_catalog.npz"),
        "ss_state_npz": str(Path(ss_layout["raw_ss"]) / "stationary_state.npz"),
        "ss_history_npz": str(Path(ss_layout["raw_ss"]) / "relaxation_history.npz"),
    }

    missing_inputs = [key for key, value in expected_inputs.items() if not Path(value).exists()]
    photon_energy_J = float(args.photon_energy_eV) * 1.602176634e-19

    placeholder_npz = raw_photon / "photon_placeholder_state.npz"
    np.savez_compressed(
        placeholder_npz,
        photon_energy_eV=np.array(float(args.photon_energy_eV)),
        photon_energy_J=np.array(photon_energy_J),
        absorption_x_fraction=np.array(float(args.absorption_x_fraction)),
        absorption_y_fraction=np.array(float(args.absorption_y_fraction)),
        eta_abs=np.array(float(args.eta_abs)),
        t_final_ps=np.array(float(args.t_final_ps)),
        dt_fs=np.array(float(args.dt_fs)),
        status=np.array("placeholder_no_coupled_transient_yet"),
    )

    summary = {
        "run_name": run_name,
        "pre_run_name": pre_name,
        "ss_run_name": ss_name,
        "backend": "photon_placeholder_v1",
        "status": "placeholder_no_coupled_transient_yet",
        "honest_scope": (
            "This pipeline does not yet evolve Te, Tph, Delta, electric potential, "
            "I_SNSPD or V_out. It only fixes the official PHOTON-run interface."
        ),
        "expected_inputs": expected_inputs,
        "missing_inputs": missing_inputs,
        "photon": {
            "energy_eV": float(args.photon_energy_eV),
            "energy_J": photon_energy_J,
            "eta_abs": float(args.eta_abs),
            "absorption_x_fraction": float(args.absorption_x_fraction),
            "absorption_y_fraction": float(args.absorption_y_fraction),
        },
        "time_axis_request": {
            "t_final_ps": float(args.t_final_ps),
            "dt_fs": float(args.dt_fs),
        },
        "outputs": {
            "placeholder_npz": str(placeholder_npz),
        },
    }

    summary_path = raw_photon / "photon_summary.yaml"
    _write_yaml(summary_path, summary)
    manifest_path = write_manifest(
        cfg,
        run_name,
        stage="photon",
        extra={
            "pipeline": "03_photon_run_template.py",
            "purpose": "Placeholder PHOTON-run interface for the future coupled transient.",
            "summary": summary,
            "outputs": summary["outputs"] | {"photon_summary": str(summary_path)},
        },
    )

    print("PHOTON-run placeholder")
    print(f"  run_name:     {run_name}")
    print(f"  pre_run_name: {pre_name}")
    print(f"  ss_run_name:  {ss_name}")
    print(f"  raw_photon:   {raw_photon}")
    print()
    if missing_inputs:
        print("Missing expected inputs")
        for key in missing_inputs:
            print(f"  {key}: {expected_inputs[key]}")
        print()
    print("Photon request")
    print(f"  energy_eV: {args.photon_energy_eV}")
    print(f"  energy_J:  {photon_energy_J:.6e}")
    print(f"  eta_abs:   {args.eta_abs}")
    print()
    print("Outputs")
    print(f"  placeholder_npz: {placeholder_npz}")
    print(f"  photon_summary:  {summary_path}")
    print(f"  photon_manifest: {manifest_path}")
    print("Status: PLACEHOLDER")
    return 0


def _write_yaml(path: str | Path, data: dict[str, Any]) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as f:
        yaml.safe_dump(data, f, sort_keys=False, allow_unicode=True)


if __name__ == "__main__":
    raise SystemExit(main())
