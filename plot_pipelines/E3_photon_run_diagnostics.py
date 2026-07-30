#!/usr/bin/env python3
"""Evaluate one completed photon run without comparing it to another run."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import yaml

from pysnspd.analysis.timing import analyze_photon_timing
from pysnspd.analysis.timing_cli import (
    add_timing_analysis_arguments,
    timing_criteria_from_args,
)
from pysnspd.config import load_config, validate_config
from pysnspd.io.manager import create_run_layout
from pysnspd.mesh.delaunay import load_mesh_npz
from pysnspd.plotting.photon_diagnostics import (
    make_photon_run_diagnostic_figures,
    nearest_unique_snapshot_indices,
)
from pysnspd.plotting.style import THESIS_DPI

HBAR_J_S = 1.054571817e-34
K_B_J_K = 1.380649e-23


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create scalar and field diagnostics for one completed photon run."
    )
    parser.add_argument("--config", required=True, help="Absolute YAML project configuration.")
    parser.add_argument("--pre-run-name", required=True, help="PRE run containing the mesh.")
    parser.add_argument("--run-name", required=True, help="Completed photon run to evaluate.")
    parser.add_argument(
        "--times-ps",
        nargs="+",
        type=float,
        default=(50.0, 55.0, 60.0, 100.0),
        help=(
            "Requested field-map times in ps. The nearest stored snapshot is used; "
            "duplicate resolved snapshots are plotted only once."
        ),
    )
    parser.add_argument(
        "--xi-nm",
        type=float,
        default=None,
        help=(
            "Reference coherence length for |q|*xi in nm. By default it is reconstructed "
            "from the PRE effective gTDGL diffusion coefficient at the bias temperature."
        ),
    )
    parser.add_argument(
        "--figures-subdir",
        default="E3_photon_diagnostics",
        help="Subdirectory below plots/<run>/figures.",
    )
    parser.add_argument("--dpi", type=int, default=THESIS_DPI)
    add_timing_analysis_arguments(parser)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    cfg = validate_config(load_config(args.config))
    pre_layout = create_run_layout(cfg, args.pre_run_name)
    run_layout = create_run_layout(cfg, args.run_name)

    raw_pre = Path(pre_layout["raw_pre"])
    raw_photon = Path(run_layout["raw_photon"])
    output_dir = Path(run_layout["plots_figures"]) / str(args.figures_subdir)
    output_dir.mkdir(parents=True, exist_ok=True)

    mesh = load_mesh_npz(_require_file(raw_pre / "mesh.npz", "PRE mesh"))
    pre_summary = _read_yaml(
        _require_file(raw_pre / "usadel_dos_summary.yaml", "Usadel summary")
    )
    history = _load_npz(
        _require_file(raw_photon / "transient_history.npz", "photon history")
    )
    snapshots = _load_npz(
        _require_file(raw_photon / "transient_snapshots.npz", "photon snapshots")
    )
    summary = _read_yaml(
        _require_file(raw_photon / "photon_summary.yaml", "photon summary")
    )

    delta0_meV = _find_numeric(pre_summary, "delta0_meV")
    if not np.isfinite(delta0_meV) or delta0_meV <= 0.0:
        raise ValueError("PRE metadata do not provide a positive delta0_meV.")
    xi_m, xi_source = _resolve_xi_m(
        cfg=cfg,
        pre_summary=pre_summary,
        snapshots=snapshots,
        override_nm=args.xi_nm,
    )

    detection_criteria, recovery_criteria = timing_criteria_from_args(args)
    timing = analyze_photon_timing(
        history,
        snapshots=snapshots,
        detection=detection_criteria,
        recovery=recovery_criteria,
    )
    saved = make_photon_run_diagnostic_figures(
        mesh=mesh,
        history=history,
        snapshots=snapshots,
        summary=summary,
        delta0_meV=float(delta0_meV),
        xi_m=xi_m,
        requested_times_ps=args.times_ps,
        output_dir=output_dir,
        dpi=int(args.dpi),
        timing=timing,
    )
    manifest_path = _write_manifest(
        args=args,
        raw_pre=raw_pre,
        raw_photon=raw_photon,
        output_dir=output_dir,
        saved=saved,
        delta0_meV=float(delta0_meV),
        xi_m=xi_m,
        xi_source=xi_source,
        snapshots=snapshots,
        timing=timing,
    )

    resolved_times = _resolved_times(snapshots, args.times_ps)
    print("E3 single photon-run diagnostics")
    print(f" run_name:       {args.run_name}")
    print(f" pre_run_name:   {args.pre_run_name}")
    print(f" raw_photon:     {raw_photon}")
    print(f" figures_dir:    {output_dir}")
    print(f" Delta_BCS(0):   {delta0_meV:.9g} meV")
    print(f" xi:             {1.0e9 * xi_m:.9g} nm ({xi_source})")
    print(
        " snapshot_times: "
        f"requested={', '.join(f'{float(t):g}' for t in args.times_ps)} ps; "
        f"resolved={', '.join(f'{t:g}' for t in resolved_times)} ps"
    )
    print(
        " timing:         "
        f"t_lat={dict(timing.get('latency', {})).get('t_lat_ps', 'censored')} ps, "
        f"t_rec={dict(dict(timing.get('recovery', {})).get('selected', {})).get('t_rec_ps', 'censored')} ps"
    )
    print("Figures")
    for key, path in saved.items():
        print(f" {key}: {path}")
    print(f" manifest: {manifest_path}")
    print("Status: OK")
    return 0


def _write_manifest(
    *,
    args: argparse.Namespace,
    raw_pre: Path,
    raw_photon: Path,
    output_dir: Path,
    saved: Mapping[str, Path],
    delta0_meV: float,
    xi_m: float,
    xi_source: str,
    snapshots: Mapping[str, Any],
    timing: Mapping[str, Any],
) -> Path:
    manifest = {
        "schema_version": 1,
        "pipeline": "plot_pipelines/E3_photon_run_diagnostics.py",
        "purpose": "Scalar and selected-field diagnostics for one completed photon run.",
        "run_name": str(args.run_name),
        "pre_run_name": str(args.pre_run_name),
        "raw_pre": str(raw_pre),
        "raw_photon": str(raw_photon),
        "output_dir": str(output_dir),
        "snapshot_selection": {
            "requested_times_ps": [float(value) for value in args.times_ps],
            "resolved_times_ps": _resolved_times(snapshots, args.times_ps),
            "policy": "nearest stored snapshot; duplicate resolved indices omitted",
        },
        "normalizations": {
            "order_parameter": "abs(Delta) / Delta_BCS(0)",
            "superfluid_momentum": "abs(q) * xi",
        },
        "delta0_meV": float(delta0_meV),
        "xi_m": float(xi_m),
        "xi_source": str(xi_source),
        "timing": dict(timing),
        "figures": {key: str(path) for key, path in saved.items()},
    }
    path = output_dir / "E3_photon_diagnostics_manifest.yaml"
    with path.open("w", encoding="utf-8") as stream:
        yaml.safe_dump(
            manifest,
            stream,
            sort_keys=False,
            allow_unicode=True,
            default_flow_style=False,
        )
    return path


def _resolve_xi_m(
    *,
    cfg: Mapping[str, Any],
    pre_summary: Mapping[str, Any],
    snapshots: Mapping[str, Any],
    override_nm: float | None,
) -> tuple[float, str]:
    if override_nm is not None:
        xi_m = float(override_nm) * 1.0e-9
        if not np.isfinite(xi_m) or xi_m <= 0.0:
            raise ValueError("--xi-nm must be finite and positive.")
        return xi_m, "command-line --xi-nm"

    stored = np.asarray(snapshots.get("stationarity_xi_m", []), dtype=float).reshape(-1)
    stored = stored[np.isfinite(stored) & (stored > 0.0)]
    if stored.size:
        return float(stored[0]), "transient_snapshots.npz:stationarity_xi_m"

    diffusion = _first_positive(
        _nested_value(pre_summary, ("usadel", "gtdgl_allmaras", "D_effective_m2_s")),
        _nested_value(pre_summary, ("metadata", "gtdgl_allmaras_D_m2_s")),
        _nested_value(pre_summary, ("usadel", "D_m2_s")),
        _nested_value(cfg, ("material", "D_m2_s")),
    )
    Tc_K = _first_positive(
        _nested_value(pre_summary, ("metadata", "Tc_K")),
        _nested_value(cfg, ("material", "Tc_K")),
    )
    bias_K = _first_positive(
        _nested_value(cfg, ("bias", "T_bias_K")),
        _initial_snapshot_temperature(snapshots),
    )
    if not all(np.isfinite(value) and value > 0.0 for value in (diffusion, Tc_K, bias_K)):
        raise ValueError(
            "Could not reconstruct xi. Supply --xi-nm or provide PRE D, Tc and bias temperature."
        )
    xi2_m2 = (
        math.pi
        * HBAR_J_S
        * diffusion
        / (4.0 * math.sqrt(2.0) * K_B_J_K * Tc_K * math.sqrt(1.0 + bias_K / Tc_K))
    )
    return (
        math.sqrt(xi2_m2),
        "PRE effective gTDGL D and project bias temperature",
    )


def _resolved_times(
    snapshots: Mapping[str, Any],
    requested_times_ps: list[float] | tuple[float, ...],
) -> list[float]:
    if "snapshot_t_ps" in snapshots:
        stored = np.asarray(snapshots["snapshot_t_ps"], dtype=float).reshape(-1)
    else:
        stored = (
            np.asarray(snapshots.get("snapshot_t_s", []), dtype=float).reshape(-1)
            / 1.0e-12
        )
    if stored.size == 0:
        return []
    indices = nearest_unique_snapshot_indices(stored, requested_times_ps)
    return [float(stored[index]) for index in indices]


def _initial_snapshot_temperature(snapshots: Mapping[str, Any]) -> float:
    values = np.asarray(snapshots.get("Te_snapshot_K", []), dtype=float)
    if values.ndim < 2 or values.shape[0] == 0:
        return np.nan
    finite = values[0][np.isfinite(values[0]) & (values[0] > 0.0)]
    return float(np.nanmedian(finite)) if finite.size else np.nan


def _first_positive(*values: Any) -> float:
    for value in values:
        try:
            numeric = float(value)
        except (TypeError, ValueError):
            continue
        if np.isfinite(numeric) and numeric > 0.0:
            return numeric
    return np.nan


def _nested_value(mapping: Mapping[str, Any], path: tuple[str, ...]) -> Any:
    value: Any = mapping
    for key in path:
        if not isinstance(value, Mapping) or key not in value:
            return None
        value = value[key]
    return value


def _find_numeric(value: Any, target: str) -> float:
    if isinstance(value, Mapping):
        if target in value:
            try:
                return float(value[target])
            except (TypeError, ValueError):
                pass
        for nested in value.values():
            found = _find_numeric(nested, target)
            if np.isfinite(found):
                return found
    return np.nan


def _load_npz(path: Path) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=True) as data:
        return {key: np.asarray(data[key]) for key in data.files}


def _read_yaml(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as stream:
        value = yaml.safe_load(stream) or {}
    return value if isinstance(value, dict) else {}


def _require_file(path: Path, label: str) -> Path:
    if not path.exists():
        raise FileNotFoundError(f"Missing {label}: {path}")
    return path


if __name__ == "__main__":
    raise SystemExit(main())
