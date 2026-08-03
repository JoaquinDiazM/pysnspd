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
from pysnspd.analysis.photon_snapshots import compute_photon_snapshot_plot_diagnostics
from pysnspd.analysis.timing_cli import (
    add_timing_analysis_arguments,
    timing_criteria_from_args,
)
from pysnspd.config import load_config, validate_config
from pysnspd.io.manager import create_run_layout
from pysnspd.mesh.delaunay import load_mesh_npz
from pysnspd.mesh.edges import load_edges_npz
from pysnspd.mesh.operators import build_fv_operators
from pysnspd.plotting.photon_diagnostics import (
    make_photon_run_diagnostic_figures,
    nearest_unique_snapshot_indices,
)
from pysnspd.plotting.style import THESIS_DPI
from pysnspd.thermal.evolution import build_central_thermal_mask

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
    edge_data = load_edges_npz(_require_file(raw_pre / "edges.npz", "PRE edges"))
    pre_summary = _read_yaml(
        _require_file(raw_pre / "usadel_dos_summary.yaml", "Usadel summary")
    )
    history = _load_npz(
        _require_file(raw_photon / "transient_history.npz", "photon history")
    )
    snapshots_path = _require_file(
        raw_photon / "transient_snapshots.npz", "photon snapshots"
    )
    power_table_path = _require_file(
        raw_pre / "power_table_catalog.npz", "PRE power table"
    )
    snapshots = _load_npz(snapshots_path)
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
    selected_indices = nearest_unique_snapshot_indices(
        _snapshot_times(snapshots), args.times_ps
    )
    transient_config = dict(summary.get("config", {}))
    thermal_enabled = bool(transient_config.get("thermal_enabled", False))
    thermal_window_m = float(transient_config.get("thermal_window_m", 100.0e-9))
    thermal_mask = (
        build_central_thermal_mask(
            np.asarray(mesh.nodes, dtype=float),
            window_m=thermal_window_m,
        )
        if thermal_enabled
        else np.zeros(np.asarray(mesh.nodes).shape[0], dtype=bool)
    )
    snapshot_diagnostics = compute_photon_snapshot_plot_diagnostics(
        snapshots=snapshots,
        selected_indices=selected_indices,
        power_table_npz=str(power_table_path),
        ops=build_fv_operators(mesh, edge_data),
        sigma_n_S_m=float(cfg["material"]["sigma_n_S_m"]),
        thermal_active_mask=thermal_mask,
        thermal_bath_K=_initial_snapshot_field_temperature(
            snapshots, "Tph_snapshot_K"
        ),
        global_limits_override=_reusable_global_limits(
            output_dir / "E3_photon_diagnostics_manifest.yaml",
            snapshots_path=snapshots_path,
            power_table_path=power_table_path,
            requested_times_ps=args.times_ps,
        ),
        progress_callback=_snapshot_diagnostic_progress(),
    )
    snapshot_diagnostics.update(
        {
            "nodes_x_nm": 1.0e9 * np.asarray(mesh.nodes, dtype=float)[:, 0],
            "nodes_y_nm": 1.0e9 * np.asarray(mesh.nodes, dtype=float)[:, 1],
            "triangles": np.asarray(mesh.triangles, dtype=np.int64),
        }
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
        snapshot_diagnostics=snapshot_diagnostics,
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
        snapshot_diagnostics=snapshot_diagnostics,
        snapshots_path=snapshots_path,
        power_table_path=power_table_path,
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
    print(
        " global_scales:  "
        + (
            "reused from matching manifest"
            if bool(snapshot_diagnostics.get("global_limits_reused", False))
            else "scanned over all persisted snapshots"
        )
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
    snapshot_diagnostics: Mapping[str, Any],
    snapshots_path: Path,
    power_table_path: Path,
) -> Path:
    manifest = {
        "schema_version": 2,
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
        "snapshot_scale_policy": str(
            snapshot_diagnostics.get("global_scale_policy", "unavailable")
        ),
        "snapshot_scale_schema": "photon_snapshot_thermodynamics_v1_multilinear",
        "snapshot_global_limits": {
            key: [float(value) for value in np.asarray(bounds, dtype=float).reshape(-1)]
            for key, bounds in dict(
                snapshot_diagnostics.get("snapshot_global_limits", {})
            ).items()
        },
        "source_fingerprints": {
            "transient_snapshots": _file_fingerprint(snapshots_path),
            "power_table_catalog": _file_fingerprint(power_table_path),
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


def _initial_snapshot_field_temperature(
    snapshots: Mapping[str, Any], key: str
) -> float:
    values = np.asarray(snapshots.get(key, []), dtype=float)
    if values.ndim != 2 or values.shape[0] == 0:
        raise ValueError(f"Photon snapshots lack initial {key} values.")
    finite = values[0][np.isfinite(values[0]) & (values[0] > 0.0)]
    if finite.size == 0:
        raise ValueError(f"Photon snapshots have no positive finite initial {key} values.")
    return float(np.nanmedian(finite))


def _snapshot_times(snapshots: Mapping[str, Any]) -> np.ndarray:
    if "snapshot_t_ps" in snapshots:
        return np.asarray(snapshots["snapshot_t_ps"], dtype=float).reshape(-1)
    return np.asarray(snapshots.get("snapshot_t_s", []), dtype=float).reshape(-1) / 1.0e-12


def _snapshot_diagnostic_progress():
    state = {"bucket": -1}

    def report(completed: int, total: int) -> None:
        if total <= 0:
            return
        bucket = min(10, int(10 * completed / total))
        if bucket > state["bucket"]:
            state["bucket"] = bucket
            print(f" snapshot thermodynamics: {10 * bucket}% ({completed}/{total})")

    return report


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


def _reusable_global_limits(
    manifest_path: Path,
    *,
    snapshots_path: Path,
    power_table_path: Path,
    requested_times_ps: list[float] | tuple[float, ...],
) -> dict[str, np.ndarray] | None:
    """Reuse a completed all-snapshot scan only when its inputs still match."""

    if not manifest_path.exists():
        return None
    try:
        manifest = _read_yaml(manifest_path)
        if int(manifest.get("schema_version", 0)) < 2:
            return None
        scale_schema = manifest.get("snapshot_scale_schema")
        if scale_schema not in (
            None,
            "photon_snapshot_thermodynamics_v1_multilinear",
        ):
            return None
        previous = np.asarray(
            dict(manifest.get("snapshot_selection", {})).get(
                "requested_times_ps", []
            ),
            dtype=float,
        )
        requested = np.asarray(requested_times_ps, dtype=float)
        if previous.shape != requested.shape or not np.allclose(
            previous, requested, rtol=0.0, atol=1.0e-12
        ):
            return None
        fingerprints = dict(manifest.get("source_fingerprints", {}))
        if fingerprints:
            if fingerprints.get("transient_snapshots") != _file_fingerprint(
                snapshots_path
            ):
                return None
            if fingerprints.get("power_table_catalog") != _file_fingerprint(
                power_table_path
            ):
                return None
        elif manifest_path.stat().st_mtime_ns < max(
            snapshots_path.stat().st_mtime_ns,
            power_table_path.stat().st_mtime_ns,
        ):
            return None
        raw_limits = dict(manifest.get("snapshot_global_limits", {}))
        required = {
            "q_abs_snapshot_m_inv",
            "joule_snapshot_W_m3",
            "P_S_snapshot_W_m3",
            "P_R_snapshot_W_m3",
            "P_total_snapshot_W_m3",
            "P_diff_snapshot_W_m3",
            "P_esc_snapshot_W_m3",
            "u_e_snapshot_J_m3",
            "u_ph_snapshot_J_m3",
            "C_e_snapshot_J_m3_K",
            "C_ph_snapshot_J_m3_K",
            "kappa_s_snapshot_W_m_K",
        }
        if not required.issubset(raw_limits):
            return None
        limits = {
            key: np.asarray(raw_limits[key], dtype=float).reshape(-1)
            for key in required
        }
        if not all(values.size == 2 and np.all(np.isfinite(values)) for values in limits.values()):
            return None
        return limits
    except (OSError, TypeError, ValueError):
        return None


def _file_fingerprint(path: Path) -> dict[str, int]:
    stat = path.stat()
    return {"size_bytes": int(stat.st_size), "mtime_ns": int(stat.st_mtime_ns)}


if __name__ == "__main__":
    raise SystemExit(main())
