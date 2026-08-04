#!/usr/bin/env python3
"""Compare two completed photon runs that differ only in impact position."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import yaml

from pysnspd.analysis.photon_snapshots import compute_photon_snapshot_plot_diagnostics
from pysnspd.config import load_config, validate_config
from pysnspd.io.manager import create_run_layout
from pysnspd.mesh.delaunay import load_mesh_npz
from pysnspd.mesh.edges import load_edges_npz
from pysnspd.mesh.operators import build_fv_operators
from pysnspd.plotting.photon_comparison import make_photon_position_figures
from pysnspd.plotting.photon_diagnostics import nearest_unique_snapshot_indices
from pysnspd.plotting.style import THESIS_DPI
from pysnspd.thermal.evolution import build_central_thermal_mask
from pysnspd.analysis.timing import analyze_photon_timing
from pysnspd.analysis.timing_cli import (
    add_timing_analysis_arguments,
    timing_criteria_from_args,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create thesis-ready PDF comparisons from two completed photon runs."
    )
    parser.add_argument("--config", required=True, help="Absolute YAML project configuration.")
    parser.add_argument("--pre-run-name", required=True, help="PRE run containing the common mesh.")
    parser.add_argument("--center-run-name", required=True, help="Completed central-impact photon run.")
    parser.add_argument("--edge-run-name", required=True, help="Completed edge-impact photon run.")
    parser.add_argument(
        "--times-ps",
        nargs="+",
        type=float,
        default=(50.0, 51.0, 52.0, 53.0),
        help=(
            "Requested matched field-map times; the nearest stored snapshot is used "
            "and duplicate resolved snapshots are omitted."
        ),
    )
    parser.add_argument(
        "--xi-nm",
        type=float,
        default=None,
        help=(
            "Reference coherence length for |q|*xi in nm. By default it is taken "
            "from the photon snapshots or reconstructed from the PRE metadata."
        ),
    )
    parser.add_argument(
        "--output-run-name",
        default="E3_photon_position_comparison",
        help="Run-like directory under big_data_root/plots.",
    )
    parser.add_argument("--figures-subdir", default="figures")
    parser.add_argument("--dpi", type=int, default=THESIS_DPI)
    add_timing_analysis_arguments(parser)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    cfg = validate_config(load_config(args.config))
    pre_layout = create_run_layout(cfg, args.pre_run_name)
    center_layout = create_run_layout(cfg, args.center_run_name)
    edge_layout = create_run_layout(cfg, args.edge_run_name)
    output_layout = create_run_layout(cfg, args.output_run_name)

    raw_pre = Path(pre_layout["raw_pre"])
    center_raw = Path(center_layout["raw_photon"])
    edge_raw = Path(edge_layout["raw_photon"])
    output_dir = Path(output_layout["plots_run"]) / str(args.figures_subdir)
    output_dir.mkdir(parents=True, exist_ok=True)

    mesh = load_mesh_npz(_require_file(raw_pre / "mesh.npz", "PRE mesh"))
    edge_data = load_edges_npz(_require_file(raw_pre / "edges.npz", "PRE edges"))
    ops = build_fv_operators(mesh, edge_data)
    power_table_path = _require_file(
        raw_pre / "power_table_catalog.npz", "PRE power table"
    )
    pre_summary = _read_yaml(
        _require_file(raw_pre / "usadel_dos_summary.yaml", "Usadel summary")
    )
    center_history = _load_npz(_require_file(center_raw / "transient_history.npz", "center history"))
    center_snapshots_path = _require_file(
        center_raw / "transient_snapshots.npz", "center snapshots"
    )
    center_snapshots = _load_npz(center_snapshots_path)
    center_summary = _read_yaml(_require_file(center_raw / "photon_summary.yaml", "center summary"))
    edge_history = _load_npz(_require_file(edge_raw / "transient_history.npz", "edge history"))
    edge_snapshots_path = _require_file(
        edge_raw / "transient_snapshots.npz", "edge snapshots"
    )
    edge_snapshots = _load_npz(edge_snapshots_path)
    edge_summary = _read_yaml(_require_file(edge_raw / "photon_summary.yaml", "edge summary"))
    delta0_meV = _read_delta0_meV(raw_pre)
    xi_m, xi_source = _resolve_xi_m(
        cfg=cfg,
        pre_summary=pre_summary,
        snapshots=center_snapshots,
        override_nm=args.xi_nm,
    )
    detection_criteria, recovery_criteria = timing_criteria_from_args(args)
    center_timing = analyze_photon_timing(
        center_history,
        snapshots=center_snapshots,
        detection=detection_criteria,
        recovery=recovery_criteria,
    )
    edge_timing = analyze_photon_timing(
        edge_history,
        snapshots=edge_snapshots,
        detection=detection_criteria,
        recovery=recovery_criteria,
    )
    center_diagnostics = _snapshot_diagnostics(
        label="center",
        cfg=cfg,
        mesh=mesh,
        ops=ops,
        snapshots=center_snapshots,
        snapshots_path=center_snapshots_path,
        summary=center_summary,
        requested_times_ps=args.times_ps,
        power_table_path=power_table_path,
        existing_manifest=(
            Path(center_layout["plots_figures"])
            / "E3_photon_diagnostics"
            / "E3_photon_diagnostics_manifest.yaml"
        ),
    )
    edge_diagnostics = _snapshot_diagnostics(
        label="edge",
        cfg=cfg,
        mesh=mesh,
        ops=ops,
        snapshots=edge_snapshots,
        snapshots_path=edge_snapshots_path,
        summary=edge_summary,
        requested_times_ps=args.times_ps,
        power_table_path=power_table_path,
        existing_manifest=(
            Path(edge_layout["plots_figures"])
            / "E3_photon_diagnostics"
            / "E3_photon_diagnostics_manifest.yaml"
        ),
    )

    saved = make_photon_position_figures(
        mesh=mesh,
        center_history=center_history,
        center_snapshots=center_snapshots,
        center_summary=center_summary,
        edge_history=edge_history,
        edge_snapshots=edge_snapshots,
        edge_summary=edge_summary,
        delta0_meV=delta0_meV,
        xi_m=xi_m,
        requested_times_ps=args.times_ps,
        output_dir=output_dir,
        dpi=int(args.dpi),
        center_timing=center_timing,
        edge_timing=edge_timing,
        center_snapshot_diagnostics=center_diagnostics,
        edge_snapshot_diagnostics=edge_diagnostics,
    )
    manifest_path = _write_manifest(
        args=args,
        raw_pre=raw_pre,
        center_raw=center_raw,
        edge_raw=edge_raw,
        output_dir=output_dir,
        saved=saved,
        delta0_meV=delta0_meV,
        xi_m=xi_m,
        xi_source=xi_source,
        center_history=center_history,
        edge_history=edge_history,
        center_snapshots=center_snapshots,
        edge_snapshots=edge_snapshots,
        center_timing=center_timing,
        edge_timing=edge_timing,
        center_snapshot_diagnostics=center_diagnostics,
        edge_snapshot_diagnostics=edge_diagnostics,
        center_snapshots_path=center_snapshots_path,
        edge_snapshots_path=edge_snapshots_path,
        power_table_path=power_table_path,
    )

    print("E3 photon-position comparison")
    print(f" center_run: {args.center_run_name}")
    print(f" edge_run:   {args.edge_run_name}")
    print(f" output_dir: {output_dir}")
    print(f" Delta_BCS(0): {delta0_meV:.9g} meV")
    print(f" xi:           {1.0e9 * xi_m:.9g} nm ({xi_source})")
    print(
        " global_scales: "
        f"center={'reused' if center_diagnostics.get('global_limits_reused') else 'scanned'}, "
        f"edge={'reused' if edge_diagnostics.get('global_limits_reused') else 'scanned'}"
    )
    print(
        " center timing: "
        f"t_lat={dict(center_timing.get('latency', {})).get('t_lat_ps', 'censored')} ps, "
        f"t_rec={dict(dict(center_timing.get('recovery', {})).get('selected', {})).get('t_rec_ps', 'censored')} ps"
    )
    print(
        " edge timing:   "
        f"t_lat={dict(edge_timing.get('latency', {})).get('t_lat_ps', 'censored')} ps, "
        f"t_rec={dict(dict(edge_timing.get('recovery', {})).get('selected', {})).get('t_rec_ps', 'censored')} ps"
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
    center_raw: Path,
    edge_raw: Path,
    output_dir: Path,
    saved: Mapping[str, Path],
    delta0_meV: float,
    xi_m: float,
    xi_source: str,
    center_history: Mapping[str, Any],
    edge_history: Mapping[str, Any],
    center_snapshots: Mapping[str, Any],
    edge_snapshots: Mapping[str, Any],
    center_timing: Mapping[str, Any],
    edge_timing: Mapping[str, Any],
    center_snapshot_diagnostics: Mapping[str, Any],
    edge_snapshot_diagnostics: Mapping[str, Any],
    center_snapshots_path: Path,
    edge_snapshots_path: Path,
    power_table_path: Path,
) -> Path:
    manifest = {
        "schema_version": 2,
        "pipeline": "plot_pipelines/E3_photon_position_comparison.py",
        "purpose": "Matched center/edge photon-impact field and circuit-response PDFs.",
        "pre_run_name": str(args.pre_run_name),
        "center_run_name": str(args.center_run_name),
        "edge_run_name": str(args.edge_run_name),
        "raw_pre": str(raw_pre),
        "center_raw_photon": str(center_raw),
        "edge_raw_photon": str(edge_raw),
        "output_dir": str(output_dir),
        "requested_times_ps": [float(value) for value in args.times_ps],
        "center_resolved_times_ps": _resolved_times(center_snapshots, args.times_ps),
        "edge_resolved_times_ps": _resolved_times(edge_snapshots, args.times_ps),
        "delta0_meV": float(delta0_meV),
        "xi_m": float(xi_m),
        "xi_source": str(xi_source),
        "normalizations": {
            "order_parameter": "abs(Delta) / Delta_BCS(0)",
            "superfluid_momentum": "abs(q) * xi",
        },
        "snapshot_scale_policy": (
            "shared exact finite extrema over all persisted snapshots in both runs"
        ),
        "snapshot_global_limits": {
            "center": _serializable_limits(center_snapshot_diagnostics),
            "edge": _serializable_limits(edge_snapshot_diagnostics),
        },
        "source_fingerprints": {
            "center_transient_snapshots": _file_fingerprint(center_snapshots_path),
            "edge_transient_snapshots": _file_fingerprint(edge_snapshots_path),
            "power_table_catalog": _file_fingerprint(power_table_path),
        },
        "center_photon_time_ps": _photon_time(center_history),
        "edge_photon_time_ps": _photon_time(edge_history),
        "center_timing": dict(center_timing),
        "edge_timing": dict(edge_timing),
        "figures": {key: str(path) for key, path in saved.items()},
    }
    path = output_dir / "E3_photon_position_manifest.yaml"
    with path.open("w", encoding="utf-8") as stream:
        yaml.safe_dump(manifest, stream, sort_keys=False, allow_unicode=True, default_flow_style=False)
    return path


def _read_delta0_meV(raw_pre: Path) -> float:
    summary = _read_yaml(_require_file(raw_pre / "usadel_dos_summary.yaml", "Usadel summary"))
    value = _find_numeric(summary, "delta0_meV")
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError("PRE metadata do not provide a positive delta0_meV.")
    return float(value)


HBAR_J_S = 1.054571817e-34
K_B_J_K = 1.380649e-23


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
        _initial_snapshot_field_temperature(snapshots, "Te_snapshot_K"),
    )
    if not all(np.isfinite(value) and value > 0.0 for value in (diffusion, Tc_K, bias_K)):
        raise ValueError("Could not reconstruct xi; supply --xi-nm or complete PRE metadata.")
    xi2_m2 = (
        math.pi * HBAR_J_S * diffusion
        / (4.0 * math.sqrt(2.0) * K_B_J_K * Tc_K * math.sqrt(1.0 + bias_K / Tc_K))
    )
    return math.sqrt(xi2_m2), "PRE effective gTDGL D and project bias temperature"


def _snapshot_diagnostics(
    *,
    label: str,
    cfg: Mapping[str, Any],
    mesh: Any,
    ops: Any,
    snapshots: Mapping[str, Any],
    snapshots_path: Path,
    summary: Mapping[str, Any],
    requested_times_ps: list[float] | tuple[float, ...],
    power_table_path: Path,
    existing_manifest: Path,
) -> dict[str, Any]:
    selected = nearest_unique_snapshot_indices(
        _snapshot_times_ps(snapshots), requested_times_ps
    )
    transient_config = dict(summary.get("config", {}))
    thermal_enabled = bool(transient_config.get("thermal_enabled", False))
    thermal_window_m = float(transient_config.get("thermal_window_m", 100.0e-9))
    nodes = np.asarray(mesh.nodes, dtype=float)
    thermal_mask = (
        build_central_thermal_mask(nodes, window_m=thermal_window_m)
        if thermal_enabled
        else np.zeros(nodes.shape[0], dtype=bool)
    )
    diagnostics = compute_photon_snapshot_plot_diagnostics(
        snapshots=snapshots,
        selected_indices=selected,
        power_table_npz=str(power_table_path),
        ops=ops,
        sigma_n_S_m=float(cfg["material"]["sigma_n_S_m"]),
        thermal_active_mask=thermal_mask,
        thermal_bath_K=_initial_snapshot_field_temperature(
            snapshots, "Tph_snapshot_K"
        ),
        global_limits_override=_reusable_global_limits(
            existing_manifest,
            snapshots_path=snapshots_path,
            power_table_path=power_table_path,
        ),
        progress_callback=_snapshot_diagnostic_progress(label),
    )
    diagnostics.update(
        {
            "nodes_x_nm": 1.0e9 * nodes[:, 0],
            "nodes_y_nm": 1.0e9 * nodes[:, 1],
            "triangles": np.asarray(mesh.triangles, dtype=np.int64),
        }
    )
    return diagnostics


def _reusable_global_limits(
    manifest_path: Path,
    *,
    snapshots_path: Path,
    power_table_path: Path,
) -> dict[str, np.ndarray] | None:
    """Reuse all-snapshot limits; selected plotting times do not affect them."""
    if not manifest_path.exists():
        return None
    try:
        manifest = _read_yaml(manifest_path)
        if int(manifest.get("schema_version", 0)) < 2:
            return None
        fingerprints = dict(manifest.get("source_fingerprints", {}))
        if fingerprints.get("transient_snapshots") != _file_fingerprint(snapshots_path):
            return None
        if fingerprints.get("power_table_catalog") != _file_fingerprint(power_table_path):
            return None
        raw_limits = dict(manifest.get("snapshot_global_limits", {}))
        required = {
            "q_abs_snapshot_m_inv", "joule_snapshot_W_m3", "P_S_snapshot_W_m3",
            "P_R_snapshot_W_m3", "P_total_snapshot_W_m3", "P_diff_snapshot_W_m3",
            "P_esc_snapshot_W_m3", "u_e_snapshot_J_m3", "u_ph_snapshot_J_m3",
            "C_e_snapshot_J_m3_K", "C_ph_snapshot_J_m3_K", "kappa_s_snapshot_W_m_K",
        }
        if not required.issubset(raw_limits):
            return None
        limits = {key: np.asarray(raw_limits[key], dtype=float).reshape(-1) for key in required}
        if not all(values.size == 2 and np.all(np.isfinite(values)) for values in limits.values()):
            return None
        return limits
    except (OSError, TypeError, ValueError):
        return None


def _snapshot_diagnostic_progress(label: str):
    state = {"bucket": -1}

    def report(completed: int, total: int) -> None:
        if total <= 0:
            return
        bucket = min(10, int(10 * completed / total))
        if bucket > state["bucket"]:
            state["bucket"] = bucket
            print(f" {label} snapshot thermodynamics: {10 * bucket}% ({completed}/{total})")

    return report


def _snapshot_times_ps(snapshots: Mapping[str, Any]) -> np.ndarray:
    if "snapshot_t_ps" in snapshots:
        return np.asarray(snapshots["snapshot_t_ps"], dtype=float).reshape(-1)
    return np.asarray(snapshots.get("snapshot_t_s", []), dtype=float).reshape(-1) / 1.0e-12


def _initial_snapshot_field_temperature(snapshots: Mapping[str, Any], key: str) -> float:
    values = np.asarray(snapshots.get(key, []), dtype=float)
    if values.ndim != 2 or values.shape[0] == 0:
        raise ValueError(f"Photon snapshots lack initial {key} values.")
    finite = values[0][np.isfinite(values[0]) & (values[0] > 0.0)]
    if finite.size == 0:
        raise ValueError(f"Photon snapshots have no positive finite initial {key} values.")
    return float(np.nanmedian(finite))


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


def _file_fingerprint(path: Path) -> dict[str, int]:
    stat = path.stat()
    return {"size_bytes": int(stat.st_size), "mtime_ns": int(stat.st_mtime_ns)}


def _serializable_limits(diagnostics: Mapping[str, Any]) -> dict[str, list[float]]:
    return {
        key: [float(value) for value in np.asarray(bounds, dtype=float).reshape(-1)]
        for key, bounds in dict(diagnostics.get("snapshot_global_limits", {})).items()
    }


def _find_numeric(value: Any, target: str) -> float:
    if isinstance(value, Mapping):
        if target in value:
            try:
                return float(value[target])
            except Exception:
                pass
        for nested in value.values():
            found = _find_numeric(nested, target)
            if np.isfinite(found):
                return found
    return np.nan


def _resolved_times(snapshots: Mapping[str, Any], requested: list[float]) -> list[float]:
    stored = np.asarray(snapshots.get("snapshot_t_ps", []), dtype=float)
    if stored.size == 0:
        return []
    return [float(stored[int(np.nanargmin(np.abs(stored - float(value))))]) for value in requested]


def _photon_time(history: Mapping[str, Any]) -> float | None:
    time = np.asarray(history.get("t_ps", []), dtype=float)
    applied = np.asarray(history.get("photon_applied", []), dtype=bool)
    if time.size == 0 or applied.size == 0:
        return None
    if applied.size != time.size:
        applied = np.resize(applied, time.size)
    indices = np.flatnonzero(applied)
    return float(time[indices[0]]) if indices.size else None


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
