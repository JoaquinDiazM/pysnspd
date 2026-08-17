#!/usr/bin/env python3
"""Diagnose P_Delta/P_q and an energy-variable update from a saved photon run."""

from __future__ import annotations

import argparse
import hashlib
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import yaml

from pysnspd.analysis.energy_projection import extract_energy_projection_diagnostics
from pysnspd.config import load_config, validate_config
from pysnspd.io.manager import create_run_layout
from pysnspd.mesh.delaunay import load_mesh_npz
from pysnspd.mesh.edges import load_edges_npz
from pysnspd.mesh.operators import build_fv_operators
from pysnspd.plotting.energy_projection import write_energy_projection_figures
from pysnspd.plotting.style import THESIS_DPI

HBAR_J_S = 1.054571817e-34
K_B_J_K = 1.380649e-23

DEFAULT_CONFIG = "configs/geminga_local_v3.yaml"
DEFAULT_PRE = "pre_oe6_v4_L360nm_mesh4p0nm_k101T121D141Q_phase200T31D41Q2400W_power200Tph_01"
DEFAULT_RUN = "photon_phasecg_I30uA_0p8eV_sigma10nm_t50ps_1500ps_stiffness_01"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Project an already-saved photon trajectory on u_e(Te,|Delta|,q), "
            "estimate P_spec/P_Delta/P_q, and create three diagnostic PDFs."
        )
    )
    parser.add_argument("--config", default=DEFAULT_CONFIG)
    parser.add_argument("--pre-run-name", default=DEFAULT_PRE)
    parser.add_argument("--run-name", default=DEFAULT_RUN)
    parser.add_argument(
        "--impact-role",
        choices=("center", "edge"),
        required=True,
        help="Role of this persisted run in the matched center/edge audit.",
    )
    parser.add_argument("--window-nm", type=float, default=100.0)
    parser.add_argument("--chunk-rows", type=int, default=32)
    parser.add_argument(
        "--snapshot-stride",
        type=int,
        default=1,
        help="Use every Nth persisted state. One is the faithful full-resolution diagnostic.",
    )
    parser.add_argument(
        "--max-snapshots",
        type=int,
        default=None,
        help=(
            "Stop after this many persisted states for a smoke run. The resulting PDFs "
            "are explicitly labelled truncated and use a separate output directory."
        ),
    )
    parser.add_argument("--figures-subdir", default=None)
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Optional explicit output directory for the compact cache, manifest, and PDFs.",
    )
    parser.add_argument("--force", action="store_true", help="Ignore a matching compact cache.")
    parser.add_argument("--dpi", type=int, default=THESIS_DPI)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.chunk_rows < 1 or args.snapshot_stride < 1:
        raise ValueError("--chunk-rows and --snapshot-stride must be positive.")
    if args.max_snapshots is not None and args.max_snapshots < 2:
        raise ValueError("--max-snapshots must be at least two.")
    if not np.isfinite(args.window_nm) or args.window_nm <= 0.0:
        raise ValueError("--window-nm must be finite and positive.")

    cfg = validate_config(load_config(args.config))
    pre_layout = create_run_layout(cfg, args.pre_run_name)
    run_layout = create_run_layout(cfg, args.run_name)
    raw_pre = Path(pre_layout["raw_pre"])
    raw_photon = Path(run_layout["raw_photon"])
    if args.figures_subdir:
        subdir = str(args.figures_subdir)
    elif args.max_snapshots is None:
        subdir = f"D3_energy_projection_{args.impact_role}"
    else:
        subdir = (
            f"D3_energy_projection_{args.impact_role}_smoke_"
            f"{int(args.max_snapshots)}"
        )
    output_dir = (
        Path(args.output_dir).expanduser().resolve()
        if args.output_dir
        else Path(run_layout["plots_figures"]) / subdir
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    snapshots_path = _require(raw_photon / "transient_snapshots.npz")
    history_path = _require(raw_photon / "transient_history.npz")
    power_table_path = _require(raw_pre / "power_table_catalog.npz")
    usadel_current_path = _require(raw_pre / "usadel_dos_catalog.npz")
    mesh_path = _require(raw_pre / "mesh.npz")
    edges_path = _require(raw_pre / "edges.npz")
    pre_summary_path = _require(raw_pre / "usadel_dos_summary.yaml")
    mesh = load_mesh_npz(mesh_path)
    edges = load_edges_npz(edges_path)
    ops = build_fv_operators(mesh, edges)
    history = _load_npz(history_path)
    pre_summary = _read_yaml(pre_summary_path)
    xi_m, xi_source = _resolve_xi_m(cfg=cfg, pre_summary=pre_summary)

    parameters = {
        "diagnostic_schema": 4,
        "impact_role": str(args.impact_role),
        "window_nm": float(args.window_nm),
        "chunk_rows": int(args.chunk_rows),
        "snapshot_stride": int(args.snapshot_stride),
        "max_snapshots": None if args.max_snapshots is None else int(args.max_snapshots),
    }
    sources = {
        "transient_snapshots": _fingerprint(snapshots_path),
        "transient_history": _fingerprint(history_path),
        "power_table_catalog": _fingerprint(power_table_path),
        "usadel_current_catalog": _fingerprint(usadel_current_path),
        "mesh": _fingerprint(mesh_path),
        "edges": _fingerprint(edges_path),
    }
    cache_path = output_dir / "D3_energy_projection_cache.npz"
    manifest_path = output_dir / "D3_energy_projection_manifest.yaml"
    cached_manifest = _read_yaml(manifest_path) if manifest_path.exists() else {}
    use_cache = (
        not bool(args.force)
        and cache_path.exists()
        and cached_manifest.get("parameters") == parameters
        and cached_manifest.get("source_fingerprints") == sources
    )
    if use_cache:
        print(f"Reusing compact D3 cache: {cache_path}")
        diagnostics = _load_npz(cache_path)
    else:
        print("D3 streaming extraction (no solver step will be executed)")
        print(f" snapshots: {snapshots_path}")
        print(f" central window: {args.window_nm:g} nm")
        diagnostics = extract_energy_projection_diagnostics(
            snapshots_npz=snapshots_path,
            power_table_npz=power_table_path,
            usadel_current_npz=usadel_current_path,
            history=history,
            nodes_m=np.asarray(mesh.nodes, dtype=float),
            triangles=np.asarray(mesh.triangles, dtype=np.int64),
            ops=ops,
            sigma_n_S_m=float(cfg["material"]["sigma_n_S_m"]),
            thickness_m=float(cfg["material"]["thickness_m"]),
            Tc_K=float(cfg["material"]["Tc_K"]),
            xi_m=float(xi_m),
            window_m=float(args.window_nm) * 1.0e-9,
            chunk_rows=int(args.chunk_rows),
            snapshot_stride=int(args.snapshot_stride),
            max_snapshots=args.max_snapshots,
            progress_callback=_progress_bar(),
        )
        np.savez_compressed(cache_path, **diagnostics)
        print(f" compact cache: {cache_path}")

    saved = write_energy_projection_figures(
        diagnostics,
        output_dir,
        dpi=int(args.dpi),
    )
    manifest = {
        "schema_version": 2,
        "pipeline": "plot_pipelines/D3_photon_energy_projection.py",
        "purpose": (
            "A-posteriori Simon-energy projection of a persisted photon run; "
            "no resimulation and no modification of the solver state."
        ),
        "run_name": str(args.run_name),
        "impact_role": str(args.impact_role),
        "pre_run_name": str(args.pre_run_name),
        "raw_photon": str(raw_photon),
        "raw_pre": str(raw_pre),
        "output_dir": str(output_dir),
        "cache": str(cache_path),
        "cache_reused": bool(use_cache),
        "parameters": parameters,
        "source_fingerprints": sources,
        "source_content_hashes": {
            "usadel_current_catalog_sha256": _sha256(usadel_current_path),
        },
        "xi_m": float(xi_m),
        "xi_source": str(xi_source),
        "selected_times_ps": [
            float(value) for value in np.asarray(diagnostics["selected_times_ps"]).reshape(-1)
        ],
        "dropped_duplicate_time_count": int(
            np.asarray(diagnostics["dropped_duplicate_time_count"]).reshape(-1)[0]
        ),
        "truncated": bool(np.asarray(diagnostics["truncated"]).reshape(-1)[0]),
        "sign_convention": {
            "P_ep": "positive means energy leaves electrons",
            "P_spec": "positive means energy is stored by the moving spectrum at fixed Te",
            "temperature_equation": "Ce*dTe/dt = Pdiff + PJ - Pep - Pspec",
            "P_q_scope": (
                "spectral dependence of the Simon u_e catalogue only; it is not asserted "
                "to be the complete superflow/circuit energy"
            ),
        },
        "finite_step_policy": {
            "P_delta": "change Delta first at fixed (Te,q_n)",
            "P_q": "then change q at fixed (Te,Delta_n+1)",
            "P_path": "difference in P_delta under the opposite q-then-Delta path",
        },
        "strict_q_supercurrent_audit": _strict_q_supercurrent_summary(diagnostics),
        "figures": {key: str(path) for key, path in saved.items()},
    }
    with manifest_path.open("w", encoding="utf-8") as stream:
        yaml.safe_dump(manifest, stream, sort_keys=False, allow_unicode=True)

    print("D3 energy projection diagnostics")
    print(f" run:             {args.run_name}")
    print(f" impact role:     {args.impact_role}")
    print(f" selected times:  {', '.join(f'{v:.6g}' for v in manifest['selected_times_ps'])} ps")
    print(f" xi:              {1.0e9 * xi_m:.6g} nm ({xi_source})")
    print(f" truncated:       {manifest['truncated']}")
    print(f" duplicate times: {manifest['dropped_duplicate_time_count']} later state kept")
    for key, path in saved.items():
        print(f" {key:16s} {path}")
    print(f" manifest:        {manifest_path}")
    if manifest["truncated"]:
        print("Full run: rerun without --max-snapshots; the compact cache makes later plotting cheap.")
    print("Status: OK")
    return 0


def _progress_bar() -> Any:
    state = {"bucket": -1}

    def report(done: int, total: int, elapsed_s: float) -> None:
        fraction = min(max(float(done) / max(float(total), 1.0), 0.0), 1.0)
        bucket = int(20 * fraction)
        if bucket <= state["bucket"] and done < total:
            return
        state["bucket"] = bucket
        width = 30
        filled = int(round(width * fraction))
        rate = float(done) / max(float(elapsed_s), 1.0e-9)
        print(
            f" [{('#' * filled).ljust(width, '-')}] {100.0 * fraction:5.1f}% "
            f"({done}/{total}; {rate:.1f} snapshots/s; {elapsed_s:.1f}s)"
        )

    return report


def _resolve_xi_m(
    *,
    cfg: Mapping[str, Any],
    pre_summary: Mapping[str, Any],
) -> tuple[float, str]:
    diffusion = _first_positive(
        _nested(pre_summary, ("usadel", "gtdgl_allmaras", "D_effective_m2_s")),
        _nested(pre_summary, ("metadata", "gtdgl_allmaras_D_m2_s")),
        _nested(pre_summary, ("usadel", "D_m2_s")),
        _nested(cfg, ("material", "D_m2_s")),
    )
    Tc_K = _first_positive(
        _nested(pre_summary, ("metadata", "Tc_K")),
        _nested(cfg, ("material", "Tc_K")),
    )
    bias_K = _first_positive(_nested(cfg, ("bias", "T_bias_K")))
    if not all(np.isfinite(value) and value > 0.0 for value in (diffusion, Tc_K, bias_K)):
        raise ValueError("Could not reconstruct xi from PRE/config metadata.")
    xi2 = math.pi * HBAR_J_S * diffusion / (
        4.0 * math.sqrt(2.0) * K_B_J_K * Tc_K * math.sqrt(1.0 + bias_K / Tc_K)
    )
    return math.sqrt(xi2), "PRE effective gTDGL D, Tc and project bias temperature"


def _strict_q_supercurrent_summary(
    diagnostics: Mapping[str, Any],
) -> dict[str, Any]:
    fraction = np.asarray(
        diagnostics["strict_q_clipped_fraction"], dtype=float
    ).reshape(-1)
    dt_ps = np.asarray(diagnostics["dt_ps"], dtype=float).reshape(-1)
    p95_ratio = np.asarray(
        diagnostics["strict_q_clipped_js_p95_over_catalog_max"], dtype=float
    ).reshape(-1)
    near_zero_fraction = np.asarray(
        diagnostics["strict_q_clipped_js_near_zero_fraction"], dtype=float
    ).reshape(-1)
    active = np.isfinite(fraction) & (fraction > 0.0)
    finite_ratio = p95_ratio[active & np.isfinite(p95_ratio)]
    finite_near_zero = near_zero_fraction[
        active & np.isfinite(near_zero_fraction)
    ]
    max_p95_ratio = (
        float(np.nanmax(finite_ratio)) if finite_ratio.size else None
    )
    minimum_near_zero_fraction = (
        float(np.nanmin(finite_near_zero)) if finite_near_zero.size else None
    )
    supports_near_zero = bool(
        max_p95_ratio is not None
        and minimum_near_zero_fraction is not None
        and max_p95_ratio <= 1.0e-3
        and minimum_near_zero_fraction >= 0.95
    )
    return {
        "q_axis_max_m_inv": float(
            np.asarray(diagnostics["strict_q_max_m_inv"], dtype=float).reshape(-1)[0]
        ),
        "catalog_js_reference_A_m2": float(
            np.asarray(
                diagnostics["strict_js_reference_A_m2"], dtype=float
            ).reshape(-1)[0]
        ),
        "maximum_clipped_node_fraction": (
            float(np.nanmax(fraction)) if fraction.size else 0.0
        ),
        "sampled_duration_with_q_clipping_ps": float(
            np.sum(dt_ps[: active.size][active[: dt_ps.size]])
        ),
        "clipping_active_sample_count": int(np.sum(active)),
        "median_temporal_p95_js_over_catalog_max": (
            float(np.nanmedian(finite_ratio)) if finite_ratio.size else None
        ),
        "temporal_p95_of_p95_js_over_catalog_max": (
            float(np.nanpercentile(finite_ratio, 95.0))
            if finite_ratio.size
            else None
        ),
        "maximum_temporal_p95_js_over_catalog_max": max_p95_ratio,
        "median_temporal_fraction_below_1e-3_catalog_max": (
            float(np.nanmedian(finite_near_zero))
            if finite_near_zero.size
            else None
        ),
        "fraction_of_clipping_active_samples_with_95pct_near_zero": (
            float(np.mean(finite_near_zero >= 0.95))
            if finite_near_zero.size
            else None
        ),
        "minimum_temporal_fraction_below_1e-3_catalog_max": minimum_near_zero_fraction,
        "near_zero_definition": "abs(js) <= 1e-3 * max(abs(js_A_m2 catalogue))",
        "supports_near_zero_claim_under_declared_definition": supports_near_zero,
        "lookup_policy": (
            "Te, Delta, and abs(q) are clipped to the strict PRE table axes before "
            "interpolation; the statistic reports the current actually returned at "
            "the table boundary, not an extrapolated current."
        ),
    }


def _first_positive(*values: Any) -> float:
    for value in values:
        try:
            number = float(value)
        except (TypeError, ValueError):
            continue
        if np.isfinite(number) and number > 0.0:
            return number
    return np.nan


def _nested(mapping: Mapping[str, Any], path: tuple[str, ...]) -> Any:
    value: Any = mapping
    for key in path:
        if not isinstance(value, Mapping) or key not in value:
            return None
        value = value[key]
    return value


def _fingerprint(path: Path) -> dict[str, int]:
    stat = path.stat()
    return {"size_bytes": int(stat.st_size), "mtime_ns": int(stat.st_mtime_ns)}


def _sha256(path: Path, *, chunk_bytes: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(chunk_bytes):
            digest.update(chunk)
    return digest.hexdigest()


def _read_yaml(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as stream:
        value = yaml.safe_load(stream) or {}
    return value if isinstance(value, dict) else {}


def _load_npz(path: Path) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as data:
        return {key: np.asarray(data[key]) for key in data.files}


def _require(path: Path) -> Path:
    if not path.exists():
        raise FileNotFoundError(path)
    return path


if __name__ == "__main__":
    raise SystemExit(main())
