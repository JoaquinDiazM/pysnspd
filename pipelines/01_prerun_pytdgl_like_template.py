#!/usr/bin/env python3
"""pyTDGL-faithful PRE-run entry point in SI units.

This wrapper delegates the heavy PRE computation to the standard pySNSPD
``01_prerun_template.py`` after removing the extra pyTDGL-only CLI flags from
``sys.argv``.  After the standard PRE-run finishes, it writes the full
pyTDGL-like finite-volume sidecar: Mesh, EdgeMesh and Voronoi arrays.

The meshing order follows pyTDGL's ``Device.make_mesh`` pattern:

    generate_mesh(...) -> Mesh.from_triangulation(..., create_submesh=False)
    -> Mesh.smooth(..., create_submesh=False)
    -> Mesh.from_triangulation(..., create_submesh=True)

Coordinates remain in SI meters; no pyTDGL coherence-length normalization is
applied here.
"""
from __future__ import annotations

import argparse
import os
import runpy
import sys
import tempfile
from pathlib import Path
from typing import Any, Mapping

import yaml

from pysnspd.config import load_config
from pysnspd.mesh.pytdgl_like import (
    build_pytdgl_like_mesh_summary,
    generate_rectangular_pytdgl_fvm_mesh_from_parameters,
    parameters_from_config,
    save_pytdgl_like_mesh_npz,
)
from pysnspd.plotting.pytdgl_mesh import plot_pytdgl_fvm_mesh


def _parse_wrapper_args(argv: list[str]) -> tuple[argparse.Namespace, list[str]]:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--config", required=True)
    parser.add_argument("--run-name", required=True)
    parser.add_argument("--jitter-fraction", type=float, default=0.0)
    parser.add_argument("--boundary-guard-layers", type=int, default=1)
    parser.add_argument("--pytdgl-min-angle", type=float, default=32.5)
    parser.add_argument("--pytdgl-smooth", type=int, default=None)
    parser.add_argument("--pytdgl-max-edge-length-m", type=float, default=None)
    parser.add_argument("--pytdgl-min-points", type=int, default=None)
    return parser.parse_known_args(argv)


def _config_section(config: Mapping[str, Any], name: str) -> Mapping[str, Any]:
    section = config.get(name, {})
    if not isinstance(section, Mapping):
        raise TypeError(f"config[{name!r}] must be a mapping.")
    return section


def _load_raw_yaml_config(path: str | Path) -> dict[str, Any]:
    data = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    if data is None:
        return {}
    if not isinstance(data, dict):
        raise TypeError(f"Expected a YAML mapping in {path!s}.")
    return data


def _write_mesh_override_config(config_path: str | Path, args: argparse.Namespace) -> Path:
    """Write a temporary config carrying pyTDGL mesh controls.

    The standard PRE pipeline only knows the old pySNSPD CLI.  Passing a
    temporary config keeps its mesh generation and this wrapper's sidecar on the
    same pyTDGL parameters without changing the standard pipeline parser.
    """
    cfg = _load_raw_yaml_config(config_path)
    mesh_cfg = dict(cfg.get("mesh", {}) or {})
    if args.pytdgl_max_edge_length_m is not None:
        mesh_cfg["pytdgl_max_edge_length_m"] = float(args.pytdgl_max_edge_length_m)
    else:
        mesh_cfg.pop("pytdgl_max_edge_length_m", None)
    if args.pytdgl_min_angle is not None:
        mesh_cfg["pytdgl_min_angle_deg"] = float(args.pytdgl_min_angle)
    if args.pytdgl_smooth is not None:
        mesh_cfg["pytdgl_smooth"] = int(args.pytdgl_smooth)
    if args.pytdgl_min_points is not None:
        mesh_cfg["pytdgl_min_points"] = int(args.pytdgl_min_points)
    cfg["mesh"] = mesh_cfg

    tmp = tempfile.NamedTemporaryFile(
        mode="w",
        suffix=".yaml",
        prefix="pysnspd_pytdgl_mesh_",
        delete=False,
        encoding="utf-8",
    )
    with tmp:
        yaml.safe_dump(cfg, tmp, sort_keys=False)
    return Path(tmp.name)


def _iter_path_values(config: Mapping[str, Any]) -> list[Path]:
    """Collect plausible data-root paths from an arbitrary config mapping."""
    out: list[Path] = []

    def visit(obj: Any, key: str = "") -> None:
        if isinstance(obj, Mapping):
            for k, v in obj.items():
                visit(v, str(k))
            return
        if not isinstance(obj, (str, os.PathLike)):
            return
        key_l = key.lower()
        value = str(obj)
        if (
            "big_data" in key_l
            or "data_root" in key_l
            or "root" == key_l
            or "big_data" in value
        ):
            out.append(Path(value).expanduser())

    visit(config)
    return out


def _run_dirs_from_config(config: Mapping[str, Any], run_name: str) -> tuple[Path, Path]:
    """Resolve PRE output directories robustly after the delegated PRE stage.

    The standard PRE pipeline has already written ``raw/<run_name>/pre`` before
    this wrapper writes the pyTDGL-like sidecar.  Therefore this resolver first
    honors explicit config roots, then falls back to the same practical roots
    used throughout the local pySNSPD workflow.  This avoids failing after a
    successful PRE-run merely because the temporary wrapper config does not
    expose ``paths.big_data_root`` under that exact key.
    """
    candidates: list[Path] = []
    paths = config.get("paths", {})
    if isinstance(paths, Mapping):
        for key in (
            "big_data_root",
            "big_data_root_path",
            "data_root",
            "root",
            "base_dir",
        ):
            value = paths.get(key)
            if value is not None:
                candidates.append(Path(str(value)).expanduser())

    candidates.extend(_iter_path_values(config))

    for env_key in (
        "PYSNSPD_BIG_DATA_ROOT",
        "PYSNSPD_DATA_ROOT",
        "BIG_DATA_ROOT",
        "SCRATCH_BIG_DATA_ROOT",
    ):
        value = os.environ.get(env_key)
        if value:
            candidates.append(Path(value).expanduser())

    candidates.extend(
        [
            Path.home() / "scratch" / "big_data",
            Path("/home/jdiaz/scratch/big_data"),
        ]
    )

    unique: list[Path] = []
    seen: set[str] = set()
    for root in candidates:
        root = root.expanduser()
        marker = str(root)
        if marker not in seen:
            unique.append(root)
            seen.add(marker)

    for root in unique:
        raw_pre = root / "raw" / run_name / "pre"
        if raw_pre.exists():
            plots_mesh = root / "plots" / run_name / "mesh"
            plots_mesh.mkdir(parents=True, exist_ok=True)
            return raw_pre, plots_mesh

    if unique:
        root = unique[0]
        raw_pre = root / "raw" / run_name / "pre"
        plots_mesh = root / "plots" / run_name / "mesh"
        raw_pre.mkdir(parents=True, exist_ok=True)
        plots_mesh.mkdir(parents=True, exist_ok=True)
        return raw_pre, plots_mesh

    raise KeyError(
        "Could not resolve big-data root. Expected an existing "
        f"raw/{run_name}/pre directory or a config/env data root."
    )


def _run_standard_prerun(
    *,
    standard_pipeline: Path,
    args: argparse.Namespace,
    passthrough: list[str],
) -> None:
    """Run the standard PRE pipeline without exposing pyTDGL-only flags to it."""
    forwarded_argv = [
        str(standard_pipeline),
        "--config",
        str(args.config),
        "--run-name",
        str(args.run_name),
        "--jitter-fraction",
        str(args.jitter_fraction),
        "--boundary-guard-layers",
        str(args.boundary_guard_layers),
        *passthrough,
    ]
    old_argv = sys.argv[:]
    try:
        sys.argv = forwarded_argv
        try:
            runpy.run_path(str(standard_pipeline), run_name="__main__")
        except SystemExit as exc:
            if exc.code not in (None, 0):
                raise
    finally:
        sys.argv = old_argv


def main() -> int:
    args, passthrough = _parse_wrapper_args(sys.argv[1:])
    standard_pipeline = Path(__file__).with_name("01_prerun_template.py")

    effective_config = _write_mesh_override_config(args.config, args)
    args_for_standard = argparse.Namespace(**vars(args))
    args_for_standard.config = str(effective_config)

    _run_standard_prerun(
        standard_pipeline=standard_pipeline,
        args=args_for_standard,
        passthrough=passthrough,
    )

    cfg = load_config(effective_config)
    raw_pre, plots_mesh = _run_dirs_from_config(cfg, args.run_name)

    params = parameters_from_config(
        cfg,
        jitter_fraction=args.jitter_fraction,
        boundary_guard_layers=args.boundary_guard_layers,
        max_edge_length_m=args.pytdgl_max_edge_length_m,
        min_angle_deg=args.pytdgl_min_angle,
        smooth=args.pytdgl_smooth,
        min_points=args.pytdgl_min_points,
    )
    fvm_mesh = generate_rectangular_pytdgl_fvm_mesh_from_parameters(params)
    sidecar = save_pytdgl_like_mesh_npz(fvm_mesh, raw_pre / "pytdgl_fvm_mesh.npz")
    summary = build_pytdgl_like_mesh_summary(fvm_mesh)
    summary.update(
        {
            "max_edge_length_m": float(params.max_edge_length_m),
            "min_angle_deg": float(params.min_angle_deg),
            "smooth": int(params.smooth),
            "min_points": None if params.min_points is None else int(params.min_points),
            "units_policy": "SI meters; no pyTDGL coherence-length normalization in PRE mesh.",
        }
    )
    summary_path = raw_pre / "pytdgl_fvm_mesh_summary.yaml"
    with summary_path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(summary, f, sort_keys=False)
    plot_path = plot_pytdgl_fvm_mesh(fvm_mesh, plots_mesh / "pytdgl_fvm_mesh.png")

    print("pyTDGL-like finite-volume sidecar")
    print(f"  pytdgl_fvm_mesh_npz    : {sidecar}")
    print(f"  pytdgl_fvm_mesh_summary: {summary_path}")
    print(f"  pytdgl_fvm_mesh_plot   : {plot_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
