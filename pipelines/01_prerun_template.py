"""Official PRE-run pipeline for pySNSPD.

This stage builds the reusable, expensive objects needed by later stationary and
photon runs:

1. the 2D nanowire mesh and boundary edge table;
2. the dirty-limit Usadel/DOS catalogue;
3. a Matsubara Usadel supercurrent-density table saved into the same NPZ;
4. the superconducting phase-space catalogue used by the kinetic layer.

The PRE-run is the only pipeline that should spend time building catalogues.
Later SS/PHOTON stages load these objects instead of recomputing them.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import numpy as np
import yaml

from pysnspd.config import load_config, validate_config
from pysnspd.io.manager import create_run_layout, write_manifest
from pysnspd.mesh.delaunay import (
    generate_rectangular_delaunay_mesh,
    mesh_summary,
    save_mesh_npz,
)
from pysnspd.mesh.edges import (
    assert_edge_data_consistent,
    build_edge_data,
    edge_summary,
    save_edges_npz,
)
from pysnspd.usadel.catalog import (
    build_usadel_catalog_from_config,
    catalog_summary,
    save_usadel_catalog_npz,
)
from pysnspd.kinetic.phase_space import (
    build_phase_space_catalog_from_usadel_catalog,
    phase_space_summary,
    save_phase_space_catalog_npz,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate PRE-run mesh, Usadel and phase-space catalogues."
    )
    parser.add_argument("--config", required=True, help="Path to YAML project config.")
    parser.add_argument(
        "--run-name",
        default=None,
        help="Optional run name. If omitted, project.default_run_name is used.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Reserved for catalogue-building workflows; recorded in the manifest.",
    )
    parser.add_argument(
        "--no-progress",
        action="store_true",
        help="Disable the PRE-run stage progress bar.",
    )

    # Usadel / OE3.
    parser.add_argument("--eta-fraction", type=float, default=1.0e-3)
    parser.add_argument("--gamma-max-fraction", type=float, default=0.80)
    parser.add_argument("--energy-max-factor", type=float, default=30.0)

    # Phase space / OE4.
    parser.add_argument("--skip-phase-space", action="store_true")
    parser.add_argument("--phase-n-Te", type=int, default=6)
    parser.add_argument("--phase-n-delta", type=int, default=6)
    parser.add_argument("--phase-n-q", type=int, default=6)
    parser.add_argument("--phase-n-omega", type=int, default=480)
    parser.add_argument("--phase-omega-max-meV", type=float, default=35.0)
    parser.add_argument("--phase-Te-min-K", type=float, default=None)
    parser.add_argument("--phase-Te-max-K", type=float, default=None)

    return parser.parse_args()


def main() -> int:
    args = parse_args()
    cfg = validate_config(load_config(args.config))
    layout = create_run_layout(cfg, args.run_name)
    run_name = layout["run_name"]
    raw_pre = Path(layout["raw_pre"])
    raw_pre.mkdir(parents=True, exist_ok=True)

    progress = _ProgressBar(total=4, enabled=not args.no_progress)

    # ------------------------------------------------------------------
    # OE2: pyTDGL-style mesh and boundary edges.
    # ------------------------------------------------------------------
    progress.begin("building pyTDGL-style mesh and edge table")
    mesh = generate_rectangular_delaunay_mesh(cfg)
    edge_data = build_edge_data(
        mesh.nodes,
        mesh.triangles,
        length_m=mesh.length_m,
        width_m=mesh.width_m,
    )
    assert_edge_data_consistent(edge_data)

    mesh_npz = save_mesh_npz(mesh, raw_pre / "mesh.npz")
    edges_npz = save_edges_npz(edge_data, raw_pre / "edges.npz")

    mesh_edge_summary = {
        "run_name": run_name,
        "mesh": mesh_summary(mesh),
        "edges": edge_summary(edge_data),
    }
    mesh_summary_path = raw_pre / "mesh_summary.yaml"
    _write_yaml(mesh_summary_path, mesh_edge_summary)
    progress.advance("mesh and edge table ready")

    # ------------------------------------------------------------------
    # OE3: dirty-limit Usadel/DOS catalogue.
    # ------------------------------------------------------------------
    progress.begin("building dirty-limit Usadel catalogue")
    usadel_catalog = build_usadel_catalog_from_config(
        cfg,
        eta_fraction=args.eta_fraction,
        gamma_max_fraction=args.gamma_max_fraction,
        energy_max_factor=args.energy_max_factor,
    )
    usadel_npz = save_usadel_catalog_npz(usadel_catalog, raw_pre / "usadel_dos_catalog.npz")
    _append_matsubara_supercurrent_table(usadel_npz, usadel_catalog)

    usadel_summary = catalog_summary(usadel_catalog)
    usadel_summary["supercurrent_table"] = {
        "stored_in": str(usadel_npz),
        "table_key": "j_s_A_m2",
        "axis_key": "q_axis_m_inv",
        "source": "Matsubara Usadel calibration sweep at T_bias.",
        "purpose": "Used by the flat gTDGL SS backend when supercurrent_law=usadel_poisson.",
    }
    usadel_summary_path = raw_pre / "usadel_dos_summary.yaml"
    _write_yaml(
        usadel_summary_path,
        {
            "run_name": run_name,
            "usadel": usadel_summary,
            "metadata": usadel_catalog.metadata,
        },
    )
    progress.advance("Usadel catalogue and Matsubara current table ready")

    outputs: dict[str, str] = {
        "mesh_npz": str(mesh_npz),
        "edges_npz": str(edges_npz),
        "mesh_summary": str(mesh_summary_path),
        "usadel_npz": str(usadel_npz),
        "usadel_summary": str(usadel_summary_path),
    }
    phase_summary_data: dict[str, Any] | None = None

    # ------------------------------------------------------------------
    # OE4: phase-space catalogue.
    # ------------------------------------------------------------------
    if args.skip_phase_space:
        progress.advance("phase-space catalogue skipped")
    else:
        progress.begin("building superconducting phase-space catalogue")
        phase_catalog = build_phase_space_catalog_from_usadel_catalog(
            usadel_catalog,
            cfg,
            n_Te=args.phase_n_Te,
            n_delta=args.phase_n_delta,
            n_q=args.phase_n_q,
            n_omega=args.phase_n_omega,
            Te_min_K=args.phase_Te_min_K,
            Te_max_K=args.phase_Te_max_K,
            omega_max_meV=args.phase_omega_max_meV,
        )
        phase_npz = save_phase_space_catalog_npz(
            phase_catalog,
            raw_pre / "phase_space_catalog.npz",
        )
        phase_summary_data = phase_space_summary(phase_catalog)
        phase_summary_path = raw_pre / "phase_space_summary.yaml"
        _write_yaml(
            phase_summary_path,
            {
                "run_name": run_name,
                "phase_space": phase_summary_data,
                "metadata": phase_catalog.metadata,
            },
        )
        outputs.update(
            {
                "phase_space_npz": str(phase_npz),
                "phase_space_summary": str(phase_summary_path),
            }
        )
        progress.advance("phase-space catalogue ready")

    progress.begin("writing PRE manifest and summary")
    manifest_path = write_manifest(
        cfg,
        run_name,
        stage="pre",
        extra={
            "pipeline": "01_prerun_template.py",
            "purpose": "Official PRE-run: pyTDGL-style mesh, dirty-limit Usadel, Matsubara current table, phase-space catalogue.",
            "workers": int(args.workers),
            "outputs": outputs,
            "mesh_edge_summary": mesh_edge_summary,
            "usadel_summary": usadel_summary,
            "phase_space_summary": phase_summary_data,
        },
    )
    progress.advance("PRE manifest written")

    print()
    print("PRE-run generation")
    print(f"  run_name: {run_name}")
    print(f"  raw_pre:  {raw_pre}")
    print()
    print("Mesh")
    _print_dict(mesh_edge_summary["mesh"])
    print()
    print("Edges")
    _print_dict(mesh_edge_summary["edges"])
    print()
    print("Usadel")
    _print_dict(usadel_summary)
    print()
    if phase_summary_data is None:
        print("Phase-space: skipped")
    else:
        print("Phase-space")
        _print_dict(phase_summary_data)
    print()
    print("Outputs")
    for key, value in outputs.items():
        print(f"  {key}: {value}")
    print(f"  pre_manifest: {manifest_path}")
    print("Status: OK")

    return 0


def _append_matsubara_supercurrent_table(npz_path: str | Path, usadel_catalog: Any) -> None:
    """Add the Matsubara Usadel supercurrent table to the PRE NPZ.

    The gTDGL Appendix-B current closure is obtained from a precalculated table.
    The current table saved here is the stable calibration sweep

        j_s(q, T_bias) = (2 pi k_B T_bias / |e| hbar) sigma_n q sum_n s_n^2,

    already computed by the Usadel layer. The flat gTDGL backend can use this
    one-dimensional table directly for stationary/frozen-temperature tests.
    """

    path = Path(npz_path)
    with np.load(path, allow_pickle=True) as data:
        arrays = {key: data[key] for key in data.files}

    arrays["j_s_A_m2"] = np.asarray(
        usadel_catalog.calibration_current_density_values_A_m2,
        dtype=float,
    )
    arrays["js_A_m2"] = arrays["j_s_A_m2"]
    arrays["q_axis_m_inv"] = np.asarray(
        usadel_catalog.calibration_q_values_m_inv,
        dtype=float,
    )
    arrays["js_table_T_K"] = np.array(float(usadel_catalog.metadata["T_bias_K"]))
    arrays["js_table_n_matsubara"] = np.array(
        int(usadel_catalog.metadata.get("n_matsubara_configured", -1)),
        dtype=np.int64,
    )

    np.savez_compressed(path, **arrays)


def _write_yaml(path: str | Path, data: MappingLike) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as f:
        yaml.safe_dump(data, f, sort_keys=False, allow_unicode=True)


def _print_dict(data: dict[str, Any]) -> None:
    for key, value in data.items():
        print(f"  {key}: {value}")


class _ProgressBar:
    """Small dependency-free stage progress bar for PRE smoke tests."""

    def __init__(self, *, total: int, enabled: bool = True, width: int = 28) -> None:
        self.total = int(total)
        self.enabled = bool(enabled)
        self.width = int(width)
        self.current = 0

    def begin(self, label: str) -> None:
        if not self.enabled:
            return
        print(f"PRE-run: {label} ...", flush=True)

    def advance(self, label: str) -> None:
        if not self.enabled:
            return
        self.current = min(self.current + 1, self.total)
        frac = self.current / self.total if self.total > 0 else 1.0
        filled = int(round(self.width * frac))
        bar = "#" * filled + "-" * (self.width - filled)
        percent = int(round(100.0 * frac))
        print(f"PRE-run [{bar}] {percent:3d}%  {label}", flush=True)


MappingLike = dict[str, Any]


if __name__ == "__main__":
    raise SystemExit(main())
