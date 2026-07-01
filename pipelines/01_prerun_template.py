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
from pysnspd.plotting.pre_diagnostics import write_pre_diagnostic_plots
from pysnspd.usadel.supercurrent_table import (
    append_supercurrent_table_3d_to_npz,
    build_matsubara_supercurrent_table_3d,
    supercurrent_table_summary,
    temperature_axis_from_request,
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
        default=None,
        help="Catalogue workers. If omitted, use parallel.workers from the config when enabled.",
    )
    parser.add_argument(
        "--no-progress",
        action="store_true",
        help="Disable the PRE-run stage progress bar.",
    )
    parser.add_argument(
        "--no-diagnostic-plots",
        action="store_true",
        help="Do not write PRE diagnostic plots.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=480,
        help="Resolution for PRE diagnostic plots.",
    )

    # Usadel / OE3.
    parser.add_argument("--eta-fraction", type=float, default=1.0e-3)
    parser.add_argument("--gamma-max-fraction", type=float, default=0.80)
    parser.add_argument("--energy-max-factor", type=float, default=30.0)


    # Strict 3D Usadel supercurrent table for SS.
    parser.add_argument(
        "--js-table-n-Te",
        type=int,
        default=3,
        help="Number of electronic-temperature points for js_A_m2[Te,delta,q]. Use 1 for a cheap smoke PRE-run.",
    )
    parser.add_argument(
        "--js-table-Te-min-K",
        type=float,
        default=None,
        help="Minimum Te for the 3D Usadel supercurrent table. Defaults to bias T.",
    )
    parser.add_argument(
        "--js-table-Te-max-K",
        type=float,
        default=None,
        help="Maximum Te for the 3D Usadel supercurrent table. Defaults near Tc when n_Te > 1.",
    )
    parser.add_argument(
        "--js-table-n-delta",
        type=int,
        default=None,
        help="Number of |Delta| points for the strict current table. Defaults to catalogs.dos.n_delta.",
    )
    parser.add_argument(
        "--js-table-n-q",
        type=int,
        default=None,
        help="Number of q points for the strict current table. Defaults to catalogs.dos.n_q.",
    )

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
    workers = _resolve_workers(cfg, args.workers)
    layout = create_run_layout(cfg, args.run_name)
    run_name = layout["run_name"]
    raw_pre = Path(layout["raw_pre"])
    raw_pre.mkdir(parents=True, exist_ok=True)

    progress = _ProgressBar(total=5, enabled=not args.no_progress)

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

    js_n_delta = int(args.js_table_n_delta or usadel_catalog.delta_values_J.size)
    js_n_q = int(args.js_table_n_q or usadel_catalog.q_values_m_inv.size)
    js_delta_axis = np.linspace(0.0, float(np.max(usadel_catalog.delta_values_J)), js_n_delta)
    js_q_axis = np.linspace(0.0, float(np.max(usadel_catalog.q_values_m_inv)), js_n_q)
    js_Te_axis = temperature_axis_from_request(
        T_bias_K=float(usadel_catalog.metadata["T_bias_K"]),
        Tc_K=float(usadel_catalog.metadata["Tc_K"]),
        n_Te=int(args.js_table_n_Te),
        Te_min_K=args.js_table_Te_min_K,
        Te_max_K=args.js_table_Te_max_K,
    )
    js_table = build_matsubara_supercurrent_table_3d(
        Te_axis_K=js_Te_axis,
        delta_axis_J=js_delta_axis,
        q_axis_m_inv=js_q_axis,
        D_m2_s=float(usadel_catalog.metadata["D_m2_s"]),
        sigma_n_S_m=float(usadel_catalog.metadata["sigma_n_S_m"]),
        n_matsubara=int(usadel_catalog.metadata.get("n_matsubara_configured", 500)),
        workers=workers,
    )
    append_supercurrent_table_3d_to_npz(str(usadel_npz), js_table)
    js_summary = supercurrent_table_summary(js_table)

    usadel_summary = catalog_summary(usadel_catalog)
    usadel_summary["supercurrent_table"] = {
        "stored_in": str(usadel_npz),
        "table_key": "js_A_m2",
        "layout": "Te,delta,q",
        "axis_keys": ["Te_axis_K", "delta_axis_J", "q_axis_m_inv"],
        "source": "Matsubara Usadel local-current table over q, |Delta| and Te.",
        "purpose": "Required by SS usadel_poisson; legacy 1D j_s(q) tables are rejected.",
        **js_summary,
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
    progress.advance("Usadel catalogue and strict 3D Matsubara current table ready")

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

    # ------------------------------------------------------------------
    # PRE diagnostic plots.
    # ------------------------------------------------------------------
    diagnostic_plot_outputs: dict[str, str] = {}
    if args.no_diagnostic_plots:
        progress.advance("PRE diagnostic plots skipped")
    else:
        progress.begin("writing PRE diagnostic plots")
        diagnostic_plot_outputs = write_pre_diagnostic_plots(
            mesh=mesh,
            edge_data=edge_data,
            usadel_catalog=usadel_catalog,
            output_dir=raw_pre / "plots_diagnostics",
            dpi=int(args.dpi),
        )
        outputs.update(diagnostic_plot_outputs)
        progress.advance("PRE diagnostic plots written")

    progress.begin("writing PRE manifest and summary")
    manifest_path = write_manifest(
        cfg,
        run_name,
        stage="pre",
        extra={
            "pipeline": "01_prerun_template.py",
            "purpose": "Official PRE-run: pyTDGL-style mesh, dirty-limit Usadel, Matsubara current table, phase-space catalogue.",
            "workers": int(workers),
            "outputs": outputs,
            "mesh_edge_summary": mesh_edge_summary,
            "usadel_summary": usadel_summary,
            "phase_space_summary": phase_summary_data,
            "diagnostic_plots": diagnostic_plot_outputs,
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


def _resolve_workers(cfg: dict[str, Any], requested: int | None) -> int:
    if requested is not None:
        return max(1, int(requested))
    parallel = cfg.get("parallel", {}) if isinstance(cfg, dict) else {}
    if bool(parallel.get("enabled", False)):
        return max(1, int(parallel.get("workers", 1)))
    return 1


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
