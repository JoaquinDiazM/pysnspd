#!/usr/bin/env python3
"""E1 PRE plotting pipeline.

This plotting-only pipeline reads an existing PRE-run Usadel catalogue and writes
E-type thesis figures in PDF format.

Default output reads only completed PRE catalogues and does not reconstruct the
expensive Delta_eq(T) figure:

- E1_usadel_supercurrent_curve.pdf
- E1_usadel_dos_curves_delta_eq.pdf
- E1_usadel_dos_curves_delta0.pdf
- E1_eliashberg_spectral_function_phdos.pdf
- E1_power_channels_Te_Tph_maps.pdf
- E1_power_exchange_vs_temperature.pdf
- E1_energy_heat_capacity_curves.pdf
- E1_electronic_thermal_conductivity_curves.pdf
- mesh_pytdgl_style.pdf

Use --with-gap-plot only when the additional Delta_eq(T) figure is needed.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any

import yaml

# Allow running as ``python plot_pipelines/E1_plot_prerun.py`` from repo root.
_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from pysnspd.config import load_config, validate_config
from pysnspd.io.manager import create_run_layout
from pysnspd.kinetic.eliashberg import load_simon_eliashberg_dat
from pysnspd.mesh.delaunay import load_mesh_npz
from pysnspd.plotting.eliashberg_spectrum import plot_eliashberg_spectrum
from pysnspd.plotting.mesh import plot_mesh_pytdgl_style
from pysnspd.plotting.pre_diagnostics import plot_usadel_supercurrent_curve
from pysnspd.plotting.power_diagnostics import (
    load_power_table_plot_catalog,
    plot_electronic_thermal_conductivity_curves,
    plot_energy_heat_capacity_curves,
    plot_power_channels_Te_Tph_maps,
    plot_power_total_Te_curves,
)
from pysnspd.plotting.usadel_dos_curves import (
    plot_usadel_dos_curves_equilibrium_gap,
    plot_usadel_dos_curves_fixed_delta0,
)
from pysnspd.plotting.usadel_gap import load_usadel_gap_catalog, plot_gap_eq_vs_temperature
from pysnspd.usadel.catalog import load_usadel_catalog_npz


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "E1 PRE plotting pipeline: write the Usadel supercurrent PDF, two DOS-curve PDFs, "
            "the Eliashberg/PhDOS spectrum PDF, and optionally the expensive Delta_eq(T) PDF."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--config",
        required=True,
        type=Path,
        help="pySNSPD YAML config used by the PRE-run.",
    )
    parser.add_argument(
        "--pre-run-name",
        required=True,
        help="Existing PRE-run name produced by 01_prerun_template.py.",
    )
    parser.add_argument(
        "--catalog-npz",
        type=Path,
        default=None,
        help=(
            "Optional direct path to the Usadel .npz catalog; defaults to "
            "raw/<pre-run-name>/pre/usadel_dos_catalog.npz."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Optional output directory; defaults to plots/<pre-run-name>/figures/E1_prerun.",
    )
    parser.add_argument(
        "--eliashberg-dat",
        default=None,
        help=(
            "Path to Simon/MIT NbN alpha2F/PhDOS data. Defaults to "
            "<project.big_data_root>/catalogs/simon_2025/nbn-a2f-ph.dat."
        ),
    )
    parser.add_argument(
        "--eliashberg-pdf-name",
        default="E1_eliashberg_spectral_function_phdos.pdf",
        help="Output PDF filename for the full Eliashberg alpha2F and PhDOS spectrum.",
    )
    parser.add_argument(
        "--skip-eliashberg-spectrum",
        action="store_true",
        help="Skip the Eliashberg alpha2F/PhDOS spectrum PDF.",
    )
    parser.add_argument(
        "--power-table-npz",
        type=Path,
        default=None,
        help=(
            "Optional direct path to the PRE power table; defaults to "
            "raw/<pre-run-name>/pre/power_table_catalog.npz."
        ),
    )
    parser.add_argument(
        "--skip-power-table-figures",
        action="store_true",
        help="Skip the four E1 thermodynamic and projected-power figures.",
    )

    parser.add_argument(
        "--supercurrent-pdf-name",
        default="E1_usadel_supercurrent_curve.pdf",
        help="Output PDF filename for the supercurrent calibration figure.",
    )
    parser.add_argument(
        "--dos-eq-pdf-name",
        default="E1_usadel_dos_curves_delta_eq.pdf",
        help="Output PDF filename for the DOS curves evaluated at Delta_eq(q).",
    )
    parser.add_argument(
        "--dos-delta0-pdf-name",
        default="E1_usadel_dos_curves_delta0.pdf",
        help="Output PDF filename for the DOS curves evaluated at fixed Delta_0.",
    )
    parser.add_argument(
        "--skip-dos-curves",
        action="store_true",
        help="Only write the supercurrent figure; skip the two DOS-curve PDFs.",
    )
    parser.add_argument(
        "--dos-current-fractions",
        nargs="+",
        type=float,
        default=(0.0, 0.25, 0.50, 0.65, 0.80, 0.95),
        help="Current fractions I_s/I_c used for both DOS-curve PDFs.",
    )
    parser.add_argument(
        "--dos-energy-max-meV",
        type=float,
        default=None,
        help="Optional upper energy limit for the DOS-curve PDFs.",
    )
    parser.add_argument(
        "--dos-energy-window",
        dest="dos_energy_window",
        action="store_true",
        default=True,
        help="Use an automatic compact energy window for the DOS-curve PDFs.",
    )
    parser.add_argument(
        "--no-dos-energy-window",
        dest="dos_energy_window",
        action="store_false",
        help="Use the full catalogue energy axis for the DOS-curve PDFs.",
    )

    parser.add_argument(
        "--with-gap-plot",
        action="store_true",
        help="Also reconstruct and write the expensive Delta_eq(T) PDF.",
    )
    parser.add_argument(
        "--gap-pdf-name",
        default="E1_usadel_gap_eq_vs_temperature.pdf",
        help="Output PDF filename for the optional Delta_eq(T) figure.",
    )
    parser.add_argument(
        "--n-curves",
        type=int,
        default=4,
        help="Number of q curves for the optional Delta_eq(T) plot, including q=0 and q=q_c.",
    )
    parser.add_argument(
        "--n-temperature",
        type=int,
        default=64,
        help=(
            "Number of temperature samples for the optional Delta_eq(T) plot. "
            "Default is a fast smoke-style resolution."
        ),
    )
    parser.add_argument(
        "--T-min-K",
        type=float,
        default=None,
        help="Optional minimum temperature for the optional Delta_eq(T) plot. Defaults to PRE metadata T_bias_K.",
    )
    parser.add_argument(
        "--T-max-K",
        type=float,
        default=None,
        help=(
            "Optional maximum temperature for the optional Delta_eq(T) plot. Values at or below Tc "
            "are automatically extended to show Tc inside the axis."
        ),
    )
    parser.add_argument(
        "--q-critical-m-inv",
        type=float,
        default=None,
        help="Optional q_c override in m^-1 for the optional Delta_eq(T) plot.",
    )
    parser.add_argument(
        "--n-matsubara",
        type=int,
        default=96,
        help=(
            "Optional Matsubara cutoff for the optional Delta_eq(T) reconstruction. "
            "Default is a fast smoke-style value."
        ),
    )
    parser.add_argument(
        "--progress",
        dest="progress",
        action="store_true",
        default=False,
        help="Show a progress bar while solving the optional Delta_eq(T) Matsubara self-consistency points.",
    )
    parser.add_argument(
        "--no-progress",
        dest="progress",
        action="store_false",
        help="Disable the progress bar for the optional Delta_eq(T) solve.",
    )

    parser.add_argument(
        "--dpi",
        type=int,
        default=480,
        help="PDF rasterization DPI for any rasterized artists.",
    )
    parser.add_argument(
        "--title",
        default=None,
        help="Optional title for the optional Delta_eq(T) plot. Default keeps the figure title-free.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    cfg = validate_config(load_config(args.config))
    layout = create_run_layout(cfg, args.pre_run_name)

    raw_pre = Path(layout["raw_pre"])
    figures_dir = (
        Path(layout["plots_figures"]) / "E1_prerun"
        if args.output_dir is None
        else args.output_dir.expanduser().resolve()
    )
    figures_dir.mkdir(parents=True, exist_ok=True)

    catalog_path = (
        args.catalog_npz.expanduser().resolve()
        if args.catalog_npz is not None
        else raw_pre / "usadel_dos_catalog.npz"
    )
    if not catalog_path.exists():
        raise FileNotFoundError(f"Usadel catalog not found: {catalog_path}")

    saved: dict[str, Path] = {}
    usadel_catalog = load_usadel_catalog_npz(catalog_path)

    mesh_path = raw_pre / "mesh.npz"
    if not mesh_path.exists():
        raise FileNotFoundError(f"PRE mesh not found: {mesh_path}")
    saved["mesh_pytdgl_style_pdf"] = plot_mesh_pytdgl_style(
        load_mesh_npz(mesh_path),
        figures_dir / "mesh_pytdgl_style.pdf",
        dpi=int(args.dpi),
    )

    supercurrent_output = plot_usadel_supercurrent_curve(
        usadel_catalog,
        figures_dir / _ensure_pdf_name(args.supercurrent_pdf_name),
        dpi=int(args.dpi),
    )
    saved["E1_usadel_supercurrent_curve_pdf"] = supercurrent_output

    if not args.skip_dos_curves:
        dos_eq_output = plot_usadel_dos_curves_equilibrium_gap(
            usadel_catalog,
            figures_dir / _ensure_pdf_name(args.dos_eq_pdf_name),
            current_fractions=tuple(float(v) for v in args.dos_current_fractions),
            dpi=int(args.dpi),
            energy_max_meV=args.dos_energy_max_meV,
            energy_window=bool(args.dos_energy_window),
        )
        saved["usadel_dos_curves_delta_eq_pdf"] = dos_eq_output

        dos_delta0_output = plot_usadel_dos_curves_fixed_delta0(
            usadel_catalog,
            figures_dir / _ensure_pdf_name(args.dos_delta0_pdf_name),
            current_fractions=tuple(float(v) for v in args.dos_current_fractions),
            dpi=int(args.dpi),
            energy_max_meV=args.dos_energy_max_meV,
            energy_window=bool(args.dos_energy_window),
        )
        saved["usadel_dos_curves_delta0_pdf"] = dos_delta0_output

    if not args.skip_eliashberg_spectrum:
        eliashberg_path = _resolve_eliashberg_path(cfg, args.eliashberg_dat)
        spectrum = load_simon_eliashberg_dat(eliashberg_path)
        eliashberg_output = plot_eliashberg_spectrum(
            spectrum,
            figures_dir / _ensure_pdf_name(args.eliashberg_pdf_name),
            dpi=int(args.dpi),
        )
        saved["eliashberg_spectral_function_phdos_pdf"] = eliashberg_output

    power_table_path = (
        args.power_table_npz.expanduser().resolve()
        if args.power_table_npz is not None
        else raw_pre / "power_table_catalog.npz"
    )
    if not args.skip_power_table_figures:
        if not power_table_path.exists():
            raise FileNotFoundError(f"PRE power table not found: {power_table_path}")
        power_catalog = load_power_table_plot_catalog(power_table_path)
        saved["power_channels_Te_Tph_pdf"] = plot_power_channels_Te_Tph_maps(
            power_catalog,
            figures_dir / "E1_power_channels_Te_Tph_maps.pdf",
            dpi=int(args.dpi),
        )
        saved["power_exchange_vs_temperature_pdf"] = plot_power_total_Te_curves(
            power_catalog,
            figures_dir / "E1_power_exchange_vs_temperature.pdf",
            dpi=int(args.dpi),
        )
        saved["energy_heat_capacity_pdf"] = plot_energy_heat_capacity_curves(
            power_catalog,
            figures_dir / "E1_energy_heat_capacity_curves.pdf",
            dpi=int(args.dpi),
        )
        saved["electronic_thermal_conductivity_pdf"] = plot_electronic_thermal_conductivity_curves(
            power_catalog,
            figures_dir / "E1_electronic_thermal_conductivity_curves.pdf",
            dpi=int(args.dpi),
        )

    gap_source = "not_requested"
    q_critical_m_inv = None
    gap_metadata: dict[str, Any] = {}
    if args.with_gap_plot:
        gap_catalog = load_usadel_gap_catalog(
            catalog_path,
            n_curves=int(args.n_curves),
            n_temperature=int(args.n_temperature),
            T_min_K=args.T_min_K,
            T_max_K=args.T_max_K,
            q_critical_m_inv=args.q_critical_m_inv,
            n_matsubara=args.n_matsubara,
            progress=bool(args.progress),
        )
        gap_output = plot_gap_eq_vs_temperature(
            gap_catalog,
            figures_dir / _ensure_pdf_name(args.gap_pdf_name),
            dpi=int(args.dpi),
            title=args.title,
        )
        saved["usadel_gap_eq_vs_temperature_pdf"] = gap_output
        gap_source = gap_catalog.source_key
        q_critical_m_inv = float(gap_catalog.q_critical_m_inv)
        gap_metadata = dict(gap_catalog.metadata)

    manifest_path = _write_manifest(
        pre_run_name=args.pre_run_name,
        raw_pre=raw_pre,
        figures_dir=figures_dir,
        catalog_path=catalog_path,
        saved=saved,
        skip_dos_curves=bool(args.skip_dos_curves),
        dos_settings={
            "current_fractions": [float(v) for v in args.dos_current_fractions],
            "energy_max_meV": args.dos_energy_max_meV,
            "energy_window": bool(args.dos_energy_window),
        },
        skip_eliashberg_spectrum=bool(args.skip_eliashberg_spectrum),
        eliashberg_settings={
            "dat_path": str(_resolve_eliashberg_path(cfg, args.eliashberg_dat))
            if not args.skip_eliashberg_spectrum
            else None,
            "pdf_name": _ensure_pdf_name(args.eliashberg_pdf_name),
        },
        power_table_settings={
            "skipped": bool(args.skip_power_table_figures),
            "catalog_path": str(power_table_path),
        },
        with_gap_plot=bool(args.with_gap_plot),
        gap_source=gap_source,
        q_critical_m_inv=q_critical_m_inv,
        gap_metadata=gap_metadata,
        gap_settings={
            "n_curves": int(args.n_curves),
            "n_temperature": int(args.n_temperature),
            "n_matsubara": int(args.n_matsubara) if args.n_matsubara is not None else None,
            "T_min_K": args.T_min_K,
            "T_max_K": args.T_max_K,
            "progress": bool(args.progress),
        },
    )

    print("E1 pre-run Usadel plots")
    print(f" pre_run_name: {args.pre_run_name}")
    print(f" raw_pre: {raw_pre}")
    print(f" figures_dir: {figures_dir}")
    print(f" catalog_npz: {catalog_path}")
    print()
    print("Figures")
    for key, path in saved.items():
        print(f" {key}: {path}")
    print(f" manifest: {manifest_path}")
    print("Status: OK")
    return 0


def _resolve_eliashberg_path(cfg: dict[str, Any], requested_path: str | None) -> Path:
    """Resolve the Simon/MIT NbN alpha2F/PhDOS data path as in 01_prerun_template.py."""
    if requested_path:
        return Path(requested_path).expanduser().resolve()
    root = Path(str(cfg["project"]["big_data_root"])).expanduser()
    return (root / "catalogs" / "simon_2025" / "nbn-a2f-ph.dat").resolve()


def _ensure_pdf_name(name: str | Path) -> str:
    raw = str(name)
    return raw if raw.lower().endswith(".pdf") else f"{raw}.pdf"


def _write_manifest(
    *,
    pre_run_name: str,
    raw_pre: Path,
    figures_dir: Path,
    catalog_path: Path,
    saved: dict[str, Path],
    skip_dos_curves: bool,
    dos_settings: dict[str, Any],
    skip_eliashberg_spectrum: bool,
    eliashberg_settings: dict[str, Any],
    power_table_settings: dict[str, Any],
    with_gap_plot: bool,
    gap_source: str,
    q_critical_m_inv: float | None,
    gap_metadata: dict[str, Any],
    gap_settings: dict[str, Any],
) -> Path:
    manifest: dict[str, Any] = {
        "schema_version": 4,
        "pipeline": "plot_pipelines/E1_plot_prerun.py",
        "purpose": "E-type PRE figures in PDF format: supercurrent curve, DOS curves, Eliashberg/PhDOS spectrum, and optional Delta_eq(T,q).",
        "pre_run_name": pre_run_name,
        "raw_pre": str(raw_pre),
        "figures_dir": str(figures_dir),
        "catalog_npz": str(catalog_path),
        "figures": {key: str(path) for key, path in saved.items()},
        "skip_dos_curves": bool(skip_dos_curves),
        "dos_settings": dos_settings,
        "skip_eliashberg_spectrum": bool(skip_eliashberg_spectrum),
        "eliashberg_settings": eliashberg_settings,
        "power_table_settings": power_table_settings,
        "with_gap_plot": bool(with_gap_plot),
        "gap_source": gap_source,
        "q_critical_m_inv": None if q_critical_m_inv is None else float(q_critical_m_inv),
        "gap_settings": gap_settings,
        "gap_metadata": gap_metadata,
    }
    path = figures_dir / "E1_plot_prerun_manifest.yaml"
    with path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(manifest, f, sort_keys=False, allow_unicode=True, default_flow_style=False)
    return path


if __name__ == "__main__":
    raise SystemExit(main())
