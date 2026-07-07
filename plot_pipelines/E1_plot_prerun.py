#!/usr/bin/env python3
"""E1 extra plot pipeline: Usadel equilibrium gap versus temperature.

Input: a normal pre-run produced by ``pipelines/01_prerun_template.py``.
Output: a PDF plot with six ``Delta_eq(T)`` curves from q=0 up to q_c.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any, Mapping

import yaml

# Allow running as ``python plot_pipelines/E1_plot_prerun.py`` from repo root.
_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from pysnspd.plotting.usadel_gap import load_usadel_gap_catalog, plot_gap_eq_vs_temperature


def _read_yaml(path: Path) -> Mapping[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")
    with path.open("r", encoding="utf-8") as fh:
        obj = yaml.safe_load(fh) or {}
    if not isinstance(obj, Mapping):
        raise ValueError(f"Config file did not parse as a mapping: {path}")
    return obj


def _get_nested(mapping: Mapping[str, Any], path: tuple[str, ...]) -> Any | None:
    current: Any = mapping
    for key in path:
        if not isinstance(current, Mapping) or key not in current:
            return None
        current = current[key]
    return current


def _find_data_root(config: Mapping[str, Any]) -> Path:
    """Find pySNSPD's data root in several historical config layouts."""

    candidates = (
        ("data_root",),
        ("data", "root"),
        ("data", "data_root"),
        ("paths", "data_root"),
        ("io", "data_root"),
        ("storage", "data_root"),
        ("run", "data_root"),
    )
    for path in candidates:
        value = _get_nested(config, path)
        if isinstance(value, str) and value.strip():
            return Path(value).expanduser()

    def walk(obj: Any) -> Path | None:
        if isinstance(obj, Mapping):
            for key, value in obj.items():
                key_l = str(key).lower()
                if key_l in {"data_root", "root_dir", "root"} and isinstance(value, str):
                    if "big_data" in value or value.startswith("/") or value.startswith("~"):
                        return Path(value).expanduser()
                found = walk(value)
                if found is not None:
                    return found
        return None

    found = walk(config)
    if found is not None:
        return found

    # Geminga default used throughout the current pySNSPD workflow.
    return Path("/home/jdiaz/scratch/big_data")


def _candidate_pre_dirs(data_root: Path, pre_run_name: str) -> list[Path]:
    return [
        data_root / "raw" / pre_run_name / "pre",
        data_root / "raw" / pre_run_name,
        data_root / pre_run_name / "pre",
        data_root / pre_run_name,
    ]


def _locate_pre_dir(data_root: Path, pre_run_name: str) -> Path:
    candidates = _candidate_pre_dirs(data_root, pre_run_name)
    for path in candidates:
        if path.exists() and path.is_dir():
            return path
    rendered = "\n".join(f"  - {path}" for path in candidates)
    raise FileNotFoundError(f"Could not locate pre-run directory. Tried:\n{rendered}")


def _score_catalog_path(path: Path) -> tuple[int, str]:
    name = path.name.lower()
    score = 0
    if "usadel" in name:
        score += 20
    if "catalog" in name:
        score += 10
    if "eq" in name or "equilibrium" in name:
        score += 6
    if "phase" in name:
        score -= 50
    if "summary" in name:
        score -= 20
    return -score, str(path)


def _locate_usadel_catalog(pre_dir: Path) -> Path:
    preferred = (
        pre_dir / "usadel_catalog.npz",
        pre_dir / "usadel_equilibrium_catalog.npz",
        pre_dir / "equilibrium_usadel_catalog.npz",
        pre_dir / "usadel_dos_catalog.npz",
    )
    for path in preferred:
        if path.exists():
            return path

    candidates = [path for path in pre_dir.rglob("*.npz") if "usadel" in path.name.lower()]
    candidates = [path for path in candidates if "phase" not in path.name.lower()]
    if not candidates:
        raise FileNotFoundError(f"No Usadel .npz catalog found under: {pre_dir}")
    return sorted(candidates, key=_score_catalog_path)[0]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="E1 extra plot: Delta_eq(T) from a pre-run Usadel catalog, saved as PDF.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--config", required=True, type=Path, help="pySNSPD YAML config used by the pre-run.")
    parser.add_argument("--pre-run-name", required=True, help="Name of the pre-run produced by 01_prerun_template.py.")
    parser.add_argument(
        "--pre-run-dir",
        type=Path,
        default=None,
        help="Optional direct path to the pre-run directory; overrides automatic raw/<pre>/pre lookup.",
    )
    parser.add_argument(
        "--catalog-npz",
        type=Path,
        default=None,
        help="Optional direct path to the Usadel catalog .npz; overrides pre-run directory search.",
    )
    parser.add_argument("--output-dir", type=Path, default=None, help="Optional output directory for the PDF.")
    parser.add_argument("--pdf-name", default="E1_usadel_gap_eq_vs_temperature.pdf", help="Output PDF filename.")
    parser.add_argument("--n-curves", type=int, default=6, help="Number of q curves, including q=0 and q=q_c.")
    parser.add_argument(
        "--q-critical-m-inv",
        type=float,
        default=None,
        help="Optional q_c override in m^-1. If omitted, it is estimated from the low-T gap branch.",
    )
    parser.add_argument("--dpi", type=int, default=480, help="PDF rasterization DPI for any rasterized artists.")
    parser.add_argument("--title", default=None, help="Optional figure title. Default keeps the memory-ready figure title-free.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = _read_yaml(args.config)
    data_root = _find_data_root(config)

    if args.catalog_npz is not None:
        catalog_path = args.catalog_npz.expanduser().resolve()
        pre_dir = args.pre_run_dir.expanduser().resolve() if args.pre_run_dir is not None else catalog_path.parent
    else:
        pre_dir = (
            args.pre_run_dir.expanduser().resolve()
            if args.pre_run_dir is not None
            else _locate_pre_dir(data_root, args.pre_run_name)
        )
        catalog_path = _locate_usadel_catalog(pre_dir)

    if args.output_dir is not None:
        output_dir = args.output_dir.expanduser().resolve()
    else:
        output_dir = data_root / "plots" / args.pre_run_name / "E1_prerun"

    pdf_name = args.pdf_name if args.pdf_name.lower().endswith(".pdf") else f"{args.pdf_name}.pdf"
    output_path = output_dir / pdf_name

    catalog = load_usadel_gap_catalog(catalog_path)
    output = plot_gap_eq_vs_temperature(
        catalog,
        output_path,
        n_curves=args.n_curves,
        q_critical_m_inv=args.q_critical_m_inv,
        dpi=args.dpi,
        title=args.title,
    )

    print("E1 pre-run Usadel gap plot")
    print(f" pre_run_name: {args.pre_run_name}")
    print(f" pre_dir:      {pre_dir}")
    print(f" catalog_npz:  {catalog_path}")
    print(f" gap_field:    {catalog.source_key}")
    print(f" output_pdf:   {output}")


if __name__ == "__main__":
    main()
