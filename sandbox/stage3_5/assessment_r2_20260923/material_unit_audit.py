"""Read-only algebraic audit of public phonon tables. No collision RHS or dynamics.

Inputs are cached, exact bytes from the commit recorded below. Download URLs and
SHA256 hashes are written in the result; complete upstream files stay under tmp/.
No new material normalization is applied or admitted by this script.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform
import time

import numpy as np

REVISION = "5b6bd747f80016da5ccd51db73c110a8ecc6abf6"
EXPECTED_NBN_SHA256 = "9aeea0948033d771deedae40da0cb4dc59fef80ac6ea41ac4d3a67b180efc610"
H_EV_PS = 4.135667696923859e-3  # exact SI h/e converted to eV ps
HBAR_EV_PS = H_EV_PS / (2 * np.pi)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path: Path) -> tuple[np.ndarray, int]:
    raw = np.loadtxt(path)
    stable = raw[np.argsort(raw[:, 0], kind="stable")]
    _, index = np.unique(stable[:, 0], return_index=True)
    return stable[index], len(raw)


def moments(x: np.ndarray, a: np.ndarray) -> dict:
    positive = x > 0
    e = H_EV_PS * x
    area = float(np.trapezoid(a, e))
    energy_moment = float(np.trapezoid(e * a, e))
    rate = 2 * np.pi / HBAR_EV_PS * area
    return {
        "lambda_if_column2_is_dimensionless_alpha2F": float(
            2 * np.trapezoid(a[positive] / x[positive], x[positive])
        ),
        "alpha2F_energy_integral_eV_if_axis_is_nu_THz": area,
        "alpha2F_first_energy_moment_eV2_if_axis_is_nu_THz": energy_moment,
        "normal_empty_band_single_particle_emission_limit_ps_inverse": rate,
        "inverse_of_that_formal_limit_ps": 1 / rate,
        "warning": (
            "Conditional quadrature only: normal empty final band, no phonons, "
            "rho=1, coherence=1, E above phonon support. This is not a measured "
            "lifetime, a cascade time, a thermalization time, or a bound for the "
            "superconducting detector. No RHS or trajectory was evaluated."
        ),
    }


def run(cache: Path, repository: Path) -> dict:
    original = cache / "nbn-a2f-ph.dat"
    if sha(original) != EXPECTED_NBN_SHA256:
        raise ValueError("NbN source bytes differ from the stage1 frozen source")
    files = {}
    for path in sorted(cache.iterdir()):
        if path.is_file():
            files[path.name] = {
                "sha256": sha(path), "bytes": path.stat().st_size,
                "url": f"https://raw.githubusercontent.com/qnngroup/proj-KE-solver/{REVISION}/{path.name}",
            }
    tables = {}
    for name in ["Al-a2f-phdos.dat", "nb-a2f-phdos.dat", "TiN-prim-a2f-phdos.dat", "nbn-a2f-ph.dat"]:
        path = cache / name
        a, nraw = rows(path)
        x = a[:, 0]
        first = path.read_text().splitlines()[0].strip()
        tables[name] = {
            "first_line": first, "has_header": first.startswith("#"),
            "rows_raw": nraw, "rows_stable_first_unique": len(a),
            "integral_column3_on_stored_axis_raw": float(np.trapezoid(a[:, 2], x)),
            "negative_column3_nodes": int(np.count_nonzero(a[:, 2] < 0)),
            "axis_min": float(x[0]), "axis_max": float(x[-1]),
            "moments": moments(x, a[:, 1]),
        }
    shape = repository / "docs/implementation/stage1_closure/material_nbn_shape_v1.csv"
    v = np.loadtxt(shape, delimiter=",", skiprows=1)
    m = float(np.trapezoid(v[:, 2], v[:, 0]))
    h_mev_ps = 1000 * H_EV_PS
    result = {
        "schema": "pysnspd.stage3_5.r2.material_unit_audit.v1",
        "status": "DOCUMENTARY_AND_ALGEBRAIC_ONLY_NOT_ABSOLUTE_DOS_ADMISSION",
        "new_trajectories": 0, "new_rhs_evaluations": 0,
        "revision": REVISION, "script_sha256": sha(Path(__file__)),
        "preparation_note": (
            "An initial guard used the older wrapped-source hash e94f1127... and stopped "
            "before writing a result. Corrected to source_numeric_payload_sha256 from "
            "stage1_closure/material_results.json. Exact upstream git blob is "
            "d479ac67a90d86ba704597ad8f8a11e673e5dd27; no data were changed."
        ),
        "runtime": {"python": platform.python_version(), "numpy": np.__version__},
        "inputs": files, "tables": tables,
        "existing_shape": {
            "path": shape.relative_to(repository).as_posix(), "sha256": sha(shape),
            "policy": "Existing stage1 stable-first common-support derivative; read without modification",
            "integral_on_stored_axis": m,
            "conditional_modes_if_column3_states_per_THz": m,
            "conditional_modes_if_column3_states_per_meV": m * h_mev_ps,
            "conversion_ratio_between_those_hypotheses": h_mev_ps,
            "neither_hypothesis_is_selected": True,
            "moments": moments(v[:, 0], v[:, 1]),
        },
        "energy_axis_convention": {
            "formula": "E=h*nu=hbar*omega; omega=2*pi*nu",
            "h_eV_ps": H_EV_PS,
            "meV_per_THz_cycles_per_second": h_mev_ps,
            "code_h_eV_ps_using_upstream_rounded_hbar": 2 * np.pi * 658.2e-6,
        },
        "BG19_cell_count_comparison_not_SI25_metadata": {
            "structure": "delta-NbN rocksalt, 2 atoms per primitive formula unit, 6 branches",
            "theoretical_volume_A3_per_formula_unit": 21.67,
            "experimental_volume_A3_per_formula_unit": 21.16,
            "theoretical_formula_units_nm_inverse3": 1000 / 21.67,
            "experimental_formula_units_nm_inverse3": 1000 / 21.16,
            "theoretical_atoms_nm_inverse3": 2000 / 21.67,
            "experimental_atoms_nm_inverse3": 2000 / 21.16,
            "source": "BG19 Table I; these volumes do not identify the upstream NbN file's DFPT cell",
        },
    }
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=Path("tmp/stage3_5_r2/material_original"))
    parser.add_argument("--output", type=Path, default=Path("docs/implementation/stage3_5/assessment_r2_20260923/material/unit_audit.json"))
    args = parser.parse_args()
    start = time.monotonic()
    root = Path(__file__).resolve().parents[3]
    result = run(args.cache, root)
    result["elapsed_seconds"] = time.monotonic() - start
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.output.exists():
        raise FileExistsError(f"Preserve prior evidence: choose a new --output: {args.output}")
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "output": str(args.output), "seconds": result["elapsed_seconds"]}))


if __name__ == "__main__":
    main()
