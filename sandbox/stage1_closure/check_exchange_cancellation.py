"""Reproduce the rejected subtractive exchange formula without changing the solver.

The historical formula is algebraically exact, but not numerically convex.
This diagnostic preserves its failure separately from the accepted cell outputs.
"""
from pathlib import Path
import hashlib
import json
import platform
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import numpy as np
from scipy.special import expit
from pysnspd.experimental.cell_validation import ElectronicCell, EnergyFacePair
from pysnspd.experimental.energy_catalog import OccupationEnergyCatalog


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def subtractive_exchange(pair, population, time):
    """Rejected mean plus/minus difference expression, with no clipping."""
    result = population.copy()
    cl, cr = pair.capacities[:, pair.active]
    left, right = population[:, pair.active]
    mean = (cl*left+cr*right)/(cl+cr)
    difference = (left-right)*np.exp(-pair.conductance[pair.active]*(1/cl+1/cr)*time)
    result[0, pair.active] = mean+cr/(cl+cr)*difference
    result[1, pair.active] = mean-cl/(cl+cr)*difference
    return result


def main():
    start = time.perf_counter()
    catalog_path = ROOT / "docs/implementation/stage1_r2/catalogs/occupation_catalog.npz"
    catalog = OccupationEnergyCatalog.load(catalog_path)
    left, right = ElectronicCell(catalog, .43, .02), ElectronicCell(catalog, 1.2, .2)
    common_upper = min(left.energies[-1], right.energies[-1])
    rows = []
    for nodes in (65, 513):
        grid = np.r_[np.linspace(1e-6, common_upper, nodes),
                     max(left.energies[-1], right.energies[-1])]
        pair = EnergyFacePair(left, right, grid)
        cases = {
            "thermal_initial_runner": np.array([expit(-grid/.4), expit(-grid/.12)]),
            "cold_regression": np.array([np.linspace(.01, .2, len(grid)),
                                          np.full_like(grid, 1e-44)]),
        }
        for name, population in cases.items():
            for elapsed in (0., 1e-14):
                old = subtractive_exchange(pair, population, elapsed)
                new = pair.exact(population, elapsed)
                rows.append({
                    "scenario": name, "common_nodes": nodes, "time": elapsed,
                    "initial_minimum": float(np.min(population)),
                    "subtractive_minimum": float(np.min(old)),
                    "subtractive_negative_entries": int(np.count_nonzero(old < 0)),
                    "subtractive_max_population_change": float(np.max(abs(old-population))),
                    "convex_minimum": float(np.min(new)),
                    "convex_negative_entries": int(np.count_nonzero(new < 0)),
                    "convex_zero_time_bitwise_identity": bool(np.array_equal(new, population)) if elapsed == 0 else None,
                })
    result = {
        "schema": "pysnspd.stage1_closure.exchange_cancellation.v1",
        "host": platform.node(), "python": platform.python_version(),
        "runtime_seconds": time.perf_counter()-start,
        "catalog_sha256": sha(catalog_path),
        "cell_code_sha256": sha(ROOT / "pysnspd/experimental/cell_validation.py"),
        "script_sha256": sha(__file__),
        "provenance": "Fresh reproduction of the rejected algebraic form using the final unchanged interface; not a recovered historical full run.",
        "regression_test": "tests/test_experimental_cell_validation.py::test_exact_exchange_does_not_cancel_a_cold_population_at_zero_time",
        "rows": rows,
    }
    output = ROOT / "docs/implementation/stage1_closure/pilot_initial/exchange_cancellation.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"output": str(output.relative_to(ROOT)),
                      "runtime_seconds": result["runtime_seconds"],
                      "subtractive_minimum": min(r["subtractive_minimum"] for r in rows),
                      "convex_negative_entries": sum(r["convex_negative_entries"] for r in rows)}))


if __name__ == "__main__":
    main()
