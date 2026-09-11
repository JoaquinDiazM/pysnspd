"""Measure loading and 1000 fixed-p queries of existing experimental catalogs.

This script never builds, modifies, or optimizes a catalog. Run only after the
final NPZ artifacts are complete. The timings describe these API calls on this
machine, not the duration of a future coupled detector transient.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import statistics
import sys
import time

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import pysnspd.experimental.energy_catalog as catalog_module
from pysnspd.experimental.energy_catalog import OccupationEnergyCatalog, UniformVacuumCatalog

QUERY_COUNT = 1000
REPETITIONS = 3
WARMUP_QUERIES = 12
SEED = 20260911


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def array_sha256(array: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(array, dtype="<f8").tobytes(order="C")).hexdigest()


def describe_path(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT).as_posix()
    except ValueError:
        return str(path.resolve())


def benchmark_load(path: Path, loader):
    durations = []
    loaded = None
    for _ in range(REPETITIONS):
        started = time.perf_counter()
        loaded = loader(path)
        durations.append(time.perf_counter() - started)
    return loaded, {
        "seconds_by_repetition": durations,
        "median_seconds": statistics.median(durations),
        "median_ms_per_load": 1000 * statistics.median(durations),
        "includes": "NPZ decompression, array validation, and construction of interpolation objects",
        "cache_context": "same process; SHA-256 reads precede loading; no cold-storage claim",
    }


def benchmark_queries(evaluate, queries: np.ndarray):
    for delta, gamma in queries[:WARMUP_QUERIES]:
        result = np.asarray(evaluate(float(delta), float(gamma)))
        if result.shape != (3,) or not np.all(np.isfinite(result)):
            raise FloatingPointError("The warmed query did not return three finite values.")
    durations = []
    checksums = []
    for _ in range(REPETITIONS):
        # Allocate outside the timed loop; include the unmodified scalar API
        # call and storage of its three outputs, without vectorization.
        results = np.empty((QUERY_COUNT, 3), dtype=float)
        started = time.perf_counter()
        for index, (delta, gamma) in enumerate(queries):
            results[index] = evaluate(float(delta), float(gamma))
        durations.append(time.perf_counter() - started)
        if not np.all(np.isfinite(results)):
            raise FloatingPointError("A timed query returned a nonfinite value.")
        checksums.append(array_sha256(results))
    if len(set(checksums)) != 1:
        raise AssertionError("Identical query batches produced different results.")
    median = statistics.median(durations)
    return {
        "queries_per_repetition": QUERY_COUNT,
        "warmup_queries": WARMUP_QUERIES,
        "seconds_by_repetition": durations,
        "median_seconds": median,
        "median_ms_per_query": 1000 * median / QUERY_COUNT,
        "result_sha256_float64_little_endian": checksums[0],
        "repeat_results_identical": True,
        "timed_operation": "scalar evaluate API returning energy, amplitude force, and Gamma conjugate",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    default_root = ROOT / "docs/implementation/stage1"
    parser.add_argument("--vacuum", type=Path, default=default_root / "catalogs/vacuum_catalog.npz")
    parser.add_argument("--occupation", type=Path, default=default_root / "catalogs/occupation_catalog.npz")
    parser.add_argument("--output", type=Path, default=default_root / "query_benchmark.json")
    args = parser.parse_args()
    started_total = time.perf_counter()
    module_path = Path(catalog_module.__file__).resolve()
    hashes_before = {"module": sha256(module_path), "vacuum": sha256(args.vacuum),
                     "occupation": sha256(args.occupation)}
    vacuum, vacuum_loading = benchmark_load(args.vacuum, UniformVacuumCatalog.load)
    occupation, occupation_loading = benchmark_load(args.occupation, OccupationEnergyCatalog.load)

    support = {
        "delta": [max(float(vacuum.delta_axis[0]), float(occupation.vacuum.delta_axis[0])),
                  min(float(vacuum.delta_axis[-1]), float(occupation.vacuum.delta_axis[-1]))],
        "gamma": [max(float(vacuum.gamma_axis[0]), float(occupation.vacuum.gamma_axis[0])),
                  min(float(vacuum.gamma_axis[-1]), float(occupation.vacuum.gamma_axis[-1]))],
    }
    if any(not lower < upper for lower, upper in support.values()):
        raise ValueError("The two final catalogs have no common two-dimensional query support.")
    rng = np.random.default_rng(SEED)
    fractions = rng.uniform(0.001, 0.999, size=(QUERY_COUNT, 2))
    low = np.array([support["delta"][0], support["gamma"][0]])
    high = np.array([support["delta"][1], support["gamma"][1]])
    queries = low + fractions * (high - low)
    count_cutoff = float(np.sum(occupation.count_weights))
    reduced_count = occupation.count_nodes / count_cutoff
    fixed_p = (0.12 * np.exp(-0.5 * ((reduced_count - 0.35) / 0.08) ** 2)
               + 0.07 * np.exp(-0.5 * ((reduced_count - 0.75) / 0.07) ** 2))
    fixed_p.setflags(write=False)

    vacuum_queries = benchmark_queries(vacuum.evaluate, queries)
    occupation_queries = benchmark_queries(
        lambda delta, gamma: occupation.evaluate(delta, gamma, fixed_p), queries
    )
    hashes_after = {"module": sha256(module_path), "vacuum": sha256(args.vacuum),
                    "occupation": sha256(args.occupation)}
    if hashes_before != hashes_after:
        raise RuntimeError("A catalog or its module changed during the benchmark; discard this timing.")
    result = {
        "schema": "pysnspd.experimental.catalog-query-benchmark.v1",
        "runtime_seconds": time.perf_counter() - started_total,
        "platform": platform.platform(), "machine": platform.machine(),
        "processor": platform.processor(), "python": platform.python_version(),
        "numpy": np.__version__, "scipy": scipy.__version__,
        "thread_environment": {name: os.environ.get(name) for name in
                               ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")},
        "seed": SEED, "repetitions": REPETITIONS,
        "query_support_dimensionless": support,
        "query_coordinates_sha256_float64_little_endian": array_sha256(queries),
        "common_queries": True,
        "fixed_occupation": {
            "description": "Two nonthermal Gaussian populations fixed at the state-count quadrature nodes.",
            "formula": "p(x)=0.12 exp[-0.5((x/X-0.35)/0.08)^2]+0.07 exp[-0.5((x/X-0.75)/0.07)^2]",
            "X_count_cutoff_over_Delta0": count_cutoff,
            "count_nodes": len(fixed_p), "minimum": float(fixed_p.min()), "maximum": float(fixed_p.max()),
            "occupation_sha256_float64_little_endian": array_sha256(fixed_p),
            "thermal_reduction": False,
        },
        "module": {"path": describe_path(module_path), "sha256": hashes_before["module"]},
        "benchmark_script_sha256": sha256(Path(__file__)),
        "catalogs": {
            "vacuum": {"path": describe_path(args.vacuum), "bytes": args.vacuum.stat().st_size,
                       "sha256": hashes_before["vacuum"], "loading": vacuum_loading, "queries": vacuum_queries},
            "occupation": {"path": describe_path(args.occupation), "bytes": args.occupation.stat().st_size,
                           "sha256": hashes_before["occupation"], "loading": occupation_loading,
                           "queries": occupation_queries},
        },
        "inputs_unchanged": True,
        "interpretation": "Measured cost of loading and scalar queries only; not a coupled transient runtime estimate.",
        "excluded_work": "catalog construction, spectral solves, kinetic collisions, spatial PDE, circuit evolution, and transient time stepping",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"runtime_seconds": result["runtime_seconds"],
                      "vacuum_median_ms_per_query": vacuum_queries["median_ms_per_query"],
                      "occupation_median_ms_per_query": occupation_queries["median_ms_per_query"],
                      "output": str(args.output)}))


if __name__ == "__main__":
    main()
