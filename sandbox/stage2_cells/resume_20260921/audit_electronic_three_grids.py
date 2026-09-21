"""Separate electronic refinement, fixed Bose interpolation and eta bias.

The archived three-grid event matrix is read, never rebuilt. Only two BCS
finite-eta partner integrals missing from that historical comparison are new.
No change of tolerance or reinterpretation of a combined error as a pure
quadrature error is made here.
"""
from pathlib import Path
import json
import platform
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "docs/implementation/stage2"
REVIEW = BASE / "review"
OUT = BASE / "resume_20260921"
sys.path[:0] = [str(ROOT), str(REVIEW)]
import weak_population_finite_eta as reference
from pysnspd.experimental.energy_catalog import OccupationEnergyCatalog
from pysnspd.experimental.refined_cells import refined_count_catalog


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def main():
    started = time.perf_counter()
    source = BASE / "complementary_event_final.json"
    matrix = read(source)
    criteria = read(BASE / "acceptance_criteria.json")
    tolerance = criteria["continuous_consistency"]["relative_error_max"]
    budget = criteria["continuous_consistency"]["reference_refinement_relative_budget_max"]
    catalog_path = ROOT / "docs/implementation/stage1_r2/catalogs/occupation_catalog.npz"
    factory_path = ROOT / "pysnspd/experimental/refined_cells.py"
    if matrix["catalog_sha256"] != reference.sha(catalog_path) or matrix["count_refinement_code_sha256"] != reference.sha(factory_path):
        raise ValueError("immutable parent or measured count factory mismatch")
    parent = OccupationEnergyCatalog.load(catalog_path)
    references = {1: read(REVIEW / "weak_population_reference_finite_eta_bcs_1.json")}
    for refinement in (2, 4):
        catalog = refined_count_catalog(parent, refinement)
        levels = []
        for outer, inner in ((4, 32), (8, 64)):
            tick = time.perf_counter()
            grid, values = reference.weak_reference(catalog.count_nodes, .72, catalog.eta, outer, inner)
            levels.append(values)
            print(json.dumps({"refinement": refinement, "outer": outer, "inner": inner,
                              "seconds": time.perf_counter()-tick}), flush=True)
        rows = []
        for ip, profile in enumerate(reference.PROFILES):
            channels = {}
            for ic, channel in enumerate(("scattering", "recombination")):
                expected = levels[-1][ip, ic]
                scale = float(np.sum(abs(expected)))
                channels[channel] = {"number_rhs_reference": expected.tolist(), "reference_L1": scale,
                    "last_refinement_relative_L1": float(np.sum(abs(expected-levels[-2][ip, ic]))/scale)}
            rows.append({"amplitude": .72, "gamma": 0., "profile": profile,
                         "energy_nodes": grid.tolist(), "channels": channels})
        data = {"schema": "pysnspd.stage2.same_eta_bcs_weak_reference.v1", "count_refinement": refinement,
            "count_nodes": catalog.count_nodes.tolist(), "count_weights": catalog.count_weights.tolist(),
            "electron_states": len(catalog.count_nodes), "eta": catalog.eta, "omega_support": [.01, 4.],
            "orders": [[4, 32], [8, 64]], "rows": rows, "catalog_sha256": reference.sha(catalog_path),
            "wrapper_sha256": reference.sha(factory_path), "reference_source_sha256": reference.sha(Path(reference.__file__)),
            "script_sha256": reference.sha(Path(__file__))}
        path = OUT / f"same_eta_bcs_reference_{refinement}.json"
        path.write_text(json.dumps(data, indent=2, allow_nan=False)+"\n", encoding="utf-8")
        references[refinement] = data
    comparisons = []
    for refinement in (1, 2, 4):
        ideal = read(REVIEW / f"weak_population_reference_final_all_{refinement}.json")
        same_eta = references[refinement]
        for actual in (r for r in matrix["rows"] if r["refinement"] == refinement):
            if actual["metadata"]["coupling_support"] != [.01, 4.]:
                raise ValueError("this audit concerns only the preserved historical .01 matrix")
            profile, a, g = actual["profile"], actual["amplitude"], actual["gamma"]
            ideal_row = next(r for r in ideal["rows"] if (r["amplitude"], r["gamma"], r["profile"]) == (a, g, profile))
            eta_row = next(r for r in same_eta["rows"] if r["profile"] == profile) if (a, g) == (.72, 0.) else ideal_row
            for channel in ("scattering", "recombination"):
                r0 = np.asarray(ideal_row["channels"][channel]["number_rhs_reference"])
                reta = np.asarray(eta_row["channels"][channel]["number_rhs_reference"])
                scale = float(np.sum(abs(reta)))
                full, pauli, point = [np.asarray(actual[rep][channel]["electronic_number_rhs"])
                    for rep in ("full_operator", "native_pauli_point_bose", "point_activities")]
                parts = {"fixed_phonon_interpolation": full-pauli,
                         "pauli_reconstruction": pauli-point,
                         "electronic_quadrature_vs_matched_reference": point-reta,
                         "explicit_bcs_regulator_bias": reta-r0}
                electronic = float(np.sum(abs(pauli-reta))/scale)
                referror = eta_row["channels"][channel]["last_refinement_relative_L1"]
                comparisons.append({"field": actual["field"], "amplitude": a, "gamma": g, "profile": profile,
                    "channel": channel, "refinement": refinement, "electron_states": actual["native_count_nodes"],
                    "reference_eta": parent.eta if (a, g) == (.72, 0.) else 0.,
                    "reference_policy": "matched finite eta BCS" if (a, g) == (.72, 0.) else "ideal normal/gapless reference; gapless error retains any finite-eta residual",
                    "electronic_point_bose_error": electronic,
                    "combined_vs_ideal_error": float(np.sum(abs(full-r0))/np.sum(abs(r0))),
                    "combined_vs_matched_reference_error": float(np.sum(abs(full-reta))/scale),
                    "component_relative_L1": {key: float(np.sum(abs(value))/scale) for key, value in parts.items()},
                    "vector_decomposition_residual": float(np.sum(abs(sum(parts.values())-(full-r0)))/scale),
                    "reference_refinement_relative_L1": referror,
                    "electronic_absolute_threshold_pass": electronic <= tolerance and referror <= budget})
    series = []
    for first in (r for r in comparisons if r["refinement"] == 1):
        matched = [next(r for r in comparisons if (r["field"], r["profile"], r["channel"], r["refinement"])
                        == (first["field"], first["profile"], first["channel"], level)) for level in (1, 2, 4)]
        values = [r["electronic_point_bose_error"] for r in matched]
        # Do not turn a reference-accuracy floor into a roundoff exception.
        # The literal monotonic test and uncertainty budgets remain separate.
        decrease = values[-1] < values[-2]
        series.append({"field": first["field"], "profile": first["profile"], "channel": first["channel"],
            "electronic_point_bose_errors": values, "last_two_decrease": decrease,
            "all_three_decrease": values[2] < values[1] < values[0],
            "reference_budgets": [r["reference_refinement_relative_L1"] for r in matched],
            "all_absolute_thresholds_pass": all(r["electronic_absolute_threshold_pass"] for r in matched)})
    all_decrease = all(r["last_two_decrease"] for r in series)
    result = {"schema": "pysnspd.stage2.electronic_three_grid_decomposition.v1", "host": platform.node(),
        "status": "PASS_SEPARATED_ELECTRONIC_REFINEMENT" if all_decrease and all(r["all_absolute_thresholds_pass"] for r in series) else "UNRESOLVED_CASES_REMAIN",
        "scope": "Historical IR=.01 electronic quadrature refined separately from fixed phonon interpolation and explicit BCS regulator bias; not a claim that combined ideal-limit errors decrease",
        "selected_configuration_limitation": "The selected IR=.005/630/1025 static comparison is separately measured. This audit does not create three electronic grids at .005.",
        "comparisons": comparisons, "series": series, "runtime_seconds": time.perf_counter()-started,
        "maxima_by_electronic_grid": [max(r["electronic_point_bose_error"] for r in comparisons if r["refinement"] == level) for level in (1, 2, 4)],
        "nondecreasing_cases": [r for r in series if not r["last_two_decrease"]],
        "hashes": {"script": reference.sha(Path(__file__)), "reference_source": reference.sha(Path(reference.__file__)),
            "measured_matrix": reference.sha(source), "criteria": reference.sha(BASE / "acceptance_criteria.json"),
            "factory": reference.sha(factory_path), "catalog": reference.sha(catalog_path),
            "existing_same_eta_630_reference": reference.sha(REVIEW / "weak_population_reference_finite_eta_bcs_1.json")}}
    path = OUT / "electronic_three_grid_decomposition.json"
    path.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in ("status", "runtime_seconds", "maxima_by_electronic_grid", "nondecreasing_cases")}), flush=True)


if __name__ == "__main__":
    main()
