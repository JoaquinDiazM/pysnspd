"""Assemble already-measured event evidence without running new physics."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "docs/implementation/stage2"


def load(name):
    return json.loads((OUT / name).read_text(encoding="utf-8"))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    comparison = []
    for refinement in (1, 2, 4):
        data = load(f"weak_reference_final_comparison_{refinement}.json")
        for representation in ("point_activities", "native_pauli_point_bose", "full_operator"):
            worst = max((r for r in data["rows"] if r["representation"] == representation),
                        key=lambda r: r["relative_L1"])
            comparison.append({"refinement": refinement, **worst})
    grid = load("projected_phonon_convergence.json")
    finer = load("projected_phonon_1025.json")
    metrics = ("power_relative_error", "gross_power_relative_error", "electronic_rhs_relative_L1",
               "phonon_17_hats_relative_L1", "phonon_33_hats_relative_L1")
    phonons = []
    for count in (129, 257, 513, 1025):
        rows = [r for r in (grid["rows"]+finer["rows"]) if r["phonon_nodes"] == count]
        maxima = {}
        for metric in metrics:
            worst = max(rows, key=lambda r: r[metric])
            maxima[metric] = {key: worst[key] for key in ("field", "profile", "channel", metric)}
        phonons.append({"nodes": count, "maxima": maxima})
    selected = load("selected_event_measurements.json")
    expensive = load("complementary_event_final.json")
    files = ["complementary_event_final.json", "projected_phonon_convergence.json",
             "projected_phonon_1025.json", "selected_event_measurements.json"]
    files += [f"weak_reference_final_comparison_{r}.json" for r in (1, 2, 4)]
    result = {
        "schema": "pysnspd.stage2.event_evidence_summary.v1",
        "status": "STATIC_EVIDENCE_ONLY_STAGE2_NOT_CLOSED",
        "limitations": [
            "The continuous electronic reference and final three-grid comparison use IR=.01; the selected dynamical candidate uses .005.",
            "Selected IR=.005 has isolated phonon interpolation and inner-quadrature checks; its exact same-cut independent continuous comparison is still pending.",
            "The weak BCS error against ideal eta0 has a small nondecreasing floor; finite-eta reference work is separate and must not be called electronic-quadrature convergence.",
            "A temporal reference elsewhere in stage2 exhausted its bounded run; no complete stage2 admission or production promotion follows from these static checks.",
        ],
        "criteria_sha256": sha(OUT / "acceptance_criteria.json"),
        "kernel_sha256": sha(ROOT / "pysnspd/experimental/kinetic_events.py"),
        "factory_sha256": sha(ROOT / "pysnspd/experimental/refined_cells.py"),
        "electronic_weak_maxima": comparison,
        "phonon_interpolation_maxima": phonons,
        "selected_inner_quadrature_maximum_relative_L1": max(r["electronic_weighted_L1_relative_error"] for r in selected["quadrature_self_convergence"]),
        "selected_costs": [{key: r[key] for key in ("field", "native_count_nodes", "active_events", "network_construction_seconds", "rhs_milliseconds_five_queries")}
            for r in selected["rows"] if r["quadrature_order"] == 2 and r["profile"] == "hot_electrons"],
        "three_grid_runtime_seconds": expensive["runtime_seconds"],
        "selected_static_runtime_seconds": selected["runtime_seconds"],
        "artifact_sha256": {name: sha(OUT / name) for name in files},
        "script_sha256": sha(__file__),
    }
    (OUT / "event_summary.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "artifacts": len(files)}))


if __name__ == "__main__":
    main()
