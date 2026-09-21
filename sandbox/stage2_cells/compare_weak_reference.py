"""Read-only comparison of event deposition against independent weak RHS data."""
from pathlib import Path
import argparse
import hashlib
import json

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "docs/implementation/stage2"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--reference", type=Path, default=OUT / "review/weak_population_reference_1.json")
    parser.add_argument("--order", type=int, default=2)
    parser.add_argument("--refinement", type=int)
    parser.add_argument("--output", type=Path, default=OUT / "weak_reference_comparison.json")
    args = parser.parse_args()
    candidate = json.loads(args.candidate.read_text(encoding="utf-8"))
    reference = json.loads(args.reference.read_text(encoding="utf-8"))
    for actual_key, reference_key in (("criteria_sha256", "criteria_sha256"),
                                     ("catalog_sha256", "catalog_sha256"),
                                     ("count_refinement_code_sha256", "wrapper_sha256")):
        if candidate[actual_key] != reference[reference_key]:
            raise ValueError(f"candidate and reference disagree on {reference_key}")
    rows = []
    for row in candidate["rows"]:
        if row["quadrature_order"] != args.order:
            continue
        if args.refinement is not None and row["refinement"] != args.refinement:
            continue
        if row["refinement"] != reference["count_refinement"]:
            raise ValueError("select the refinement matching the independent reference")
        if row["metadata"]["coupling_support"] != reference["omega_support"]:
            raise ValueError("candidate and reference must use the same phonon support")
        matches = [r for r in reference["rows"] if r["amplitude"] == row["amplitude"]
                   and r.get("gamma", 0.) == row["gamma"] and r["profile"] == row["profile"]]
        if not matches:
            continue
        ref = matches[0]
        for kind in ("scattering", "recombination"):
            for representation in ("point_activities", "native_pauli_point_bose", "full_operator"):
                observed = row[representation][kind]["electronic_number_rhs"]
                expected = ref["channels"][kind]["number_rhs_reference"]
                if len(observed) != len(expected):
                    raise ValueError("unmatched native grids cannot be compared")
                difference = sum(abs(a-b) for a, b in zip(observed, expected))
                rows.append({"field": row["field"], "profile": row["profile"], "channel": kind,
                             "representation": representation,
                             "relative_L1": difference/ref["channels"][kind]["reference_L1"]})
    result = {"schema": "pysnspd.stage2.weak_reference_comparison.v1",
              "candidate_sha256": sha(args.candidate), "reference_sha256": sha(args.reference),
              "script_sha256": sha(__file__), "order": args.order, "refinement": args.refinement, "rows": rows}
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    worst = max((r["relative_L1"], r["field"], r["profile"], r["channel"]) for r in rows if r["representation"] == "full_operator")
    print(json.dumps({"worst_full_operator_relative_L1": worst, "rows": len(rows)}))


if __name__ == "__main__":
    main()
