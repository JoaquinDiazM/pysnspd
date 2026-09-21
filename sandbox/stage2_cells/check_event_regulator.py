"""Hold count quadrature fixed while reducing the causal spectral regulator."""
import json
import platform
import time
import numpy as np
from event_diagnostics import (ROOT, OUT, CATALOG, CRITERIA, FIELDS, PROFILES, sha,
                               phonon_grid, point_rates, power_summary)
from pysnspd.experimental.energy_catalog import OccupationEnergyCatalog, energy_at_count_batch, retarded_spectrum
from pysnspd.experimental.kinetic_events import ElectronPhononEvents


def main():
    if not CRITERIA.exists():
        raise SystemExit("Frozen criteria are required")
    start = time.perf_counter()
    catalog = OccupationEnergyCatalog.load(CATALOG)
    rows = []
    for label, amplitude, gamma in FIELDS:
        for eta in (1e-7, 1e-8, 1e-9):
            energy = energy_at_count_batch(catalog.count_nodes, delta=amplitude, gamma=gamma, eta=eta)
            c, s = retarded_spectrum(energy, delta=amplitude, gamma=gamma, eta=eta)
            network = ElectronPhononEvents.from_arrays(energy, catalog.count_weights, s.imag/c.real,
                phonon_grid(65), rate_prefactor=1., label=f"direct causal spectrum eta={eta}; unchanged180countnodes")
            for profile in PROFILES:
                rates, gross = point_rates(network, profile)
                rows.append({"field": label, "profile": profile, "eta": eta,
                             "powers": power_summary(network, rates, gross)})
    changes = []
    for row in rows:
        if row["eta"] != 1e-8:
            continue
        refined = next(r for r in rows if r["field"] == row["field"] and r["profile"] == row["profile"] and r["eta"] == 1e-9)
        for kind in ("scattering", "recombination"):
            a, b = row["powers"][kind], refined["powers"][kind]
            changes.append({"field": row["field"], "profile": row["profile"], "channel": kind,
                            "net_relative_change": abs(a["net"]-b["net"])/abs(b["net"]),
                            "gross_relative_change": abs(a["gross"]-b["gross"])/b["gross"]})
    result = {"schema": "pysnspd.stage2.event_regulator.v1", "host": platform.node(),
              "criteria_sha256": sha(CRITERIA), "catalog_sha256": sha(CATALOG),
              "runner_sha256": sha(__file__), "kernel_sha256": sha(ROOT / "pysnspd/experimental/kinetic_events.py"),
              "scope": "Direct spectra at the same native count nodes; point Bose occupations. Isolates regulator sensitivity, not electronic quadrature or field-table interpolation.",
              "rows": rows, "eta8_to_eta9_changes": changes,
              "runtime_seconds": time.perf_counter()-start}
    output = OUT / "event_regulator.json"
    output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print(json.dumps({"runtime_seconds": result["runtime_seconds"], "max_net_relative_change": max(r["net_relative_change"] for r in changes)}))


if __name__ == "__main__":
    main()
