"""Isolate fixed-phonon interpolation from electronic event quadrature.

The direct comparator retains exactly the same native Pauli representation and
events, replacing only Bose activities by their analytic value at event Omega.
Phonon RHS is compared as a measure on fixed hats, not as nodal values.
"""
from pathlib import Path
import argparse
import json
import platform
import time
import numpy as np
from scipy.special import expit

from event_diagnostics import ROOT, OUT, CATALOG, CRITERIA, PROFILES, sha, populations, project_measure, relative_error, bose
from pysnspd.experimental.cell_validation import ElectronicCell
from pysnspd.experimental.cell_closures import CellScales, DebyePhonons
from pysnspd.experimental.energy_catalog import OccupationEnergyCatalog
from pysnspd.experimental.refined_cells import refined_count_catalog
from pysnspd.experimental.kinetic_events import PhononGrid, ProjectedElectronPhononEvents, _signed_exponential_difference


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--nodes", type=int, nargs="+", default=[129, 257, 513])
    parser.add_argument("--infrared", type=float, default=.005)
    parser.add_argument("--output", type=Path, default=OUT / "projected_phonon_convergence.json")
    args = parser.parse_args()
    if not CRITERIA.exists():
        raise SystemExit("Frozen acceptance criteria are required")
    started = time.perf_counter()
    base = OccupationEnergyCatalog.load(CATALOG)
    catalog = refined_count_catalog(base, 1)
    # Unit scales do not alter the normalized test: n_atom/(N0 Delta0)=10.
    scales = CellScales(1e-22, 1e47, 10., 2., 1.)
    material = DebyePhonons(scales, 4., 10*scales.N0_per_J_m3*scales.delta0_J,
        .03, "algorithm-only synthetic Debye; not admitted NbN", args.infrared)
    rows, algebra = [], []
    for name, amplitude, gamma in (("normal", 0., 0.), ("gapped", .72, 0.), ("depaired", .72, .2), ("gapless", .35, .65)):
        cell = ElectronicCell(catalog, amplitude, gamma)
        for count in args.nodes:
            quadrature = material.quadrature(count)
            grid = PhononGrid(quadrature.energies, quadrature.capacities, material.alpha2F,
                material.coupling_support, material.synthetic_label)
            network = ProjectedElectronPhononEvents.from_cell(cell, grid, rate_prefactor=1.,
                quadrature_order=2, integration_layout="resolved_panels", outer_order=2,
                max_energy_panel=.125)
            for profile in PROFILES:
                p, n = populations(profile, cell.energies, grid.energies)
                _, n_exact = populations(profile, cell.energies, network.omega)
                logf, logr = network.log_activities(p, n)
                direct_logf = logf-network._barycentric_log(np.log1p(n))+np.log1p(n_exact)
                direct_logr = logr-network._barycentric_log(np.log(n))+np.log(n_exact)
                direct = network.coefficients*_signed_exponential_difference(direct_logf, direct_logr)
                direct_gross = network.coefficients*(np.exp(direct_logf)+np.exp(direct_logr))
                actual = network.rates(p, n)
                actual_gross = network.coefficients*(np.exp(logf)+np.exp(logr))
                for kind, mask in (("scattering", ~network.recombination), ("recombination", network.recombination)):
                    dp, dn = network.rhs_from_rates(np.where(mask, actual, 0.))
                    direct_dp, _ = network.rhs_from_rates(np.where(mask, direct, 0.))
                    reference_power = np.dot(network.omega[mask], direct[mask])
                    gross_power = np.dot(network.omega[mask], direct_gross[mask])
                    row = {"field": name, "amplitude": amplitude, "gamma": gamma,
                        "profile": profile, "channel": kind, "phonon_nodes": count,
                        "electron_nodes": len(cell.energies), "active_events": len(network.omega),
                        "direct_power": float(reference_power), "actual_power": float(np.dot(network.omega[mask], actual[mask])),
                        "power_relative_error": relative_error(np.array([np.dot(network.omega[mask], actual[mask])]), np.array([reference_power]), np.array([gross_power])),
                        "gross_power_relative_error": float(abs(np.dot(network.omega[mask], actual_gross[mask])-gross_power)/gross_power),
                        "electronic_rhs_relative_L1": relative_error(network.electron_capacities*dp, network.electron_capacities*direct_dp, direct_gross[mask]),
                        "energy_residual": network.energy_rate(dp, dn)}
                    for hats in (17, 33):
                        expected = project_measure(network.omega[mask], direct[mask], hats)
                        observed = project_measure(grid.energies, grid.capacities*dn, hats)
                        row[f"phonon_{hats}_hats_relative_L1"] = relative_error(observed, expected, direct_gross[mask])
                        row[f"phonon_{hats}_hats_reference"] = expected.tolist()
                        row[f"phonon_{hats}_hats_actual"] = observed.tolist()
                    rows.append(row)
            if count == args.nodes[-1]:
                entry = {"field": name, "equilibria": []}
                for temperature in (.08, .2, .5):
                    p, n = expit(-cell.energies/temperature), bose(grid.energies, temperature)
                    lf, lr = network.log_activities(p, n)
                    rates = network.rates(p, n)
                    gross = network.coefficients*(np.exp(lf)+np.exp(lr))
                    represented = gross > 1e-100
                    dp, dn = network.rhs_from_rates(rates)
                    entry["equilibria"].append({"temperature": temperature,
                        "maximum_relative_event_imbalance": float(np.max(abs(rates[represented])/gross[represented])),
                        "weighted_rhs_L1": float(np.dot(network.electron_capacities, abs(dp))+np.dot(grid.capacities, abs(dn)))})
                rng = np.random.default_rng(9721)
                p = rng.choice([0., .25, 1.], len(cell.energies))
                n = rng.choice([0., .3], len(grid.energies))
                dp, dn = network.rhs(p, n)
                entry["faces_inward"] = bool(np.all(dp[p == 0] >= 0) and np.all(dp[p == 1] <= 0) and np.all(dn[n == 0] >= 0))
                p = rng.uniform(.01, .9, len(cell.energies)); n = rng.uniform(.01, 1., len(grid.energies))
                lf, lr = network.log_activities(p, n)
                entry["entropy_production"] = float(np.dot(network.rates(p, n), lf-lr))
                entry["metadata"] = network.metadata
                algebra.append(entry)
    result = {"schema": "pysnspd.stage2.projected_phonon_convergence.v1", "host": platform.node(),
        "criteria_sha256": sha(CRITERIA), "catalog_sha256": sha(CATALOG), "runner_sha256": sha(__file__),
        "sources": {name: sha(ROOT / "pysnspd/experimental" / name) for name in ("kinetic_events.py", "refined_cells.py", "cell_closures.py")},
        "phonon_rule": "Gauss-Lobatto with actual Debye capacities used by coupled driver",
        "infrared": args.infrared, "nodes": args.nodes, "material": material.metadata(),
        "rows": rows, "algebra": algebra, "runtime_seconds": time.perf_counter()-started}
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    metrics = ("power_relative_error", "gross_power_relative_error", "electronic_rhs_relative_L1", "phonon_17_hats_relative_L1", "phonon_33_hats_relative_L1")
    print(json.dumps({"runtime_seconds": result["runtime_seconds"], "finest_maxima": {key: max(r[key] for r in rows if r["phonon_nodes"] == args.nodes[-1]) for key in metrics}}))


if __name__ == "__main__":
    main()
