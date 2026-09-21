"""Separate native reconstruction, nested quadrature and activity errors."""
from pathlib import Path
import argparse
import json
import platform
import time
import numpy as np

from event_diagnostics import (ROOT, OUT, CATALOG, CRITERIA, PROFILES, sha,
                               populations, phonon_grid, project_measure, relative_error)
from pysnspd.experimental.cell_validation import ElectronicCell
from pysnspd.experimental.energy_catalog import OccupationEnergyCatalog
from pysnspd.experimental.refined_cells import refined_count_catalog
from pysnspd.experimental.cell_closures import CellScales, DebyePhonons
from pysnspd.experimental.kinetic_events import ProjectedElectronPhononEvents, HybridElectronPhononEvents, _signed_exponential_difference
from pysnspd.experimental.kinetic_events import PhononGrid

REFERENCE = OUT / "review/continuous_reactions.json"
PROJECTION_REFERENCE = OUT / "review/continuous_projection_check.json"


def channel_observables(network, rates, gross):
    result = {}
    for kind, mask in (("scattering", ~network.recombination), ("recombination", network.recombination)):
        second_moment = np.where(network.recombination,
                                 -network.target_energy_i**2-network.target_energy_j**2,
                                 network.target_energy_i**2-network.target_energy_j**2)
        dp, dn = network.rhs_from_rates(np.where(mask, rates, 0.))
        result[kind] = {"observables": [float(np.sum(network.omega[mask]*rates[mask])),
                                         float(np.sum(network.omega[mask]*gross[mask])),
                                         float(np.sum(rates[mask])), float(np.sum(gross[mask])),
                                         float(np.sum(second_moment[mask]*rates[mask]))],
                        "native_second_moment": float(np.dot(network.electron_capacities*network.electron_energies**2, dp)),
                        "electronic_number_rhs": (network.electron_capacities*dp).tolist(),
                        "target_projected_17": project_measure(network.omega[mask], rates[mask]).tolist(),
                        "projected_17": project_measure(network.phonons.energies, network.phonons.capacities*dn).tolist(),
                        "projected_33": project_measure(network.phonons.energies, network.phonons.capacities*dn, 33).tolist()}
    return result


def run(orders, phonon_nodes, method, refinement=None, integration_layout="global", max_count_panel=.5,
        phonon_rule="uniform", infrared=.01, reference_path=REFERENCE, projection_path=PROJECTION_REFERENCE,
        measurement_only=False):
    base = OccupationEnergyCatalog.load(CATALOG)
    catalog = base if refinement is None else refined_count_catalog(base, refinement)
    reference = json.loads(reference_path.read_text(encoding="utf-8"))
    projection_reference = json.loads(projection_path.read_text(encoding="utf-8"))
    if not measurement_only and reference["omega_support"] != [infrared, 4.]:
        raise ValueError("continuous reference must have exactly the measured coupling support")
    if phonon_rule == "lobatto":
        scales = CellScales(1e-22, 1e47, 10., 2., 1.)
        material = DebyePhonons(scales, 4., 10*scales.N0_per_J_m3*scales.delta0_J,
            .03, "algorithm-only synthetic Debye; not admitted NbN", infrared)
        quadrature = material.quadrature(phonon_nodes)
        grid = PhononGrid(quadrature.energies, quadrature.capacities, material.alpha2F,
            material.coupling_support, material.synthetic_label)
    else:
        grid = phonon_grid(phonon_nodes, infrared)
    rows = []
    for label, amplitude, gamma in (("normal", 0., 0.), ("gapped", .72, 0.), ("gapless", .35, .65)):
        cell_started = time.perf_counter()
        cell = ElectronicCell(catalog, amplitude, gamma)
        cell_seconds = time.perf_counter()-cell_started
        for order in orders:
            constructor = HybridElectronPhononEvents if method == "hybrid" else ProjectedElectronPhononEvents
            network_started = time.perf_counter()
            network = constructor.from_cell(cell, grid,
                rate_prefactor=1., quadrature_order=order, integration_layout=integration_layout,
                max_count_panel=max_count_panel)
            network_seconds = time.perf_counter()-network_started
            for profile in PROFILES:
                p_native, n_grid = populations(profile, cell.energies, network.phonons.energies)
                pi, n_exact = populations(profile, network.target_energy_i, network.omega)
                pj, _ = populations(profile, network.target_energy_j, network.omega)
                forward = np.where(network.recombination, pi*pj, (1-pi)*pj)*(1+n_exact)
                reverse = np.where(network.recombination, (1-pi)*(1-pj), pi*(1-pj))*n_exact
                point_rate = network.coefficients*(forward-reverse)
                point_gross = network.coefficients*(forward+reverse)
                # Remove only the grid Bose activity from the full log rate,
                # replacing it with the direct Bose value at the event energy.
                logf, logr = network.log_activities(p_native, n_grid)
                direct_logf = logf-network._barycentric_log(np.log1p(n_grid))+np.log1p(n_exact)
                direct_logr = logr-network._barycentric_log(np.log(n_grid))+np.log(n_exact)
                pauli_rate = network.coefficients*_signed_exponential_difference(direct_logf, direct_logr)
                pauli_gross = network.coefficients*(np.exp(direct_logf)+np.exp(direct_logr))
                actual_rate = network.rates(p_native, n_grid)
                actual_gross = network.coefficients*(np.exp(logf)+np.exp(logr))
                query_started = time.perf_counter()
                for _ in range(5):
                    network.rhs(p_native, n_grid)
                rhs_milliseconds = (time.perf_counter()-query_started)*1000/5
                row = {"field": label, "amplitude": amplitude, "gamma": gamma,
                       "profile": profile, "quadrature_order": order, "phonon_nodes": phonon_nodes,
                       "phonon_rule": phonon_rule, "infrared": infrared,
                       "refinement": refinement, "native_count_nodes": len(catalog.count_nodes),
                       "active_events": len(network.omega), "cell_construction_seconds": cell_seconds,
                       "network_construction_seconds": network_seconds,
                       "rhs_milliseconds_five_queries": rhs_milliseconds,
                       "metadata": network.metadata,
                       "point_activities": channel_observables(network, point_rate, point_gross),
                       "native_pauli_point_bose": channel_observables(network, pauli_rate, pauli_gross),
                       "full_operator": channel_observables(network, actual_rate, actual_gross)}
                if measurement_only:
                    row["reference_status"] = "NOT_COMPARED: same-cut independent reference is still required"
                    dp, dn = network.rhs_from_rates(actual_rate)
                    row["instantaneous_energy_residual"] = network.energy_rate(dp, dn)
                    rows.append(row)
                    continue
                ref = next(r for r in reference["rows"] if r["amplitude"] == amplitude and r["gamma"] == gamma and r["profile"] == profile)
                errors = {}
                for representation in ("point_activities", "native_pauli_point_bose", "full_operator"):
                    errors[representation] = {}
                    for kind in ("scattering", "recombination"):
                        observed = np.array(row[representation][kind]["observables"])
                        expected = np.array(ref["channels"][kind]["reference"])
                        errors[representation][kind] = {"relative_observables": (abs(observed-expected)/np.maximum(abs(expected), 1e-300)).tolist()}
                        for hats in (17, 33):
                            projected_ref = next(r for r in projection_reference["rows"] if r["amplitude"] == amplitude
                                and r["gamma"] == gamma and r["profile"] == profile and r["hat_nodes"] == hats)
                            errors[representation][kind][f"projected_{hats}_relative_L1"] = relative_error(
                                np.array(row[representation][kind][f"projected_{hats}"]),
                                np.array(projected_ref["channels"][kind]["reference"]), np.array([observed[3]]))
                row["errors_against_independent_eta_zero"] = errors
                dp, dn = network.rhs_from_rates(actual_rate)
                row["instantaneous_energy_residual"] = network.energy_rate(dp, dn)
                rows.append(row)
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--orders", type=int, nargs="+", default=[24, 48, 96])
    parser.add_argument("--phonon-nodes", type=int, default=513)
    parser.add_argument("--phonon-rule", choices=["uniform", "lobatto"], default="uniform")
    parser.add_argument("--infrared", type=float, default=.01)
    parser.add_argument("--reference", type=Path, default=REFERENCE)
    parser.add_argument("--projection-reference", type=Path, default=PROJECTION_REFERENCE)
    parser.add_argument("--measurement-only", action="store_true",
        help="store measured operators and self-convergence without comparing unlike reference supports")
    parser.add_argument("--method", choices=["projected", "hybrid"], default="projected")
    parser.add_argument("--refinements", type=int, nargs="+")
    parser.add_argument("--integration-layout", choices=["global", "native_intervals", "energy_panels", "resolved_panels"], default="global")
    parser.add_argument("--max-count-panel", type=float, default=.5)
    parser.add_argument("--output", type=Path, default=OUT / "projected_event_pilot_01.json")
    args = parser.parse_args()
    if not CRITERIA.exists() or not args.reference.exists() or not args.projection_reference.exists():
        raise SystemExit("Frozen criteria and independent reference are required")
    started = time.perf_counter()
    rows = []
    for refinement in args.refinements or [None]:
        rows.extend(run(args.orders, args.phonon_nodes, args.method, refinement,
                        args.integration_layout, args.max_count_panel, args.phonon_rule,
                        args.infrared, args.reference, args.projection_reference, args.measurement_only))
    convergence = []
    if len(args.orders) > 1:
        for row in rows:
            if row["quadrature_order"] == args.orders[-1]:
                continue
            ref = next(r for r in rows if r["field"] == row["field"] and r["profile"] == row["profile"]
                       and r["refinement"] == row["refinement"] and r["quadrature_order"] == args.orders[-1])
            for kind in ("scattering", "recombination"):
                actual = np.array(row["full_operator"][kind]["electronic_number_rhs"])
                expected = np.array(ref["full_operator"][kind]["electronic_number_rhs"])
                convergence.append({"field": row["field"], "profile": row["profile"], "channel": kind,
                                    "refinement": row["refinement"], "order": row["quadrature_order"],
                                    "reference_order": args.orders[-1],
                                    "electronic_weighted_L1_relative_error": relative_error(actual, expected, expected)})
    result = {"schema": "pysnspd.stage2.projected_event_diagnostics.v1", "host": platform.node(),
              "criteria_sha256": sha(CRITERIA), "catalog_sha256": sha(CATALOG), "reference_sha256": sha(args.reference),
              "projection_reference_sha256": sha(args.projection_reference),
              "count_refinement_code_sha256": sha(ROOT / "pysnspd/experimental/refined_cells.py"),
              "event_code_sha256": sha(ROOT / "pysnspd/experimental/kinetic_events.py"), "runner_sha256": sha(__file__),
              "rows": rows, "runtime_seconds": time.perf_counter()-started,
              "quadrature_self_convergence": convergence,
              "method": args.method,
              "measurement_only": args.measurement_only,
              "status": "PENDING_CONTINUOUS_AND_DYNAMICAL_ADJUDICATION"}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    worst = max(((max(r["errors_against_independent_eta_zero"]["full_operator"][kind]["relative_observables"][:2]),
                 r["field"], r["profile"], kind) for r in rows if r["quadrature_order"] == args.orders[-1]
                 and "errors_against_independent_eta_zero" in r
                 for kind in ("scattering", "recombination")), default=None)
    print(json.dumps({"runtime_seconds": result["runtime_seconds"], "worst_power_error": worst,
                      "worst_electronic_quadrature_relative_L1": max((r["electronic_weighted_L1_relative_error"] for r in convergence), default=None),
                      "output": str(args.output)}))


if __name__ == "__main__":
    main()
