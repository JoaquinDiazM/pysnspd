"""Bounded fresh diagnostics of electron--phonon event discretization.

No complete transient or NbN calibration is implied. The continuous-reference
comparison is kept separate from the algebraic conservation of each event.
"""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import numpy as np
from scipy.special import expit
from numpy.polynomial.legendre import leggauss
from pysnspd.experimental.cell_validation import ElectronicCell
from pysnspd.experimental.energy_catalog import OccupationEnergyCatalog, energy_at_count_batch, retarded_spectrum
from pysnspd.experimental.kinetic_events import ElectronPhononEvents, PhononGrid

OUT = ROOT / "docs/implementation/stage2"
CATALOG = ROOT / "docs/implementation/stage1_r2/catalogs/occupation_catalog.npz"
CRITERIA = OUT / "acceptance_criteria.json"
FIELDS = [("gapped", .72, 0.), ("depaired", .72, .2), ("gapless", .35, .65)]
PROFILES = ("hot_electrons", "phonon_bubble", "nonthermal")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def bose(omega, temperature):
    z = np.asarray(omega)/temperature
    result = np.zeros_like(z)
    active = z < 700
    result[active] = 1/np.expm1(z[active])
    return result


def populations(profile, energy, omega):
    if profile == "hot_electrons":
        return expit(-energy/.4), bose(omega, .15)
    if profile == "phonon_bubble":
        return expit(-energy/.12), bose(omega, .12)+.3*np.exp(-((omega-2)/.35)**2)
    if profile == "nonthermal":
        p = .12*np.exp(-energy/.25)+.1*np.exp(-((energy-1.3)/.25)**2)
        return p, bose(omega, .18)+.1*np.exp(-((omega-1.7)/.4)**2)
    raise ValueError(profile)


def phonon_grid(nodes, infrared=.01):
    energies = np.linspace(infrared, 4., nodes)
    weights = np.r_[np.diff(energies)[0]/2, (energies[2:]-energies[:-2])/2, np.diff(energies)[-1]/2]
    capacities = weights*90*energies**2/4**3
    return PhononGrid(energies, capacities, lambda omega: .03*(omega/4)**2,
                      (infrared, 4.), "synthetic Debye: 30 modes in N0*Delta0 units; explicit IR cutoff")


def point_rates(network, profile):
    """Direct D.14 at event Omega, without any phonon-grid interpolation."""
    p, n = populations(profile, network.electron_energies, network.omega)
    i, j, rec = network.index_i, network.index_j, network.recombination
    forward = np.where(rec, p[i]*p[j], (1-p[i])*p[j])*(1+n)
    reverse = np.where(rec, (1-p[i])*(1-p[j]), p[i]*(1-p[j]))*n
    return network.coefficients*(forward-reverse), network.coefficients*(forward+reverse)


def project_measure(energies, quantities, nodes=17):
    """Positive projection on fixed hat bands, preserving count and energy."""
    band_nodes = np.linspace(0., 4., nodes)
    upper = np.searchsorted(band_nodes, energies, side="left")
    lower = np.maximum(upper-1, 0)
    span = band_nodes[upper]-band_nodes[lower]
    beta = np.divide(energies-band_nodes[lower], span, out=np.zeros_like(energies), where=span > 0)
    return (np.bincount(lower, weights=(1-beta)*quantities, minlength=len(band_nodes))
            +np.bincount(upper, weights=beta*quantities, minlength=len(band_nodes)))


def relative_error(actual, reference, gross):
    denominator = max(float(np.sum(abs(reference))), 1e-12*float(np.sum(abs(gross))), 1e-300)
    return float(np.sum(abs(actual-reference))/denominator)


def power_summary(network, rates, gross):
    return {name: {"net": float(np.sum(network.omega[mask]*rates[mask])),
                   "gross": float(np.sum(network.omega[mask]*gross[mask]))}
            for name, mask in (("scattering", ~network.recombination), ("recombination", network.recombination))}


def algebraic_checks(catalog):
    rows = []
    rng = np.random.default_rng(482103)
    for label, delta, gamma in [("normal", 0., 0.), *FIELDS]:
        network = ElectronPhononEvents.from_cell(ElectronicCell(catalog, delta, gamma), phonon_grid(129), rate_prefactor=1.)
        recovered = network.beta_lower*network.phonons.energies[network.phonon_lower]
        recovered += network.beta_upper*network.phonons.energies[network.phonon_upper]
        p = rng.uniform(.01, .8, len(network.electron_energies))
        n = rng.uniform(.01, 2, len(network.phonons.energies))
        dp, dn = network.rhs(p, n)
        e_power = float(np.dot(network.electron_capacities*network.electron_energies, dp))
        ph_power = float(np.dot(network.phonons.capacities*network.phonons.energies, dn))
        row = {"field": label, "amplitude": delta, "gamma": gamma,
               "event_energy_relative_error": float(np.max(abs(recovered-network.omega)/network.omega)),
               "rhs_energy_normalized_residual": abs(e_power+ph_power)/max(abs(e_power)+abs(ph_power), 1e-300),
               "metadata": network.metadata, "thermal_equilibria": []}
        logf, logr = network.log_activities(p, n)
        row["interior_collision_entropy_production"] = float(np.dot(network.rates(p, n), logf-logr))
        row["negative_control_reversed_phonon_energy_sign_residual"] = abs(e_power-ph_power)/max(1., abs(e_power)+abs(ph_power))
        for temperature in (.08, .2, .5):
            peq, neq = expit(-network.electron_energies/temperature), bose(network.phonons.energies, temperature)
            logf, logr = network.log_activities(peq, neq)
            gross = network.coefficients*np.exp(np.maximum(logf, logr))
            rates = network.rates(peq, neq)
            representable = gross > 1e-100
            dp_eq, dn_eq = network.rhs(peq, neq)
            row["thermal_equilibria"].append({"temperature": temperature,
                "maximum_event_relative_residual": float(np.max(abs(rates[representable])/gross[representable])),
                "weighted_absolute_residual": float(np.sum(abs(rates))),
                "weighted_equilibrium_RHS_absolute": float(np.dot(network.electron_capacities, abs(dp_eq))
                    +np.dot(network.phonons.capacities, abs(dn_eq))),
                "represented_events": int(np.count_nonzero(representable))})
        p = rng.choice([0., .2, 1.], len(network.electron_energies))
        n = rng.choice([0., .3], len(network.phonons.energies))
        dp, dn = network.rhs(p, n)
        row["physical_faces"] = {"electron_zero_min_rhs": float(np.min(dp[p == 0])),
                                  "electron_one_max_rhs": float(np.max(dp[p == 1])),
                                  "phonon_zero_min_rhs": float(np.min(dn[n == 0]))}
        rows.append(row)
    return rows


def phonon_convergence(catalog, resolutions):
    rows = []
    for label, delta, gamma in FIELDS:
        cell = ElectronicCell(catalog, delta, gamma)
        for nodes in resolutions:
            network = ElectronPhononEvents.from_cell(cell, phonon_grid(nodes), rate_prefactor=1.)
            for profile in PROFILES:
                p, n = populations(profile, network.electron_energies, network.phonons.energies)
                actual = network.rates(p, n)
                reference, gross = point_rates(network, profile)
                dp, dn = network.rhs_from_rates(actual)
                dp_ref, _ = network.rhs_from_rates(reference)
                band_actual = project_measure(network.phonons.energies, network.phonons.capacities*dn)
                band_reference = project_measure(network.omega, reference)
                band_gross = project_measure(network.omega, gross)
                finer_actual = project_measure(network.phonons.energies, network.phonons.capacities*dn, 33)
                finer_reference = project_measure(network.omega, reference, 33)
                finer_gross = project_measure(network.omega, gross, 33)
                row = {"field": label, "profile": profile, "phonon_nodes": nodes,
                       "electronic_weighted_L1_relative_error": relative_error(network.electron_capacities*dp,
                           network.electron_capacities*dp_ref, np.array([2*np.sum(gross)])),
                       "phonon_projected_L1_relative_error": relative_error(band_actual, band_reference, band_gross),
                       "phonon_projected_33_L1_relative_error": relative_error(finer_actual, finer_reference, finer_gross),
                       "point_powers": power_summary(network, reference, gross),
                       "interpolated_powers": power_summary(network, actual, gross)}
                for kind in ("scattering", "recombination"):
                    r = row["point_powers"][kind]
                    row[kind+"_power_relative_error"] = relative_error(
                        np.array([row["interpolated_powers"][kind]["net"]]), np.array([r["net"]]), np.array([r["gross"]]))
                rows.append(row)
    return rows


def electronic_quadrature(orders):
    """Independent x grids; point Bose factors isolate electronic quadrature."""
    rows = []
    for label, delta, gamma in [("normal", 0., 0.), FIELDS[0], FIELDS[2]]:
        for order in orders:
            node, weight = leggauss(order)
            counts, weights = 6*(node+1), 6*weight
            energies = energy_at_count_batch(counts, delta=delta, gamma=gamma, eta=1e-8)
            c, s = retarded_spectrum(energies, delta=delta, gamma=gamma, eta=1e-8)
            network = ElectronPhononEvents.from_arrays(energies, weights, s.imag/c.real,
                phonon_grid(65), rate_prefactor=1., label=f"direct causal x Gauss-Legendre order {order}")
            for profile in PROFILES:
                rates, gross = point_rates(network, profile)
                rows.append({"field": label, "amplitude": delta, "gamma": gamma, "eta": 1e-8,
                             "profile": profile, "count_order": order, "count_cutoff": 12.,
                             "powers": power_summary(network, rates, gross)})
    return rows


def native_electronic_quadrature(catalog):
    rows = []
    for label, delta, gamma in [("normal", 0., 0.), *FIELDS]:
        network = ElectronPhononEvents.from_cell(ElectronicCell(catalog, delta, gamma),
                                                 phonon_grid(65), rate_prefactor=1.)
        for profile in PROFILES:
            rates, gross = point_rates(network, profile)
            rows.append({"field": label, "amplitude": delta, "gamma": gamma, "eta": catalog.eta,
                         "profile": profile, "count_nodes": len(catalog.count_nodes),
                         "last_count_node": float(catalog.count_nodes[-1]),
                         "powers": power_summary(network, rates, gross)})
    return rows


def cutoff_checks(catalog):
    rows = []
    for label, delta, gamma in FIELDS:
        cell = ElectronicCell(catalog, delta, gamma)
        for infrared in (0., .02, .01, .005):
            # A positive marker below every existing event is sufficient here:
            # the rate is evaluated at Omega directly, not at the grid markers.
            grid = PhononGrid(np.array([1e-16, 4.]), np.ones(2), lambda omega: .03*(omega/4)**2,
                              (infrared, 4.), "point-rate-only compact cutoff diagnostic")
            network = ElectronPhononEvents.from_cell(cell, grid, rate_prefactor=1.)
            for profile in PROFILES:
                rates, gross = point_rates(network, profile)
                rows.append({"field": label, "profile": profile, "infrared_cutoff": infrared,
                             "powers": power_summary(network, rates, gross),
                             "metadata": network.metadata})
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUT / "event_diagnostics.json")
    parser.add_argument("--phonon-nodes", type=int, nargs="+", default=[65, 129, 257])
    parser.add_argument("--electron-orders", type=int, nargs="+", default=[24, 48, 96])
    args = parser.parse_args()
    if not CRITERIA.exists():
        raise SystemExit("Freeze stage2 acceptance_criteria.json before producing numerical results")
    started = time.perf_counter()
    catalog = OccupationEnergyCatalog.load(CATALOG)
    result = {"schema": "pysnspd.stage2_cells.event_diagnostics.v1", "host": platform.node(),
              "python": platform.python_version(), "criteria_sha256": sha(CRITERIA),
              "catalog_sha256": sha(CATALOG), "runner_sha256": sha(__file__),
              "event_code_sha256": sha(ROOT / "pysnspd/experimental/kinetic_events.py"),
              "coupling": "alpha2F=.03*(Omega/4)^2 on the declared [IR,4] support; zero outside",
              "rate_prefactor": 1., "scale_status": "synthetic normalization, no absolute material time",
              "projection": "17 fixed hat-band nodes 0,.25,...4; positive projection preserves count and first energy moment",
              "algebraic": algebraic_checks(catalog),
              "phonon_convergence": phonon_convergence(catalog, args.phonon_nodes),
              "electronic_quadrature": electronic_quadrature(args.electron_orders),
              "native_electronic_quadrature": native_electronic_quadrature(catalog),
              "infrared_cutoff": cutoff_checks(catalog),
              "status": "PENDING_INDEPENDENT_CONTINUOUS_REFERENCE_AND_DYNAMICAL_ADJUDICATION"}
    result["runtime_seconds"] = time.perf_counter()-started
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print(json.dumps({"runtime_seconds": result["runtime_seconds"], "output": str(args.output),
                      "maximum_phonon_L1_error": max(r["phonon_projected_L1_relative_error"] for r in result["phonon_convergence"]
                                                      if r["phonon_nodes"] == args.phonon_nodes[-1])}))


if __name__ == "__main__":
    main()
