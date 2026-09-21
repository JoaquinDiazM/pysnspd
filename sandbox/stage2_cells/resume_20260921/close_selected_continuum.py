"""One bounded independent continuum check at the selected IR=.005 cut.

No event list or trajectory is evaluated. Existing measured native RHS values
are compared with B.7 partner integrals and independent (E,Omega) moments.
The earlier .01-cut references and their results remain unchanged.
"""
from pathlib import Path
import json
import platform
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
REVIEW = ROOT / "docs/implementation/stage2/review"
OUT = ROOT / "docs/implementation/stage2/resume_20260921"
sys.path[:0] = [str(ROOT), str(REVIEW)]
import continuous_reactions as continuous
import weak_population_reference as ideal
import weak_population_finite_eta as finite
from pysnspd.experimental.energy_catalog import OccupationEnergyCatalog
from pysnspd.experimental.refined_cells import refined_count_catalog

CUT = .005
UPPER = 4.
ORDER_PAIRS = ((2, 16), (4, 32), (8, 64))
PROFILE_NAMES = continuous.PROFILES
CHANNELS = ("scattering", "recombination")
FIELDS = ((0., 0.), (.72, 0.), (.35, .65))
CRITERIA = ROOT / "docs/implementation/stage2/acceptance_criteria.json"
MEASURED = ROOT / "docs/implementation/stage2/selected_event_measurements.json"
CATALOG = ROOT / "docs/implementation/stage1_r2/catalogs/occupation_catalog.npz"


def partner_rhs(coordinate, a, g, emax, order, eta=None):
    """B.7 at fixed energy; partner DOS removed analytically when Gamma=0.

    This is the previously audited independent formula with an explicit IR
    argument. Its integration coordinates are unrelated to the event mesh.
    With eta supplied, coordinate is x and both BCS ratios are analytic.
    """
    energy = coordinate if eta is None else finite.energy_at_x(coordinate, a, eta)
    z, w = np.polynomial.legendre.leggauss(order)
    result = np.zeros((3, 2, len(energy)))
    gap = a if g == 0 and eta is None else 0.
    if eta is not None:
        ratio = finite.ratio_at_x(coordinate, a, eta)
    elif g:
        ni, ri = ideal.gapless_values(energy, a, g)
        ratio = ri/ni
    breaks = np.array([0., .025, .05, .1, .2, .4, .65, .8, 1., 1.2,
                       1.5, 2., 2.5, 3., 4., 6., 9., 13.])
    if eta is not None:
        breaks = np.unique(np.r_[breaks, np.geomspace(1e-10, .01, 18)])
    for kind in ("upper", "lower", "recombination"):
        if kind == "upper":
            lower, upper = energy+CUT, np.minimum(emax, energy+UPPER)
        elif kind == "lower":
            lower, upper = np.maximum(gap, energy-UPPER), energy-CUT
        else:
            lower, upper = np.maximum(gap, CUT-energy), np.minimum(emax, UPPER-energy)
        valid = upper > lower
        if eta is not None:
            xmin = finite.count_at_energy(np.maximum(0., lower), a, eta)
            xmax = finite.count_at_energy(np.maximum(0., upper), a, eta)
        elif g == 0:
            xmin = np.sqrt(np.maximum(0., lower*lower-a*a))
            xmax = np.sqrt(np.maximum(0., upper*upper-a*a))
        else:
            xmin, xmax = lower, upper
        for lo, hi in zip(breaks[:-1], breaks[1:]):
            begin, end = np.maximum(xmin, lo), np.minimum(xmax, hi)
            ids = np.flatnonzero(valid & (end > begin))
            if not len(ids):
                continue
            x = begin[ids, None]+(end-begin)[ids, None]*(z[None, :]+1)/2
            weights = (end-begin)[ids, None]*w[None, :]/2
            ei = energy[ids, None]
            ej = finite.energy_at_x(x, a, eta) if eta is not None else np.hypot(a, x) if g == 0 else x
            omega = ej-ei if kind == "upper" else ei-ej if kind == "lower" else ei+ej
            sign_coherence = 1 if kind == "recombination" else -1
            if eta is not None:
                coherence = 1+sign_coherence*ratio[ids, None]*finite.ratio_at_x(x, a, eta)
            elif g:
                nj, rj = ideal.gapless_values(ej, a, g)
                coherence = nj+sign_coherence*ratio[ids, None]*rj
            else:
                coherence = 1+sign_coherence*a*a/(ei*ej)
            factor = .03*(omega/4.)**2*coherence*weights/4
            for ip, profile in enumerate(PROFILE_NAMES):
                pi, n = continuous.populations(profile, ei, omega)
                pj, _ = continuous.populations(profile, ej, omega)
                if kind == "upper":
                    net, sign = (1-pi)*pj*(1+n)-pi*(1-pj)*n, 1
                elif kind == "lower":
                    net, sign = (1-pj)*pi*(1+n)-pj*(1-pi)*n, -1
                else:
                    net, sign = pi*pj*(1+n)-(1-pi)*(1-pj)*n, -1
                result[ip, int(kind == "recombination"), ids] += sign*np.sum(factor*net, axis=1)
    return result


def weak_reference(counts, grid, a, g, outer_order, inner_order, eta=None):
    z, w = np.polynomial.legendre.leggauss(outer_order)
    if eta is not None:
        x = counts[:-1, None]+np.diff(counts)[:, None]*(z[None, :]+1)/2
        measure = np.diff(counts)[:, None]*w[None, :]/2
        energy = finite.energy_at_x(x, a, eta)
        ids = np.arange(len(counts)-1)
        coordinates = x
    else:
        low, high = np.maximum(grid[:-1], a if g == 0 else 0.), grid[1:]
        valid = high > low
        xl = np.sqrt(np.maximum(0., low[valid]**2-a*a)) if g == 0 else low[valid]
        xh = np.sqrt(high[valid]**2-a*a) if g == 0 else high[valid]
        x = xl[:, None]+(xh-xl)[:, None]*(z[None, :]+1)/2
        measure = (xh-xl)[:, None]*w[None, :]/2
        energy = np.hypot(a, x) if g == 0 else x
        if g:
            ni, _ = ideal.gapless_values(energy, a, g)
            measure *= ni
        ids = np.flatnonzero(valid)
        coordinates = energy
    beta = (energy-grid[ids, None])/(grid[ids+1, None]-grid[ids, None])
    rhs = partner_rhs(coordinates.ravel(), a, g, float(grid[-1]), inner_order, eta).reshape(3, 2, len(ids), outer_order)
    values = np.zeros((3, 2, len(grid)))
    values[:, :, ids] += 4*np.sum(rhs*measure[None, None, :, :]*(1-beta)[None, None, :, :], axis=-1)
    values[:, :, ids+1] += 4*np.sum(rhs*measure[None, None, :, :]*beta[None, None, :, :], axis=-1)
    return values


def project_33_to_17(values):
    """Exact nested-hat identity; no reinterpolation of a phonon density."""
    result = values[::2].copy()
    result[:-1] += .5*values[1::2]
    result[1:] += .5*values[1::2]
    return result


def main():
    started = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    measured = json.loads(MEASURED.read_text(encoding="utf-8"))
    criteria = json.loads(CRITERIA.read_text(encoding="utf-8"))
    tolerance = criteria["continuous_consistency"]["relative_error_max"]
    budget = criteria["continuous_consistency"]["reference_refinement_relative_budget_max"]
    if measured["event_code_sha256"] != continuous.sha(ROOT / "pysnspd/experimental/kinetic_events.py"):
        raise ValueError("measured operator does not match the frozen kernel")
    if measured["count_refinement_code_sha256"] != continuous.sha(ROOT / "pysnspd/experimental/refined_cells.py"):
        raise ValueError("measured count grid does not match the frozen factory")
    if measured["catalog_sha256"] != continuous.sha(CATALOG) or measured["criteria_sha256"] != continuous.sha(CRITERIA):
        raise ValueError("catalogue or acceptance criteria hash changed")
    selected = [r for r in measured["rows"] if r["quadrature_order"] == 2 and r["refinement"] == 1]
    if len(selected) != 9 or any(r["phonon_nodes"] != 1025 or r["metadata"]["coupling_support"] != [CUT, UPPER] for r in selected):
        raise ValueError("expected exactly the selected 630/1025 IR=.005 nine cases")
    catalog = refined_count_catalog(OccupationEnergyCatalog.load(CATALOG), 1)
    references, comparisons = [], []
    for a, g in FIELDS:
        grid = catalog.energy_kernel(a, g)[0]
        for eta in ((None, catalog.eta) if (a, g) == (.72, 0.) else (None,)):
            levels = []
            for outer, inner in ORDER_PAIRS:
                tick = time.perf_counter()
                levels.append(weak_reference(catalog.count_nodes, grid, a, g, outer, inner, eta))
                print(json.dumps({"kind": "weak", "field": [a, g], "eta": eta, "order": [outer, inner], "seconds": time.perf_counter()-tick}), flush=True)
            for ip, profile in enumerate(PROFILE_NAMES):
                observed = next(r for r in selected if r["amplitude"] == a and r["gamma"] == g and r["profile"] == profile)
                for ic, channel in enumerate(CHANNELS):
                    expected = levels[-1][ip, ic]
                    scale = float(np.sum(abs(expected)))
                    refinement_error = float(np.sum(abs(expected-levels[-2][ip, ic]))/scale)
                    reference_row = {"amplitude": a, "gamma": g, "eta": 0. if eta is None else eta,
                        "profile": profile, "channel": channel, "number_rhs_reference": expected.tolist(),
                        "reference_L1": scale, "last_refinement_relative_L1": refinement_error,
                        "energy_nodes": grid.tolist()}
                    references.append(reference_row)
                    for representation in ("point_activities", "native_pauli_point_bose", "full_operator"):
                        actual = np.asarray(observed[representation][channel]["electronic_number_rhs"])
                        difference = float(np.sum(abs(actual-expected))/scale)
                        comparisons.append({key: reference_row[key] for key in ("amplitude", "gamma", "eta", "profile", "channel")}
                            | {"representation": representation, "relative_L1": difference,
                               "reference_refinement_relative_L1": refinement_error,
                               "pass": difference <= tolerance and refinement_error <= budget})
    # The independent E/Omega reference integrates on every 33-hat boundary.
    # Its 17-hat values follow an exact nested-basis identity, including both
    # zeroth and first moments. No experimental event routine is imported.
    axis = np.linspace(0., 4., 33)
    original_segments, original_project = continuous.segments, continuous.project_hats
    def aligned_segments(bounds, order, **kwargs):
        cuts = sorted(set(list(bounds)+[float(v) for v in axis if bounds[0] < v < bounds[-1]]))
        return original_segments(cuts, order, **kwargs)
    continuous.segments = aligned_segments
    continuous.project_hats = lambda omega, rates, weights, unused: original_project(omega, rates, weights, axis)
    moment_rows = []
    try:
        for a, g in FIELDS:
            coarse = continuous.evaluate_field(a, g, 32, omega_min=CUT)
            fine = continuous.evaluate_field(a, g, 64, omega_min=CUT)
            for profile in PROFILE_NAMES:
                observed = next(r for r in selected if r["amplitude"] == a and r["gamma"] == g and r["profile"] == profile)
                for channel in CHANNELS:
                    expected = fine[0][profile][channel]
                    previous = coarse[0][profile][channel]
                    actual = np.asarray(observed["full_operator"][channel]["observables"])
                    reference_error = abs(expected-previous)/np.maximum(abs(expected), 1e-100)
                    error = abs(actual-expected)/np.maximum(abs(expected), 1e-100)
                    hats = {}
                    for size in (17, 33):
                        ref = fine[1][profile][channel]
                        old = coarse[1][profile][channel]
                        if size == 17:
                            ref, old = project_33_to_17(ref), project_33_to_17(old)
                        measured_hats = np.asarray(observed["full_operator"][channel][f"projected_{size}"])
                        scale = float(np.sum(abs(ref)))
                        err, referr = float(np.sum(abs(measured_hats-ref))/scale), float(np.sum(abs(ref-old))/scale)
                        hats[str(size)] = {"reference": ref.tolist(), "actual": measured_hats.tolist(),
                            "relative_L1": err, "reference_refinement_relative_L1": referr,
                            "pass": err <= tolerance and referr <= budget}
                    moment_rows.append({"amplitude": a, "gamma": g, "profile": profile, "channel": channel,
                        "reference": expected.tolist(), "actual": actual.tolist(),
                        "relative_errors": error.tolist(), "reference_refinement_relative_errors": reference_error.tolist(),
                        "power_and_count_pass": bool(np.all(error[:4] <= tolerance) and np.all(reference_error[:4] <= budget)),
                        "hats": hats})
    finally:
        continuous.segments, continuous.project_hats = original_segments, original_project
    passed = (all(r["pass"] for r in comparisons if r["representation"] == "full_operator")
              and all(r["power_and_count_pass"] and all(h["pass"] for h in r["hats"].values()) for r in moment_rows))
    data = {"schema": "pysnspd.stage2.selected_continuum_resume.v1", "host": platform.node(),
        "status": "PASS_STATIC_SELECTED_CONTINUUM" if passed else "FAIL_STATIC_SELECTED_CONTINUUM",
        "scope": "Existing selected630/1025 event measurements, IR=.005, three fields/profiles; no temporal or global-stage verdict",
        "eta_policy": "Ideal eta0 reference for all fields plus separate analytic BCS reference at eta1e-8; finite-eta bias is not quadrature error",
        "omega_support": [CUT, UPPER], "event_orders": [2, 2], "weak_reference_orders": ORDER_PAIRS,
        "moment_reference_orders": [32, 64], "rate_prefactor": 1., "alpha2F": ".03*(Omega/4)^2",
        "tolerance": tolerance, "reference_budget": budget,
        "count_nodes": catalog.count_nodes.tolist(), "count_weights": catalog.count_weights.tolist(),
        "weak_references": references, "weak_comparisons": comparisons, "moments_and_hats": moment_rows,
        "runtime_seconds": time.perf_counter()-started,
        "hashes": {"criteria": continuous.sha(CRITERIA), "measured": continuous.sha(MEASURED),
            "catalog": continuous.sha(CATALOG), "script": continuous.sha(Path(__file__)),
            "kernel": continuous.sha(ROOT / "pysnspd/experimental/kinetic_events.py"),
            "factory": continuous.sha(ROOT / "pysnspd/experimental/refined_cells.py"),
            **{name: continuous.sha(REVIEW / name) for name in ("continuous_reactions.py", "weak_population_reference.py", "weak_population_finite_eta.py")}}}
    target = OUT / "selected_continuum.json"
    target.write_text(json.dumps(data, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print(json.dumps({"status": data["status"], "runtime_seconds": data["runtime_seconds"],
        "weak_maximum": max(r["relative_L1"] for r in comparisons if r["representation"] == "full_operator"),
        "power_maximum": max(max(r["relative_errors"][:2]) for r in moment_rows),
        "hats_maximum": max(h["relative_L1"] for r in moment_rows for h in r["hats"].values())}), flush=True)


if __name__ == "__main__":
    main()
