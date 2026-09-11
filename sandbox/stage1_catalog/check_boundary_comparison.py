"""Compare the final Hermite and previous spline at identical stored nodes.

This lightweight diagnostic uses the saved occupation pilot, never a transient.
The previous interpolant is reconstructed on final nodal energies to isolate
the interpolation change from spectral-solver roundoff. The direct reference
uses the same finite eta, count nodes, quadrature weights and fixed occupation.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import platform
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import numpy as np
from scipy.interpolate import RectBivariateSpline

from pysnspd.experimental.energy_catalog import (
    OccupationEnergyCatalog,
    energy_at_count,
    retarded_spectrum,
    vacuum_state,
)


def main() -> None:
    started = time.perf_counter()
    npz = ROOT / "docs/implementation/stage1/catalogs/occupation_catalog.npz"
    catalog = OccupationEnergyCatalog.load(npz)
    source = ROOT / "pysnspd/experimental/energy_catalog.py"
    previous_source = ROOT / "docs/implementation/stage1/history/before_gamma_hermite/energy_catalog.py"
    increments = np.diff(
        catalog.excitation_energies,
        axis=-1,
        prepend=np.zeros((*catalog.excitation_energies.shape[:-1], 1)),
    )
    previous = [
        RectBivariateSpline(
            catalog.vacuum.delta_axis,
            catalog.vacuum.gamma_axis,
            np.log(increments[:, :, k]),
            s=0,
        )
        for k in range(len(catalog.count_nodes))
    ]
    occupation = .2 * np.exp(-catalog.count_nodes / .25)
    rows = []
    for delta, gamma in ((.17, 0), (.72, 0), (1.2, .0001), (.72, .002), (.43, .02)):
        energy = energy_at_count(catalog.count_nodes, delta=delta, gamma=gamma, eta=catalog.eta)
        c, s = retarded_spectrum(energy, delta=delta, gamma=gamma, eta=catalog.eta)
        reference = np.asarray(vacuum_state(delta, gamma)) + 4 * np.sum(
            np.array([energy, s.imag / c.real, -s.real * s.imag / c.real])
            * (occupation * catalog.count_weights)[None, :],
            axis=1,
        )
        current = np.asarray(catalog.evaluate(delta, gamma, occupation))
        previous_increments = np.exp([float(spline.ev(delta, gamma)) for spline in previous])
        previous_slope = np.cumsum(
            previous_increments * np.array([float(spline.ev(delta, gamma, dy=1)) for spline in previous])
        )
        previous_gamma = catalog.vacuum.evaluate(delta, gamma)[2] + 4 * np.dot(
            previous_slope * occupation, catalog.count_weights
        )
        rows.append({
            "delta": delta,
            "gamma": gamma,
            "eta": catalog.eta,
            "current_slope_reference_same_eta_same_quadrature": float(reference[2]),
            "new_gamma_conjugate": float(current[2]),
            "new_signed_difference": float(current[2] - reference[2]),
            "new_relative_error": float(abs((current[2] - reference[2]) / reference[2])),
            "previous_interpolant_same_final_nodes": float(previous_gamma),
            "previous_relative_error": float(abs((previous_gamma - reference[2]) / reference[2])),
        })
    report = {
        "host": platform.node(),
        "runtime_seconds": time.perf_counter() - started,
        "occupation": "0.2 exp(-x/0.25)",
        "normalization": "conjugate to Gamma divided by N0 Delta0; proportional to j/q for fixed scales",
        "scope": (
            "independent direct causal reference versus final Hermite interpolation; "
            "previous interpolation reconstructed on identical final nodal energies"
        ),
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "catalog_sha256": hashlib.sha256(npz.read_bytes()).hexdigest(),
        "previous_interpolant_source_sha256": hashlib.sha256(previous_source.read_bytes()).hexdigest(),
        "rows": rows,
        "caution": (
            "This isolates interpolation error. Eta=1e-3 still biases the ideal-limit slope "
            "by 8.61% in the delta=.72 reference; q=0 current itself is zero."
        ),
    }
    output = ROOT / "docs/implementation/stage1/electronic/boundary_comparison.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
