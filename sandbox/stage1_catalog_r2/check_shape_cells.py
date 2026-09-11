"""Check cubic energy ordering over every ratio interval at sampled amplitudes.

Bernstein controls give a sufficient bound; stationary points of the same cubic
give its minimum without a dense ratio grid. This is not a proof uniform in the
amplitude coordinate. No spectrum is solved and no table data are repaired.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import numpy as np
from pysnspd.experimental.energy_catalog import OccupationEnergyCatalog


def cubic_minimum(controls, upper):
    b0, b1, b2, b3 = controls
    a = -b0 + 3*b1 - 3*b2 + b3
    b = 3*b0 - 6*b1 + 3*b2
    c = -3*b0 + 3*b1
    def value(t):
        return ((a*t + b)*t + c)*t + b0
    minimum = np.minimum(b0, value(upper))
    # Roots of 3*a*t**2 + 2*b*t + c; the stable quadratic formula avoids
    # cancellation when one stationary point is close to a cell boundary.
    disc = b*b - 3*a*c
    root = np.sqrt(np.maximum(disc, 0))
    q = -b - np.copysign(root, b)
    with np.errstate(divide='ignore', invalid='ignore'):
        candidates = (q/(3*a), c/q, -c/(2*b))
    for index, t in enumerate(candidates):
        valid = (t > 0) & (t < upper) & np.isfinite(t)
        valid &= ((a == 0) & (b != 0)) if index == 2 else ((a != 0) & (disc >= 0))
        safe_t = np.where(valid, t, 0.)
        minimum = np.minimum(minimum, np.where(valid, value(safe_t), np.inf))
    return minimum


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('catalog', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    started = time.perf_counter()
    cat = OccupationEnergyCatalog.load(args.catalog)
    axes, ratios = cat.vacuum.delta_axis, cat.gamma_ratio_axis
    y = np.log1p(ratios/cat.gamma_coordinate_scale)
    h = np.diff(y)
    fractions = np.array([.21132486540518713, .5, .7886751345948129])
    deltas = np.sort(np.r_[axes, (axes[:-1,None]+np.diff(axes)[:,None]*fractions).ravel()])
    rows, flagged, crossing = [], set(), set()
    endpoints_bad = cells = 0
    global_min = np.inf
    for delta in deltas:
        i = min(np.searchsorted(axes, delta, side='right')-1, len(axes)-2)
        offset = delta-axes[i]
        def poly(coefficients):
            v = coefficients[:,i,:,:]
            return ((v[0]*offset+v[1])*offset+v[2])*offset+v[3]
        base = cat.count_nodes*np.sqrt(1+delta*delta/(cat.count_nodes**2+cat.eta**2))
        energy = base+poly(cat._log_coefficients)
        slope = poly(cat._gamma_log_coefficients)
        de, dm = np.diff(energy,axis=1), np.diff(slope,axis=1)
        endpoints_bad += int(np.sum(de[ratios*delta <= 1.2] <= 0))
        b0, b3 = de[:-1], de[1:]
        controls = np.stack((b0, b0+h[:,None]*dm[:-1]/3,
                             b3-h[:,None]*dm[1:]/3, b3))
        for j in np.flatnonzero(ratios[:-1]*delta < 1.2):
            upper = min(1., (np.log1p(1.2/delta/cat.gamma_coordinate_scale)-y[j])/h[j])
            actual = float(cubic_minimum(controls[:,j], upper).min())
            lower = float(controls[:,j].min())
            global_min = min(global_min, actual)
            cells += 1
            if lower <= 0:
                flagged.add(int(j))
                rows.append(dict(delta=float(delta), index=int(j), r_low=float(ratios[j]),
                                 r_high=float(ratios[j+1]), accessible_fraction=upper,
                                 minimum_Bernstein_control_increment=lower,
                                 minimum_cubic_increment=actual))
            if actual <= 0:
                crossing.add(int(j))
    def midpoints(indices):
        return [float(cat.gamma_coordinate_scale*np.expm1((y[j]+y[j+1])/2)) for j in sorted(indices)]
    digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    record = dict(schema='pysnspd.stage1_r2.shape_cells.v1',
        status='PASS' if not crossing and endpoints_bad == 0 else 'FAIL',
        method='Minimum of each exact Hermite energy-difference cubic over the physical ratio interval, evaluated at nodal amplitudes and three interior amplitude fractions. Bernstein bounds are also reported; a negative control alone is not a failed query.',
        coverage_limit='Amplitude is sampled, not certified uniformly; strictly positive energies are checked independently by the query assessment.',
        runtime_seconds=time.perf_counter()-started, catalog=str(args.catalog),
        catalog_sha256=digest(args.catalog), source_sha256=digest(ROOT/'pysnspd/experimental/energy_catalog.py'),
        script_sha256=digest(Path(__file__)), amplitudes=len(deltas), amplitude_fractions=fractions.tolist(),
        checked_cells=cells, endpoint_nonpositive_increments=endpoints_bad,
        minimum_cubic_increment=global_min, flagged_Bernstein_intervals=len(flagged),
        crossing_intervals=len(crossing), ratio_midpoints_proposed=midpoints(crossing),
        Bernstein_midpoints=midpoints(flagged), rows=rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record,indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k:v for k,v in record.items() if k not in ('rows','method','coverage_limit')},indent=2))


if __name__ == '__main__':
    main()
