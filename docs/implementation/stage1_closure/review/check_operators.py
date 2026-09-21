"""Independent instantaneous checks of the reduced-cell operators.

This does not call the implementer's diagnostic runner. The reference FD
inversion uses bisection and the heat reference uses a direct normalized
formula where it is well conditioned. These test numerical contracts, not
microscopic accuracy of the effective BGK/heating laws.
"""
from pathlib import Path
import hashlib
import json
import sys
import time

import numpy as np
from scipy.special import expit

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
from pysnspd.experimental.energy_catalog import OccupationEnergyCatalog
from pysnspd.experimental.cell_validation import ElectronicCell


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def reference_temperature(energy, weights, p):
    target = 4*np.dot(weights*energy, p)
    if target == 0:
        return 0.
    lo, hi = 0., 1.
    while 4*np.dot(weights*energy, expit(-energy/hi)) < target:
        hi *= 2
    for _ in range(120):
        mid = (lo+hi)/2
        if 4*np.dot(weights*energy, expit(-energy/mid)) < target:
            lo = mid
        else:
            hi = mid
    return (lo+hi)/2


def main():
    started = time.perf_counter()
    catalogue_path = ROOT/'docs/implementation/stage1_r2/catalogs/occupation_catalog.npz'
    catalog = OccupationEnergyCatalog.load(catalogue_path)
    x, weights = catalog.count_nodes, catalog.count_weights
    rows = []
    fields = [(1., 0.), (.72, .002), (.72, .2), (.35, .65), (.16, 1.1), (1.45, 1.15)]
    for a, g in fields:
        cell = ElectronicCell(catalog, a, g)
        energy = cell.energies
        for name, p in [
            ('vacuum', np.zeros_like(x)),
            ('low', .2*np.exp(-x/.25)),
            ('band', .16*np.exp(-((x-.95)/.3)**2)),
            ('edge', .2*np.exp(-((x-.035)/.01)**2)),
            ('thermal', expit(-energy/.25))]:
            reference_t = reference_temperature(energy, weights, p)
            rhs, obtained_t = cell.bgk(p, .7)
            expected_rhs = (np.zeros_like(p) if reference_t == 0 else expit(-energy/reference_t))-p
            expected_rhs /= .7
            available = (p > 0) & (p < 1)
            entropy_rate = float(-4*np.dot(weights[available]*np.log(p[available]/(1-p[available])), rhs[available]))
            power, bath = .03, .08
            heat = cell.heating(p, power, bath)
            direct_shape = energy*expit(-energy/max(reference_t,bath))*(1-p)
            direct_heat = power*direct_shape/(4*np.dot(weights*energy, direct_shape))
            row = dict(amplitude=a, gamma=g, profile=name,
                       temperature_abs_error=abs(obtained_t-reference_t),
                       bgk_reference_abs_error=float(np.max(abs(rhs-expected_rhs))),
                       bgk_energy_moment=abs(float(4*np.dot(weights*energy,rhs))),
                       entropy_rate=entropy_rate,
                       source_reference_abs_error=float(np.max(abs(heat-direct_heat))),
                       source_energy_moment_error=abs(float(4*np.dot(weights*energy,heat))-power),
                       source_minimum=float(np.min(heat)))
            # dU/dt = X*a_dot + <E,p_dot>; same U supplies X, but
            # independently assembled scalar cancellation detects wrong factors.
            _, force, _ = cell.moments(p)
            mobility=.13
            adot=-mobility*force
            dissipated=mobility*force*force
            deposited=cell.heating(p,dissipated,bath)
            row['self_consistent_instantaneous_balance']=abs(float(force*adot+4*np.dot(weights*energy,deposited+rhs)))
            row['dissipation']=dissipated
            rows.append(row)

    cell = ElectronicCell(catalog,.72,.2)
    p = np.zeros_like(x)
    # Block one occupied low-energy state while other states remain available.
    p[0] = 1
    filled_source = cell.heating(p,.03,.08)
    cold_source = cell.heating(np.zeros_like(x),.03,1e-5)
    failures={}
    for name, operation in [
        ('tau_zero',lambda:cell.bgk(np.zeros_like(x),0)),
        ('tau_nan',lambda:cell.bgk(np.zeros_like(x),float('nan'))),
        ('saturated_heat',lambda:cell.heating(np.ones_like(x),.03,.08)),
        ('energy_above_FD_range',lambda:cell.equivalent_temperature(np.full_like(x,.8))),
        ('negative_power',lambda:cell.heating(np.zeros_like(x),-.03,.08)),
        ('zero_bath',lambda:cell.heating(np.zeros_like(x),.03,0))]:
        try:
            operation()
        except (ValueError,FloatingPointError):
            failures[name]=True
        else:
            failures[name]=False
    summary={key:float(max(abs(row[key]) for row in rows)) for key in
             ('temperature_abs_error','bgk_reference_abs_error','bgk_energy_moment',
              'source_reference_abs_error','source_energy_moment_error','self_consistent_instantaneous_balance')}
    summary.update(entropy_rate_minimum=float(min(row['entropy_rate'] for row in rows)),
                   source_minimum=float(min(row['source_minimum'] for row in rows)),
                   filled_state_source=float(filled_source[0]),
                   cold_source_energy_moment_error=abs(float(4*np.dot(cell.weights*cell.energies,cold_source))-.03))
    passed=(summary['bgk_reference_abs_error']<=1e-10 and summary['bgk_energy_moment']<=1e-11
            and summary['source_reference_abs_error']<=1e-10 and summary['source_energy_moment_error']<=1e-11
            and summary['self_consistent_instantaneous_balance']<=1e-10
            and summary['entropy_rate_minimum']>=-1e-11 and summary['source_minimum']>=0
            and summary['filled_state_source']==0 and summary['cold_source_energy_moment_error']<=1e-11
            and all(failures.values()))
    output=dict(schema='pysnspd.stage1_closure.independent_operators.v1',
                status='PASS' if passed else 'FAIL',
                scope='Instantaneous electronic operators; no trajectory, transport or phonon validation inferred.',
                criteria_sha256=sha(ROOT/'docs/implementation/stage1_closure/acceptance_criteria.json'),
                catalog_sha256=sha(catalogue_path),
                query_source_sha256=sha(ROOT/'pysnspd/experimental/energy_catalog.py'),
                cell_source_sha256=sha(ROOT/'pysnspd/experimental/cell_validation.py'),
                reviewer_source_sha256=sha(Path(__file__)),rows=rows,summary=summary,
                rejected_invalid_inputs=failures,runtime_seconds=time.perf_counter()-started)
    destination=Path(__file__).with_name('operator_checks.json')
    destination.write_text(json.dumps(output,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({key:output[key] for key in ('status','summary','runtime_seconds')},indent=2))


if __name__=='__main__':
    main()
