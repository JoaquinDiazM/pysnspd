"""Retarded spectral oracle on the completed thermal core, with matched exterior.

Each energy/eta/case continues independently from the lowest positive Matsubara
frequency. The R=12 prescribed-tanh exterior BVP is continued with it; neither
gap nor spectral boundary is replaced by a homogeneous-contact approximation.
No physical time, population update, circuit or photon is simulated here.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import multiprocessing as mp
import os
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
from parallel_runtime import initialize_affinity, limit_thread_environment, linux_resources, resource_budget
limit_thread_environment()
sys.path.insert(0, str(ROOT))
import numpy as np
from scipy.integrate import solve_bvp
from scipy.interpolate import CubicHermiteSpline
from pysnspd.experimental import thermal_spatial_usadel as thermal
from pysnspd.experimental import retarded_spatial_usadel as retarded

DEFAULT_PLAN = ROOT/'docs/implementation/stage4/self_consistent_review_20260924/next_retarded_plan.json'
RADIAL = ROOT/'docs/implementation/stage4/followup_20260923/raw/stage4_radial_reference_20260923'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    with Path(path).open('x', encoding='utf8') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')


def outer_angle(z, gap, radius, previous, tolerance=2e-12):
    """Continue the same depairing boundary root; complex Newton, no clipping."""
    d, gamma, theta = gap*np.tanh(radius), 1/radius**2, complex(previous)
    def equation(value):
        return z*np.sin(value)-d*np.cos(value)+gamma*np.sin(value)*np.cos(value)
    for iteration in range(40):
        residual = equation(theta)
        if abs(residual) <= tolerance:
            return theta
        derivative = z*np.cos(theta)+d*np.sin(theta)+gamma*np.cos(2*theta)
        step = -residual/derivative
        for reduction in range(24):
            trial = theta+2.**(-reduction)*step
            if abs(equation(trial)) < abs(residual):
                theta = trial
                break
        else:
            raise RuntimeError('Exterior depairing root line search failed')
    raise RuntimeError('Exterior depairing root iteration budget exhausted')


def radial_step(z, gap, mesh, values, previous_outer, *, tolerance, max_nodes):
    """Same radial BVP as the frozen thermal reference, analytically continued.

    Split real and imaginary components explicitly so a real-valued BVP solver
    uses the full Cauchy-Riemann Jacobian. No complex conjugate enters the ODE.
    """
    inner, radius = float(mesh[0]), float(mesh[-1])
    outer = outer_angle(z, gap, radius, previous_outer)
    guess = np.vstack((values[0].real, values[0].imag, values[1].real, values[1].imag))
    def ode(x, y):
        theta, derivative = y[0]+1j*y[1], y[2]+1j*y[3]
        second = -derivative/x+np.sin(theta)*np.cos(theta)/(x*x)+z*np.sin(theta)-gap*np.tanh(x)*np.cos(theta)
        return np.vstack((derivative.real, derivative.imag, second.real, second.imag))
    def jacobian(x, y):
        theta = y[0]+1j*y[1]
        coefficient = np.cos(2*theta)/(x*x)+z*np.cos(theta)+gap*np.tanh(x)*np.sin(theta)
        result = np.zeros((4, 4, len(x)))
        result[0, 2] = 1.; result[1, 3] = 1.
        result[2, 0] = coefficient.real; result[2, 1] = -coefficient.imag
        result[3, 0] = coefficient.imag; result[3, 1] = coefficient.real
        result[2, 2] = -1/x; result[3, 3] = -1/x
        return result
    def boundary(left, right):
        return np.array([left[2]-left[0]/inner, left[3]-left[1]/inner,
                         right[0]-outer.real, right[1]-outer.imag])
    def boundary_jacobian(left, right):
        lhs, rhs = np.zeros((4, 4)), np.zeros((4, 4))
        lhs[0, 0], lhs[0, 2] = -1/inner, 1.
        lhs[1, 1], lhs[1, 3] = -1/inner, 1.
        rhs[2, 0], rhs[3, 1] = 1., 1.
        return lhs, rhs
    result = solve_bvp(ode, boundary, mesh, guess, fun_jac=jacobian, bc_jac=boundary_jacobian,
                       tol=tolerance, max_nodes=max_nodes)
    if not result.success:
        raise RuntimeError('Analytic exterior BVP failed: '+result.message)
    values = np.vstack((result.y[0]+1j*result.y[1], result.y[2]+1j*result.y[3]))
    boundary_error = float(np.max(abs(boundary(result.y[:, 0], result.y[:, -1]))))
    if np.min(np.cos(values[0]).real) < -1e-7:
        raise RuntimeError('Exterior solution has negative DOS; no branch repair')
    return result.x, values, outer, dict(maximum_rms_residual=float(np.max(result.rms_residuals)),
        boundary_residual=boundary_error, nodes=len(result.x))


def job_id(case, eta, energy):
    return f"{case['id']}_eta{eta:.4f}_E{energy:.4f}".replace('.', 'p')


def worker(job):
    started = time.monotonic()
    plan, case = job['plan'], job['case']
    with np.load(ROOT/case['fields_path']) as data:
        d, xy = data['d'], data['coordinates_bar']
        graph = thermal.ThermalGraph(data['area_weights'], data['edges'], data['conductance'], xy, data['boundary_nodes'])
        f, g = data['f_lowest'], data['g_lowest']
    radius = np.linalg.norm(xy, axis=1)
    direction = np.divide(xy[:, 0]+1j*xy[:, 1], radius, out=np.zeros(len(d), complex), where=radius>0)
    a, b = f/(1+g), np.conj(f)/(1+g)
    gap, anchor = plan['gap_reference_kBTc'], np.pi*plan['T_K']/plan['Tc_K']
    initial = retarded.evaluate(graph, d, anchor, a, b)
    thermal_value = thermal.spectral_energy_gradient(graph, d, anchor, f/g)
    bridge = dict(action_absolute_difference=float(abs(initial.action-thermal_value.energy)),
        current_maximum_difference=float(np.max(abs(initial.link_derivative_alpha-thermal_value.link_derivative_alpha))),
        g_maximum_difference=float(np.max(abs(initial.g-g))),
        f_maximum_difference=float(np.max(abs(initial.f-f))))
    if max(bridge.values()) > plan['matsubara_bridge_tolerance']:
        raise RuntimeError('Exact thermal-to-stereographic coordinate bridge failed')
    with np.load(RADIAL/'tasks/R12_n0000.npz') as reference:
        mesh = reference['radius']
        values = np.vstack((reference['theta'], reference['derivative'])).astype(complex)
        outer = values[0, -1]
    target = gap*(job['eta_relative']-1j*job['energy_relative'])
    path = []
    for index, fraction in enumerate(np.linspace(0, 1, plan['continuation_steps']+1)[1:], start=1):
        z = anchor+(target-anchor)*fraction
        mesh, values, outer, exterior = radial_step(z, gap, mesh, values, outer,
            tolerance=plan['radial_tolerance'], max_nodes=plan['radial_max_nodes'])
        interpolation = CubicHermiteSpline(mesh, values[0], values[1])
        boundary = graph.boundary_nodes
        theta = interpolation(radius[boundary])
        spectral = np.tan(theta/2)
        solution = retarded.solve(graph, d, z, initial_a=a, initial_b=b,
            fixed_nodes=boundary, fixed_a=spectral*direction[boundary],
            fixed_b=spectral*np.conj(direction[boundary]), tolerance=plan['spectral_tolerance'],
            max_iterations=plan['maximum_newton_iterations'], causal_tolerance=plan['causal_tolerance'])
        a, b = solution.a, solution.b
        path.append(dict(step=index, z_real=float(z.real), z_imag=float(z.imag),
            root_residual=solution.residual, root_iterations=solution.iterations,
            normalization_residual=solution.normalization_residual, minimum_DOS=solution.minimum_DOS,
            maximum_backtracks=max(solution.backtracks, default=0), exterior=exterior))
        if index == 1 or index % 6 == 0 or index == plan['continuation_steps']:
            elapsed = time.monotonic()-started
            print(json.dumps(dict(event='SPECTRAL_PATH', job_id=job['id'], step=index,
                total=plan['continuation_steps'], elapsed_seconds=elapsed,
                eta_seconds=elapsed*(plan['continuation_steps']-index)/index,
                minimum_DOS=solution.minimum_DOS, root_residual=solution.residual)), flush=True)
    final = retarded.evaluate(graph, d, target, a, b)
    core = radius <= plan['core_radius_ell0']
    center = int(np.argmin(radius))
    normal = graph.area_weights[core].sum()
    metadata = dict(id=job['id'], case_id=case['id'], energy_relative=job['energy_relative'],
        eta_relative=job['eta_relative'], z_real=float(target.real), z_imag=float(target.imag),
        node_count=graph.n_nodes, worker_seconds=time.monotonic()-started, worker_pid=os.getpid(),
        center_DOS=float(final.g[center].real), minimum_DOS=float(np.min(final.g.real)),
        maximum_DOS=float(np.max(final.g.real)),
        core_area_mean_DOS=float(np.dot(graph.area_weights[core], final.g[core].real)/normal),
        normalization_residual=final.normalization_residual,
        root_residual=solution.residual, bridge=bridge, continuation=path,
        scope='Causal spectrum of the fixed completed thermal gap; finite positive eta, same prescribed radial exterior; no kinetic populations')
    arrays = dict(coordinates_bar=xy, d=d, a=a, b=b, f=final.f, f_tilde=final.f_tilde,
        g=final.g, DOS=final.g.real, area_weights=graph.area_weights,
        spectral_link_derivative_alpha=final.link_derivative_alpha)
    return metadata, arrays


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, default=DEFAULT_PLAN)
    parser.add_argument('--output-root', type=Path, required=True)
    parser.add_argument('--pilot', action='store_true')
    parser.add_argument('--max-workers', type=int)
    parser.add_argument('--execute', action='store_true')
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text(encoding='utf8'))
    if plan['schema'] != 'pysnspd.stage4.retarded_oracle.v1':
        raise ValueError('Unknown retarded oracle plan')
    input_hashes = {}
    for case in plan['cases']:
        path = ROOT/case['fields_path']
        actual = sha(path)
        if actual != case['fields_sha256']:
            raise ValueError('Completed thermal core fields changed: '+str(path))
        input_hashes[case['fields_path']] = actual
    for name in ('identity.json', 'tasks/R12_n0000.npz'):
        path = RADIAL/name
        input_hashes[path.relative_to(ROOT).as_posix()] = sha(path)
    jobs = []
    for energy in plan['energies_relative']:
        for eta in plan['etas_relative']:
            for case in plan['cases']:
                if args.pilot and (case['id'] not in plan['pilot']['cases'] or energy not in plan['pilot']['energies_relative'] or eta not in plan['pilot']['etas_relative']):
                    continue
                jobs.append(dict(id=job_id(case, eta, energy), plan=plan, case=case, energy_relative=energy, eta_relative=eta))
    if not jobs:
        raise ValueError('No registered jobs selected')
    resources = linux_resources()
    budget = resource_budget(resources, .9, .9, min(len(jobs), args.max_workers or len(jobs)))
    paths = [Path(__file__), ROOT/'pysnspd/experimental/retarded_spatial_usadel.py',
             ROOT/'pysnspd/experimental/thermal_spatial_usadel.py', HERE/'parallel_runtime.py']
    identity = dict(plan_sha256=sha(args.plan), sources={p.relative_to(ROOT).as_posix():sha(p) for p in paths},
        inputs=input_hashes, resources=resources, budget=budget, pilot=args.pilot,
        jobs=[job['id'] for job in jobs], physical_time_steps=0, production_changed=False,
        scope='Retarded spectral oracle only; not a complete nonequilibrium closure or completed stage4')
    if not args.execute:
        print(json.dumps(dict(status='DRY_RUN', **identity), indent=2)); return
    if args.output_root.exists():
        raise FileExistsError('Choose a fresh output directory; no overwrite or implicit resume')
    args.output_root.mkdir(parents=True); (args.output_root/'spectra').mkdir()
    write_json(args.output_root/'identity.json', identity)
    started, records = time.monotonic(), []
    affinity = os.sched_getaffinity(0)
    log = (args.output_root/'progress.jsonl').open('x', encoding='utf8', buffering=1)
    def event(name, **values):
        line = json.dumps(dict(event=name, elapsed_seconds=time.monotonic()-started, **values))
        print(line, flush=True); log.write(line+'\n')
    try:
        event('START', jobs=len(jobs), budget=budget)
        os.sched_setaffinity(0, {budget['coordinator_cpu']})
        context = mp.get_context('spawn'); counter = context.Value('i', 0)
        with ProcessPoolExecutor(max_workers=budget['workers'], mp_context=context,
            initializer=initialize_affinity, initargs=(budget['worker_affinity_cpus'], counter)) as pool:
            futures = {pool.submit(worker, job):job for job in jobs}
            try:
                for future in as_completed(futures):
                    record, arrays = future.result()
                    path = args.output_root/'spectra'/(record['id']+'.npz')
                    np.savez_compressed(path, **arrays)
                    record.update(fields_path=path.relative_to(args.output_root).as_posix(), fields_sha256=sha(path))
                    write_json(path.with_suffix('.json'), record)
                    records.append(record)
                    elapsed = time.monotonic()-started
                    fraction = len(records)/len(jobs)
                    eta = elapsed*(1-fraction)/fraction
                    bar = '#'*int(24*fraction)+'-'*(24-int(24*fraction))
                    print(f'[{bar}] espectros retardados {len(records)}/{len(jobs)} | {elapsed/60:.1f} min | ETA {eta/60:.1f} min', flush=True)
                    event('COMPLETED', id=record['id'], complete=len(records), total=len(jobs), eta_seconds=eta)
            except BaseException:
                for future in futures:
                    future.cancel()
                raise
        summary = dict(status='RETARDED_PILOT_COMPLETE' if args.pilot else 'RETARDED_ORACLE_COMPLETE',
            runtime_seconds=time.monotonic()-started, records=records, stage4_complete=False,
            nonequilibrium_population_closure_admitted=False, production_changed=False, physical_time_steps=0)
        write_json(args.output_root/'summary.json', summary)
        event('COMPLETE', status=summary['status'], seconds=summary['runtime_seconds'])
    except BaseException as exc:
        event('FAILED', error=type(exc).__name__, reason=str(exc))
        write_json(args.output_root/'failure.json', dict(error=type(exc).__name__, reason=str(exc), completed=[x['id'] for x in records]))
        raise
    finally:
        os.sched_setaffinity(0, affinity); log.close()


if __name__ == '__main__':
    main()
