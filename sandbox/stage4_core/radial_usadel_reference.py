"""Prescribed thermal radial vortex: independent spatial Usadel force reference.

No production solver, photon, stationary-gap/barrier search or time trajectory.
Matsubara radial BVPs share one bounded process pool across outer-radius cases.
The raw finite sums and a separately labelled leading omitted-tail estimate are
both reported. No unresolved tail is silently treated as physical precision.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import math
import multiprocessing as mp
import os
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
from parallel_runtime import (THREAD_VARIABLES, initialize_affinity,
                              limit_thread_environment, linux_resources, resource_budget)
limit_thread_environment()
import numpy as np
from scipy.constants import hbar, Boltzmann
from scipy.integrate import solve_bvp, trapezoid
from scipy.optimize import brentq
from scipy.special import digamma, zeta

DEFAULT_PLAN = ROOT/'docs/implementation/stage4/followup_20260923/radial_reference_plan.json'
GAP_ZERO = np.pi*np.exp(-np.euler_gamma)
K0 = np.pi/4


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    path = Path(path)
    with path.open('x', encoding='utf8') as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write('\n')


def thermal_gap(t):
    """BCS equilibrium gap d=|Delta|/(kBTc), with analytic high-frequency tail."""
    if not 0 < t < 1:
        raise ValueError('This reference requires 0<T/Tc<1')
    count = 4096
    eps = 2*np.pi*t*(np.arange(count)+.5)
    def residual(d):
        en = np.hypot(eps, d)
        total = np.sum(d*d/(eps*en*(eps+en)))
        total += sum(c*zeta(power, count+.5)/(2*np.pi*t)**power for power, c in
                     ((3, d*d/2), (5, -3*d**4/8), (7, 5*d**6/16)))
        return float(np.log(t)+2*np.pi*t*total)
    gap = float(brentq(residual, .001, 4., xtol=2e-13))
    return gap, residual(gap)


def profile(radius, amplitude, kind='tanh'):
    """Amplitude, first/second derivatives and L1 d; r is in ell0 units."""
    r = np.asarray(radius, float)
    if np.any(r <= 0):
        raise ValueError('Use a positive inner radius; the regular limit is supplied analytically')
    if kind == 'uniform':
        return (np.full_like(r, amplitude), np.zeros_like(r), np.zeros_like(r),
                -amplitude/(r*r))
    if kind == 'linear':
        return amplitude*r, np.full_like(r, amplitude), np.zeros_like(r), np.zeros_like(r)
    if kind != 'tanh':
        raise ValueError('Unknown prescribed profile')
    th = np.tanh(r)
    sech2 = 1-th*th
    d, dp, dpp = amplitude*th, amplitude*sech2, -2*amplitude*th*sech2
    lap = dpp+dp/r-d/(r*r)
    small = r < .01
    rr = r[small]
    lap[small] = amplitude*(-8*rr/3+16*rr**3/5-272*rr**5/105+992*rr**7/567)
    return d, dp, dpp, lap


def uniform_angle(amplitude, eps, gamma):
    """Unique positive thermal root, bracketed in tan(theta), with no clipping.

    d=eps*tan(theta)+gamma*tan(theta)/sqrt(1+tan(theta)^2) is strictly
    increasing for eps>0. Thus the physical root lies in [0,d/eps].
    """
    if amplitude < 0 or eps <= 0 or gamma < 0:
        raise ValueError('Nonnegative amplitude/depairing and positive Matsubara frequency required')
    if amplitude == 0:
        return 0.
    if gamma == 0:
        return float(np.arctan(amplitude/eps))
    root = brentq(lambda value: eps*value+gamma*value/np.hypot(1., value)-amplitude,
                  0., amplitude/eps, xtol=1e-14, rtol=1e-13)
    return float(np.arctan(root))


def output_grid(radius, rmin=.001):
    if radius < 4 or rmin <= 0 or rmin >= .25:
        raise ValueError('Require R>=4 and 0<rmin<.25 for the common inner-core grid')
    return np.r_[np.geomspace(rmin, .25, 180, endpoint=False),
                 np.linspace(.25, 4., 301, endpoint=False),
                 (np.array([4.]) if radius == 4 else np.linspace(4., radius, 181))]


def solve_mode(eps, amplitude, radius, *, rmin=.001, winding=1, kind='tanh',
               tolerance=1e-5, max_nodes=12000, grid=None):
    """One spectral boundary-value problem on the regular Matsubara branch.

    theta''+theta'/r-m^2 sin(theta)cos(theta)/r^2-eps sin(theta)+d cos(theta)=0.
    For m=1, theta'(rmin)=theta(rmin)/rmin is the leading regular Robin limit.
    m=0, uniform d is exposed only for an independent exact-solution test.
    """
    if eps <= 0 or amplitude < 0 or radius <= rmin or winding not in (0, 1):
        raise ValueError('Invalid radial reference parameters')
    r = output_grid(radius, rmin) if grid is None else np.asarray(grid, float)
    if np.any(np.diff(r) <= 0) or r[0] != rmin or r[-1] != radius:
        raise ValueError('Output grid must increase exactly from rmin to R')
    d_outer = float(profile(np.array([radius]), amplitude, kind)[0][0])
    outer = uniform_angle(d_outer, eps, winding*winding/(radius*radius))
    initial = np.r_[np.geomspace(rmin, .2, 80, endpoint=False), np.linspace(.2, radius, 181)]
    if winding == 0 and kind == 'uniform':
        guess = np.vstack((np.full_like(initial, outer), np.zeros_like(initial)))
    else:
        guess = np.vstack((outer*np.tanh(initial)/np.tanh(radius),
                           outer*(1-np.tanh(initial)**2)/np.tanh(radius)))
    def ode(x, y):
        d = profile(x, amplitude, kind)[0]
        sine, cosine = np.sin(y[0]), np.cos(y[0])
        return np.vstack((y[1], -y[1]/x+winding*winding*sine*cosine/(x*x)+eps*sine-d*cosine))
    def jac(x, y):
        d = profile(x, amplitude, kind)[0]
        value = np.zeros((2, 2, len(x)))
        value[0, 1] = 1.
        value[1, 0] = winding*winding*np.cos(2*y[0])/(x*x)+eps*np.cos(y[0])+d*np.sin(y[0])
        value[1, 1] = -1/x
        return value
    def boundary(left, right):
        return np.array([left[1]-winding*left[0]/rmin, right[0]-outer])
    def boundary_jac(left, right):
        return np.array([[-winding/rmin, 1.], [0., 0.]]), np.array([[0., 0.], [1., 0.]])
    start = time.monotonic()
    solution = solve_bvp(ode, boundary, initial, guess, fun_jac=jac,
                         bc_jac=boundary_jac, tol=tolerance, max_nodes=max_nodes)
    if not solution.success:
        raise RuntimeError(f'Radial BVP failed: eps={eps:g}, R={radius:g}: {solution.message}')
    theta, derivative = solution.sol(r)
    if np.any(~np.isfinite(theta)) or np.min(theta) < 0 or np.max(theta) >= np.pi/2:
        raise RuntimeError('Spectral angle left the positive physical Matsubara branch; no clipping')
    return dict(radius=r, theta=theta, derivative=derivative,
        eps=float(eps), outer_angle=outer, mesh_nodes=int(len(solution.x)),
        maximum_rms_residual=float(max(solution.rms_residuals)),
        boundary_residual=boundary(solution.y[:, 0], solution.y[:, -1]).tolist(),
        runtime_seconds=time.monotonic()-start)


def uniform_moments(amplitude, q, count, t):
    """Thermal amplitude and q derivatives, historical C.7 digamma subtraction.

    The physical roots use a monotone vector bisection (tan theta), not clipped
    Newton iterations. The linear-in-d large-Gamma tail is summed exactly;
    the remaining leading nonlinear tail has eps+Gamma in its denominator.
    """
    d, qq = np.broadcast_arrays(np.atleast_1d(amplitude).astype(float), np.atleast_1d(q).astype(float))
    if np.any(d < 0) or np.any(qq < 0) or count < 1 or t <= 0:
        raise ValueError('Invalid uniform thermal moments')
    gamma = qq*qq
    eps = 2*np.pi*t*(np.arange(count)[:, None]+.5)
    lower = np.zeros((count, len(d)))
    upper = d[None, :]/eps
    # Exactly bracketed monotone roots. 48 halvings give a narrow relative
    # bracket even when Gamma is large; no iterate is projected onto a branch.
    for _ in range(48):
        middle = (lower+upper)/2
        positive = eps*middle+gamma[None, :]*middle/np.sqrt(1+middle*middle) > d[None, :]
        upper = np.where(positive, middle, upper)
        lower = np.where(positive, lower, middle)
    tangent = (lower+upper)/2
    sine = tangent/np.sqrt(1+tangent*tangent)
    residual = eps*tangent+gamma[None, :]*sine-d[None, :]
    if np.max(abs(residual)) > 1e-9*max(1., float(np.max(d))):
        raise RuntimeError('Uniform bracket residual unresolved; no fallback')
    den = eps+gamma[None, :]
    lin = np.log(t)+digamma(.5+gamma/(2*np.pi*t))-digamma(.5)
    tail = lambda power: zeta(power, count+.5+gamma/(2*np.pi*t))/(2*np.pi*t)**power
    force = 2*d*lin+4*np.pi*t*(np.sum(d[None, :]/den-sine, axis=0)
                            +d**3/2*(tail(3)-gamma*tail(4)))
    moment = 4*np.pi*t*qq*(np.sum(sine*sine, axis=0)+d*d*tail(2)
                           -d**4*(tail(4)-gamma*tail(5)))
    return force, moment, float(np.max(abs(residual)))


def candidate_force(radius, amplitude, fraction, t, count=600):
    """Thermal C.9/C.7.2 candidate on the same prescribed radial vortex."""
    r = np.asarray(radius, float)
    d, dp, dpp, lap = profile(r, amplitude)
    delta = fraction*GAP_ZERO
    rho = d*d
    denominator = rho+delta*delta
    q = rho/(denominator*r)
    x, moment, residual = uniform_moments(d, q, count, t)
    dq_da = 2*d*delta*delta/(denominator*denominator*r)
    force = x-2*K0*lap-2*K0*d*q*q+(moment-2*K0*rho*q)*dq_da
    return dict(radius=r, amplitude=d, force=force, q_delta=q,
                uniform_force=x, uniform_q_moment=moment, root_residual_max=residual)


def _worker(job):
    started = time.monotonic()
    if job['kind'] == 'mode':
        result = solve_mode(job['eps'], job['gap'], job['R'], rmin=job['rmin'],
                            tolerance=job['tolerance'], max_nodes=job['max_nodes'])
    else:
        result = candidate_force(output_grid(job['R'], job['rmin']), job['gap'],
                                 job['fraction'], job['t'], job['count'])
    return dict(job=job, values=result, worker_seconds=time.monotonic()-started, pid=os.getpid())


def force_sum(radius, amplitude, theta, t):
    eps = 2*np.pi*t*(np.arange(len(theta))+.5)
    d, _, _, lap = profile(radius, amplitude)
    raw = 2*d*np.log(t)+4*np.pi*t*np.sum(d[None, :]/eps[:, None]-np.sin(theta), axis=0)
    leading_tail = -4*np.pi*t*lap*zeta(2, len(theta)+.5)/(2*np.pi*t)**2
    return raw, leading_tail


def norm(radius, values):
    return float(np.sqrt(trapezoid(2*np.pi*radius*values*values, radius)))


def compare_force(radius, reference, candidate):
    ref_norm = norm(radius, reference)
    difference = candidate-reference
    return dict(reference_L2=ref_norm, candidate_L2=norm(radius, candidate),
        difference_L2=norm(radius, difference),
        relative_L2_difference=norm(radius, difference)/ref_norm if ref_norm else None,
        reference_peak_absolute=float(max(abs(reference))), candidate_peak_absolute=float(max(abs(candidate))),
        maximum_absolute_difference=float(max(abs(difference))),
        reference_peak_radius=float(radius[np.argmax(abs(reference))]),
        candidate_peak_radius=float(radius[np.argmax(abs(candidate))]))


def tasks(plan, gap, *, pilot=False):
    common = dict(gap=gap, t=plan['T_K']/plan['Tc_K'], rmin=plan['rmin_ell0'],
                  tolerance=plan['bvp_relative_tolerance'], max_nodes=plan['bvp_max_nodes'])
    jobs = []
    if pilot:
        for n in plan['pilot_modes']:
            jobs.append(dict(**common, kind='mode', id=f'R8_n{n:04d}', R=8., n=n,
                             eps=2*np.pi*common['t']*(n+.5)))
        return jobs
    # One fair queue; the R=8 first128 modes are reused by its N=128 sum.
    radii = list(dict.fromkeys(case['R_ell0'] for case in plan['cases']))
    for radius in radii:
        for fraction in plan['candidate_delta_fractions']:
            jobs.append(dict(**common, kind='candidate', id=f'R{radius:g}_delta{fraction:g}',
                R=radius, fraction=fraction, count=plan['candidate_matsubara_count']))
    maximum = {radius:max(case['matsubara_count'] for case in plan['cases'] if case['R_ell0']==radius) for radius in radii}
    for n in range(max(maximum.values())):
        for radius in radii:
            if n < maximum[radius]:
                jobs.append(dict(**common, kind='mode', id=f'R{radius:g}_n{n:04d}', R=radius,
                    n=n, eps=2*np.pi*common['t']*(n+.5)))
    return jobs


def summarize(plan, gap, stored, output):
    t = plan['T_K']/plan['Tc_K']
    reports = {}
    curves = {}
    for case in plan['cases']:
        radius, count = case['R_ell0'], case['matsubara_count']
        r = output_grid(radius, plan['rmin_ell0'])
        theta = np.array([stored[f'R{radius:g}_n{n:04d}']['theta'] for n in range(count)])
        raw, tail = force_sum(r, gap, theta, t)
        core = r <= plan['comparison_window_ell0']
        rc = r[core]
        entry = dict(R_ell0=radius, matsubara_count=count, core_force_L2=norm(rc, raw[core]),
            leading_tail_estimate_L2=norm(rc, tail[core]), candidate_comparisons={})
        arrays = dict(radius=r, theta=theta, force_raw=raw, leading_tail_estimate=tail,
                      force_with_leading_tail_estimate=raw+tail, core_mask=core)
        for fraction in plan['candidate_delta_fractions']:
            candidate = stored[f'R{radius:g}_delta{fraction:g}']['force']
            entry['candidate_comparisons'][str(fraction)] = dict(
                versus_raw_finite_sum=compare_force(rc, raw[core], candidate[core]),
                versus_leading_tail_estimate=compare_force(rc, (raw+tail)[core], candidate[core]))
            arrays[f'candidate_force_delta_{fraction:g}'] = candidate
        case_id = case['id']
        np.savez_compressed(output/(case_id+'.npz'), **arrays)
        entry['fields_file'] = case_id+'.npz'
        entry['fields_sha256'] = sha(output/entry['fields_file'])
        reports[case_id] = entry
        curves[case_id] = (rc, raw[core], tail[core])
    comparisons = {}
    for name, first, second in [('cutoff_R8_N128_to_N256','R8_N128','R8_N256'),
                               ('outer_R8_to_R12_N256','R8_N256','R12_N256')]:
        ra, a, ta = curves[first]
        rb, b, tb = curves[second]
        if not np.array_equal(ra, rb):
            raise ValueError('Core comparison coordinates must be identical')
        comparisons[name] = dict(raw=compare_force(rb, b, a),
            with_leading_tail_estimate=compare_force(rb, b+tb, a+ta),
            interpretation='Observed change, not a rigorous tail bound or physical confidence interval')
    return dict(cases=reports, comparisons=comparisons,
        status='THERMAL_FORCE_REFERENCE_COMPUTED_NOT_PHYSICAL_CORE_ADMISSION',
        interpretation=plan['interpretation'], stage4_complete=False, physical_core_admitted=False,
        physical_time_steps=0, photon=False, energy_or_barrier_validation=False,
        primary_quantity='Raw finite Matsubara force; leading omitted-tail estimate reported separately',
        tail_formula='DeltaX_leading=-4*pi*(T/Tc)*L1[d]*sum_{n>=N} eps_n^-2; L1=d2/dr2+(1/r)d/dr-1/r2',
        source_equations='C.7.2 candidate; spatial Usadel specialization of A and D.4.4 physical reference',
        core_window_ell0=plan['comparison_window_ell0'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, default=DEFAULT_PLAN)
    parser.add_argument('--output-root', type=Path, required=True)
    parser.add_argument('--max-workers', type=int)
    parser.add_argument('--cpu-fraction', type=float, default=.9)
    parser.add_argument('--memory-fraction', type=float, default=.9)
    parser.add_argument('--pilot', action='store_true')
    parser.add_argument('--execute', action='store_true')
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text(encoding='utf8'))
    if plan['schema'] != 'pysnspd.stage4.radial_thermal_reference.v1':
        raise ValueError('Unknown radial reference plan')
    gap, gap_residual = thermal_gap(plan['T_K']/plan['Tc_K'])
    jobs = tasks(plan, gap, pilot=args.pilot)
    linux = hasattr(os, 'sched_getaffinity')
    if not linux and not args.pilot:
        raise RuntimeError('Full campaign requires Linux affinity and resource budget')
    budget = (resource_budget(linux_resources(), args.cpu_fraction, args.memory_fraction, args.max_workers)
              if linux else dict(workers=1, total_processes=1, scope='Serial bounded local pilot only'))
    source_paths = [Path(__file__), HERE/'parallel_runtime.py']
    identity = dict(plan_sha256=sha(args.plan), sources={p.relative_to(ROOT).as_posix():sha(p) for p in source_paths},
        budget=budget, jobs=len(jobs), pilot=args.pilot, gap_kBTc=gap, gap_residual=gap_residual,
        ell0_nm=float(np.sqrt(hbar*plan['D_m2_s']/(2*Boltzmann*plan['Tc_K']))*1e9),
        restrictions=plan['interpretation'])
    if not args.execute:
        print(json.dumps(dict(status='DRY_RUN_NO_RADIAL_SOLVES', **identity), indent=2))
        return
    if args.output_root.exists():
        raise FileExistsError('Choose a new output root; previous runs are preserved')
    args.output_root.mkdir(parents=True)
    (args.output_root/'tasks').mkdir()
    write_json(args.output_root/'identity.json', identity)
    started = time.monotonic()
    log = (args.output_root/'progress.jsonl').open('x', encoding='utf8', buffering=1)
    stored, records = {}, []
    last_print = 0.
    affinity = os.sched_getaffinity(0) if linux else None
    def event(name, **fields):
        record = dict(event=name, elapsed_seconds=time.monotonic()-started, **fields)
        line = json.dumps(record, allow_nan=False)
        log.write(line+'\n')
        print(line, flush=True)
    def accept(result):
        nonlocal last_print
        job, values = result['job'], result['values']
        stored[job['id']] = values
        arrays = {key:value for key,value in values.items() if isinstance(value, np.ndarray)}
        meta = {key:value for key,value in values.items() if not isinstance(value, np.ndarray)}
        file = args.output_root/'tasks'/(job['id']+'.npz')
        np.savez_compressed(file, **arrays)
        record = dict(id=job['id'], kind=job['kind'], job=job, diagnostics=meta,
            worker_seconds=result['worker_seconds'], worker_pid=result['pid'],
            arrays=file.relative_to(args.output_root).as_posix(), sha256=sha(file))
        records.append(record)
        now = time.monotonic()
        if len(records)==len(jobs) or now-last_print > 3:
            last_print = now
            fraction = len(records)/len(jobs)
            remaining = (now-started)*(len(jobs)-len(records))/len(records)
            bar = '#'*int(24*fraction)+'-'*(24-int(24*fraction))
            print(f'[{bar}] referencia {len(records)}/{len(jobs)} | transcurrido {(now-started)/60:.1f} min | ETA {remaining/60:.1f} min', flush=True)
            event('PROGRESS', complete=len(records), total=len(jobs), eta_seconds=remaining)
    try:
        event('START', pilot=args.pilot, jobs=len(jobs), budget=budget)
        if linux:
            os.sched_setaffinity(0, {budget['coordinator_cpu']})
            context = mp.get_context('spawn')
            counter = context.Value('i', 0)
            with ProcessPoolExecutor(max_workers=budget['workers'], mp_context=context,
                initializer=initialize_affinity, initargs=(budget['worker_affinity_cpus'], counter)) as pool:
                futures = {pool.submit(_worker, job):job for job in jobs}
                try:
                    for future in as_completed(futures):
                        accept(future.result())
                except BaseException:
                    for future in futures:
                        future.cancel()
                    raise
        else:
            for job in jobs:
                accept(_worker(job))
        if args.pilot:
            summary = dict(status='BOUNDED_PILOT_ONLY_NO_FORCE_COMPARISON',
                mode_count=len(records), stage4_complete=False, physical_core_admitted=False,
                summed_worker_seconds=sum(row['worker_seconds'] for row in records),
                maximum_mode_seconds=max(row['worker_seconds'] for row in records),
                note='Pilot samples selected modes; it is not a truncated physical-force result')
        else:
            summary = summarize(plan, gap, stored, args.output_root)
        summary.update(runtime_seconds=time.monotonic()-started, task_records=records)
        write_json(args.output_root/'summary.json', summary)
        event('COMPLETE', output=str(args.output_root/'summary.json'))
    except BaseException as exc:
        write_json(args.output_root/'failure.json', dict(status='STOPPED', reason=str(exc),
            exception=type(exc).__name__, completed=records,
            policy='No retry, clipping, automatic fallback, overwritten output or continuation'))
        event('STOPPED', reason=str(exc))
        raise
    finally:
        log.close()
        if affinity is not None:
            os.sched_setaffinity(0, affinity)


if __name__ == '__main__':
    main()
