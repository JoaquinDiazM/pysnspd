"""Finite-frequency admission of the missing dynamic coupling, without a photon.

Uniform equilibrium permits an exact decomposition in the eigenmodes of the
actual dual-mesh Laplacian. Energy, charge and anomalous moments are evaluated
from the SAME Keldysh response. This is a weak-response experiment, not a new
nonlinear production solver or a fitted quasiparticle lifetime.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor, FIRST_COMPLETED, wait
import hashlib
import json
import math
import multiprocessing as mp
import os
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from sandbox.stage4_core.parallel_runtime import (limit_thread_environment,
    linux_resources, resource_budget, initialize_affinity, THREAD_VARIABLES)
from sandbox.stage4_core.coupled_campaign import atomic_json, sha
limit_thread_environment()
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.constants import hbar, k as KB, e as EC
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import eigsh, spsolve
from sandbox.stage4_core.dual_kwt_time import load_graph, uniform_reference

I = np.eye(2, dtype=complex)
Z = np.diag([1., -1.]).astype(complex)
DX = np.array([[0, 1], [1, 0]], complex)
DY = np.array([[0, 1j], [-1j, 0]], complex)
FORCING = np.stack((DX, DY, 1j*I))


def encode(value):
    value = np.asarray(value)
    return {'real': value.real.tolist(), 'imag': value.imag.tolist()}


def energy_rule(d, omega, eta, order, cutoff, tail_order=0):
    """Resolve shifted edges, optionally integrating both infinite tails.

    The reciprocal coordinate x=cutoff/|E| maps each tail to (0,1), with
    d|E|=cutoff dx/x**2. Every tail sample evaluates the SAME physical
    integrand: there is no fitted counterterm, artificial relaxation or DOS
    renormalization. tail_order=0 preserves the archived finite-cutoff rule.
    """
    points = [-cutoff, 0., cutoff]
    for center in (-d-omega/2, -d+omega/2, d-omega/2, d+omega/2, -d, d):
        points.extend(center+eta*np.array([-64, -16, -4, -1, 0, 1, 4, 16, 64]))
    # These fixed thermal panels resolve tanh independently of the gap edges.
    points.extend(np.linspace(-2., 2., 17))
    if cutoff > 2:
        outer = np.geomspace(2., cutoff, 12)
        points.extend(outer); points.extend(-outer)
    points = np.unique(np.clip(points, -cutoff, cutoff))
    x, w = leggauss(order)
    widths = np.diff(points)/2
    energy = ((points[:-1, None]+points[1:, None])/2+widths[:, None]*x).ravel()
    weights = (widths[:, None]*w).ravel()
    if type(tail_order) is not int or tail_order < 0:
        raise ValueError('tail_order must be a nonnegative integer')
    if tail_order:
        if not np.isfinite(cutoff) or cutoff <= 0:
            raise ValueError('Positive finite transition energy required for reciprocal tails')
        nodes, quadrature = leggauss(tail_order)
        coordinate, measure = (nodes+1)/2, quadrature/2
        tail_energy = cutoff/coordinate
        tail_weights = cutoff*measure/coordinate**2
        energy = np.r_[-tail_energy, energy, tail_energy[::-1]]
        weights = np.r_[tail_weights, weights, tail_weights[::-1]]
    return energy, weights


def backgrounds(d, energy, eta):
    z = eta-1j*np.asarray(energy)
    root = np.sqrt(z*z+d*d)
    root = np.where(root.real < 0, -root, root)
    B = d*DX+z[:, None, None]*Z
    R = B/root[:, None, None]
    A = -Z@R.conj().transpose(0, 2, 1)@Z
    return R, A, root, d*DX-1j*np.asarray(energy)[:, None, None]*Z


def spectral_response(Rp, Rm, Ap, Am, bp, bm, lam, forcing=FORCING):
    dR = (forcing-Rp[:, None]@forcing@Rm[:, None])/(bp+bm+2*lam)[:, None, None, None]
    dA = (forcing-Ap[:, None]@forcing@Am[:, None])/(bp.conj()+bm.conj()+2*lam)[:, None, None, None]
    return dR, dA


def projection(matrix):
    return np.stack(((matrix[..., 0, 0]+matrix[..., 1, 1])/4,
                     (matrix[..., 0, 0]-matrix[..., 1, 1])/4), axis=-1)


def gap_moments(matrix):
    return np.stack((.5j*(matrix[..., 0, 1]+matrix[..., 1, 0]),
                     .5*(matrix[..., 0, 1]-matrix[..., 1, 0])), axis=-1)


def kinetic_modal(Rp, Rm, Ap, Am, Bp, Bm, dR, dA, hp, hm, lam):
    """Three independent gap/potential forcings; no eta collision is added."""
    Kp, Km = hp[:, None, None]*(Rp-Ap), hm[:, None, None]*(Rm-Am)
    basis = np.stack((I, Z))
    Kbasis = Rp[:, None]@basis-basis@Am[:, None]
    source = dR*hm[:, None, None, None]-hp[:, None, None, None]*dA
    def local(K):
        return -.5*(Bp[:, None]@K-K@Bm[:, None])
    def divergence(K, dr, da):
        return -lam/2*(Rp[:, None]@K-dr@Km[:, None]+Kp[:, None]@da-K@Am[:, None])
    coefficient = projection(local(Kbasis)+divergence(Kbasis, 0.*I, 0.*I)).transpose(0, 2, 1)
    rhs_source = projection(local(source)+divergence(source, dR, dA)
        -.5*(FORCING@Km[:, None]-Kp[:, None]@FORCING)).transpose(0, 2, 1)
    distribution = np.linalg.solve(coefficient, -rhs_source)
    K = source+np.einsum('ecf,ecij->efij', distribution, Kbasis)
    residual = coefficient@distribution+rhs_source
    return K, distribution, float(np.max(abs(residual))), coefficient


def thermal_stiffness(d, t, epsilon, lam):
    g = epsilon/np.hypot(epsilon, d)
    radial = 2*np.log(t)+4*np.pi*t*np.sum(1/epsilon-g**3/(epsilon+g*lam))
    phase = 4*np.pi*t*np.sum(g*g*lam/(epsilon*(epsilon+g*lam)))
    return float(radial), float(phase)


def kernel(task):
    started = time.monotonic()
    d, t, omega, eta, lam = [task[key] for key in ('d', 't', 'omega', 'eta', 'lambda')]
    energy, weights = energy_rule(d, omega, eta, task['order'], task['cutoff'], task.get('tail_order', 0))
    Rp, Ap, bp, Bp = backgrounds(d, energy+omega/2, eta)
    Rm, Am, bm, Bm = backgrounds(d, energy-omega/2, eta)
    hp, hm = np.tanh((energy+omega/2)/(2*t)), np.tanh((energy-omega/2)/(2*t))
    dR, dA = spectral_response(Rp, Rm, Ap, Am, bp, bm, lam)
    K, distribution, kinetic_error, coefficient = kinetic_modal(Rp, Rm, Ap, Am, Bp, Bm, dR, dA, hp, hm, lam)
    R0, A0, b0, _ = backgrounds(d, energy, eta)
    sR, sA = spectral_response(R0, R0, A0, A0, b0, b0, lam)
    static = np.tanh(energy/(2*t))[:, None, None, None]*(sR-sA)
    # A static scalar potential is not part of the equilibrium gap Hessian.
    static[:, 2] = 0.
    correction = np.einsum('e,efc->cf', weights/2, gap_moments(K-static))
    charge = np.einsum('e,ef->f', weights/2, projection(K)[..., 0])
    epsilon = 2*np.pi*t*(np.arange(task['matsubara_count'])+.5)
    hx, hy = thermal_stiffness(d, t, epsilon, lam)
    matrix = np.vstack((correction, charge))
    matrix += np.diag([hx, hy, 1.])
    # These are measured kernels; no frequency-by-frequency mobility subtraction
    # is silently installed as a causal time-domain material law.
    reduced = matrix[0, 0]-matrix[0, 1:]@np.linalg.solve(matrix[1:, 1:], matrix[1:, 0])
    gamma_resolved = -2*reduced.imag/omega
    tau = 1/(t/task['tau_ee_Tc_ps']+t**3/task['tau_ep_Tc_ps'])
    kappa = (tau/task['tD_ps'])**2
    gamma_kwt = 2*np.sqrt((1+t)/2)*(np.pi/4)*np.sqrt(1+kappa*d*d)
    hybrid = matrix.copy()
    hybrid[0, 0] -= .5j*omega*gamma_kwt
    forcing = np.diag([hx, max(hy, 1e-12), 1.])
    response = np.linalg.solve(hybrid, forcing)
    # Charge and anomalous populations follow this *same* coupled solution.
    h_response = distribution@response
    return dict(kind='mode', mode=task['mode'], lambda_value=lam, omega=omega,
        eta=eta, quadrature_order=task['order'], energy_cutoff=task['cutoff'], energy_nodes=len(energy),
        tail_order=task.get('tail_order', 0),
        kernel=encode(matrix), candidate_response=encode(response),
        gamma_resolved=float(gamma_resolved), gamma_kwt=float(gamma_kwt),
        resolved_over_kwt=float(gamma_resolved/gamma_kwt),
        kinetic_residual=kinetic_error,
        charge_solve_residual=float(np.max(abs((hybrid@response-forcing)[2]))),
        maximum_distribution_response=float(np.max(abs(h_response))),
        maximum_kinetic_condition=float(np.max(np.linalg.cond(coefficient))),
        status='MEASURED_NOT_A_NONLINEAR_CLOSURE', elapsed_seconds=time.monotonic()-started)


def film_port_coefficients(task):
    """Intrinsic film coefficients, before assigning exterior circuit elements.

    Keeping this separate permits the normal-state conductivity check without
    pretending that a normal metal has a finite superconducting inductance.
    """
    d, t, omega, eta = [task[key] for key in ('d', 't', 'omega', 'eta')]
    energy, weights = energy_rule(d, omega, eta, task['order'], task['cutoff'], task.get('tail_order', 0))
    Rp, Ap, bp, Bp = backgrounds(d, energy+omega/2, eta)
    Rm, Am, bm, Bm = backgrounds(d, energy-omega/2, eta)
    hp, hm = np.tanh((energy+omega/2)/(2*t)), np.tanh((energy-omega/2)/(2*t))
    Kp, Km = hp[:, None, None]*(Rp-Ap), hm[:, None, None]*(Rm-Am)
    dr, da = (Z@Rm-Rp@Z)/omega, (Z@Am-Ap@Z)/omega
    dk = (Z@Km-Kp@Z)/omega
    flux = .5*(Rp@dk-dr@Km+Kp@da-dk@Am)
    current = np.dot(weights, projection(flux)[:, 1])
    R0, A0, b0, _ = backgrounds(d, energy, eta)
    h0 = np.tanh(energy/(2*t))
    K0 = h0[:, None, None]*(R0-A0)
    dB0 = (-2j*d/omega)*DY
    r0, a0 = spectral_response(R0, R0, A0, A0, b0, b0, 0., dB0[None])
    r0, a0 = r0[:, 0], a0[:, 0]
    k0 = h0[:, None, None]*(r0-a0)
    flux0 = .5*(R0@k0-r0@K0+K0@a0-k0@A0)
    current0 = np.dot(weights, projection(flux0)[:, 1])
    ep = 2*np.pi*t*(np.arange(task['matsubara_count'])+.5)
    stiffness = 4*np.pi*t*np.sum(d*d/(ep*ep+d*d))
    G = task['graph_conductance']
    ybar = 2j*stiffness*G/omega-G*(current-current0)
    admittance = ybar/(2*task['sheet_resistance_ohm'])
    neutral = 1+np.dot(weights/2, projection(dk)[:, 0])
    return dict(admittance=admittance, neutral=neutral, current=current,
        current0=current0, stiffness=stiffness, energy_nodes=len(energy))


def port(task):
    """Voltage-driven film admittance; superconducting contacts obey Josephson."""
    started = time.monotonic()
    omega, eta = task['omega'], task['eta']
    film = film_port_coefficients(task)
    admittance, neutral = film['admittance'], film['neutral']
    current0, stiffness = film['current0'], film['stiffness']
    G = task['graph_conductance']
    from pysnspd.experimental.electrical_ports import ThesisCircuitParameters
    from pysnspd.experimental.thesis_circuit_harmonic import solve_harmonic
    E0 = KB*task['Tc_K']
    omega_si = omega*E0/hbar
    resolved_L = hbar*task['sheet_resistance_ohm']/(E0*stiffness*G)
    external_L = task['total_inductance_H']-resolved_L
    if external_L <= 0:
        raise ValueError('Resolved inductance exceeds the thesis total; no negative exterior allowed')
    parameters = ThesisCircuitParameters(Lk_ext_H=external_L)
    # Circuit helper originally uses exp(+i omega t): an explicit conjugation
    # preserves its convention and converts back on output.
    circuit = solve_harmonic(parameters, omega_si, np.conj(1/admittance), bias_voltage_peak_V=1e-6)
    return dict(kind='port', omega=omega, eta=eta, quadrature_order=task['order'],
        energy_cutoff=task['cutoff'], energy_nodes=film['energy_nodes'], tail_order=task.get('tail_order', 0),
        admittance_S=encode(admittance), impedance_ohm=encode(1/admittance),
        neutral_gauge_residual=encode(neutral),
        static_retarded_current_coefficient=encode(current0),
        thermal_stiffness=float(stiffness), resolved_inductance_H=float(resolved_L),
        external_inductance_H=float(external_L),
        bias_voltage_peak_V=1e-6, time_convention='exp(-iwt)',
        circuit_state_peak=encode(circuit.state_peak.conj()),
        readout_peak_V=encode(np.conj(circuit.readout_voltage_peak_V)),
        circuit_power_residual_W=float(circuit.average_power_residual_W),
        film_average_power_per_voltage_squared_W=float(admittance.real/2),
        status='WEAK_PORT_RESPONSE', elapsed_seconds=time.monotonic()-started)


def mesh_modes(mesh, count):
    graph, profile = load_graph(mesh)
    tail, head = graph.edges.T
    c = graph.conductance
    lap = coo_matrix((np.r_[c, c, -c, -c], (np.r_[tail, head, tail, head], np.r_[tail, head, head, tail])),
        shape=(graph.n_nodes, graph.n_nodes)).tocsc()
    free = np.ones(graph.n_nodes, bool); free[graph.boundary_nodes] = False
    values, vectors = eigsh(lap[free][:, free], k=count, M=diags(graph.area_weights[free]), sigma=0., which='LM')
    indices = np.argsort(values); values, vectors = values[indices], vectors[:, indices]
    full = np.zeros((graph.n_nodes, count)); full[free] = vectors
    coords = graph.coordinates_bar
    xmin, xmax = coords[:, 0].min(), coords[:, 0].max()
    lift = np.zeros(graph.n_nodes)
    lift[graph.boundary_nodes] = np.where(coords[graph.boundary_nodes, 0] < (xmin+xmax)/2, .5, -.5)
    lift[free] = spsolve(lap[free][:, free], -(lap@lift)[free])
    conductance = float(lift@(lap@lift))
    error = lap@full-graph.area_weights[:, None]*full*values
    return graph, values, full, lift, conductance, float(np.max(abs(error[free])))


def execute_task(task):
    value = port(task) if task['kind'] == 'port' else kernel(task)
    value['worker'] = runtime_telemetry()
    return value


def runtime_telemetry():
    import resource
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return dict(pid=os.getpid(), affinity_cpus=sorted(os.sched_getaffinity(0)),
        peak_rss_bytes=int(usage.ru_maxrss)*1024,
        user_cpu_seconds=usage.ru_utime, system_cpu_seconds=usage.ru_stime)


def admit_runtime(workers, resources, inherited=None):
    """Use one shared 90% admission; never reserve another 10% inside a case."""
    if type(workers) is not int or workers < 1:
        raise ValueError('workers must be a positive integer')
    if inherited is None:
        budget = resource_budget(resources, max_workers=workers)
        if budget['workers'] < workers:
            raise ValueError('Requested workers exceed the fresh standalone CPU/RAM budget')
        return dict(source='standalone_fresh_admission', **budget)
    cpus = inherited['cpus']
    shared = inherited['shared_budget']
    if (inherited.get('parent_pid') != os.getppid() or len(set(cpus)) != len(cpus)
            or set(cpus) != set(resources['logical_cpus'])
            or not set(cpus) <= set(shared['selected_affinity_cpus'])
            or inherited['workers'] != len(cpus)-1 or workers > inherited['workers']):
        raise ValueError('Inherited campaign allocation does not match the actual child affinity')
    ceiling = inherited['total_process_cpu_ceiling']
    if ceiling > shared['logical_budget'] or ceiling > len(shared['selected_affinity_cpus']):
        raise ValueError('Inherited campaign exceeds its shared CPU budget')
    quota = resources.get('cpu_quota_units')
    if quota is not None and ceiling > math.floor(.9*quota):
        raise ValueError('The current cgroup CPU quota is tighter than campaign admission')
    reserved = shared['coordinator_reserve_bytes']+workers*shared['worker_reserve_bytes']
    if reserved > math.floor(.9*resources['available_memory_bytes']):
        raise ValueError('Available RAM no longer admits this case coordinator and workers')
    return dict(source='inherited_shared_campaign_admission', parent_allocation=inherited,
        workers=workers, coordinator_cpu=cpus[-1], worker_affinity_cpus=cpus[:workers],
        selected_affinity_cpus=cpus, threads_per_process=1, total_processes=workers+1,
        estimated_memory_reservation_bytes=reserved)


def source_hashes():
    names = ['sandbox/stage4_core/coupled_response.py', 'sandbox/stage4_core/coupled_campaign.py',
        'sandbox/stage4_core/parallel_runtime.py', 'sandbox/stage4_core/dual_kwt_time.py',
        'pysnspd/experimental/thesis_circuit_harmonic.py', 'pysnspd/experimental/electrical_ports.py',
        'pysnspd/experimental/thermal_spatial_usadel.py']
    return {name: sha(ROOT/name) for name in names if (ROOT/name).is_file()}


def validate_plan(plan):
    if plan.get('schema') != 'pysnspd.stage4.coupled_response.v1':
        raise ValueError('Unsupported plan')
    for name in ('omegas', 'eta_over_gap', 'orders', 'cutoffs'):
        values = plan.get(name)
        if not isinstance(values, list) or not values or any(
                isinstance(x, bool) or not isinstance(x, (float, int)) or not math.isfinite(x) or x <= 0
                for x in values):
            raise ValueError(name+' must contain positive finite values')
    for name in ('mode_count', 'matsubara_count'):
        if type(plan.get(name)) is not int or plan[name] < 1:
            raise ValueError(name+' must be a positive integer')
    if any(type(x) is not int or x < 1 for x in plan['orders']):
        raise ValueError('Quadrature orders must be positive integers')
    if not isinstance(plan.get('mode_indices'), list) or any(
            type(x) is not int or not 0 <= x < plan['mode_count'] for x in plan['mode_indices']):
        raise ValueError('Mode indices must address admitted eigenmodes')
    if not 0 < plan['T_K'] < plan['Tc_K']:
        raise ValueError('This uniform superconducting reference requires 0 < T < Tc')
    return plan


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=1)
    args = parser.parse_args()
    plan = validate_plan(json.loads(args.plan.read_text(encoding='utf8')))
    if args.output.exists():
        raise FileExistsError('Distinct output required; no overwrite')
    mesh = ROOT/plan['mesh']
    if hashlib.sha256(mesh.read_bytes()).hexdigest() != plan['mesh_sha256']:
        raise ValueError('Mesh hash mismatch')
    resources = linux_resources()
    parent = os.environ.get('PYSNSPD_SHARED_ALLOCATION')
    budget = admit_runtime(args.workers, resources, json.loads(parent) if parent else None)
    previous_affinity = os.sched_getaffinity(0)
    sources = source_hashes()
    if any(sources.get(name) != digest for name, digest in plan.get('source_sha256', {}).items()):
        raise ValueError('A plan-pinned source hash does not match the current implementation')
    args.output.mkdir(parents=True)
    started = time.monotonic()
    manifest = dict(status='RUNNING', plan=plan, plan_sha256=sha(args.plan), sources=sources,
        resources=resources, budget=budget, thread_environment={name:os.environ[name] for name in THREAD_VARIABLES},
        retry_policy='No automatic retries; independent energy queries continue after a failed query')
    atomic_json(args.output/'manifest.json', manifest)
    try:
        os.sched_setaffinity(0, {budget['coordinator_cpu']})
        return run_admitted(args, plan, mesh, budget, manifest, started)
    except BaseException as error:
        manifest.update(status='FAILED', exception=type(error).__name__, reason=str(error),
            elapsed_seconds=time.monotonic()-started, coordinator=runtime_telemetry())
        atomic_json(args.output/'manifest.json', manifest)
        raise
    finally:
        os.sched_setaffinity(0, previous_affinity)


def run_admitted(args, plan, mesh, budget, manifest, started):
    graph, modes, vectors, lift, conductance, eig_error = mesh_modes(mesh, plan['mode_count'])
    gap, _, _ = uniform_reference(graph, plan['T_K']/plan['Tc_K'], plan['matsubara_count'])
    np.savez_compressed(args.output/'spatial_modes.npz', modes=modes, vectors=vectors, lift=lift,
        coordinates=graph.coordinates_bar, area_weights=graph.area_weights)
    common = dict(plan, d=float(gap[0].real), t=plan['T_K']/plan['Tc_K'],
        tD_ps=hbar/(2*KB*plan['Tc_K'])*1e12, graph_conductance=conductance)
    tasks = []
    for omega in plan['omegas']:
        for eta_ratio in plan['eta_over_gap']:
            for order in plan['orders']:
                for cutoff in plan['cutoffs']:
                    base = dict(common, omega=omega, eta=eta_ratio*common['d'], order=order, cutoff=cutoff)
                    tasks.append(dict(base, kind='port'))
                    for number in plan['mode_indices']:
                        tasks.append(dict(base, kind='mode', mode=number, **{'lambda':float(modes[number])}))
    records = []
    workers = min(args.workers, len(tasks))
    manifest.update(d=common['d'], graph_conductance=conductance,
        eigen_residual=eig_error, tasks=len(tasks), workers=workers)
    atomic_json(args.output/'manifest.json', manifest)
    print(json.dumps(dict(event='START', tasks=len(tasks), workers=workers, nodes=graph.n_nodes)), flush=True)
    failures = []
    context = mp.get_context('spawn')
    counter = context.Value('i', 0)
    with ProcessPoolExecutor(max_workers=workers, mp_context=context,
            initializer=initialize_affinity, initargs=(budget['worker_affinity_cpus'], counter)) as pool:
        pending = {pool.submit(execute_task, task):index for index, task in enumerate(tasks)}
        while pending:
            ready, _ = wait(pending, timeout=5., return_when=FIRST_COMPLETED)
            for future in ready:
                index = pending.pop(future)
                try:
                    record = future.result()
                    atomic_json(args.output/f'query_{index:05d}.json', record)
                    records.append(dict(index=index, **record))
                except Exception as error:
                    failures.append(dict(index=index, reason=repr(error)))
                    atomic_json(args.output/f'query_{index:05d}.failed.json', failures[-1])
            completed = len(records)+len(failures)
            elapsed = time.monotonic()-started
            print(json.dumps(dict(event='PROGRESS', completed=completed, total=len(tasks),
                fraction=completed/len(tasks), elapsed_seconds=elapsed,
                eta_seconds=elapsed*(len(tasks)-completed)/completed if completed else None,
                failures=len(failures))), flush=True)
            manifest.update(completed_queries=completed, failed_queries=len(failures), elapsed_seconds=elapsed)
            atomic_json(args.output/'manifest.json', manifest)
    status = 'COMPLETED' if not failures else 'COMPLETED_WITH_FAILURES'
    atomic_json(args.output/'results.json', dict(status=status,
        records=sorted(records, key=lambda row:row['index']), failures=failures,
        elapsed_seconds=time.monotonic()-started, scope=plan['scope']))
    manifest.update(status=status, elapsed_seconds=time.monotonic()-started, coordinator=runtime_telemetry(),
        outputs={name:sha(args.output/name) for name in ('spatial_modes.npz', 'results.json')})
    atomic_json(args.output/'manifest.json', manifest)
    print(json.dumps(dict(event='COMPLETE', status=status, fraction=1., eta_seconds=0.)), flush=True)
    return 1 if failures else 0


if __name__ == '__main__':
    raise SystemExit(main())
