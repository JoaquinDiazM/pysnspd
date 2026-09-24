"""User-run weak dynamic coupling on a self-consistent biased 2D strip.

The gap and potential responses use a declared spatial Galerkin basis; every
spectral and kinetic solve retains all mesh nodes. Basis, quadrature, cutoff
and regulator comparisons are outputs, not silently repaired acceptance.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor, FIRST_COMPLETED, wait
from contextlib import contextmanager
from functools import lru_cache
import hashlib
import json
import multiprocessing as mp
import os
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from sandbox.stage4_core.parallel_runtime import limit_thread_environment, linux_resources, initialize_affinity
limit_thread_environment()
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.constants import hbar, k as KB, e as EC
from pysnspd.experimental.thermal_spatial_usadel import ThermalGraph
from pysnspd.experimental.electrical_ports import ThesisCircuitParameters
from pysnspd.experimental.thesis_circuit_harmonic import solve_harmonic
from sandbox.stage4_core.coupled_response import atomic_json, encode, mesh_modes, admit_runtime, runtime_telemetry


@contextmanager
def durable_run(output, manifest, budget):
    """All preparation, worker and postprocessing failures receive a final status."""
    original = os.sched_getaffinity(0)
    started = time.monotonic()
    try:
        os.sched_setaffinity(0, {budget['coordinator_cpu']})
        yield
    except BaseException as error:
        receipt = dict(status='FAILED', exception=type(error).__name__, reason=str(error),
            phase=manifest.get('phase'), elapsed_seconds=time.monotonic()-started,
            policy='Completed evidence preserved; no automatic retry')
        if 'failure' not in manifest:
            manifest['failure'] = receipt
            atomic_json(Path(output)/'failure.json', receipt)
        manifest.update(status='FAILED', elapsed_seconds=time.monotonic()-started)
        atomic_json(Path(output)/'manifest.json', manifest)
        print(json.dumps(dict(event='FAILED', **receipt)), flush=True)
        raise
    finally:
        os.sched_setaffinity(0, original)


class PhaseProgress:
    """Query progress with phase-specific cost estimates and bounded output."""
    def __init__(self, phase_counts, emit=print, clock=time.monotonic, interval=5.):
        self.counts = dict(phase_counts)
        self.done = {name:0 for name in self.counts}
        self.emit, self.clock, self.interval = emit, clock, interval
        self.started = clock()
        self.phase_started = None
        self.phase = None
        self.last_print = -float('inf')

    def begin(self, phase):
        self.phase, self.phase_started = phase, self.clock()
        self.report(force=True)

    def completed(self, count=1):
        self.done[self.phase] += count
        self.report(force=self.done[self.phase] == self.counts[self.phase])

    def report(self, *, force=False):
        now = self.clock()
        if not force and now-self.last_print < self.interval:
            return
        done, total = self.done[self.phase], self.counts[self.phase]
        phase_elapsed = now-self.phase_started
        phase_eta = phase_elapsed*(total-done)/done if done else None
        final_phase = self.phase == next(reversed(self.counts))
        payload = dict(event='PROGRESS', phase=self.phase, phase_completed=done,
            phase_total=total, completed=sum(self.done.values()), total=sum(self.counts.values()),
            fraction=sum(self.done.values())/sum(self.counts.values()),
            fraction_meaning='Fraction of declared queries, not fraction of runtime',
            elapsed_seconds=now-self.started, phase_elapsed_seconds=phase_elapsed,
            phase_eta_seconds=phase_eta, eta_seconds=phase_eta if final_phase else None,
            eta_basis='Observed mean cost within the current phase; energy-query cost can vary with energy')
        self.emit(json.dumps(payload, allow_nan=False), flush=True)
        self.last_print = now


def completed_futures(jobs, progress):
    """Keep printing a bounded heartbeat while a costly query is still running."""
    pending = set(jobs)
    while pending:
        ready, pending = wait(pending, timeout=5., return_when=FIRST_COMPLETED)
        for future in ready:
            yield future
        if not ready:
            progress.report()


def admitted_reference(path):
    """A saved final iterate is not necessarily an admitted stationary branch."""
    path = Path(path)
    receipt_path = path.parent/'manifest.json'
    receipt = json.loads(receipt_path.read_text(encoding='utf8'))
    if receipt.get('status') != 'FINITE_SUM_STATIONARY':
        raise ValueError('The stationary reference did not reach its declared gap criterion')
    final = receipt.get('latest_checkpoint') or {}
    if final.get('sha256') != hashlib.sha256(path.read_bytes()).hexdigest():
        raise ValueError('Reference bytes differ from the stationary acceptance receipt')
    if (path.parent/final.get('path','')).resolve() != path.resolve():
        raise ValueError('The reference path differs from the final admitted checkpoint')
    return receipt


@contextmanager
def fail_fast_pool(failure_directory, manifest, **executor_arguments):
    """Stop this campaign's workers after failure instead of draining the queue.

    A standard ProcessPoolExecutor context waits for all submitted jobs on
    exception, including a large queued energy campaign. Only child processes
    owned by this executor are terminated here; no external screen or process
    name is targeted. Completed files are preserved and no task is retried.
    """
    pool = ProcessPoolExecutor(**executor_arguments)
    try:
        yield pool
    except BaseException as error:
        # Older supported Python versions lack terminate_workers(). Capture
        # owned processes before shutdown clears the executor's process map.
        owned = list((getattr(pool, '_processes', None) or {}).values())
        terminated = []
        for process in owned:
            if process.is_alive():
                terminated.append(process.pid)
                process.terminate()
        deadline = time.monotonic()+5.
        for process in owned:
            process.join(timeout=max(0., deadline-time.monotonic()))
        killed = []
        for process in owned:
            if process.is_alive():
                killed.append(process.pid)
                process.kill()
        pool.shutdown(wait=False, cancel_futures=True)
        receipt = dict(status='FAILED', exception=type(error).__name__, reason=str(error),
            active_query=manifest.get('active_query'), terminated_owned_worker_pids=terminated,
            killed_owned_worker_pids=killed,
            policy='Queued work cancelled; owned workers terminated; completed evidence preserved; no retry')
        manifest.update(status='FAILED', failure=receipt)
        atomic_json(Path(failure_directory)/'failure.json', receipt)
        atomic_json(Path(failure_directory)/'manifest.json', manifest)
        print(json.dumps(dict(event='FAILED', **receipt)), flush=True)
        raise
    else:
        pool.shutdown(wait=True)


@lru_cache(maxsize=2)
def reference(path):
    with np.load(path) as data:
        arrays = {key:data[key].copy() for key in data.files}
    d = arrays.get('d', arrays.get('delta_bar'))
    epsilon = arrays.get('epsilon', arrays.get('epsilon_bar'))
    graph = ThermalGraph(arrays['area_weights'], arrays['edges'], arrays['conductance'],
        arrays['coordinates_bar'], arrays['boundary_nodes'])
    return arrays, graph, d, epsilon


@lru_cache(maxsize=2)
def load_modes(path):
    with np.load(path) as data:
        return data['modes'].copy(), data['lift'].copy()


def thermal_query(payload):
    from pysnspd.experimental.thermal_contact_tangent import ThermalContactTangent
    arrays, graph, d, epsilon = reference(payload['reference_npz'])
    modes, lift = load_modes(payload['modes_npz'])
    count = modes.shape[1]
    phase = d/abs(d)
    directions = np.column_stack((phase[:, None]*modes, 1j*phase[:, None]*modes, 1j*d*lift))
    index = payload['n']
    factor = ThermalContactTangent(graph, d, epsilon[index], arrays['u'][index], alpha=arrays.get('alpha'))
    force = np.empty(directions.shape, complex)
    current = np.empty((len(graph.edges), 2*count+1))
    residual = 0.
    for j in range(directions.shape[1]):
        response = factor.apply(directions[:, j])
        force[:, j] = 4*np.pi*payload['t']*graph.area_weights*(directions[:, j]/epsilon[index]-response.df)
        current[:, j] = 2*np.pi*payload['t']*response.current_derivative
        residual = max(residual, response.equation_residual)
    rotated = phase.conj()[:, None]*force
    return np.stack((rotated.real, rotated.imag), axis=1), current, residual


def energy_rule(plan, d):
    """Dense physical-energy panels also cover the depairing-shifted gap edge."""
    cutoff = float(plan['energy_cutoff'])
    core = max(4., 2*float(np.max(abs(d)))+plan['omega'])
    if cutoff <= core:
        raise ValueError('Energy cutoff must exceed the resolved core interval')
    points = np.r_[np.linspace(-core, core, int(np.ceil(2*core/plan['maximum_energy_panel']))+1),
        -np.geomspace(core, cutoff, 14), np.geomspace(core, cutoff, 14)]
    points = np.unique(points)
    x, w = leggauss(plan['quadrature_order'])
    scale = np.diff(points)/2
    energy = ((points[:-1, None]+points[1:, None])/2+scale[:, None]*x).ravel()
    weights = (scale[:, None]*w).ravel()
    tail_order = int(plan.get('energy_tail_order', 0))
    if tail_order:
        # Exact variable substitution E=C/x on (0,1), not an inferred balance
        # correction. Both signed semi-infinite intervals are integrated.
        tx, tw = leggauss(tail_order)
        tx, tw = (tx+1)/2, tw/2
        tail_energy, tail_weights = cutoff/tx, cutoff*tw/(tx*tx)
        energy = np.r_[energy, tail_energy, -tail_energy]
        weights = np.r_[weights, tail_weights, tail_weights]
    return energy, weights


def port_sum(graph, currents):
    coordinates = graph.coordinates_bar
    midpoint = .5*(coordinates[:, 0].min()+coordinates[:, 0].max())
    left = np.zeros(graph.n_nodes, bool)
    left[graph.boundary_nodes] = coordinates[graph.boundary_nodes, 0] < midpoint
    tail, head = graph.edges.T
    return (left[tail].astype(int)-left[head].astype(int))@currents


def divergence(graph, currents):
    result = np.zeros((graph.n_nodes, currents.shape[1]), complex)
    np.add.at(result, graph.edges[:, 0], currents)
    np.add.at(result, graph.edges[:, 1], -currents)
    return result


def summarize(plan, graph, d, modes, lift, thermal_force, thermal_current,
              anomaly, charge, current_correction, output):
    k, n = modes.shape[1], graph.n_nodes
    mass = graph.area_weights
    t = plan['T_K']/plan['Tc_K']
    omega = plan['omega']
    # Thermal arrays contain physical real Cartesian directions; the port
    # phase amplitude is -2i/Î© times the last real direction i*d*lift.
    tf = np.zeros((n, 2, 3*k+1), complex)
    ti = np.zeros((len(graph.edges), 3*k+1), complex)
    tf[:, :, :2*k], ti[:, :2*k] = thermal_force[:, :, :2*k], thermal_current[:, :2*k]
    tf[:, :, -1], ti[:, -1] = (-2j/omega)*thermal_force[:, :, -1], (-2j/omega)*thermal_current[:, -1]
    total_force, total_current = tf+anomaly, ti+current_correction
    potential = np.column_stack((np.zeros((n, 2*k)), modes, lift))
    neutral = charge+potential
    matrix = np.vstack((modes.T@total_force[:, 0], modes.T@total_force[:, 1],
                        modes.T@(mass[:, None]*neutral)))
    tD_ps = hbar/(2*KB*plan['Tc_K'])*1e12
    tau = 1/(t/plan['tau_ee_Tc_ps']+t**3/plan['tau_ep_Tc_ps'])
    gamma = 2*np.sqrt((1+t)/2)*np.pi/4*np.sqrt(1+(tau/tD_ps)**2*abs(d)**2)
    gamma_projected = modes.T@((mass*gamma)[:, None]*modes)
    microscopic = matrix[:, :-1].copy()
    # Diagnose the radial kernel only AFTER eliminating phase and charge.
    schur = microscopic[:k, :k]-microscopic[:k, k:]@np.linalg.solve(microscopic[k:, k:], microscopic[k:, :k])
    gamma_resolved = (schur.conj().T-schur)/(1j*omega)
    gamma_residual = gamma_projected-gamma_resolved
    # Candidate adds the inherited radial dissipation once. It is explicitly
    # provisional until the resolved/residual mobility comparison is admitted.
    candidate = microscopic.copy()
    candidate[:k, :k] -= .5j*omega*gamma_projected
    coefficients = np.linalg.solve(candidate, -matrix[:, -1])
    all_coefficients = np.r_[coefficients, 1.]
    radial = modes@coefficients[:k]
    angular = modes@coefficients[k:2*k]+(-2j/omega)*abs(d)*lift
    voltage = potential@all_coefficients
    neutral_value = neutral@all_coefficients
    force_value = np.einsum('ncf,f->nc', total_force, all_coefficients)
    force_value[:, 0] -= .5j*omega*mass*gamma*radial
    currents = total_current@all_coefficients
    div = divergence(graph, currents[:, None])[:, 0]
    free = np.ones(n, bool); free[graph.boundary_nodes] = False
    def rms(value):
        return float(np.sqrt(np.sum(mass[free]*abs(value[free])**2)/mass[free].sum()))
    current_scale = max(float(np.max(abs(currents))), 1e-12)
    full_projection = dict(
        radial_force_rms=rms(force_value[:, 0]/mass),
        angular_force_rms=rms(force_value[:, 1]/mass),
        neutral_rms=rms(neutral_value),
        voltage_rms=rms(voltage),
        neutral_relative_to_voltage=rms(neutral_value)/max(rms(voltage), 1e-12),
        current_divergence_rms=rms(div/mass),
        current_divergence_over_link_current=rms(div/mass)/current_scale,
        projected_equation_max=float(np.max(abs(candidate@coefficients+matrix[:, -1]))))
    Ibar = port_sum(graph, currents)
    E0 = KB*plan['Tc_K']
    current_unit = E0/(2*EC*plan['sheet_resistance_ohm'])
    admittance = (EC/E0)*current_unit*Ibar
    # Equilibrated differential inductance at the same branch and port planes.
    hs = np.vstack((modes.T@thermal_force[:, 0], modes.T@thermal_force[:, 1]))
    static_coefficients = np.linalg.solve(hs[:, :-1], -hs[:, -1])
    static_current = thermal_current@np.r_[static_coefficients, 1.]
    derivative = -current_unit*port_sum(graph, static_current)  # lift: thetaR-thetaL=-1
    resolved_L = hbar/(2*EC)/derivative
    external_L = plan['total_inductance_H']-resolved_L
    if not np.isfinite(resolved_L) or resolved_L <= 0 or external_L <= 0:
        raise ValueError('No passive positive exterior inductance at this branch; no clipping')
    dc_current = plan.get('reference_dc_current_A')
    dc_supply = 0. if dc_current is None else 1e4*dc_current
    circuit = solve_harmonic(ThesisCircuitParameters(Lk_ext_H=float(external_L), V_bias_V=dc_supply),
        omega*E0/hbar, 1/admittance, bias_voltage_peak_V=plan['bias_voltage_peak_V'], time_convention='exp(-iwt)')
    minimum_stiffness = float(np.linalg.eigvalsh((hs[:, :-1]+hs[:, :-1].T)/2).min())
    np.savez_compressed(output/'coupled_fields.npz', gap_reference=d, modes=modes,
        radial_response=radial, angular_response=angular, potential_response=voltage,
        neutral_residual=neutral_value, current_response=currents,
        force_residual=force_value, coefficients=coefficients, microscopic_matrix=microscopic,
        candidate_matrix=candidate, thermal_hessian=hs[:, :-1], gamma_resolved=gamma_resolved,
        gamma_residual=gamma_residual, coordinates_bar=graph.coordinates_bar,
        area_weights=mass, edges=graph.edges, lift=lift)
    radial_heat_per_voltage_squared = omega**2/(16*plan['sheet_resistance_ohm'])*float(np.sum(mass*gamma*abs(radial)**2))
    port_power_per_voltage_squared = float(admittance.real/2)
    return dict(status='COMPUTED_REQUIRES_COMPARISON', scope=plan['scope'],
        omega=omega, angular_frequency_per_ps=omega/(2*tD_ps),
        mode_count=k, minimum_static_projected_stiffness=minimum_stiffness,
        minimum_residual_mobility_eigenvalue=float(np.linalg.eigvalsh(gamma_residual).min()),
        resolved_mobility_norm_over_kwt=float(np.linalg.norm(gamma_resolved, 2)/np.linalg.norm(gamma_projected, 2)),
        full_field_residuals=full_projection, admittance_S=encode(admittance),
        impedance_ohm=encode(1/admittance), resolved_inductance_H=float(resolved_L),
        external_inductance_H=float(external_L), bias_voltage_peak_V=plan['bias_voltage_peak_V'],
        circuit_state_peak=encode(circuit.state_peak), readout_peak_V=encode(circuit.readout_voltage_peak_V),
        circuit_power_residual_W=circuit.average_power_residual_W,
        reference_dc_current_A=dc_current, reference_dc_supply_voltage_V=dc_supply,
        film_port_power_per_voltage_squared_W=port_power_per_voltage_squared,
        candidate_radial_heat_per_voltage_squared_W=radial_heat_per_voltage_squared,
        port_power_minus_radial_heat_per_voltage_squared_W=port_power_per_voltage_squared-radial_heat_per_voltage_squared,
        energy_scope='AC port power and positive radial KWT heat are separate observables. Their difference is not an independently computed reservoir heat flux; this is not a nonlinear energy certificate.',
        circuit_time_convention='exp(-iwt)',
        acceptance='Compare independent refinements with the declared 2% observable margin; no automatic nonlinear closure from this result.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--workers', type=int, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text(encoding='utf8'))
    if plan.get('schema') != 'pysnspd.stage4.biased_coupled_response.v1':
        raise ValueError('Unsupported plan')
    if args.output.exists():
        raise FileExistsError('A fresh output directory is required')
    inherited = os.environ.get('PYSNSPD_SHARED_ALLOCATION')
    budget = admit_runtime(args.workers, linux_resources(), json.loads(inherited) if inherited else None)
    mesh_path = ROOT/plan['mesh']
    if hashlib.sha256(mesh_path.read_bytes()).hexdigest() != plan['mesh_sha256']:
        raise ValueError('Mesh hash changed')
    reference_receipt = admitted_reference(args.reference)
    for key in ('T_K', 'Tc_K'):
        if reference_receipt['plan'][key] != plan[key]:
            raise ValueError('Reference and response must share '+key)
    args.output.mkdir(parents=True)
    reference_path = str(args.reference.resolve())
    manifest = dict(status='RUNNING', plan=plan, reference=reference_path,
        reference_sha256=hashlib.sha256(args.reference.read_bytes()).hexdigest(),
        workers=args.workers, budget=budget, phase='preparation',
        reference_receipt_sha256=hashlib.sha256((args.reference.parent/'manifest.json').read_bytes()).hexdigest(),
        source_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()
            for path in [Path(__file__).resolve(), ROOT/'sandbox/stage4_core/biased_harmonic_query.py',
                ROOT/'pysnspd/experimental/harmonic_usadel.py', ROOT/'pysnspd/experimental/harmonic_kinetic_usadel.py',
                ROOT/'pysnspd/experimental/thermal_contact_tangent.py']})
    atomic_json(args.output/'manifest.json', manifest)
    started = time.monotonic()
    with durable_run(args.output, manifest, budget):
        execute_admitted(args, plan, budget, manifest, started, mesh_path, reference_path)


def execute_admitted(args, plan, budget, manifest, started, mesh_path, reference_path):
    arrays, graph, d, epsilon = reference(reference_path)
    if 'current_bar' in arrays:
        plan['reference_dc_current_A'] = float(KB*plan['Tc_K']/(2*EC*plan['sheet_resistance_ohm'])*port_sum(graph, arrays['current_bar']))
    mesh_graph, eigenvalues, modes, lift, _, eigen_residual = mesh_modes(mesh_path, plan['mode_count'])
    if (mesh_graph.signature != graph.signature
            or not np.array_equal(mesh_graph.boundary_nodes,graph.boundary_nodes)
            or not np.array_equal(mesh_graph.coordinates_bar,graph.coordinates_bar)):
        raise ValueError('Stationary reference and response must use the same graph')
    mode_path = str((args.output/'modes.npz').resolve())
    np.savez_compressed(mode_path, modes=modes, lift=lift, eigenvalues=eigenvalues)
    workers = args.workers
    # Admission checked the full assigned set before pinning this coordinator.
    manifest.update(eigen_residual=eigen_residual, phase='thermal_tangent')
    atomic_json(args.output/'manifest.json', manifest)
    t = plan['T_K']/plan['Tc_K']
    phase = d/abs(d)
    directions = np.column_stack((phase[:, None]*modes, 1j*phase[:, None]*modes, 1j*d*lift))
    base_force = 2*graph.area_weights[:, None]*directions*np.log(t)
    rotated = phase.conj()[:, None]*base_force
    thermal_force = np.stack((rotated.real, rotated.imag), axis=1)
    thermal_current = np.zeros((len(graph.edges), directions.shape[1]))
    common = dict(reference_npz=reference_path, modes_npz=mode_path, t=t)
    thermal_errors = []
    from sandbox.stage4_core.biased_harmonic_query import energy_query
    energy, weights = energy_rule(plan, d)
    completed = 0
    progress = PhaseProgress({'thermal_tangent':len(epsilon), 'coupled_spectral_kinetic':len(energy)})
    anomaly = np.zeros((graph.n_nodes, 2, 3*plan['mode_count']+1), complex)
    charge = np.zeros((graph.n_nodes, 3*plan['mode_count']+1), complex)
    current = np.zeros((len(graph.edges), 3*plan['mode_count']+1), complex)
    context = mp.get_context('spawn')
    counter = context.Value('i', 0)
    with fail_fast_pool(args.output, manifest, max_workers=workers, mp_context=context,
            initializer=initialize_affinity, initargs=(budget['worker_affinity_cpus'], counter)) as pool:
        progress.begin('thermal_tangent')
        jobs = {pool.submit(thermal_query, dict(common, n=n)):n for n in range(len(epsilon))}
        for future in completed_futures(jobs, progress):
            manifest['active_query'] = dict(phase='thermal_tangent', index=jobs.pop(future))
            force, flux, error = future.result()
            thermal_force += force; thermal_current += flux; thermal_errors.append(error)
            completed += 1; progress.completed()
        np.savez_compressed(args.output/'thermal_tangent.npz', force=thermal_force, current=thermal_current)
        manifest.update(phase='coupled_spectral_kinetic', active_query=None)
        atomic_json(args.output/'manifest.json', manifest)
        progress.begin('coupled_spectral_kinetic')
        jobs = {}
        for index, value in enumerate(energy):
            payload = dict(common, E=float(value), energy=float(value), omega=plan['omega'],
                eta=plan['eta_over_gap']*float(np.max(abs(d))),
                spectral_tolerance=plan['spectral_tolerance'], continuation_steps=plan['continuation_steps'])
            jobs[pool.submit(energy_query, payload)] = index
        # Every completed query is independently checkpointed. No retry or
        # accidental continuation after a failed physical query is performed.
        for future in completed_futures(jobs, progress):
            index = jobs.pop(future)
            manifest['active_query'] = dict(phase='coupled_spectral_kinetic', index=index, energy=float(energy[index]))
            result = future.result()
            weight = weights[index]/2
            anomaly += weight*result['anomaly_correction']
            charge += weight*result['charge']
            current += weight*result['current_correction']
            metadata = {**result['metadata'], 'index':index, 'energy':float(energy[index]), 'weight':float(weights[index])}
            atomic_json(args.output/f'energy_{index:05d}.json', metadata)
            # Aggregates suffice for reconstruction of this declared quadrature.
            # Individual spectra are not claimed as restart checkpoints.
            completed += 1; progress.completed()
            if completed % 32 == 0:
                np.savez_compressed(args.output/'partial_moments.npz', anomaly=anomaly, charge=charge,
                    current=current, completed=completed)
    manifest.update(phase='postprocessing', active_query=None)
    atomic_json(args.output/'manifest.json', manifest)
    np.savez_compressed(args.output/'integrated_moments.npz', anomaly=anomaly, charge=charge, current=current)
    result = summarize(plan, graph, d, modes, lift, thermal_force, thermal_current, anomaly, charge, current, args.output)
    result.update(elapsed_seconds=time.monotonic()-started, energy_nodes=len(energy),
        maximum_thermal_tangent_residual=max(thermal_errors))
    atomic_json(args.output/'results.json', result)
    manifest.update(status='COMPLETED', phase='complete', elapsed_seconds=result['elapsed_seconds'], coordinator=runtime_telemetry(),
        outputs={name:hashlib.sha256((args.output/name).read_bytes()).hexdigest()
            for name in ('results.json','coupled_fields.npz','integrated_moments.npz','thermal_tangent.npz','modes.npz')})
    atomic_json(args.output/'manifest.json', manifest)
    print(json.dumps(dict(event='COMPLETE', elapsed_seconds=result['elapsed_seconds'], status=result['status'], fraction=1., eta_seconds=0.)), flush=True)


if __name__ == '__main__':
    main()
