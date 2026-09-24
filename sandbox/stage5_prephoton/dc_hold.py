"""Preserve a current-carrying bulk-matched thermal state before any photon.

Every accepted step uses the inherited first-order local KWT update and the
finite Matsubara force/current on every free mesh node.  The initial force is
never subtracted.  Spectral boundary values come from the supplied biased
reference, not from the zero-current expression u=Delta/epsilon.

This is a DC preservation diagnostic, not a nonequilibrium photon solver.
Thermal populations are prescribed equilibrium distributions; the normal
potential uses the inherited Ohmic closure and equipotential current ports.
The full thesis circuit supplies its current and receives the passive strip
voltage; its external inductance is constant. Bulk amplitudes remain fixed at
the two end planes while their phases evolve by the Josephson relation.
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
from pathlib import Path
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from sandbox.stage4_core.parallel_runtime import limit_thread_environment, linux_resources
limit_thread_environment()
import numpy as np
from scipy.constants import Boltzmann, elementary_charge, hbar
from scipy.sparse.linalg import splu
from pysnspd.experimental import thermal_spatial_usadel as thermal
from pysnspd.experimental.thermal_stable_newton import solve_frequency
from pysnspd.experimental.stable_kwt_euler import stable_inherited_kwt_step
from pysnspd.experimental.electrical_ports import ThesisCircuitParameters
from pysnspd.experimental.dc_current_ports import (DCCurrentPorts,
    CurrentDrivenThermalKWT, advance_bulk_contact_phases)
from sandbox.stage4_core.coupled_campaign import atomic_json, sha
from sandbox.stage4_core.coupled_response import admit_runtime
from sandbox.stage5_prephoton.provenance import source_sha

SCHEMA = 'pysnspd.stage5.prephoton_dc_hold.v1'


def validate_plan(plan):
    if plan.get('schema') != SCHEMA:
        raise ValueError('Unsupported DC hold schema')
    for name in ('T_K', 'Tc_K', 'tau_ee_Tc_ps', 'tau_ep_Tc_ps', 'dt_ps',
                 'duration_ps', 'spectral_tolerance', 'delta0_over_kBTc',
                 'sheet_resistance_ohm', 'gap_drift_relative_limit',
                 'gap_ripple_relative_limit', 'current_cut_relative_limit',
                 'current_absolute_floor_A', 'maximum_relative_displacement'):
        value = plan.get(name)
        if isinstance(value, bool) or not isinstance(value, (float, int)) or not np.isfinite(value) or value <= 0:
            raise ValueError(name+' must be positive and finite')
    if not plan['T_K'] < plan['Tc_K']:
        raise ValueError('Require T < Tc')
    for name in ('maximum_newton_iterations', 'observation_cadence'):
        if type(plan.get(name)) is not int or plan[name] < 1:
            raise ValueError(name+' must be a positive integer')
    if not np.isfinite(plan.get('bias_current_A', np.nan)):
        raise ValueError('An explicit finite DC bias current in ampere is required')
    circuit = plan.get('circuit', {})
    for key in ('R_load_ohm', 'R_bias_ohm', 'L_bias_H', 'Lk_ext_H', 'C_couple_F'):
        if key not in circuit:
            raise ValueError('Explicit constant circuit component required: '+key)
    ThesisCircuitParameters(**circuit)
    for name, default in (('dark_readout_fraction_limit',1e-3),('energy_balance_relative_limit',1e-2)):
        value=plan.get(name,default)
        if isinstance(value,bool) or not isinstance(value,(int,float)) or not np.isfinite(value) or value<=0:
            raise ValueError(name+' must be positive and finite')
    if 'V_bias_V' in circuit:
        raise ValueError('Do not set V_bias_V: the runner derives R_bias times the initial discrete strip current')
    return plan


def load_reference(path):
    with np.load(path, allow_pickle=False) as data:
        fixed_key = 'boundary_nodes' if 'boundary_nodes' in data else 'fixed_nodes'
        graph = thermal.ThermalGraph(data['area_weights'], data['edges'],
            data['conductance'], data['coordinates_bar'], data[fixed_key])
        gap, epsilon, u = (data[key].copy() for key in ('delta_bar', 'epsilon_bar', 'u'))
        if 'alpha' in data and np.any(data['alpha'] != 0):
            raise ValueError('Stable-energy Newton backend requires alpha=0; use the physical q*x phase gauge')
    if gap.shape != (graph.n_nodes,) or u.shape != (len(epsilon), graph.n_nodes):
        raise ValueError('Reference gap and frequency-resolved u shapes do not match the graph')
    if any(np.any(~np.isfinite(value)) for value in (gap, epsilon, u)):
        raise ValueError('Reference fields must be finite')
    x = graph.coordinates_bar[:, 0]
    ends = np.flatnonzero(np.isclose(x, x.min(), atol=1e-10, rtol=0) |
                         np.isclose(x, x.max(), atol=1e-10, rtol=0))
    if not np.array_equal(np.sort(graph.boundary_nodes), ends):
        raise ValueError('Only longitudinal ends may be fixed; lateral edges must remain natural')
    return graph, gap, epsilon, u


def cut_currents(graph, current, count=9):
    """Oriented flux through nine interior x cuts; no nodal interpolation."""
    x = graph.coordinates_bar[:, 0]
    locations = np.linspace(x.min(), x.max(), count+2)[1:-1]
    tail, head = graph.edges.T
    values = []
    for cut in locations:
        signs = ((x[tail] < cut) & (x[head] >= cut)).astype(float)
        signs -= ((x[head] < cut) & (x[tail] >= cut))
        values.append(float(np.dot(signs, current)))
    return locations, np.asarray(values)


def _mode_contribution(graph, gap, epsilon, u, fixed_u, plan, factor_cache=None):
    """One exact spectral contribution, with an optional resident chord factor.

    A previous sparse Jacobian is only a predictor/preconditioner. The proposed
    corrected field is admitted by the *nonlinear* residual at the new gap and
    new reservoir phases, using the same tolerance as full Newton. A rejected
    predictor invokes the unchanged Newton solver and refreshes the factor.
    No spectral value, force or energy is interpolated or subtracted.
    """
    fixed = graph.boundary_nodes
    free = np.ones(graph.n_nodes, bool); free[fixed] = False
    guess = u.copy(); guess[fixed] = fixed_u
    value = thermal.spectral_energy_gradient(graph, gap, epsilon, guess)
    scale=graph.area_weights[free]*np.maximum(1., abs(gap[free]))
    def norm(evaluation):
        return float(np.max(abs(evaluation.residual[free])/scale,initial=0.))
    residual=norm(value)
    reused=residual<=plan['spectral_tolerance']
    corrections=factorizations=newton_solves=0
    components=np.flatnonzero(np.repeat(free,2))

    def refresh(at):
        _,jacobian=thermal.spectral_residual_jacobian(graph,gap,epsilon,at)
        factor_cache['factor']=splu(jacobian[components][:,components].tocsc())

    if not reused and factor_cache is not None and len(components):
        if 'factor' not in factor_cache:
            refresh(guess);factorizations+=1
        rhs=np.column_stack((value.residual.real,value.residual.imag)).ravel()
        step=factor_cache['factor'].solve(-rhs[components])
        trial=guess.copy()
        trial[free]+=step[::2]+1j*step[1::2]
        if np.all(np.isfinite(trial)):
            proposed=thermal.spectral_energy_gradient(graph,gap,epsilon,trial)
            proposed_residual=norm(proposed)
            if proposed_residual<=plan['spectral_tolerance']:
                guess,value,residual=trial,proposed,proposed_residual
                corrections=1
    if residual>plan['spectral_tolerance']:
        solution = solve_frequency(graph, gap, epsilon, initial_u=guess,
            fixed_nodes=fixed, fixed_u=fixed_u, tol=plan['spectral_tolerance'],
            max_iterations=plan['maximum_newton_iterations'])
        newton_solves=1
        guess, residual = solution.u, solution.residual
        value = thermal.spectral_energy_gradient(graph, gap, epsilon, guess)
        if factor_cache is not None and len(components):
            refresh(guess);factorizations+=1
    mass = graph.area_weights
    return guess, dict(force=2*mass*(gap/epsilon-value.f),
        current=-value.link_derivative_alpha,
        energy=float(value.energy+np.dot(mass, abs(gap)**2)/epsilon),
        spectral_torque=np.imag(np.conj(guess)*value.gradient_u),
        residual=residual, newton_solves=newton_solves, reuses=int(reused),
        chord_corrections=corrections,spectral_factorizations=factorizations)


class SerialSpectralEvaluator:
    """Small reproducible control; production campaign uses resident workers."""
    def __init__(self, graph, epsilon, initial_u, plan):
        self.graph, self.epsilon, self.plan = graph, epsilon, plan
        self.u = initial_u.copy()
        self.fixed_u = initial_u[:, graph.boundary_nodes].copy()
        self.fixed_phase = self.fixed_u[0]/abs(self.fixed_u[0])
        self.newton_solves = self.reuses = self.batches = 0
        self.chord_corrections=self.spectral_factorizations=0
        self.factor_cache=[{} for _ in epsilon]
        self.maximum_residual = 0.

    def __call__(self, gap):
        parts = []
        phase_ratio = gap[self.graph.boundary_nodes]/abs(gap[self.graph.boundary_nodes])/self.fixed_phase
        for n, epsilon in enumerate(self.epsilon):
            self.u[n], part = _mode_contribution(self.graph, gap, epsilon,
                self.u[n], self.fixed_u[n]*phase_ratio, self.plan,self.factor_cache[n])
            parts.append(part)
        return _aggregate(self, gap, parts)


def _aggregate(evaluator, gap, parts):
    graph, plan = evaluator.graph, evaluator.plan
    ratio = plan['T_K']/plan['Tc_K']; weight = 2*np.pi*ratio
    force = 2*graph.area_weights*gap*np.log(ratio)
    current = np.zeros(len(graph.edges))
    torque = np.zeros(graph.n_nodes)
    energy = float(np.dot(graph.area_weights, abs(gap)**2)*np.log(ratio))
    for part in parts:
        force += weight*part['force']; current += weight*part['current']
        torque += weight*part['spectral_torque']
        energy += weight*part['energy']
        evaluator.newton_solves += part['newton_solves']
        evaluator.reuses += part['reuses']
        evaluator.chord_corrections+=part['chord_corrections']
        evaluator.spectral_factorizations+=part['spectral_factorizations']
        evaluator.maximum_residual = max(evaluator.maximum_residual, part['residual'])
    evaluator.batches += 1
    return dict(gradient=force, current=current, energy=energy, spectral_torque=torque)


def _actor(pipe, cpu, graph, epsilon, initial_u, indices, plan):
    try:
        os.sched_setaffinity(0, {cpu})
        cache = {n: initial_u[n].copy() for n in indices}
        factor_cache={n:{} for n in indices}
        fixed_u = {n: initial_u[n, graph.boundary_nodes].copy() for n in indices}
        fixed_phase = initial_u[0, graph.boundary_nodes]/abs(initial_u[0, graph.boundary_nodes])
        pipe.send(dict(event='READY', cpu=cpu, frequencies=len(indices)))
        while True:
            command, gap = pipe.recv()
            if command == 'STOP':
                break
            if command != 'EVALUATE':
                raise ValueError('Unexpected resident spectral command')
            total = dict(force=np.zeros(graph.n_nodes, complex),
                current=np.zeros(len(graph.edges)), energy=0., residual=0.,
                spectral_torque=np.zeros(graph.n_nodes),
                newton_solves=0, reuses=0,chord_corrections=0,spectral_factorizations=0)
            phase_ratio = gap[graph.boundary_nodes]/abs(gap[graph.boundary_nodes])/fixed_phase
            for n in indices:
                cache[n], part = _mode_contribution(graph, gap, epsilon[n],
                    cache[n], fixed_u[n]*phase_ratio, plan,factor_cache[n])
                for name in ('force', 'current', 'energy', 'spectral_torque', 'newton_solves',
                             'reuses','chord_corrections','spectral_factorizations'):
                    total[name] += part[name]
                total['residual'] = max(total['residual'], part['residual'])
            pipe.send(total)
    except BaseException as exc:
        try:
            pipe.send(dict(event='ERROR', reason=str(exc), traceback=traceback.format_exc()))
        except (BrokenPipeError, EOFError):
            pass
    finally:
        pipe.close()


class ResidentSpectralPool:
    def __init__(self, graph, epsilon, initial_u, plan, budget, event):
        self.graph, self.plan = graph, plan
        self.newton_solves = self.reuses = self.batches = 0
        self.chord_corrections=self.spectral_factorizations=0
        self.maximum_residual = 0.
        self.processes, self.pipes = [], []
        context = mp.get_context('spawn')
        cpus = budget['worker_affinity_cpus'][:min(len(epsilon), budget['workers'])]
        try:
            for index, cpu in enumerate(cpus):
                parent, child = context.Pipe()
                worker = context.Process(target=_actor, args=(child, cpu, graph,
                    epsilon, initial_u, list(range(index, len(epsilon), len(cpus))), plan))
                worker.start(); child.close()
                self.processes.append(worker); self.pipes.append(parent)
            for pipe in self.pipes:
                response = pipe.recv()
                if response.get('event') != 'READY':
                    raise RuntimeError(str(response))
                event('WORKER_READY', **{k: v for k, v in response.items() if k != 'event'})
        except BaseException:
            self.close(); raise

    def __call__(self, gap):
        for pipe in self.pipes:
            pipe.send(('EVALUATE', gap))
        parts = [pipe.recv() for pipe in self.pipes]
        for part in parts:
            if part.get('event') == 'ERROR':
                raise RuntimeError(part['reason']+'\n'+part['traceback'])
        return _aggregate(self, gap, parts)

    def close(self):
        for process, pipe in zip(self.processes, self.pipes):
            if process.is_alive():
                try:
                    pipe.send(('STOP', None))
                except (BrokenPipeError, EOFError):
                    pass
        for process in self.processes:
            process.join(timeout=3)
            if process.is_alive():
                process.terminate(); process.join(timeout=3)
        for pipe in self.pipes:
            pipe.close()


def execute_hold(graph, initial_gap, epsilon, plan, evaluate, *, output=None,
                 event=lambda *args, **kwargs: None):
    validate_plan(plan)
    expected = 2*np.pi*plan['T_K']/plan['Tc_K']*(np.arange(len(epsilon))+.5)
    if not np.allclose(epsilon, expected, rtol=1e-12, atol=0):
        raise ValueError('Reference frequencies do not match the declared temperature')
    material = {k: plan[k] for k in ('T_K', 'Tc_K', 'tau_ee_Tc_ps', 'tau_ep_Tc_ps')}
    gap = initial_gap.copy(); base = abs(initial_gap)
    baseline_scale = float(np.average(base, weights=graph.area_weights))
    if baseline_scale <= 0:
        raise ValueError('Nonzero superconducting reference required')
    I0 = Boltzmann*plan['Tc_K']/(2*elementary_charge*plan['sheet_resistance_ohm'])
    U0 = hbar*I0/(2*elementary_charge)
    ports = DCCurrentPorts(graph, Tc_K=plan['Tc_K'], sheet_resistance_ohm=plan['sheet_resistance_ohm'])
    first_value = evaluate(gap)
    divergence = np.zeros(graph.n_nodes)
    np.add.at(divergence, graph.edges[:, 0], first_value['current'])
    np.add.at(divergence, graph.edges[:, 1], -first_value['current'])
    initial_current_A = float(.5*I0*(divergence[ports.left].sum()-divergence[ports.right].sum()))
    state = np.array([initial_current_A, initial_current_A, 0.])
    circuit = ThesisCircuitParameters(**plan['circuit'],
        V_bias_V=plan['circuit']['R_bias_ohm']*initial_current_A)
    intervals = int(np.ceil(plan['duration_ps']/plan['dt_ps']))
    dt = plan['duration_ps']/intervals
    history, snapshots = [], []
    initial_energy = None; integrated_loss = 0.; previous_loss = None
    integrated_port_work = 0.; previous_port = None
    started = time.monotonic(); last_progress = -float('inf')
    maximum = dict(gap_drift_relative=0., gap_ripple_relative=0.,
        current_cut_relative=0., current_bias_mismatch_relative=0.,
        force_density_rms=0., interior_potential_absolute_V=0.,
        dark_readout_absolute_V=0.,device_voltage_absolute_V=0., energy_balance_absolute_J=0.,
        circuit_current_drift_relative=0., complete_power_residual_W=0.)
    condensation_scale = float(U0*np.dot(graph.area_weights, base**2))
    for index in range(intervals+1):
        value = first_value if index == 0 else evaluate(gap)
        model = CurrentDrivenThermalKWT(graph, gap, **material,
            ports=ports, imposed_current_A=state[1])
        response = model.response(value['gradient'], value['current'])
        voltage = model.last_voltage
        free = model.free
        cut_x, currents_bar = cut_currents(graph, value['current']+response['normal_current'])
        currents = I0*currents_bar
        scale = max(abs(plan['bias_current_A']), plan['current_absolute_floor_A'])
        loss = U0*(response['kwt_loss']+response['normal_loss'])/model.tD_ps
        port_per_ps = voltage.port_power_W*1e-12
        if initial_energy is None:
            initial_energy = value['energy']
        elif previous_loss is not None:
            integrated_loss += .5*dt*(loss+previous_loss)
            integrated_port_work += .5*dt*(port_per_ps+previous_port)
        previous_loss = loss
        previous_port = port_per_ps
        drift = float(np.max(abs(abs(gap)-base))/baseline_scale)
        if np.max(abs(gap-initial_gap))/baseline_scale > plan['maximum_relative_displacement']:
            raise RuntimeError('DC trajectory left the admitted neighborhood; no clipping or baseline-force subtraction')
        potential_V = voltage.phi_V
        # The bulk spectral phases move with their contact condensate phases.
        # Their conjugate work is required by the common-action Noether identity.
        spectral_boundary_rate = -.5*np.dot(value['spectral_torque'][graph.boundary_nodes],
            response['potential_v'][graph.boundary_nodes])
        complete_rate_W = float(U0*(response['free_energy_rate']+spectral_boundary_rate)/model.tD_ps*1e12)
        total_loss_W = loss*1e12
        circuit_outputs = circuit.outputs(state, voltage.Vdev_V)
        row = dict(step=index, time_ps=index*dt,
            gap_drift_relative=drift,
            gap_ripple_relative=float(np.ptp(abs(gap))/baseline_scale),
            current_cut_relative=float(np.ptp(currents)/scale),
            current_bias_mismatch_relative=float(np.max(abs(currents-state[1]))/scale),
            circuit_current_drift_relative=float(abs(state[1]-initial_current_A)/scale),
            current_cut_mean_A=float(currents.mean()),
            force_density_rms=float(np.sqrt(np.sum(abs(value['gradient'][free])**2/
                graph.area_weights[free])/graph.area_weights[free].sum())),
            interior_potential_absolute_V=float(np.max(abs(potential_V[free]), initial=0.)),
            dark_readout_absolute_V=abs(circuit_outputs.Vout_V),
            device_voltage_absolute_V=abs(voltage.Vdev_V),
            V_out_V=circuit_outputs.Vout_V,
            circuit_I_b_A=float(state[0]), circuit_I_s_A=float(state[1]),
            circuit_v_c_V=float(state[2]), terminal_voltage_V=voltage.Vdev_V,
            free_energy_difference_J=U0*(value['energy']-initial_energy),
            integrated_dissipation_J=integrated_loss,
            integrated_port_work_J=integrated_port_work,
            energy_balance_absolute_J=abs(U0*(value['energy']-initial_energy)+integrated_loss-integrated_port_work),
            device_port_power_W=voltage.port_power_W,
            free_energy_rate_with_boundary_work_W=complete_rate_W,
            spectral_boundary_work_W=float(U0*spectral_boundary_rate/model.tD_ps*1e12),
            complete_power_residual_W=abs(complete_rate_W+total_loss_W-voltage.port_power_W),
            circuit_power_residual_W=circuit.power_balance(state,voltage.Vdev_V).residual_W,
            normal_dissipation_W=U0*response['normal_loss']/model.tD_ps*1e12,
            kwt_dissipation_W=U0*response['kwt_loss']/model.tD_ps*1e12,
            continuity_max=float(response['continuity_max']),
            maximum_spectral_residual=evaluate.maximum_residual)
        for name in maximum:
            maximum[name] = max(maximum[name], row[name])
        observe = index % plan['observation_cadence'] == 0 or index == intervals
        if observe:
            history.append(row)
            snapshots.append(dict(time_ps=index*dt, gap=gap.copy(),
                current_bar=value['current'].copy(), normal_current_bar=response['normal_current'].copy(),
                potential_V=potential_V.copy(), force_density=value['gradient']/graph.area_weights))
            if output is not None:
                atomic_json(Path(output)/'history.json', history)
                np.savez_compressed(Path(output)/'latest_checkpoint.npz',
                    time_ps=index*dt, delta_bar=gap, current_bar=value['current'],
                    potential_V=potential_V, circuit_state_A_A_V=state)
        if observe or time.monotonic()-last_progress >= 5.:
            last_progress=time.monotonic()
            fraction = index/intervals
            eta = None if index == 0 else (time.monotonic()-started)*(intervals-index)/index
            event('DC_HOLD', fraction=fraction, eta_seconds=eta, **row)
            print(f"[{'#'*int(24*fraction):24s}] DC hold {index}/{intervals} | "
                  f"t={index*dt:.4g} ps | ETA {eta if eta is not None else 'estimating'} s", flush=True)
        if index == intervals:
            break
        update = stable_inherited_kwt_step(model, value['gradient'], value['current'],
            dt_ps=dt, delta0_over_kBTc=plan['delta0_over_kBTc'])
        if np.any(~np.isfinite(update.gap)) or not np.array_equal(
                update.gap[graph.boundary_nodes], gap[graph.boundary_nodes]):
            raise RuntimeError('Invalid field or local KWT unexpectedly moved bulk boundary')
        gap = advance_bulk_contact_phases(update.gap, response['potential_v'],
            graph.boundary_nodes, dt_ps=dt, tD_ps=model.tD_ps)
        if not np.allclose(abs(gap[graph.boundary_nodes]), base[graph.boundary_nodes],rtol=1e-11,atol=0):
            raise RuntimeError('Bulk boundary amplitude drifted')
        state = state+dt*1e-12*circuit.rhs(state, voltage.Vdev_V)
        if np.any(~np.isfinite(state)):
            raise RuntimeError('Nonfinite full-circuit state')
    limits = dict(gap_drift_relative=plan['gap_drift_relative_limit'],
        gap_ripple_relative=plan['gap_ripple_relative_limit'],
        current_cut_relative=plan['current_cut_relative_limit'],
        current_bias_mismatch_relative=plan['current_cut_relative_limit'],
        circuit_current_drift_relative=plan['current_cut_relative_limit'])
    checks = {name: bool(maximum[name] <= limit) for name, limit in limits.items()}
    requested_mismatch = abs(initial_current_A-plan['bias_current_A'])/max(
        abs(plan['bias_current_A']), plan['current_absolute_floor_A'])
    checks['reference_vs_requested_current'] = bool(requested_mismatch <= plan['current_cut_relative_limit'])
    dark_scale = circuit.R_load_ohm*max(abs(initial_current_A), plan['current_absolute_floor_A'])
    dark_ratio = maximum['dark_readout_absolute_V']/dark_scale
    energy_ratio = maximum['energy_balance_absolute_J']/condensation_scale
    limits['dark_readout_fraction_of_Rload_Iref'] = plan.get('dark_readout_fraction_limit',1e-3)
    limits['thermal_energy_balance_relative_to_condensation_scale'] = plan.get('energy_balance_relative_limit',1e-2)
    checks['dark_readout_fraction_of_Rload_Iref'] = bool(dark_ratio <= limits['dark_readout_fraction_of_Rload_Iref'])
    checks['thermal_energy_balance_relative_to_condensation_scale'] = bool(energy_ratio <= limits['thermal_energy_balance_relative_to_condensation_scale'])
    rhs = circuit.rhs(state, voltage.Vdev_V)
    summary = dict(status='COMPLETED', dc_hold_passed=bool(all(checks.values())),
        checks=checks, limits=limits, maxima=maximum, accepted_steps=intervals,
        dt_ps=dt, duration_ps=plan['duration_ps'], time_order=1,
        initial=history[0], final=history[-1], circuit=dict(vars(circuit)),
        history_sampling='Initial, final and observation cadence; maxima and energy integrals evaluate every accepted step',
        final_circuit_rhs_A_s_A_s_V_s=rhs.tolist(),
        requested_current_A=plan['bias_current_A'], initial_discrete_current_A=float(initial_current_A),
        reference_vs_requested_current_relative=float(requested_mismatch),
        energy_balance_relative_to_condensation_scale=energy_ratio,
        dark_readout_fraction_of_Rload_Iref=dark_ratio,
        dark_readout_scale_V=float(dark_scale),
        dark_readout_scale_meaning='R_load times initial DC current is a pulse voltage scale, not an experimental photon trigger threshold',
        condensation_scale_J=float(condensation_scale), current_scale_A=float(I0),
        condensation_scale_meaning='U0 times sum(node_area*abs(initial_gap_bar)**2); reference energy scale, not a calculated F_normal minus F_superconducting',
        circuit_inductance_policy='One fixed external series L; no time-dependent kinetic inductance or dynamic subtraction',
        circuit_scope='Full three-state CM Euler step; Is drives equipotential strip ports, Vdev drives the circuit, and contact phases obey Josephson',
        population_scope='Prescribed equilibrium hL=tanh(E/(2kBT)), hT=0 and Bose occupation; no kinetic population evolution',
        potential_scope='Inherited Ohmic normal-current closure on an equipotential-end quotient graph, with imposed circuit current; right port is voltage gauge zero',
        energy_scope='Isothermal finite-Matsubara free energy and inherited KWT dissipation; not nonlinear internal-heat closure',
        boundary_scope='Both longitudinal amplitudes and spectral magnitudes fixed to the biased bulk reference; end phases evolve with port voltages; natural lateral boundary',
        photon=False, ac_source=False, nonlinear_heat_closed=False,
        production_changed=False, source_force_subtracted=False,
        local_euler_algebra='Inherited KWT quadratic solved for delta(|psi|^2), algebraically equivalent and first order; avoids stationary large-gamma cancellation',
        numerical_worker_solves=evaluate.newton_solves, spectral_stationary_reuses=evaluate.reuses,
        spectral_chord_corrections=evaluate.chord_corrections,
        spectral_factorizations=evaluate.spectral_factorizations,
        spectral_cache_policy='Resident sparse Jacobian predicts one chord correction; every accepted spectrum satisfies the unchanged nonlinear residual tolerance; otherwise full Newton refreshes the factor',
        spectral_batches=evaluate.batches, maximum_spectral_residual=evaluate.maximum_residual)
    if output is not None:
        output = Path(output)
        atomic_json(output/'summary.json', summary)
        arrays = {key: np.asarray([snapshot[key] for snapshot in snapshots]) for key in snapshots[0]}
        np.savez_compressed(output/'fields.npz', **arrays, initial_gap=initial_gap,
            coordinates_bar=graph.coordinates_bar, area_weights=graph.area_weights,
            edges=graph.edges, conductance=graph.conductance, boundary_nodes=graph.boundary_nodes,
            cut_positions_bar=cut_x, epsilon_bar=epsilon)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('plan', 'reference', 'output'):
        parser.add_argument('--'+key, type=Path, required=True)
    parser.add_argument('--workers', type=int, default=1)
    args = parser.parse_args()
    plan = validate_plan(json.loads(args.plan.read_text(encoding='utf8')))
    if args.output.exists():
        raise FileExistsError('Use a fresh output directory; no overwrite')
    graph, gap, epsilon, initial_u = load_reference(args.reference)
    if plan.get('reference_sha256') is not None and sha(args.reference) != plan['reference_sha256']:
        raise ValueError('Reference hash mismatch')
    for path, digest in plan.get('source_sha256', {}).items():
        digest_function=source_sha if plan.get('source_hash_normalization')=='CRLF-to-LF' else sha
        if digest_function(ROOT/path) != digest:
            raise ValueError('Pinned source mismatch: '+path)
    resources = linux_resources()
    inherited = os.environ.get('PYSNSPD_SHARED_ALLOCATION')
    budget = admit_runtime(args.workers, resources, json.loads(inherited) if inherited else None)
    args.output.mkdir(parents=True)
    manifest = dict(status='RUNNING', plan=plan, resources=resources, budget=budget,
        plan_sha256=sha(args.plan), reference_sha256=sha(args.reference),
        reference_path=str(args.reference.resolve()), source_sha256=sha(__file__))
    atomic_json(args.output/'manifest.json', manifest)
    started = time.monotonic(); previous_affinity = os.sched_getaffinity(0)
    log = (args.output/'progress.jsonl').open('x', encoding='utf8', buffering=1)
    def event(name, **values):
        record = dict(event=name, elapsed_seconds=time.monotonic()-started, **values)
        line = json.dumps(record, allow_nan=False)
        log.write(line+'\n'); print(line, flush=True)
    pool = None
    try:
        os.sched_setaffinity(0, {budget['coordinator_cpu']})
        event('START', workers=budget['workers'], nodes=graph.n_nodes, frequencies=len(epsilon))
        pool = ResidentSpectralPool(graph, epsilon, initial_u, plan, budget, event)
        summary = execute_hold(graph, gap, epsilon, plan, pool, output=args.output, event=event)
        manifest.update(status='COMPLETED', elapsed_seconds=time.monotonic()-started,
            dc_hold_passed=summary['dc_hold_passed'])
        atomic_json(args.output/'manifest.json', manifest)
        event('COMPLETE', dc_hold_passed=summary['dc_hold_passed'])
    except BaseException as exc:
        manifest.update(status='FAILED', reason=str(exc), exception=type(exc).__name__,
            elapsed_seconds=time.monotonic()-started)
        atomic_json(args.output/'manifest.json', manifest)
        event('FAILED', reason=str(exc)); raise
    finally:
        if pool is not None:
            pool.close()
        log.close(); os.sched_setaffinity(0, previous_affinity)


if __name__ == '__main__':
    main()
