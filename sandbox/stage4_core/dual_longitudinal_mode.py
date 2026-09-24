"""Exact invariant-mode control of weak longitudinal gap/population exchange.

At a uniform, real equilibrium every radial Matsubara response is a function
of the same mass-weighted graph Laplacian. A single spatial eigenmode is an
invariant subspace of the linearized equations, not a detector dimensional
reduction. Euler advances its (1 + number of energies) amplitudes; a small
matrix exponential is only an independent temporal reference. No collisions,
photon, charge imbalance, physical material rate or total internal energy is
supplied by this control. The main nonlinear trajectory uses inherited KWT.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if __name__ == '__main__':
    from sandbox.stage4_core.parallel_runtime import limit_thread_environment
    limit_thread_environment()

import numpy as np
from scipy.linalg import expm
from scipy.optimize import brentq
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import eigsh

from pysnspd.experimental import frozen_kinetic_usadel as kinetic
from pysnspd.experimental import retarded_spatial_usadel as retarded
from pysnspd.experimental import thermal_spatial_usadel as thermal
from pysnspd.experimental.longitudinal_reciprocal import (
    LongitudinalReciprocalBridge, thermal_susceptibility)
from pysnspd.experimental.thermal_weak_response import ThermalKWTNormal


def uniform_equilibrium(temperature_ratio, matsubara_count):
    eps = 2*np.pi*temperature_ratio*(np.arange(matsubara_count)+.5)
    gap = brentq(lambda d: np.log(temperature_ratio)+2*np.pi*temperature_ratio*
        np.sum(1/eps-1/np.sqrt(eps*eps+d*d)), .01, 4.)
    return gap, eps


def lowest_mode(graph):
    """Full dual graph, fixed end contacts, natural side boundaries."""
    free = np.ones(graph.n_nodes, bool)
    free[graph.boundary_nodes] = False
    tail, head = graph.edges.T
    c = graph.conductance
    lap = coo_matrix((np.r_[c, c, -c, -c],
        (np.r_[tail, head, tail, head], np.r_[tail, head, head, tail])),
        shape=(graph.n_nodes, graph.n_nodes)).tocsr()
    mass = diags(graph.area_weights[free])
    values, vectors = eigsh(lap[free][:, free], k=1, M=mass, sigma=0.,
        which='LM', v0=np.ones(free.sum()), tol=1e-10)
    vector = np.zeros(graph.n_nodes)
    vector[free] = vectors[:, 0]
    if np.sum(vector) < 0:
        vector *= -1
    vector /= np.max(abs(vector))
    residual = lap@vector-values[0]*graph.area_weights*vector
    scale = np.linalg.norm((values[0]*graph.area_weights*vector)[free])
    return float(values[0]), vector, float(np.linalg.norm(residual[free])/scale)


def radial_hessian_symbol(eigenvalue, gap, eps, temperature_ratio):
    """H v = M v h(lambda), from J_n = eps_n M + g_n L.

    At the real uniform state df_n = g_n**3 du_n. Thus a Laplacian
    eigenvector satisfies df_n/v = g_n**3/(eps_n + g_n lambda).
    This identity is checked against the existing exact SpectralTangent.
    """
    g = eps/np.sqrt(eps*eps+gap*gap)
    return float(2*np.log(temperature_ratio)+4*np.pi*temperature_ratio*
        np.sum(1/eps-g**3/(eps+g*eigenvalue)))


def make_control(graph, *, matsubara_count=256, temperature_ratio=.9/8.65):
    gap, eps = uniform_equilibrium(temperature_ratio, matsubara_count)
    eigenvalue, mode, eigen_residual = lowest_mode(graph)
    hessian = radial_hessian_symbol(eigenvalue, gap, eps, temperature_ratio)
    if hessian <= 0:
        raise ValueError('Control reference is not stable in its chosen mode')
    # Finite, explicitly prescribed diagnostic quadrature, not a converged
    # physical continuum or an electron/phonon collision quadrature.
    energy = np.array([.2, .45, .7, 1., 1.4, 1.8, 2.4])
    weight = np.array([.2, .25, .25, .3, .4, .4, .6])
    eta = .01
    d = np.full(graph.n_nodes, gap, complex)
    operators = []
    for value in energy:
        a, b = retarded.uniform_fields(d, eta-1j*value)
        f, ft, g = retarded.fields(a, b)
        operators.append(kinetic.assemble(graph, d, g, f, ft))
    chi = thermal_susceptibility(energy, temperature_ratio)
    rho = np.array([op.R[0, 0, 0].real for op in operators])
    kernel = np.array([-2*op.R[0, 0, 1].imag for op in operators])
    active_edge = int(np.flatnonzero(graph.conductance > 0)[0])
    diffusion = np.array([op.edge_head_blocks[active_edge, 0, 0]/
        graph.conductance[active_edge] for op in operators])
    model = ThermalKWTNormal(graph, d, Tc_K=8.65, T_K=temperature_ratio*8.65)
    free = model.free
    radial_inverse = 1/((model.denominator/graph.area_weights)*model.stretch)
    mobility = float(radial_inverse[free][0])
    size = len(energy)+1
    matrix = np.zeros((size, size))
    matrix[0, 0] = -mobility*hessian
    matrix[0, 1:] = -mobility*weight*chi*kernel
    matrix[1:] = .5*(kernel/rho)[:, None]*matrix[0]
    matrix[1:, 1:] -= np.diag(diffusion*eigenvalue/rho)
    # Coefficients above use tau = t/tD. Every output time is physical ps.
    matrix /= model.tD_ps
    initial = np.r_[-.002*gap, -.002*np.exp(-((energy-.7)/.3)**2)]
    norm_mass = float(np.dot(graph.area_weights, mode*mode))
    metric = norm_mass*np.r_[.5*hessian, weight*chi*rho]

    def restricted_hessian(argument):
        coefficient = np.dot(graph.area_weights*mode, argument)/norm_mass
        if np.linalg.norm(argument-coefficient*mode) > 1e-10*max(np.linalg.norm(argument), 1e-15):
            raise ValueError('Analytic callback is valid only in the declared invariant mode')
        return hessian*graph.area_weights*argument

    bridge = LongitudinalReciprocalBridge(graph, operators, energy, weight,
        temperature_ratio, restricted_hessian, lambda f:model.local_inverse(f).real,
        spectral_z=eta-1j*energy)
    x = initial[0]*mode
    y = initial[1:, None]*mode
    response = bridge.response(x, y)
    full_velocity = np.vstack((response.gap_velocity,
        response.population_velocity))/model.tD_ps
    embedded_velocity = (matrix@initial)[:, None]*mode
    embedding_error = float(np.linalg.norm(full_velocity-embedded_velocity)/
        np.linalg.norm(embedded_velocity))
    return dict(graph=graph, gap=gap, eps=eps, eigenvalue=eigenvalue,
        mode=mode, eigen_residual=eigen_residual, hessian=hessian,
        energy=energy, weight=weight, chi=chi, rho=rho, kernel=kernel,
        diffusion=diffusion, mobility=mobility, tD_ps=model.tD_ps,
        matrix=matrix, initial=initial, metric=metric, norm_mass=norm_mass,
        embedding_error=embedding_error, reciprocal_power_residual=
        response.reciprocal_power_residual, eta=eta,
        temperature_ratio=temperature_ratio)


def losses(control, state):
    x, y = state[0], state[1:]
    force = control['hessian']*x+np.dot(control['weight']*control['chi']*
        control['kernel'], y)
    kwt = control['norm_mass']*control['mobility']*force*force/control['tD_ps']
    transport = 2*control['norm_mass']*control['eigenvalue']*np.dot(
        control['weight']*control['chi']*control['diffusion'], y*y)/control['tD_ps']
    return np.array([kwt, transport])


def evolve(control, *, steps, duration_ps):
    if steps < 1 or duration_ps <= 0:
        raise ValueError('Positive duration and step count required')
    matrix = control['matrix']
    dt = duration_ps/steps
    if np.max(abs(np.linalg.eigvals(np.eye(len(matrix))+dt*matrix))) > 1+1e-12:
        raise ValueError('Euler step exceeds stability interval; choose a smaller step explicitly')
    states = np.empty((steps+1, len(matrix)))
    states[0] = control['initial']
    heat = np.zeros((steps+1, 2))
    for index in range(steps):
        states[index+1] = states[index]+dt*(matrix@states[index])
        # Trapezoidal observation of the two continuous loss functions. This
        # does not turn the Euler state update into a second-order integrator.
        heat[index+1] = heat[index]+.5*dt*(losses(control, states[index])+
            losses(control, states[index+1]))
    times = np.linspace(0., duration_ps, steps+1)
    available = (states*states)@control['metric']
    exact = expm(duration_ps*matrix)@states[0]
    metric_norm = lambda z:float(np.sqrt(np.dot(control['metric'], z*z)))
    error = metric_norm(states[-1]-exact)/metric_norm(states[0])
    initial_h = np.tanh(control['energy']/(2*control['temperature_ratio']))
    # An eigenvector of this connected positive graph has 0 <= mode <= 1;
    # compute extrema using both extrema rather than assume the sign.
    values = (initial_h[None, :, None]+control['chi'][None, :, None]*
        states[:, 1:, None]*np.array([control['mode'].min(), control['mode'].max()]))
    occupation = .5*(1-values)
    return dict(times=times, states=states, integrated_loss=heat,
        availability=available, exact_final=exact,
        metric_error_over_initial=error,
        integrated_balance_over_initial=float(abs(available[-1]-available[0]+
            heat[-1].sum())/available[0]),
        minimum_occupation=float(occupation.min()), maximum_occupation=float(occupation.max()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh', type=Path, required=True)
    parser.add_argument('--output-root', type=Path, required=True)
    parser.add_argument('--steps', type=int, default=512)
    parser.add_argument('--duration-ps', type=float, default=1.)
    args = parser.parse_args()
    if args.output_root.exists():
        raise SystemExit('Preserving existing control output')
    from sandbox.stage4_core.parallel_runtime import (linux_resources,
        resource_budget, limit_thread_environment)
    limit_thread_environment()
    resources = linux_resources()
    budget = resource_budget(resources, max_workers=1)
    original = os.sched_getaffinity(0)
    os.sched_setaffinity(0, {budget['worker_affinity_cpus'][0]})
    started = time.monotonic()
    try:
        print(json.dumps(dict(event='START', workers=1, processes=1,
            threads_per_process=1, scope='Uniform invariant-mode longitudinal control')), flush=True)
        with np.load(args.mesh) as data:
            graph = thermal.ThermalGraph(data['area_weights'], data['edges'],
                data['conductance'], data['coordinates_bar'], data['fixed_nodes'])
        control = make_control(graph)
        print(json.dumps(dict(event='PROGRESS', fraction=1/3,
            elapsed_seconds=time.monotonic()-started,
            eta_seconds=2*(time.monotonic()-started),
            meaning='Exact graph mode and analytic uniform spectra ready')), flush=True)
        primary = evolve(control, steps=args.steps, duration_ps=args.duration_ps)
        refined = evolve(control, steps=2*args.steps, duration_ps=args.duration_ps)
        passed = bool(control['embedding_error'] < 1e-7 and control['eigen_residual'] < 1e-7
            and primary['metric_error_over_initial'] < .02
            and primary['integrated_balance_over_initial'] < .02
            and primary['minimum_occupation'] >= 0 and primary['maximum_occupation'] <= 1
            and refined['metric_error_over_initial'] < primary['metric_error_over_initial'])
        args.output_root.mkdir(parents=True)
        np.savez_compressed(args.output_root/'evolution.npz',
            coordinates_bar=graph.coordinates_bar, mode=control['mode'],
            times_ps=primary['times'], amplitudes=primary['states'],
            availability=primary['availability'], integrated_losses=primary['integrated_loss'],
            refined_times_ps=refined['times'], refined_amplitudes=refined['states'],
            exact_final=primary['exact_final'], energy_over_kBTc=control['energy'])
        def summary(result):
            return {key:result[key] for key in ('metric_error_over_initial',
                'integrated_balance_over_initial', 'minimum_occupation', 'maximum_occupation')}
        root = Path(__file__).resolve().parents[2]
        paths = ['sandbox/stage4_core/dual_longitudinal_mode.py',
            'pysnspd/experimental/longitudinal_reciprocal.py',
            'pysnspd/experimental/frozen_kinetic_usadel.py',
            'pysnspd/experimental/thermal_weak_response.py',
            'tests/test_dual_longitudinal_mode.py']
        receipt = dict(schema='pysnspd.stage4.dual_longitudinal_invariant_mode.v1',
            status='PASSED_CONTROL' if passed else 'FAILED_CONTROL',
            scope='Exact invariant sector of weak real-gap longitudinal equations at uniform equilibrium',
            excluded=['Nonlinear nonthermal detector closure', 'Internal energy total',
                'Physical material tau_kin', 'Phonon collisions or heat deposition',
                'Charge/phase/circuit trajectory', 'Photon or hotbelt'],
            nodes=graph.n_nodes, free_nodes=graph.n_nodes-len(graph.boundary_nodes),
            T_K=.9, Tc_K=8.65, matsubara_count=256, gap_over_kBTc=control['gap'],
            energies_over_kBTc=control['energy'].tolist(), quadrature_weights=control['weight'].tolist(),
            quadrature_scope='Prescribed finite diagnostic quadrature, not continuum convergence',
            eta=control['eta'], eta_is_physical_bath=False,
            eigenvalue=control['eigenvalue'], eigenmode_residual=control['eigen_residual'],
            full_operator_embedding_relative_error=control['embedding_error'],
            instantaneous_reciprocal_power_residual=control['reciprocal_power_residual'],
            initial_amplitudes=control['initial'].tolist(),
            preparation='0.2 percent gap suppression plus independent nonthermal energy profile; no photon',
            population_definition='delta hL(E,X)=chi(E)*y(E)*mode(X); occupation=(1-hL)/2',
            gap_definition='delta Delta/(kBTc)=x*mode(X)',
            availability_definition='Quadratic free-energy availability of perturbation, not total internal energy',
            loss_columns=['KWT availability loss', 'Longitudinal transport availability loss'],
            steps=args.steps, refined_steps=2*args.steps, duration_ps=args.duration_ps,
            numerical_method='Forward Euler of the same radial KWT tangent; matrix exponential only reference',
            primary=summary(primary), refined=summary(refined),
            error_reduction_ratio=primary['metric_error_over_initial']/refined['metric_error_over_initial'],
            final_loss_over_initial=(primary['integrated_loss'][-1]/primary['availability'][0]).tolist(),
            final_availability_over_initial=float(primary['availability'][-1]/primary['availability'][0]),
            tolerance_policy='2 percent of initial availability metric for temporal observables and integrated balance; support and sign retained',
            runtime_seconds=time.monotonic()-started,
            actual_resources=dict(workers=1, processes=1, blas_threads=1,
                pinned_cpus=sorted(os.sched_getaffinity(0)), detected=resources,
                admission_budget=budget),
            mesh_sha256=hashlib.sha256(args.mesh.read_bytes()).hexdigest(),
            source_sha256={path:hashlib.sha256((root/path).read_bytes()).hexdigest() for path in paths})
        (args.output_root/'receipt.json').write_text(json.dumps(receipt, indent=2)+'\n', encoding='utf-8')
        print(json.dumps(dict(event='COMPLETE', status=receipt['status'], fraction=1.,
            eta_seconds=0., runtime_seconds=receipt['runtime_seconds'], primary=receipt['primary'])), flush=True)
        if not passed:
            raise SystemExit(1)
    finally:
        os.sched_setaffinity(0, original)


if __name__ == '__main__':
    main()
