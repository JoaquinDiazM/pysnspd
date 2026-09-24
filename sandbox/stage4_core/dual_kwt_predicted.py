"""The frozen dual KWT trajectory with a verified uniform spectral predictor.

Only the worker backend changes. A previously implemented SpectralTangent
factorization predicts each new spectral state. The exact nonlinear residual
must satisfy the SAME tolerance; otherwise the original Newton solver runs.
The inherited Euler integrator, forces, currents, energies and boundary
conditions are imported unchanged from dual_kwt_time.
"""
from __future__ import annotations
import multiprocessing as mp
import os
from pathlib import Path
import sys
import traceback
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from sandbox.stage4_core import dual_kwt_time as base
import numpy as np
from pysnspd.experimental.thermal_weak_response import SpectralTangent


class UniformPredictorMode:
    """One exact uniform Jacobian shared by all trajectories of one frequency."""
    def __init__(self, graph, d0, epsilon, plan):
        self.graph, self.d0, self.epsilon, self.plan = graph, d0, epsilon, plan
        self.free = np.ones(graph.n_nodes, bool)
        self.free[graph.boundary_nodes] = False
        self.u0 = d0/epsilon
        value = base.thermal.spectral_energy_gradient(graph, d0, epsilon, self.u0)
        residual = self.residual(value, d0)
        if residual > plan['spectral_tolerance']:
            raise RuntimeError('Analytic uniform reference is not spectrally stationary')
        scale = graph.area_weights*np.maximum(1., abs(d0))
        gradient_residual = float(np.max(abs(value.gradient_u[self.free])/scale[self.free]))
        solution = base.thermal.SpectralSolution(float(epsilon), self.u0, value.f,
            value.g, value.energy, residual, gradient_residual, 0, tuple(),
            graph.boundary_nodes, base.thermal._signature(graph, d0, np.zeros(len(graph.edges))))
        self.tangent = SpectralTangent(graph, d0, solution)
        self.cache = {}

    def residual(self, evaluated, gap):
        scale = self.graph.area_weights*np.maximum(1., abs(gap))
        return float(np.max(abs(evaluated.residual[self.free])/scale[self.free]))

    def evaluate(self, name, gap):
        previous_gap, previous_u = self.cache.get(name, (self.d0, self.u0))
        increment = self.tangent.apply(gap-previous_gap).du
        predictor = previous_u+increment
        predictor[self.graph.boundary_nodes] = self.u0[self.graph.boundary_nodes]
        value = base.thermal.spectral_energy_gradient(self.graph, gap, self.epsilon, predictor)
        predicted_residual = self.residual(value, gap)
        if predicted_residual <= self.plan['spectral_tolerance']:
            u, residual, iterations, solved = predictor, predicted_residual, 0, False
        else:
            solution = base.solve_frequency(self.graph, gap, self.epsilon,
                initial_u=predictor, fixed_nodes=self.graph.boundary_nodes,
                fixed_u=self.u0[self.graph.boundary_nodes],
                tol=self.plan['spectral_tolerance'],
                max_iterations=self.plan['maximum_newton_iterations'])
            u, residual, iterations, solved = (solution.u, solution.residual,
                solution.iterations, True)
        self.cache[name] = (gap.copy(), u.copy())
        return dict(u=u, residual=residual, predicted_residual=predicted_residual,
            iterations=iterations, newton_solved=solved)


def predicted_actor(pipe, cpu, graph, frequencies, names, d0, epsilon, plan):
    try:
        os.sched_setaffinity(0, {cpu})
        modes = {n: UniformPredictorMode(graph, d0, epsilon[n], plan) for n in frequencies}
        weight = 2*np.pi*plan['T_K']/plan['Tc_K']
        pipe.send(dict(event='READY', cpu=cpu, uniform_factorizations=len(modes),
            task_count=len(modes)*len(names), predictor='incremental exact uniform SpectralTangent; nonlinear residual verified'))
        while True:
            command, fields = pipe.recv()
            if command == 'STOP':
                break
            if command != 'EVALUATE':
                raise ValueError('Unknown predicted spectral request')
            output = {name: [np.zeros(graph.n_nodes, complex),
                np.zeros(len(graph.edges)), 0.] for name in fields}
            maximum = 0.
            solved = predicted = iterations = 0
            for n, mode in modes.items():
                for name, gap in fields.items():
                    result = mode.evaluate(name, gap)
                    solved += int(result['newton_solved'])
                    predicted += int(not result['newton_solved'])
                    iterations += result['iterations']
                    diff = base.spectral_difference(graph, d0, mode.u0,
                        gap, result['u'], epsilon[n])
                    output[name][0] += weight*diff['gap_gradient_difference']
                    output[name][1] += weight*diff['current_difference']
                    output[name][2] += weight*diff['renormalized_energy_difference']
                    maximum = max(maximum, result['residual'])
            pipe.send(dict(output=output, maximum_residual=maximum,
                newton_solves=solved, stationary_reuses=predicted, newton_iterations=iterations))
    except BaseException as error:
        try:
            pipe.send(dict(event='ERROR', reason=str(error), traceback=traceback.format_exc()))
        except (BrokenPipeError, EOFError):
            pass
    finally:
        pipe.close()


class UniformPredictorPool(base.SpectralPool):
    """Reuse the frozen pool's aggregation and closure; replace worker setup."""
    def __init__(self, graph, names, d0, epsilon, plan, budget, event):
        self.graph, self.d0 = graph, d0
        self.logt = np.log(plan['T_K']/plan['Tc_K'])
        self.processes, self.pipes = [], []
        self.batches = self.newton_solves = self.stationary_reuses = self.iterations = 0
        self.maximum_residual = 0.
        context = mp.get_context('spawn')
        frequencies = list(range(len(epsilon)))
        try:
            for index, cpu in enumerate(budget['worker_affinity_cpus']):
                parent, child = context.Pipe()
                process = context.Process(target=predicted_actor, args=(child, cpu,
                    graph, frequencies[index::budget['workers']], names, d0, epsilon, plan))
                process.start(); child.close()
                self.processes.append(process); self.pipes.append(parent)
            for index, pipe in enumerate(self.pipes):
                result = pipe.recv()
                if result.get('event') != 'READY':
                    raise RuntimeError(str(result))
                event('WORKER_READY', worker=index, **{k:v for k,v in result.items() if k != 'event'})
            event('PREDICTOR_READY', uniform_factorizations=len(epsilon),
                stationary_reuses_meaning='Spectral predictions accepted only after the exact nonlinear residual meets the registered tolerance',
                fallback='Original stable Newton with the predicted state as initial guess')
        except BaseException:
            self.close()
            raise


def main():
    # Explicit backend selection for this separate entry point. All integration
    # and output code remains the byte-identical previously executed runner.
    base.SpectralPool = UniformPredictorPool
    base.main()


if __name__ == '__main__':
    main()
