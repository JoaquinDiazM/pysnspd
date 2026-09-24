"""Lightweight unit-conversion check against the actual inherited KWT update."""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import time
from unittest.mock import patch

import numpy as np
import scipy

from pysnspd.experimental import thermal_spatial_usadel as thermal
from pysnspd.experimental import thermal_stable_newton as stable
from pysnspd.experimental.heredado_kwt_bridge import inherited_kwt_step
from pysnspd.experimental.thermal_weak_response import ThermalKWTNormal
from pysnspd.solver.core import TDGLSolver


def bridge_diagnostic():
    start = time.monotonic()
    graph = thermal.rectangular_graph(9, 7, 6., 4.)
    x, y = graph.coordinates_bar.T
    envelope = np.sin(np.pi*x/6) * np.cos(np.pi*y/4)
    gap = (1.55 + .08*envelope) * np.exp(.12j*envelope)
    temperature, critical = .9, 8.65
    ratio = temperature/critical
    spectral = []
    for index in range(16):
        epsilon = np.pi*ratio*(2*index+1)
        spectral.append(stable.solve_frequency(graph, gap, epsilon,
            fixed_nodes=graph.boundary_nodes,
            fixed_u=gap[graph.boundary_nodes]/epsilon, tol=1e-7))
    evaluation = thermal.evaluate_thermal(graph, gap, ratio, spectral)
    model = ThermalKWTNormal(graph, gap, Tc_K=critical, T_K=temperature,
                            tau_ee_Tc_ps=.5, tau_ep_Tc_ps=2.47)
    response = model.response(evaluation.gap_gradient, evaluation.current_bar)
    exact = response['velocity'] / model.tD_ps
    phase_axis = gap/np.abs(gap)
    def norm(value):
        return float(np.sqrt(np.dot(graph.area_weights, np.abs(value)**2)))
    material_velocity = response['material_velocity']/model.tD_ps
    gauge_velocity = -.5j*response['potential_v']*gap/model.tD_ps
    records = []
    original = TDGLSolver.solve_for_psi_squared
    with patch.object(TDGLSolver, 'solve_for_psi_squared', wraps=original) as call:
        for dt_ps in (1e-4, 5e-5, 2.5e-5):
            result = inherited_kwt_step(model, evaluation.gap_gradient,
                evaluation.current_bar, dt_ps=dt_ps, delta0_over_kBTc=1.764)
            velocity = (result.gap-gap)/dt_ps
            error = velocity-exact
            records.append(dict(dt_ps=dt_ps,
                velocity_relative_error=norm(error)/norm(exact),
                amplitude_velocity_absolute_error=norm(np.real(np.conj(phase_axis)*error)),
                phase_velocity_absolute_error=norm(np.imag(np.conj(phase_axis)*error)),
                fixed_contact_change=float(np.max(abs(result.gap[~model.free]-gap[~model.free]))),
                backend_dt=result.backend_dt,
                backend_tau0_ps=result.backend_tau0_ps,
                backend_gamma=result.backend_gamma))
        direct_calls = call.call_count
    errors = [item['velocity_relative_error'] for item in records]
    ratios = [errors[i]/errors[i+1] for i in range(len(errors)-1)]
    nonzero_fields = dict(
        minimum_gap=float(np.min(abs(gap))),
        integrated_force_norm=norm(evaluation.gap_gradient),
        potential_max=float(np.max(abs(response['potential_v']))),
        material_velocity_norm_per_ps=norm(material_velocity),
        gauge_velocity_norm_per_ps=norm(gauge_velocity),
        total_velocity_norm_per_ps=norm(exact))
    passed = bool(direct_calls == 3 and
        all(1.8 < value < 2.2 for value in ratios) and
        errors[-1] < 5e-4 and
        all(item['fixed_contact_change'] == 0 for item in records) and
        nonzero_fields['minimum_gap'] > 1 and
        nonzero_fields['potential_max'] > 1e-3 and
        nonzero_fields['material_velocity_norm_per_ps'] > 1e-3 and
        nonzero_fields['gauge_velocity_norm_per_ps'] > 1e-3)
    return dict(status='PASSED' if passed else 'FAILED',
        runtime_seconds=time.monotonic()-start,
        nodes=graph.n_nodes, free_nodes=int(np.count_nonzero(model.free)),
        temperature_K=temperature, critical_temperature_K=critical,
        gap_scale_Delta0_over_kBTc=1.764,
        spectral_frequencies=16, spectral_tolerance=1e-7,
        maximum_spectral_residual=max(item.residual for item in spectral),
        direct_inherited_solver_calls=direct_calls,
        reference='ThermalKWTNormal.response velocity divided by tD_ps',
        error_norm='sqrt(sum_i m_i |delta velocity_i|^2), integrated adimensional area',
        nonzero_fields=nonzero_fields, refinements=records,
        first_order_ratios=ratios,
        scope='Local KWT time and force conversion on a smooth nonzero-gap state. '
              'No finite-horizon trajectory, GL current, circuit or photon validation.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit('Refusing to overwrite a diagnostic receipt')
    result = bridge_diagnostic()
    root = Path(__file__).resolve().parents[2]
    paths = [
        'pysnspd/experimental/heredado_kwt_bridge.py',
        'pysnspd/experimental/thermal_weak_response.py',
        'pysnspd/experimental/thermal_spatial_usadel.py',
        'pysnspd/experimental/thermal_stable_newton.py',
        'pysnspd/solver/core.py',
        'sandbox/stage4_core/kwt_bridge_check.py',
        'tests/test_heredado_kwt_bridge.py',
    ]
    result['source_sha256'] = {path: hashlib.sha256((root/path).read_bytes()).hexdigest()
                               for path in paths}
    result['runtime'] = dict(python=platform.python_version(), numpy=np.__version__,
                             scipy=scipy.__version__)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, indent=2))
    if result['status'] != 'PASSED':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
