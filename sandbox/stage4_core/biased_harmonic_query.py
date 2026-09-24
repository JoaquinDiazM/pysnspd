"""One independently parallelizable energy query on a biased Usadel graph.

No integration over energy is performed here. The worker preserves signed
energy, complex Fourier amplitudes, explicit reservoir response and the full
spectral/kinetic graph. It is intended for the user-run long campaign.
"""
from __future__ import annotations
import argparse
from functools import lru_cache
import hashlib
import json
import os
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from sandbox.stage4_core.parallel_runtime import limit_thread_environment
limit_thread_environment()
import numpy as np
from scipy.sparse.linalg import splu
from pysnspd.experimental import retarded_spatial_usadel as spectral
from pysnspd.experimental.thermal_spatial_usadel import ThermalGraph
from pysnspd.experimental.harmonic_usadel import HarmonicUsadelFactor
from pysnspd.experimental import harmonic_kinetic_usadel as kinetic

I, Z = kinetic.I2, kinetic.TAU3


def _first(arr, n):
    arr = np.asarray(arr)
    if arr.shape == (n,):
        return arr
    if arr.ndim == 2 and arr.shape[1] == n:
        return arr[0]
    raise ValueError('Expected node fields or Matsubara-by-node fields')


@lru_cache(maxsize=4)
def _load_reference(path_string, mtime_ns, byte_count):
    path = Path(path_string)
    with np.load(path, allow_pickle=False) as archive:
        data = {name: archive[name].copy() for name in archive.files}
    d = np.asarray(data.get('delta_bar', data.get('d')), complex)
    n = len(d)
    graph = ThermalGraph(data['area_weights'], data['edges'], data['conductance'],
        data['coordinates_bar'], data['boundary_nodes'])
    if d.shape != (graph.n_nodes,) or np.any(~np.isfinite(d)):
        raise ValueError('Reference must contain a finite complex gap per node')
    epsilon = np.asarray(data.get('epsilon_bar', data.get('epsilon')), float).reshape(-1)
    if not len(epsilon) or epsilon[0] <= 0:
        raise ValueError('Positive lowest Matsubara frequency required')
    if 'f' in data and 'g' in data:
        f, g = _first(data['f'], n), _first(data['g'], n)
        a, b = f/(1+g), np.conj(f)/(1+g)
    elif 'u' in data:
        # The thermal module uses u=f/g, not retarded stereographic a=f/(1+g).
        u = np.asarray(_first(data['u'], n), complex)
        g = 1/np.sqrt(1+abs(u)**2)
        f = u*g
        a, b = f/(1+g), np.conj(f)/(1+g)
    else:
        raise ValueError('Reference requires lowest Matsubara f/g or u')
    alpha = np.asarray(data.get('alpha', np.zeros(len(graph.edges))), float)
    return graph, d, float(epsilon[0]), a, b, alpha, data, hashlib.sha256(path.read_bytes()).hexdigest()


def _reference(payload):
    path = Path(payload['reference_npz']).resolve()
    stat = path.stat()
    result = _load_reference(str(path), stat.st_mtime_ns, stat.st_size)
    if payload.get('reference_sha256') and payload['reference_sha256'] != result[-1]:
        raise ValueError('Reference hash mismatch')
    return result


def _response_basis(payload, archive, graph):
    source = archive
    if payload.get('modes_npz'):
        with np.load(payload['modes_npz'], allow_pickle=False) as loaded:
            source = {**archive, **{name: loaded[name].copy() for name in loaded.files}}
    modes = np.asarray(payload.get('modes', source.get('modes', source.get('vectors'))))
    if modes.ndim != 2 or modes.shape[0] != graph.n_nodes or np.iscomplexobj(modes) or np.any(~np.isfinite(modes)):
        raise ValueError('Finite real node-by-mode basis required')
    if 'mode_indices' in payload:
        modes = modes[:, np.asarray(payload['mode_indices'], int)]
    lift = np.asarray(payload.get('lift', source.get('lift')), float)
    if lift.shape != (graph.n_nodes,) or np.any(~np.isfinite(lift)):
        raise ValueError('Finite harmonic voltage lift required')
    if np.max(abs(modes[graph.boundary_nodes]), initial=0) > 1e-10:
        raise ValueError('Internal response modes must vanish at the reservoirs')
    gram = modes.T@(graph.area_weights[:, None]*modes)
    error = float(np.max(abs(gram-np.eye(modes.shape[1])), initial=0))
    if error > 1e-7:
        raise ValueError(f'Response modes are not M-orthonormal: {error}')
    return modes, lift, error


def _matrix_fields(solution):
    n = len(solution.g)
    R = np.empty((n, 2, 2), complex)
    R[:, 0, 0], R[:, 1, 1] = solution.g, -solution.g
    R[:, 0, 1], R[:, 1, 0] = solution.f, solution.f_tilde
    A = -Z@R.conj().transpose(0, 2, 1)@Z
    return R, A


def _generators(d, energy, eta):
    B = kinetic.gap_vertex(d.real, d.imag)-1j*energy*Z
    return B+eta*Z, B-eta*Z, B


def _continued_spectrum(graph, d, epsilon, a_anchor, b_anchor, alpha, energy,
                        eta, tolerance, steps):
    """One declared straight path; no retry, alternate contour or clipping."""
    target = complex(eta, -energy)
    a, b = a_anchor.copy(), b_anchor.copy()
    history = []
    fixed = graph.boundary_nodes
    for index, z in enumerate(np.linspace(complex(epsilon), target, steps+1)[1:], 1):
        fa, fb = spectral.uniform_fields(d[fixed], z)
        solution = spectral.solve(graph, d, z, alpha, initial_a=a, initial_b=b,
            fixed_nodes=fixed, fixed_a=fa, fixed_b=fb, tolerance=tolerance,
            causal_tolerance=max(1e-9, tolerance))
        if solution.normalization_residual > max(1e-8, 10*tolerance):
            raise RuntimeError('Retarded continuation lost normalization')
        history.append({'step': index, 'z_real': float(z.real), 'z_imag': float(z.imag),
            'iterations': solution.iterations, 'residual': solution.residual,
            'normalization_residual': solution.normalization_residual,
            'minimum_DOS': solution.minimum_DOS})
        a, b = solution.a, solution.b
    R, A = _matrix_fields(solution)
    return R, A, history


def _column_fields(d, modes, lift, omega):
    phase = np.angle(d)
    c, s = np.cos(phase), np.sin(phase)
    count = modes.shape[1]
    columns = 3*count+1
    x, y, voltage = (np.zeros((len(d), columns), complex) for _ in range(3))
    x[:, :count], y[:, :count] = c[:, None]*modes, s[:, None]*modes
    x[:, count:2*count], y[:, count:2*count] = -s[:, None]*modes, c[:, None]*modes
    voltage[:, 2*count:3*count] = modes
    theta = -2j*lift/omega
    x[:, -1], y[:, -1], voltage[:, -1] = -d.imag*theta, d.real*theta, lift
    labels = ([f'radial_mode_{i}' for i in range(count)]
        +[f'phase_mode_{i}' for i in range(count)]
        +[f'voltage_mode_{i}' for i in range(count)]+['voltage_port_lift'])
    return x, y, voltage, theta, c, s, labels


def energy_query(payload):
    """Return per-energy arrays; the caller performs signed-energy quadrature.

    anomaly_correction has shape (N,2,F), includes nodal mass, and is rotated
    into local radial/tangent quadratures. charge_density=(Tr deltaK)/4 excludes
    mass. current_correction=2*delta_flux_T includes edge conductance. Integrate
    these three with signed energy weight/2 to match existing conventions.
    """
    start = time.monotonic()
    graph, d, epsilon, a_anchor, b_anchor, alpha, archive, sha = _reference(payload)
    modes, lift, gram_error = _response_basis(payload, archive, graph)
    energy = float(payload.get('energy', payload.get('E')))
    omega, eta, temperature = [float(payload[name]) for name in ('omega', 'eta', 't')]
    tolerance = float(payload.get('spectral_tolerance', payload.get('spectraltol', 1e-8)))
    steps = int(payload.get('continuation_steps', 12))
    if not np.isfinite([energy, omega, eta, temperature, tolerance]).all() or min(omega, eta, temperature, tolerance) <= 0 or steps < 1:
        raise ValueError('Finite signed energy and positive frequency, eta, temperature, tolerance, steps required')
    fixed = graph.boundary_nodes
    if fixed is None or not len(fixed):
        raise ValueError('Explicit BCS reservoir nodes required')
    fields, continuation = {}, {}
    for label, value in [('plus', energy+omega/2), ('minus', energy-omega/2), ('static', energy)]:
        R, A, history = _continued_spectrum(graph, d, epsilon, a_anchor, b_anchor,
            alpha, value, eta, tolerance, steps)
        fields[label] = (R, A, *_generators(d, value, eta))
        continuation[label] = history
    Rp, Ap, BRp, BAp, Bp = fields['plus']
    Rm, Am, BRm, BAm, Bm = fields['minus']
    R0, A0, BR0, BA0, B0 = fields['static']
    options = dict(alpha=alpha, fixed_nodes=fixed,
        background_tolerance=max(1e-7, 20*tolerance), normalization_tolerance=1e-8)
    factors = [HarmonicUsadelFactor(graph, left, right, Bleft, Bright, **options)
        for left, right, Bleft, Bright in ((Rp, Rm, BRp, BRm), (Ap, Am, BAp, BAm),
                                         (R0, R0, BR0, BR0), (A0, A0, BA0, BA0))]
    hp, hm, h0 = [float(np.tanh(value/(2*temperature))) for value in (energy+omega/2, energy-omega/2, energy)]
    x, y, voltage, theta, cosine, sine, labels = _column_fields(d, modes, lift, omega)
    n, columns = x.shape
    zero = np.zeros((n, 2, 2), complex)
    base_kinetic = kinetic.assemble(graph, Rp, Rm, Ap, Am, Bp, Bm, zero, zero,
        zero, hp, hm, alpha=alpha)
    free_nodes = np.ones(n, bool); free_nodes[fixed] = False
    free_components = np.repeat(free_nodes, 2)
    kinetic_lu = splu(base_kinetic.matrix[free_components][:, free_components].tocsc())
    anomalous = np.empty((n, 2, columns), complex)
    charge = np.empty((n, columns), complex)
    current = np.empty((len(graph.edges), columns), complex)
    gap_dynamic = np.empty_like(anomalous)
    gap_static = np.empty_like(anomalous)
    current_dynamic, current_static = np.empty_like(current), np.empty_like(current)
    distributions = np.empty((n, 2, columns), complex)
    metrics = []
    response_tolerance = max(1e-7, 20*tolerance)
    for column in range(columns):
        dynamic_B = kinetic.gap_vertex(x[:, column], y[:, column], voltage[:, column])
        static_B = kinetic.gap_vertex(x[:, column], y[:, column])
        is_port = column == columns-1
        fixed_R = lift[fixed, None, None]/omega*(Z@Rm[fixed]-Rp[fixed]@Z) if is_port else 0.
        fixed_A = lift[fixed, None, None]/omega*(Z@Am[fixed]-Ap[fixed]@Z) if is_port else 0.
        fixed_R0 = .5j*theta[fixed, None, None]*(Z@R0[fixed]-R0[fixed]@Z) if is_port else 0.
        fixed_A0 = .5j*theta[fixed, None, None]*(Z@A0[fixed]-A0[fixed]@Z) if is_port else 0.
        R_response = factors[0].solve(dynamic_B, fixed_R, residual_tolerance=response_tolerance)
        A_response = factors[1].solve(dynamic_B, fixed_A, residual_tolerance=response_tolerance)
        if np.max(abs(static_B), initial=0) == 0:
            static_R = static_A = zero
            static_metrics = []
        else:
            R_static_response = factors[2].solve(static_B, fixed_R0, residual_tolerance=response_tolerance)
            A_static_response = factors[3].solve(static_B, fixed_A0, residual_tolerance=response_tolerance)
            static_R, static_A = R_static_response.delta_R, A_static_response.delta_R
            static_metrics = [R_static_response.metrics, A_static_response.metrics]
        dynamic_operator = kinetic.assemble(graph, Rp, Rm, Ap, Am, Bp, Bm,
            R_response.delta_R, A_response.delta_R, dynamic_B, hp, hm, alpha=alpha)
        h = np.zeros((n, 2), complex)
        if is_port:
            h[fixed, 1] = -lift[fixed]*(hp-hm)/omega
        rhs = -dynamic_operator.source.ravel()-base_kinetic.matrix@h.ravel()
        h.ravel()[free_components] = kinetic_lu.solve(rhs[free_components])
        if np.any(~np.isfinite(h)):
            raise RuntimeError('Nonfinite biased kinetic response; no cutoff or retry applied')
        dynamic = dynamic_operator.evaluate(h)
        scaled_residual = float(np.max(abs(dynamic.residual[free_nodes]/graph.area_weights[free_nodes, None]), initial=0))
        source_scale = max(1., float(np.max(abs(dynamic_operator.source/graph.area_weights[:, None]), initial=0)))
        if scaled_residual > response_tolerance*source_scale:
            raise RuntimeError(f'Biased kinetic residual {scaled_residual} exceeds {response_tolerance*source_scale}')
        static_operator = kinetic.assemble(graph, R0, R0, A0, A0, B0, B0,
            static_R, static_A, static_B, h0, h0, alpha=alpha)
        static = static_operator.evaluate(np.zeros((n, 2), complex))
        for destination, observation in ((gap_dynamic, dynamic), (gap_static, static)):
            gx, gy = observation.gap_cartesian.T
            destination[:, 0, column] = cosine*gx+sine*gy
            destination[:, 1, column] = -sine*gx+cosine*gy
        anomalous[:, :, column] = gap_dynamic[:, :, column]-gap_static[:, :, column]
        charge[:, column] = kinetic.project(dynamic.delta_K)[:, 0]
        current_dynamic[:, column], current_static[:, column] = dynamic.charge_current, static.charge_current
        current[:, column] = dynamic.charge_current-static.charge_current
        distributions[:, :, column] = h
        spectral_metrics = [R_response.metrics, A_response.metrics, *static_metrics]
        metrics.append({'column': column, 'label': labels[column],
            'kinetic_scaled_residual': scaled_residual, 'kinetic_source_scale': source_scale,
            'maximum_distribution_response': float(np.max(abs(h))),
            'full_spectral_scaled_max': max(item['full_spectral_scaled_max'] for item in spectral_metrics),
            'normalization_absolute_max': max(item['normalization_absolute_max'] for item in spectral_metrics),
            'static_kinetic_projected_residual': float(np.max(abs(static.residual[free_nodes]/graph.area_weights[free_nodes, None]), initial=0))})
    metadata = dict(schema='pysnspd.stage4.biased_harmonic_energy.v1',
        reference_sha256=sha, graph_signature=graph.signature, energy=energy, omega=omega,
        eta=eta, temperature_over_Tc=temperature, nodes=n, edges=len(graph.edges),
        modes=modes.shape[1], columns=labels, continuation=continuation,
        response_metrics=metrics, mode_gram_error=gram_error,
        elapsed_seconds=time.monotonic()-start,
        time_convention='exp(-i omega t); Omega=hbar*omega/(kB*Tc)',
        quadrature_convention='Signed energy weight/2 for anomaly/current arrays and charge_density=Tr(K)/4',
        scope='Linear response on supplied biased stationary gap, no imposed heat bath, no energy integration in this worker')
    return dict(anomaly_correction=anomalous,
        anomaly_density_correction=anomalous/graph.area_weights[:, None, None],
        charge_density=charge, charge=charge, current_correction=current,
        anomaly_dynamic=gap_dynamic, anomaly_static=gap_static,
        current_dynamic=current_dynamic, current_static=current_static,
        distributions=distributions, metadata=metadata)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--payload', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Distinct output required; no overwrite')
    result = energy_query(json.loads(args.payload.read_text(encoding='utf8')))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output, **{key: value for key, value in result.items() if key != 'metadata'})
    args.output.with_suffix('.json').write_text(json.dumps(result['metadata'], indent=2)+'\n', encoding='utf8')
    print(json.dumps(result['metadata'], indent=2), flush=True)


if __name__ == '__main__':
    main()
