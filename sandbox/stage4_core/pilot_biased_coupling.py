"""Bounded full-mesh API/runtime pilot; its reference is deliberately UNIFORM.

This does not accept a biased physical branch or an energy quadrature. It runs
one complete energy worker on the admitted 1712-node mesh, with four spatial
modes, and one thermal tangent query. Use an external 60--120 s timeout.
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
sys.path.insert(0, str(ROOT))
from sandbox.stage4_core.parallel_runtime import limit_thread_environment, linux_resources, resource_budget
limit_thread_environment()
import numpy as np
from sandbox.stage4_core.coupled_response import mesh_modes, atomic_json
from sandbox.stage4_core.dual_kwt_time import uniform_reference
from sandbox.stage4_core.biased_harmonic_query import energy_query
from sandbox.stage4_core.biased_coupled_response import thermal_query
from pysnspd.experimental import thermal_spatial_usadel as thermal

DEFAULT_MESH = 'docs/implementation/stage4/practical_time_review_20260924/dual_mesh/resampled/mesh.npz'


def run_pilot(output, mesh=ROOT/DEFAULT_MESH, *, tiny=False, continuation_steps=6):
    output = Path(output)
    if output.exists():
        raise FileExistsError('Fresh pilot output required; no overwrite')
    output.mkdir(parents=True)
    started = time.monotonic()
    if hasattr(os, 'sched_getaffinity'):
        resources = linux_resources()
        budget = resource_budget(resources, max_workers=1)
        os.sched_setaffinity(0, {budget['coordinator_cpu']})
        runtime = dict(resources=resources, budget=budget, actual_affinity=sorted(os.sched_getaffinity(0)),
            numerical_workers=1, numerical_library_threads=1)
    else:
        runtime = dict(platform=sys.platform, numerical_workers=1,
            numerical_library_threads=1, affinity_inventory='Linux inventory unavailable on this local platform')
    mesh = Path(mesh)
    if tiny:
        base = thermal.rectangular_graph(5, 3, 2., 1.)
        x = base.coordinates_bar[:, 0]
        fixed = np.flatnonzero((x == x.min()) | (x == x.max()))
        mesh = output/'tiny_mesh.npz'
        np.savez_compressed(mesh, area_weights=base.area_weights, edges=base.edges,
            conductance=base.conductance, coordinates_bar=base.coordinates_bar,
            fixed_nodes=fixed, smooth_profile=np.sin(np.pi*x/2))
    graph, values, modes, lift, conductance, eigen_error = mesh_modes(mesh, 4)
    if not tiny and graph.n_nodes != 1712:
        raise ValueError('Full pilot must use the admitted 1712-node mesh')
    t = .9/8.65
    d, force, epsilon = uniform_reference(graph, t, 256)
    u = d[None]/epsilon[:, None]
    g = 1/np.sqrt(1+abs(u)**2)
    f = u*g
    reference_path, modes_path = output/'uniform_reference.npz', output/'modes.npz'
    np.savez_compressed(reference_path, delta_bar=d, epsilon_bar=epsilon, u=u, f=f, g=g,
        area_weights=graph.area_weights, edges=graph.edges, conductance=graph.conductance,
        coordinates_bar=graph.coordinates_bar, boundary_nodes=graph.boundary_nodes,
        alpha=np.zeros(len(graph.edges)))
    np.savez_compressed(modes_path, modes=modes, lift=lift, eigenvalues=values)
    sources = ['sandbox/stage4_core/pilot_biased_coupling.py',
        'sandbox/stage4_core/biased_harmonic_query.py',
        'sandbox/stage4_core/biased_coupled_response.py',
        'pysnspd/experimental/harmonic_usadel.py',
        'pysnspd/experimental/harmonic_kinetic_usadel.py',
        'pysnspd/experimental/thermal_contact_tangent.py',
        'pysnspd/experimental/retarded_spatial_usadel.py']
    payload = dict(reference_npz=str(reference_path.resolve()), modes_npz=str(modes_path.resolve()),
        reference_sha256=hashlib.sha256(reference_path.read_bytes()).hexdigest(),
        energy=.5, omega=.2, eta=.01*float(d[0].real), t=t,
        spectral_tolerance=1e-8, continuation_steps=continuation_steps)
    manifest = dict(schema='pysnspd.stage4.biased_api_pilot.v1',
        status='RUNNING', scope='Uniform reference; full-mesh API/performance pilot, NOT biased physical acceptance',
        mesh=str(mesh), mesh_sha256=hashlib.sha256(mesh.read_bytes()).hexdigest(),
        nodes=graph.n_nodes, modes=4, matsubara_count=256,
        graph_conductance=conductance, eigen_residual=eigen_error,
        maximum_uniform_force=float(np.max(abs(force))), runtime=runtime,
        source_sha256={path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sources},
        payload=payload)
    atomic_json(output/'manifest.json', manifest)
    print(json.dumps(dict(event='PILOT_START', nodes=graph.n_nodes, columns=13,
        energy_queries=1, thermal_queries=1, scope=manifest['scope'])), flush=True)
    thermal_start = time.monotonic()
    thermal_force, thermal_current, thermal_residual = thermal_query(dict(payload, n=0))
    thermal_seconds = time.monotonic()-thermal_start
    print(json.dumps(dict(event='THERMAL_QUERY_COMPLETE', elapsed_seconds=thermal_seconds,
        maximum_equation_residual=thermal_residual)), flush=True)
    result = energy_query(payload)
    np.savez_compressed(output/'energy_response.npz',
        **{name:value for name,value in result.items() if name != 'metadata'})
    atomic_json(output/'energy_metadata.json', result['metadata'])
    np.savez_compressed(output/'thermal_response.npz', force=thermal_force, current=thermal_current)
    summary = dict(status='PILOT_COMPLETE', scope=manifest['scope'], nodes=graph.n_nodes,
        modes=4, columns=13, thermal_query_seconds=thermal_seconds,
        thermal_query_residual=thermal_residual,
        energy_query_seconds=result['metadata']['elapsed_seconds'],
        maximum_spectral_residual=max(row['full_spectral_scaled_max'] for row in result['metadata']['response_metrics']),
        maximum_kinetic_residual=max(row['kinetic_scaled_residual'] for row in result['metadata']['response_metrics']),
        elapsed_seconds=time.monotonic()-started)
    atomic_json(output/'summary.json', summary)
    manifest.update(status='PILOT_COMPLETE', elapsed_seconds=summary['elapsed_seconds'])
    atomic_json(output/'manifest.json', manifest)
    print(json.dumps(summary, indent=2), flush=True)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--mesh', type=Path, default=ROOT/DEFAULT_MESH)
    parser.add_argument('--tiny', action='store_true', help='Local 15-node API smoke only')
    parser.add_argument('--continuation-steps', type=int, default=6)
    args = parser.parse_args()
    run_pilot(args.output, args.mesh, tiny=args.tiny, continuation_steps=args.continuation_steps)


if __name__ == '__main__':
    main()
