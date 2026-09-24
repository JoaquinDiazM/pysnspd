"""Certify and compact a completed thermal trajectory, retaining its failed gate.

No numerical physics is executed. Every observation and accepted checkpoint is
verified against the execution records; the original status is never rewritten.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
from parallel_runtime import limit_thread_environment, linux_resources, resource_budget
limit_thread_environment()
import numpy as np


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024**2), b''):
            digest.update(block)
    return digest.hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf8'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    if args.output.exists():
        raise FileExistsError('Fresh review export required; no overwrite')
    identity, summary, plan = [read(args.raw/n) for n in ('identity.json', 'summary.json', 'executed_plan.json')]
    assert sha(args.raw/'executed_plan.json') == identity['plan_sha256']
    assert summary['horizon_ps'] == plan['observation_times_ps'][-1]
    assert summary['refinement'] == read(args.raw/'refinement.json')
    for name, digest in identity['sources'].items():
        assert sha(ROOT/name) == digest, ('source changed', name)
    operator = Path(identity['operator_root'])
    for name, digest in identity['operator_inputs'].items():
        assert sha(operator/name) == digest, ('operator changed', name)

    passes = {name: read(args.raw/name/'summary.json') for name in plan['passes']}
    expected, accepted = {}, []
    for name, data in passes.items():
        assert [row['time_ps'] for row in data['observations']] == plan['observation_times_ps']
        for row in data['observations']:
            expected[row['fields_path']] = row['fields_sha256']
        rows = [row for row in data['history'] if row['accepted']]
        assert len(list((args.raw/name).glob('accepted_*.npz'))) == len(rows)
        for row in rows:
            path = args.raw/name/row['fields_path']
            expected[path.relative_to(args.raw).as_posix()] = row['fields_sha256']
            with np.load(path) as arrays:
                assert float(arrays['time_ps']) == row['time_ps']
                maximum = float(np.max(abs(arrays['displacements']))/plan['gap_reference_kBTc'])
                assert maximum == row['maximum_relative_displacement']
                assert maximum <= plan['maximum_relative_displacement']
            accepted.append(dict(pass_name=name, **row))
    budget = resource_budget(linux_resources(), .9, .9, max_workers=4)
    def verify(path):
        name = path.relative_to(args.raw).as_posix()
        assert name in expected, ('unregistered checkpoint', name)
        actual = sha(path)
        assert actual == expected[name], ('checkpoint changed', name)
        with np.load(path) as arrays:
            assert all(np.all(np.isfinite(arrays[key])) for key in arrays.files)
            shapes = {key: list(arrays[key].shape) for key in arrays.files}
        return dict(path=name, sha256=actual, bytes=path.stat().st_size,
                    recorded_hash_verified=True, arrays_finite=True, shapes=shapes)
    affinity = os.sched_getaffinity(0)
    os.sched_setaffinity(0, set(budget['worker_affinity_cpus']+[budget['coordinator_cpu']]))
    try:
        with ThreadPoolExecutor(max_workers=budget['workers']) as pool:
            certificate = list(pool.map(verify, sorted(args.raw.rglob('*.npz'))))
    finally:
        os.sched_setaffinity(0, affinity)
    common = ('d0', 'G0', 'coordinates_bar', 'area_weights', 'edges', 'conductance', 'boundary_nodes')
    first = args.raw/passes['primary']['observations'][0]['fields_path']
    with np.load(first) as arrays:
        compact = {key: arrays[key].copy() for key in common}
        fields = sorted(set(arrays.files)-set(common))
    for name, data in passes.items():
        all_fields = {key: [] for key in fields}
        for row in data['observations']:
            with np.load(args.raw/row['fields_path']) as arrays:
                assert all(np.array_equal(arrays[key], compact[key]) for key in common)
                for key in fields:
                    all_fields[key].append(arrays[key].copy())
        compact.update({name+'_'+key: np.asarray(values) for key, values in all_fields.items()})
        compact[name+'_times_ps'] = np.asarray([row['time_ps'] for row in data['observations']])
    args.output.mkdir(parents=True)
    for name in ('identity.json','summary.json','executed_plan.json','refinement.json','progress.jsonl','failure.json'):
        if (args.raw/name).exists():
            shutil.copyfile(args.raw/name, args.output/name)
    for name in passes:
        (args.output/name).mkdir()
        shutil.copyfile(args.raw/name/'summary.json', args.output/name/'summary.json')
    np.savez_compressed(args.output/'trajectory_compact.npz', **compact)
    copied = {p.relative_to(args.output).as_posix(): sha(p) for p in args.output.rglob('*') if p.is_file()}
    receipt = dict(status='CHECKPOINTS_VERIFIED_ORIGINAL_GATE_RETAINED',
        original_status=summary['status'], raw_directory=str(args.raw.resolve()),
        runtime_seconds=time.monotonic()-started, script_sha256=sha(__file__),
        sources_matching_identity=identity['sources'], operator_inputs_verified=len(identity['operator_inputs']),
        budget=budget, checkpoint_certificate=certificate, copied_file_sha256=copied,
        accepted_steps=accepted, observed_fields=fields,
        maximum_accepted_relative_displacement=max(row['maximum_relative_displacement'] for row in accepted),
        new_physics_solves=0, originals_retained_remotely=True)
    with (args.output/'extraction_receipt.json').open('x', encoding='utf8') as stream:
        json.dump(receipt, stream, indent=2, allow_nan=False); stream.write('\n')
    print(json.dumps(dict(status=receipt['status'], original_status=receipt['original_status'],
        checkpoints=len(certificate), accepted_steps=len(accepted),
        compact_bytes=(args.output/'trajectory_compact.npz').stat().st_size,
        maximum_displacement_percent=100*receipt['maximum_accepted_relative_displacement'],
        runtime_seconds=receipt['runtime_seconds'])))


if __name__ == '__main__':
    main()
