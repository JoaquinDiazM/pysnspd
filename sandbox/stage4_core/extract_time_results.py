"""Certify completed affine thermal trajectories without executing any physics.

Original checkpoints stay on Geminga. A compact archive removes repeated mesh
data while retaining every observed field for independent analysis and figures.
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


def write(path, data):
    with Path(path).open('x', encoding='utf8') as stream:
        json.dump(data, stream, indent=2, allow_nan=False)
        stream.write('\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    if args.output.exists():
        raise FileExistsError('Fresh review export required; no overwrite')
    identity, summary, plan = [read(args.raw / n) for n in ('identity.json', 'summary.json', 'executed_plan.json')]
    assert summary['status'] == 'FULL_NODE_AFFINE_THERMAL_TIME_COMPLETE'
    assert sha(args.raw / 'executed_plan.json') == identity['plan_sha256']
    assert summary['horizon_ps'] == plan['observation_times_ps'][-1]
    assert summary['refinement'] == read(args.raw / 'refinement.json')
    assert summary['refinement']['all_comparisons_met']
    assert not (args.raw / 'failure.json').exists()
    for name, digest in identity['sources'].items():
        assert sha(ROOT / name) == digest, ('source changed', name)
    operator = Path(identity['operator_root'])
    for name, digest in identity['inputs'].items():
        assert sha(operator / name) == digest, ('operator changed', name)
    for name in ('summary', 'identity'):
        assert sha(operator / (name + '.json')) == plan['operator_' + name + '_sha256']

    passes = {name: read(args.raw / name / 'summary.json') for name in plan['passes']}
    jobs = []
    expected = {}
    for name, data in passes.items():
        assert [row['time_ps'] for row in data['observations']] == plan['observation_times_ps']
        assert all(row['admitted'] for row in summary['refinement']['records'])
        for row in data['observations']:
            expected[row['fields_path']] = row['fields_sha256']
    all_npz = sorted(args.raw.rglob('*.npz'))
    budget = resource_budget(linux_resources(), .9, .9, max_workers=4)

    def verify(path):
        name = path.relative_to(args.raw).as_posix()
        actual = sha(path)
        if name in expected:
            assert actual == expected[name], ('observation hash mismatch', name)
        with np.load(path) as arrays:
            assert all(np.all(np.isfinite(arrays[key])) for key in arrays.files), name
            shapes = {key: list(arrays[key].shape) for key in arrays.files}
        return dict(path=name, sha256=actual, bytes=path.stat().st_size,
                    recorded_observation_hash_verified=name in expected,
                    arrays_finite=True, shapes=shapes)

    previous_affinity = os.sched_getaffinity(0)
    os.sched_setaffinity(0, set(budget['worker_affinity_cpus'] + [budget['coordinator_cpu']]))
    try:
        with ThreadPoolExecutor(max_workers=budget['workers']) as pool:
            certificate = list(pool.map(verify, all_npz))
    finally:
        os.sched_setaffinity(0, previous_affinity)

    common = ('d0', 'G0', 'coordinates_bar', 'area_weights', 'edges', 'conductance', 'boundary_nodes')
    compact = {}
    first = args.raw / passes['primary']['observations'][0]['fields_path']
    with np.load(first) as arrays:
        compact.update({key: arrays[key].copy() for key in common})
        fields = sorted(set(arrays.files) - set(common))
    free = np.ones(len(compact['d0']), bool)
    free[compact['boundary_nodes']] = False
    scale = np.sqrt(compact['area_weights'][free]) / plan['gap_reference_kBTc']
    size = 2 * np.count_nonzero(free)
    accepted_records = []
    for name, data in passes.items():
        for key in fields:
            values = []
            for row in data['observations']:
                with np.load(args.raw / row['fields_path']) as arrays:
                    assert all(np.array_equal(arrays[k], compact[k]) for k in common)
                    values.append(arrays[key].copy())
            compact[name + '_' + key] = np.asarray(values)
        compact[name + '_times_ps'] = np.asarray([row['time_ps'] for row in data['observations']])
        accepted = [(i, row) for i, row in enumerate(data['history'], 1) if row['accepted']]
        assert len(list((args.raw / name).glob('accepted_*.npz'))) == len(accepted)
        for index, row in accepted:
            path = args.raw / name / f'accepted_{index:05d}.npz'
            with np.load(path) as arrays:
                vector = arrays['affine_state']
                assert abs(float(arrays['time_ps']) - row['time_ps']) < 1e-13
                # Arnoldi approximates the augmented exponential too. Report
                # drift of its constant coordinate; do not invent a new gate.
                affine_constant_defect = float(vector[-1] - 1)
                raw_blocks = vector[:-1].reshape(3, size)
                blocks = np.zeros((3, len(free)), complex)
                for j, (block, reference) in enumerate(zip(raw_blocks, data['scaled_reference_norms'])):
                    blocks[j, free] = (block[::2] + 1j*block[1::2]) * reference / scale
                states = np.array([blocks[0], blocks[0]+blocks[1], blocks[0]+blocks[2]])
                maximum = float(np.max(abs(states)) / plan['gap_reference_kBTc'])
                assert maximum <= plan['maximum_relative_displacement']
                assert row['integrated_defect'] <= row['tolerance']
                accepted_records.append(dict(pass_name=name, path=path.relative_to(args.raw).as_posix(),
                    time_ps=row['time_ps'], maximum_relative_displacement=maximum,
                    integrated_defect=row['integrated_defect'], tolerance=row['tolerance'],
                    affine_constant_defect=affine_constant_defect))

    args.output.mkdir(parents=True)
    for name in ('identity.json', 'summary.json', 'executed_plan.json', 'refinement.json', 'progress.jsonl'):
        shutil.copyfile(args.raw / name, args.output / name)
    for name in passes:
        (args.output / name).mkdir()
        shutil.copyfile(args.raw / name / 'summary.json', args.output / name / 'summary.json')
    np.savez_compressed(args.output / 'trajectory_compact.npz', **compact)
    copied = {p.relative_to(args.output).as_posix(): sha(p) for p in args.output.rglob('*') if p.is_file()}
    receipt = dict(status='ALL_THERMAL_TIME_CHECKPOINTS_VERIFIED_NO_SOLVES',
        raw_directory=str(args.raw.resolve()), runtime_seconds=time.monotonic()-started,
        script_sha256=sha(__file__), sources_matching_identity=identity['sources'],
        operator_inputs_verified=len(identity['inputs']), budget=budget,
        checkpoint_certificate=certificate, copied_file_sha256=copied,
        accepted_steps=accepted_records, observed_fields=fields,
        maximum_accepted_relative_displacement=max(r['maximum_relative_displacement'] for r in accepted_records),
        maximum_affine_constant_defect=max(abs(r['affine_constant_defect']) for r in accepted_records),
        new_physics_solves=0, originals_retained_remotely=True)
    write(args.output / 'extraction_receipt.json', receipt)
    print(json.dumps(dict(status=receipt['status'], checkpoints=len(certificate),
        accepted_steps=len(accepted_records), raw_bytes=sum(r['bytes'] for r in certificate),
        compact_bytes=(args.output/'trajectory_compact.npz').stat().st_size,
        maximum_displacement_percent=100*receipt['maximum_accepted_relative_displacement'],
        runtime_seconds=receipt['runtime_seconds'])))


if __name__ == '__main__':
    main()
