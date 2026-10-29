"""Append an authenticated electrical-only 0D phase to a native 1D checkpoint.

Prepare-only is the default. No historical field, population, source or ledger
is rewritten. The full native 2D/1D producers remain source-frozen.
"""
import argparse
from dataclasses import asdict
import json
import os
from pathlib import Path
import signal
import sys

for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
            'NUMEXPR_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[key] = '1'
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
from pysnspd.experimental.electrical_0d import AdmissionConfig
from pysnspd.experimental.electrical_continuation import prepare, execute, scratch, source_lock


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--from-run', type=Path, required=True)
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--final-time-ns', type=float, default=100.)
    p.add_argument('--sample-ps', type=float, default=10.)
    p.add_argument('--execute', action='store_true')
    p.add_argument('--resume', action='store_true')
    p.add_argument('--cpu-list', type=str)
    for key, value in asdict(AdmissionConfig()).items():
        p.add_argument('--admit-'+key.replace('_', '-'), type=float, default=value)
    return p


def main(argv=None):
    p = parser(); a = p.parse_args(argv)
    if a.execute and not a.cpu_list:
        p.error('Execution requires an admitted --cpu-list; the pair coordinator allocates it automatically')
    if a.cpu_list:
        cpus = set(map(int, a.cpu_list.split(',')))
        if not cpus or not cpus <= os.sched_getaffinity(0):
            raise ValueError('CPU mask outside admitted affinity')
        os.sched_setaffinity(0, cpus)
    config = AdmissionConfig(**{k: getattr(a, 'admit_'+k) for k in asdict(AdmissionConfig())})
    source = scratch(a.from_run)
    with source_lock(source):
        plan = prepare(source, a.output_dir, a.final_time_ns, a.sample_ps, config)
    print(json.dumps(dict(status='ELECTRICAL_0D_ADMITTED' if plan['admission']['admitted']
                         else 'ELECTRICAL_0D_NOT_ADMITTED', admission=plan['admission'])), flush=True)
    if not a.execute:
        return 0 if plan['admission']['admitted'] else 3
    stopped = {'value': False}
    def stop(sig, frame):
        stopped['value'] = True
    signal.signal(signal.SIGINT, stop); signal.signal(signal.SIGTERM, stop)
    m = execute(a.output_dir, resume=a.resume, stop_requested=lambda: stopped['value'])
    print(json.dumps(dict(status=m['status'], phase='0d', time_ps=m['time_ps'],
                         output_dir=str(a.output_dir))), flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
