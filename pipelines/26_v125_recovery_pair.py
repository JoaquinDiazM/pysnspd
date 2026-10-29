"""Two independent one-thread electrical continuations with shared CPU leases.

The 0D phases are cheap, but retain independent logs, checkpoints and energy
ledgers. --execute is required; default behavior is authenticated preparation.
"""
import argparse
import contextlib
from concurrent.futures import ThreadPoolExecutor
import fcntl
import importlib.util
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
SCRATCH = Path('/home/jdiaz/scratch')
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
            'NUMEXPR_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[key] = '1'
from pysnspd.experimental.electrical_continuation import prepare, scratch, sha, write_json, source_lock
from sandbox.stage4_core.parallel_runtime import linux_resources, resource_budget

NAMES = ('rs608_bias15p5_775', 'rs608_bias15p5_1550')


def pid_start(pid):
    try:
        return Path(f'/proc/{pid}/stat').read_text().rsplit(')', 1)[1].split()[19]
    except OSError:
        return None


@contextlib.contextmanager
def reservation(campaign):
    lockpath = SCRATCH/'v124_resource_admission.lock'
    registry = SCRATCH/'v124_resource_leases.json'
    def live():
        rows = json.loads(registry.read_text()) if registry.exists() else []
        return [r for r in rows if pid_start(r['pid']) == r.get('start_ticks')]
    with lockpath.open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        rows = live(); res = linux_resources(); budget = resource_budget(res, .9, .9, 1)
        occupied = {c for r in rows for c in r['all_cpus']}
        allowed = set(budget['selected_affinity_cpus'])
        # Legacy independent workers without leases are counted too.
        spec = importlib.util.spec_from_file_location('_v125_resource_scan', ROOT/'pipelines/15_calibration_basecases.py')
        util = importlib.util.module_from_spec(spec); spec.loader.exec_module(util)
        active = util.active_simulations()
        for r in active:
            occupied.update(r['affinity_cpus'])
        if not occupied <= allowed:
            raise ValueError('Other work occupies CPUs outside the admitted ceiling; preserved')
        free = [list(g) for g in res['physical_core_groups'] if set(g) <= allowed and not occupied.intersection(g)]
        if len(free) < 3:
            raise ValueError('Need two whole worker cores and one coordinator core; jobs preserved')
        chosen = free[:3]; mine_cpus = sorted(c for g in chosen for c in g)
        remaining = [g for g in res['physical_core_groups'] if not set(g).intersection(occupied | set(mine_cpus))]
        memory = 3*1024**3
        if (len(remaining) < 2 or len(occupied | set(mine_cpus)) > budget['logical_budget'] or
            memory+sum(r.get('memory_reservation_bytes', 0) for r in rows) > budget['memory_limit_bytes']):
            raise ValueError('Shared CPU/two-whole-core/RAM budget unavailable')
        mine = dict(pid=os.getpid(), start_ticks=pid_start(os.getpid()), campaign_dir=str(campaign),
                    all_cpus=mine_cpus, memory_reservation_bytes=memory)
        write_json(registry, rows+[mine])
        affinity = set(os.sched_getaffinity(0)); os.sched_setaffinity(0, set(chosen[2]))
        allocation = dict(inventory=res, budget=budget, worker_cpus=[g[0] for g in chosen[:2]],
            coordinator_cpus=chosen[2], reserved_whole_core_groups=remaining, worker_threads=2,
            combined_logical_cpus=sorted(occupied | set(mine_cpus)), prior_leases=rows)
    try:
        yield allocation
    finally:
        os.sched_setaffinity(0, affinity)
        with lockpath.open('a+') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            write_json(registry, [r for r in live() if (r['pid'], r['start_ticks']) != (mine['pid'], mine['start_ticks'])])


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-campaign', type=Path, required=True)
    p.add_argument('--campaign-dir', type=Path, required=True)
    p.add_argument('--final-time-ns', type=float, default=100.)
    p.add_argument('--execute', action='store_true'); p.add_argument('--resume', action='store_true')
    a = p.parse_args(argv); source = scratch(a.source_campaign); out = scratch(a.campaign_dir)
    if source == out or out.is_relative_to(source):
        raise ValueError('A separate 0D coordinator directory is required')
    out.mkdir(parents=True, exist_ok=True)
    with (out/'.execution.lock').open('a+') as local:
        try:
            fcntl.flock(local, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise ValueError('This campaign already has a coordinator')
        plans = []
        for name in NAMES:
            origin = scratch(source/'runs'/name)
            with source_lock(origin):
                plans.append(prepare(origin, out/'runs'/name, a.final_time_ns))
        plan = dict(model_version='1.2.5', source_campaign=str(source), final_time_ns=a.final_time_ns,
                    cases=list(NAMES), source_pins_sha256={'pipelines/26_v125_recovery_pair.py': sha(__file__)})
        old = out/'campaign_plan.json'
        if old.exists() and json.loads(old.read_text()) != plan:
            raise ValueError('Campaign identity differs; use a new directory')
        write_json(old, plan)
        for name, phase in zip(NAMES, plans):
            ad = phase['admission']
            print(f'{name}: '+('ADMITTED' if ad['admitted'] else 'NOT ADMITTED')+
                  f" | native {phase['initial_time_ps']/1000:.6f}ns | overlap {ad['checks']}checks | "+
                  f"max output error {ad['errors']['voltage_error_V']*1e6:.6f}uV", flush=True)
        if not all(v['admission']['admitted'] for v in plans):
            return 3
        if not a.execute:
            print('PREPARED; no electrical continuation executed.'); return 0
        if (out/'campaign_execution.json').exists() and not a.resume:
            raise ValueError('Existing execution: use --resume')
        processes = []; stopped = {'value': False}
        def stop(sig, frame):
            stopped['value'] = True
            for proc in processes:
                if proc.poll() is None:
                    proc.send_signal(signal.SIGTERM)
        signal.signal(signal.SIGINT, stop); signal.signal(signal.SIGTERM, stop)
        with reservation(out) as alloc:
            write_json(out/'execution_resources.json', alloc)
            print('PYSNSPD1.2.5 |2 electrical 0D cases |1 thread/case |100ns global horizon', flush=True)
            logs = out/'logs'; logs.mkdir(exist_ok=True)
            def worker(name, cpu):
                cmd = [sys.executable, '-u', str(ROOT/'pipelines/25_v125_electrical_phase.py'),
                    '--from-run', str(source/'runs'/name), '--output-dir', str(out/'runs'/name),
                    '--final-time-ns', str(a.final_time_ns), '--execute', '--cpu-list', str(cpu)]
                if a.resume:
                    cmd.append('--resume')
                env = dict(os.environ, PYTHONPYCACHEPREFIX=str(SCRATCH/'tmp/pysnspd/pycache'),
                    NUMBA_CACHE_DIR=str(SCRATCH/'tmp/pysnspd/numba_cache'), TMPDIR=str(SCRATCH/'tmp/pysnspd'))
                with (logs/(name+'.log')).open('a') as log:
                    proc = subprocess.Popen(['taskset', '-c', str(cpu)]+cmd, cwd=ROOT, env=env,
                                            stdout=log, stderr=subprocess.STDOUT)
                    processes.append(proc)
                    if stopped['value']:
                        proc.send_signal(signal.SIGTERM)
                    code = proc.wait()
                path = out/'runs'/name/'manifest.json'
                m = json.loads(path.read_text()) if path.exists() else {}
                return dict(case=name, exit_code=code, status=m.get('status', 'FAILED'),
                            time_ps=m.get('time_ps'), run_complete=m.get('run_complete', False))
            with ThreadPoolExecutor(max_workers=2) as pool:
                futures = [pool.submit(worker, name, cpu) for name, cpu in zip(NAMES, alloc['worker_cpus'])]
                while not all(f.done() for f in futures):
                    rows = []
                    for name in NAMES:
                        path = out/'runs'/name/'manifest.json'
                        try:
                            m = json.loads(path.read_text())
                            rows.append(f"{name}: {m['status']} phase={m['phase']} t={m['time_ps']/1000:.6f}ns")
                        except (FileNotFoundError, json.JSONDecodeError):
                            rows.append(name+': authenticating native checkpoint')
                    print('\n'.join(rows), flush=True)
                    time.sleep(2.)
                results = [f.result() for f in futures]
            write_json(out/'campaign_execution.json', dict(total=2, finished=2, results=results))
            for r in results:
                print(json.dumps(r), flush=True)
            return 0 if all(r['exit_code'] == 0 for r in results) else 2


if __name__ == '__main__':
    raise SystemExit(main())
