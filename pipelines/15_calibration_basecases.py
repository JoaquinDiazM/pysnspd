"""Eight foreground v1.2.3 latency cases, with a human-readable terminal board.

Prepare only by default. --execute starts all eight cases after resource admission.
No equations, tolerances or historical data are changed by this coordinator.
"""
from __future__ import annotations

import argparse
from collections import deque
from concurrent.futures import ThreadPoolExecutor
import fcntl
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import signal
import shutil
import subprocess
import sys
import threading
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
SCRATCH = Path('/home/jdiaz/scratch')
PYTHON = '/home/jdiaz/.conda/envs/snspd/bin/python'
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
from sandbox.stage4_core.parallel_runtime import linux_resources, resource_budget


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    temp = path.with_name(path.name + '.tmp')
    temp.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    os.replace(temp, path)


def catalog():
    return [dict(name=f'rs{rs}_bias{str(current).replace(".", "p")}_{wavelength}',
                 sheet_resistance_ohm=rs, current_uA=current, wavelength_nm=wavelength,
                 threads=4 if current == 15.5 and wavelength == 1550 else 3)
            for rs in (425, 608) for current in (15.5, 21.5) for wavelength in (775, 1550)]


def active_simulations():
    """Read existing producer/coordinator affinity without disturbing any process."""
    names = {'05_effective_photon_run.py', '08_effective_basecase_sweep.py',
             '09_thermal_stiffness_photon_run.py', '10_thermal_stiffness_basecases.py',
             '11_bias_escape_controls.py', '12_low_bias_latency_controls.py',
             '13_dimensional_photon_run.py', '14_dimensional_basecases.py', Path(__file__).name}
    found = []
    for p in Path('/proc').iterdir():
        if not p.name.isdigit() or int(p.name) == os.getpid():
            continue
        try:
            if p.stat().st_uid != os.getuid():
                continue
            words = (p/'cmdline').read_bytes().decode(errors='replace').split('\0')
            if any(Path(word).name in names for word in words):
                found.append(dict(pid=int(p.name), affinity_cpus=sorted(os.sched_getaffinity(int(p.name)))))
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            continue
    return found


def allocation(resources=None, active=None):
    res = linux_resources() if resources is None else resources
    budget = resource_budget(res, .9, .9, 1)
    active = active_simulations() if active is None else active
    if active:
        raise ValueError('Active simulation detected; all eight cases require a fresh shared budget. Existing jobs preserved.')
    groups = res['physical_core_groups']
    allowed = set(budget['selected_affinity_cpus'])
    admitted = [list(g) for g in groups if set(g) <= allowed]
    if len(admitted) < 14 or any(len(g) != 2 for g in admitted[:14]):
        raise ValueError('This eight-case layout requires fourteen admitted SMT2 physical cores and two reserved cores.')
    # Two difficult low-current infrared cases receive complete pairs of cores.
    slots = {}
    difficult = [c for c in catalog() if c['threads'] == 4]
    for i, case in enumerate(difficult):
        slots[case['name']] = admitted[2*i] + admitted[2*i+1]
    ordinary = [c for c in catalog() if c['threads'] == 3]
    for pair in range(3):
        a, b, c = admitted[4+3*pair:7+3*pair]
        slots[ordinary[2*pair]['name']] = a + [b[0]]
        slots[ordinary[2*pair+1]['name']] = [b[1]] + c
    coordinator = admitted[13]
    worker = [v for cpus in slots.values() for v in cpus]
    if len(worker) != 26 or len(set(worker)) != 26 or set(worker) & set(coordinator):
        raise ValueError('Overlapping or incomplete logical CPU allocation')
    used = set(worker) | set(coordinator)
    reserved = [list(g) for g in groups if not set(g) & used]
    total = budget['effective_cpu_quota_units']
    if len(reserved) < 2 or len(used) > budget['logical_budget'] or len(used)/total > .9:
        raise ValueError('Global CPU ceiling or two-core reserve violated')
    if len(worker)/total < .8:
        raise ValueError('This layout does not allocate the requested 80% worker capacity on this host')
    memory = (3 + 8*5)*1024**3
    if memory > budget['memory_limit_bytes']:
        raise ValueError('Insufficient RAM for eight conservative 5 GiB worker reservations')
    sharing = []
    for group in groups:
        owners = [name for name, cpus in slots.items() if set(group) & set(cpus)]
        if len(owners) > 1:
            sharing.append(dict(physical_core_cpus=list(group), cases=owners))
    return dict(inventory=res, budget=budget, slots=slots, coordinator_cpus=coordinator,
                reserved_whole_core_groups=reserved, worker_threads=26,
                worker_capacity_fraction=26/total, admitted_capacity_fraction=len(used)/total,
                nested_blas_threads_per_worker=1, shared_SMT_cores=sharing,
                memory_reservation_bytes=memory, active_simulations=active,
                scope='Disjoint logical CPUs; three ordinary pairs share SMT siblings. Allocation is not measured CPU usage or speedup.')


def command(case, out, cpus, args, resume=False):
    words = [PYTHON, '-u', str(ROOT/'pipelines/13_dimensional_photon_run.py'),
             '--output-dir', str(out), '--cpu-list', ','.join(map(str, cpus)),
             '--threads', str(case['threads']), '--physical-cores', '2', '--latency-only',
             '--current-uA', str(case['current_uA']), '--wavelength-nm', str(case['wavelength_nm']),
             '--sheet-resistance-ohm', str(case['sheet_resistance_ohm']),
             '--retained-fraction', '0.9', '--tau-ee-Tc-ps', '5', '--tau-ep-Tc-ps', '24.7',
             '--tau-escape-ps', '20', '--target-depairing-uA', str(38.8*80/120),
             '--stiffness-temperature', 'electron', '--gradient-closure', 'thermal-depairing',
             '--length-nm', '300', '--spacing-nm', '5', '--sigma-nm', '7.549834435',
             '--bath-K', '0.9', '--threshold-mV', '0.1', '--confirmation-ps', '0.5',
             '--extra-time-ps', '2', '--ceiling-ps', str(args.ceiling_ps)]
    if resume:
        words.append('--resume')
    return words


def hours(seconds):
    if seconds is None or not math.isfinite(seconds) or seconds < 0:
        return '--'
    minutes = round(seconds/60)
    return f'{minutes//60}h{minutes%60:02d}'


class Board:
    """Only the main thread writes the terminal; workers keep complete raw logs."""
    def __init__(self, cases, out, ceiling, ansi):
        self.out, self.ceiling, self.ansi = out, ceiling, ansi
        self.lock = threading.Lock()
        self.rows = {c['name']:dict(case=c, status='WAIT', samples=deque(maxlen=8),
                                  t=None, V=None, step=None, steps=None, rejects=None,
                                  latency=None, target=None, pid=None, proc_sample=None, cpu=None,
                                  rss=None, last_message='') for c in cases}
        self.lines = 0
        self.started = time.monotonic()

    def update(self, name, **values):
        with self.lock:
            self.rows[name].update(values)

    def ingest(self, name, line):
        try:
            data = json.loads(line)
        except (ValueError, TypeError):
            with self.lock:
                if line.strip():
                    self.rows[name]['last_message'] = line.strip()[-240:]
            return
        if not isinstance(data, dict):
            return
        with self.lock:
            row = self.rows[name]
            for key, source in [('t', 'time_ps'), ('V', 'Vout_mV'), ('step', 'step_ps'),
                                ('steps', 'accepted_steps'), ('rejects', 'rejections')]:
                if data.get(source) is not None:
                    row[key] = data[source]
            if data.get('status'):
                row['status'] = data['status']
            t = data.get('time_ps')
            if t is not None and (not row['samples'] or t > row['samples'][-1][1]):
                row['samples'].append((time.monotonic(), t))

    def refresh_files(self):
        for name in self.rows:
            try:
                path = self.out/'runs'/name/'manifest.json'
                data = json.loads(path.read_text())
                with self.lock:
                    row = self.rows[name]
                    if row['status'] not in ('STARTING', 'WAIT') and row.get('exit_code') is not None:
                        continue
                    row['status'] = data.get('status', row['status'])
                    if data.get('time_ps') is not None:
                        row['t'] = max(data['time_ps'], row['t'] or 0.)
                    row['steps'] = max(data.get('accepted_steps', 0), row['steps'] or 0)
                    row['rejects'] = data.get('rejections', row['rejects'])
                    obs = data.get('observation', {})
                    monitor = data.get('latest_monitor', {})
                    if monitor.get('Vout_mV') is not None and (row['t'] or 0.) <= data.get('time_ps', 0.):
                        row['V'] = monitor['Vout_mV']
                    crossing = obs.get('crossing_time_s')
                    if crossing is not None:
                        row['latency'] = crossing*1e12
                        row['target'] = row['latency'] + 2.
                    if data.get('error'):
                        row['last_message'] = data['error'][-240:]
            except (FileNotFoundError, json.JSONDecodeError):
                pass
        now = time.monotonic()
        hz = os.sysconf(os.sysconf_names['SC_CLK_TCK'])
        with self.lock:
            for row in self.rows.values():
                if row['pid'] is None or row.get('exit_code') is not None:
                    continue
                try:
                    stat = Path(f'/proc/{row["pid"]}/stat').read_text().rsplit(')', 1)[1].split()
                    total = (int(stat[11]) + int(stat[12]))/hz
                    rss = int(stat[21])*os.sysconf('SC_PAGE_SIZE')/1024**3
                    previous = row['proc_sample']
                    if previous:
                        row['cpu'] = 100*(total-previous[1])/max(now-previous[0], .001)
                    row['proc_sample'] = (now, total)
                    row['rss'] = rss
                except (FileNotFoundError, ProcessLookupError, ValueError):
                    pass

    def render(self, final=False):
        self.refresh_files()
        with self.lock:
            width = shutil.get_terminal_size((140, 40)).columns if self.ansi else 160
            compact = width < 132
            text = [f'PYSNSPD 1.2.3 | 8 latency cases | eta=0.90, tau_ee=5 ps, escape=20 ps | wall {hours(time.monotonic()-self.started)}',
                    ('Case                         State        t[ps]  Vout[mV]    ETA*   CPU%' if compact else
                     'Case                         State        t[ps]  Vout[mV]  dt[ps]  accept/reject  ps/h   ETA*    CPU%   RAM'),
                    '-'*123]
            for row in self.rows.values():
                status = row['status']
                short = {'PREPARATION':'DC/source', 'STARTING':'DC/source', 'RUNNING':'RUN',
                         'TRIGGER_MARGIN_COMPLETE':'DETECTED', 'TIME_CENSORED_NO_TRIGGER':'CENSORED',
                         'INTERRUPTED_CHECKPOINTED':'PAUSED', 'FAILED_BEFORE_MANIFEST':'FAILED'}.get(status, status[:11])
                if row.get('exit_code') not in (None, 0) and status != 'NOT_STARTED':
                    short = 'FAILED'
                speed = None
                if len(row['samples']) >= 2:
                    a, b = row['samples'][0], row['samples'][-1]
                    speed = 3600*(b[1]-a[1])/max(b[0]-a[0], .001)
                target = row['target'] if row['target'] is not None else self.ceiling
                eta = 3600*max(0., target-(row['t'] or 0.))/speed if speed and speed > 0 else None
                if row.get('exit_code') is not None:
                    eta = 0.
                fmt = lambda value, pattern: '--' if value is None else format(value, pattern)
                prefix = (f'{row["case"]["name"]:<28} {short:<11} {fmt(row["t"],"7.2f"):>7} '
                          f'{fmt(row["V"],"8.5f"):>9}')
                if compact:
                    text.append(prefix + f' {hours(eta):>7} {fmt(row["cpu"],"5.0f"):>6}')
                else:
                    text.append(prefix + f' {fmt(row["step"],".1e"):>7} '
                                f'{str(row["steps"] or 0)+"/"+str(row["rejects"] or 0):>13} '
                                f'{fmt(speed,"5.1f"):>5} {hours(eta):>7} '
                                f'{fmt(row["cpu"],"5.0f"):>6} {fmt(row["rss"],".1f"):>5}G')
            text.extend(['-'*123,
                         f'* ETA is rolling time to {self.ceiling:g} ps cap before trigger, then crossing + 2 ps; not a detection forecast.',
                         'CPU% per process: 100% = one logical CPU. 26 Numba threads + coordinator; 2 whole cores reserved.',
                         f'Raw logs and checkpoints: {self.out}'])
            errors = [f'{name}: {row["last_message"]}' for name, row in self.rows.items()
                      if row.get('exit_code') not in (None, 0) and row['last_message']]
            text.extend(errors[:8])
            if self.ansi:
                text = [line[:max(20, width-1)] for line in text]
            if self.ansi and self.lines:
                sys.stdout.write(f'\x1b[{self.lines}A\x1b[J')
            sys.stdout.write('\n'.join(text)+'\n')
            sys.stdout.flush()
            self.lines = len(text)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--campaign-dir', type=Path, required=True)
    mode = p.add_mutually_exclusive_group()
    mode.add_argument('--prepare-only', action='store_true')
    mode.add_argument('--execute', action='store_true')
    p.add_argument('--resume', action='store_true')
    p.add_argument('--ceiling-ps', type=float, default=100.)
    p.add_argument('--refresh-seconds', type=float, default=10.)
    p.add_argument('--plain', action='store_true', help='Print periodic complete tables without ANSI cursor movement')
    args = p.parse_args(argv)
    if not math.isfinite(args.ceiling_ps) or args.ceiling_ps <= 0 or not 1 <= args.refresh_seconds <= 60:
        p.error('Positive finite ceiling and refresh interval in [1,60] seconds required')
    out = args.campaign_dir.expanduser().resolve()
    if not out.is_relative_to(SCRATCH.resolve()):
        p.error('All campaign products must stay under scratch')
    cases = catalog()
    resources = allocation()
    spec = importlib.util.spec_from_file_location('_calibration_v123', ROOT/'pipelines/13_dimensional_photon_run.py')
    run = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(run)
    pins = run.source_pins()
    pins[str(Path(__file__).relative_to(ROOT))] = digest(__file__)
    pins['sandbox/stage4_core/parallel_runtime.py'] = digest(ROOT/'sandbox/stage4_core/parallel_runtime.py')
    definition = dict(schema='pysnspd.calibration_latency_campaign.v1', model_version='1.2.3',
                      cases=cases, retained_fraction=.9, tau_ee_Tc_ps=5., tau_ep_Tc_ps=24.7,
                      tau_escape_ps=20., ceiling_ps=args.ceiling_ps, confirmation_ps=.5,
                      extra_time_ps=2., latency_only=True, source_pins_sha256=pins,
                      material_control='Rs425 versus Rs608; recalibrate D to the same 25.8666667 uA depairing target; derive N0, ell0, gamma, geometry and DC inductance coherently',
                      commands={c['name']:command(c, out/'runs'/c['name'], resources['slots'][c['name']], args)
                                for c in cases},
                      hypotheses='Effective postcascade Gaussian with common initial FD/BE; sigma remains independent of calibrated D; one common material per four-case set; no quantitative fit claimed',
                      scheduling=dict(simultaneous_cases=8, worker_threads=26, threads_per_case=[c['threads'] for c in cases],
                                      SMT_sharing_declared=True, minimum_reserved_whole_cores=2,
                                      coordinator_on_separate_core=True))
    # Runtime affinity is recorded separately; a resumed run may use a new admitted placement.
    stable = dict(definition)
    stable.pop('commands')
    out.mkdir(parents=True, exist_ok=True)
    plan = out/'campaign_plan.json'
    if plan.exists():
        previous = json.loads(plan.read_text())
        previous.pop('commands', None)
        if previous != stable:
            raise ValueError('Existing plan sources/parameters differ; preserve it and choose a new campaign directory')
        if (out/'campaign_execution.json').exists() and not (args.resume or not args.execute):
            raise FileExistsError('Execution exists; use --resume or a fresh directory')
    write_json(plan, definition)
    write_json(out/'preparation_resources.json', resources)
    print('Prepared eight cases: Rs425/Rs608 x 15.5/21.5 uA x 775/1550 nm; eta=0.9, tau_ee=5 ps, escape=20 ps.')
    print(f'Workers: 26 Numba threads ({100*resources["worker_capacity_fraction"]:.2f}% worker capacity); one separate coordinator core; two whole cores reserved.')
    print(f'Stop: confirmed +0.1 mV crossing + 2 ps; no trigger by {args.ceiling_ps:g} ps is censored. No 1D recovery evolution.')
    if not args.execute:
        print(f'PREPARED ONLY. No simulation started. Plan: {plan}')
        return 0
    with (out/'.execution.lock').open('a+') as local_lock, (SCRATCH/'thermal_candidate_runtime.lock').open('a+') as global_lock:
        try:
            fcntl.flock(local_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            fcntl.flock(global_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError('Another campaign owns the shared runtime lock; existing jobs preserved')
        # Re-inventory immediately before launching all workers.
        resources = allocation()
        write_json(out/'execution_resources.json', resources)
        original_affinity = set(os.sched_getaffinity(0))
        os.sched_setaffinity(0, set(resources['coordinator_cpus']))
        ansi = sys.stdout.isatty() and os.environ.get('TERM') != 'dumb' and not args.plain
        board = Board(cases, out, args.ceiling_ps, ansi)
        stopped = threading.Event()
        processes = {}
        process_lock = threading.Lock()
        def request_stop(sig, frame):
            stopped.set()
            with process_lock:
                for process in processes.values():
                    if process.poll() is None:
                        process.send_signal(signal.SIGTERM)
        signal.signal(signal.SIGTERM, request_stop)
        signal.signal(signal.SIGINT, request_stop)
        def worker(case):
            name = case['name']
            dest = out/'runs'/name
            dest.parent.mkdir(exist_ok=True)
            if (dest/'manifest.json').exists():
                previous = json.loads((dest/'manifest.json').read_text())
                if previous.get('source_pins_sha256') != run.source_pins():
                    raise ValueError(f'{name}: saved producer source identity differs')
                if previous.get('run_complete'):
                    board.update(name, status=previous['status'], t=previous['time_ps'], exit_code=0)
                    return dict(case=name, status=previous['status'], run_complete=True, reused_complete=True)
                if not args.resume:
                    raise FileExistsError(f'{name}: existing result preserved; use --resume')
                if previous.get('status') != 'INTERRUPTED_CHECKPOINTED':
                    board.update(name, status=previous['status'], t=previous.get('time_ps'), exit_code=2,
                                 last_message='Diagnostic state preserved; not automatically resumed')
                    return dict(case=name, status=previous['status'], run_complete=False, diagnostic_preserved=True)
            if stopped.is_set():
                board.update(name, status='NOT_STARTED', exit_code=130)
                return dict(case=name, status='NOT_STARTED', run_complete=False)
            words = command(case, dest, resources['slots'][name], args, resume=(dest/'manifest.json').exists())
            env = dict(os.environ, NUMBA_NUM_THREADS=str(case['threads']),
                       PYTHONPYCACHEPREFIX=str(SCRATCH/'tmp/pysnspd/pycache'),
                       MPLCONFIGDIR=str(SCRATCH/'tmp/pysnspd/matplotlib'),
                       NUMBA_CACHE_DIR=str(SCRATCH/'tmp/pysnspd/numba_cache'))
            board.update(name, status='STARTING')
            logs = out/'logs'
            logs.mkdir(exist_ok=True)
            with (logs/(name+'.log')).open('a', buffering=1) as log:
                process = subprocess.Popen(['taskset', '-c', ','.join(map(str, sorted(original_affinity)))]+words,
                                           cwd=ROOT, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                           text=True, bufsize=1)
                with process_lock:
                    processes[name] = process
                    if stopped.is_set():
                        process.send_signal(signal.SIGTERM)
                board.update(name, pid=process.pid)
                for line in process.stdout:
                    log.write(line)
                    board.ingest(name, line)
                code = process.wait()
            saved = dest/'manifest.json'
            m = json.loads(saved.read_text()) if saved.exists() else {}
            status = m.get('status', 'FAILED_BEFORE_MANIFEST')
            if code != 0 and status in ('RUNNING', 'PREPARATION'):
                status = 'FAILED_EXIT'
            good = code == 0 and m.get('run_complete', False)
            board.update(name, status=status, exit_code=code, t=m.get('time_ps'),
                         steps=m.get('accepted_steps'), rejects=m.get('rejections'),
                         V=m.get('latest_monitor', {}).get('Vout_mV'),
                         latency=(m.get('observation', {}).get('crossing_time_s') or 0.)*1e12 or None,
                         cpu=None, rss=None,
                         last_message=m.get('error') or m.get('monitor_reason') or (status if code else ''))
            return dict(case=name, status=status, exit_code=code, time_ps=m.get('time_ps'),
                        run_complete=good, detected=m.get('observation', {}).get('detected'),
                        crossing_time_s=m.get('observation', {}).get('crossing_time_s'),
                        error=m.get('error'))
        results = []
        board.render()
        with ThreadPoolExecutor(max_workers=8) as pool:
            futures = {pool.submit(worker, case):case for case in cases}
            while futures:
                completed = [future for future in futures if future.done()]
                for future in completed:
                    case = futures.pop(future)
                    try:
                        result = future.result()
                    except Exception as exc:
                        result = dict(case=case['name'], status='COORDINATOR_WORKER_ERROR',
                                      run_complete=False, error=str(exc))
                        board.update(case['name'], status=result['status'], exit_code=2, last_message=str(exc))
                    results.append(result)
                    write_json(out/'campaign_execution.json', dict(results=results, finished=len(results), total=8))
                board.render()
                if futures:
                    # No long blocking wait: bounded refresh also handles quiet preparation phases.
                    time.sleep(args.refresh_seconds)
        board.render(final=True)
        good = all(result.get('run_complete') for result in results)
        print(f'Campaign {"COMPLETE" if good else "STOPPED WITH CHECKPOINTS OR DIAGNOSTICS"}: {len(results)}/8 cases finished.')
        for rs in (425, 608):
            for current in (15.5, 21.5):
                pair = [r for r in results if r['case'].startswith(f'rs{rs}_bias{str(current).replace(".", "p")}_')]
                times = {int(r['case'].rsplit('_', 1)[1]):r.get('crossing_time_s')
                         for r in pair if r.get('detected') and r.get('crossing_time_s') is not None}
                if len(times) == 2:
                    print(f'Rs={rs} ohm/sq, I={current:g} uA: relative latency (1550-775) = {(times[1550]-times[775])*1e12:.4f} ps')
                else:
                    print(f'Rs={rs} ohm/sq, I={current:g} uA: relative latency unavailable (incomplete/censored pair).')
        return 0 if good else 2


if __name__ == '__main__':
    raise SystemExit(main())
