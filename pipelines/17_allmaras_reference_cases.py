"""Four isolated historical-parameter controls alongside an identified live campaign.

Only prepare by default. --execute starts four cases if all resource checks pass.
The current campaign and its global exclusive lock are not altered or signaled.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
import contextlib
import fcntl
import importlib.util
import io
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import threading
import time

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
SCRATCH=Path('/home/jdiaz/scratch');PYTHON='/home/jdiaz/.conda/envs/snspd/bin/python'
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module
C=load('_allmaras_terminal_and_inventory',ROOT/'pipelines/15_calibration_basecases.py')
H=load('_allmaras_reference_definition',ROOT/'pipelines/16_allmaras_parameter_control.py')

def cases():
    return [dict(name=f'allmaras2020_bias{str(i).replace(".","p")}_{w}',
        current_uA=i,wavelength_nm=w,threads=3,sheet_resistance_ohm=608.)
        for i in (15.5,21.5) for w in (775,1550)]

def scratch(value):
    value=Path(value).expanduser().resolve()
    if not value.is_relative_to(SCRATCH.resolve()):raise ValueError('Campaign products must stay under scratch')
    return value

def inventory(coexist,resources=None,records=None):
    res=C.linux_resources() if resources is None else resources
    budget=C.resource_budget(res,.9,.9,1);allowed=set(budget['selected_affinity_cpus'])
    if records is None:
        records=[]
        for row in C.active_simulations():
            try:words=Path(f'/proc/{row["pid"]}/cmdline').read_bytes().decode(errors='replace').split('\0')
            except FileNotFoundError:continue
            records.append(dict(row,words=words))
    busy=set();old_coordinator=None
    for row in records:
        words=row['words'];names={Path(word).name for word in words}
        def option(flag):
            return words[words.index(flag)+1] if flag in words else None
        output=option('--output-dir');campaign=option('--campaign-dir')
        producer=('13_dimensional_photon_run.py' in names and output is not None
                  and Path(output).resolve().is_relative_to(coexist/'runs'))
        coordinator=('15_calibration_basecases.py' in names and campaign is not None
                     and Path(campaign).resolve()==coexist)
        if not producer and not coordinator:
            raise ValueError(f'Unidentified active simulation/coordinator PID{row["pid"]}; existing job preserved')
        busy.update(row['affinity_cpus'])
        if coordinator:old_coordinator=list(row['affinity_cpus'])
    if not busy<=allowed:raise ValueError('Existing campaign exceeds the currently available shared CPU budget')
    admitted=[list(g) for g in res['physical_core_groups'] if set(g)<=allowed]
    free=[g for g in admitted if not set(g)&busy]
    # Two I/O coordinators share the already admitted coordination core, never
    # an occupied numerical-worker core. If the old coordinator has exited, use
    # one free complete core instead.
    if old_coordinator:
        coord=next((g for g in admitted if set(old_coordinator)<=set(g)),None)
        if coord is None:raise ValueError('Existing coordinator does not occupy one complete admitted core')
    else:
        if not free:raise ValueError('No complete coordinator core free')
        coord=free.pop()
    if len(free)<6:
        raise ValueError('Need six free whole SMT2 worker cores for four simultaneous three-thread cases; retry after Rs608 finishes')
    selected=free[:6]
    if any(len(g)!=2 for g in selected):raise ValueError('This layout requires SMT2 core groups')
    slots={}
    for pair in range(2):
        a,b,c=selected[3*pair:3*pair+3]
        slots[cases()[2*pair]['name']]=a+[b[0]]
        slots[cases()[2*pair+1]['name']]=[b[1]]+c
    worker=[v for values in slots.values() for v in values]
    if len(worker)!=12 or len(set(worker))!=12 or set(worker)&busy:
        raise ValueError('New and existing worker CPU masks overlap')
    occupied=busy|set(worker)|set(coord)
    reserved=[list(g) for g in res['physical_core_groups'] if not set(g)&occupied]
    if len(occupied)>budget['logical_budget'] or len(reserved)<2:
        raise ValueError('Combined campaigns violate the global CPU ceiling/two-whole-core reserve')
    memory=23*1024**3
    if memory>budget['memory_limit_bytes']:raise ValueError('Insufficient available RAM for four new conservative worker reservations')
    return dict(inventory=res,budget=budget,external_processes=records,slots=slots,
        coordinator_cpus=coord,coordinator_shared_with_existing=bool(old_coordinator),
        reserved_whole_core_groups=reserved,worker_threads=12,existing_admitted_logical_cpus=sorted(busy),
        combined_admitted_logical_cpus=sorted(occupied),memory_reservation_bytes=memory,
        numerical_worker_masks_disjoint=True,scope='All four start together; SMT shared only among the new numerical cases; existing campaign untouched')

def command(case,out,cpus,args,resume=False):
    words=[PYTHON,'-u',str(ROOT/'pipelines/16_allmaras_parameter_control.py'),
        '--source-runner',str(args.source_runner),'--output-dir',str(out),
        '--threads','3','--cpu-list',','.join(map(str,cpus)),
        '--current-uA',str(case['current_uA']),'--wavelength-nm',str(case['wavelength_nm']),
        '--length-nm','300','--spacing-nm','5','--threshold-mV','0.1',
        '--confirmation-ps','0.5','--extra-time-ps','2','--ceiling-ps','100']
    if resume:words.append('--resume')
    return words

class Board(C.Board):
    def __init__(self,out,ansi):
        super().__init__(cases(),out,100.,ansi);self.previous_lines=0
    def render(self,final=False):
        self.lines=0;capture=io.StringIO()
        with contextlib.redirect_stdout(capture):super().render(final)
        text=capture.getvalue().replace('8 latency cases | eta=0.90, tau_ee=5 ps, escape=20 ps',
            '4 HISTORICAL controls | eta=.667, ee=6ps, gamma=60, fixedD=.5, Tb=2K')
        text=text.replace('26 Numba threads + coordinator; 2 whole cores reserved.',
            '12 new Numba threads; combined budget checked; 2 whole cores reserved.')
        if self.ansi and self.previous_lines:sys.stdout.write(f'\x1b[{self.previous_lines}A\x1b[J')
        sys.stdout.write(text);sys.stdout.flush();self.previous_lines=len(text.splitlines())

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--campaign-dir',type=Path,required=True)
    p.add_argument('--coexist-with-campaign',type=Path,required=True)
    p.add_argument('--source-runner',type=Path,required=True,
        help='Current native producer source; explicit argument also exposes this control to the unchanged existing resource scanner')
    mode=p.add_mutually_exclusive_group();mode.add_argument('--execute',action='store_true');mode.add_argument('--prepare-only',action='store_true')
    p.add_argument('--resume',action='store_true');p.add_argument('--refresh-seconds',type=float,default=10.)
    p.add_argument('--plain',action='store_true');a=p.parse_args(argv)
    a.source_runner=a.source_runner.resolve()
    if a.source_runner!=(ROOT/'pipelines/13_dimensional_photon_run.py').resolve():p.error('Current native producer required')
    if not 1<=a.refresh_seconds<=60:p.error('Refresh interval must be in[1,60]s')
    out=scratch(a.campaign_dir);coexist=scratch(a.coexist_with_campaign)
    if out==coexist or out.is_relative_to(coexist):p.error('Use a separate campaign directory')
    old=json.loads((coexist/'campaign_plan.json').read_text())
    if old.get('model_version')!='1.2.3':p.error('Expected the identified current v1.2.3 campaign')
    native=load('_reference_native_identity',a.source_runner)
    pins=native.source_pins()
    for path in [Path(__file__),ROOT/'pipelines/16_allmaras_parameter_control.py',ROOT/'pipelines/15_calibration_basecases.py']:
        pins[str(path.relative_to(ROOT))]=C.digest(path)
    definition=dict(schema='pysnspd.allmaras_historical_control.v1',model_version='1.2.3',
        historical_parameter_control=True,reference=H.REFERENCE,cases=cases(),
        coexist_with_campaign=str(coexist),source_pins_sha256=pins,
        horizon='confirmed positive0.1mV crossing+2ps, confirmation0.5ps; no trigger at100ps is censored',
        physical_scope='Current native2D+actual circuit;40nm transverse hotbelt; fixed D and independent gamma; no reimplementation of historical1D/no stochastic IRF fit',
        resources='four simultaneous3-thread cases; runtime combined admission with existing campaign and at least two whole cores reserved')
    out.mkdir(parents=True,exist_ok=True);plan=out/'campaign_plan.json'
    if plan.exists() and json.loads(plan.read_text())!=definition:raise ValueError('Existing historical plan changed; preserve it and use a new directory')
    if (out/'campaign_execution.json').exists() and a.execute and not a.resume:raise FileExistsError('Use --resume or a fresh campaign directory')
    C.write_json(plan,definition)
    if not a.execute:
        try:res=inventory(coexist);status='PREPARED'
        except ValueError as exc:res=dict(reason=str(exc));status='PREPARED_WAITING_CAPACITY'
        C.write_json(out/'preparation_resources.json',res)
        print(f'{status}: four historical controls, no simulations started. {res.get("reason", "Capacity available now; rechecked at execution.")}')
        return 0
    # Dedicated concurrent-campaign lease; never take, release or overwrite the
    # exclusive lock held by the unchanged main campaign. Explicit source-runner
    # argv markers make both new coordinator and children visible to its scanner.
    with (out/'.execution.lock').open('a+') as local,(SCRATCH/'concurrent_reference_resource.lock').open('a+') as lease:
        try:
            fcntl.flock(local,fcntl.LOCK_EX|fcntl.LOCK_NB);fcntl.flock(lease,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:raise RuntimeError('A reference campaign is already active; existing jobs preserved')
        resources=inventory(coexist);C.write_json(out/'execution_resources.json',resources)
        affinity=sorted(os.sched_getaffinity(0));os.sched_setaffinity(0,set(resources['coordinator_cpus']))
        print('Starting four historical controls alongside the identified campaign. Existing workers are not repinned or signaled.')
        board=Board(out,sys.stdout.isatty() and os.environ.get('TERM')!='dumb' and not a.plain)
        stop=threading.Event();processes={};lock=threading.Lock()
        def request(sig,frame):
            stop.set()
            with lock:
                for process in processes.values():
                    if process.poll() is None:process.send_signal(signal.SIGTERM)
        signal.signal(signal.SIGINT,request);signal.signal(signal.SIGTERM,request)
        def worker(case):
            name=case['name'];dest=out/'runs'/name;dest.parent.mkdir(exist_ok=True)
            if (dest/'manifest.json').exists():
                m=json.loads((dest/'manifest.json').read_text())
                expected={k:v for k,v in pins.items() if k not in
                    ('pipelines/17_allmaras_reference_cases.py','pipelines/15_calibration_basecases.py')}
                if m.get('source_pins_sha256')!=expected:raise ValueError(name+': saved source identity changed')
                if m.get('run_complete'):
                    board.update(name,status=m['status'],t=m['time_ps'],exit_code=0)
                    return dict(case=name,current_uA=case['current_uA'],wavelength_nm=case['wavelength_nm'],
                        run_complete=True,status=m['status'],detected=m['observation'].get('detected'),
                        crossing_time_s=m['observation'].get('crossing_time_s'))
                if not a.resume or m.get('status')!='INTERRUPTED_CHECKPOINTED':
                    raise ValueError(name+': diagnostic/existing state preserved; only an interrupted own checkpoint can resume')
            if stop.is_set():return dict(case=name,run_complete=False,status='NOT_STARTED')
            cmd=command(case,dest,resources['slots'][name],a,(dest/'manifest.json').exists())
            env=dict(os.environ,NUMBA_NUM_THREADS='3',PYTHONPYCACHEPREFIX=str(SCRATCH/'tmp/pysnspd/pycache'),
                MPLCONFIGDIR=str(SCRATCH/'tmp/pysnspd/matplotlib'),NUMBA_CACHE_DIR=str(SCRATCH/'tmp/pysnspd/numba_cache'))
            logs=out/'logs';logs.mkdir(exist_ok=True);board.update(name,status='STARTING')
            with (logs/(name+'.log')).open('a',buffering=1) as log:
                proc=subprocess.Popen(['taskset','-c',','.join(map(str,affinity))]+cmd,cwd=ROOT,env=env,
                    stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,bufsize=1)
                with lock:
                    processes[name]=proc
                    if stop.is_set():proc.send_signal(signal.SIGTERM)
                board.update(name,pid=proc.pid)
                for line in proc.stdout:log.write(line);board.ingest(name,line)
                code=proc.wait()
            m=json.loads((dest/'manifest.json').read_text()) if (dest/'manifest.json').exists() else {}
            status=m.get('status','FAILED_BEFORE_MANIFEST');good=code==0 and m.get('run_complete',False)
            board.update(name,status=status,exit_code=code,t=m.get('time_ps'),cpu=None,rss=None,
                V=m.get('latest_monitor',{}).get('Vout_mV'),last_message=m.get('error') or (status if code else ''))
            return dict(case=name,current_uA=case['current_uA'],wavelength_nm=case['wavelength_nm'],
                status=status,run_complete=good,exit_code=code,detected=m.get('observation',{}).get('detected'),
                crossing_time_s=m.get('observation',{}).get('crossing_time_s'),error=m.get('error'))
        results=[];board.render()
        with ThreadPoolExecutor(max_workers=4) as pool:
            pending={pool.submit(worker,c):c for c in cases()}
            while pending:
                for future in [v for v in pending if v.done()]:
                    c=pending.pop(future)
                    try:result=future.result()
                    except Exception as exc:
                        result=dict(case=c['name'],run_complete=False,status='WORKER_ERROR',error=str(exc))
                        board.update(c['name'],status='WORKER_ERROR',exit_code=2,last_message=str(exc))
                    results.append(result);C.write_json(out/'campaign_execution.json',dict(results=results,finished=len(results),total=4))
                board.render()
                if pending:time.sleep(a.refresh_seconds)
        for current in (15.5,21.5):
            pair={v['wavelength_nm']:v.get('crossing_time_s') for v in results
                if v.get('current_uA')==current and v.get('detected') and v.get('crossing_time_s') is not None}
            if len(pair)==2:print(f'Historical control I={current:g}uA:1550-minus775={(pair[1550]-pair[775])*1e12:.4f}ps')
            else:print(f'Historical control I={current:g}uA: relative latency unavailable; incomplete/censored pair')
        good=all(v.get('run_complete') for v in results)
        print('Reference campaign '+('COMPLETE' if good else 'STOPPED WITH CHECKPOINTS OR DIAGNOSTICS'))
        return 0 if good else 2

if __name__=='__main__':raise SystemExit(main())
