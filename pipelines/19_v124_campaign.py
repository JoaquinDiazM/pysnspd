"""Four spatial latency cases OR four authenticated baseline1D continuations.

Two independently checkpointed screens can coexist through atomic shared CPU
leases. Each campaign runs4x3 Numba threads, a separate coordination core,
BLAS1, and preserves at least two whole physical cores globally.
"""
from __future__ import annotations
import argparse,contextlib,fcntl,importlib.util,json,os,shutil,signal,subprocess,sys,threading,time
from concurrent.futures import ThreadPoolExecutor
from collections import deque
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
SCRATCH=Path('/home/jdiaz/scratch');PYTHON='/home/jdiaz/.conda/envs/snspd/bin/python'
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
os.environ.setdefault('NUMBA_NUM_THREADS','3')
def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
C=load('_v124_existing_progress_utilities',ROOT/'pipelines/15_calibration_basecases.py')

def catalog(mode):
    if mode=='recovery':return [dict(name=f'rs608_bias{str(i).replace(".","p")}_{w}',
        current_uA=i,wavelength_nm=w,impact_y_nm=0.,deposited_energy_eV=None)
        for i in (15.5,21.5) for w in (775,1550)]
    return [dict(name='edge35_I15p5_1550',current_uA=15.5,wavelength_nm=1550,impact_y_nm=35.,deposited_energy_eV=None),
        dict(name='edge35_I21p5_1550',current_uA=21.5,wavelength_nm=1550,impact_y_nm=35.,deposited_energy_eV=None),
        dict(name='center_I18p5_1550',current_uA=18.5,wavelength_nm=1550,impact_y_nm=0.,deposited_energy_eV=None),
        dict(name='center_I15p5_E096_775',current_uA=15.5,wavelength_nm=775,impact_y_nm=0.,deposited_energy_eV=.96)]

def scratch(path):
    value=Path(path).expanduser().resolve()
    if not value.is_relative_to(SCRATCH):raise ValueError('All campaign products must remain in scratch')
    return value

def pid_start(pid):
    try:return Path(f'/proc/{pid}/stat').read_text().rsplit(')',1)[1].split()[19]
    except OSError:return None

def active_workers(leases=()):
    records=[]
    for row in C.active_simulations():
        try:words=Path(f'/proc/{row["pid"]}/cmdline').read_bytes().decode(errors='replace').split('\0')
        except FileNotFoundError:continue
        # Peer19 coordinators may still have their inherited full affinity
        # before admission. Their complete reservation is protected by THIS
        # same atomic registry; count leases, not their transient startup mask.
        if any(Path(w).name=='19_v124_campaign.py' for w in words) and '--output-dir' not in words:continue
        # A newly exec'd18 worker inventories the full inherited topology before
        # its native runner pins the assigned mask. Its parent's live atomic
        # lease already reserves that complete allocation during this bootstrap.
        if any(Path(w).name=='18_spatial_photon_run.py' for w in words) and '--cpu-list' in words:
            try:
                stat=Path(f'/proc/{row["pid"]}/stat').read_text().rsplit(')',1)[1].split()
                parent=int(stat[1]);assigned=set(map(int,words[words.index('--cpu-list')+1].split(',')))
                owner=next((v for v in leases if v['pid']==parent and pid_start(parent)==v['start_ticks']),None)
                if owner and assigned and assigned<=set(owner['all_cpus']):continue
            except (OSError,ValueError,IndexError):pass
        records.append(row)
    return records

def allocation(resources,active,leases):
    budget=C.resource_budget(resources,.9,.9,1);allowed=set(budget['selected_affinity_cpus'])
    occupied={c for record in active for c in record['affinity_cpus']}
    occupied.update(c for record in leases for c in record['all_cpus'])
    if not occupied<=allowed:raise ValueError('Active jobs occupy CPUs outside the shared admitted budget; jobs preserved')
    free=[list(g) for g in resources['physical_core_groups'] if set(g)<=allowed and not set(g)&occupied]
    if len(free)<7:raise ValueError('Four3-thread cases need six FREE SMT2 worker cores plus one coordinator core; existing jobs preserved')
    chosen=free[:6];coordinator=free[6]
    if any(len(g)!=2 for g in chosen):raise ValueError('This layout requires SMT2 physical cores')
    slots=[]
    for i in (0,3):
        a,b,c=chosen[i:i+3];slots.extend([a+[b[0]],[b[1]]+c])
    all_cpus=sorted(set(coordinator)|{c for row in slots for c in row})
    if len(all_cpus)!=14 or set(all_cpus)&occupied:raise ValueError('CPU allocation overlaps')
    reserved=[list(g) for g in resources['physical_core_groups'] if not set(g)&(occupied|set(all_cpus))]
    if len(reserved)<2 or len(occupied|set(all_cpus))>budget['logical_budget']:raise ValueError('Global90%/two-whole-core reservation violated')
    memory=23*1024**3
    if memory+sum(v.get('memory_reservation_bytes',0) for v in leases)>budget['memory_limit_bytes']:
        raise ValueError('Insufficient conservatively reserved RAM for concurrent campaigns')
    return dict(slots=slots,coordinator_cpus=coordinator,all_cpus=all_cpus,
        reserved_whole_core_groups=reserved,worker_threads=12,memory_reservation_bytes=memory,
        external_processes=active,prior_leases=leases,combined_logical_cpus=sorted(occupied|set(all_cpus)),
        SMT_sharing_only_within_new_case_pairs=True,numerical_worker_masks_disjoint=True)

@contextlib.contextmanager
def reservation(campaign):
    lockpath=SCRATCH/'v124_resource_admission.lock';registry=SCRATCH/'v124_resource_leases.json'
    def read():
        records=json.loads(registry.read_text()) if registry.exists() else []
        return [v for v in records if pid_start(v['pid'])==v.get('start_ticks')]
    with lockpath.open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        records=read();resources=allocation(C.linux_resources(),active_workers(records),records)
        mine=dict(pid=os.getpid(),start_ticks=pid_start(os.getpid()),campaign_dir=str(campaign),
            all_cpus=resources['all_cpus'],memory_reservation_bytes=resources['memory_reservation_bytes'])
        C.write_json(registry,records+[mine]);os.sched_setaffinity(0,set(resources['coordinator_cpus']))
    try:yield resources
    finally:
        with lockpath.open('a+') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            C.write_json(registry,[v for v in read() if (v['pid'],v['start_ticks'])!=(mine['pid'],mine['start_ticks'])])

def command(case,dest,cpus,args,resume=False):
    cmd=[PYTHON,'-u',str(ROOT/'pipelines/18_spatial_photon_run.py'),'--source-runner',str(args.source_runner),
        '--output-dir',str(dest),'--threads','3','--cpu-list',','.join(map(str,cpus))]
    if args.mode=='recovery':
        cmd+=['--max-1d-wall-hours',str(args.max_1d_wall_hours),'--final-time-ns','15']
        if not resume:cmd+=['--recover-from-run',str(args.base_campaign/'runs'/case['name'])]
    else:
        cmd+=['--latency-only','--current-uA',str(case['current_uA']),'--wavelength-nm',str(case['wavelength_nm']),
            '--impact-y-nm',str(case['impact_y_nm']),'--sheet-resistance-ohm','608','--retained-fraction','0.9',
            '--tau-ee-Tc-ps','5','--tau-ep-Tc-ps','24.7','--tau-escape-ps','20','--bath-K','0.9',
            '--length-nm','300','--spacing-nm','5','--sigma-nm','7.549834435',
            '--threshold-mV','0.1','--confirmation-ps','0.5','--extra-time-ps','2','--ceiling-ps','100']
        if case['deposited_energy_eV'] is not None:cmd+=['--deposited-energy-eV',str(case['deposited_energy_eV'])]
    if resume:cmd.append('--resume')
    return cmd

class Board:
    def __init__(self,cases,out,mode,ansi,wall):
        self.out,self.mode,self.ansi,self.wall=out,mode,ansi,wall;self.lines=0;self.started=time.monotonic();self.lock=threading.Lock()
        self.rows={v['name']:dict(status='WAIT',phase='--',t=None,V=None,step=None,steps=0,rejects=0,
            pid=None,exit_code=None,cpu=None,rss=None,proc_sample=None,samples=deque(maxlen=8),target=100.,budget=None,error='') for v in cases}
    def update(self,name,**values):
        with self.lock:self.rows[name].update(values)
    def ingest(self,name,line):
        try:data=json.loads(line)
        except (ValueError,TypeError):return
        if not isinstance(data,dict):return
        with self.lock:
            r=self.rows[name];status=data.get('status',r['status']);phase=data.get('phase',r['phase'])
            if status in ('REPLAYING_1D_OVERLAP','1D_OVERLAP_REPLAY_COMPLETE'):phase='overlap'
            if status=='HANDOFF_ADMITTED':phase='1d'
            if phase!=r['phase']:r['samples'].clear();r['phase']=phase
            r['status']=status
            for k,v in [('t','time_ps'),('V','Vout_mV'),('step','step_ps'),('steps','accepted_steps'),('rejects','rejections')]:
                if data.get(v) is not None:r[k]=data[v]
            if data.get('time_ps') is not None:r['samples'].append((time.monotonic(),float(data['time_ps'])))
    def refresh(self):
        for name in self.rows:
            try:m=json.loads((self.out/'runs'/name/'manifest.json').read_text())
            except (FileNotFoundError,json.JSONDecodeError):m={}
            with self.lock:
                r=self.rows[name]
                if r['exit_code'] is None:
                    # Saved final2D timestamps must not overwrite the independent
                    # replay clock while it starts before that final timestamp.
                    if r['phase']!='overlap' or m.get('phase')=='1d':
                        phase=m.get('phase',r['phase'])
                        if phase!=r['phase']:r['samples'].clear();r['phase']=phase
                        r['status']=m.get('status',r['status'])
                        if m.get('time_ps') is not None and (r['t'] is None or m['time_ps']>=r['t']):r['t']=m['time_ps']
                    mon=m.get('latest_monitor') or {}
                    if r['phase']=='overlap' and m.get('phase')!='1d':r['V']=None
                    elif mon.get('Vout_mV') is not None:r['V']=mon['Vout_mV']
                    crossing=m.get('observation',{}).get('crossing_time_s')
                    r['target']=15000. if self.mode=='recovery' else 100. if crossing is None else crossing*1e12+2
                    r['budget']=max(0,self.wall*3600-m.get('runtime_1d_session_wall_s',0)) if self.mode=='recovery' else None
                    r['error']=m.get('error') or m.get('monitor_reason') or ''
                if r['pid'] is not None and r['exit_code'] is None:
                    try:
                        s=Path(f'/proc/{r["pid"]}/stat').read_text().rsplit(')',1)[1].split()
                        total=(int(s[11])+int(s[12]))/os.sysconf('SC_CLK_TCK');now=time.monotonic()
                        if r['proc_sample']:r['cpu']=100*(total-r['proc_sample'][1])/max(now-r['proc_sample'][0],.001)
                        r['proc_sample']=(now,total);r['rss']=int(s[21])*os.sysconf('SC_PAGE_SIZE')/1024**3
                    except (OSError,ValueError):pass
    def render(self):
        self.refresh();width=shutil.get_terminal_size((160,40)).columns if self.ansi else 160;compact=width<145
        with self.lock:
            text=[f'PYSNSPD1.2.4 |4 {self.mode} cases |3 threads/case |wall {C.hours(time.monotonic()-self.started)}',
                'Case                       State        Phase     t[ps] Vout[mV]  ps/h    ETA*  WallLeft' + ('' if compact else '  dt[ps]  accept/reject CPU% RAM'),'-'*min(width-1,158)]
            for name,r in self.rows.items():
                status={'TRIGGER_MARGIN_COMPLETE':'DETECTED','TIME_CENSORED_NO_TRIGGER':'CENSORED',
                    'WALL_BUDGET_CHECKPOINTED':'WALL-PAUSE','INTERRUPTED_CHECKPOINTED':'PAUSED',
                    'PREPARATION':'DC/source','RECOVERY_PREPARED':'forked','RUNNING':'RUN',
                    'REPLAYING_1D_OVERLAP':'replay1D','POSSIBLE_DIMENSIONAL_ERROR':'JOIN-STOP'}.get(r['status'],r['status'][:11])
                speed=None
                if len(r['samples'])>1:
                    a,b=r['samples'][0],r['samples'][-1];speed=3600*(b[1]-a[1])/max(b[0]-a[0],.001)
                eta=None if not speed or speed<=0 else 3600*max(0,r['target']-(r['t'] or 0))/speed
                fmt=lambda v,f:'--' if v is None else format(v,f)
                line=f'{name:<26} {status:<11} {r["phase"]:<7} {fmt(r["t"],"8.2f"):>8} {fmt(r["V"],"8.5f"):>8} {fmt(speed,"5.1f"):>5} {C.hours(eta):>7} {C.hours(r["budget"]):>8}'
                if not compact:line+=f' {fmt(r["step"],".1e"):>7} {str(r["steps"])+"/"+str(r["rejects"]):>13} {fmt(r["cpu"],"4.0f"):>4} {fmt(r["rss"],".1f"):>4}G'
                text.append(line)
            text+=['-'*min(width-1,158),
                '* ETA:15ns for recovery (guards may stop earlier); latency100ps cap then crossing+2ps. Not a detection forecast.',
                f'Recovery session budget:{self.wall:g}h native1D;12worker threads; two whole cores globally reserved.',str(self.out)]
            text += [name+': '+r['error'][-180:] for name,r in self.rows.items() if r.get('exit_code') not in (None,0) and r['error']]
            if self.ansi:
                text=[v[:max(width-1,20)] for v in text]
                if self.lines:sys.stdout.write(f'\x1b[{self.lines}A\x1b[J')
            sys.stdout.write('\n'.join(text)+'\n');sys.stdout.flush();self.lines=len(text)

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mode',choices=('spatial','recovery'),required=True)
    p.add_argument('--campaign-dir',type=Path,required=True)
    p.add_argument('--base-campaign',type=Path,default=SCRATCH/'photon_runs/calibration_eta090_ee5_rs425_rs608_20261021_01')
    p.add_argument('--source-runner',type=Path,required=True,help='Native source and marker for legacy CPU scanners')
    p.add_argument('--execute',action='store_true');p.add_argument('--prepare-only',action='store_true')
    p.add_argument('--resume',action='store_true');p.add_argument('--refresh-seconds',type=float,default=10.)
    p.add_argument('--max-1d-wall-hours',type=float,default=8.);p.add_argument('--plain',action='store_true')
    a=p.parse_args(argv)
    if a.execute and a.prepare_only:p.error('Choose execute OR prepare-only')
    if not 1<=a.refresh_seconds<=60 or not 0<a.max_1d_wall_hours<float('inf'):p.error('Invalid refresh/wall budget')
    out=scratch(a.campaign_dir);a.base_campaign=scratch(a.base_campaign)
    if out==a.base_campaign or out.is_relative_to(a.base_campaign):p.error('Use a separate fresh campaign directory')
    adapter=load('_v124_single_identity',ROOT/'pipelines/18_spatial_photon_run.py');R,_=adapter.configure(a.source_runner,3)
    producer_pins=R.source_pins();pins=dict(producer_pins);pins['pipelines/19_v124_campaign.py']=C.digest(__file__)
    pins['pipelines/15_calibration_basecases.py']=C.digest(ROOT/'pipelines/15_calibration_basecases.py')
    cases=catalog(a.mode)
    definition=dict(schema='pysnspd.v124_campaign.v1',model_version='1.2.4',mode=a.mode,cases=cases,
        base_campaign=str(a.base_campaign),source_pins_sha256=pins,threads_per_case=3,simultaneous_cases=4,
        max_1d_wall_hours=a.max_1d_wall_hours,final_time_ns=15.,coexisting_campaigns=True,
        source_policy='Reflecting transverse Gaussian images; one retained-energy inversion' if a.mode=='spatial' else 'Original source preserved; replay saved overlap then guarded1D')
    out.mkdir(parents=True,exist_ok=True);plan=out/'campaign_plan.json'
    if plan.exists() and json.loads(plan.read_text())!=definition:raise ValueError('Plan changed; preserve old campaign and use a new path')
    if (out/'campaign_execution.json').exists() and a.execute and not a.resume:raise ValueError('Existing execution: use --resume')
    C.write_json(plan,definition)
    if not a.execute:
        try:res=allocation(C.linux_resources(),active_workers(),[]);status='PREPARED'
        except ValueError as exc:res=dict(reason=str(exc));status='PREPARED_WAITING_CAPACITY'
        C.write_json(out/'preparation_resources.json',res)
        print(f'{status}:4 {a.mode} cases; no simulations started. '+res.get('reason','Runtime atomic leases are rechecked at execution.'))
        return 0
    with (out/'.execution.lock').open('a+') as local:
        try:fcntl.flock(local,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:raise RuntimeError('This campaign already runs; preserved')
        original_affinity=sorted(os.sched_getaffinity(0))
        with reservation(out) as resources:
            C.write_json(out/'execution_resources.json',resources)
            board=Board(cases,out,a.mode,sys.stdout.isatty() and os.environ.get('TERM')!='dumb' and not a.plain,a.max_1d_wall_hours)
            stopped=threading.Event();processes={};lock=threading.Lock()
            def stop(sig,frame):
                stopped.set()
                with lock:
                    for proc in processes.values():
                        if proc.poll() is None:proc.send_signal(signal.SIGTERM)
            signal.signal(signal.SIGINT,stop);signal.signal(signal.SIGTERM,stop)
            def worker(case,cpus):
                name=case['name'];dest=out/'runs'/name;dest.parent.mkdir(exist_ok=True)
                exists=(dest/'manifest.json').exists()
                if exists:
                    saved=json.loads((dest/'manifest.json').read_text())
                    if saved['source_pins_sha256']!=producer_pins:raise ValueError('Saved source identity changed')
                    if saved.get('run_complete'):
                        board.update(name,status=saved['status'],phase=saved.get('phase','--'),t=saved.get('time_ps'),exit_code=0)
                        return dict(case=name,status=saved['status'],run_complete=True,reused=True)
                    if not a.resume or saved.get('status') not in ('INTERRUPTED_CHECKPOINTED','WALL_BUDGET_CHECKPOINTED','DEVELOPMENT_DURATION_REACHED'):
                        raise ValueError('Existing diagnostic preserved; only an own checkpoint pause can resume')
                if stopped.is_set():return dict(case=name,status='NOT_STARTED',run_complete=False)
                cmd=command(case,dest,cpus,a,exists);logs=out/'logs';logs.mkdir(exist_ok=True)
                env=dict(os.environ,NUMBA_NUM_THREADS='3',PYTHONPYCACHEPREFIX=str(SCRATCH/'tmp/pysnspd/pycache'),
                    MPLCONFIGDIR=str(SCRATCH/'tmp/pysnspd/matplotlib'),NUMBA_CACHE_DIR=str(SCRATCH/'tmp/pysnspd/numba_cache'))
                board.update(name,status='STARTING')
                with (logs/(name+'.log')).open('a',buffering=1) as log:
                    proc=subprocess.Popen(['taskset','-c',','.join(map(str,original_affinity))]+cmd,cwd=ROOT,env=env,
                        stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,bufsize=1)
                    with lock:
                        processes[name]=proc
                        if stopped.is_set():proc.send_signal(signal.SIGTERM)
                    board.update(name,pid=proc.pid)
                    for line in proc.stdout:log.write(line);board.ingest(name,line)
                    code=proc.wait()
                path=dest/'manifest.json';m=json.loads(path.read_text()) if path.exists() else {}
                board.update(name,status=m.get('status','FAILED'),phase=m.get('phase','--'),t=m.get('time_ps'),
                    exit_code=code,cpu=None,rss=None,error=m.get('error') or m.get('monitor_reason') or '')
                return dict(case=name,status=m.get('status','FAILED'),exit_code=code,time_ps=m.get('time_ps'),
                    phase=m.get('phase'),run_complete=bool(m.get('run_complete')),error=m.get('error') or m.get('monitor_reason'))
            results=[];board.render()
            with ThreadPoolExecutor(max_workers=4) as pool:
                pending={pool.submit(worker,c,s):c for c,s in zip(cases,resources['slots'])}
                while pending:
                    for fut in [f for f in pending if f.done()]:
                        case=pending.pop(fut)
                        try:result=fut.result()
                        except Exception as exc:
                            result=dict(case=case['name'],status='WORKER_ERROR',error=str(exc),run_complete=False)
                            board.update(case['name'],status='WORKER_ERROR',exit_code=2,error=str(exc))
                        results.append(result);C.write_json(out/'campaign_execution.json',dict(results=results,total=4,finished=len(results)))
                    board.render()
                    if pending:time.sleep(a.refresh_seconds)
            paused=all(v.get('run_complete') or v['status'] in ('WALL_BUDGET_CHECKPOINTED','INTERRUPTED_CHECKPOINTED') for v in results)
            print('Campaign '+('COMPLETE_OR_CHECKPOINTED' if paused else 'STOPPED_WITH_DIAGNOSTICS')+f':{out}')
            return 0 if paused else 2

if __name__=='__main__':raise SystemExit(main())
