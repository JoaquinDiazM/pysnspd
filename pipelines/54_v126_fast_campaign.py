"""v1.2.6 independent fast cases; fresh fast-mode origin and strict own resumes."""
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

CATALOG_PATH=SCRATCH/'v126_fast_red_20261116_01/batch_catalog.json'
CATALOG={};BATCH_SPECS={}

def set_catalog(path):
    global CATALOG_PATH,CATALOG,BATCH_SPECS
    CATALOG_PATH=scratch(path);CATALOG=json.loads(CATALOG_PATH.read_text())
    if CATALOG.get('schema')!='pysnspd.v126_fast_batches.v1':raise ValueError('Unrecognized immutable fast batch catalog')
    BATCH_SPECS={k:len(v) for k,v in CATALOG['batches'].items()}

def expansion_catalog():return [dict(c) for rows in CATALOG['batches'].values() for c in rows]
def catalog(mode,batch):
    if mode!='spatial' or batch not in BATCH_SPECS:raise ValueError('Unknown complete fast-mode batch')
    return [dict(c) for c in CATALOG['batches'][batch]]



def prepare_recovery_links(out,args,cases,producer_pins):
    """Reference exact existing 1D runs; never copy or reinitialize their state."""
    source_root=scratch(args.recovery_campaign)
    if out==source_root or out.is_relative_to(source_root):
        raise ValueError('Use a separate coordinator directory; original campaign preserved')
    rows=[]
    for case in cases:
        origin=source_root/'runs'/case['name'];path=origin/'manifest.json'
        m=json.loads(path.read_text())
        if not m.get('run_complete') and m.get('status') not in ('WALL_BUDGET_CHECKPOINTED','INTERRUPTED_CHECKPOINTED'):
            raise ValueError('Only paused native 1D checkpoints may continue')
        if m.get('phase')!='1d' or m['parameter_identity']['current_uA']!=15.5:
            raise ValueError('Exact low-current 1D continuation required')
        if m['source_pins_sha256']!=producer_pins:raise ValueError('Native source identity changed')
        for file,key in [('checkpoint_latest.npz','checkpoint_sha256'),('dc1d.npz','dc1d_sha256')]:
            if C.digest(origin/file)!=m[key]:raise ValueError('Native checkpoint/DC checksum mismatch')
        rows.append((case,origin.resolve()))
    (out/'runs').mkdir(exist_ok=True)
    for case,origin in rows:
        link=out/'runs'/case['name']
        if link.is_symlink():
            if link.resolve()!=origin:raise ValueError('Existing run link differs')
        elif link.exists():raise ValueError('Existing run directory preserved')
        else:link.symlink_to(origin,target_is_directory=True)
    return {case['name']:str(origin) for case,origin in rows}


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
        if any(Path(w).name in ('19_v124_campaign.py','20_low_bias_spatial_campaign.py','21_low_bias_followup.py','22_latency_surface_refinement.py','23_extended_recovery.py','24_midbias_latency_anchors.py','26_v125_recovery_pair.py','27_midbias_energy_position8.py','28_midbias_followup3.py','29_midbias_final6_shared_host.py','30_highbias_anchors_energy8.py','31_highbias_energy_position8.py','32_highbias_energy_position6_shared_host.py','33_highbias_final7.py','34_energy_expansion_centers6.py','36_energy_expansion_centers6_fixed.py','37_energy_expansion_batches4.py','39_energy_remaining_pairs.py') for w in words) and '--output-dir' not in words:continue
        # A newly exec'd18 worker inventories the full inherited topology before
        # its native runner pins the assigned mask. Its parent's live atomic
        # lease already reserves that complete allocation during this bootstrap.
        if any(Path(w).name in ('18_spatial_photon_run.py','35_spatial_energy_entry.py','53_v126_fast_photon.py') for w in words) and '--cpu-list' in words:
            try:
                stat=Path(f'/proc/{row["pid"]}/stat').read_text().rsplit(')',1)[1].split()
                parent=int(stat[1]);assigned=set(map(int,words[words.index('--cpu-list')+1].split(',')))
                owner=next((v for v in leases if v['pid']==parent and pid_start(parent)==v['start_ticks']),None)
                if owner and assigned and assigned<=set(owner['all_cpus']):continue
            except (OSError,ValueError,IndexError):pass
        records.append(row)
    return records

def foreign_cpu_sample(seconds=2.):
    """Read-only per-process/per-thread load sample; never inspect private argv.

    Foreign user compute tasks above25% CPU count toward a conservative whole
    core capacity budget. Sleeping BLAS pools are not confused with active work.
    Affinity is reported, not changed. Unpinned jobs can migrate after admission.
    """
    import math
    def take():
        rows={};stamp=time.monotonic()
        for p in Path('/proc').iterdir():
            if not p.name.isdigit():continue
            try:
                uid=p.stat().st_uid
                if uid==os.getuid() or uid<1000:continue
                st=(p/'stat').read_text().rsplit(')',1)[1].split()
                threads={}
                for task in (p/'task').iterdir():
                    try:
                        a=(task/'stat').read_text().rsplit(')',1)[1].split()
                        threads[int(task.name)]=(int(a[11])+int(a[12]),int(a[36]))
                    except OSError:continue
                rows[int(p.name)]=dict(uid=uid,start=st[19],ticks=int(st[11])+int(st[12]),
                    threads=threads,rss_bytes=int(st[21])*os.sysconf('SC_PAGE_SIZE'),
                    comm=(p/'comm').read_text().strip())
            except OSError:continue
        return stamp,rows
    first=take();time.sleep(seconds);last=take();dt=last[0]-first[0]
    hz=os.sysconf('SC_CLK_TCK');rows=[];cores=0;hot={}
    for pid,b in last[1].items():
        a=first[1].get(pid)
        if not a or a['start']!=b['start']:continue
        percent=100*(b['ticks']-a['ticks'])/(hz*dt)
        if percent<25.:continue
        active=[]
        for tid,(ticks,cpu) in b['threads'].items():
            before=a['threads'].get(tid,(0,cpu))[0]
            part=100*max(0,ticks-before)/(hz*dt)
            if part>=5.:
                active.append(dict(tid=tid,cpu_percent=part,last_cpu=cpu))
                hot[cpu]=hot.get(cpu,0.)+part
        # Ignore at most5% of one CPU of sampling jitter when rounding units;
        # the number of active threads is an independent conservative bound.
        reserve=max(len(active),int(math.ceil(percent/100.-.05)),1)
        cores+=reserve
        try:affinity=sorted(os.sched_getaffinity(pid))
        except OSError:affinity=[]
        rows.append(dict(pid=pid,uid=b['uid'],comm=b['comm'],start_ticks=b['start'],
            cpu_percent=percent,total_threads=len(b['threads']),active_threads=active,
            reserved_core_units=reserve,affinity_cpus=affinity,rss_bytes=b['rss_bytes']))
    return dict(sample_seconds=dt,processes=rows,reserved_core_units=cores,
        measured_cpu_units=sum(r['cpu_percent'] for r in rows)/100.,hot_cpus=hot,
        policy='Capacity margin based on sampled active work; foreign tasks remain untouched and may migrate')

def allocation(resources,active,leases,case_count,foreign):
    if case_count not in (1,2,3,4,6,7,8):raise ValueError('Case counts1/2/3/4/6/7/8 required')
    budget=C.resource_budget(resources,.9,.9,1);allowed=set(budget['selected_affinity_cpus'])
    groups=resources['physical_core_groups']
    occupied={c for r in active for c in r['affinity_cpus']}
    occupied.update(c for r in leases for c in r['all_cpus'])
    if not occupied<=set(resources['logical_cpus']):raise ValueError('Other jobs outside visible host affinity; preserved')
    free_groups=[list(g) for g in groups if not set(g)&occupied]
    def score(group):return sum(foreign['hot_cpus'].get(c,foreign['hot_cpus'].get(str(c),0.)) for c in group)
    foreign_groups=sorted(free_groups,key=lambda g:(-score(g),g[0]))[:foreign['reserved_core_units']]
    if len(foreign_groups)!=foreign['reserved_core_units']:
        raise ValueError('Foreign compute load exhausts available whole cores')
    external={c for g in foreign_groups for c in g}
    free=[list(g) for g in groups if set(g)<=allowed and not set(g)&(occupied|external)]
    worker_cores=(3*case_count+1)//2
    if len(free)<worker_cores+1:
        raise ValueError(f'{case_count} cases do not fit after reserving {len(foreign_groups)} external core units; use fewer cases or wait')
    chosen=free[:worker_cores];coordinator=free[worker_cores]
    if any(len(g)!=2 for g in chosen):raise ValueError('SMT2 topology required')
    # Preserve whole-core ownership for odd case counts. The final unused SMT
    # sibling stays leased; it must not become a peer campaign's worker slot.
    worker_pool=[cpu for group in chosen for cpu in group]
    slots=[worker_pool[i:i+3] for i in range(0,3*case_count,3)]
    idle_worker_cpus=worker_pool[3*case_count:]
    own=set(coordinator)|set(worker_pool)
    accounted=own|external
    reserved=[list(g) for g in groups if not set(g)&(occupied|accounted)]
    if (len(own)!=2*worker_cores+2 or own&(occupied|external) or len(reserved)<2 or
        len(occupied|accounted)>budget['logical_budget']):
        raise ValueError('Global90%/two-whole-core ceiling exceeded including foreign capacity margin')
    memory=(3+5*case_count)*1024**3
    if memory+sum(r.get('memory_reservation_bytes',0) for r in leases)>budget['memory_limit_bytes']:
        raise ValueError('Insufficient90% available-RAM reservation')
    return dict(slots=slots,idle_worker_cpus=idle_worker_cpus,coordinator_cpus=coordinator,own_cpus=sorted(own),
        all_cpus=sorted(accounted),external_capacity_cpus=sorted(external),
        external_capacity_core_groups=foreign_groups,foreign_load_sample=foreign,
        worker_threads=3*case_count,reserved_whole_core_groups=reserved,
        combined_logical_cpus=sorted(occupied|accounted),memory_reservation_bytes=memory,
        prior_leases=leases,external_processes=active,numerical_worker_masks_disjoint=True,
        foreign_affinity_unchanged=True,exclusive_foreign_core_isolation_guaranteed=False)

def admission(case_count):
    registry=SCRATCH/'v124_resource_leases.json'
    records=json.loads(registry.read_text()) if registry.exists() else []
    records=[r for r in records if pid_start(r['pid'])==r.get('start_ticks')]
    foreign=foreign_cpu_sample()
    return records,allocation(C.linux_resources(),active_workers(records),records,case_count,foreign)

def preview(case_count):
    with (SCRATCH/'v124_resource_admission.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        return admission(case_count)[1]

@contextlib.contextmanager
def reservation(campaign,case_count):
    lockpath=SCRATCH/'v124_resource_admission.lock';registry=SCRATCH/'v124_resource_leases.json'
    def read():
        records=json.loads(registry.read_text()) if registry.exists() else []
        return [r for r in records if pid_start(r['pid'])==r.get('start_ticks')]
    with lockpath.open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        records,resources=admission(case_count)
        mine=dict(pid=os.getpid(),start_ticks=pid_start(os.getpid()),campaign_dir=str(campaign),
            all_cpus=resources['all_cpus'],memory_reservation_bytes=resources['memory_reservation_bytes'],
            external_capacity_cpus=resources['external_capacity_cpus'])
        # The lease includes the foreign CAPACITY margin so peer coordinators
        # cannot consume it; none of the foreign processes is altered or owned.
        C.write_json(registry,records+[mine])
        try:os.sched_setaffinity(0,set(resources['coordinator_cpus']))
        except Exception:
            C.write_json(registry,records);raise
    try:yield resources
    finally:
        with lockpath.open('a+') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            C.write_json(registry,[r for r in read() if (r['pid'],r['start_ticks'])!=(mine['pid'],mine['start_ticks'])])

def case_source_pins(case,pins):
    # Low-energy trajectories continue through the exact original entrypoint.
    # New high-energy trajectories additionally authenticate the CLI extension.
    return dict(pins)


def command(case,dest,cpus,args,resume=False):
    cmd=[PYTHON,'-u',str(ROOT/'pipelines/53_v126_fast_photon.py'),'--source-runner',str(args.source_runner),
        '--output-dir',str(dest),'--threads','3','--cpu-list',','.join(map(str,cpus))]
    if args.mode=='recovery':
        cmd+=['--max-1d-wall-hours',str(args.max_1d_wall_hours),'--final-time-ns','15']
        if not resume:raise ValueError('Recovery only resumes authenticated existing native1D states')
    else:
        cmd+=['--latency-only','--current-uA',str(case['current_uA']),'--wavelength-nm',str(case['wavelength_nm']),
            '--impact-y-nm',str(case['impact_y_nm']),'--sheet-resistance-ohm','608','--retained-fraction','0.9',
            '--tau-ee-Tc-ps','5','--tau-ep-Tc-ps','24.7','--tau-escape-ps','20','--bath-K','0.9',
            '--length-nm','300','--spacing-nm','5','--sigma-nm','7.549834435',
            '--threshold-mV','0.1','--confirmation-ps','0.5','--extra-time-ps','2','--ceiling-ps','100']
        if case['deposited_energy_eV'] is not None:cmd+=['--deposited-energy-eV',str(case['deposited_energy_eV'])]
    if resume:cmd.append('--resume')
    if getattr(args,'native_duration_ps',None) is not None:cmd+=['--duration-ps',str(args.native_duration_ps)]
    return cmd

def physical_input_live(case):
    from pysnspd.analysis.simulation_inventory import live_workers,energy_of,key
    wanted=key(case['current_uA'],case['deposited_energy_eV'],case['impact_y_nm'])
    for run in live_workers():
        try:
            m=json.loads((Path(run)/'manifest.json').read_text());q=m['parameter_identity']
            if m.get('phase') in ('2d',None,'PREPARATION') and q.get('sheet_resistance_ohm')==608 and key(q['current_uA'],energy_of(m),q.get('impact_y_nm',0))==wanted:return True
        except (OSError,ValueError,KeyError):continue
    return False


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
            text=[f'PYSNSPD1.2.6 FAST |native1.2.4 |{len(self.rows)} {self.mode} cases |3 threads/case |wall {C.hours(time.monotonic()-self.started)}',
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
                '* ETA:100ps cap before trigger, then confirmed crossing+2ps. Not a detection forecast; no1D recovery.',
                f'{3*len(self.rows)}worker threads (3/case); confirmed +0.1mV for0.5ps, then crossing+2ps; two whole cores globally reserved.',str(self.out)]
            text += [name+': '+r['error'][-180:] for name,r in self.rows.items() if r.get('exit_code') not in (None,0) and r['error']]
            if self.ansi:
                text=[v[:max(width-1,20)] for v in text]
                if self.lines:sys.stdout.write(f'\x1b[{self.lines}A\x1b[J')
            sys.stdout.write('\n'.join(text)+'\n');sys.stdout.flush();self.lines=len(text)

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mode',choices=('spatial',),required=True)
    p.add_argument('--batch',choices=tuple(BATCH_SPECS),required=True)
    p.add_argument('--campaign-dir',type=Path,required=True)
    p.add_argument('--base-campaign',type=Path,default=SCRATCH/'photon_runs/calibration_eta090_ee5_rs425_rs608_20261021_01')
    p.add_argument('--source-runner',type=Path,required=True,help='Native source and marker for legacy CPU scanners')
    p.add_argument('--execute',action='store_true');p.add_argument('--prepare-only',action='store_true')
    p.add_argument('--resume',action='store_true');p.add_argument('--refresh-seconds',type=float,default=10.)
    p.add_argument('--native-duration-ps',type=float,help='Development-only native fast checkpoint probe; same-mode resume')
    p.add_argument('--max-1d-wall-hours',type=float,default=18.);p.add_argument('--plain',action='store_true')
    p.add_argument('--recovery-campaign',type=Path,default=SCRATCH/'photon_runs/v124_base_recovery_20261023_01')
    a=p.parse_args(argv)
    if a.execute and a.prepare_only:p.error('Choose execute OR prepare-only')
    if not 1<=a.refresh_seconds<=60 or not 0<a.max_1d_wall_hours<float('inf'):p.error('Invalid refresh/wall budget')
    out=scratch(a.campaign_dir);a.base_campaign=scratch(a.base_campaign)
    if out==a.base_campaign or out.is_relative_to(a.base_campaign):p.error('Use a separate fresh campaign directory')
    adapter=load('_v126_fast_identity',ROOT/'pipelines/53_v126_fast_photon.py');R,_=adapter.configure(a.source_runner,3)
    producer_pins=R.source_pins();pins=dict(producer_pins);pins['pipelines/54_v126_fast_campaign.py']=C.digest(__file__)
    pins['pipelines/15_calibration_basecases.py']=C.digest(ROOT/'pipelines/15_calibration_basecases.py')
    cases=catalog(a.mode,a.batch)
    if hasattr(adapter,'validate_arguments'):
        for case in cases:
            adapter.validate_arguments(command(case,out/'runs'/case['name'],[0,16,1],a,False)[3:])
    out.mkdir(parents=True,exist_ok=True)
    run_refs=prepare_recovery_links(out,a,cases,producer_pins) if a.mode=='recovery' else {}
    definition=dict(schema='pysnspd.v124_campaign.v1',model_version='1.2.6',native_producer_version='1.2.4',execution_mode='fast',batch_catalog_sha256=C.digest(CATALOG_PATH),mode=a.mode,batch=a.batch,cases=cases,
        base_campaign=str(a.base_campaign),source_pins_sha256=pins,threads_per_case=3,simultaneous_cases=len(cases),run_references=run_refs,
        max_1d_wall_hours=a.max_1d_wall_hours,final_time_ns=15.,coexisting_campaigns=True,
        sampling_scope='Every red physical input in the saved global diagnostic except live cases; full-fast fromt0, sparse history, frozen HB-tail estimate; 0.48eV isolated',
        source_policy='Reflecting transverse Gaussian images; one retained-energy inversion' if a.mode=='spatial' else 'Original source preserved; replay saved overlap then guarded1D')
    out.mkdir(parents=True,exist_ok=True);plan=out/'campaign_plan.json'
    if plan.exists() and json.loads(plan.read_text())!=definition:raise ValueError('Plan changed; preserve old campaign and use a new path')
    if (out/'campaign_execution.json').exists() and a.execute and not a.resume:raise ValueError('Existing execution: use --resume')
    C.write_json(plan,definition)
    if not a.execute:
        try:res=preview(len(cases));status='PREPARED'
        except ValueError as exc:res=dict(reason=str(exc));status='PREPARED_WAITING_CAPACITY'
        C.write_json(out/'preparation_resources.json',res)
        print(f'{status}:{len(cases)} {a.mode} cases; no simulations started. '+res.get('reason','Runtime atomic leases are rechecked at execution.'))
        return 0
    with (out/'.execution.lock').open('a+') as local:
        try:fcntl.flock(local,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:raise RuntimeError('This campaign already runs; preserved')
        original_affinity=sorted(os.sched_getaffinity(0))
        with contextlib.ExitStack() as shared_locks:
            if a.mode=='recovery':
                origin_lock=shared_locks.enter_context((scratch(a.recovery_campaign)/'.execution.lock').open('a+'))
                try:fcntl.flock(origin_lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
                except BlockingIOError:raise RuntimeError('Original recovery campaign already has a writer; preserved')
            resources=shared_locks.enter_context(reservation(out,len(cases)))
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
                    if saved.get('execution_mode')!='fast':raise ValueError('Normal/fast checkpoint mixing is forbidden')
                    if saved['source_pins_sha256']!=case_source_pins(case,producer_pins):raise ValueError('Saved source identity changed')
                    if saved.get('run_complete'):
                        board.update(name,status=saved['status'],phase=saved.get('phase','--'),t=saved.get('time_ps'),exit_code=0)
                        return dict(case=name,status=saved['status'],run_complete=True,reused=True)
                    if not a.resume or saved.get('status') not in ('INTERRUPTED_CHECKPOINTED','WALL_BUDGET_CHECKPOINTED','DEVELOPMENT_DURATION_REACHED'):
                        raise ValueError('Existing diagnostic preserved; only an own checkpoint pause can resume')
                if stopped.is_set():return dict(case=name,status='NOT_STARTED',run_complete=False)
                if physical_input_live(case):return dict(case=name,status='DEFERRED_NATIVE_ACTIVE',run_complete=False)
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
                if code and not m.get('error'):
                    tail=(logs/(name+'.log')).read_text(errors='replace').splitlines()
                    m['error']=next((line for line in reversed(tail) if line.strip()),'Worker exited before manifest')
                board.update(name,status=m.get('status','FAILED'),phase=m.get('phase','--'),t=m.get('time_ps'),
                    exit_code=code,cpu=None,rss=None,error=m.get('error') or m.get('monitor_reason') or '')
                return dict(case=name,status=m.get('status','FAILED'),exit_code=code,time_ps=m.get('time_ps'),
                    phase=m.get('phase'),run_complete=bool(m.get('run_complete')),error=m.get('error') or m.get('monitor_reason'))
            results=[];board.render()
            with ThreadPoolExecutor(max_workers=len(cases)) as pool:
                pending={pool.submit(worker,c,s):c for c,s in zip(cases,resources['slots'])}
                while pending:
                    for fut in [f for f in pending if f.done()]:
                        case=pending.pop(fut)
                        try:result=fut.result()
                        except Exception as exc:
                            result=dict(case=case['name'],status='WORKER_ERROR',error=str(exc),run_complete=False)
                            board.update(case['name'],status='WORKER_ERROR',exit_code=2,error=str(exc))
                        results.append(result);C.write_json(out/'campaign_execution.json',dict(results=results,total=len(cases),finished=len(results)))
                    board.render()
                    if pending:time.sleep(a.refresh_seconds)
            paused=all(v.get('run_complete') or v['status'] in ('WALL_BUDGET_CHECKPOINTED','INTERRUPTED_CHECKPOINTED','DEVELOPMENT_DURATION_REACHED','DEFERRED_NATIVE_ACTIVE') for v in results)
            print('Campaign '+('COMPLETE_OR_CHECKPOINTED' if paused else 'STOPPED_WITH_DIAGNOSTICS')+f':{out}')
            return 0 if paused else 2

if __name__=='__main__':
    if CATALOG_PATH.exists():set_catalog(CATALOG_PATH)
    raise SystemExit(main())
