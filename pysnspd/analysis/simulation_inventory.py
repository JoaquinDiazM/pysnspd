"""Audited global coverage grid, distinct from a publication fit or PCR claim."""
from __future__ import annotations
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from zoneinfo import ZoneInfo
import hashlib, json, math, os, statistics

CURRENTS = (15.5,18.5,21.5)
POSITIONS = tuple(range(0,41,5))
EV_J = 1.602176634e-19
NONDETECTION_PS = 150.

def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024**2),b''): h.update(b)
    return h.hexdigest()

def finite(value):
    return isinstance(value,(int,float)) and math.isfinite(value)

def key(current,energy,y): return (float(current),round(float(energy),9),round(abs(float(y)),7))

def energy_of(m):
    p=m['parameter_identity']; source=m.get('source',{})
    if finite(source.get('retained_energy_J')): return source['retained_energy_J']/EV_J
    if finite(p.get('deposited_energy_eV')): return float(p['deposited_energy_eV'])
    return p['retained_fraction']*6.62607015e-34*299792458/(p['wavelength_nm']*1e-9)/EV_J

def read_json(path):
    raw=Path(path).read_bytes(); value=json.loads(raw)
    if Path(path).read_bytes()!=raw: raise ValueError('File changed during read')
    return value,hashlib.sha256(raw).hexdigest()

def live_workers():
    """Inspect only own processes; authenticate output paths from actual argv."""
    result={}
    for proc in Path('/proc').iterdir():
        if not proc.name.isdigit(): continue
        try:
            if proc.stat().st_uid!=os.getuid(): continue
            argv=proc.joinpath('cmdline').read_bytes().decode(errors='replace').split(chr(0))
            if not argv or not Path(argv[0]).name.startswith(('python','pypy')): continue
            if '--output-dir' not in argv or not any(Path(w).name in ('18_spatial_photon_run.py','35_spatial_energy_entry.py','13_dimensional_photon_run.py','53_v126_fast_photon.py') for w in argv): continue
            run=Path(argv[argv.index('--output-dir')+1])
            if not run.is_absolute(): run=proc.joinpath('cwd').resolve()/run
            stat=proc.joinpath('stat').read_text().rsplit(')',1)[1].split()
            result[str(run.resolve())]=dict(pid=int(proc.name),start_ticks=stat[19],affinity_cpus=sorted(os.sched_getaffinity(int(proc.name))))
        except (OSError,ValueError,IndexError): continue
    return result

def recent_speed(run):
    """Physical ps per worker wall-hour from recent native JSON progress lines."""
    log=run.parent.parent/'logs'/(run.name+'.log')
    try:
        with log.open('rb') as f:
            f.seek(max(0,log.stat().st_size-1024**2)); lines=f.read().decode(errors='replace').splitlines()
        rows=[]
        for line in lines:
            try:
                v=json.loads(line)
                if v.get('phase')=='2d' and finite(v.get('time_ps')) and finite(v.get('elapsed_session_wall_s')):
                    if rows and v['elapsed_session_wall_s']<rows[-1][0]: rows=[]
                    rows.append((v['elapsed_session_wall_s'],v['time_ps']))
            except (ValueError,TypeError): continue
        if len(rows)<2: return None
        tail=rows[-120:]; dt=tail[-1][0]-tail[0][0]; advance=tail[-1][1]-tail[0][1]
        return 3600*advance/dt if dt>=60 and advance>0 else None
    except OSError: return None

def resolve_cell(current,energy,y,accepted,candidates,planned):
    """Known detection wins; otherwise an active continuation remains red."""
    base=dict(current_uA=current,energy_eV=energy,display_energy_eV=f'{energy:.2f}',y_nm=y,
              tau_ps=None,elapsed_ps=None,value_ps=None,run_dir='',planned=planned,note='')
    if energy==0:
        return dict(base,color='yellow',state='analytic_zero',note='E=0: prior knowledge; no simulation; tau=infinity')
    detected=[r for r in accepted if not r['censored']]
    active=[r for r in candidates if r.get('live') and r['phase'] in ('2d',None,'PREPARATION')]
    if detected:
        r=min(detected,key=lambda v:v['tau_d_ps'])
        # Conflicting accepted detections are an audit error, never silently averaged.
        if max(v['tau_d_ps'] for v in detected)-r['tau_d_ps']>1e-6:
            return dict(base,color='red',state='conflict',note='Conflicting admitted delays')
        return dict(base,color='green',state='detected',tau_ps=r['tau_d_ps'],value_ps=r['tau_d_ps'],run_dir=r['run_dir'],note='Confirmed native output-threshold latency')
    if active:
        r=max(active,key=lambda v:v.get('time_ps') or 0.)
        return dict(base,color='red',state='running',elapsed_ps=r.get('time_ps'),value_ps=r.get('time_ps'),run_dir=r['run_dir'],note='Saved native time; unfinished; not a latency',speed_ps_h=r.get('speed_ps_h'))
    censored=[r for r in accepted if r['censored']]
    if censored:
        r=max(censored,key=lambda v:v['censored_lower_bound_ps'])
        t=r['censored_lower_bound_ps']; no_detection=t>NONDETECTION_PS
        return dict(base,color='yellow' if no_detection else 'gray',state='no_detection' if no_detection else 'censored',
                    elapsed_ps=t,value_ps=t,run_dir=r['run_dir'],note='Operational no-detection label (>150ps); not mathematical proof of infinite latency' if no_detection else 'Finite observation without confirmed detection')
    if candidates:
        r=max(candidates,key=lambda v:(v.get('time_ps') or 0.,v.get('manifest_mtime',0)))
        state='stale' if r['status']=='RUNNING' and not r.get('live') else 'paused'
        if r['status'] in ('INTEGRATION_FAILURE','NUMERICAL_RESOLUTION_LIMIT','FAILED'): state='failed'
        return dict(base,color='red',state=state,elapsed_ps=r.get('time_ps'),value_ps=r.get('time_ps'),run_dir=r['run_dir'],note='No admitted completed observation; status='+str(r['status']))
    excluded=current==15.5 and energy<=.2 and not planned
    return dict(base,color='red',state='excluded' if excluded else ('pending' if planned else 'unplanned'),note='Excluded from previous launches; not an assumed physical non-detection' if excluded else ('Prepared, not started' if planned else 'Grid input without a prepared campaign'))

def estimate(cells,positive_energies):
    """Transparent finite-horizon workload scenario, never an infinity ETA."""
    speeds=[r['speed_ps_h'] for r in cells if r['state']=='running' and finite(r.get('speed_ps_h')) and r['speed_ps_h']>0]
    red=[r for r in cells if r['energy_eV']>0 and r['color']=='red']
    cap_work=sum(max(0.,100.-(r.get('elapsed_ps') or 0.)) for r in red)
    result=dict(speed_samples=len(speeds),remaining_positive_inputs=len(red),remaining_ps_to_100ps=cap_work,
                scope='Scenario to100ps per unfinished input; future triggers may stop earlier; >150ps no-detection certification is separate and not scheduled; wall speed may vary')
    if not speeds: return dict(result,available=False)
    lo,hi=min(speeds),max(speeds);median=statistics.median(speeds)
    current_parallel=len(speeds)
    observed=[r for r in cells if r['energy_eV']>0 and r['color']!='red']
    def predicted_work(r):
        peers=[v for v in observed if key(v['current_uA'],v['energy_eV'],0)==key(r['current_uA'],r['energy_eV'],0)]
        peer=min(peers,key=lambda v:abs(v['y_nm']-r['y_nm'])) if peers else None
        stop=min(100.,peer['tau_ps']+2.) if peer and peer['state']=='detected' else 100.
        return max(0.,stop-(r.get('elapsed_ps') or 0.))
    neighbor_work=sum(predicted_work(r) for r in red)
    planned_work=sum(predicted_work(r) for r in red if r['planned'])
    return dict(result,available=True,observed_speed_min_ps_h=lo,observed_speed_median_ps_h=median,observed_speed_max_ps_h=hi,
                worker_hours_fast=cap_work/hi,worker_hours_slow=cap_work/lo,
                calendar_days_six_workers_fast=cap_work/hi/6/24,calendar_days_six_workers_slow=cap_work/lo/6/24,
                calendar_days_at_measured_workers_fast=cap_work/hi/current_parallel/24,
                calendar_days_at_measured_workers_slow=cap_work/lo/current_parallel/24,
                assumed_continuous_workers=6,nearest_remaining_physical_ps=neighbor_work,planned_remaining_physical_ps=planned_work,
                nearest_calendar_days_fast=neighbor_work/hi/6/24,nearest_calendar_days_slow=neighbor_work/lo/6/24,
                planned_calendar_days_fast=planned_work/hi/6/24,planned_calendar_days_slow=planned_work/lo/6/24,
                nearest_stop_scope='Same bias/exact retained energy; closest observed position, stop=tau+2ps or100ps if censored/unavailable. Heuristic only; no simulated or fitted delays added to cells')

def collect(database,source_root,reuse_catalog=None):
    from pysnspd.analysis.latency_catalog import catalog,family_mismatches,BASE_CAMPAIGN,BASE_NAMES
    started=datetime.now(timezone.utc); database=Path(database).resolve(); source_root=Path(source_root).resolve()
    prior=Path(reuse_catalog) if reuse_catalog else None
    if prior is None:
        pointer=Path('/home/jdiaz/scratch/latency_response/latest.json')
        if pointer.exists():
            latest,_=read_json(pointer);prior=Path(latest['artifact']).with_name('latency_catalog.csv')
    # Exact immutable catalog receipts are reused; new/changed cases undergo native replay.
    accepted,audit=catalog(database,source_root,reuse_catalog=prior,currents=CURRENTS)
    reference=json.loads((database/BASE_CAMPAIGN/'runs'/BASE_NAMES[0]/'manifest.json').read_text())
    workers=live_workers(); planned=set(); candidates=defaultdict(list); energies={r['retained_energy_eV'] for r in accepted}; observations=[]; excluded=[]; seen=set()
    # Discover manifests only; large field snapshots are not traversed or loaded.
    for path in sorted(database.glob('*/runs/*/manifest.json')):
        try:
            run=path.parent.resolve()
            if run in seen: continue
            seen.add(run);m,digest=read_json(path); p=m.get('parameter_identity',{})
            if family_mismatches(m,reference,currents=CURRENTS): continue
            if m.get('phase')=='0d' or (m.get('continuation') and not m.get('latency_observation_extension')): continue
            I=p['current_uA']; E=energy_of(m); y=abs(p.get('impact_y_nm',0.)); energies.add(E)
            item=dict(run_dir=str(run),current_uA=I,energy_eV=E,y_nm=y,status=m.get('status'),phase=m.get('phase'),
                      time_ps=m.get('time_ps') if m.get('phase')=='2d' else None,manifest_sha256=digest,
                      manifest_mtime=path.stat().st_mtime,live=str(run) in workers)
            if item['live'] and item['phase']=='2d': item['speed_ps_h']=recent_speed(run)
            candidates[key(I,E,y)].append(item); observations.append(item)
        except (OSError,ValueError,TypeError,KeyError) as exc: excluded.append(dict(path=str(path),reason=str(exc)))
    source_cache={}
    for path in sorted(database.glob('*/campaign_plan.json')):
        try:
            plan,_=read_json(path)
            # Only authenticated campaign provenance, or the original calibration plan.
            pins=plan.get('source_pins_sha256',{})
            if not pins or not all(pins.get(rel)==h for rel,h in reference['source_pins_sha256'].items()): continue
            for rel in pins:
                if rel not in source_cache: source_cache[rel]=sha256(source_root/rel)
            if not all(source_cache[rel]==h for rel,h in pins.items()): continue
            for c in plan.get('cases',[]):
                I=c.get('current_uA'); E=c.get('deposited_energy_eV'); y=c.get('impact_y_nm',0.)
                if I in CURRENTS and finite(E) and y in POSITIONS: planned.add(key(I,E,y))
        except (OSError,ValueError,TypeError,KeyError): continue
    # Map floating variants by tight physical identity; display rounding never changes keys.
    energy_groups={round(E,9):E for E in energies if E>0 and finite(E)}
    positive=sorted(energy_groups.values()); grouped=defaultdict(list)
    for r in accepted: grouped[key(r['current_uA'],r['retained_energy_eV'],r['y0_nm'])].append(r)
    # Fast estimates never enter the normal publication catalog/replay table.
    from pysnspd.analysis.fast_inventory import discover as discover_fast
    fast_accepted,fast_attempts,fast_errors=discover_fast(database,source_root,reference)
    cells=[]
    for I in CURRENTS:
        for E in [0.]+positive:
            for y in POSITIONS:
                k=key(I,E,y); cell=resolve_cell(I,E,y,grouped[k],candidates[k],k in planned or bool(grouped[k]) or bool(candidates[k]))
                cell['execution_mode']='normal'
                normal=dict(cell);modes={}
                if E>0 and (grouped[k] or candidates[k]):modes['normal']=normal
                if E>0 and (fast_accepted.get(k) or fast_attempts.get(k)):
                    rapid=resolve_cell(I,E,y,fast_accepted.get(k,[]),fast_attempts.get(k,[]),True)
                    rapid['execution_mode']='fast';rapid['estimated']=False
                    if rapid['state']=='detected':
                        selected=next(r for r in fast_accepted[k] if r['run_dir']==rapid['run_dir'])
                        rapid['estimated']=selected['estimated']
                        if selected['estimated']:
                            rapid['state']='estimated';rapid['empirical_error_allowance_ps']=selected['empirical_error_allowance_ps']
                    modes['fast']=rapid
                    if normal['color']=='red':cell=dict(rapid)
                cell['modes']=modes
                cell['both_modes_complete']=len(modes)==2 and all(r['color']!='red' for r in modes.values())
                cell['id']=len(cells)+1; cells.append(cell)
    pos=[r for r in cells if r['energy_eV']>0]; counts=dict(Counter(r['state'] for r in pos))
    counts.update(positive_grid=len(pos),analytic_zero_cells=sum(r['energy_eV']==0 for r in cells),total_cells=len(cells),
                  completed=sum(r['color']!='red' for r in pos),not_computed=sum(r['color']=='red' for r in pos),
                  planned_positive=sum(r['planned'] for r in pos),completion_fraction=sum(r['color']!='red' for r in pos)/len(pos),
                  normal_completed=sum(r.get('modes',{}).get('normal',{}).get('color') in ('green','gray','yellow') for r in pos),
                  fast_completed=sum(r.get('modes',{}).get('fast',{}).get('color') in ('green','gray','yellow') for r in pos),
                  fast_estimated=sum(r.get('modes',{}).get('fast',{}).get('estimated',False) for r in pos),
                  both_modes_completed=sum(r['both_modes_complete'] for r in pos))
    ended=datetime.now(timezone.utc)
    return dict(schema='pysnspd.simulation_inventory.v1',generated_utc=ended.isoformat(),generated_local=ended.astimezone(ZoneInfo('America/Santiago')).isoformat(),
                scan_started_utc=started.isoformat(),scan_finished_utc=ended.isoformat(),scan_wall_seconds=(ended-started).total_seconds(),
                currents_uA=list(CURRENTS),positive_energies_eV=positive,positions_nm=list(POSITIONS),non_detection_threshold_ps=NONDETECTION_PS,
                threshold_policy='Strictly more than150ps without confirmed detection; operational label only, not proof of mathematical divergence',
                display_energy_policy='Two decimals in reports; exact physical energy retained in JSON/CSV and source simulations',
                verification_scope=audit['verification_scope'],catalog_audit=audit,candidate_read_errors=excluded,observations=observations,live_workers=workers,
                cells=cells,counts=counts,fast_authentication_errors=fast_errors,eta=estimate(cells,positive),database=str(database),source_root=str(source_root))
