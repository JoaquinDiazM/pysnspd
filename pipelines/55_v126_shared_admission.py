"""Position-triplet admission compatible with source-frozen active v1/v2/v3 jobs.

Only scheduling metadata changes. Authenticated one-shot coordinators never
respawn completed tasks. Unknown owners remain fully reserved. Native sources,
worker affinity and numerical trajectories are never modified here.
"""
from __future__ import annotations
import contextlib,copy,fcntl,json,os,time
from pathlib import Path

LEGACY_SOURCES={
    '54_v126_fast_campaign.py':'e18cba4b2b8f5a1ae92a3be96ceae2063a6a0ba83c231fefb4f4ba92b11b47c1',
    '49_position_triplets.py':'f1cdea8723df6008e696ceec6c96f913b7b1d819086931747dcf81634911b6b3',
    '37_energy_expansion_batches4.py':'8d90f3e597718c4c1d05d9034e79fb2123a7c5aab4c72867781d567ac1c1002d',
    '39_energy_remaining_pairs.py':'3d2abf09b4b198b689a8dcb0f836db95a2800d1d3b262471bab0d1c7c491d397',
    '42_lowenergy_pairs.py':'9183dfbdc7bb832f41b03b3234b2770b94f508bff6f5cf38c6acab2ae720f124',
    '45_ultralow_energy_pairs.py':'d6ad6b0d58a9cf97f4403a2fce8fd1672c5cf00e47b71b18fd550e34ccb2c740'}
DRIVER='56_v126_fast_batches.py'
POLICY='pysnspd.shared_fast_admission.v5'
DRIVER_SPECS={
    '41_energy_pairs_shared.py':dict(coordinator='39_energy_remaining_pairs.py',
        policy='pysnspd.shared_pair_admission.v1',pins={
            'pipelines/40_shared_pair_admission.py':'e1697277947f1d3c9e68c54117b606412fe8afaa97ca5e6d10b2586b67825914',
            'pipelines/41_energy_pairs_shared.py':'4dfdca7d3e95bd3f2d940471702af982884521843e80c70a96cf7d36d096bd24'}),
    '44_lowenergy_pairs_shared.py':dict(coordinator='42_lowenergy_pairs.py',
        policy='pysnspd.shared_pair_admission.v2',pins={
            'pipelines/43_shared_energy_admission.py':'3eb06ba3cd3ac58c13a66b5c73f539834eb2ebfc2b39ccfb19bce0af88911b61',
            'pipelines/44_lowenergy_pairs_shared.py':'6d356a59c49cb78bb2c5e2f7ae3ca8e925814b9a9a9c007c86a0c0e045e2e546'}),
    '47_ultralow_energy_pairs_shared.py':dict(coordinator='45_ultralow_energy_pairs.py',
        policy='pysnspd.shared_pair_admission.v3',pins={
            'pipelines/46_shared_energy_admission.py':'b1bd074520c5383ca9a5fda0cbab8e293d4ffecd063f1e8631933d6952ee9920',
            'pipelines/47_ultralow_energy_pairs_shared.py':'e4d4b7e47c032996e73ad3053ff914dbe1322321cdb7625b7cd572b07d3a275c'}),
    '51_position_triplets_shared.py':dict(coordinator='49_position_triplets.py',policy='pysnspd.shared_position_admission.v4',pins={'pipelines/50_shared_position_admission.py':'211247fe77ffef26b21e1421e90d686492a3f63e4b8efc486f62be2d43cc71d5','pipelines/51_position_triplets_shared.py':'3d3a42663dbf54149244be2e53dd6fed9e3220221a5886f87cb7d3276868b506'}),
    DRIVER:dict(coordinator='54_v126_fast_campaign.py',policy=POLICY,pins=None)}

def process_words(pid):
    try:
        path=Path(f'/proc/{pid}')
        if path.stat().st_uid!=os.getuid():return []
        return (path/'cmdline').read_bytes().decode(errors='replace').split(chr(0))
    except OSError:return []

def normalize_leases(records,proofs,occupied,groups):
    """Pure proposal: retain active/future masks and RAM, transfer old margins."""
    rows=copy.deepcopy(records);audit=[]
    for row in rows:
        proof=proofs.get(row['pid'])
        if proof is None:continue
        before=set(row['all_cpus']);after=set(before)
        coord=set(proof['coordinator_cpus']);margin=set(row.get('external_capacity_cpus',[]))
        if not coord<=before or margin&coord or margin&occupied:
            raise ValueError('Owner coordinator/foreign margin overlaps live numerical work')
        if not margin<=before or any(set(g)&margin and not set(g)<=margin for g in groups):
            raise ValueError('Invalid foreign capacity provenance')
        released=[]
        for item in proof['finished']:
            cpus=set(item['cpus'])
            if len(cpus)!=3 or cpus&coord:raise ValueError('Invalid completed worker mask')
            if cpus<=after and not cpus&occupied:
                after-=cpus;released.append(item)
        after-=margin
        if occupied&(before-after):raise ValueError('Attempt to release live numerical CPU')
        if before==after:continue
        item=dict(pid=row['pid'],campaign_dir=row['campaign_dir'],before_cpus=sorted(before),
            after_cpus=sorted(after),completed_slots=released,released_foreign_margin=sorted(margin),
            proof='Authenticated no-respawn owner, terminal receipt plus manifest, no live native occupant')
        row['all_cpus']=sorted(after);row['external_capacity_cpus']=[]
        row.setdefault('shared_admission_history',[]).append(item)
        audit.append(item)
    return rows,audit

class Admission:
    def __init__(self,source,driver):
        self.S=source;self.root=source.ROOT;self.driver=Path(driver).resolve()
        self.scheduler_pins={str(Path(__file__).resolve().relative_to(self.root)):source.C.digest(__file__),
            str(self.driver.relative_to(self.root)):source.C.digest(self.driver)}

    def live(self):
        path=self.S.SCRATCH/'v124_resource_leases.json'
        rows=json.loads(path.read_text()) if path.exists() else []
        return [r for r in rows if self.S.pid_start(r['pid'])==r.get('start_ticks')]

    def process_path(self,word,pid=None):
        path=Path(word)
        if not path.is_absolute():
            base=Path.cwd() if pid is None else Path(f'/proc/{pid}/cwd').resolve(strict=True)
            path=base/path
        return path.resolve()

    def scheduler_identity(self,kind):
        spec=DRIVER_SPECS[kind]
        pins=self.scheduler_pins if kind==self.driver.name else spec['pins']
        return spec['policy'],pins

    def coordinator_kind(self,argv,pid=None):
        for word in argv:
            raw=Path(word)
            if raw.name not in LEGACY_SOURCES and raw.name not in DRIVER_SPECS:continue
            try:path=self.process_path(word,pid)
            except OSError:continue
            if path.name in LEGACY_SOURCES and path==(self.root/'pipelines'/path.name).resolve():
                if self.S.C.digest(path)==LEGACY_SOURCES[path.name]:return path.name
            if path.name in DRIVER_SPECS and path==(self.root/'pipelines'/path.name).resolve():
                policy,pins=self.scheduler_identity(path.name)
                if pins and all(self.S.C.digest(self.root/k)==v for k,v in pins.items()):return path.name
        return None

    def owner_proof(self,row):
        if self.S.pid_start(row['pid'])!=row.get('start_ticks'):raise ValueError('PID/start identity changed')
        argv=process_words(row['pid']);kind=self.coordinator_kind(argv,row['pid'])
        if kind is None:raise ValueError('Unknown owner preserved')
        camp=self.S.scratch(row['campaign_dir'])
        if '--campaign-dir' not in argv or self.process_path(argv[argv.index('--campaign-dir')+1],row['pid'])!=camp:
            raise ValueError('Owner command/campaign mismatch')
        plan=json.loads((camp/'campaign_plan.json').read_text())
        # New leases carry their original map atomically, closing the small
        # window before the coordinator writes execution_resources.json.
        res=row['admission_resource_map'] if kind in DRIVER_SPECS else json.loads((camp/'execution_resources.json').read_text())
        native_coordinator=DRIVER_SPECS[kind]['coordinator'] if kind in DRIVER_SPECS else kind
        if plan['source_pins_sha256'].get('pipelines/'+native_coordinator)!=LEGACY_SOURCES[native_coordinator]:
            raise ValueError('Owner source provenance mismatch')
        if kind in DRIVER_SPECS:
            policy,pins=self.scheduler_identity(kind)
            if (row.get('scheduler_policy')!=policy or
                row.get('scheduler_source_pins_sha256')!=pins or
                res.get('scheduler_source_pins_sha256')!=pins):
                raise ValueError('Shared-admission owner provenance mismatch')
        if '--batch' not in argv or argv[argv.index('--batch')+1]!=plan['batch']:
            raise ValueError('Owner batch identity mismatch')
        if (plan['mode']!='spatial' or plan['threads_per_case']!=3 or
            plan['simultaneous_cases']!=len(plan['cases']) or len(res['slots'])!=len(plan['cases']) or
            len(res['coordinator_cpus'])!=2):raise ValueError('Unsupported one-shot plan')
        coord=set(res['coordinator_cpus']);owned=set(coord)
        for mask in res['slots']:
            if len(mask)!=3 or len(set(mask))!=3 or owned&set(mask):raise ValueError('Overlapping original masks')
            owned.update(mask)
        owned.update(res.get('idle_worker_cpus',[]))
        margin=set(res.get('external_capacity_cpus',[]))
        if owned&margin or owned|margin!=set(res['all_cpus']):raise ValueError('Incomplete original allocation')
        if not set(row['all_cpus'])<=set(res['all_cpus']):raise ValueError('Lease exceeds original allocation')
        if not set(row.get('external_capacity_cpus',[]))<=margin:raise ValueError('Foreign margin changed')
        done_path=camp/'campaign_execution.json'
        done=json.loads(done_path.read_text()).get('results',[]) if done_path.exists() else []
        receipts={r['case']:r for r in done if r.get('run_complete') and
            (r.get('exit_code')==0 or r.get('reused') is True) and
            r.get('status') in ('TRIGGER_MARGIN_COMPLETE','TIME_CENSORED_NO_TRIGGER','FAST_DELAY_ESTIMATED')}
        finished=[]
        for case,mask in zip(plan['cases'],res['slots']):
            receipt=receipts.get(case['name'])
            if receipt is None:continue
            m=json.loads((camp/'runs'/case['name']/'manifest.json').read_text())
            if (not m.get('run_complete') or m.get('status')!=receipt['status'] or m.get('error') or
                m.get('observation',{}).get('status')!=receipt['status']):continue
            if kind==DRIVER and (m.get('execution_mode')!='fast' or not m.get('fast_protocol')):raise ValueError('Fast owner has a normal result')
            if any(m['parameter_identity'].get(k)!=case[k] for k in ('current_uA','deposited_energy_eV','impact_y_nm')):
                raise ValueError('Completed physical input mismatch')
            if any(plan['source_pins_sha256'].get(k)!=v for k,v in m['source_pins_sha256'].items()):
                raise ValueError('Completed native source mismatch')
            finished.append(dict(case=case['name'],cpus=list(mask)))
        return dict(coordinator_cpus=list(coord),finished=finished)

    def activity(self,records):
        """Exclude known startup coordinators; map protected bootstrap workers."""
        active=[]
        for row in self.S.C.active_simulations():
            if row['pid']==os.getpid():continue
            argv=process_words(row['pid'])
            if not argv or not Path(argv[0]).name.startswith(('python','pypy')):continue
            if '--output-dir' not in argv and self.coordinator_kind(argv,row['pid']) is not None:continue
            affinity=list(row['affinity_cpus'])
            if '--cpu-list' in argv and any(Path(w).name in ('18_spatial_photon_run.py','35_spatial_energy_entry.py','53_v126_fast_photon.py') for w in argv):
                try:
                    stat=Path(f'/proc/{row["pid"]}/stat').read_text().rsplit(')',1)[1].split()
                    parent=int(stat[1]);assigned=set(map(int,argv[argv.index('--cpu-list')+1].split(',')))
                    owner=next((r for r in records if r['pid']==parent),None)
                    if owner and assigned<=set(owner['all_cpus']) and assigned:affinity=sorted(assigned)
                except (OSError,ValueError,IndexError):pass
            active.append(dict(pid=row['pid'],affinity_cpus=affinity))
        return active

    def propose(self,count):
        before=self.live();resources=self.S.C.linux_resources();activity=self.activity(before)
        occupied={c for r in activity for c in r['affinity_cpus']};proofs={};preserved=[]
        for row in before:
            try:proofs[row['pid']]=self.owner_proof(row)
            except (OSError,ValueError,KeyError,IndexError,TypeError) as exc:
                preserved.append(dict(pid=row['pid'],reason=str(exc)))
        rows,audit=normalize_leases(before,proofs,occupied,resources['physical_core_groups'])
        # All understood old margins were removed in the pure proposal. Reserve
        # the complete current aggregate foreign load once, in the new lease.
        foreign=self.S.foreign_cpu_sample()
        result=self.S.allocation(resources,activity,rows,count,foreign)
        for row in rows:
            if self.S.pid_start(row['pid'])!=row.get('start_ticks'):
                raise ValueError('Owner changed during admission; retry')
        result.update(shared_admission_audit=audit,preserved_unknown_owners=preserved,
            scheduler_policy=POLICY,scheduler_source_pins_sha256=self.scheduler_pins,
            foreign_load_counted_once_for_authenticated_leases=True)
        return before,rows,result

    def preview(self,count):
        with (self.S.SCRATCH/'v124_resource_admission.lock').open('a+') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            return self.propose(count)[2]

    @contextlib.contextmanager
    def reservation(self,campaign,count):
        registry=self.S.SCRATCH/'v124_resource_leases.json';lockpath=self.S.SCRATCH/'v124_resource_admission.lock'
        with lockpath.open('a+') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            before,rows,result=self.propose(count)
            mine=dict(pid=os.getpid(),start_ticks=self.S.pid_start(os.getpid()),campaign_dir=str(campaign),
                all_cpus=result['all_cpus'],memory_reservation_bytes=result['memory_reservation_bytes'],
                external_capacity_cpus=result['external_capacity_cpus'],scheduler_policy=POLICY,
                scheduler_source_pins_sha256=self.scheduler_pins,
                admission_resource_map={k:result[k] for k in ('slots','coordinator_cpus','own_cpus','all_cpus',
                    'external_capacity_cpus','idle_worker_cpus','scheduler_source_pins_sha256')})
            self.S.C.write_json(registry,rows+[mine])
            try:os.sched_setaffinity(0,set(result['coordinator_cpus']))
            except Exception:
                self.S.C.write_json(registry,before);raise
        try:yield result
        finally:
            with lockpath.open('a+') as lock:
                fcntl.flock(lock,fcntl.LOCK_EX)
                self.S.C.write_json(registry,[r for r in self.live() if
                    (r['pid'],r['start_ticks'])!=(mine['pid'],mine['start_ticks'])])

    def cleanup_only(self,count,audit_dir):
        """Apply verified cleanup now, without launching a worker or reserving new CPUs."""
        registry=self.S.SCRATCH/'v124_resource_leases.json'
        with (self.S.SCRATCH/'v124_resource_admission.lock').open('a+') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            before,rows,result=self.propose(count)
            path=self.S.scratch(audit_dir);path.mkdir(parents=True,exist_ok=True)
            receipt=path/'cleanup_receipt.json'
            if receipt.exists():raise ValueError('Existing cleanup receipt preserved')
            self.S.C.write_json(receipt,dict(registry_before=before,registry_after=rows,proposed_admission=result,
                physical_processes_unchanged=True))
            self.S.C.write_json(registry,rows)
            return result
