"""Authenticated fresh-directory continuation of completed v1.2.3 latency cases."""
import hashlib,json,os,shutil,subprocess
from pathlib import Path
import h5py
import numpy as np

SCRATCH=Path('/home/jdiaz/scratch')
NEW_SOURCES={'pipelines/18_spatial_photon_run.py',
    'pysnspd/experimental/reflected_photon.py','pysnspd/experimental/latency_continuation.py'}
SCHEDULING={'threads','physical_cores','latency_only','impact_y_nm','deposited_energy_eV','source_boundary'}

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(2**20),b''):h.update(block)
    return h.hexdigest()

def physical_identity(identity):return {k:v for k,v in identity.items() if k not in SCHEDULING}

def authenticate(source,pins,identity):
    source=Path(source).resolve();raw=(source/'manifest.json').read_bytes();manifest=json.loads(raw)
    expected={k:v for k,v in pins.items() if k not in NEW_SOURCES}
    if manifest['source_pins_sha256']!=expected:raise ValueError('Native evolution source identity changed; latency continuation refused')
    if physical_identity(manifest['parameter_identity'])!=physical_identity(identity):
        raise ValueError('Physical/numerical inputs differ from the completed base case')
    if (manifest.get('model_version')!='1.2.3' or not manifest.get('run_complete')
        or manifest.get('status')!='TRIGGER_MARGIN_COMPLETE' or manifest.get('phase')!='2d'
        or not manifest['parameter_identity'].get('latency_only')
        or not manifest.get('observation',{}).get('margin_complete')):
        raise ValueError('Only completed native v1.2.3 latency-only trigger+margin states are admitted')
    for proc in Path('/proc').iterdir():
        if not proc.name.isdigit() or int(proc.name)==os.getpid():continue
        try:
            words=(proc/'cmdline').read_bytes().decode(errors='replace').split('\0')
            if '--output-dir' in words and Path(words[words.index('--output-dir')+1]).resolve()==source:
                raise ValueError('Original base case is still active')
        except OSError:pass
    if sha(source/'checkpoint_latest.npz')!=manifest['checkpoint_sha256'] or sha(source/'dc1d.npz')!=manifest['dc1d_sha256']:
        raise ValueError('Original checkpoint/DC integrity mismatch')
    with np.load(source/'checkpoint_latest.npz',allow_pickle=False) as z:
        arrays={k:z[k].copy() for k in z.files}
    meta=json.loads(str(arrays['metadata_json']))
    if (meta['source_pins_sha256']!=expected or meta['parameter_identity']!=manifest['parameter_identity']
        or meta['phase']!='2d' or meta.get('terminal_diagnostic')
        or abs(meta['time_ps']-manifest['time_ps'])>1e-9):raise ValueError('Original manifest/checkpoint mismatch')
    counts=manifest['archive_counts']
    if counts['1d']['trace'] or counts['shadow1d']['trace']:raise ValueError('Expected a latency archive without prior 1D evolution')
    with h5py.File(source/'data.h5','r',swmr=True) as f:
        for phase,sizes in counts.items():
            for key,n in sizes.items():
                if len(f[phase][key])<n:raise ValueError('Archive shorter than committed counts')
        q=f['2d'];last=counts['2d']['snapshot_time_ps']-1
        if abs(float(q['snapshot_time_ps'][last])-meta['time_ps'])>1e-8:raise ValueError('Final saved 2D state absent')
        for key in ['delta_bar','p','n']:
            if not np.array_equal(q[key][last],arrays[key]):raise ValueError('Final archive/checkpoint state differs: '+key)
        trace=q['trace'][:counts['2d']['trace']];j=int(np.argmin(abs(trace['time_ps']-meta['time_ps'])))
        if abs(float(trace['time_ps'][j])-meta['time_ps'])>1e-8:raise ValueError('Final accepted circuit absent')
        if not np.array_equal(np.array([trace[k][j] for k in ('Ib_A','Is_A','vc_V')]),arrays['circuit_SI']):
            raise ValueError('Final accepted circuit/checkpoint mismatch')
        times=q['snapshot_time_ps'][:counts['2d']['snapshot_time_ps']]
        cross=manifest['observation']['crossing_time_s']*1e12
        begin=max(0,int(np.searchsorted(times,cross,side='right'))-1)
        if cross-times[begin]>.11 or len(times)-begin<3:raise ValueError('Saved overlap has insufficient time resolution')
    if (source/'manifest.json').read_bytes()!=raw:raise ValueError('Original changed during authentication')
    return manifest,meta,arrays

def copy_file(source,dest):
    # CoW when supported; ordinary independent copy otherwise. Never hard-link.
    subprocess.run(['cp','--reflink=auto','--',str(source),str(dest)],check=True)
    shutil.copystat(source,dest);return str(dest)

def prepare_latency_continuation(source,destination,pins,identity):
    source=Path(source).resolve();destination=Path(destination).resolve()
    if not source.is_relative_to(SCRATCH) or not destination.is_relative_to(SCRATCH):raise ValueError('Scratch paths required')
    if destination.exists():raise FileExistsError('Original/new data preserved; destination exists')
    manifest,meta,arrays=authenticate(source,pins,identity)
    print(json.dumps(dict(status='FORKING_SAVED_2D',time_ps=meta['time_ps'],source_run=str(source))),flush=True)
    archive_sha=sha(source/'data.h5')
    shutil.copytree(source,destination,copy_function=copy_file)
    if sha(destination/'data.h5')!=archive_sha:raise ValueError('Copied archive hash differs')
    (destination/'pre_update_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    shutil.copy2(source/'checkpoint_latest.npz',destination/'pre_update_checkpoint.npz')
    provenance=dict(schema='pysnspd.v124_latency_continuation.v1',source_run=str(source),
        source_checkpoint_sha256=manifest['checkpoint_sha256'],source_archive_sha256=archive_sha,
        source_pins_sha256=manifest['source_pins_sha256'],reused_until_ps=meta['time_ps'],reused_phase='2d',
        original_2d_physics_changed=False,physical_kernel_changed=False,photon_redeposited=False,
        overlap_policy='Replay native 1D only from the last saved 2D at/before trigger through trigger+2ps; compare exact simultaneous saved states',
        original_2d_trajectory_reintegrated=False,source_profile_preserved=True)
    meta.update(source_pins_sha256=pins,parameter_identity=identity,terminal_diagnostic=False,
        continuation=provenance,replay_saved_overlap=True,reaudit_saved_overlap=False,
        shadow_active=False,shadow_checks=[],shadow_error=None,projection=None)
    arrays={k:v for k,v in arrays.items() if not k.startswith('shadow_')}
    arrays['metadata_json']=np.array(json.dumps(meta,allow_nan=False))
    temp=destination/'checkpoint_recovery.tmp'
    with temp.open('wb') as f:np.savez_compressed(f,**arrays)
    os.replace(temp,destination/'checkpoint_latest.npz')
    manifest.update(model_version='1.2.4',run_complete=False,status='RECOVERY_PREPARED',error=None,
        source_pins_sha256=pins,parameter_identity=identity,continuation=provenance,
        checkpoint_sha256=sha(destination/'checkpoint_latest.npz'),
        observation_scope=dict(mode='authenticated full 1D continuation',one_d_evolution_requested=True,full_reset_requested=False))
    (destination/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(dict(status='RECOVERY_PREPARED',time_ps=meta['time_ps'],photon_redeposited=False)),flush=True)
