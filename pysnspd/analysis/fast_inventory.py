"""Authenticate fast-mode results separately from normal native observations."""
from __future__ import annotations
from pathlib import Path
import copy,json,hashlib
import h5py,numpy as np
from pysnspd.experimental.fast_latency import PROTOCOL,ESTIMATED,COMPLETE_STATUSES,FAST_PARAMETER_KEYS,TailModel,HotbeltDetector,digest

def validate_run(run,source_root):
    run=Path(run);source_root=Path(source_root);raw=(run/'manifest.json').read_bytes();m=json.loads(raw)
    if m.get('execution_mode')!='fast' or m.get('fast_protocol')!=PROTOCOL or m.get('model_version')!='1.2.6':raise ValueError('Not an immutable v1.2.6 fast trajectory')
    if m.get('error') or not m.get('run_complete') or m.get('status') not in COMPLETE_STATUSES:raise ValueError('Fast observation incomplete or failed')
    p=m['parameter_identity']
    if p.get('execution_mode')!='fast' or p.get('fast_protocol')!=PROTOCOL:raise ValueError('Fast identity missing from native parameters')
    model=TailModel(source_root/'pysnspd/data/v126_fast_tail.json')
    required=set(model.data['baseline_source_pins'])|{'pipelines/53_v126_fast_photon.py','pipelines/18_spatial_photon_run.py','pipelines/35_spatial_energy_entry.py',
                                                   'pysnspd/experimental/fast_latency.py','pysnspd/data/v126_fast_tail.json'}
    if not required<=set(m['source_pins_sha256']):raise ValueError('Required native/fast source provenance missing')
    if m.get('fast_model_sha256')!=model.sha256 or p.get('fast_model_sha256')!=model.sha256:raise ValueError('Frozen model changed')
    for rel,h in m['source_pins_sha256'].items():
        path=(source_root/rel).resolve()
        if not path.is_relative_to(source_root.resolve()) or digest(path)!=h:raise ValueError('Fast source changed: '+rel)
    if digest(run/'checkpoint_latest.npz')!=m['checkpoint_sha256']:raise ValueError('Fast checkpoint changed')
    if digest(run/'dc1d.npz')!=m['dc1d_sha256']:raise ValueError('Fast reference changed')
    with np.load(run/'checkpoint_latest.npz',allow_pickle=False) as z:
        meta=json.loads(str(z['metadata_json'].item()))
        if meta['parameter_identity']!=p or meta['source_pins_sha256']!=m['source_pins_sha256'] or meta.get('terminal_diagnostic'):raise ValueError('Fast checkpoint identity/guard differs')
        if meta['controller'].get('fast_mode_state',{}).get('execution_mode')!='fast':raise ValueError('Fast checkpoint contains a normal controller')
        for k in ('delta_bar','p','n','circuit_SI'):
            if not np.isfinite(z[k]).all():raise ValueError('Nonfinite fast state')
    if not np.isfinite(m['max_global_energy_fraction']) or m['max_global_energy_fraction']>p['global_energy_tolerance']:raise ValueError('Fast global energy audit failed')
    if m['domain_fraction']>p['guard_cut_energy_relative']:raise ValueError('Fast domain guard failed')
    with h5py.File(run/'data.h5','r',swmr=True) as h:
        if h.attrs.get('execution_mode')!='fast':raise ValueError('Archive mode differs')
        q=h['2d'];n=m['archive_counts']['2d']['trace'];trace=q['trace'][:n]
        if len(trace)!=n or not len(trace) or not np.isfinite(trace['time_ps']).all() or np.any(np.diff(trace['time_ps'])<=0):raise ValueError('Invalid native fast trace')
        if trace['time_ps'][0]!=0 or abs(float(trace['time_ps'][-1])-m['time_ps'])>1e-8:raise ValueError('Fast trace does not span its own t0/checkpoint')
        estimated=m['status']==ESTIMATED
        if estimated:
            decision=m['fast_result']['decision'];hb=decision['HB_time_ps'];j=int(np.argmin(abs(trace['time_ps']-hb)))
            if not 0<j<len(trace)-1 or abs(trace['time_ps'][j]-hb)>1e-8:raise ValueError('Missing native HB derivative bracket')
            names=trace.dtype.names
            row=lambda i:{k:float(trace[k][i]) for k in names}
            recalculated=model.predict(p['current_uA'],p['deposited_energy_eV'],p['impact_y_nm'],row(j),row(j-1),row(j+1),ceiling_ps=p['ceiling_ps'],confirmation_ps=p['confirmation_ps'],stop_ps=m['time_ps'])
            if not recalculated.get('admitted') or abs(recalculated['estimated_latency_ps']-decision['estimated_latency_ps'])>1e-8:raise ValueError('Tail estimate cannot be independently reproduced')
            times=q['snapshot_time_ps'][:m['archive_counts']['2d']['snapshot_time_ps']];k=int(np.argmin(abs(times-hb)))
            if abs(times[k]-hb)>1e-8 or not HotbeltDetector(q['dc_delta_bar'][:],q['edges'][:],q['coordinates_m'][:])(q['delta_bar'][k]):raise ValueError('Missing independently observed HB map')
            tau=decision['estimated_latency_ps'];censored=False;bound=None
        else:
            from pysnspd.analysis.latency_catalog import replay_observation
            observation=replay_observation(trace,m['observation']);censored=observation['censored']
            tau=None if censored else observation['latency_s']*1e12;bound=observation['censored_lower_bound_s']*1e12 if censored else None
    if (run/'manifest.json').read_bytes()!=raw:raise ValueError('Fast manifest changed during extraction')
    return dict(case=run.name,run_dir=str(run),execution_mode='fast',current_uA=p['current_uA'],retained_energy_eV=p['deposited_energy_eV'],y0_nm=abs(p['impact_y_nm']),
                censored=censored,censored_lower_bound_ps=bound,tau_d_ps=tau,estimated=estimated,status=m['status'],manifest_sha256=hashlib.sha256(raw).hexdigest(),
                empirical_error_allowance_ps=model.allowance_ps if estimated else None)

def discover(database,source_root,reference):
    from pysnspd.analysis.latency_catalog import family_mismatches
    from pysnspd.analysis.simulation_inventory import key,live_workers,recent_speed
    accepted={};attempts={};errors=[];live=live_workers()
    for path in sorted(Path(database).glob('*/runs/*/manifest.json')):
        try:
            m=json.loads(path.read_text())
            if m.get('execution_mode')!='fast':continue
            p=m['parameter_identity'];clean=copy.deepcopy(m)
            for name in FAST_PARAMETER_KEYS:clean['parameter_identity'].pop(name,None)
            if family_mismatches(clean,reference):raise ValueError('Fast input belongs to another physical family')
            k=key(p['current_uA'],p['deposited_energy_eV'],p['impact_y_nm']);run=path.parent.resolve()
            item=dict(run_dir=str(run),status=m.get('status'),phase=m.get('phase'),time_ps=m.get('time_ps'),manifest_mtime=path.stat().st_mtime,
                      live=str(run) in live,execution_mode='fast')
            if item['live']:item['speed_ps_h']=recent_speed(run)
            attempts.setdefault(k,[]).append(item)
            if m.get('run_complete'):accepted.setdefault(k,[]).append(validate_run(run,source_root))
        except (OSError,ValueError,KeyError,TypeError) as exc:errors.append(dict(run_dir=str(path.parent),reason=str(exc)))
    return accepted,attempts,errors
