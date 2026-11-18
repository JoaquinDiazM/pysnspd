"""v1.2.6 delay-only observation over unchanged native evolution.

The run is fast from t=0, including cases whose uncertainty requires completing
the native observation window. No normal checkpoint can enter this protocol.
"""
from __future__ import annotations
import hashlib, json, math
from pathlib import Path
import numpy as np

PROTOCOL='pysnspd.fast_latency.v1'
ESTIMATED='FAST_DELAY_ESTIMATED'
COMPLETE_STATUSES=(ESTIMATED,'TRIGGER_MARGIN_COMPLETE','TIME_CENSORED_NO_TRIGGER')
FAST_PARAMETER_KEYS=('execution_mode','fast_protocol','fast_model_sha256')

def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024**2),b''):h.update(b)
    return h.hexdigest()

class TailModel:
    def __init__(self,path):
        self.path=Path(path);self.sha256=digest(path);self.data=json.loads(self.path.read_text())
        if self.data.get('schema')!='pysnspd.fast_tail.v1':raise ValueError('Unsupported frozen tail calibration')
        self.mean=np.array(self.data['mean'],float);self.scale=np.array(self.data['scale'],float);self.beta=np.array(self.data['beta'],float)
        if self.mean.shape!=(6,) or self.scale.shape!=(6,) or self.beta.shape!=(7,):raise ValueError('Invalid tail feature dimensions')
        if not all(np.isfinite(a).all() for a in (self.mean,self.scale,self.beta)) or np.any(self.scale<=0):raise ValueError('Invalid frozen coefficients')
        self.allowance_ps=float(self.data['empirical_error_allowance_ps'])
        if not 0<self.allowance_ps<=2:raise ValueError('Invalid calibrated error allowance')

    def predict(self,current,energy,y,hb,previous,following,*,ceiling_ps=100.,confirmation_ps=.5,stop_ps=None):
        """Central derivative uses exactly one accepted step after the HB sample."""
        if not previous['time_ps']<hb['time_ps']<following['time_ps']:return dict(admitted=False,reason='missing causal derivative bracket')
        slope=(following['Vout_V']-previous['Vout_V'])*1e3/(following['time_ps']-previous['time_ps'])
        raw=np.array([float(current==I) for I in (15.5,18.5,21.5)]+[hb['Vout_V']*1e3,hb['Vdev_V']*1e3,slope])
        if not np.isfinite(raw).all():return dict(admitted=False,reason='nonfinite electrical features')
        support=self.data['support'].get(str(float(current)))
        if support is None:return dict(admitted=False,reason='uncalibrated bias')
        if not support['energy_min_eV']-1e-9<=energy<=support['energy_max_eV']+1e-9 or not 0<=abs(y)<=40:
            return dict(admitted=False,reason='outside calibrated energy/position extent')
        lower=np.array(support['feature_min']);upper=np.array(support['feature_max']);pad=.05*np.maximum(upper-lower,1e-9)
        if np.any(raw[3:]<lower-pad) or np.any(raw[3:]>upper+pad):return dict(admitted=False,reason='outside calibrated HB electrical envelope')
        lag=float(np.r_[1.,(raw-self.mean)/self.scale]@self.beta);tau=hb['time_ps']+lag;eps=self.allowance_ps
        result=dict(admitted=False,reason='',model_sha256=self.sha256,HB_time_ps=float(hb['time_ps']),estimated_tail_ps=lag,
                    estimated_latency_ps=tau,empirical_error_allowance_ps=eps,features=raw.tolist(),
                    derivative_bracket_ps=[previous['time_ps'],following['time_ps']],
                    estimate_kind='calibrated estimate; not an integrated or confirmed trigger')
        if lag<=0 or tau-eps<=(following['time_ps'] if stop_ps is None else stop_ps):result['reason']='noncausal or too-near prediction'
        elif tau+eps+confirmation_ps>ceiling_ps:result['reason']='prediction uncertainty overlaps observation/confirmation horizon'
        else:result.update(admitted=True,reason='within frozen support and horizon allowance')
        return result

class HotbeltDetector:
    """The inherited <=0.2DC node-path criterion, with static adjacency reuse."""
    def __init__(self,dc,edges,xy):
        self.dc=np.abs(np.asarray(dc));xy=np.asarray(xy);self.adjacency=[[] for _ in self.dc]
        if np.any(self.dc<=0):raise ValueError('Nonzero DC gap required')
        for a,b in np.asarray(edges,int):self.adjacency[a].append(b);self.adjacency[b].append(a)
        y=xy[:,1];tol=1e-12*max(float(np.ptp(y)),1e-9)
        self.bottom=np.flatnonzero(abs(y-y.min())<=tol);self.top=np.flatnonzero(abs(y-y.max())<=tol)
    def __call__(self,delta):
        low=np.abs(delta)<=.2*self.dc;seeds=self.bottom[low[self.bottom]];targets=set(self.top[low[self.top]])
        if not len(seeds) or not targets:return False
        seen=set(map(int,seeds));todo=list(seen)
        while todo:
            i=todo.pop()
            if i in targets:return True
            for j in self.adjacency[i]:
                if low[j] and j not in seen:seen.add(j);todo.append(j)
        return False

class FastObserver:
    def __init__(self,model):
        self.model=model;self.args=None;self.controller=None;self.detector=None;self.state=None;self.rhs=None;self.reference=None
        self.previous=None;self.pending=None;self.next_sample_ps=0.;self.decision=None;self.saw_hb=False;self.monitor_stopped=False
    def bind(self,system,geometry,dc):
        self.detector=HotbeltDetector(dc.delta_bar,geometry.graph.edges,geometry.coordinates_m)
        self.reference=float(system.kinetic_residual_energy_scale_J)
    def process(self,row):
        # Both independent native ledger checks must pass at the stopping state.
        valid=(abs(row['energy_defect_J'])/(self.reference+row['absolute_external_J'])<=self.args.global_energy_tolerance
               and row['domain_fraction']<=self.args.guard_cut_energy_relative)
        if self.pending is not None and row['time_ps']>self.pending['hb']['time_ps']:
            p=self.pending;self.pending=None
            if valid:
                self.decision=self.model.predict(self.args.current_uA,self.args.deposited_energy_eV,abs(self.args.impact_y_nm),p['hb'],p['previous'],row,
                    ceiling_ps=self.args.ceiling_ps,confirmation_ps=self.args.confirmation_ps)
            else:self.decision=dict(admitted=False,reason='native independent energy/domain gate')
        if self.state is not None and not self.saw_hb and row['time_ps']>=self.next_sample_ps-1e-9:
            self.next_sample_ps=row['time_ps']+self.args.snapshot_every_ps
            if self.detector(self.state.delta_bar):
                self.saw_hb=True
                if self.previous is not None:self.pending=dict(hb=dict(row),previous=dict(self.previous))
                else:self.decision=dict(admitted=False,reason='initial HB lacks derivative history')
        self.previous=dict(row)
    def checkpoint(self):
        return dict(protocol=PROTOCOL,execution_mode='fast',model_sha256=self.model.sha256,previous=self.previous,pending=self.pending,
                    next_sample_ps=self.next_sample_ps,decision=self.decision,saw_hb=self.saw_hb)
    def restore(self,value):
        if value.get('protocol')!=PROTOCOL or value.get('execution_mode')!='fast' or value.get('model_sha256')!=self.model.sha256:
            raise ValueError('Checkpoint is not the identical fast protocol/calibration')
        for k in ('previous','pending','next_sample_ps','decision','saw_hb'):setattr(self,k,value[k])

def install(R,model,entry_relative):
    """Install observation/storage hooks; checked_step and every RHS stay intact."""
    observer=FastObserver(model);original_run=R.run;original_pins=R.source_pins;OriginalWindow=R.B.TriggerWindow
    original_metrics=R.monitor_metrics;original_monitor=R.monitor_observe;OriginalArchive=R.PhaseArchive
    raw_write=R.OLD._ORIGINAL_WRITE_JSON
    def pins():
        result=original_pins()
        for name in (entry_relative,'pysnspd/experimental/fast_latency.py','pysnspd/data/v126_fast_tail.json'):result[name]=digest(R.ROOT/name)
        return result
    def validate(a):
        if not a.latency_only or a.recover_from_run is not None or a.dc_only:raise ValueError('Fast mode requires a fresh latency-only photon, or its own resume')
        if not math.isclose(a.length_nm,300.,abs_tol=1e-12):raise ValueError('v1.2.6 fast calibration requires300nm')
        if a.deposited_energy_eV is None or math.isclose(a.deposited_energy_eV,.48,abs_tol=1e-9):raise ValueError('Explicit active retained energy required;0.48eV is isolated')
        allowed=set(model.data['allowed_parameter_variations'])|set(FAST_PARAMETER_KEYS)
        for k,v in model.data['reference_parameter_identity'].items():
            if k not in allowed and getattr(a,k,None)!=v:raise ValueError('Fast calibration physical/numerical input differs: '+k)
        if a.resume:
            saved=json.loads((Path(a.output_dir)/'manifest.json').read_text())
            if saved.get('execution_mode')!='fast' or saved.get('fast_protocol')!=PROTOCOL:raise ValueError('Normal-to-fast continuation is forbidden')
        a.execution_mode='fast';a.fast_protocol=PROTOCOL;a.fast_model_sha256=model.sha256
        observer.args=a
        return a
    class FastWindow(OriginalWindow):
        _terminal_statuses=OriginalWindow._terminal_statuses|{ESTIMATED}
        def __init__(self,config):super().__init__(config);observer.controller=self
        def observe_accepted(self,t,v):
            result=super().observe_accepted(t,v)
            if not self.terminal and self.crossing_time_s is None and self.candidate_crossing_s is None and observer.decision and observer.decision.get('admitted'):
                self.status=ESTIMATED;self.actual_stop_s=self.last_time_s
                return self.summary()
            return result
        def checkpoint(self):
            result=super().checkpoint();result['fast_mode_state']=observer.checkpoint();return result
        @classmethod
        def from_checkpoint(cls,value):
            if 'fast_mode_state' not in value:raise ValueError('Normal-to-fast controller continuation is forbidden')
            observer.restore(value['fast_mode_state']);return super().from_checkpoint(value)
    class SparseArchive(OriginalArchive):
        def __init__(self,path,systems,geometries,dcs,**kwargs):
            super().__init__(path,systems,geometries,dcs,**kwargs)
            self.fast_system=systems['2d']
            observer.bind(systems['2d'],geometries['2d'],dcs['2d'])
            self.file.attrs['execution_mode']='fast';self.file.attrs['fast_protocol']=PROTOCOL;self.file.attrs.modify('model_version','1.2.6')
            self.file.attrs['field_history']='initial and terminal only; resumable full state also in native checkpoints'
        def append(self,name,row):
            super().append(name,row)
            if name=='2d':
                observer.process(row)
                if observer.pending is not None and abs(observer.pending['hb']['time_ps']-row['time_ps'])<1e-12:
                    super().snapshot(name,self.fast_system,observer.state,observer.rhs,row['time_ps'])
        def snapshot(self,name,system,state,rhs,t):
            if t==0:
                observer.state=state
                if observer.detector(state.delta_bar):
                    observer.saw_hb=True;observer.decision=dict(admitted=False,reason='initial HB lacks derivative history')
                observer.next_sample_ps=observer.args.snapshot_every_ps
            at_hb=observer.pending is not None and abs(observer.pending['hb']['time_ps']-t)<1e-9
            if t==0 or at_hb or observer.monitor_stopped or (observer.controller is not None and observer.controller.terminal):
                super().snapshot(name,system,state,rhs,t)
    def metrics(system,geometry,state,rhs,*args,**kwargs):
        observer.state=state;observer.rhs=rhs
        return original_metrics(system,geometry,state,rhs,*args,**kwargs)
    def monitor(*args,**kwargs):
        result=original_monitor(*args,**kwargs)
        if result[0]:observer.monitor_stopped=True
        return result
    def write(path,value):
        if Path(path).name=='manifest.json':
            value.update(model_version='1.2.6',execution_mode='fast',fast_protocol=PROTOCOL,fast_model_sha256=model.sha256,
                         mode_locked_from_time_ps=0.,fast_origin='fresh photon or same-fast checkpoint; never a normal continuation')
            value['fast_result']=dict(decision=observer.decision,HB_observed=observer.saw_hb,estimated_trigger=value.get('status')==ESTIMATED,
                                      native_evolution_end_ps=value.get('time_ps'),field_history='sparse; full-state checkpoints retained')
            if value.get('status')==ESTIMATED and not value.get('error'):
                if not observer.decision or not observer.decision.get('admitted'):raise ValueError('Fast completion lacks admitted prediction')
                value['run_complete']=True
                value['scope']='Native300nm toHB plus one accepted derivative step; calibrated electrical-tail estimate, not native trigger confirmation'
        return raw_write(path,value)
    def run(a):
        validate(a)
        for rel,h in model.data['baseline_source_pins'].items():
            if pins().get(rel)!=h:raise ValueError('Fast calibration native source changed: '+rel)
        code=original_run(a)
        m=json.loads((Path(a.output_dir)/'manifest.json').read_text())
        return 0 if m.get('execution_mode')=='fast' and m.get('status') in COMPLETE_STATUSES and m.get('run_complete') and not m.get('error') else code
    R.source_pins=pins;R.B.TriggerWindow=FastWindow;R.PhaseArchive=SparseArchive;R.monitor_metrics=metrics;R.monitor_observe=monitor;R.write_json=write;R.run=run
    R.fast_validate=validate;R.fast_observer=observer
    return observer
