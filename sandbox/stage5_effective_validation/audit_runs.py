"""External saved-data audit for effective BCS/KWT photon runs.

Interpolation is confined to comparing saved observables. Electrical detection
and latency are read only from the accepted-step trigger controller.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,os,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
 os.environ[key]="1"
import numpy as np
from scipy.constants import elementary_charge as E_CHARGE
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
SCRATCH=Path("/home/jdiaz/scratch")
MEV_J=1e-3*E_CHARGE
KINDS={"temporal":.02,"mesh":.05,"domain":.05,"spectral":.02}
RUNTIME_KEYS={"output_dir","resume","cpu","cpu_start","threads","snapshot_every_ps",
 "checkpoint_every_ps","checkpoint_wall_s","duration_ps","max_iterations","max_rejections"}
REFINEMENT_KEYS={"temporal":{"dt_ps"},"mesh":{"spacing_nm"},"domain":{"length_nm"},
 "spectral":{"xi_cutoff","xi_order","phonon_nodes","transport_order"}}


TEMPORAL_ENERGY_GUARD_CONTROL=dict(reference_dt_cap_ps=.02,candidate_dt_cap_ps=.01,
 reference_local_energy_tolerance=1e-5,candidate_local_energy_tolerance=1.25e-6,
 declaration_utc="2026-10-11T15:12:53.580984+00:00",
 rationale="Local midpoint energy defect is approximately cubic in step;8x stricter guard targets2x smaller accepted critical steps. This is a declared temporal numerical control, not a physical change.")
def declared_temporal_guard_control(kind,ai,bi):
 if kind!="temporal" or ai.get("local_energy_tolerance")==bi.get("local_energy_tolerance"):return None
 c=TEMPORAL_ENERGY_GUARD_CONTROL
 pairs=((ai.get("dt_ps"),c["reference_dt_cap_ps"]),(bi.get("dt_ps"),c["candidate_dt_cap_ps"]),
  (ai.get("local_energy_tolerance"),c["reference_local_energy_tolerance"]),
  (bi.get("local_energy_tolerance"),c["candidate_local_energy_tolerance"]),
  (ai.get("step_tolerance"),2e-7),(bi.get("step_tolerance"),2e-7))
 if any(actual!=expected for actual,expected in pairs):return None
 return dict(c)
def accepted_step_statistics(run,start_s,end_s):
 history=run["history"]
 if "accepted_step_ps" not in history:return dict(available=False,reason="Accepted-step history unavailable")
 t=np.asarray(history["time_s"]);h=np.asarray(history["accepted_step_ps"])
 clock_roundoff=32*np.finfo(float).eps*max(abs(start_s),abs(end_s),1e-12)
 mask=(t>start_s+clock_roundoff)&(t<=end_s+clock_roundoff)&(h>0)
 values=h[mask]
 if not len(values):return dict(available=False,reason="No positive accepted steps in common interval")
 if not np.isfinite(values).all():raise ValueError("Nonfinite saved accepted steps")
 return dict(available=True,count=len(values),minimum_ps=float(values.min()),
  median_ps=float(np.median(values)),maximum_ps=float(values.max()),
  scope="Steps ending in common saved trajectory interval; actual accepted h, diagnostic not nominal cap")


def guarded_common_step_statistics(reference,candidate,start_s,end_s):
 anchors=[]
 for run in (reference,candidate):
  guard=run["summary"].get("local_energy_guard",{})
  if not guard.get("applicable") or "guard_active_since_time_ps" not in guard:
   return dict(available=False,reason="Both saved local-guard anchors are required")
  anchors.append(float(guard["guard_active_since_time_ps"])*1e-12)
 start=max(start_s,*anchors)
 clock_roundoff=32*np.finfo(float).eps*max(abs(start),abs(end_s),1e-12)
 if end_s-start<=clock_roundoff:
  return dict(available=False,reason="No nonzero common post-anchor accepted interval within binary64 clock roundoff",
   maximum_guard_anchor_ps=max(anchors)*1e12)
 return dict(available=True,interval_ps=[start*1e12,end_s*1e12],
  reference=accepted_step_statistics(reference,start,end_s),
  candidate=accepted_step_statistics(candidate,start,end_s),
  scope="Common interval after both guard anchors; no automatic2x ratio of every adaptive step is inferred")

def sha256(path):
 h=hashlib.sha256()
 with Path(path).open("rb") as stream:
  for block in iter(lambda:stream.read(1024*1024),b""):h.update(block)
 return h.hexdigest()
def json_ready(value):
 if isinstance(value,np.ndarray):return json_ready(value.tolist())
 if isinstance(value,np.generic):return json_ready(value.item())
 if isinstance(value,Path):return str(value)
 if isinstance(value,float) and not math.isfinite(value):return None
 if isinstance(value,dict):return {str(k):json_ready(v) for k,v in value.items()}
 if isinstance(value,(list,tuple)):return [json_ready(v) for v in value]
 return value
def load_npz(path):
 with np.load(path,allow_pickle=False) as saved:return {key:saved[key].copy() for key in saved.files}
def finite(name,value):
 a=np.asarray(value)
 if np.any(~np.isfinite(a)):raise ValueError(name+" contains nonfinite values")
 return a
def increasing_time(value,name):
 t=finite(name,value).astype(float)
 if t.ndim!=1 or len(t)<1 or np.any(np.diff(t)<=0):raise ValueError(name+" must increase strictly")
 return t
def require_series(data,key,t):
 value=finite(key,data[key]).astype(float)
 if value.shape!=t.shape:raise ValueError(key+" has inconsistent trajectory shape")
 return value
def checkpoint_support(path):
 arrays=load_npz(path);p=finite("checkpoint.p",arrays["p"]);n=finite("checkpoint.n",arrays["n"])
 delta=finite("checkpoint.delta_bar",arrays["delta_bar"]);cm=finite("checkpoint.circuit_SI",arrays["circuit_SI"])
 meta=json.loads(str(arrays["metadata_json"].item()))
 return dict(p_min=float(p.min()),p_max=float(p.max()),n_min=float(n.min()),
  gap_min_bar=float(abs(delta).min()),Ib_A=float(cm[0]),Is_A=float(cm[1]),
  physical=bool(np.all((p>=0)&(p<=1)) and np.all(n>=0) and np.all(cm[:2]>0)),
  checkpoint_time_s=float(meta["time_ps"])*1e-12,metadata=meta)
def event_record(manifest,final_time_s):
 obs=manifest.get("observation",{});confirmed=obs.get("trigger_confirmed") is True
 result=dict(confirmed=confirmed,censored=bool(obs.get("censored",False)),
  status=obs.get("status",manifest.get("status")),latency_s=None,crossing_time_s=None,
  confirmation_time_s=None,source="saved accepted-step trigger controller")
 if not confirmed:return result
 cross=obs.get("crossing_time_s",obs.get("first_crossing_time_s"))
 confirmation=obs.get("confirmation_time_s",obs.get("trigger_time_s"))
 latency=obs.get("latency_s")
 if latency is None and cross is not None and "event_time_s" in obs.get("config",{}):
  latency=float(cross)-float(obs["config"]["event_time_s"])
 if any(v is None or not np.isfinite(float(v)) for v in (cross,confirmation,latency)):
  raise ValueError("Confirmed event has no complete saved crossing/confirmation/latency")
 if not 0<=float(latency) or not float(cross)<=float(confirmation)<=final_time_s+1e-24:
  raise ValueError("Saved trigger contradicts its accepted observation interval")
 result.update(latency_s=float(latency),crossing_time_s=float(cross),confirmation_time_s=float(confirmation))
 return result
def node_areas(manifest,snap):
 geometry=manifest.get("geometry",{});path=geometry.get("source_path")
 if path and Path(path).is_file():
  if geometry.get("mesh_sha256") and sha256(path)!=geometry["mesh_sha256"]:
   raise ValueError("Pinned geometry changed after the run")
  with np.load(path,allow_pickle=False) as mesh:
   if "area_weights" in mesh and "ell0_m" in mesh:
    areas=mesh["area_weights"]*float(mesh["ell0_m"])**2
    if len(areas)!=len(snap["coordinates_m"]):raise ValueError("Mesh areas do not match saved nodes")
    return areas,"original Delaunay-Voronoi dual areas"
 triangles=np.asarray(snap.get("triangles",[]),int);xy=snap["coordinates_m"]
 if triangles.ndim!=2 or triangles.shape[1]!=3 or not len(triangles):
  return None,"unavailable; no area is invented"
 points=xy[triangles];u=points[:,1]-points[:,0];v=points[:,2]-points[:,0]
 cross=u[:,0]*v[:,1]-u[:,1]*v[:,0]
 area=abs(cross)/2;areas=np.zeros(len(xy))
 np.add.at(areas,triangles.ravel(),np.repeat(area/3,3))
 return areas,"triangle barycentric areas, analysis fallback"
def hotbelt_metrics(snap,areas,*,fraction=.2):
 t=increasing_time(snap["time_s"],"snapshots.time_s");xy=finite("coordinates_m",snap["coordinates_m"])
 gap=abs(finite("delta_J",snap["delta_J"]))/MEV_J
 Te=finite("Te_K",snap["Te_K"]);Tph=finite("Tph_K",snap["Tph_K"])
 if xy.ndim!=2 or xy.shape[1]!=2 or gap.shape!=(len(t),len(xy)) or Te.shape!=gap.shape or Tph.shape!=gap.shape:
  raise ValueError("Snapshot geometry or physical fields have inconsistent shapes")
 if np.any(Te<0) or np.any(Tph<0):raise ValueError("Negative equivalent temperature")
 initial=gap[0]
 if np.any(initial<=0):raise ValueError("Positive fresh-DC gap reference required")
 edges=np.asarray(snap.get("edges",[]),int)
 if edges.size==0:
  tri=np.asarray(snap.get("triangles",[]),int)
  edges=np.unique(np.sort(np.vstack((tri[:,[0,1]],tri[:,[1,2]],tri[:,[2,0]])),axis=1),axis=0) if tri.size else np.empty((0,2),int)
 if edges.ndim!=2 or edges.shape[1]!=2 or np.any(edges<0) or np.any(edges>=len(xy)):
  raise ValueError("Invalid snapshot graph edges")
 y=xy[:,1];tol=max(np.ptp(y)*1e-10,1e-18)
 low=np.flatnonzero(abs(y-y.min())<=tol);high=np.flatnonzero(abs(y-y.max())<=tol)
 center=np.array([(xy[:,0].min()+xy[:,0].max())/2,(y.min()+y.max())/2])
 central=int(np.argmin(np.sum((xy-center)**2,axis=1)))
 rows=[]
 for index,time_s in enumerate(t):
  mask=gap[index]<=fraction*initial;selected=edges[mask[edges].all(axis=1)]
  spanning=False;extent=None;component_area=0.
  if len(selected):
   adj=coo_matrix((np.ones(2*len(selected)),(np.r_[selected[:,0],selected[:,1]],
    np.r_[selected[:,1],selected[:,0]])),shape=(len(xy),len(xy))).tocsr()
   count,labels=connected_components(adj,directed=False)
   shared=set(labels[low[mask[low]]])&set(labels[high[mask[high]]])
   if shared:
    spanning=True;members=np.flatnonzero(mask&np.isin(labels,list(shared)))
    extent=float(np.ptp(xy[members,0])*1e9)
    component_area=None if areas is None else float(areas[members].sum()*1e18)
  rows.append(dict(time_s=float(time_s),time_ps=float(time_s*1e12),
   gap_min_meV=float(gap[index].min()),gap_center_meV=float(gap[index,central]),
   gap_max_meV=float(gap[index].max()),Te_max_K=float(Te[index].max()),Tph_max_K=float(Tph[index].max()),
   suppressed_area_nm2=None if areas is None else float(areas[mask].sum()*1e18),
   gap_deficit_meV_nm2=None if areas is None else float(np.sum(areas*np.maximum(initial-gap[index],0))*1e18),
   width_spanning=spanning,spanning_component_x_extent_nm=extent,
   spanning_component_area_nm2=component_area))
 return rows


def audit_run(name,path,*,energy_limit=.01,hotbelt_fraction=.2,dc_reference_energy_J=None):
 path=Path(path).resolve();manifest=json.loads((path/"manifest.json").read_text())
 history=load_npz(path/"trajectory.npz");snap=load_npz(path/"snapshots.npz")
 t=increasing_time(history["time_s"],"trajectory.time_s")
 for key in ("Vout_V","Vdev_V","Ib_A","Is_A","energy_total_J","energy_defect_J",
  "external_energy_integral_J","instant_power_residual_W","escape_power_W","reservoir_power_out_W"):
  history[key]=require_series(history,key,t)
 E0=float(manifest.get("energy_reference",{}).get("initial_total_energy_J",history["energy_total_J"][0]))
 recomputed=history["energy_total_J"]-E0-history["external_energy_integral_J"]
 mismatch=float(np.max(abs(recomputed-history["energy_defect_J"])))
 roundoff=64*np.finfo(float).eps*max(float(max(abs(history["energy_total_J"]))),1e-300)
 retained=float(manifest.get("source",{}).get("retained_energy_J",0.))
 if retained<0:raise ValueError("Negative retained source energy")
 reference=retained if retained>0 else dc_reference_energy_J
 absolute_work=np.r_[abs(history["external_energy_integral_J"][0]),
  abs(np.diff(history["external_energy_integral_J"]))].cumsum()
 energy_applicable=reference is not None and reference>0
 relative=abs(recomputed)/(float(reference)+absolute_work) if energy_applicable else np.full(len(t),np.nan)
 support=checkpoint_support(path/"checkpoint_latest.npz")
 checkpoint_time_consistent=abs(support["checkpoint_time_s"]-t[-1])<=max(1e-24,64*np.finfo(float).eps*abs(t[-1]))
 checksums={}
 for file,expected in manifest.get("output_sha256",{}).items():
  checksums[file]=sha256(path/file)==expected
 areas,area_method=node_areas(manifest,snap);hotbelt=hotbelt_metrics(snap,areas,fraction=hotbelt_fraction)
 snapshot_t=np.array([row["time_s"] for row in hotbelt])
 if snapshot_t[0]<t[0]-1e-24 or snapshot_t[-1]>t[-1]+1e-24:raise ValueError("Snapshot times escape saved trajectory")
 initial_source=None;source_ok=True
 if (path/"source_initial_state.npz").is_file():
  source=load_npz(path/"source_initial_state.npz")
  pp=finite("initial.p",source["p"]);nn=finite("initial.n",source["n"])
  target=finite("node_deposition_J",source["node_deposition_J"])
  physical=bool(np.all((pp>=0)&(pp<=1)) and np.all(nn>=0) and np.all(target>=0))
  total=float(target.sum());defect=total-retained
  source_ok=physical and abs(defect)<=1e-11*max(retained,1e-300)
  initial_source=dict(physical=physical,p_min=float(pp.min()),p_max=float(pp.max()),n_min=float(nn.min()),
   node_deposition_sum_J=total,declared_retained_energy_J=retained,energy_defect_J=defect)
 initial_density_energy=None
 initial_density_ok=retained==0
 if retained>0 and areas is not None and "u_excess_J_m3" in snap:
  density=finite("u_excess_J_m3",snap["u_excess_J_m3"])
  thickness=float(manifest.get("physicalparams",{}).get("film_thickness_m",0.))
  if density.shape!=(len(snapshot_t),len(areas)) or thickness<=0:
   raise ValueError("Physical energy density/film thickness does not match snapshot geometry")
  integrated=float(np.sum(density[0]*areas)*thickness)
  density_defect=integrated-retained
  initial_density_ok=abs(density_defect)<=1e-9*retained
  initial_density_energy=dict(integral_J=integrated,retained_energy_J=retained,
   defect_J=density_defect,relative_defect=abs(density_defect)/retained,
   method="Saved initial full-functional excess-energy density integrated over original nodal film volumes")
 limiting=None
 if "limited_event_energy_J" in history and "totalabs_event_energy_J" in history:
  limited=require_series(history,"limited_event_energy_J",t)
  totalabs=require_series(history,"totalabs_event_energy_J",t)
  limiting=dict(cumulative_limited_event_energy_J=float(limited[-1]),
   cumulative_absolute_event_energy_J=float(totalabs[-1]),
   absolute_event_energy_fraction=float(limited[-1]/totalabs[-1]) if totalabs[-1]>0 else 0.,
   limited_energy_per_retained_source=None if retained==0 else float(limited[-1]/retained),
   acceptance="Diagnostic only; its effect on observables is tested by temporal refinement")
 event=event_record(manifest,float(t[-1]))
 event_consistent=True
 if event["confirmed"]:
  obs=manifest["observation"];config=obs.get("config",{})
  threshold=float(obs.get("threshold_voltage_V",config.get("threshold_V",0.)))
  baseline=float(obs.get("baseline_Vout_V",config.get("baseline_V",0.)))
  sign=1. if config.get("polarity","positive")=="positive" else -1.
  excursion=sign*(float(np.interp(event["crossing_time_s"],t,history["Vout_V"]))-baseline)
  event["saved_crossing_voltage_difference_V"]=excursion-threshold
  event_consistent=threshold>0 and abs(excursion-threshold)<=max(1e-14,threshold*1e-10)
 circuit=manifest.get("physicalparams",{}).get("circuit",{})
 Rb=float(circuit["R_bias_ohm"]);Rl=float(circuit["R_load_ohm"]);Vb=float(circuit["V_bias_V"])
 if not all(np.isfinite(v) for v in (Rb,Rl,Vb)) or min(Rb,Rl)<=0:
  raise ValueError("Finite physical circuit and positive resistances required for instantaneous power audit")
 power_device=history["Is_A"]*history["Vdev_V"]
 power_source=Vb*history["Ib_A"]
 power_bias=Rb*history["Ib_A"]**2
 power_load=Rl*(history["Ib_A"]-history["Is_A"])**2
 power_circuit_net=power_source-power_bias-power_load
 power_scale=(abs(power_device)+abs(history["escape_power_W"])+
  abs(history["reservoir_power_out_W"])+abs(power_circuit_net))
 power_tolerance=1e-20+1e-10*power_scale
 power_absres=abs(history["instant_power_residual_W"])
 power_ratio=power_absres/power_tolerance
 worst_power_index=int(np.argmax(power_ratio))
 instant_power_ok=bool(np.all(power_absres<=power_tolerance))
 instant_power=dict(absolute_tolerance_W=1e-20,relative_tolerance=1e-10,
  normalization="abs(Is*Vdev)+abs(Pescape)+abs(Pcuts)+abs(Vb*Ib-Rb*Ib^2-Rl*(Ib-Is)^2)",
  maximum_residual_over_tolerance=float(power_ratio[worst_power_index]),
  minimum_allowed_minus_residual_W=float(np.min(power_tolerance-power_absres)),
  worst_ratio_time_ps=float(t[worst_power_index]*1e12),
  scale_at_worst_ratio_W=float(power_scale[worst_power_index]),
  tolerance_at_worst_ratio_W=float(power_tolerance[worst_power_index]),
  residual_at_worst_ratio_W=float(history["instant_power_residual_W"][worst_power_index]),
  maximum_power_scale_W=float(power_scale.max()),
  maximum_circuit_net_power_W=float(np.max(abs(power_circuit_net))),
  scope="Saved independent full film+CM instantaneous power residual; DC source/bias cancellation is not used to inflate the scale")
 energy_pass=bool(energy_applicable and float(np.max(relative))<=energy_limit)
 checks=dict(complete_declared_window=manifest.get("run_complete") is True,
  successful_saved_window_status=manifest.get("status") in {"TRIGGER_MARGIN_COMPLETE",
   "TIME_CENSORED_NO_TRIGGER","DEVELOPMENT_DURATION_REACHED","DC_HOLD_COMPLETED"},
  DC_baseline_within_nV=retained>0 or bool(np.max(abs(history["Vout_V"]))<=1e-9),
  energy_ledger_recomputation=mismatch<=roundoff,energy_relative_defect=energy_pass,
  instantaneous_energy_balance=instant_power_ok,
  checkpoint_support=support["physical"],checkpoint_time_consistent=checkpoint_time_consistent,
  positive_global_currents=bool(np.all(history["Ib_A"]>0) and np.all(history["Is_A"]>0)),
  positive_snapshot_temperatures=True,source_initial_support_and_deposition=source_ok,
  saved_output_checksums={"trajectory.npz","snapshots.npz","checkpoint_latest.npz"}.issubset(checksums) and all(checksums.values()),
  independent_initial_source_energy=initial_density_ok,
  initial_photon_state_present=retained==0 or initial_source is not None,
  saved_electrical_event_consistent=event_consistent)
 guard=audit_local_energy_guard(manifest,history,support["metadata"],retained)
 checks["independent_local_energy_guard"]=guard["passed"]
 spatial=None
 if manifest.get("source",{}).get("center_m") is not None:
  spatial=spatial_snapshot_diagnostics(manifest,snap)
 spans=[index for index,row in enumerate(hotbelt) if row["width_spanning"]]
 first=None if not spans else spans[0]
 belt_record=None if first is None else dict(first_saved_time_s=hotbelt[first]["time_s"],
  bracket_s=[hotbelt[max(0,first-1)]["time_s"],hotbelt[first]["time_s"]])
 baseline=float(manifest.get("observation",{}).get("baseline_Vout_V",0.))
 pulse=history["Vout_V"]-baseline
 summary=dict(name=name,path=str(path),passed=all(checks.values()),checks=checks,
  stage5_admission_inferred=False,status=manifest.get("status"),declared_mode=manifest.get("observation",{}).get("mode"),
  final_time_ps=float(t[-1]*1e12),energy=dict(retained_energy_J=retained,reference_energy_J=reference,
   normalization="Eretained+sum(abs(accepted midpoint external-energy increments))" if retained>0 else "explicit DC reference energy plus absolute external work",
   applicable=energy_applicable,relative_limit=energy_limit,maximum_relative_defect=None if not energy_applicable else float(np.max(relative)),
   maximum_absolute_defect_J=float(np.max(abs(recomputed))),absolute_external_work_J=float(absolute_work[-1]),
   saved_ledger_mismatch_J=mismatch,ledger_roundoff_bound_J=roundoff,
   maximum_instant_power_residual_W=float(np.max(abs(history["instant_power_residual_W"])))),
  metrics=dict(Vout_peak_mV=float(pulse.max()*1e3),Vout_min_mV=float(pulse.min()*1e3),
   Vdev_peak_mV=float(history["Vdev_V"].max()*1e3),Ib_min_uA=float(history["Ib_A"].min()*1e6),
   Is_min_uA=float(history["Is_A"].min()*1e6),gap_min_meV=float(min(row["gap_min_meV"] for row in hotbelt)),
   Te_max_K=float(max(row["Te_max_K"] for row in hotbelt)),Tph_max_K=float(max(row["Tph_max_K"] for row in hotbelt))),
  local_energy_guard=guard,spatial_diagnostics=None if spatial is None else {key:value for key,value in spatial.items() if key not in {"time_s","deficit_J_m2"}},
  event=event,checkpoint_support={key:v for key,v in support.items() if key!="metadata"},
  instantaneous_power=instant_power,initial_source=initial_source,initial_source_density_integral=initial_density_energy,output_checksums=checksums,numerical_event_limiting=limiting,
  positivity_scope="Initial/full final distributions directly verified; intermediate accepted distributions are enforced by the frozen integrator, not reconstructed from sparse snapshots",
  hotbelt=dict(gap_fraction_of_initial_local_DC=hotbelt_fraction,
   definition="Connected suppressed-gap component touching both transverse strip boundaries; diagnostic only, never an electrical detection event",
   first_saved_width_spanning=belt_record,area_method=area_method))
 traces=dict(time_s=t,energy_relative_defect=relative,energy_recomputed_defect_J=recomputed,
  absolute_external_work_J=absolute_work,Vout_relative_V=pulse,Ib_A=history["Ib_A"],Is_A=history["Is_A"],
  instant_power_residual_W=history["instant_power_residual_W"],instant_power_scale_W=power_scale,
  instant_power_tolerance_W=power_tolerance,instant_power_residual_over_tolerance=power_ratio)
 return dict(summary=summary,manifest=manifest,history=history,snap=snap,hotbelt=hotbelt,traces=traces,spatial=spatial)


def audit_local_energy_guard(manifest,history,metadata,retained_J):
 """Independent endpoint/work reconstruction; no retrospective legacy guarding."""
 identity=manifest.get("parameter_identity",{});numerical=manifest.get("numerical",{})
 declared=("local_energy_tolerance" in identity or "local_energy_tolerance" in numerical
  or numerical.get("local_energy_guard",{}).get("enabled") is True)
 if not declared:return dict(applicable=False,passed=True,scope="Legacy data: no local guard is certified")
 required=("absolute_external_energy_integral_J","sum_abs_local_energy_defect_J",
  "max_accepted_local_energy_fraction","max_guarded_accepted_local_energy_fraction",
  "guard_active_since_time_ps","local_energy_reference_energy_J","local_energy_rejections","energy_guard_ledger_origin")
 if not all(key in metadata for key in required):
  return dict(applicable=True,passed=False,error="Incomplete checkpoint energy-guard ledger")
 if "local_energy_tolerance" not in identity:
  return dict(applicable=True,passed=False,error="Energy guard missing from strict parameter identity")
 tol=float(identity["local_energy_tolerance"])
 saved_tol=float(numerical.get("local_energy_tolerance",np.nan))
 vals=[float(metadata[key]) for key in required[:-1]]
 if not np.isfinite(vals).all() or min(vals)<0 or not np.isfinite(tol) or tol<=0:
  return dict(applicable=True,passed=False,error="Invalid nonfinite/negative energy-guard metadata")
 reference=float(metadata["local_energy_reference_energy_J"])
 if reference<=0:return dict(applicable=True,passed=False,error="Positive guard energy reference required")
 t=np.asarray(history["time_s"]);energy=np.asarray(history["energy_total_J"])
 work=np.asarray(history["external_energy_integral_J"]);dw=np.diff(work);de=np.diff(energy)
 defect=de-dw
 absolute=np.cumsum(abs(dw));denominator=reference+absolute;fractions=abs(defect)/denominator
 anchor=float(metadata["guard_active_since_time_ps"])*1e-12
 time_roundoff=32*np.finfo(float).eps*max(float(abs(t).max()),1e-12)
 anchor_indices=np.flatnonzero(abs(t-anchor)<=time_roundoff)
 post=t[:-1]>=anchor-time_roundoff
 # Persisted energies are binary64. This bound explicitly allows only arithmetic
 # rounding of endpoint subtraction and cumulative-work differencing.
 eps=np.finfo(float).eps
 roundoff_J=16*eps*(abs(energy[:-1])+abs(energy[1:])+abs(work[:-1])+abs(work[1:]))
 fraction_roundoff=roundoff_J/denominator
 global_max=float(fractions.max()) if len(fractions) else 0.
 guarded_max=float(fractions[post].max()) if np.any(post) else 0.
 max_roundoff=float(fraction_roundoff.max()) if len(fractions) else 0.
 guarded_roundoff=float(fraction_roundoff[post].max()) if np.any(post) else 0.
 abs_work=math.fsum(abs(float(v)) for v in dw)
 sum_defect=math.fsum(abs(float(v)) for v in defect)
 # Same endpoints are used by runner and reader; sum mismatch arises only from
 # subtracting persisted cumulative W instead of original midpoint receipts.
 work_sum_roundoff_J=16*eps*math.fsum(float(v) for v in abs(work[:-1])+abs(work[1:]))
 max_roundoff_W=max(work_sum_roundoff_J,32*eps*max(reference,abs_work))
 saved=manifest.get("local_energy_control",{})
 matching_manifest=all(key in saved and saved[key]==metadata[key] for key in required)
 count=float(metadata["local_energy_rejections"])
 checks=dict(strict_parameter_and_manifest_tolerance=bool(np.isfinite(saved_tol) and saved_tol==tol),
  positive_photon_reference_matches_retained=retained_J<=0 or math.isclose(reference,retained_J,rel_tol=1e-13,abs_tol=0.),
  anchor_is_saved_accepted_time=bool(len(anchor_indices)==1 and anchor>=-time_roundoff and anchor<=t[-1]+time_roundoff),
  post_anchor_local_defects_bounded=bool(np.all(fractions[post]<=tol+fraction_roundoff[post])),
  absolute_work_ledger=abs(float(metadata["absolute_external_energy_integral_J"])-abs_work)<=max_roundoff_W,
  absolute_defect_ledger=abs(float(metadata["sum_abs_local_energy_defect_J"])-sum_defect)<=max_roundoff_W,
  maximum_all_ledger=abs(float(metadata["max_accepted_local_energy_fraction"])-global_max)<=max_roundoff,
  maximum_guarded_ledger=abs(float(metadata["max_guarded_accepted_local_energy_fraction"])-guarded_max)<=guarded_roundoff,
  rejection_count_valid=count==int(count),
  checkpoint_manifest_guard_ledger_matches=matching_manifest)
 return dict(applicable=True,passed=all(checks.values()),checks=checks,
  tolerance=tol,reference_energy_J=reference,guard_active_since_time_ps=anchor*1e12,
  post_anchor_accepted_steps=int(post.sum()),legacy_accepted_steps=int((~post).sum()),
  maximum_all_reconstructed_fraction=global_max,maximum_guarded_reconstructed_fraction=guarded_max,
  maximum_post_anchor_excess_fraction=float(np.max(fractions[post]-tol-fraction_roundoff[post])) if np.any(post) else None,
  absolute_work_reconstructed_J=abs_work,sum_abs_local_defect_reconstructed_J=sum_defect,
  ledger_sum_roundoff_bound_J=max_roundoff_W,maximum_fraction_roundoff=max_roundoff,
  maximum_guarded_fraction_roundoff=guarded_roundoff,
  roundoff_definition="Per step16eps*(abs(Eold)+abs(Enew)+abs(Wold)+abs(Wnew)); ledger sum16eps*sum(abs(Wold)+abs(Wnew)), with32eps*max(Eref,sumabsdW) lower arithmetic bound",
  denominator="Ereference+cumsum(abs(diff(saved cumulative midpoint external work)))",
  scope="Only accepted steps starting at/after the saved guard anchor are tested against local tolerance; legacy maxima remain visible; 1% global energy gate remains independent")

def spatial_snapshot_diagnostics(manifest,snap):
 """Nodal deficit and complex-interpolant topology; no electrical event inferred."""
 from sandbox.stage5_effective_validation import vortex_diagnostics as vd
 xy=np.asarray(snap["coordinates_m"],float);tri,_=vd.ccw_triangles(xy,snap["triangles"])
 center=manifest.get("source",{}).get("center_m")
 if center is None:raise ValueError("Saved physical source.center_m required for robust spatial comparison")
 center=finite("source.center_m",center)
 if np.asarray(center).shape!=(2,):raise ValueError("Planar photon source center required")
 if abs(xy[:,1].min()+40e-9)>1e-15 or abs(xy[:,1].max()-40e-9)>1e-15:
  raise ValueError("The fixed physical robust-comparison ROI requires strip y in [-40,40]nm")
 weights=vd.central_linear_weights(xy,tri,center)
 expected=120e-9*80e-9
 if abs(float(weights.sum())-expected)>1e-10*expected:
  raise ValueError("Central physical120x80nm analysis region is not fully inside the saved domain")
 d=finite("snapshots.delta_J",snap["delta_J"]);d0=abs(d[0])
 deficit=np.maximum(d0[None,:]-abs(d),0)@weights
 rows=[]
 for index,time_s in enumerate(snap["time_s"]):
  minimum,_,_=vd.complex_triangle_minimum(d[index],tri)
  phase_defined=bool(np.all(abs(d[index])>0))
  positive=negative=total=None
  if phase_defined:
   charges,_,_,_=vd.triangle_topology(d[index],xy,tri)
   positive=int((charges>0).sum());negative=int((charges<0).sum());total=int(charges.sum())
  rows.append(dict(time_s=float(time_s),time_ps=float(time_s*1e12),
   nodal_minimum_meV=float(abs(d[index]).min()/MEV_J),
   complex_linear_minimum_meV=float(minimum.min()/MEV_J),
   central_nodal_deficit_meV_nm2=float(deficit[index]/MEV_J*1e18),
   phase_winding_defined=phase_defined,n_positive=positive,n_negative=negative,total_charge=total))
 return dict(rows=rows,time_s=np.asarray(snap["time_s"]),deficit_J_m2=deficit,
  central_origin_x_m=float(center[0]),physical_y_origin_m=0.,
  region_relative_nm=[-60,60,-40,40],central_area_m2=float(weights.sum()),
  DC_gap_scale_J=float(weights@d0/weights.sum()),
  DC_gap_area_scale_J_m2=float(weights@d0),
  deficit_definition="D_h,nodal=sum exact clipped linear-nodal weights * positive nodal gap deficit; spatial approximation to continuum, not integral of magnitude of complex interpolant",
  topology_scope="A=0,CCW winding of actual saved states only; zeros of interpolant do not prove resolved core or electrical detection")

def complex_ROI_sampler(snap,origin_x_m,spacing_nm=.5):
 """Fixed physical quadrature, saved triangulation, complex interpolation first."""
 import matplotlib.tri as mtri
 xy=np.asarray(snap["coordinates_m"],float).copy();xy[:,0]-=origin_x_m
 triangles=np.asarray(snap["triangles"],int)
 nx=int(np.ceil(120/spacing_nm));ny=int(np.ceil(80/spacing_nm))
 x=(-60+(np.arange(nx)+.5)*120/nx)*1e-9
 y=(-40+(np.arange(ny)+.5)*80/ny)*1e-9
 xx,yy=np.meshgrid(x,y);points=np.column_stack((xx.ravel(),yy.ravel()))
 mesh=mtri.Triangulation(xy[:,0],xy[:,1],triangles)
 labels=np.asarray(mesh.get_trifinder()(points[:,0],points[:,1]),int)
 if np.any(labels<0):raise ValueError("Common physical field quadrature extends outside a saved domain")
 nodes=triangles[labels];p=xy[nodes]
 a=p[:,1]-p[:,0];b=p[:,2]-p[:,0];r=points-p[:,0]
 cross=lambda u,v:u[:,0]*v[:,1]-u[:,1]*v[:,0]
 determinant=cross(a,b)
 if np.any(determinant==0):raise ValueError("Degenerate interpolation triangle")
 l1=cross(r,b)/determinant;l2=cross(a,r)/determinant
 bary=np.column_stack((1-l1-l2,l1,l2))
 return dict(nodes=nodes,barycentric=bary,sample_spacing_actual_nm=[120/nx,80/ny],
  samples=len(nodes),rule="Uniform midpoint spatial quadrature on fixed physical120x80nm ROI")

def complex_saved_frame_at(snap,time_s):
 """Linear comparison in time only; never used for winding or trigger."""
 t=np.asarray(snap["time_s"]);j=int(np.searchsorted(t,time_s))
 if j==0:return snap["delta_J"][0]
 if j==len(t):return snap["delta_J"][-1]
 fraction=(time_s-t[j-1])/(t[j]-t[j-1])
 return (1-fraction)*snap["delta_J"][j-1]+fraction*snap["delta_J"][j]

def compare_robust_spatial(reference,candidate,grid):
 a,b=reference["spatial"],candidate["spatial"]
 norm=float(a["DC_gap_scale_J"])
 if norm<=0:raise ValueError("Positive initial DC gap required for field normalization")
 ra=np.interp(grid,a["time_s"],a["deficit_J_m2"])
 cb=np.interp(grid,b["time_s"],b["deficit_J_m2"])
 difference=float(abs(cb-ra).max());peak=float(ra.max())
 peak_error=(0. if difference==0 else np.inf) if peak==0 else difference/peak
 sampler_a=complex_ROI_sampler(reference["snap"],a["central_origin_x_m"])
 sampler_b=complex_ROI_sampler(candidate["snap"],b["central_origin_x_m"])
 l2=[]
 for time_s in grid:
  fields=[]
  for run,sampler in ((reference,sampler_a),(candidate,sampler_b)):
   frame=complex_saved_frame_at(run["snap"],time_s)
   fields.append(abs(np.sum(frame[sampler["nodes"]]*sampler["barycentric"],axis=1)))
  l2.append(float(np.sqrt(np.mean((fields[1]-fields[0])**2))/norm))
 l2=np.asarray(l2)
 metric=dict(region_relative_nm=[-60,60,-40,40],x_center="saved source.center_m[0]",y_coordinate="physical",
  complex_interpolation="Saved linear complex triangulation first, magnitude afterwards; no vortex alignment",
  time_interpolation="Piecewise linear complex saved-field comparison only; winding and event are never computed on interpolated states",
  quadrature=sampler_a["rule"],sample_spacing_actual_nm=sampler_a["sample_spacing_actual_nm"],
  L2_field=dict(maximum_relative_error=float(l2.max()),DC_gap_scale_meV=norm/MEV_J),
  deficit=dict(maximum_absolute_difference_meV_nm2=difference/MEV_J*1e18,
   reference_actual_peak_meV_nm2=peak/MEV_J*1e18,
   relative_error_actual_peak=float(peak_error),
   relative_error_DCarea=difference/float(a["DC_gap_area_scale_J_m2"]),
   zero_peak_rule="Zero difference/zero peak ->0; nonzero difference/zero peak ->failure, no arbitrary floor",
   definition=a["deficit_definition"]))
 curves=dict(field_L2_relative=l2,
  reference_central_nodal_deficit_meV_nm2=ra/MEV_J*1e18,
  candidate_central_nodal_deficit_meV_nm2=cb/MEV_J*1e18)
 return metric,curves


def common_grid(a,b):
 start=max(float(a[0]),float(b[0]));end=min(float(a[-1]),float(b[-1]))
 if end<=start:raise ValueError("No nonzero common saved time interval")
 return np.unique(np.r_[start,a[(a>=start)&(a<=end)],b[(b>=start)&(b<=end)],end])
def compare_runs(kind,reference,candidate):
 if kind not in KINDS:raise ValueError("Unknown comparison kind: "+kind)
 ar,br=reference["summary"],candidate["summary"];am,bm=reference["manifest"],candidate["manifest"]
 ai,bi=am.get("parameter_identity",{}),bm.get("parameter_identity",{})
 changed={key:{"reference":ai.get(key),"candidate":bi.get(key)}
  for key in set(ai)|set(bi) if ai.get(key)!=bi.get(key)}
 forbidden={key:value for key,value in changed.items() if key not in RUNTIME_KEYS|REFINEMENT_KEYS[kind]}
 declared_control=declared_temporal_guard_control(kind,ai,bi)
 if declared_control is not None:forbidden.pop("local_energy_tolerance",None)
 source_match=bool(am.get("source_pins_sha256")) and am.get("source_pins_sha256")==bm.get("source_pins_sha256")
 tg=common_grid(reference["traces"]["time_s"],candidate["traces"]["time_s"])
 rv=np.interp(tg,reference["traces"]["time_s"],reference["traces"]["Vout_relative_V"])
 cv=np.interp(tg,candidate["traces"]["time_s"],candidate["traces"]["Vout_relative_V"])
 peak=float(max(abs(rv)));threshold=float(am.get("observation",{}).get("threshold_voltage_V",0.))
 voltage_scale=max(peak,threshold)
 if voltage_scale<=0:raise ValueError("A positive saved voltage threshold or pulse scale is required")
 vdiff=float(max(abs(rv-cv)));voltage_error=vdiff/voltage_scale
 ast=np.array([r["time_s"] for r in reference["hotbelt"]]);bst=np.array([r["time_s"] for r in candidate["hotbelt"]])
 gt=common_grid(ast,bst);ag=np.array([r["gap_min_meV"] for r in reference["hotbelt"]])
 bg=np.array([r["gap_min_meV"] for r in candidate["hotbelt"]])
 rg=np.interp(gt,ast,ag);cg=np.interp(gt,bst,bg)
 gap_scale=float(np.median(abs(reference["snap"]["delta_J"][0]))/MEV_J)
 gdiff=float(max(abs(rg-cg)));gap_error=gdiff/gap_scale;limit=KINDS[kind]
 ae,be=ar["event"],br["event"];latency=dict(available=False,passed=None,
  reason="No pair of saved confirmed electrical events; no latency is assigned")
 latency_check=True;latency_complete=False
 if ae["confirmed"] and be["confirmed"]:
  delta=abs(ae["latency_s"]-be["latency_s"])
  tolerance=max(limit*abs(ae["latency_s"]),.1e-12 if kind=="temporal" else 0.)
  latency=dict(available=True,passed=delta<=tolerance,reference_ps=ae["latency_s"]*1e12,
   candidate_ps=be["latency_s"]*1e12,absolute_difference_ps=delta*1e12,tolerance_ps=tolerance*1e12)
  latency_check=latency["passed"];latency_complete=True
 elif ae["confirmed"]!=be["confirmed"]:
  observed=ae if ae["confirmed"] else be
  other_end=br["final_time_ps"]*1e-12 if ae["confirmed"] else ar["final_time_ps"]*1e-12
  if other_end>=observed["confirmation_time_s"]:
   latency_check=False
   latency.update(passed=False,reason="Only one run confirms detection despite both covering the confirmation time")
  else:latency.update(reason="The shorter saved window ends before the other confirmed event; comparison pending")
 declared_refinement_ok=True
 if kind in {"mesh","domain"}:
  key="spacing_nm" if kind=="mesh" else "length_nm"
  av=float(ai.get(key,np.nan));bv=float(bi.get(key,np.nan))
  declared_refinement_ok=bool(np.isfinite(av) and np.isfinite(bv) and av>0 and bv>0 and
   (bv<av if kind=="mesh" else bv>av))
 elif kind=="spectral":
  pairs=[(float(ai[key]),float(bi[key])) for key in REFINEMENT_KEYS[kind] if key in ai and key in bi]
  declared_refinement_ok=bool(pairs and all(np.isfinite(av) and np.isfinite(bv) and av>0 and bv>=av
   for av,bv in pairs) and any(bv>av for av,bv in pairs))
 temporal_ratio=None;refinement_ok=True
 if kind=="temporal":
  temporal_ratio=float(ai.get("dt_ps",np.nan))/float(bi.get("dt_ps",np.nan))
  refinement_ok=bool(np.isfinite(temporal_ratio) and temporal_ratio>=2-1e-12)
 spatial=None;spatial_curves={}
 spatial_available=reference.get("spatial") is not None and candidate.get("spatial") is not None
 if spatial_available:spatial,spatial_curves=compare_robust_spatial(reference,candidate,gt)
 checks=dict(reference_saved_window_passed=ar["passed"],candidate_saved_window_passed=br["passed"],
  same_frozen_model_sources=source_match,same_nonrefined_parameters=not forbidden,
  actual_declared_refinement=declared_refinement_ok,
  voltage_waveform=voltage_error<=limit,
  electrical_latency_if_observed=latency_check,temporal_refinement_at_least_two=refinement_ok)
 if kind in {"mesh","domain"}:
  checks["common_physical_ROI_available"]=spatial_available
  checks["central_gap_magnitude_L2"]=spatial_available and spatial["L2_field"]["maximum_relative_error"]<=limit
  checks["central_deficit_over_reference_peak"]=spatial_available and spatial["deficit"]["relative_error_actual_peak"]<=limit
 else:checks["gap_min_waveform"]=gap_error<=limit
 legacy_checks={key:value for key,value in checks.items() if key not in
  {"common_physical_ROI_available","central_gap_magnitude_L2","central_deficit_over_reference_peak"}}
 legacy_checks["gap_min_waveform"]=gap_error<=limit
 result=dict(kind=kind,reference=ar["name"],candidate=br["name"],passed=all(checks.values()),
  checks=checks,relative_limit=limit,stage5_admission_inferred=False,
  policy_version="spatial_observables_v2",legacy_nodal_policy_passed=all(legacy_checks.values()),
  legacy_nodal_policy_checks=legacy_checks,robust_spatial=spatial,
  common_trajectory_interval_ps=[float(tg[0]*1e12),float(tg[-1]*1e12)],
  common_gap_snapshot_interval_ps=[float(gt[0]*1e12),float(gt[-1]*1e12)],
  interpolation="Piecewise linear saved-data comparison only; no physical state or trigger is interpolated by this auditor",
  changed_parameters=changed,forbidden_parameter_changes=forbidden,step_refinement_ratio=temporal_ratio,
  declared_numerical_control=declared_control,
  accepted_step_statistics_after_both_guard_anchors=guarded_common_step_statistics(reference,candidate,tg[0],tg[-1]),
  accepted_step_statistics=dict(reference=accepted_step_statistics(reference,tg[0],tg[-1]),
   candidate=accepted_step_statistics(candidate,tg[0],tg[-1])) ,
  Vout=dict(maximum_absolute_difference_mV=vdiff*1e3,reference_peak_mV=peak*1e3,
   saved_threshold_mV=threshold*1e3,admission_scale_mV=voltage_scale*1e3,
   relative_error=voltage_error,pulse_scaled_error=None if peak==0 else vdiff/peak),
  gap_min=dict(maximum_absolute_difference_meV=gdiff,DC_reference_scale_meV=gap_scale,relative_error=gap_error,diagnostic_only=kind in {"mesh","domain"},former_gate_passed=gap_error<=limit),
  latency=latency,latency_comparison_complete=latency_complete)
 curves=dict(time_ps=tg*1e12,reference_Vout_mV=rv*1e3,candidate_Vout_mV=cv*1e3,
  difference_Vout_mV=(cv-rv)*1e3)
 gapcurves=dict(time_ps=gt*1e12,reference_gap_min_meV=rg,candidate_gap_min_meV=cg,difference_gap_min_meV=cg-rg)
 gapcurves.update(spatial_curves)
 return result,curves,gapcurves


def write_table(path,data):
 if isinstance(data,dict):
  keys=list(data);rows=zip(*(data[key] for key in keys))
 else:
  keys=list(data[0]) if data else [];rows=([row[key] for key in keys] for row in data)
 with Path(path).open("w",newline="") as stream:
  writer=csv.writer(stream);writer.writerow(keys)
  for row in rows:writer.writerow(["" if value is None else value for value in row])
def main(argv=None):
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument("--run",action="append",required=True,metavar="NAME=DIR")
 parser.add_argument("--compare",action="append",nargs=3,default=[],metavar=("KIND","REFERENCE","CANDIDATE"))
 parser.add_argument("--output-dir",type=Path,required=True)
 parser.add_argument("--energy-relative-limit",type=float,default=.01)
 parser.add_argument("--hotbelt-gap-fraction",type=float,default=.2)
 parser.add_argument("--dc-reference-energy-J",type=float,default=None)
 a=parser.parse_args(argv);output=a.output_dir.expanduser().resolve()
 if not output.is_relative_to(SCRATCH.resolve()):raise ValueError("Audit outputs must remain under scratch")
 if output.exists():raise FileExistsError("Existing audit evidence preserved; choose another output directory")
 if not 0<a.energy_relative_limit<1 or not 0<a.hotbelt_gap_fraction<1:raise ValueError("Invalid declared audit fraction")
 if a.dc_reference_energy_J is not None and a.dc_reference_energy_J<=0:raise ValueError("DC reference energy must be positive")
 runs={};errors={}
 for item in a.run:
  name,separator,path=item.partition("=")
  if not separator or not name or not name.replace("_","").replace("-","").isalnum():raise ValueError("Run requires a simple NAME=DIR")
  if name in runs or name in errors:raise ValueError("Duplicate run name")
  try:runs[name]=audit_run(name,path,energy_limit=a.energy_relative_limit,
   hotbelt_fraction=a.hotbelt_gap_fraction,dc_reference_energy_J=a.dc_reference_energy_J)
  except (ValueError,KeyError,OSError) as error:errors[name]=repr(error)
 comparisons=[]
 for kind,reference,candidate in a.compare:
  if reference not in runs or candidate not in runs:
   comparisons.append(dict(kind=kind,reference=reference,candidate=candidate,passed=False,
    error="Required run could not be audited"));continue
  try:
   result,curves,gapcurves=compare_runs(kind,runs[reference],runs[candidate])
   comparisons.append(result)
   result["_curves"]=curves;result["_gapcurves"]=gapcurves
  except (ValueError,KeyError) as error:
   comparisons.append(dict(kind=kind,reference=reference,candidate=candidate,passed=False,error=repr(error)))
 output.mkdir(parents=True)
 for name,run in runs.items():
  write_table(output/("trajectory_metrics_"+name+".csv"),run["traces"])
  write_table(output/("hotbelt_snapshots_"+name+".csv"),run["hotbelt"])
 for result in comparisons:
  prefix=result["kind"]+"_"+result["reference"]+"_"+result["candidate"]
  if "_curves" in result:
   write_table(output/(prefix+"_Vout.csv"),result.pop("_curves"))
   write_table(output/(prefix+"_gap.csv"),result.pop("_gapcurves"))
 report=dict(schema="pysnspd.effective_saved_data_audit.v2",stage5_admission_inferred=False,
  all_requested_checks_passed=not errors and all(run["summary"]["passed"] for run in runs.values())
   and all(result["passed"] for result in comparisons),
  criteria=dict(energy_relative_limit=a.energy_relative_limit,relative_limits=KINDS,
   instantaneous_energy_absolute_tolerance_W=1e-20,instantaneous_energy_relative_tolerance=1e-10,
   temporal_latency_absolute_margin_ps=.1,hotbelt_gap_fraction=a.hotbelt_gap_fraction,
   DC_reference_energy_J=a.dc_reference_energy_J,
   spatial_policy="Mesh/domain:Vout+latency+central magnitude L2/DC+central nodal-deficit/reference actual peak; temporal/spectral retain nodal-min gate",
   spatial_ROI_nm=[-60,60,-40,40],field_quadrature_spacing_nm=.5),
  scope="Saved finite windows; passing this audit does not establish experiment calibration, full reset, jitter statistics, or stage5 by itself",
  runs={name:run["summary"] for name,run in runs.items()},comparisons=comparisons,errors=errors,
  auditor_source_sha256=sha256(__file__),spatial_helper_sha256=sha256(Path(__file__).with_name("vortex_diagnostics.py")))
 (output/"audit.json").write_text(json.dumps(json_ready(report),indent=2,allow_nan=False)+"\n")
 print(json.dumps(dict(output_dir=str(output),all_requested_checks_passed=report["all_requested_checks_passed"],
  stage5_admission_inferred=False,run_count=len(runs),comparison_count=len(comparisons),errors=errors)))
 return report
if __name__=="__main__":main()
