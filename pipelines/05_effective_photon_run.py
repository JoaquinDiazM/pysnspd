"""Unified effective BCS/KWT/native-kinetic/CM photon pipeline in SI."""
from __future__ import annotations
import argparse
from dataclasses import asdict
import hashlib,json,math,os,signal,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from sandbox.stage4_core.parallel_runtime import limit_thread_environment
limit_thread_environment()
_early_parser=argparse.ArgumentParser(add_help=False)
_early_parser.add_argument("--threads",type=int,default=int(os.environ.get("NUMBA_NUM_THREADS","1")))
_early_threads,_early_unknown=_early_parser.parse_known_args()
if _early_threads.threads<1:raise ValueError("Positive numerical thread count required")
os.environ["NUMBA_NUM_THREADS"]=str(_early_threads.threads)
os.environ["NUMBA_CACHE_DIR"]="/home/jdiaz/scratch/tmp/pysnspd/numba_cache"
import numpy as np
from scipy.constants import Boltzmann as kB,hbar,elementary_charge as e,Planck,speed_of_light
from scipy.special import zeta
from pysnspd.experimental.effective_bcs_functional import EffectiveBCSFunctional,EffectiveBCSBulk,bcs_count_quadrature
from pysnspd.experimental.effective_bcs_preparation import get_effective_geometry,prepare_effective_dc,save_effective_dc,EffectiveBCSDC
from pysnspd.experimental.effective_bcs_kinetics import EffectiveBCSKinetics,acoustic_phonon_grid
from pysnspd.experimental.effective_bcs_ports import EffectiveNeutralPorts,freeze_exterior_inductance
from pysnspd.experimental.effective_material import allmaras_acoustic_material
from pysnspd.experimental.electrical_ports import ThesisCircuitParameters
from pysnspd.experimental.effective_coupled import EffectiveState,EffectiveCoupledSystem
from pysnspd.experimental.trigger_window import TriggerWindow,TriggerWindowConfig
SCRATCH=Path("/home/jdiaz/scratch");SCHEMA="pysnspd.effective_photon_run.v1"
HISTORY_KEYS=("time_s","Vout_V","Vdev_V","Ib_A","Is_A","vc_V","energy_total_J",
 "energy_defect_J","escape_power_W","domain_energy_J","instant_power_residual_W",
 "reservoir_power_out_W","external_energy_integral_J","limited_event_energy_J",
 "totalabs_event_energy_J","accepted_step_ps")
SNAPSHOT_KEYS=("time_s","delta_J","Te_K","Tph_K","u_excess_J_m3","electron_energy_J_m3","phonon_energy_J_m3")

def sha256(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def safe_json(v):
 if isinstance(v,Path):return str(v)
 if isinstance(v,np.ndarray):return safe_json(v.tolist())
 if isinstance(v,np.generic):return safe_json(v.item())
 if isinstance(v,float) and not math.isfinite(v):
  if math.isnan(v):raise ValueError("NaN forbidden in metadata")
  return "Infinity" if v>0 else "-Infinity"
 if isinstance(v,dict):return {str(k):safe_json(x) for k,x in v.items()}
 if isinstance(v,(list,tuple)):return [safe_json(x) for x in v]
 return v
def write_json(path,v):
 p=Path(path);tmp=p.with_name(p.name+".part")
 tmp.write_text(json.dumps(safe_json(v),indent=2,allow_nan=False)+"\n");os.replace(tmp,p)
def atomic_npz(path,arrays):
 p=Path(path);tmp=p.with_name(p.name+".part")
 with tmp.open("wb") as s:np.savez_compressed(s,**arrays)
 os.replace(tmp,p)
def scratch_path(value):
 p=Path(value).expanduser().resolve()
 if not p.is_relative_to(SCRATCH.resolve()):raise ValueError("Run data must remain under scratch")
 return p
def source_pins():
 names=["effective_bcs_functional.py","effective_bcs_preparation.py","effective_bcs_kinetics.py",
 "effective_bcs_ports.py","effective_coupled.py","effective_material.py","electrical_ports.py",
 "thermal_spatial_usadel.py","kinetic_events.py","trigger_window.py","effective_bcs_kinetics_numba.py"]
 paths=[ROOT/"pysnspd/experimental"/n for n in names]+[Path(__file__),ROOT/"sandbox/stage4_core/parallel_runtime.py"]
 return {str(p.relative_to(ROOT)):sha256(p) for p in paths}

class _LocalEnergyStepRejected(ArithmeticError):
 """Numerical trial rejection; no accepted physical state has been changed."""

def _local_energy_step(before_J,after_J,external_increment_J,reference_J,absolute_external_J):
 values=(before_J,after_J,external_increment_J,reference_J,absolute_external_J)
 if not all(math.isfinite(float(v)) for v in values):
  raise ArithmeticError("Nonfinite local-energy acceptance data")
 if reference_J<=0 or absolute_external_J<0:
  raise ValueError("Positive reference and nonnegative absolute work required")
 defect=float(after_J)-float(before_J)-float(external_increment_J)
 denominator=float(reference_J)+float(absolute_external_J)+abs(float(external_increment_J))
 return dict(local_energy_defect_J=defect,local_energy_denominator_J=denominator,
  local_energy_fraction=abs(defect)/denominator)

def _energy_guard_from_history(history,reference_J):
 """Measure legacy accepted increments; this does not certify past guarding."""
 energy=np.asarray(history["energy_total_J"],float)
 external=np.asarray(history["external_energy_integral_J"],float)
 if energy.ndim!=1 or external.shape!=energy.shape or not len(energy):
  raise ValueError("Coherent nonempty energy/work history required")
 if not np.all(np.isfinite(energy)) or not np.all(np.isfinite(external)):
  raise ValueError("Finite energy/work history required")
 if not math.isfinite(float(reference_J)) or reference_J<=0:
  raise ValueError("Positive finite local-energy reference required")
 dw=np.diff(external);de=np.diff(energy);defect=de-dw
 denominators=reference_J+np.cumsum(np.abs(dw))
 fractions=np.abs(defect)/denominators
 return dict(absolute_external_energy_integral_J=math.fsum(abs(float(v)) for v in dw),
  sum_abs_local_energy_defect_J=math.fsum(abs(float(v)) for v in defect),
  max_accepted_local_energy_fraction=float(np.max(fractions)) if len(fractions) else 0.)

def _energy_guard_ledger(meta,history,reference_J,time_ps):
 keys=("absolute_external_energy_integral_J","sum_abs_local_energy_defect_J",
  "local_energy_rejections","max_accepted_local_energy_fraction",
  "max_guarded_accepted_local_energy_fraction","guard_active_since_time_ps",
  "local_energy_reference_energy_J","energy_guard_ledger_origin")
 if meta is not None and any(key in meta for key in keys) and not all(key in meta for key in keys):
  raise ValueError("Partial saved energy-guard ledger")
 if meta is not None and all(key in meta for key in keys):
  result={key:meta[key] for key in keys}
  for key in keys[:-1]:
   value=float(result[key])
   if not math.isfinite(value) or value<0:raise ValueError("Invalid saved energy-guard ledger: "+key)
  if int(result["local_energy_rejections"])!=result["local_energy_rejections"]:
   raise ValueError("Invalid saved local-energy rejection count")
  if not math.isclose(result["local_energy_reference_energy_J"],reference_J,rel_tol=1e-13,abs_tol=0.):
   raise ValueError("Local-energy reference changed")
  if result["guard_active_since_time_ps"]>time_ps:
   raise ValueError("Energy guard starts after checkpoint time")
  return result
 result=_energy_guard_from_history(history,reference_J)
 result.update(local_energy_rejections=0,max_guarded_accepted_local_energy_fraction=0.,
  guard_active_since_time_ps=float(time_ps),local_energy_reference_energy_J=float(reference_J),
  energy_guard_ledger_origin=("reconstructed_from_coherent_pre_guard_history; accepted_midpoint_receipts_after_anchor"
   if meta is not None else "accepted_midpoint_receipts_since_time_zero"))
 return result

def validate_args(a):
 for n in ("length_nm","spacing_nm","xi_cutoff","sigma_nm","dt_ps","ceiling_ps","threshold_mV",
  "tau_ee_Tc_ps","tau_ep_Tc_ps","current_uA","bath_K","snapshot_every_ps","checkpoint_every_ps","checkpoint_wall_s","step_tolerance","local_energy_tolerance"):
  v=float(getattr(a,n))
  if not np.isfinite(v) or v<=0:raise ValueError(n+" must be finite positive")
 if not 0<float(a.retained_fraction)<=1:raise ValueError("retained_fraction must lie in (0,1]")
 if a.wavelength_nm not in (775.,1550.):raise ValueError("wavelength_nm must be775 or1550")
 if not np.isfinite(a.confirmation_ps) or a.confirmation_ps<0 or not np.isfinite(a.extra_time_ps) or a.extra_time_ps<0:
  raise ValueError("Nonnegative confirmation and extra time required")
 if np.isnan(a.tau_escape_ps) or a.tau_escape_ps<=0:raise ValueError("Positive escape time or infinity required")
 if a.duration_ps is not None and (not np.isfinite(a.duration_ps) or a.duration_ps<=0):raise ValueError("Positive development duration required")
 if a.length_nm>5000:raise ValueError("Resolved length exceeds5um wire")
 if a.threads<1 or a.cpu_start<0:raise ValueError("Positive thread count and nonnegative physical-core start required")
 if a.xi_order<2 or a.phonon_nodes<24 or a.transport_order<2 or a.max_rejections<1 or a.max_iterations<1:raise ValueError("Invalid quadrature/solver counts")
 return a
def parameter_identity(a):return safe_json({k:v for k,v in vars(a).items() if k not in ("output_dir","resume","cpu_start")})
def bose_matrix(omega,Eu,temp):
 temp=np.asarray(temp,float);x=omega[None,:]*Eu/(kB*temp[:,None])
 return np.exp(-x)/(-np.expm1(-x))
def prepare_photon(f,k,d,p0,g,*,wavelength_nm,sigma_nm,retained_fraction,bath_temperature_K):
 """Deposit once on the NEW native BCS/phonon states; leave DC gap unchanged."""
 bath=np.full(f.graph.n_nodes,bath_temperature_K);n0=bose_matrix(k.phonons.energies,f.Eu_J,bath)
 xy=g.coordinates_m;center=np.array([.5*(xy[:,0].min()+xy[:,0].max()),0.])
 shape=np.exp(-np.sum((xy-center)**2,axis=1)/(2*(sigma_nm*1e-9)**2))
 cuts=np.r_[g.left_nodes,g.right_nodes];shape[cuts]=0.
 weighted=shape*f.volumes_m3
 if weighted.sum()<=0:raise ValueError("Photon profile lacks interior support")
 photon=Planck*speed_of_light/(wavelength_nm*1e-9);retained=photon*retained_fraction
 target=retained*weighted/weighted.sum();en=f.energies_bar(d)
 ew=4*f.node_energy_scale_J[:,None]*f.weights_bar*en;pw=k.phonon_energy_weights_J
 def values(temp):
  p=f.fermi_dirac(d,temp);n=bose_matrix(k.phonons.energies,f.Eu_J,temp)
  return np.sum(ew*(p-p0),axis=1)+np.sum(pw*(n-n0),axis=1),p,n
 lo=bath.copy();hi=np.maximum(2*bath,2.);active=target>0
 for _ in range(32):
  added,_,_=values(hi);need=(added<target)&active
  if not np.any(need):break
  hi[need]*=2
 else:raise ArithmeticError("No initial FD/BE energy bracket")
 hi[~active]=bath[~active]
 for _ in range(60):
  mid=(lo+hi)/2;added,_,_=values(mid);below=added<target
  lo=np.where(below,mid,lo);hi=np.where(below,hi,mid)
 temp=np.where(active,(lo+hi)/2,bath);added,p,n=values(temp);defect=added-target
 receipt=dict(wavelength_nm=wavelength_nm,photon_energy_J=photon,retained_fraction=retained_fraction,
  retained_energy_J=retained,sigma_nm=sigma_nm,center_m=center.tolist(),
  initialization="Common FD/BE transfer on NEW native BCS/phonon quadratures; fresh DC gap unchanged",
  time_origin="t=0 is effective post-cascade transfer; no additional cascade delay simulated",
  cascade_reference="Exploratory sigma=sqrt(2*D*0.57ps), freely varied; not a measured775/1550nm width",
  energy_deposited_J=float(added.sum()),energy_defect_J=float(defect.sum()),
  maximum_nodal_energy_defect_J=float(max(abs(defect))),relative_total_energy_defect=float(abs(defect.sum())/retained),
  maximum_initial_temperature_K=float(temp.max()),cut_deposition_J=float(target[cuts].sum()),
  source_must_not_be_applied_again=True,old_prepared_photon_state_reused=False)
 if abs(defect.sum())>1e-11*retained:raise ArithmeticError("Initial retained photon energy not conserved")
 return p,n,temp,target,receipt
def phonon_equivalent_temperature_K(k,n):
 target=np.sum(k.phonons.capacities*k.phonons.energies*n,axis=1)
 lo=np.zeros(len(n));hi=np.full(len(n),max(k.bath_temperature_K,1.))
 weights=k.phonons.capacities*k.phonons.energies
 def moment(temp):return np.sum(weights*bose_matrix(k.phonons.energies,k.Eu_J,temp),axis=1)
 for _ in range(32):
  below=moment(hi)<target
  if not np.any(below):break
  hi[below]*=2
 else:raise ArithmeticError("No finite phonon temperature bracket")
 for _ in range(52):
  mid=(lo+hi)/2;below=moment(mid)<target
  lo=np.where(below,mid,lo);hi=np.where(below,hi,mid)
 return np.where(target==0,0.,(lo+hi)/2)


def build_system(a,output,*,resume=False):
 Tc=8.65;D=5e-5;Rs=608.;thickness=7e-9;Eu=kB*Tc
 N0=1/(2*e**2*D*Rs*thickness)
 g=get_effective_geometry(a.length_nm*1e-9,a.spacing_nm*1e-9,output_root=output/"meshes")
 xi,w=bcs_count_quadrature(cutoff_bar=a.xi_cutoff,order=a.xi_order)
 f=EffectiveBCSFunctional(g.graph,energy_scale_J=Eu,N0_per_J_m3=N0,
  ell0_m=g.ell0_m,thickness_m=thickness,count_nodes=xi,count_weights=w,
  gradient_law="bcs_radial",bath_temperature_K=a.bath_K,diffusion_m2_s=D)
 dcpath=output/"dc/effective_dc.npz"
 if resume:
  receipt=json.loads((output/"dc/effective_dc_receipt.json").read_text())
  if sha256(dcpath)!=receipt["fields_sha256"]:raise ValueError("Saved DC fields changed")
  if receipt["geometry"]["mesh_sha256"]!=g.receipt["mesh_sha256"]:raise ValueError("Resume mesh changed")
  with np.load(dcpath,allow_pickle=False) as d:
   dc=EffectiveBCSDC(d["delta_bar"].copy(),d["p"].copy(),d["phi_V"].copy(),
    EffectiveBCSBulk(**receipt["bulk"]),float(receipt["phase_left_rad"]),receipt)
 else:
  dc=prepare_effective_dc(f,left_nodes=g.left_nodes,right_nodes=g.right_nodes,temperature_K=a.bath_K,
   current_A=a.current_uA*1e-6,width_m=g.width_m,sheet_resistance_ohm=Rs)
  save_effective_dc(output/"dc",dc,g)
 tau0=720*float(zeta(5))/np.pi**2*a.tau_ep_Tc_ps*1e-12
 material=allmaras_acoustic_material(Tc_K=Tc,tau0_s=tau0,diffusivity_m2_s=D,
  sheet_resistance_ohm=Rs,thickness_m=thickness,nodes=1025)
 if not resume:material.save(output/"material.npz",output/"material.json")
 phonons=acoustic_phonon_grid(material,Eu,nodes=a.phonon_nodes)
 k=EffectiveBCSKinetics(xi_bar=xi,weights_bar=w,phonons=phonons,Eu_J=Eu,
  N0_per_J_m3=N0,D_m2_s=D,volumes_m3=f.volumes_m3,edges=g.graph.edges,
  face_areas_m2=g.face_areas_m2(thickness),edge_lengths_m=g.edge_lengths_m,
  bath_temperature_K=a.bath_K,tau_ee_Tc_s=a.tau_ee_Tc_ps*1e-12,
  tau_escape_s=a.tau_escape_ps*1e-12,transport_order=a.transport_order,numba_threads=a.threads)
 ports=EffectiveNeutralPorts(g.graph,sheet_resistance_ohm=Rs,
  left_nodes=g.left_nodes,right_nodes=g.right_nodes)
 Lperlength=hbar/(2*e)*dc.bulk.q_derivative_m_inv_per_A
 ext=freeze_exterior_inductance(Lperlength,g.length_m)
 circuit=ThesisCircuitParameters(Lk_ext_H=ext.total_H,V_bias_V=1e4*a.current_uA*1e-6)
 system=EffectiveCoupledSystem(f,k,ports,circuit,tau_ee_Tc_ps=a.tau_ee_Tc_ps,
  tau_ep_Tc_ps=a.tau_ep_Tc_ps,current_scale_A=a.current_uA*1e-6,
  kinetic_residual_energy_scale_J=None if a.dc_only else
   Planck*speed_of_light/(a.wavelength_nm*1e-9)*a.retained_fraction)
 params=dict(Tc_K=Tc,bath_temperature_K=a.bath_K,diffusion_m2_s=D,
  sheet_resistance_ohm=Rs,film_thickness_m=thickness,N0_single_spin_per_J_m3=N0,
  energy_scale_J=Eu,current_A=a.current_uA*1e-6,
  tau_ee_Tc_s=a.tau_ee_Tc_ps*1e-12,tau_ep_Tc_s=a.tau_ep_Tc_ps*1e-12,
  tau0_s=tau0,tau_escape_s=a.tau_escape_ps*1e-12,
  calibrators="Allmaras2020 nomenclature; freely adjustable, not measured NbN constants",
  exterior_inductance=asdict(ext),circuit=asdict(circuit))
 return system,dc,g,material,params


def _session_eta_s(elapsed_session_wall_s,session_start_time_ps,time_ps,remaining_ps):
 # Absolute photon time includes earlier sessions; only newly simulated time
 # measures throughput of this process. No old wall duration is inferred.
 advance_ps=time_ps-session_start_time_ps
 return None if advance_ps<=0 else elapsed_session_wall_s/advance_ps*remaining_ps

def _session_wall_info(elapsed_session_wall_s,session_start_time_ps):
 return dict(elapsed_wall_s=elapsed_session_wall_s,
  elapsed_session_wall_s=elapsed_session_wall_s,
  session_start_time_ps=session_start_time_ps,
  elapsed_wall_scope="current_process_session; elapsed_wall_s is a compatibility alias")

def history_row(system,state,rhs,t_ps,E0,external,limited,totalabs,h_ps):
 Ib,Is,vc=state.circuit_SI
 return dict(time_s=t_ps*1e-12,Vout_V=system.circuit.R_load_ohm*(Ib-Is),
  Vdev_V=rhs.electrical.Vdev_V,Ib_A=Ib,Is_A=Is,vc_V=vc,
  energy_total_J=rhs.energy_J,energy_defect_J=rhs.energy_J-E0-external,
  escape_power_W=rhs.escape_power_W,domain_energy_J=rhs.domain_energy_J,
  instant_power_residual_W=rhs.instant_power_residual_W,
  reservoir_power_out_W=rhs.reservoir_power_out_W,external_energy_integral_J=external,
  limited_event_energy_J=limited,totalabs_event_energy_J=totalabs,accepted_step_ps=h_ps)
def append_history(history,row):
 for key in HISTORY_KEYS:history[key].append(float(row[key]))
def append_snapshot(snapshots,system,state,rhs,t_ps,reference_node_J):
 f,k=system.functional,system.kinetics
 phnode=np.sum(k.phonon_energy_weights_J*state.n,axis=1)
 en=f.energies_bar(state.delta_bar)
 electron=4*f.N0_per_J_m3*f.Eu_J**2*np.sum(f.weights_bar*en*state.p,axis=1)
 values=dict(time_s=t_ps*1e-12,delta_J=state.delta_bar*f.Eu_J,Te_K=rhs.temperature_K,
  Tph_K=phonon_equivalent_temperature_K(k,state.n),
  u_excess_J_m3=(rhs.evaluation.node_energy_J+phnode-reference_node_J)/f.volumes_m3,
  electron_energy_J_m3=electron,phonon_energy_J_m3=phnode/f.volumes_m3)
 for key in SNAPSHOT_KEYS:snapshots[key].append(np.array(values[key],copy=True))
def observation_manifest(controller,a,status):
 result={} if controller is None else controller.summary()
 result.update(status=status,trigger_confirmed=bool(result.get("detected",False)),
  first_crossing_time_s=result.get("crossing_time_s"),
  trigger_time_s=result.get("confirmation_time_s"),baseline_Vout_V=0.,
  threshold_voltage_V=a.threshold_mV*1e-3,confirmation_time_s_required=a.confirmation_ps*1e-12,
  extra_time_s=a.extra_time_ps*1e-12,no_trigger_ceiling_s=a.ceiling_ps*1e-12,
  development_duration_s=None if a.duration_ps is None else a.duration_ps*1e-12,
  full_reset_requested=False,mode="DC" if a.dc_only else "photon",
  observable="Vout=Rload*(Ib-Is), positive excursion from fresh DC baseline")
 return result
def save_checkpoint(output,state,metadata,history,snapshots,g):
 arrays=dict(delta_bar=state.delta_bar,p=state.p,n=state.n,circuit_SI=state.circuit_SI,
  metadata_json=np.array(json.dumps(safe_json(metadata),allow_nan=False)))
 arrays.update({"trajectory__"+key:np.asarray(value) for key,value in history.items()})
 arrays.update({"snapshots__"+key:np.asarray(value) for key,value in snapshots.items()})
 atomic_npz(output/"checkpoint_latest.npz",arrays)
 atomic_npz(output/"trajectory.npz",{key:np.asarray(v) for key,v in history.items()})
 snapdata={key:np.asarray(v) for key,v in snapshots.items()}
 snapdata.update(coordinates_m=g.coordinates_m,triangles=g.triangles,edges=g.graph.edges)
 atomic_npz(output/"snapshots.npz",snapdata)


def run(a, *, stop_requested=None):
 validate_args(a);output=scratch_path(a.output_dir);pins=source_pins();identity=parameter_identity(a)
 if a.resume:
  if not output.is_dir():raise ValueError("Resume directory missing")
  prior=json.loads((output/"manifest.json").read_text())
  if prior.get("run_complete"):raise ValueError("Completed run is preserved; choose a new directory")
  if prior["source_pins_sha256"]!=pins:raise ValueError("Source changed; exact resume refused")
  if prior["parameter_identity"]!=identity:raise ValueError("Parameters changed; exact resume refused")
 else:
  if output.exists():raise FileExistsError("Existing data preserved; use --resume or new output directory")
  output.mkdir(parents=True)
 from sandbox.stage4_core.parallel_runtime import linux_resources,resource_budget
 resources=linux_resources();budget=resource_budget(resources,.9,.9,max_workers=1)
 original_affinity=set(os.sched_getaffinity(0))
 allowed=set(budget["selected_affinity_cpus"])
 admitted_cores=[group for group in resources["physical_core_groups"] if set(group)<=allowed]
 selected_cores=admitted_cores[a.cpu_start:a.cpu_start+a.threads]
 if len(selected_cores)!=a.threads:raise ValueError("Requested distinct physical cores exceed admitted90% budget")
 cpus=[group[0] for group in selected_cores]
 os.sched_setaffinity(0,set(cpus));started=time.monotonic()
 status="PREPARATION";manifest=None;state=None;committed=None
 try:
  system,dc,g,material,params=build_system(a,output,resume=a.resume)
  f,k=system.functional,system.kinetics
  dc_n=np.broadcast_to(k.bath_phonons,(g.graph.n_nodes,k.phonon_size)).copy()
  reference_node=f.evaluate(dc.delta_bar,dc.p).node_energy_J+np.sum(k.phonon_energy_weights_J*dc_n,axis=1)
  history={key:[] for key in HISTORY_KEYS};snapshots={key:[] for key in SNAPSHOT_KEYS}
  if a.resume:
   with np.load(output/"checkpoint_latest.npz",allow_pickle=False) as saved:
    meta=json.loads(str(saved["metadata_json"].item()))
    if meta["source_pins_sha256"]!=pins or meta["parameter_identity"]!=identity:
     raise ValueError("Checkpoint source or parameters differ")
    state=EffectiveState(*(saved[key].copy() for key in ("delta_bar","p","n","circuit_SI")))
    for key in HISTORY_KEYS:history[key]=list(saved["trajectory__"+key])
    for key in SNAPSHOT_KEYS:snapshots[key]=list(saved["snapshots__"+key])
   source=meta["source"];t_ps=float(meta["time_ps"]);E0=float(meta["initial_total_energy_J"])
   external=float(meta["external_energy_integral_J"]);limited=float(meta["limited_event_energy_J"])
   totalabs=float(meta["totalabs_event_energy_J"]);steps=int(meta["accepted_steps"])
   next_h=float(meta["next_step_ps"]);rejections=int(meta["total_rejections"])
   stable=int(meta.get("stable_accepted_steps",0));consecutive=int(meta.get("consecutive_rejections",0))
   controller=None if meta["controller"] is None else TriggerWindow.from_checkpoint(meta["controller"])
   if controller is not None and controller.status=="NUMERICAL_FAILURE":
    raise ValueError("Numerical-failure state is retained for diagnosis; exact continuation requires a valid controller")
  else:
   state=EffectiveState(dc.delta_bar.copy(),dc.p.copy(),dc_n,
    np.array([a.current_uA*1e-6,a.current_uA*1e-6,0.]))
   source=dict(mode="DC",retained_energy_J=0.,energy_deposited_J=0.,source_must_not_be_applied_again=True)
   if not a.dc_only:
    state.p,state.n,initial_temp,target,source=prepare_photon(f,k,state.delta_bar,dc.p,g,
     wavelength_nm=a.wavelength_nm,sigma_nm=a.sigma_nm,retained_fraction=a.retained_fraction,
     bath_temperature_K=a.bath_K)
    atomic_npz(output/"source_initial_state.npz",dict(delta_bar=state.delta_bar,p=state.p,n=state.n,
     temperature_K=initial_temp,node_deposition_J=target,coordinates_m=g.coordinates_m))
   t_ps=0.;external=0.;limited=0.;totalabs=0.;steps=0;next_h=a.dt_ps;rejections=0;stable=0;consecutive=0
   controller=None if a.dc_only else TriggerWindow(TriggerWindowConfig(
    event_time_s=0.,baseline_V=0.,baseline_reference="fresh effective DC Ib=Is, vc=0",
    threshold_V=a.threshold_mV*1e-3,polarity="positive",confirmation_s=a.confirmation_ps*1e-12,
    hysteresis_fraction=.05,post_trigger_margin_s=a.extra_time_ps*1e-12,
    no_trigger_time_ceiling_s=a.ceiling_ps*1e-12))
   rhs=system.evaluate(state);E0=rhs.energy_J
   append_history(history,history_row(system,state,rhs,t_ps,E0,external,limited,totalabs,0.))
   append_snapshot(snapshots,system,state,rhs,t_ps,reference_node)
   if controller is not None:controller.observe_accepted(0.,history["Vout_V"][-1])
  session_start_time_ps=t_ps
  rhs=system.evaluate(state)
  guard_reference=float(system.kinetic_residual_energy_scale_J)
  if not a.dc_only and not math.isclose(guard_reference,float(source["retained_energy_J"]),rel_tol=1e-13,abs_tol=0.):
   raise ValueError("Photon local-energy reference must equal the once-retained source energy")
  guard=_energy_guard_ledger(meta if a.resume else None,history,guard_reference,t_ps)
  engine_threads=k.numba_threads if k.accelerated else 1
  manifest=dict(schema=SCHEMA,run_complete=False,status="RUNNING",
   parameter_identity=identity,source_pins_sha256=pins,physicalparams=params,geometry=g.receipt,
   source=source,material=material.manifest,kinetics=k.metadata(),
   closure=dict(name="effective neutral native BCS/KWT/kinetic/three-state CM",
    gap_energy="One local BCS plus radial-rigidity graph functional, shared force/current/energy",
    gradient_rigidity="BCS radial coefficient at fixed bath temperature; mesoscopic hypothesis",
    electrical="Quasineutral equipotential full-width cuts; fixed normal sheet conductance",
    cuts="Rigid superconducting fresh-DC amplitude; D_tDelta=0; Josephson phase; thermal reservoirs",
    kinetic_energy="Native quasiparticle count with full acoustic Debye phonon band",
    microscopic_exact_Usadel=False,calibrated_experimental_model=False,
    exterior_inductance="Frozen once from continuum fresh-DC dq/dI;96nH plus unresolved5um wire"),
   admission=dict(stage51=False,stage52=False,stage53=False,stage5=False,
    policy="Execution completion does not imply physical/numerical admission"),
   resources=dict(runtime_inventory=resources,scheduling_budget=budget,actual_processes=1,
    actual_threads=engine_threads,requested_threads=a.threads,reserved_physical_core_count=len(selected_cores),
    actual_affinity_cpus=cpus,selected_physical_core_groups=selected_cores,
    threading=("One Numba thread per distinct physical core; BLAS/OpenMP remain1" if k.accelerated
     else "Serial NumPy fallback; requested cores remain reserved, BLAS/OpenMP remain1"),global_resource_ceiling_fraction=.9),
   session_timing=dict(start_time_ps=session_start_time_ps,wall_duration_scope="Current process session only; earlier durations were not saved or accumulated"),
   local_energy_control=guard.copy(),
   numerical=dict(method="Coupled implicit midpoint with shared conservative event extent limits",
    requested_step_ps=a.dt_ps,step_tolerance=a.step_tolerance,max_iterations=a.max_iterations,
    convergence_norm="Gap/CM maximum scaled residual; unsigned electron/phonon residual energy divided by retained photon energy (DC: sum N0 Eu^2 volume); induced gap-force correction guarded at the same tolerance",
    kinetic_reference_energy_J=system.kinetic_residual_energy_scale_J,
    local_energy_tolerance=a.local_energy_tolerance,
    local_energy_guard=dict(enabled=True,reference_energy_J=guard_reference,
     reference_kind="sum N0 Eu^2 nodal volume" if a.dc_only else "retained photon energy",
     acceptance="abs(Enew-Eold-dWmid)<=local_energy_tolerance*(reference+sum(abs(accepted dWmid))+abs(trial dWmid))",
     cumulative_global_energy_audit="Independent; local acceptance does not establish global admission",
     pre_guard_history="Measured accepted increments only; never retroactively certified"),
    no_population_clipping=True,no_raw_forward_euler_time_restriction=True,
    checkpoint_policy="Atomic latest full state plus complete trajectory and snapshots, in scratch",
    resume_requires_identical_sources_parameters_DC_and_mesh=True),
   energy_reference=dict(initial_total_energy_J=E0,
    defect="Efilm(t)+ECM(t)-Efilm(0)-ECM(0)-integral(Psource-Pbias-Pload-Pescape-Pcuts)dt",
    u_excess="Local full BCS functional plus phonon energy minus fresh DC, divided by nodal volume",
    Te_Tph="Energy-equivalent temperatures of evolving nonequilibrium p and n; not forced thermal states"),
   observation=observation_manifest(controller,a,"RUNNING"))
  duration=a.duration_ps if a.duration_ps is not None else (10. if a.dc_only else np.inf)
  next_snapshot=(math.floor(t_ps/a.snapshot_every_ps+1e-8)+1)*a.snapshot_every_ps
  next_checkpoint=(math.floor(t_ps/a.checkpoint_every_ps+1e-8)+1)*a.checkpoint_every_ps
  last_checkpoint=time.monotonic();last_progress=last_checkpoint
  def metadata():
   return dict(schema=SCHEMA,time_ps=t_ps,initial_total_energy_J=E0,
    external_energy_integral_J=external,limited_event_energy_J=limited,
    totalabs_event_energy_J=totalabs,accepted_steps=steps,total_rejections=rejections,
    next_step_ps=next_h,stable_accepted_steps=stable,consecutive_rejections=consecutive,
    source=source,source_pins_sha256=pins,parameter_identity=identity,
    controller=None if controller is None else controller.checkpoint(),
    session_start_time_ps=session_start_time_ps,
    session_wall_time_scope="Descriptive current process session only; earlier durations not accumulated",**guard)
  def capture_committed():
   # State arrays are immutable between accepted assignments. Histories append
   # only: retain their coherent lengths, not another heavy state copy.
   return dict(state=state,rhs=rhs,time_ps=t_ps,external=external,limited=limited,
    totalabs=totalabs,steps=steps,next_h=next_h,rejections=rejections,
    stable=stable,consecutive=consecutive,guard=guard.copy(),controller=None if controller is None else controller.checkpoint(),
    history_size=len(history["time_s"]),snapshot_size=len(snapshots["time_s"]))
  def restore_committed(saved):
   nonlocal state,rhs,t_ps,external,limited,totalabs,steps,next_h,rejections,stable,consecutive,controller,guard
   state=saved["state"];rhs=saved["rhs"];t_ps=saved["time_ps"]
   external=saved["external"];limited=saved["limited"];totalabs=saved["totalabs"]
   steps=saved["steps"];next_h=saved["next_h"];rejections=saved["rejections"]
   stable=saved["stable"];consecutive=saved["consecutive"];guard=saved["guard"].copy()
   controller=None if saved["controller"] is None else TriggerWindow.from_checkpoint(saved["controller"])
   for key in HISTORY_KEYS:del history[key][saved["history_size"]:]
   for key in SNAPSHOT_KEYS:del snapshots[key][saved["snapshot_size"]:]
  committed=capture_committed()
  write_json(output/"manifest.json",manifest)
  save_checkpoint(output,state,metadata(),history,snapshots,g)
  print(json.dumps(safe_json(dict(status="STARTED",output_dir=output,time_ps=t_ps,
   nodes=g.graph.n_nodes,electron_states=k.electron_size,phonon_modes=k.phonon_size,
   current_uA=a.current_uA,dt_ps=a.dt_ps,actual_processes=1,actual_threads=engine_threads,
   requested_threads=a.threads,reserved_physical_core_count=len(selected_cores),
   affinity_cpus=cpus,source=source))),flush=True)

  while True:
   if stop_requested is not None and stop_requested():
    raise KeyboardInterrupt("SIGTERM: checkpoint at last complete accepted state")
   if controller is not None and controller.terminal:
    status=controller.status;break
   if np.isfinite(duration) and t_ps>=duration-64*np.finfo(float).eps*max(1.,duration):
    status="DC_HOLD_COMPLETED" if a.dc_only else "DEVELOPMENT_DURATION_REACHED";break
   horizon=duration
   if controller is not None:horizon=min(horizon,controller.next_required_time_s()*1e12)
   h_ps=min(next_h,horizon-t_ps,next_snapshot-t_ps,next_checkpoint-t_ps)
   if h_ps<=0 or t_ps+h_ps==t_ps:raise ArithmeticError("No representable positive accepted time increment")
   try:
    candidate,receipt=system.midpoint_step(state,h_ps,tolerance=a.step_tolerance,max_iterations=a.max_iterations)
    # Audit the endpoint before the accepted-state transaction and reuse its RHS.
    candidate_rhs=system.evaluate(candidate)
    local_energy=_local_energy_step(rhs.energy_J,candidate_rhs.energy_J,
     receipt["external_energy_increment_J"],guard_reference,guard["absolute_external_energy_integral_J"])
    if local_energy["local_energy_fraction"]>a.local_energy_tolerance:
     raise _LocalEnergyStepRejected("Local energy fraction "+str(local_energy["local_energy_fraction"])+
      " exceeds "+str(a.local_energy_tolerance))
   except (ArithmeticError,ValueError,RuntimeError) as error:
    consecutive+=1;rejections+=1;stable=0;next_h=h_ps/2
    if isinstance(error,_LocalEnergyStepRejected):guard["local_energy_rejections"]+=1
    committed=capture_committed()
    minimum_adaptive_step=getattr(a,"minimum_adaptive_step_ps",0.)
    if minimum_adaptive_step and next_h<minimum_adaptive_step:
     raise ArithmeticError("Declared adaptive-step resolution reached: time_ps="+str(t_ps)+
      ", next_step_ps="+str(next_h)+", minimum_adaptive_step_ps="+str(minimum_adaptive_step)+
      "; last rejected step: "+str(error)) from error
    if consecutive>=a.max_rejections:raise ArithmeticError("Repeated coupled step rejection: "+str(error)) from error
    if time.monotonic()-last_progress>=60:
     print(json.dumps(dict(status="STEP_REJECTED",time_ps=t_ps,next_step_ps=next_h,
      consecutive_rejections=consecutive,error=str(error),**_session_wall_info(time.monotonic()-started,session_start_time_ps))),flush=True)
     last_progress=time.monotonic()
    continue
   state=candidate;t_ps+=h_ps;steps+=1;consecutive=0;stable+=1
   external+=float(receipt["external_energy_increment_J"])
   limited+=float(receipt.get("limited_event_energy_J",0.))
   totalabs+=float(receipt.get("totalabs_event_energy_J",0.))
   guard["absolute_external_energy_integral_J"]+=abs(float(receipt["external_energy_increment_J"]))
   guard["sum_abs_local_energy_defect_J"]+=abs(local_energy["local_energy_defect_J"])
   guard["max_accepted_local_energy_fraction"]=max(guard["max_accepted_local_energy_fraction"],local_energy["local_energy_fraction"])
   guard["max_guarded_accepted_local_energy_fraction"]=max(guard["max_guarded_accepted_local_energy_fraction"],local_energy["local_energy_fraction"])
   rhs=candidate_rhs
   append_history(history,history_row(system,state,rhs,t_ps,E0,external,limited,totalabs,h_ps))
   if controller is not None:controller.observe_accepted(t_ps*1e-12,history["Vout_V"][-1])
   if stable>=8:next_h=min(a.dt_ps,max(next_h,h_ps)*2);stable=0
   committed=capture_committed()
   eps=64*np.finfo(float).eps*max(1.,t_ps)
   if t_ps>=next_snapshot-eps:
    append_snapshot(snapshots,system,state,rhs,t_ps,reference_node);next_snapshot+=a.snapshot_every_ps
   now=time.monotonic()
   if t_ps>=next_checkpoint-eps or now-last_checkpoint>=a.checkpoint_wall_s:
    save_checkpoint(output,state,metadata(),history,snapshots,g);last_checkpoint=now
    if t_ps>=next_checkpoint-eps:next_checkpoint+=a.checkpoint_every_ps
    manifest.update(local_energy_control=guard.copy(),observation=observation_manifest(controller,a,"RUNNING"),
     progress=dict(time_ps=t_ps,accepted_steps=steps,total_rejections=rejections,
      **_session_wall_info(now-started,session_start_time_ps),energy_defect_J=history["energy_defect_J"][-1]))
    write_json(output/"manifest.json",manifest)
   if now-last_progress>=60:
    remaining=max(0.,horizon-t_ps);elapsed=now-started
    print(json.dumps(safe_json(dict(status="RUNNING",time_ps=t_ps,accepted_steps=steps,
     **_session_wall_info(elapsed,session_start_time_ps),remaining_to_current_observation_target_ps=remaining,
     eta_to_current_target_s=_session_eta_s(elapsed,session_start_time_ps,t_ps,remaining),
     step_ps=h_ps,Vout_mV=history["Vout_V"][-1]*1e3,Is_uA=state.circuit_SI[1]*1e6,
     energy_defect_J=history["energy_defect_J"][-1],limited_event_energy_J=limited,
     totalabs_event_energy_J=totalabs,total_rejections=rejections))),flush=True)
    last_progress=now
  if not snapshots["time_s"] or snapshots["time_s"][-1]!=t_ps*1e-12:
   append_snapshot(snapshots,system,state,rhs,t_ps,reference_node)
  save_checkpoint(output,state,metadata(),history,snapshots,g)
  energydef=np.asarray(history["energy_defect_J"]);powerdef=np.asarray(history["instant_power_residual_W"])
  manifest.update(run_complete=True,status=status,local_energy_control=guard.copy(),observation=observation_manifest(controller,a,status),
   execution=dict(accepted_steps=steps,total_rejections=rejections,final_time_s=t_ps*1e-12,
    **_session_wall_info(time.monotonic()-started,session_start_time_ps),maximum_absolute_energy_defect_J=float(max(abs(energydef))),
    maximum_absolute_instant_power_residual_W=float(max(abs(powerdef))),
    limited_event_energy_J=limited,totalabs_event_energy_J=totalabs,
    limited_absolute_event_energy_fraction=limited/totalabs if totalabs else 0.,
    final_Ib_A=float(state.circuit_SI[0]),final_Is_A=float(state.circuit_SI[1]),
    minimum_Ib_A=float(min(history["Ib_A"])),minimum_Is_A=float(min(history["Is_A"])),
    full_reset_simulated=False,**guard),
   output_sha256={name:sha256(output/name) for name in ("trajectory.npz","snapshots.npz",
    "checkpoint_latest.npz","dc/effective_dc.npz","material.npz","material.json")})
  write_json(output/"manifest.json",manifest)
  print(json.dumps(safe_json(dict(status=status,output_dir=output,final_time_ps=t_ps,
   accepted_steps=steps,**_session_wall_info(time.monotonic()-started,session_start_time_ps),admission=False,
   observation=manifest["observation"],energy=manifest["execution"]))),flush=True)
  return manifest
 except BaseException as error:
  if committed is not None:restore_committed(committed)
  if manifest is not None and state is not None and "metadata" in locals():
   status="INTERRUPTED" if isinstance(error,KeyboardInterrupt) else "INTEGRATION_FAILURE"
   if controller is not None and not controller.terminal and status!="INTERRUPTED":
    controller.report_numerical_failure(str(error))
   if snapshots["time_s"][-1]!=t_ps*1e-12:
    append_snapshot(snapshots,system,state,rhs,t_ps,reference_node)
   save_checkpoint(output,state,metadata(),history,snapshots,g)
   manifest.update(status=status,run_complete=False,error=repr(error),local_energy_control=guard.copy(),
    observation=observation_manifest(controller,a,status),
    progress=dict(time_ps=t_ps,accepted_steps=steps,total_rejections=rejections,**_session_wall_info(time.monotonic()-started,session_start_time_ps)))
   write_json(output/"manifest.json",manifest)
  elif output.is_dir():
   write_json(output/"preparation_failure.json",dict(status="PREPARATION_FAILURE",error=repr(error),
    source_pins_sha256=pins,parameter_identity=identity))
  raise
 finally:os.sched_setaffinity(0,original_affinity)

def parser():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument("--output-dir",required=True,type=Path)
 p.add_argument("--resume",action="store_true")
 p.add_argument("--length-nm",type=float,default=160.)
 p.add_argument("--spacing-nm",type=float,default=5.)
 p.add_argument("--xi-cutoff",type=float,default=32.)
 p.add_argument("--xi-order",type=int,default=6)
 p.add_argument("--phonon-nodes",type=int,default=96)
 p.add_argument("--transport-order",type=int,default=48)
 p.add_argument("--wavelength-nm",type=float,choices=(775.,1550.),default=775.)
 p.add_argument("--sigma-nm",type=float,default=7.549834435)
 p.add_argument("--retained-fraction",type=float,default=.667)
 p.add_argument("--dt-ps",type=float,default=.02)
 p.add_argument("--ceiling-ps",type=float,default=200.)
 p.add_argument("--threshold-mV",type=float,default=.1)
 p.add_argument("--confirmation-ps",type=float,default=.5)
 p.add_argument("--extra-time-ps",type=float,default=5.)
 p.add_argument("--tau-ee-Tc-ps",dest="tau_ee_Tc_ps",type=float,default=5.)
 p.add_argument("--tau-ep-Tc-ps",dest="tau_ep_Tc_ps",type=float,default=24.7)
 p.add_argument("--tau-escape-ps",type=float,default=15.)
 p.add_argument("--current-uA",type=float,default=21.29720614489798)
 p.add_argument("--bath-K",dest="bath_K",type=float,default=.9)
 p.add_argument("--dc-only",action="store_true")
 p.add_argument("--duration-ps",type=float,default=None,
  help="Finite development/DC hold; not a latency or stage5 admission")
 p.add_argument("--snapshot-every-ps",type=float,default=.1)
 p.add_argument("--checkpoint-every-ps",type=float,default=1.)
 p.add_argument("--checkpoint-wall-s",type=float,default=60.)
 p.add_argument("--step-tolerance",type=float,default=2e-7)
 p.add_argument("--local-energy-tolerance",type=float,default=1e-5,
  help="Accepted-step energy defect divided by retained energy plus cumulative absolute external work")
 p.add_argument("--max-iterations",type=int,default=40)
 p.add_argument("--max-rejections",type=int,default=50)
 p.add_argument("--threads",type=int,default=_early_threads.threads,help="Numba threads, each pinned to a distinct physical core; BLAS remains1")
 p.add_argument("--cpu-start",type=int,default=0,help="First index of admitted complete physical-core groups; offset independent concurrent runs")
 return p
def main(argv=None):
 # A signal may arrive during any Python assignment or array write. Defer it
 # until the numerical stage and accepted-step ledger transaction are complete.
 requested={"signal":None}
 def terminate(signum,frame):requested["signal"]=int(signum)
 previous=signal.signal(signal.SIGTERM,terminate)
 try:
  return run(parser().parse_args(argv),stop_requested=lambda:requested["signal"] is not None)
 finally:
  signal.signal(signal.SIGTERM,previous)
if __name__=="__main__":main()
