"""Frozen-checkpoint diagnostic: no source injection, no changes to live runs."""
from __future__ import annotations
import argparse,contextlib,hashlib,importlib.util,io,json,math,os,shutil,sys,time
from pathlib import Path
import numpy as np
REPO=Path("/home/jdiaz/pysnspd")
SCRATCH=Path("/home/jdiaz/scratch")
ORIGINALS={"B":SCRATCH/"stage5_effective_L320_h35_775_dt002_v3_20261011",
           "A":SCRATCH/"stage5_effective_L320_h5_775_dt002_v3_20261011"}
def digest(data):return hashlib.sha256(data).hexdigest()
def writej(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False)+"\n")
def freeze(root):
 if root.exists():raise FileExistsError(root)
 root.mkdir(parents=True);raws={};metas={}
 for key,original in ORIGINALS.items():
  raw=(original/"checkpoint_latest.npz").read_bytes()
  with np.load(io.BytesIO(raw),allow_pickle=False) as z:meta=json.loads(str(z["metadata_json"].item()))
  mismatches=[q for q,h in meta["source_pins_sha256"].items() if digest((REPO/q).read_bytes())!=h]
  if mismatches:raise ValueError("Live source/checkpoint mismatch: "+repr(mismatches))
  target=root/("frozen_"+key);target.mkdir()
  (target/"checkpoint_latest.npz").write_bytes(raw)
  for dirname in ("dc","meshes"):
   if (original/dirname).exists():shutil.copytree(original/dirname,target/dirname)
  for filename in ("manifest.json","material.npz","material.json","source_initial_state.npz"):
   if (original/filename).is_file():shutil.copy2(original/filename,target/filename)
  writej(target/"capture_receipt.json",dict(original_directory=str(original),checkpoint_sha256=digest(raw),
    checkpoint_time_ps=meta["time_ps"],captured_unix_s=time.time(),atomic_raw_bytes=True,
    state_keys=["delta_bar","p","n","circuit_SI"],complete_history_controller_and_ledger=True))
  metas[key]=meta;raws[key]=digest(raw)
 frozen=root/"frozen_code"
 shutil.copytree(REPO/"pysnspd",frozen/"pysnspd",ignore=shutil.ignore_patterns("__pycache__","*.pyc"))
 for key in ("pipelines/05_effective_photon_run.py","sandbox/stage4_core/parallel_runtime.py"):
  target=frozen/key;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(REPO/key,target)
 for rel in ("docs/implementation/stage4/practical_time_review_20260924/dual_mesh/resampled",
             "docs/implementation/stage5/prephoton_dc_20260924/meshes/l160_h35",
             "docs/implementation/stage5/prephoton_dc_20260924/meshes/l320_h5"):
  shutil.copytree(REPO/rel,frozen/rel)
 for key,meta in metas.items():
  mismatch=[q for q,h in meta["source_pins_sha256"].items() if digest((frozen/q).read_bytes())!=h]
  if mismatch:raise ValueError("Frozen source mismatch: "+repr(mismatch))
 shutil.copy2(Path(__file__),root/Path(__file__).name)
 writej(root/"experiment_receipt.json",dict(schema="pysnspd.checkpoint_solver_probe.v1",
  originals={k:str(v) for k,v in ORIGINALS.items()},checkpoint_sha256=raws,
  source_pins_sha256=metas["A"]["source_pins_sha256"],B_source_pins_sha256=metas["B"]["source_pins_sha256"],
  helper_sha256=digest(Path(__file__).read_bytes()),affinity=list(os.sched_getaffinity(0)),
  threads=1,BLAS_threads=1,no_live_file_mutation=True,no_photon_reinjection=True,
  no_admission=True,numerical_changes_only=dict(B_tolerances=[2e-9,2e-8,2e-7],B_h_ps=.02,
  B_max_iterations=120,A_hmax_ps=[.02,.05,.1],A_duration_ps=2.,A_tolerance=2e-9)))
 print(json.dumps(dict(status="FROZEN",root=str(root),times_ps={k:m["time_ps"] for k,m in metas.items()})),flush=True)
def load_pipeline(root):
 frozen=root/"frozen_code";sys.path.insert(0,str(frozen))
 spec=importlib.util.spec_from_file_location("frozen_effective_pipeline",frozen/"pipelines/05_effective_photon_run.py")
 mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
def load_state(pipeline,directory):
 with np.load(directory/"checkpoint_latest.npz",allow_pickle=False) as z:
  meta=json.loads(str(z["metadata_json"].item()))
  state=pipeline.EffectiveState(*(z[key].copy() for key in ("delta_bar","p","n","circuit_SI")))
  history={k:list(z["trajectory__"+k]) for k in pipeline.HISTORY_KEYS}
  snapshots={k:list(z["snapshots__"+k]) for k in pipeline.SNAPSHOT_KEYS}
 a=argparse.Namespace(**meta["parameter_identity"]);a.threads=1;a.output_dir=directory;a.resume=True;a.cpu_start=13
 system,dc,g,material,params=pipeline.build_system(a,directory,resume=True)
 if pipeline.source_pins()!=meta["source_pins_sha256"]:raise ValueError("Frozen loaded source differs")
 return system,state,meta,a,history,snapshots,dc,g
def physical(system,state,rhs,meta):
 f=system.functional;Ib,Is,vc=state.circuit_SI;Eret=float(meta["source"]["retained_energy_J"])
 return dict(gap_min_meV=float(np.min(abs(state.delta_bar))*f.Eu_J/1.602176634e-22),
  gap_max_meV=float(np.max(abs(state.delta_bar))*f.Eu_J/1.602176634e-22),
  Vout_V=float(system.circuit.R_load_ohm*(Ib-Is)),Vdev_V=float(rhs.electrical.Vdev_V),
  Ib_A=float(Ib),Is_A=float(Is),vc_V=float(vc),energy_total_J=float(rhs.energy_J),
  domain_energy_J=float(rhs.domain_energy_J),instant_power_residual_W=float(rhs.instant_power_residual_W),
  p_min=float(state.p.min()),p_max=float(state.p.max()),n_min=float(state.n.min()),
  gap_force_max_J_per_gapbar=float(np.max(abs(rhs.evaluation.gap_force_J_per_gapbar))),
  electron_temperature_max_K=float(rhs.temperature_K.max()),retained_energy_J=Eret)
def compare_states(system,ref,other,h_ps):
 f=system.functional;rref=system.evaluate(ref);rother=system.evaluate(other)
 en=f.energies_bar(ref.delta_bar);dp=other.p-ref.p
 popforce=2*f.node_energy_scale_J*ref.delta_bar*np.sum(f.weights_bar*dp/en,axis=1)
 cart=np.column_stack((popforce.real,popforce.imag)).ravel()
 induced=h_ps*(system.mobility_blocks(ref.delta_bar,rref.temperature_K)@cart)
 Eref=system.kinetic_residual_energy_scale_J
 return dict(maximum_complex_gap_difference_bar=float(np.max(abs(other.delta_bar-ref.delta_bar))),
  maximum_gap_magnitude_difference_meV=float(np.max(abs(abs(other.delta_bar)-abs(ref.delta_bar)))*f.Eu_J/1.602176634e-22),
  Vout_difference_V=float(system.circuit.R_load_ohm*((other.circuit_SI[0]-other.circuit_SI[1])-(ref.circuit_SI[0]-ref.circuit_SI[1]))),
  Vdev_difference_V=float(rother.electrical.Vdev_V-rref.electrical.Vdev_V),
  domain_energy_difference_J=float(rother.domain_energy_J-rref.domain_energy_J),
  total_energy_difference_J=float(rother.energy_J-rref.energy_J),
  unsigned_population_difference_energy_J=float(np.sum(4*f.node_energy_scale_J[:,None]*f.weights_bar*en*abs(dp))),
  unsigned_phonon_difference_energy_J=float(np.sum(system.kinetics.phonon_energy_weights_J*abs(other.n-ref.n))),
  maximum_population_difference=float(abs(dp).max()),
  maximum_population_force_difference_J_per_gapbar=float(abs(popforce).max()),
  maximum_induced_gap_difference_bar=float(np.max(abs(induced[system._free_cart]))),
  maximum_total_force_difference_J_per_gapbar=float(np.max(abs(rother.evaluation.gap_force_J_per_gapbar-rref.evaluation.gap_force_J_per_gapbar))),
  normalization_retained_energy_J=Eref)
def run_B(root,pipeline):
 system,initial,meta,a,history,snapshots,dc,g=load_state(pipeline,root/"frozen_B")
 out=root/"B_single_steps";out.mkdir(exist_ok=True);records=[];states={}
 initial_rhs=system.evaluate(initial);initialphysical=physical(system,initial,initial_rhs,meta)
 for tol in (2e-9,2e-8,2e-7):
  tag=format(tol,".0e").replace("-","m");started=time.monotonic();trace=out/(tag+".jsonl")
  if (out/(tag+"_receipt.json")).exists():
   record=json.loads((out/(tag+"_receipt.json")).read_text());records.append(record)
   if record["status"]=="CONVERGED":
    with np.load(out/(tag+"_state.npz"),allow_pickle=False) as z:
     states[tag]=pipeline.EffectiveState(*(z[key].copy() for key in ("delta_bar","p","n","circuit_SI")))
   print(json.dumps(dict(record,event="B_CASE_REUSED_WITHOUT_RECOMPUTE")),flush=True);continue
  record=dict(tolerance=tol,h_ps=.02,max_iterations=120,start_time_ps=meta["time_ps"])
  os.environ["PYSNSPD_MIDPOINT_DEBUG"]="1"
  try:
   with trace.open("w") as stream,contextlib.redirect_stdout(stream):
    candidate,receipt=system.midpoint_step(initial.copy(),.02,tolerance=tol,max_iterations=120)
   rhs=system.evaluate(candidate);states[tag]=candidate
   pipeline.atomic_npz(out/(tag+"_state.npz"),dict(delta_bar=candidate.delta_bar,p=candidate.p,n=candidate.n,circuit_SI=candidate.circuit_SI))
   record.update(status="CONVERGED",receipt=receipt,physical=physical(system,candidate,rhs,meta),
    integrated_step_energy_defect_J=float(rhs.energy_J-initial_rhs.energy_J-receipt["external_energy_increment_J"]),
    inherited_plus_step_energy_defect_J=float(rhs.energy_J-meta["initial_total_energy_J"]-meta["external_energy_integral_J"]-receipt["external_energy_increment_J"]))
  except (ArithmeticError,ValueError,RuntimeError) as error:record.update(status="FAILED",error=repr(error))
  finally:os.environ.pop("PYSNSPD_MIDPOINT_DEBUG",None)
  rows=[json.loads(line) for line in trace.read_text().splitlines() if line.startswith("{")]
  record.update(elapsed_wall_s=time.monotonic()-started,trace_iterations=len(rows),last_trace=rows[-1] if rows else None,
   trace_sha256=digest(trace.read_bytes()))
  writej(out/(tag+"_receipt.json"),record);records.append(record)
  print(json.dumps(dict(record,event="B_CASE_FINISHED")),flush=True)
 strict=states.get("2em09");comparisons={}
 if strict is not None:
  for tag,state in states.items():
   if tag!="2em09":comparisons[tag]=compare_states(system,strict,state,.02)
 writej(out/"summary.json",dict(status="DIAGNOSTIC_ONLY",initial=initialphysical,cases=records,comparisons_to_2e9=comparisons,
  no_admission=True,no_RHS_changes=True,no_photon_reinjection=True))
def run_A(root,pipeline,*,tolerance=2e-9,folder="A_continuations",caps=(.02,.05,.1)):
 base=root/"frozen_A";results=[]
 for cap in caps:
  tag=("hmax_"+str(cap).replace(".","p"));out=root/folder/tag
  out.mkdir(parents=True)
  for dirname in ("dc","meshes"):
   if (base/dirname).exists():shutil.copytree(base/dirname,out/dirname)
  system,state,meta,a,history,snapshots,dc,g=load_state(pipeline,base)
  # Reuse the same full state, ledger, source and observation-controller checkpoint.
  # Reset adaptation counters exactly as the production resume entry does.
  controller=None if meta["controller"] is None else pipeline.TriggerWindow.from_checkpoint(meta["controller"])
  start_t=t=float(meta["time_ps"]);target=start_t+2.
  E0=float(meta["initial_total_energy_J"]);external=float(meta["external_energy_integral_J"])
  limited=float(meta["limited_event_energy_J"]);totalabs=float(meta["totalabs_event_energy_J"])
  steps=int(meta["accepted_steps"]);rejections=int(meta["total_rejections"]);next_h=float(meta["next_step_ps"])
  next_snapshot=(math.floor(t/a.snapshot_every_ps+1e-8)+1)*a.snapshot_every_ps
  next_checkpoint=(math.floor(t/a.checkpoint_every_ps+1e-8)+1)*a.checkpoint_every_ps
  dc_n=np.broadcast_to(system.kinetics.bath_phonons,state.n.shape)
  reference_node=system.functional.evaluate(dc.delta_bar,dc.p).node_energy_J+np.sum(system.kinetics.phonon_energy_weights_J*dc_n,axis=1)
  initial_steps=steps;initial_rejections=rejections;initial_external=external;initial_limited=limited
  stable=consecutive=0;started=last_progress=last_checkpoint=time.monotonic();step_receipts=[]
  def metadata():
   return dict(meta,time_ps=t,external_energy_integral_J=external,limited_event_energy_J=limited,
    totalabs_event_energy_J=totalabs,accepted_steps=steps,total_rejections=rejections,next_step_ps=next_h,
    controller=None if controller is None else controller.checkpoint(),
    diagnostic_override=dict(hmax_ps=cap,original_hmax_ps=a.dt_ps,step_tolerance=tolerance,original_step_tolerance=a.step_tolerance,duration_after_checkpoint_ps=2.))
  rhs=system.evaluate(state);status="DIAGNOSTIC_DURATION_REACHED";error=None
  try:
   while t<target-64*np.finfo(float).eps*max(1.,target):
    if controller is not None and controller.terminal:status=controller.status;break
    horizon=min(target,controller.next_required_time_s()*1e12) if controller is not None else target
    h=min(next_h,horizon-t,next_snapshot-t,next_checkpoint-t)
    if h<=0:raise ArithmeticError("Nonpositive diagnostic step")
    try:candidate,receipt=system.midpoint_step(state,h,tolerance=tolerance,max_iterations=a.max_iterations)
    except (ArithmeticError,ValueError,RuntimeError) as exc:
     consecutive+=1;rejections+=1;stable=0;next_h=h/2.
     if consecutive>=a.max_rejections:raise ArithmeticError("Repeated step rejection: "+str(exc)) from exc
     if time.monotonic()-last_progress>=60:
      print(json.dumps(dict(status="A_STEP_REJECTED",hmax_ps=cap,time_ps=t,next_step_ps=next_h,error=str(exc))),flush=True);last_progress=time.monotonic()
     continue
    state=candidate;t+=h;steps+=1;consecutive=0;stable+=1
    external+=float(receipt["external_energy_increment_J"]);limited+=float(receipt.get("limited_event_energy_J",0.))
    totalabs+=float(receipt.get("totalabs_event_energy_J",0.))
    step_receipts.append(dict(time_ps=t,h_ps=h,**receipt))
    rhs=system.evaluate(state)
    pipeline.append_history(history,pipeline.history_row(system,state,rhs,t,E0,external,limited,totalabs,h))
    if controller is not None:controller.observe_accepted(t*1e-12,history["Vout_V"][-1])
    if stable>=8:next_h=min(cap,max(next_h,h)*2);stable=0
    eps=64*np.finfo(float).eps*max(1.,t)
    if t>=next_snapshot-eps:
     pipeline.append_snapshot(snapshots,system,state,rhs,t,reference_node);next_snapshot+=a.snapshot_every_ps
    now=time.monotonic()
    if t>=next_checkpoint-eps or now-last_checkpoint>=60:
     pipeline.save_checkpoint(out,state,metadata(),history,snapshots,g);last_checkpoint=now
     if t>=next_checkpoint-eps:next_checkpoint+=a.checkpoint_every_ps
    if now-last_progress>=60:
     print(json.dumps(dict(status="A_RUNNING",hmax_ps=cap,time_ps=t,target_ps=target,new_steps=steps-initial_steps,
      new_rejections=rejections-initial_rejections,step_ps=h,Vout_V=history["Vout_V"][-1],
      elapsed_wall_s=now-started)),flush=True);last_progress=now
  except (ArithmeticError,ValueError,RuntimeError) as exc:status="DIAGNOSTIC_FAILURE";error=repr(exc)
  if not snapshots["time_s"] or snapshots["time_s"][-1]!=t*1e-12:
   pipeline.append_snapshot(snapshots,system,state,rhs,t,reference_node)
  pipeline.save_checkpoint(out,state,metadata(),history,snapshots,g)
  pipeline.write_json(out/"step_receipts.json",step_receipts)
  record=dict(hmax_ps=cap,step_tolerance=tolerance,status=status,error=error,start_time_ps=start_t,final_time_ps=t,
   new_steps=steps-initial_steps,new_rejections=rejections-initial_rejections,elapsed_wall_s=time.monotonic()-started,
   step_range_ps=[min([r["h_ps"] for r in step_receipts],default=0.),max([r["h_ps"] for r in step_receipts],default=0.)],
   external_increment_J=external-initial_external,limited_event_increment_J=limited-initial_limited,
   physical=physical(system,state,rhs,meta),initial_checkpoint_sha256=digest((base/"checkpoint_latest.npz").read_bytes()),
   source_pins_sha256=meta["source_pins_sha256"],controller=None if controller is None else controller.checkpoint(),
   maximum_energy_defect_J=max(abs(np.asarray(history["energy_defect_J"])[len(history["time_s"])-(steps-initial_steps):]),default=0.))
  pipeline.write_json(out/"receipt.json",record);results.append(record)
  print(json.dumps(dict(record,event="A_CASE_FINISHED")),flush=True)
 writej(root/folder/"summary.json",dict(cases=results,status="DIAGNOSTIC_ONLY",no_admission=True))
def main():
 p=argparse.ArgumentParser();p.add_argument("--output",type=Path,required=True);p.add_argument("--phase",choices=("freeze","B","A","practical","all"),default="all")
 a=p.parse_args();root=a.output.resolve()
 if not root.is_relative_to(SCRATCH):raise ValueError("scratch required")
 if a.phase in ("freeze","all"):freeze(root)
 if a.phase!="freeze":
  pipeline=load_pipeline(root)
  if a.phase in ("B","all"):run_B(root,pipeline)
  if a.phase in ("A","all"):run_A(root,pipeline)
  if a.phase=="practical":run_A(root,pipeline,tolerance=2e-7,folder="A_practical_norm_continuations")
if __name__=="__main__":main()
