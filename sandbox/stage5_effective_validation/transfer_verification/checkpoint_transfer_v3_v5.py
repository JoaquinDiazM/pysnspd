"""Explicit, non-evolving checkpoint transfer with optional v5 energy guarding; never signals jobs.

Commands:
 capture --run <live-or-stopped-v3> --anchor <new-scratch-dir>
 transfer --anchor <immutable-anchor> --output <new-campaign> --dt-ps .05
          --step-tolerance 2e-8 --threads 4 --cpu-start 0
Transfer prepares an exact --resume of the NEW capsule. It does not run it.
Only dt_ps, step_tolerance, explicit local_energy_tolerance and threads may change. An explicit --duration-ps
extension is additionally allowed ONLY after a completed development duration
with an ACTIVE observation controller; terminal events are never reopened. Sources may change only at
the two reviewed frontend paths. The physical clock/source/ledger are immutable.
"""
from __future__ import annotations
import argparse,copy,hashlib,importlib.util,io,json,math,os,shlex,shutil,sys
from pathlib import Path
import numpy as np

SCRATCH=Path("/home/jdiaz/scratch").resolve()
REPO=Path("/home/jdiaz/pysnspd")
FRONTEND={"pipelines/05_effective_photon_run.py","pysnspd/experimental/effective_bcs_preparation.py"}
STATE=("delta_bar","p","n","circuit_SI")
CONSERVED=("time_ps","initial_total_energy_J","external_energy_integral_J",
 "limited_event_energy_J","totalabs_event_energy_J","accepted_steps","total_rejections",
 "source","controller")
LEDGERS=("external_energy_integral_J","limited_event_energy_J","totalabs_event_energy_J")
SCHEMA="pysnspd.explicit_checkpoint_transfer.v2"
PARENT_HELPER_SHA256="0c70590e5f37f877f6fc988348bd45c631bad55de55f55823d9d4e8ed99111cd"
GUARD_KEYS=("absolute_external_energy_integral_J","sum_abs_local_energy_defect_J",
 "local_energy_rejections","max_accepted_local_energy_fraction",
 "max_guarded_accepted_local_energy_fraction","guard_active_since_time_ps",
 "local_energy_reference_energy_J","energy_guard_ledger_origin")

def digest(data):return hashlib.sha256(data).hexdigest()
def filehash(path):return digest(Path(path).read_bytes())
def arrayhash(a):
 a=np.ascontiguousarray(a)
 return digest(json.dumps({"dtype":a.dtype.str,"shape":a.shape},sort_keys=True).encode()+a.tobytes())
def scratch(path):
 p=Path(path).expanduser().resolve()
 if not p.is_relative_to(SCRATCH) or p==SCRATCH:raise ValueError("Path must be below scratch")
 return p
def newdir(path):
 p=scratch(path)
 if p.exists():raise FileExistsError("Original/new directory already exists; refusing overwrite: "+str(p))
 p.mkdir(parents=True);return p
def writejson(path,data):
 p=Path(path);p.parent.mkdir(parents=True,exist_ok=True)
 p.write_text(json.dumps(data,indent=2,allow_nan=False)+"\n")
def decode(payload):
 with np.load(io.BytesIO(payload),allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}
def metadata(arrays):return json.loads(str(arrays["metadata_json"].item()))
def close_scalar(x,y,label):
 # Comparison only; no physical state or ledger is modified.
 if x!=y and abs(x-y)>64*max(math.ulp(float(x)),math.ulp(float(y))):
  raise ValueError("Incoherent checkpoint "+label+": "+repr((x,y)))
def audit_arrays(arrays,manifest):
 m=metadata(arrays)
 for k in ("time_ps","initial_total_energy_J",*LEDGERS,"next_step_ps"):
  if not np.isfinite(m[k]):raise ValueError("Nonfinite metadata "+k)
 if m["time_ps"]<0 or m["next_step_ps"]<=0:raise ValueError("Invalid clock/next step")
 if m["source_pins_sha256"]!=manifest["source_pins_sha256"]:raise ValueError("Checkpoint/manifest pins differ")
 if m["parameter_identity"]!=manifest["parameter_identity"]:raise ValueError("Checkpoint/manifest parameters differ")
 if m["schema"]!=manifest["schema"]:raise ValueError("Checkpoint/manifest schemas differ")
 for k in STATE:
  if k not in arrays or not np.all(np.isfinite(arrays[k])):raise ValueError("Nonfinite/missing state "+k)
 if np.any((arrays["p"]<0)|(arrays["p"]>1)) or np.any(arrays["n"]<0):raise ValueError("Pauli/Bose support violated")
 if arrays["circuit_SI"].shape!=(3,):raise ValueError("Circuit shape")
 if np.any(arrays["circuit_SI"][:2]<0):raise ValueError("Global current polarization reversal")
 nnode=len(arrays["delta_bar"])
 if arrays["delta_bar"].shape!=(nnode,) or arrays["p"].shape[0]!=nnode or arrays["n"].shape[0]!=nnode:
  raise ValueError("Node shapes differ")
 ts=arrays["trajectory__time_s"]
 if ts.ndim!=1 or len(ts)<1 or ts[0]!=0 or np.any(np.diff(ts)<=0):raise ValueError("History clock invalid")
 if len(ts)!=int(m["accepted_steps"])+1:raise ValueError("Accepted counter/history length differ")
 close_scalar(float(ts[-1]),float(m["time_ps"])*1e-12,"physical clock")
 for k,v in arrays.items():
  if k.startswith("trajectory__") and (v.shape!=(len(ts),) or not np.all(np.isfinite(v))):
   raise ValueError("History incomplete/nonfinite: "+k)
 for k in LEDGERS:close_scalar(float(arrays["trajectory__"+k][-1]),float(m[k]),k)
 for k,i in (("Ib_A",0),("Is_A",1),("vc_V",2)):
  close_scalar(float(arrays["trajectory__"+k][-1]),float(arrays["circuit_SI"][i]),k)
 R=manifest["physicalparams"]["circuit"]["R_load_ohm"]
 close_scalar(float(arrays["trajectory__Vout_V"][-1]),R*(arrays["circuit_SI"][0]-arrays["circuit_SI"][1]),"Vout")
 close_scalar(float(arrays["trajectory__energy_defect_J"][-1]),
  float(arrays["trajectory__energy_total_J"][-1])-m["initial_total_energy_J"]-m["external_energy_integral_J"],"energy ledger")
 st=arrays["snapshots__time_s"]
 if st.ndim!=1 or len(st)<1 or st[0]!=0 or np.any(np.diff(st)<=0) or st[-1]>ts[-1]:
  raise ValueError("Snapshot clocks invalid")
 for k,v in arrays.items():
  if k.startswith("snapshots__") and (len(v)!=len(st) or not np.all(np.isfinite(v))):
   raise ValueError("Snapshot incomplete/nonfinite: "+k)
 c=m["controller"]
 if c is not None:
  if c["state"]["status"] in ("NUMERICAL_FAILURE","INVALID_INITIAL_EVENT"):raise ValueError("Invalid controller")
  close_scalar(c["state"]["last_time_s"],float(ts[-1]),"controller clock")
  if int(c["state"]["sample_count"])!=len(ts):raise ValueError("Controller/history count differ")
  cfg=c["config"];sgn=1 if cfg["polarity"]=="positive" else -1
  close_scalar(c["state"]["last_excursion_V"],sgn*(float(arrays["trajectory__Vout_V"][-1])-cfg["baseline_V"]),"controller voltage")
 if m["source"].get("source_must_not_be_applied_again") is not True:raise ValueError("Missing one-time deposition contract")
 return m

def audit_runtime(runtime,pins):
 archive_path=runtime/"archive_manifest.json"
 a=json.loads(archive_path.read_text())
 if a["runtime_pins_sha256"]!=pins:raise ValueError("Original runtime archive pins differ")
 for rel,want in pins.items():
  if filehash(runtime/rel)!=want:raise ValueError("Original runtime pin mismatch "+rel)
 closure=a["package_sources_and_small_geometry"]
 for rel,item in closure.items():
  want=item["sha256"] if isinstance(item,dict) else item
  if filehash(runtime/rel)!=want:raise ValueError("Transitive runtime source mismatch "+rel)
 return {"manifest_sha256":filehash(archive_path),"pin_count":len(pins),"transitive_file_count":len(closure)}

def capture(run,anchor):
 run=scratch(run);anchor=newdir(anchor)
 # One read of the atomic NPZ. This byte string is the unique physical anchor.
 payload=(run/"checkpoint_latest.npz").read_bytes()
 manifest_payload=(run/"manifest.json").read_bytes()
 arrays=decode(payload);manifest=json.loads(manifest_payload)
 m=audit_arrays(arrays,manifest)
 original_runtime=run/"runtime_sources"
 archive=audit_runtime(original_runtime,m["source_pins_sha256"])
 (anchor/"checkpoint_original.npz").write_bytes(payload)
 (anchor/"manifest_original.json").write_bytes(manifest_payload)
 shutil.copytree(original_runtime,anchor/"runtime_sources")
 # trajectory and snapshot arrays come from the SAME checkpoint, not from
 # independently changing output files.
 np.savez_compressed(anchor/"trajectory.npz",**{k[12:]:v for k,v in arrays.items() if k.startswith("trajectory__")})
 np.savez_compressed(anchor/"snapshots_fields.npz",**{k[11:]:v for k,v in arrays.items() if k.startswith("snapshots__")})
 for rel in ("dc","meshes"):
  if (run/rel).exists():shutil.copytree(run/rel,anchor/rel)
 for rel in ("material.npz","material.json","source_initial_state.npz"):
  if not (run/rel).is_file():raise ValueError("Missing immutable physical data "+rel)
  shutil.copy2(run/rel,anchor/rel)
 dc=json.loads((anchor/"dc/effective_dc_receipt.json").read_text())
 if filehash(anchor/"dc/effective_dc.npz")!=dc["fields_sha256"]:raise ValueError("DC data hash differs")
 # Diagnostic parents retain their exact helper/resource override separately
 # from the ancestral runner identity; preserve that evidence without editing NPZ.
 if (run/"parent_preparation_receipt.json").exists():
  shutil.copy2(run/"parent_preparation_receipt.json",anchor/"parent_preparation_receipt.json")
 if (run/"diagnostic_provenance").is_dir():
  shutil.copytree(run/"diagnostic_provenance",anchor/"diagnostic_provenance")
 # Preserve a prior continuation receipt so numerical segments remain rooted
 # in their original clocks instead of rewriting all prior history as one solver.
 prior_receipt=None
 if manifest.get("explicit_transfer_receipt"):
  prior_path=(run/manifest["explicit_transfer_receipt"]).resolve()
  if not prior_path.is_relative_to(run):raise ValueError("Prior transfer receipt escapes parent")
  prior_payload=prior_path.read_bytes();prior_receipt=json.loads(prior_payload)
  if prior_receipt["pins_after"]!=m["source_pins_sha256"]:raise ValueError("Prior transfer receipt pins differ")
  target=anchor/"provenance/parent_transfer_receipt.json";target.parent.mkdir(exist_ok=True)
  target.write_bytes(prior_payload)
 physicalfiles={str(p.relative_to(anchor)):filehash(p) for p in anchor.rglob("*") if p.is_file() and "runtime_sources" not in p.parts}
 receipt=dict(schema=SCHEMA,action="CAPTURE_ONLY_NO_SIGNAL_NO_EVOLUTION",parent_output_dir=str(run),
  checkpoint_bytes_sha256=digest(payload),manifest_sha256=digest(manifest_payload),runtime_archive=archive,
  physical_file_sha256=physicalfiles,state_arrays_sha256={k:arrayhash(arrays[k]) for k in STATE},
  conserved_metadata={k:m[k] for k in CONSERVED},parameter_identity=m["parameter_identity"],
  source_pins_sha256=m["source_pins_sha256"],live_source_checkpoint_may_advance=True,
  immutable_anchor=str(anchor))
 writejson(anchor/"capture_receipt.json",receipt)
 if filehash(anchor/"checkpoint_original.npz")!=digest(payload):raise ValueError("Capture bytes differ")
 return receipt

def load_runner(root):
 path=root/"pipelines/05_effective_photon_run.py"
 spec=importlib.util.spec_from_file_location("transfer_runner",path)
 module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module
 old=sys.argv;sys.argv=["checkpoint_transfer"]
 try:spec.loader.exec_module(module)
 finally:sys.argv=old
 return module

def resume_argv(identity,output,cpu_start):
 result=["--output-dir",str(output),"--resume","--cpu-start",str(cpu_start)]
 for k,v in identity.items():
  if v is None:continue
  flag="--"+k.replace("_","-")
  if isinstance(v,bool):
   if v:result.append(flag)
  else:result.extend((flag,str(v)))
 return result


def duration_argument(value):
 if value.lower()=="unchanged":return "unchanged"
 if value.lower()=="none":return None
 return float(value)

def extend_control_horizon(parent,oldmeta,requested):
 controller=oldmeta["controller"]
 if controller is not None and controller["state"]["status"] not in ("ACTIVE","TRIGGER_CONFIRMED"):
  raise ValueError("Terminal event may never be reopened")
 if requested=="unchanged":
  if parent.get("run_complete"):raise ValueError("Completed control requires an explicit eligible horizon extension")
  return oldmeta["parameter_identity"].get("duration_ps"),None
 if (parent.get("run_complete") is not True or parent.get("status")!="DEVELOPMENT_DURATION_REACHED"
     or controller is None or controller["state"]["status"]!="ACTIVE"):
  raise ValueError("Only completed DEVELOPMENT_DURATION_REACHED with ACTIVE controller may extend duration")
 if any(controller["state"].get(k) is not None for k in ("crossing_time_s","confirmation_time_s","stop_target_s","actual_stop_s")):
  raise ValueError("ACTIVE completed control contains a confirmed or terminal event")
 prior=oldmeta["parameter_identity"].get("duration_ps")
 if prior is None or isinstance(prior,bool) or not np.isfinite(prior) or prior<=0:
  raise ValueError("Completed control lacks a finite original duration")
 close_scalar(float(oldmeta["time_ps"]),float(prior),"completed control duration clock")
 if requested is not None:
  if isinstance(requested,bool) or not np.isfinite(requested):
   raise ValueError("Extended duration must be finite or explicit none")
  requested=float(requested)
  margin=128*np.finfo(float).eps*max(abs(prior),abs(requested))
  if requested<=prior+margin:raise ValueError("Extended duration must exceed the completed duration")
 change=dict(runtime_key="duration_ps",before_ps=prior,after_ps=requested,
  effective_time_ps=oldmeta["time_ps"],parent_status=parent["status"],
  parent_run_complete=True,controller_status_preserved="ACTIVE",
  physical_clock_preserved=True,no_event_reopened=True,
  unchanged_no_trigger_ceiling_ps=oldmeta["parameter_identity"].get("ceiling_ps"))
 return requested,change


def energy_guard_transition(oldmeta,arrays,parent,tolerance):
 """Explicit numerical policy only; never changes state, source or signed work."""
 if tolerance is None:return None,None
 if not math.isfinite(float(tolerance)) or tolerance<=0:raise ValueError("Positive finite local-energy tolerance required")
 oldtol=oldmeta["parameter_identity"].get("local_energy_tolerance")
 if oldtol is not None and (not math.isfinite(float(oldtol)) or oldtol<=0):
  raise ValueError("Invalid parent local-energy tolerance")
 ref=float(oldmeta["source"].get("retained_energy_J",0.))
 kind="retained photon energy"
 if ref<=0:
  ref=float(parent["numerical"]["kinetic_reference_energy_J"])
  kind="sum N0 Eu^2 nodal volume"
 if not math.isfinite(ref) or ref<=0:raise ValueError("Positive finite local-energy reference required")
 energy=np.asarray(arrays["trajectory__energy_total_J"],float)
 external=np.asarray(arrays["trajectory__external_energy_integral_J"],float)
 dw=np.diff(external);defect=np.diff(energy)-dw
 fractions=np.abs(defect)/(ref+np.cumsum(np.abs(dw)))
 measured=dict(absolute_external_energy_integral_J=math.fsum(abs(float(v)) for v in dw),
  sum_abs_local_energy_defect_J=math.fsum(abs(float(v)) for v in defect),
  max_accepted_local_energy_fraction=float(np.max(fractions)) if len(fractions) else 0.)
 now=float(oldmeta["time_ps"])
 if oldtol is None:
  if any(k in oldmeta for k in GUARD_KEYS):raise ValueError("Unguarded parent contains ambiguous guard metadata")
  ledger=dict(measured,local_energy_rejections=0,max_guarded_accepted_local_energy_fraction=0.,
   guard_active_since_time_ps=now,local_energy_reference_energy_J=ref,
   energy_guard_ledger_origin="reconstructed_from_coherent_pre_guard_history; accepted_midpoint_receipts_after_anchor")
 else:
  if not all(k in oldmeta for k in GUARD_KEYS):raise ValueError("Guarded parent lacks complete guard ledger")
  ledger={k:oldmeta[k] for k in GUARD_KEYS}
  for k in GUARD_KEYS[:-1]:
   if not math.isfinite(float(ledger[k])) or ledger[k]<0:raise ValueError("Invalid parent guard ledger")
  if not math.isclose(ledger["local_energy_reference_energy_J"],ref,rel_tol=1e-13,abs_tol=0.):
   raise ValueError("Parent guard reference differs")
  if oldtol!=tolerance:
   ledger.update(guard_active_since_time_ps=now,max_guarded_accepted_local_energy_fraction=0.,
    energy_guard_ledger_origin=str(ledger["energy_guard_ledger_origin"])+"; tolerance_changed_at_explicit_transfer")
 change=dict(enabled_before=oldtol is not None,enabled_after=True,
  local_energy_tolerance_before=oldtol,local_energy_tolerance_after=float(tolerance),
  effective_time_ps=now,reference_energy_J=ref,reference_kind=kind,
  denominator="reference + cumulative sum(abs(accepted external_energy_increment_J)) + abs(trial external_energy_increment_J)",
  pre_guard_history_certified=False,history_diagnostic_measured=measured,
  ledger_after=copy.deepcopy(ledger),
  old_rejection_count_policy="Preserved if recorded under a guard; zero counts only guarded rejections since activation for legacy",
  state_clock_source_signed_ledger_unchanged=True)
 return ledger,change

def numerical_segments_for_transfer(anchor,oldmeta,identity,horizon_change,guard_change):
 now=float(oldmeta["time_ps"])
 prior_path=anchor/"provenance/parent_transfer_receipt.json"
 if prior_path.exists():
  prior=json.loads(prior_path.read_text());segments=copy.deepcopy(prior["numerical_segments"])
  if not segments or segments[-1]["end_time_ps"] is not None:
   raise ValueError("Prior numerical segments lack open continuation")
  if segments[-1]["start_time_ps"]>now:raise ValueError("Prior numerical segment begins after anchor")
  for k in ("dt_ps","step_tolerance","threads"):
   if segments[-1][k]!=oldmeta["parameter_identity"][k]:raise ValueError("Prior numerical segment differs: "+k)
  if segments[-1].get("local_energy_tolerance")!=oldmeta["parameter_identity"].get("local_energy_tolerance"):
   raise ValueError("Prior numerical guard segment differs")
  segments[-1]["end_time_ps"]=now
 else:
  segments=[{"start_time_ps":0.,"end_time_ps":now,"origin":"v3",
   **{k:oldmeta["parameter_identity"][k] for k in ("dt_ps","step_tolerance","threads")}}]
  diagnostic_path=anchor/"parent_preparation_receipt.json"
  if diagnostic_path.exists():
   diagnostic=json.loads(diagnostic_path.read_text())
   if diagnostic["checkpoint_sha256"]!=filehash(anchor/"checkpoint_original.npz"):
    raise ValueError("Diagnostic parent checkpoint provenance differs")
   settings=diagnostic["executed_A_numeric_settings"]
   close_scalar(float(settings["end_time_ps"]),now,"diagnostic endpoint")
   start=float(settings["start_time_ps"])
   if not 0<=start<=now:raise ValueError("Diagnostic segment start outside history")
   if settings["max_step_ps"]!=oldmeta["parameter_identity"]["dt_ps"] or settings["step_tolerance"]!=oldmeta["parameter_identity"]["step_tolerance"]:
    raise ValueError("Diagnostic solver identity differs")
   if settings["max_iterations"]!=oldmeta["parameter_identity"]["max_iterations"]:
    raise ValueError("Diagnostic iteration cap differs")
   segments[0]["end_time_ps"]=start
   segments.append(dict(start_time_ps=start,end_time_ps=now,origin="frozen_diagnostic_A",
    dt_ps=settings["max_step_ps"],step_tolerance=settings["step_tolerance"],
    threads=settings["threads"],actual_affinity_cpus=[settings["one_physical_CPU"]],
    BLAS_threads=settings["BLAS_threads"],max_iterations=settings["max_iterations"],
    parent_preparation_receipt_sha256=filehash(diagnostic_path),
    identity_note=diagnostic["metadata_runtime_note"]))
 new={"start_time_ps":now,"end_time_ps":None,"origin":"v4",
  **{k:identity[k] for k in ("dt_ps","step_tolerance","threads")}}
 segments.append(new)
 if horizon_change is not None:
  segments[-2].update(origin="completed_control",duration_ps=horizon_change["before_ps"])
  new.update(origin="control_horizon_extension",duration_ps=horizon_change["after_ps"])
 if guard_change is not None:
  for prior in segments[:-1]:
   prior.setdefault("energy_guard_enabled",prior.get("local_energy_tolerance") is not None)
   prior.setdefault("local_energy_tolerance",None)
  new.update(origin="v5_local_energy_guard",energy_guard_enabled=True,
   local_energy_tolerance=identity["local_energy_tolerance"],
   local_energy_reference_energy_J=guard_change["reference_energy_J"],
   guard_active_since_time_ps=guard_change["ledger_after"]["guard_active_since_time_ps"])
 return segments

def transfer(anchor,output,repo,dt_ps,tolerance,threads,cpu_start,*,duration_ps_override='unchanged',local_energy_tolerance=None):
 anchor=scratch(anchor);output=scratch(output);repo=Path(repo).resolve()
 # Never reinterpret a partial directory as an admitted anchor.
 cr=json.loads((anchor/"capture_receipt.json").read_text())
 payload=(anchor/"checkpoint_original.npz").read_bytes()
 if digest(payload)!=cr["checkpoint_bytes_sha256"]:raise ValueError("Anchored checkpoint changed")
 for rel,want in cr["physical_file_sha256"].items():
  if filehash(anchor/rel)!=want:raise ValueError("Anchored physical file changed "+rel)
 arrays=decode(payload);parent=json.loads((anchor/"manifest_original.json").read_text())
 oldmeta=audit_arrays(arrays,parent);oldpins=oldmeta["source_pins_sha256"]
 new_duration,horizon_change=extend_control_horizon(parent,oldmeta,duration_ps_override)
 oldaudit=audit_runtime(anchor/"runtime_sources",oldpins)
 newpins={rel:filehash(repo/rel) for rel in oldpins}
 changes={rel:{"old":oldpins[rel],"new":newpins[rel],"old_path":str(anchor/"runtime_sources"/rel),
  "new_path":str(repo/rel)} for rel in oldpins if oldpins[rel]!=newpins[rel]}
 if set(changes)-FRONTEND:raise ValueError("Unapproved RHS/source changes: "+repr(set(changes)-FRONTEND))
 if not all(np.isfinite(v) and v>0 for v in (dt_ps,tolerance)) or threads<1 or cpu_start<0:
  raise ValueError("Invalid declared numerical/resource parameters")
 guard_capable="--local-energy-tolerance" in (repo/"pipelines/05_effective_photon_run.py").read_text()
 if guard_capable and local_energy_tolerance is None:
  raise ValueError("New guarded runner requires explicit --local-energy-tolerance; absence does not certify old history")
 if not guard_capable and local_energy_tolerance is not None:
  raise ValueError("New runner does not support local-energy guarding")
 guard_ledger,guard_change=energy_guard_transition(oldmeta,arrays,parent,local_energy_tolerance)
 identity=copy.deepcopy(oldmeta["parameter_identity"])
 identity.update(dt_ps=float(dt_ps),step_tolerance=float(tolerance),threads=int(threads))
 if local_energy_tolerance is not None:identity["local_energy_tolerance"]=float(local_energy_tolerance)
 if horizon_change is not None:identity["duration_ps"]=new_duration
 output=newdir(output);prov=output/"provenance";prov.mkdir()
 shutil.copy2(anchor/"checkpoint_original.npz",prov/"parent_checkpoint.npz")
 shutil.copy2(anchor/"manifest_original.json",prov/"parent_manifest.json")
 shutil.copy2(anchor/"capture_receipt.json",prov/"capture_receipt.json")
 if (anchor/"parent_preparation_receipt.json").exists():
  shutil.copy2(anchor/"parent_preparation_receipt.json",prov/"parent_preparation_receipt.json")
 if (anchor/"diagnostic_provenance").is_dir():
  shutil.copytree(anchor/"diagnostic_provenance",prov/"diagnostic_provenance")
 shutil.copytree(anchor/"runtime_sources",output/"runtime_sources")
 rt=output/"runtime_sources"
 for rel in FRONTEND:
  if rel in oldpins:shutil.copy2(repo/rel,rt/rel)
 shutil.copy2(rt/"archive_manifest.json",prov/"parent_runtime_archive_manifest.json")
 if (rt/"verification_receipt.json").exists():
  shutil.move(rt/"verification_receipt.json",prov/"parent_runtime_verification_receipt.json")
 for rel in ("dc","meshes"):
  if (anchor/rel).exists():shutil.copytree(anchor/rel,output/rel)
 for rel in ("material.npz","material.json","source_initial_state.npz"):shutil.copy2(anchor/rel,output/rel)
 newmeta=copy.deepcopy(oldmeta)
 newmeta.update(source_pins_sha256=newpins,parameter_identity=identity,
  next_step_ps=min(float(oldmeta["next_step_ps"]),dt_ps))
 if guard_ledger is not None:newmeta.update(guard_ledger)
 missing=[]
 for k in ("stable_accepted_steps","consecutive_rejections"):
  if k not in newmeta:missing.append(k);newmeta[k]=0
 # Complete histories and physical arrays are preserved byte-for-byte as arrays;
 # only explicit metadata above is repacked.
 arrays["metadata_json"]=np.array(json.dumps(newmeta,allow_nan=False))
 with (output/"checkpoint_latest.npz").open("wb") as s:np.savez_compressed(s,**arrays)
 np.savez_compressed(output/"trajectory.npz",**{k[12:]:v for k,v in arrays.items() if k.startswith("trajectory__")})
 newmanifest=copy.deepcopy(parent)
 newmanifest.update(source_pins_sha256=newpins,parameter_identity=identity,status="EXPLICIT_TRANSFER_PREPARED",
  run_complete=False,explicit_transfer_receipt="provenance/transfer_receipt.json")
 # Completed execution/output hashes refer to the parent, never to this campaign.
 for k in ("execution","output_sha256","error"):newmanifest.pop(k,None)
 newmanifest["numerical"].update(requested_step_ps=dt_ps,step_tolerance=tolerance)
 if guard_change is not None:
  newmanifest["local_energy_control"]=copy.deepcopy(guard_ledger)
  newmanifest["numerical"].update(local_energy_tolerance=float(local_energy_tolerance),
   local_energy_guard=dict(enabled=True,reference_energy_J=guard_change["reference_energy_J"],
    reference_kind=guard_change["reference_kind"],acceptance=guard_change["denominator"],
    pre_guard_history="Measured accepted increments only; never retroactively certified"))

 if horizon_change is not None:
  newmanifest["explicit_horizon_change"]=horizon_change
  newmanifest.setdefault("observation",{}).update(status="EXPLICIT_TRANSFER_PREPARED",development_duration_s=None if new_duration is None else new_duration*1e-12)
 writejson(output/"manifest.json",newmanifest)
 # Verify physical reconstruction on the NEW FROZEN capsule, no evolution,
 # no initial photon preparation. Limit verification to one physical CPU.
 runner=load_runner(rt)
 argv=resume_argv(identity,output,cpu_start)
 a=runner.parser().parse_args(argv);a.threads=1
 system,dc,g,material,params=runner.build_system(a,output,resume=True)
 if params!=parent["physicalparams"]:raise ValueError("Physical parameter identity changed")
 if material.manifest!=parent["material"]:raise ValueError("Effective material identity changed")
 if newmeta["controller"] is not None:runner.TriggerWindow.from_checkpoint(newmeta["controller"])
 if g.receipt["mesh_sha256"]!=parent["geometry"]["mesh_sha256"]:raise ValueError("Mesh changed")
 state=runner.EffectiveState(*(arrays[k].copy() for k in STATE))
 rhs=system.evaluate(state)
 if guard_change is not None:
  if not math.isclose(float(system.kinetic_residual_energy_scale_J),guard_change["reference_energy_J"],rel_tol=1e-13,abs_tol=0.):
   raise ValueError("Reconstructed system guard reference differs")
 close_scalar(float(rhs.energy_J),float(arrays["trajectory__energy_total_J"][-1]),"reconstructed total energy")
 snap={k[11:]:v for k,v in arrays.items() if k.startswith("snapshots__")}
 snap.update(coordinates_m=g.coordinates_m,triangles=g.triangles,edges=g.graph.edges)
 np.savez_compressed(output/"snapshots.npz",**snap)
 check=decode((output/"checkpoint_latest.npz").read_bytes())
 newmeta=audit_arrays(check,newmanifest)
 for k in arrays:
  if k!="metadata_json" and (arrays[k].dtype!=check[k].dtype or not np.array_equal(arrays[k],check[k])):
   raise ValueError("Transfer changed array "+k)
 for k in CONSERVED:
  if oldmeta[k]!=newmeta[k]:raise ValueError("Transfer changed conserved metadata "+k)
 for rel in oldpins:
  if filehash(rt/rel)!=newpins[rel]:raise ValueError("New runtime pin mismatch "+rel)
 oldarchive=json.loads((anchor/"runtime_sources/archive_manifest.json").read_text())
 newarchive=copy.deepcopy(oldarchive)
 newarchive.update(archive_root=str(rt),run_dir=str(output),runtime_pins_sha256=newpins,
  scope="Complete frozen source capsule for explicit frontend/solver continuation; parent is independently preserved")
 for rel in FRONTEND:
  if rel in newarchive["package_sources_and_small_geometry"]:
   newarchive["package_sources_and_small_geometry"][rel]={"sha256":filehash(rt/rel),"bytes":(rt/rel).stat().st_size,"origin":str(repo/rel)}
 newarchive["data_references"]={"parent_capture":str(anchor),"physical_files_in_new_run":True}
 writejson(rt/"archive_manifest.json",newarchive)
 audit_runtime(rt,newpins)
 writejson(rt/"verification_receipt.json",{"schema":"pysnspd.runtime_transfer_verification.v1","source_pins_match":True,"physical_parameters_equal":True,"DC_loaded_not_reprepared":True,"mesh_sha256":g.receipt["mesh_sha256"],"photon_reapplied":False,"evolution_executed":False})
 segments=numerical_segments_for_transfer(anchor,oldmeta,identity,horizon_change,guard_change)
 receipt=dict(schema=SCHEMA,action="EXPLICIT_TRANSFER_PREPARED_NO_EVOLUTION_NO_SIGNAL",
  parent={"output_dir":cr["parent_output_dir"],"anchor":str(anchor),
   "checkpoint_bytes_sha256":digest(payload),"manifest_sha256":cr["manifest_sha256"],
   "runtime_archive_manifest_sha256":oldaudit["manifest_sha256"]},
  output_dir=str(output),state_arrays_sha256_before={k:arrayhash(arrays[k]) for k in STATE},
  state_arrays_sha256_after={k:arrayhash(check[k]) for k in STATE},
  conserved_metadata_before={k:oldmeta[k] for k in CONSERVED},
  conserved_metadata_after={k:newmeta[k] for k in CONSERVED},
  pins_before=oldpins,pins_after=newpins,frontend_changes=changes,RHS_sources_identical=True,
  numerical_segments=segments,horizon_change=horizon_change,energy_guard_change=guard_change,
  parent_helper_sha256=PARENT_HELPER_SHA256,missing_old_adaptive_counters=missing,
  adaptive_counter_policy="Preserve recorded counters exactly; v3 unrecorded counters explicitly restart at zero",
  next_step_ps_before=oldmeta["next_step_ps"],next_step_ps_after=newmeta["next_step_ps"],
  preserved_history_sha256={k[12:]:arrayhash(v) for k,v in arrays.items() if k.startswith("trajectory__")},
  preserved_snapshots_sha256={k[11:]:arrayhash(v) for k,v in arrays.items() if k.startswith("snapshots__")},
  physical_identity_verified={"parameters_equal":True,"mesh_sha256":g.receipt["mesh_sha256"],
   "DC_sha256":filehash(output/"dc/effective_dc.npz"),"material_sha256":filehash(output/"material.npz"),
   "reconstructed_total_energy_J":float(rhs.energy_J),"photon_source_reapplied":False},
  new_checkpoint_bytes_sha256=filehash(output/"checkpoint_latest.npz"),
  new_runtime_archive_manifest_sha256=filehash(rt/"archive_manifest.json"),
  resume_argv=argv,admission=False)
 writejson(prov/"transfer_receipt.json",receipt)
 command="env PYTHONDONTWRITEBYTECODE=1 PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/jdiaz/.conda/envs/snspd/bin/python "+shlex.quote(str(rt/"pipelines/05_effective_photon_run.py"))+" "+shlex.join(argv)
 (prov/"resume_command.txt").write_text(command+"\n")
 return receipt

def main():
 p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest="command",required=True)
 c=sub.add_parser("capture");c.add_argument("--run",required=True);c.add_argument("--anchor",required=True)
 t=sub.add_parser("transfer");t.add_argument("--anchor",required=True);t.add_argument("--output",required=True)
 t.add_argument("--repo",default=str(REPO));t.add_argument("--dt-ps",required=True,type=float)
 t.add_argument("--step-tolerance",required=True,type=float);t.add_argument("--threads",required=True,type=int)
 t.add_argument("--cpu-start",required=True,type=int)
 t.add_argument("--local-energy-tolerance",type=float,default=None,
  help="Explicit new guard tolerance; old absence means unguarded history, never retroactive certification")
 t.add_argument("--duration-ps",default="unchanged",type=duration_argument,help="Explicit completed ACTIVE control extension: larger physical duration or none; existing trigger/ceiling unchanged")
 a=p.parse_args()
 if a.command=="capture":r=capture(a.run,a.anchor)
 else:r=transfer(a.anchor,a.output,a.repo,a.dt_ps,a.step_tolerance,a.threads,a.cpu_start,duration_ps_override=a.duration_ps,local_energy_tolerance=a.local_energy_tolerance)
 print(json.dumps({"status":r["action"],"receipt":str(Path(a.anchor if a.command=="capture" else a.output)/("capture_receipt.json" if a.command=="capture" else "provenance/transfer_receipt.json")),"time_ps":r.get("conserved_metadata",r.get("conserved_metadata_after",{})).get("time_ps")},indent=2))
if __name__=="__main__":main()
