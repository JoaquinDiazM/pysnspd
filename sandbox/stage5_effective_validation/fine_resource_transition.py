"""Review-only fine control resource transition; no signals or time evolution.
The sole state change is parameter_identity.threads:2->4. The old guard window,
all adaptive counters, accepted data and energy origin remain exact.
"""
from __future__ import annotations
import argparse,ast,copy,hashlib,importlib.util,json,os,shutil,sys
from pathlib import Path
import numpy as np
SCRATCH=Path("/home/jdiaz/scratch")
ROOT=Path("/home/jdiaz/pysnspd")
PARENT=SCRATCH/"stage5_effective_L320_h35_775_adaptive_v6_guard2b_20261011"
NOMINAL=SCRATCH/"stage5_effective_L320_h5_775_adaptive_v6_guard2b_20261011"
CHILD=SCRATCH/"stage5_effective_L320_h35_775_adaptive_v6_guard2b_resources4_20261011"
RUNNER_SHA="70b7c5b5d39c74f19a35ad816f413ecc47451c46651817708549b8e9fc8e97aa"
STATE_KEYS=("delta_bar","p","n","circuit_SI")
GUARD_HELPER=SCRATCH/"stage5_effective_controls_20261011/tighten_local_energy_guard_v6_guard2b.py"
GUARD_HELPER_SHA="567d90a84c9c89083b37d759e8a07ee04377c69d57810cc19f402850ddf3de41"
BOUNDED_SHA="b718d5a41874771e5afcd706ce4a4bd730fea44cf3f1074315b15878d9e8e5a6"
BASELINE_SHA="3c08c2534efd2ba8d5803842a9e9cf347b6befc46e68ebddc15adb716156c27a"
CERTIFICATE=SCRATCH/"stage5_node_fusion_candidate_20261011/equivalence_certificate.json"
CERTIFICATE_SHA="a7e0170f5c4b295a89ede835d7c234f4ebda8919f4271d2098b60fee638f937f"
APPLICATION=CERTIFICATE.with_name("root_application_receipt.json")
APPLICATION_SHA="14390a5da13607439036268e17a6be73a0912e838374d8b1400628ecc5dcb4e2"
ACTION="EXPLICIT_FINE_RESOURCE_ONLY_TRANSITION_PREPARED_NO_EVOLUTION_NO_SIGNAL"
ORIGIN="explicit_fine_resource_only_transition"
EVENTS=("candidate_crossing_s","crossing_time_s","confirmation_time_s","stop_target_s","actual_stop_s","failure_reason")
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def entry(p):return dict(path=str(Path(p).resolve()),sha256=digest(p))
def owned(p):
 p=Path(p).resolve()
 if os.getuid()!=1011 or not p.is_relative_to(SCRATCH) or p==SCRATCH or p.stat().st_uid!=1011:
  raise ValueError("jdiaz-owned scratch evidence required")
 return p
def bound(e):
 if not isinstance(e,dict) or set(e)!={"path","sha256"}:raise ValueError("Explicit path/SHA required")
 p=owned(e["path"])
 if digest(p)!=e["sha256"]:raise ValueError("Bound evidence changed")
 return p
def load(path,name,sha):
 if digest(path)!=sha:raise ValueError("Reviewed dependency changed")
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);sys.modules[name]=m
 spec.loader.exec_module(m);return m
def dependencies():
 g=load(GUARD_HELPER,"_fine_resource_guard_helper",GUARD_HELPER_SHA);g.helper_integrity()
 from sandbox.stage5_effective_validation import bounded_guard_protocol as guard,audit_runs as audit
 if digest(guard.__file__)!=BOUNDED_SHA or digest(audit.__file__)!=BASELINE_SHA:raise ValueError("Baseline/guard changed")
 return g,guard,audit
def certificate(pins):
 if digest(CERTIFICATE)!=CERTIFICATE_SHA or digest(APPLICATION)!=APPLICATION_SHA:raise ValueError("Existing certificate changed")
 c=read(CERTIFICATE);a=read(APPLICATION)
 if a.get("certificate_sha256")!=CERTIFICATE_SHA or a.get("applied_to_main") is not True or a["sources_after"]!=c["sources_after"]:
  raise ValueError("Root-applied storage optimization proof required")
 if any(pins.get(k)!=v for k,v in c["sources_after"].items()):raise ValueError("Certified runtime sources differ")
 return dict(certificate=entry(CERTIFICATE),root_application=entry(APPLICATION),
  scope="Existing tested1/2/4-thread algebraic evidence only; future trajectories and cross-thread scalar reductions are not asserted bitwise identical")
def checkpoint_namespace():
 """The pinned runner writes27numeric arrays plus metadata_json:28total keys."""
 runner=ROOT/"pipelines/05_effective_photon_run.py"
 if digest(runner)!=RUNNER_SHA:raise ValueError("Checkpoint producer source changed")
 fields={}
 for node in ast.parse(runner.read_text()).body:
  if isinstance(node,ast.Assign):
   for target in node.targets:
    if isinstance(target,ast.Name) and target.id in {"HISTORY_KEYS","SNAPSHOT_KEYS"}:fields[target.id]=ast.literal_eval(node.value)
 if set(fields)!={"HISTORY_KEYS","SNAPSHOT_KEYS"}:raise ValueError("Pinned producer schema missing")
 names=set(STATE_KEYS)|{"metadata_json"}|{"trajectory__"+k for k in fields["HISTORY_KEYS"]}|{"snapshots__"+k for k in fields["SNAPSHOT_KEYS"]}
 if len(fields["HISTORY_KEYS"])!=16 or len(fields["SNAPSHOT_KEYS"])!=7 or len(names)!=28:raise ValueError("Pinned runner schema differs from reviewed28-key producer")
 return frozenset(names)
def parent_contract(manifest,m,arrays):
 if manifest.get("status")!="INTERRUPTED" or manifest.get("run_complete") is not False:
  raise ValueError("Stable INTERRUPTED unfinished fine parent required")
 i=m["parameter_identity"]
 expected=dict(threads=2,local_energy_tolerance=2.5e-6,step_tolerance=2e-7,dt_ps=.02,
  length_nm=320.,spacing_nm=3.5,wavelength_nm=775.,xi_cutoff=32.,xi_order=6,phonon_nodes=96,transport_order=48,duration_ps=None)
 if any(isinstance(i.get(k),bool) or i.get(k)!=v for k,v in expected.items()):raise ValueError("Fixed fine guard2b case/solver/resources required")
 if manifest["parameter_identity"]!=i or set(arrays)!=checkpoint_namespace() or len(m["source_pins_sha256"])!=13:
  raise ValueError("Exact28-key schema(27numeric arrays+metadata)/13-source namespace required")
 s=m["controller"]["state"]
 if s.get("status")!="ACTIVE" or any(s.get(k) is not None for k in EVENTS) or s.get("rejected_candidate_count")!=0:
  raise ValueError("Only an ACTIVE pre-trigger parent without prior/pending events is eligible")
 if not 22.185<=m["guard_active_since_time_ps"]<=22.435 or m["time_ps"]<m["guard_active_since_time_ps"]:
  raise ValueError("Actual approved fine guard2b anchor required")
 r=manifest["resources"]
 if r.get("actual_threads")!=2 or r.get("requested_threads")!=2 or r.get("actual_affinity_cpus")!=[12,13]:
  raise ValueError("Observed Num2 fine affinity12,13 required")
 return True
def same_bytes(a,b):
 return a.dtype==b.dtype and a.shape==b.shape and a.tobytes(order="C")==b.tobytes(order="C")
def exact_json(a,b):
 return json.dumps(a,sort_keys=True,allow_nan=False)==json.dumps(b,sort_keys=True,allow_nan=False)
def true_check(value):
 """Only actual Python/NumPy boolean truth may certify a saved audit check."""
 return isinstance(value,(bool,np.bool_)) and bool(value)
AUDIT_OUTPUTS=("trajectory.npz","snapshots.npz","checkpoint_latest.npz")
def audit_manifest_projection(original,actual):
 """Fill only native INTERRUPTED missing output hashes; never repair existing ones."""
 if not isinstance(original,dict) or not isinstance(actual,dict):
  raise ValueError("Original manifest and independently recomputed hashes required")
 required=set(AUDIT_OUTPUTS)
 if not required<=set(actual) or any(not isinstance(v,str) or len(v)!=64 or any(c not in "0123456789abcdef" for c in v) for v in actual.values()):
  raise ValueError("Exact independently recomputed output SHA256 values required")
 expected=copy.deepcopy(original)
 if "output_sha256" in original:
  advertised=original["output_sha256"]
  if not isinstance(advertised,dict) or not required<=set(advertised) or any(k not in actual or v!=actual[k] for k,v in advertised.items()):
   raise ValueError("Existing saved output checksums are invalid; they must never be repaired")
  return expected,"Existing output_sha256 authenticated unchanged"
 if original.get("status")!="INTERRUPTED" or original.get("run_complete") is not False:
  raise ValueError("Only native incomplete INTERRUPTED manifests may lack output_sha256")
 expected["output_sha256"]={k:actual[k] for k in AUDIT_OUTPUTS}
 return expected,"Native INTERRUPTED runner omits output_sha256; derived audit view adds only recomputed saved-output hashes"
def audit_manifest_view(anchor,original,*,write_json=None):
 """Authenticate original bytes separately from their exact derived audit view."""
 anchor=owned(anchor);original_path=owned(anchor/"manifest_original.json")
 original_bytes=original_path.read_bytes()
 if not exact_json(json.loads(original_bytes),original):raise ValueError("Original captured manifest differs")
 advertised=original.get("output_sha256",{})
 if not isinstance(advertised,dict):raise ValueError("Existing output_sha256 must be a dictionary")
 actual={}
 for name in sorted(set(AUDIT_OUTPUTS)|set(advertised)):
  if not isinstance(name,str) or Path(name).name!=name:raise ValueError("Saved output must be an owned captured file")
  actual[name]=digest(owned(anchor/name))
 expected,reason=audit_manifest_projection(original,actual)
 derived_path=anchor/"manifest.json"
 if write_json is not None:write_json(derived_path,expected)
 if not exact_json(read(owned(derived_path)),expected):raise ValueError("Derived audit manifest must equal original plus only authenticated output_sha256")
 if original_path.read_bytes()!=original_bytes:raise ValueError("Original manifest bytes changed")
 augmentation=dict(original_manifest_sha256=hashlib.sha256(original_bytes).hexdigest(),
  derived_manifest_sha256=digest(derived_path),recomputed_output_sha256=actual,
  output_sha256_added="output_sha256" not in original,reason=reason,
  original_manifest_bytes_unchanged=True,physical_files_modified=False)
 return expected,augmentation
def transition_arrays(before,initial):
 if set(before)!=checkpoint_namespace() or set(initial)!=checkpoint_namespace():raise ValueError("Exact28-key schema(27numeric arrays+metadata) required")
 bm=json.loads(str(before["metadata_json"].item()));im=json.loads(str(initial["metadata_json"].item()))
 expected=copy.deepcopy(bm);expected["parameter_identity"]["threads"]=4
 if not exact_json(im,expected):raise ValueError("Metadata may change ONLY threads2->4")
 if bm["parameter_identity"].get("threads")!=2:raise ValueError("Parent threads must be2")
 for k in before:
  if k!="metadata_json" and not same_bytes(before[k],initial[k]):
   raise ValueError("Physical/history array changed: "+k)
 return bm,im
def live_prefix(initial,live):
 if set(initial)!=checkpoint_namespace() or set(live)!=checkpoint_namespace():raise ValueError("Advanced28-key namespace changed")
 old=json.loads(str(initial["metadata_json"].item()));now=json.loads(str(live["metadata_json"].item()))
 for k in ("parameter_identity","source_pins_sha256","source","initial_total_energy_J",
           "guard_active_since_time_ps","local_energy_reference_energy_J","energy_guard_ledger_origin"):
  if not exact_json(old[k],now[k]):raise ValueError("Advanced identity/global origin/guard window changed: "+k)
 for k,v in initial.items():
  if k.startswith(("trajectory__","snapshots__")):
   w=live[k]
   if v.dtype!=w.dtype or v.ndim!=w.ndim or v.shape[1:]!=w.shape[1:] or len(w)<len(v) or not same_bytes(v,w[:len(v)]):
    raise ValueError("Advanced accepted prefix changed: "+k)
 if now["time_ps"]<old["time_ps"] or now["accepted_steps"]<old["accepted_steps"]:raise ValueError("Advanced clock/counter reversed")
 for key in ("total_rejections","local_energy_rejections","limited_event_energy_J","totalabs_event_energy_J",
  "absolute_external_energy_integral_J","sum_abs_local_energy_defect_J","max_accepted_local_energy_fraction","max_guarded_accepted_local_energy_fraction"):
  if now[key]<old[key]:raise ValueError("Cumulative resource-prefix counter/ledger reset: "+key)
 if now["time_ps"]==old["time_ps"]:
  for key,value in initial.items():
   if key!="metadata_json" and not same_bytes(value,live[key]):
    raise ValueError("Same-clock state changed: "+key)
 return now
def resource_segments(segments,m):
 rows=copy.deepcopy(segments);now=m["time_ps"]
 if rows[-1]["threads"]!=2 or rows[-1]["local_energy_tolerance"]!=2.5e-6 or rows[-1]["guard_active_since_time_ps"]!=m["guard_active_since_time_ps"]:
  raise ValueError("Parent numerical/resource timeline differs")
 nxt=copy.deepcopy(rows[-1]);rows[-1]["end_time_ps"]=now
 nxt.update(start_time_ps=now,end_time_ps=None,threads=4,origin=ORIGIN,
  resource_affinity_cpus=[0,1,2,3],resource_cpu_start=0,
  guard_window_preserved=True)
 rows.append(nxt);return rows
def guard_projection(segments,resource_time):
 """Explicit guard-policy projection only; full resource chronology is retained."""
 rows=copy.deepcopy(segments)
 if len(rows)<2 or rows[-1].get("origin")!=ORIGIN or rows[-1]["start_time_ps"]!=resource_time:
  raise ValueError("Exactly one final resource-only tail required")
 last=rows.pop();prev=rows[-1]
 expected=copy.deepcopy(prev);expected.update(start_time_ps=last["start_time_ps"],end_time_ps=None,threads=4,origin=ORIGIN,
  resource_affinity_cpus=[0,1,2,3],resource_cpu_start=0,guard_window_preserved=True)
 if expected!=last or prev["end_time_ps"]!=resource_time:raise ValueError("Resource tail hides numerical/guard changes")
 prev["end_time_ps"]=None
 for row in rows:
  for key in ("threads","resource_affinity_cpus","resource_cpu_start"):row.pop(key,None)
 return rows
def process_owners():
 found=[]
 for p in Path("/proc").glob("[0-9]*"):
  try:
   if p.stat().st_uid!=1011:continue
   argv=[x.decode() for x in (p/"cmdline").read_bytes().split(b"\0") if x]
   if not any(x.endswith("/05_effective_photon_run.py") for x in argv) or "--output-dir" not in argv:continue
   found.append(dict(pid=int(p.name),output_dir=str(Path(argv[argv.index("--output-dir")+1]).resolve()),
    affinity=sorted(os.sched_getaffinity(int(p.name)))))
  except (OSError,ValueError,UnicodeError):continue
 return found
def resource_preflight(nominal_manifest,owners,inventory,budget):
 if nominal_manifest.get("run_complete") is not True or nominal_manifest.get("status") not in ("TRIGGER_MARGIN_COMPLETE","TIME_CENSORED_NO_TRIGGER"):
  raise ValueError("Nominal must actually finish its terminal observation first")
 if any(x["output_dir"] in {str(PARENT),str(NOMINAL),str(CHILD)} for x in owners):raise ValueError("Parent/nominal/child still has a running owner")
 groups=inventory["physical_core_groups"];allowed=set(budget["selected_affinity_cpus"])
 admitted=[g for g in groups if set(g)<=allowed]
 if [g[0] for g in admitted[:4]]!=[0,1,2,3] or len(admitted)>14 or any(set(g)&{14,15,30,31} for g in admitted):
  raise ValueError("Full host admission must reserve physical14/15 and select target0..3")
 occupied={j for x in owners for j,g in enumerate(groups) if set(g)&set(x["affinity"])}
 target={j for j,g in enumerate(groups) if set(g)&{0,1,2,3}}
 if occupied&{14,15} or occupied&target or len(occupied|target)>14 or len(occupied|target)>.9*len(groups):
  raise ValueError("Target cores occupied or shared90percent global budget exceeded")
 return dict(target_affinity_cpus=[0,1,2,3],cpu_start=0,requested_threads=4,
  active_pipeline_owners=owners,physical_cores_after=len(occupied|target),
  reserved_physical_cores=[14,15],global_fraction_ceiling=.9,
  scope="Point-in-time pipeline-owner inventory; recheck immediately before launch; no reservation or launch performed")
def preflight():
 g,guard,audit=dependencies()
 from sandbox.stage4_core.parallel_runtime import linux_resources,resource_budget
 inv=linux_resources();budget=resource_budget(inv,.9,.9,max_workers=1)
 proof=resource_preflight(read(owned(NOMINAL/"manifest.json")),process_owners(),inv,budget)
 proof.update(nominal_terminal_manifest=entry(owned(NOMINAL/"manifest.json")),runtime_inventory=inv,scheduling_budget=budget)
 return proof
def prepared_manifest(parent,identity,plan):
 out=copy.deepcopy(parent);out.update(parameter_identity=identity,status="EXPLICIT_RESOURCE_TRANSITION_PREPARED",
  explicit_transfer_receipt="provenance/transfer_receipt.json",explicit_resource_transition=plan)
 # resources retains the last actually observed parent execution until runner STARTED.
 return out
def prepare(anchor,output,protocol_path):
 g,guard,audit=dependencies();h=g.h
 if Path(output).resolve()!=CHILD:raise ValueError("Only the declared fine Num4 child is in scope")
 anchor=owned(anchor);cr=read(anchor/"capture_receipt.json")
 if cr["parent_output_dir"]!=str(PARENT):raise ValueError("Wrong fixed fine parent")
 parent_payload=(anchor/"checkpoint_original.npz").read_bytes();manifest_payload=(anchor/"manifest_original.json").read_bytes()
 if hashlib.sha256(parent_payload).hexdigest()!=cr["checkpoint_bytes_sha256"] or hashlib.sha256(manifest_payload).hexdigest()!=cr["manifest_sha256"]:
  raise ValueError("Immutable capture changed")
 for rel,want in cr["physical_file_sha256"].items():
  if digest(owned(anchor/rel))!=want:raise ValueError("Captured evidence changed: "+rel)
 if (owned(PARENT/"checkpoint_latest.npz").read_bytes()!=parent_payload or owned(PARENT/"manifest.json").read_bytes()!=manifest_payload):
  raise ValueError("Stopped parent changed since capture")
 before=h.decode(parent_payload);parent=json.loads(manifest_payload);m=h.audit_arrays(before,parent)
 parent_contract(parent,m,before);g.eligible(parent,m,before);cpuprooff=preflight()
 os.sched_setaffinity(0,{4})  # Metadata/audit work only, after checking the full host budget.
 runtime=h.audit_runtime(anchor/"runtime_sources",m["source_pins_sha256"]);cert=certificate(m["source_pins_sha256"])
 if m["source_pins_sha256"]!=read(NOMINAL/"manifest.json")["source_pins_sha256"]:raise ValueError("Nominal/fine must retain identical13 frozen runtime sources")
 for rel in m["source_pins_sha256"]:owned(anchor/"runtime_sources"/rel)
 # Keep original captured bytes; derive only the hashes omitted by native INTERRUPTED.
 shutil.copy2(anchor/"manifest_original.json",anchor/"manifest.json")
 shutil.copy2(anchor/"checkpoint_original.npz",anchor/"checkpoint_latest.npz")
 for name in ("trajectory.npz","snapshots.npz"):shutil.copy2(owned(PARENT/name),anchor/name)
 _,audit_view=audit_manifest_view(anchor,parent,write_json=h.writejson)
 parentrun=audit.audit_run("resource_parent",anchor)
 context=guard.load_protocol(protocol_path);cfg=context[0]
 pp=guard.verify_run(parentrun,cfg["controls"]["mesh"],"mesh",context)
 if any(not true_check(v) for k,v in parentrun["summary"]["checks"].items() if k not in {"complete_declared_window","successful_saved_window_status"}):
  raise ValueError("Parent independent physical/guard audit failed")
 if (PARENT/"checkpoint_latest.npz").read_bytes()!=parent_payload or (PARENT/"manifest.json").read_bytes()!=manifest_payload:
  raise ValueError("Stopped parent changed during preparation")
 identity=copy.deepcopy(m["parameter_identity"]);identity["threads"]=4
 numerical=resource_segments(pp["segments"],m);h.lineage.segments_valid(numerical)
 plan=dict(requested_threads=4,target_affinity_cpus=[0,1,2,3],cpu_start=0,
  preparation_only=True,execution_started=False,previous_observed_resources=copy.deepcopy(parent["resources"]))
 output=h.newdir(output);prov=output/"provenance";prov.mkdir()
 for rel in ("runtime_sources","dc","meshes"):shutil.copytree(anchor/rel,output/rel)
 for rel in ("material.npz","material.json","source_initial_state.npz","trajectory.npz","snapshots.npz"):shutil.copy2(anchor/rel,output/rel)
 new=copy.deepcopy(m);new["parameter_identity"]=identity;initial=copy.deepcopy(before);initial["metadata_json"]=np.array(json.dumps(new,allow_nan=False))
 np.savez_compressed(output/"checkpoint_latest.npz",**initial)
 manifest=prepared_manifest(parent,identity,plan);h.writejson(output/"manifest.json",manifest)
 transition_arrays(before,h.decode((output/"checkpoint_latest.npz").read_bytes()))
 h.audit_arrays(initial,manifest)
 for src,dst in [("checkpoint_original.npz","parent_checkpoint.npz"),("manifest_original.json","parent_manifest.json"),("capture_receipt.json","capture_receipt.json")]:
  shutil.copy2(anchor/src,prov/dst)
 shutil.copy2(anchor/"provenance/parent_transfer_receipt.json",prov/"parent_transfer_receipt.json")
 shutil.copy2(output/"checkpoint_latest.npz",prov/"initial_resource_checkpoint.npz")
 shutil.copy2(output/"manifest.json",prov/"initial_resource_manifest.json")
 argv=h.resume_argv(identity,output,0)
 r=dict(schema="pysnspd.explicit_checkpoint_transfer.v3",action=ACTION,helper_sha256=digest(__file__),
  parent=dict(output_dir=str(PARENT),anchor=str(anchor),checkpoint_bytes_sha256=cr["checkpoint_bytes_sha256"],
   manifest_sha256=cr["manifest_sha256"],runtime_archive_manifest_sha256=runtime["manifest_sha256"]),
  output_dir=str(output),pins_before=m["source_pins_sha256"],pins_after=m["source_pins_sha256"],
  state_arrays_sha256_before={k:h.arrayhash(before[k]) for k in h.STATE},state_arrays_sha256_after={k:h.arrayhash(initial[k]) for k in h.STATE},
  conserved_metadata_before={k:m[k] for k in h.CONSERVED},conserved_metadata_after={k:new[k] for k in h.CONSERVED},
  preserved_array_hashes={k:h.arrayhash(v) for k,v in before.items() if k!="metadata_json"},
  checkpoint_total_keys=28,checkpoint_numeric_arrays=27,checkpoint_namespace=sorted(checkpoint_namespace()),all_metadata_except_threads_preserved=True,
  frontend_changes={},implementation_changes={},horizon_change=None,energy_guard_change=None,RHS_sources_identical=True,
  guard_ledger_before={k:m[k] for k in h.GUARD_KEYS},guard_ledger_after={k:new[k] for k in h.GUARD_KEYS},
  next_step_ps_before=m["next_step_ps"],next_step_ps_after=new["next_step_ps"],
  global_energy_error_before=g.global_error_at_anchor(m,before),global_energy_error_after=g.global_error_at_anchor(new,initial),
  resources_before=parent["resources"],resources_after=plan,numerical_segments=numerical,
  guard_policy_projection=guard_projection(numerical,m["time_ps"]),resource_preflight=cpuprooff,
  parent_guard_entry=copy.deepcopy(cfg["controls"]["mesh"]),parent_guard_verification=pp,
  parent_guard_protocol=entry(protocol_path),existing_thread_evidence=cert,
  initial_checkpoint=entry(prov/"initial_resource_checkpoint.npz"),initial_manifest=entry(prov/"initial_resource_manifest.json"),
  parent_audit_view_manifest=audit_view,
  parent_original_audit_files={name:entry(anchor/name) for name in ("manifest.json","checkpoint_latest.npz","trajectory.npz","snapshots.npz")},
  parent_runtime_archive_manifest=entry(anchor/"runtime_sources/archive_manifest.json"),
  new_checkpoint_bytes_sha256=digest(output/"checkpoint_latest.npz"),resume_argv=argv,
  state_or_heat_correction=False,global_energy_origin_reset=False,guard_window_reset=False,
  future_cross_thread_bitwise_identity_claimed=False,admission=False)
 h.writejson(prov/"transfer_receipt.json",r);return r
def main():
 p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest="command",required=True)
 sub.add_parser("preflight")
 t=sub.add_parser("prepare");t.add_argument("--anchor",required=True);t.add_argument("--output",required=True);t.add_argument("--guard-protocol",required=True)
 a=p.parse_args();result=preflight() if a.command=="preflight" else prepare(a.anchor,a.output,a.guard_protocol)
 print(json.dumps(result,indent=2,allow_nan=False))
if __name__=="__main__":main()
