"""Closed resource-only continuations for domain/spectral; no signals or evolution."""
from __future__ import annotations
import argparse,copy,hashlib,importlib.util,json,os,shutil,sys
from pathlib import Path
import numpy as np
ROOT=Path("/home/jdiaz/pysnspd");S=Path("/home/jdiaz/scratch")
FINE_SHA="51ad9818845fb89586e7d0fac2f00798fe7d3affa528ee32f3ce1c0f4477e938"
FINE_ADAPTER_SHA="283eb0067e881d586751b922fa9cc38f6849b658f5d91174d9d58aef98d3dc12"
ACTION="EXPLICIT_DECLARED_CONTROL_RESOURCE_ONLY_TRANSITION_PREPARED_NO_EVOLUTION_NO_SIGNAL"
ORIGIN="explicit_declared_control_resource_only_transition"
DOM="stage5_effective_L480_h5_775_control_trigger_v6_guard2b_20261011"
SPEC="stage5_effective_L320_h5_775_spectral_trigger_v6_guard2b_20261011"
FINE="stage5_effective_L320_h35_775_adaptive_v6_guard2b_resources4_20261011"
N1550="stage5_effective_L320_h5_1550_adaptive_v6_20261011"
WIDE="stage5_effective_L320_h5_775_sigma_wide_v6_control20_data_20261011"
TEMP="stage5_effective_L320_h5_775_temporal_adaptive_v6_guard2b_20261011"
CASES={
 "domain_2to4":dict(kind="domain",parent=DOM,child=DOM.replace("_20261011","_resources4_20261011"),old=2,new=4,old_cpus=[10,11],cpus=[10,11,12,13],cpu_start=10,prior=None,prerequisites={"fine_moved":FINE}),
 "spectral_2to4":dict(kind="spectral",parent=SPEC,child=SPEC.replace("_20261011","_resources4_20261011"),old=2,new=4,old_cpus=[6,7],cpus=[4,5,6,7],cpu_start=4,prior=None,prerequisites={"1550_terminal":N1550,"wide_terminal":WIDE}),
 "spectral_4to8":dict(kind="spectral",parent=SPEC.replace("_20261011","_resources4_20261011"),child=SPEC.replace("_20261011","_resources8_20261011"),old=4,new=8,old_cpus=[4,5,6,7],cpus=list(range(4,12)),cpu_start=4,prior="spectral_2to4",prerequisites={"1550_terminal":N1550,"wide_terminal":WIDE,"temporal_terminal":TEMP,"domain_terminal":DOM.replace("_20261011","_resources4_20261011")})
}
def case(name):
 if name not in CASES:raise ValueError("Only the three explicitly declared control transitions are supported")
 return copy.deepcopy(CASES[name])
def dependencies():
 from sandbox.stage5_effective_validation import fine_resource_transition as f,fine_resource_guard_adapter as fa
 if f.digest(f.__file__)!=FINE_SHA or f.digest(fa.__file__)!=FINE_ADAPTER_SHA:raise ValueError("Reviewed fine bridge changed")
 g,guard,audit=f.dependencies();return f,fa,g,guard,audit
def parent_contract(name,manifest,m,arrays):
 c=case(name);f,_,_,_,_=dependencies()
 if manifest.get("status")!="INTERRUPTED" or manifest.get("run_complete") is not False:raise ValueError("Stable interrupted unfinished parent required")
 expected=dict(threads=c["old"],local_energy_tolerance=2.5e-6,step_tolerance=2e-7,dt_ps=.02,
  wavelength_nm=775.,spacing_nm=5.,duration_ps=None,
  length_nm=480. if c["kind"]=="domain" else 320.,xi_cutoff=32. if c["kind"]=="domain" else 64.,
  xi_order=6 if c["kind"]=="domain" else 8,phonon_nodes=96 if c["kind"]=="domain" else 144,
  transport_order=48 if c["kind"]=="domain" else 96)
 i=m["parameter_identity"]
 if any(isinstance(i.get(k),bool) or i.get(k)!=v for k,v in expected.items()):raise ValueError("Declared physical, solver, guard or parent resource identity changed")
 if not f.exact_json(manifest["parameter_identity"],i) or set(arrays)!=f.checkpoint_namespace() or len(m["source_pins_sha256"])!=13:raise ValueError("Actual28-key/27-numeric schema and13 frozen sources required")
 state=m["controller"]["state"]
 if state.get("status")!="ACTIVE" or any(state.get(k) is not None for k in f.EVENTS) or state.get("rejected_candidate_count")!=0:raise ValueError("Only an active pre-trigger parent is eligible")
 anchor=m["guard_active_since_time_ps"]
 if not 22.185<=anchor<=22.435 or m["time_ps"]<anchor:raise ValueError("Actual bounded guard2b anchor required")
 r=manifest["resources"]
 if any(r.get(k)!=v for k,v in dict(actual_threads=c["old"],requested_threads=c["old"],actual_affinity_cpus=c["old_cpus"]).items()):raise ValueError("Observed parent resources differ from the closed plan")
 return True
def transition_arrays(name,before,after):
 c=case(name);f,*_=dependencies()
 if set(before)!=f.checkpoint_namespace() or set(after)!=f.checkpoint_namespace():raise ValueError("Exact producer28-key namespace required")
 old=json.loads(str(before["metadata_json"].item()));new=json.loads(str(after["metadata_json"].item()))
 wanted=copy.deepcopy(old);wanted["parameter_identity"]["threads"]=c["new"]
 if old["parameter_identity"].get("threads")!=c["old"] or not f.exact_json(new,wanted):raise ValueError("Only the explicitly declared thread count may change")
 if any(not f.same_bytes(v,after[k]) for k,v in before.items() if k!="metadata_json"):raise ValueError("One of27numeric arrays changed")
 return old,new
def resource_segments(name,parent_proof,m):
 c=case(name);rows=copy.deepcopy(parent_proof.get("complete_resource_and_numerical_segments",parent_proof["segments"]))
 old=rows[-1]
 if old["threads"]!=c["old"] or old["local_energy_tolerance"]!=2.5e-6 or old["guard_active_since_time_ps"]!=m["guard_active_since_time_ps"]:raise ValueError("Parent authenticated resource/guard timeline differs")
 nxt=copy.deepcopy(old);old["end_time_ps"]=m["time_ps"]
 nxt.update(start_time_ps=m["time_ps"],end_time_ps=None,threads=c["new"],origin=ORIGIN,
  resource_case=name,resource_affinity_cpus=c["cpus"],resource_cpu_start=c["cpu_start"],guard_window_preserved=True)
 rows.append(nxt);return rows
def guard_projection(parent_proof):
 # Project ONLY an already authenticated original guard law. Full resource rows stay separately visible.
 rows=copy.deepcopy(parent_proof["segments"])
 for row in rows:
  for key in ("threads","resource_affinity_cpus","resource_cpu_start"):row.pop(key,None)
 return rows
def prerequisite_contract(name,manifests,owners):
 c=case(name)
 if set(manifests)!=set(c["prerequisites"]):raise ValueError("All explicit prerequisite manifests required")
 for key,m in manifests.items():
  path=str(S/c["prerequisites"][key])
  if any(x["output_dir"]==path for x in owners) and key!="fine_moved":raise ValueError("A required completed case still has a process owner")
  if key=="fine_moved":
   r=m["resources"]
   if m.get("status")=="EXPLICIT_RESOURCE_TRANSITION_PREPARED" or any(r.get(k)!=v for k,v in dict(actual_threads=4,requested_threads=4,actual_affinity_cpus=[0,1,2,3]).items()):raise ValueError("Fine must actually have moved toNum4 on0..3")
  else:
   allowed=("DEVELOPMENT_DURATION_REACHED",) if key=="wide_terminal" else ("TRIGGER_MARGIN_COMPLETE","TIME_CENSORED_NO_TRIGGER")
   if m.get("run_complete") is not True or m.get("status") not in allowed:raise ValueError("Required case has not actually completed its observation")
 return True
def resource_preflight(name,manifests,owners,inventory,budget):
 c=case(name);prerequisite_contract(name,manifests,owners)
 if any(x["output_dir"] in {str(S/c["parent"]),str(S/c["child"])} for x in owners):raise ValueError("Parent or child still owns an active process")
 groups=inventory["physical_core_groups"];allowed=set(budget["selected_affinity_cpus"])
 admitted=[g for g in groups if set(g)<=allowed]
 if len(admitted)!=14 or any(set(g)&{14,15,30,31} for g in admitted):raise ValueError("Full host budget must reserve complete physical14/15")
 selected=[g[0] for g in admitted[c["cpu_start"]:c["cpu_start"]+c["new"]]]
 if selected!=c["cpus"]:raise ValueError("Actual runner cpu-start selection differs")
 occupied={j for x in owners for j,g in enumerate(groups) if set(g)&set(x["affinity"])}
 target={j for j,g in enumerate(groups) if set(g)&set(c["cpus"])}
 if len(target)!=c["new"] or occupied&target or occupied&{14,15} or len(occupied|target)>14 or len(occupied|target)>.9*len(groups):raise ValueError("Target is occupied or global90percent budget/reserve fails")
 return dict(case=name,target_affinity_cpus=c["cpus"],cpu_start=c["cpu_start"],requested_threads=c["new"],
  active_pipeline_owners=owners,physical_cores_after=len(occupied|target),reserved_physical_cores=[14,15],
  scope="Point-in-time jdiaz pipeline-owner inventory only; root must recheck ALL host workloads immediately before launch; no CPU reservation")
def preflight(name,context,audit):
 f,fa,*_=dependencies()
 from sandbox.stage4_core.parallel_runtime import linux_resources,resource_budget
 c=case(name);entries={k:f.entry(f.owned(S/v/"manifest.json")) for k,v in c["prerequisites"].items()}
 manifests={k:f.read(f.bound(e)) for k,e in entries.items()}
 inv=linux_resources();budget=resource_budget(inv,.9,.9,max_workers=1)
 proof=resource_preflight(name,manifests,f.process_owners(),inv,budget)
 proof.update(prerequisite_manifests=entries,observed_prerequisite_source_paths={k:str(S/v/"manifest.json") for k,v in c["prerequisites"].items()},runtime_inventory=inv,scheduling_budget=budget)
 if name=="domain_2to4":
  rp=S/FINE/"provenance/transfer_receipt.json";r=f.read(f.owned(rp));h=dependencies()[2].h
  if r.get("helper_sha256")!=FINE_SHA or r.get("action")!=f.ACTION or r["output_dir"]!=str(S/FINE):raise ValueError("Reviewed fine resource-only transfer required")
  h.lineage.verify_receipt(rp,r["pins_after"]);_,proof["fine_resource_ancestry"]=h.lineage.resolve_receipt(rp,r["pins_after"])
  proof["fine_resource_receipt"]=f.entry(rp)
 return proof
NUM8_HELPER_SHA="85fab053549194a658835a2933416dcecac335048842d0d8ff7d7d07d9fc1fd6"
def thread_evidence(name,pins,num8=None,num8_approval=None):
 c=case(name);f,*_=dependencies();base=f.certificate(pins)
 if c["new"]==4:
  if num8 is not None or num8_approval is not None:raise ValueError("Num8 evidence is not an alias for aNum4 transition")
  return base
 if num8 is None or num8_approval is None:raise ValueError("Num8 remains blocked until real brief oracle evidence and separate root approval exist")
 p=Path(__file__).with_name("num8_resource_equivalence.py")
 if f.digest(p)!=NUM8_HELPER_SHA:raise ValueError("Reviewed separateNum8 verifier changed")
 spec=importlib.util.spec_from_file_location("_reviewed_num8_resource_equivalence",p);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
 return dict(existing_num124=base,num8=mod.verify(num8,pins,num8_approval),future_cross_thread_bitwise_identity_claimed=False)
def parent_proof(name,run,context,audit,prior_receipt=None):
 c=case(name);_,_,_,guard,_=dependencies()
 if c["prior"] is None:
  if prior_receipt is not None:raise ValueError("First resource transition must start from original bounded-guard parent")
  return guard.verify_run(run,context[0]["controls"][c["kind"]],c["kind"],context)
 if prior_receipt is None:raise ValueError("Second spectral resource transition requires its authenticatedNum2->4 parent receipt")
 return verify_resource_run(c["prior"],run,prior_receipt,context,audit)
def prepared_manifest(parent,identity,plan):
 out=copy.deepcopy(parent);out.update(parameter_identity=identity,status="EXPLICIT_RESOURCE_TRANSITION_PREPARED",
  explicit_transfer_receipt="provenance/transfer_receipt.json",explicit_resource_transition=plan);return out
def copy_prepared_state_trees(anchor,output,geometry,parent_dir):
 """Copy mandatory state trees; omit generated meshes only for authenticated archived DC."""
 f,fa,_,_,_=dependencies()
 anchor=Path(anchor);output=Path(output);parent_dir=Path(parent_dir)
 for rel in ("runtime_sources","dc"):
  source=anchor/rel
  if not source.is_dir() or source.is_symlink():raise ValueError("Mandatory captured state tree missing or linked: "+rel)
  shutil.copytree(source,output/rel)
 meshes=anchor/"meshes"
 if meshes.exists() or meshes.is_symlink():
  if not meshes.is_dir() or meshes.is_symlink():raise ValueError("Captured generated meshes must be a real directory")
  shutil.copytree(meshes,output/"meshes")
  return
 source=fa.mesh_path_without_symlinks(geometry.get("source_path",""))
 parent_dir=fa.mesh_path_without_symlinks(parent_dir)
 if not source.is_relative_to(parent_dir/fa.ARCHIVED_DC_MESH_ROOT):
  raise ValueError("Missing meshes is permitted only for the exact archived DC mesh root")
 relative=source.relative_to(parent_dir);captured=fa.mesh_path_without_symlinks(anchor/relative)
 declared=geometry.get("mesh_sha256")
 archive=f.read(anchor/"runtime_sources/archive_manifest.json")
 entry=archive.get("package_sources_and_small_geometry",{}).get(str(relative.relative_to("runtime_sources")))
 if not isinstance(declared,str) or not isinstance(entry,dict) or entry.get("sha256")!=declared or f.digest(captured)!=declared:
  raise ValueError("Captured archived mesh must match its geometry and archive-manifest SHA256")
 fa.geometry_relocation(geometry,geometry,parent_dir,output,prepared=True)
def prepare(name,anchor,output,protocol_path,*,prior_receipt=None,num8=None,num8_approval=None):
 c=case(name);f,fa,g,guard,audit=dependencies();h=g.h
 if Path(output).resolve()!=S/c["child"]:raise ValueError("Only the closed declared child path is allowed")
 anchor=f.owned(anchor);cr=f.read(anchor/"capture_receipt.json")
 if cr["parent_output_dir"]!=str(S/c["parent"]):raise ValueError("Wrong fixed parent")
 payload=(anchor/"checkpoint_original.npz").read_bytes();mb=(anchor/"manifest_original.json").read_bytes()
 if hashlib.sha256(payload).hexdigest()!=cr["checkpoint_bytes_sha256"] or hashlib.sha256(mb).hexdigest()!=cr["manifest_sha256"]:raise ValueError("Immutable captured bytes changed")
 for rel,want in cr["physical_file_sha256"].items():
  if f.digest(f.owned(anchor/rel))!=want:raise ValueError("Captured physical evidence changed")
 parentdir=S/c["parent"]
 if (parentdir/"checkpoint_latest.npz").read_bytes()!=payload or (parentdir/"manifest.json").read_bytes()!=mb:raise ValueError("Stopped parent changed since capture")
 before=h.decode(payload);parent=json.loads(mb);m=h.audit_arrays(before,parent);parent_contract(name,parent,m,before);g.eligible(parent,m,before)
 ctx=guard.load_protocol(protocol_path);resources=preflight(name,ctx,audit)
 if any(f.read(f.bound(e))["source_pins_sha256"]!=m["source_pins_sha256"] for e in resources["prerequisite_manifests"].values()):raise ValueError("Prerequisite frozen13runtime pins differ")
 os.sched_setaffinity(0,{4})
 runtime=h.audit_runtime(anchor/"runtime_sources",m["source_pins_sha256"]);threads=thread_evidence(name,m["source_pins_sha256"],num8,num8_approval)
 for old,new in (("manifest_original.json","manifest.json"),("checkpoint_original.npz","checkpoint_latest.npz")):shutil.copy2(anchor/old,anchor/new)
 for rel in ("trajectory.npz","snapshots.npz"):shutil.copy2(f.owned(parentdir/rel),anchor/rel)
 _,audit_view=f.audit_manifest_view(anchor,parent,write_json=h.writejson)
 original=audit.audit_run("resource_parent",anchor)
 if not f.exact_json(original["manifest"]["geometry"],parent["geometry"]):raise ValueError("Audited geometry differs from authenticated original manifest")
 pp=parent_proof(name,original,ctx,audit,prior_receipt)
 if any(not f.true_check(v) for k,v in original["summary"]["checks"].items() if k not in {"complete_declared_window","successful_saved_window_status"}):raise ValueError("Original physical audit failed")
 if (parentdir/"checkpoint_latest.npz").read_bytes()!=payload or (parentdir/"manifest.json").read_bytes()!=mb:raise ValueError("Stopped original changed during preparation")
 identity=copy.deepcopy(m["parameter_identity"]);identity["threads"]=c["new"]
 segments=resource_segments(name,pp,m);h.lineage.segments_valid(segments)
 plan=dict(requested_threads=c["new"],target_affinity_cpus=c["cpus"],cpu_start=c["cpu_start"],preparation_only=True,execution_started=False,previous_observed_resources=parent["resources"])
 output=h.newdir(output);prov=output/"provenance";prov.mkdir()
 for key,e in resources["prerequisite_manifests"].items():
  data=f.bound(e).read_bytes();saved=prov/("prerequisite_"+key+"_manifest.json");saved.write_bytes(data)
  resources["prerequisite_manifests"][key]=f.entry(saved)
 copy_prepared_state_trees(anchor,output,original["manifest"]["geometry"],parentdir)
 for rel in ("material.npz","material.json","source_initial_state.npz","trajectory.npz","snapshots.npz"):shutil.copy2(anchor/rel,output/rel)
 new=copy.deepcopy(m);new["parameter_identity"]=identity;initial=copy.deepcopy(before);initial["metadata_json"]=np.array(json.dumps(new,allow_nan=False))
 np.savez_compressed(output/"checkpoint_latest.npz",**initial);manifest=prepared_manifest(parent,identity,plan);h.writejson(output/"manifest.json",manifest)
 transition_arrays(name,before,h.decode((output/"checkpoint_latest.npz").read_bytes()));h.audit_arrays(initial,manifest)
 for old,newfile in (("checkpoint_original.npz","parent_checkpoint.npz"),("manifest_original.json","parent_manifest.json"),("capture_receipt.json","capture_receipt.json")):shutil.copy2(anchor/old,prov/newfile)
 shutil.copy2(anchor/"provenance/parent_transfer_receipt.json",prov/"parent_transfer_receipt.json")
 for old,newfile in (("checkpoint_latest.npz","initial_resource_checkpoint.npz"),("manifest.json","initial_resource_manifest.json")):shutil.copy2(output/old,prov/newfile)
 r=dict(schema="pysnspd.explicit_checkpoint_transfer.v3",action=ACTION,case=name,helper_sha256=f.digest(__file__),admission=False,
  parent=dict(output_dir=str(parentdir),anchor=str(anchor),checkpoint_bytes_sha256=cr["checkpoint_bytes_sha256"],manifest_sha256=cr["manifest_sha256"],runtime_archive_manifest_sha256=runtime["manifest_sha256"]),
  output_dir=str(output),pins_before=m["source_pins_sha256"],pins_after=m["source_pins_sha256"],
  state_arrays_sha256_before={k:h.arrayhash(before[k]) for k in h.STATE},state_arrays_sha256_after={k:h.arrayhash(initial[k]) for k in h.STATE},
  conserved_metadata_before={k:m[k] for k in h.CONSERVED},conserved_metadata_after={k:new[k] for k in h.CONSERVED},
  preserved_array_hashes={k:h.arrayhash(v) for k,v in before.items() if k!="metadata_json"},
  checkpoint_total_keys=28,checkpoint_numeric_arrays=27,checkpoint_namespace=sorted(f.checkpoint_namespace()),all_metadata_except_threads_preserved=True,
  frontend_changes={},implementation_changes={},horizon_change=None,energy_guard_change=None,RHS_sources_identical=True,
  guard_ledger_before={k:m[k] for k in h.GUARD_KEYS},guard_ledger_after={k:new[k] for k in h.GUARD_KEYS},
  next_step_ps_before=m["next_step_ps"],next_step_ps_after=new["next_step_ps"],
  global_energy_error_before=g.global_error_at_anchor(m,before),global_energy_error_after=g.global_error_at_anchor(new,initial),
  resources_before=parent["resources"],resources_after=plan,numerical_segments=segments,guard_policy_projection=guard_projection(pp),
  resource_preflight=resources,parent_guard_entry=ctx[0]["controls"][c["kind"]],parent_guard_verification=pp,parent_guard_protocol=f.entry(protocol_path),
  prior_resource_receipt=None if prior_receipt is None else f.entry(prior_receipt),existing_thread_evidence=threads,
  initial_checkpoint=f.entry(prov/"initial_resource_checkpoint.npz"),initial_manifest=f.entry(prov/"initial_resource_manifest.json"),
  parent_audit_view_manifest=audit_view,
  parent_original_audit_files={rel:f.entry(anchor/rel) for rel in ("manifest.json","checkpoint_latest.npz","trajectory.npz","snapshots.npz")},
  parent_runtime_archive_manifest=f.entry(anchor/"runtime_sources/archive_manifest.json"),new_checkpoint_bytes_sha256=f.digest(output/"checkpoint_latest.npz"),
  resume_argv=h.resume_argv(identity,output,c["cpu_start"]),state_or_heat_correction=False,global_energy_origin_reset=False,guard_window_reset=False,future_cross_thread_bitwise_identity_claimed=False)
 h.writejson(prov/"transfer_receipt.json",r);h.lineage.verify_receipt(prov/"transfer_receipt.json",r["pins_after"])
 resolved,_=h.lineage.resolve_receipt(prov/"transfer_receipt.json",r["pins_after"])
 if resolved!=segments:raise ValueError("Complete original recursive lineage differs")
 return r
HISTORICAL_DOMAIN_RECEIPT= S/"stage5_effective_L480_h5_775_control_trigger_v6_guard2b_resources4_20261011/provenance/transfer_receipt.json"
HISTORICAL_DOMAIN_RECEIPT_SHA="29d117bf39fb7b13d52f4a1ade6ba0bc2483e8d76eaa98f1942e7d9857c3cd48"
HISTORICAL_DOMAIN_HELPER_SHA="0458f0141d7005578d0ab98fb2e1b963f79b0c34769197a585d5651c4ebf5118"
HISTORICAL_DOMAIN_HELPER_SOURCE=Path(__file__).parent/"resource_transition_history/20261011_control_before_archived_mesh.py"
def verify_receipt_helper_origin(name,receipt_path,receipt):
 """Authenticate current origin or the one exact pre-fix Domain4 receipt.

 The historical receipt/source remain unchanged; this is an audit compatibility
 transformation, not a claim that the new helper prepared that earlier state.
 """
 digest=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
 if receipt.get("helper_sha256")==digest(__file__):return True
 if (name!="domain_2to4" or Path(receipt_path)!=HISTORICAL_DOMAIN_RECEIPT
  or receipt.get("helper_sha256")!=HISTORICAL_DOMAIN_HELPER_SHA
  or digest(receipt_path)!=HISTORICAL_DOMAIN_RECEIPT_SHA
  or digest(HISTORICAL_DOMAIN_HELPER_SOURCE)!=HISTORICAL_DOMAIN_HELPER_SHA):
  raise ValueError("Only the exact authenticated historical Domain4 helper/receipt is compatible")
 return True
def verify_resource_run(name,run,receipt_path,context,audit,*,prepared=False):
 c=case(name);f,fa,g,guard,_=dependencies();h=g.h
 rp=f.owned(receipt_path);r=f.read(rp)
 if r.get("action")!=ACTION or r.get("case")!=name or r.get("admission") is not False:raise ValueError("Explicit reviewed closed resource receipt required")
 verify_receipt_helper_origin(name,rp,r)
 if r["output_dir"]!=str(S/c["child"]) or r["parent"]["output_dir"]!=str(S/c["parent"]):raise ValueError("Closed parent/child paths differ")
 original_context=guard.load_protocol(f.bound(r["parent_guard_protocol"]))
 if r["parent_guard_entry"]!=context[0]["controls"][c["kind"]] or original_context[0]["controls"][c["kind"]]!=context[0]["controls"][c["kind"]] or original_context[0]["nominal"]!=context[0]["nominal"]:raise ValueError("Original relevant guard entries changed")
 if any(r[k]!={} for k in ("frontend_changes","implementation_changes")) or any(r[k] is not None for k in ("horizon_change","energy_guard_change")):raise ValueError("Resource tail hides another change")
 if any(r[k] is not False for k in ("state_or_heat_correction","global_energy_origin_reset","guard_window_reset","future_cross_thread_bitwise_identity_claimed")):raise ValueError("State, ledger, guard reset or future bitwise claim prohibited")
 anchor=f.owned(r["parent"]["anchor"])
 if set(r["parent_original_audit_files"])!={"manifest.json","checkpoint_latest.npz","trajectory.npz","snapshots.npz"}:raise ValueError("Complete original evidence required")
 for rel,e in r["parent_original_audit_files"].items():
  if f.bound(e)!=anchor/rel:raise ValueError("Immutable original path differs")
 if f.bound(r["parent_runtime_archive_manifest"])!=anchor/"runtime_sources/archive_manifest.json":raise ValueError("Archive path differs")
 original_manifest=f.read(anchor/"manifest_original.json")
 _,audit_view=f.audit_manifest_view(anchor,original_manifest)
 if not f.exact_json(audit_view,r.get("parent_audit_view_manifest")):raise ValueError("Authenticated derived parent audit view differs")
 payload=(anchor/"checkpoint_latest.npz").read_bytes();before=h.decode(payload);parent=audit.audit_run("original_resource_parent",anchor);pm=h.audit_arrays(before,parent["manifest"])
 if hashlib.sha256(payload).hexdigest()!=r["parent"]["checkpoint_bytes_sha256"] or f.digest(anchor/"manifest_original.json")!=r["parent"]["manifest_sha256"]:raise ValueError("Original full parent bytes changed")
 parent_contract(name,parent["manifest"],pm,before)
 if any(not f.true_check(v) for k,v in parent["summary"]["checks"].items() if k not in {"complete_declared_window","successful_saved_window_status"}):raise ValueError("Independent original parent physical/guard audit failed")
 prior=None if r["prior_resource_receipt"] is None else f.bound(r["prior_resource_receipt"])
 pp=parent_proof(name,parent,context,audit,prior)
 if pp!=r["parent_guard_verification"]:raise ValueError("Original recursive guard proof changed")
 initial=h.decode(f.bound(r["initial_checkpoint"]).read_bytes());im=f.read(f.bound(r["initial_manifest"]));old,new=transition_arrays(name,before,initial)
 if im!=prepared_manifest(original_manifest,new["parameter_identity"],r["resources_after"]) or r["new_checkpoint_bytes_sha256"]!=r["initial_checkpoint"]["sha256"]:raise ValueError("Prepared initial manifest/checkpoint changed")
 if r["pins_before"]!=r["pins_after"] or r["pins_after"]!=pm["source_pins_sha256"]:raise ValueError("Frozen13runtime pins changed")
 num8=None if c["new"]==4 else f.bound(r["existing_thread_evidence"]["num8"]["certificate"])
 approval=None if c["new"]==4 else f.bound(r["existing_thread_evidence"]["num8"]["root_approval"])
 if thread_evidence(name,r["pins_after"],num8,approval)!=r["existing_thread_evidence"]:raise ValueError("Thread evidence changed")
 expected_plan=dict(requested_threads=c["new"],target_affinity_cpus=c["cpus"],cpu_start=c["cpu_start"],preparation_only=True,execution_started=False,previous_observed_resources=parent["manifest"]["resources"])
 if r["resources_after"]!=expected_plan or r["resources_before"]!=parent["manifest"]["resources"]:raise ValueError("Resources differ from closed plan")
 pf=r["resource_preflight"];mans={k:f.read(f.bound(e)) for k,e in pf["prerequisite_manifests"].items()}
 expected=resource_preflight(name,mans,pf["active_pipeline_owners"],pf["runtime_inventory"],pf["scheduling_budget"])
 expected.update(prerequisite_manifests=pf["prerequisite_manifests"],observed_prerequisite_source_paths={k:str(S/v/"manifest.json") for k,v in c["prerequisites"].items()},runtime_inventory=pf["runtime_inventory"],scheduling_budget=pf["scheduling_budget"])
 for key,e in pf["prerequisite_manifests"].items():
  if f.bound(e)!=S/c["child"]/"provenance"/("prerequisite_"+key+"_manifest.json"):raise ValueError("Prerequisite manifest must remain the immutable observed copy")
  if mans[key]["source_pins_sha256"]!=r["pins_after"]:raise ValueError("Prerequisite13runtime pins differ")
 if name=="domain_2to4":
  fr=f.bound(pf["fine_resource_receipt"]);fine=f.read(fr)
  if fine.get("helper_sha256")!=FINE_SHA or fine.get("action")!=f.ACTION or fine["output_dir"]!=str(S/FINE):raise ValueError("Reviewed fine transition required")
  h.lineage.verify_receipt(fr,fine["pins_after"]);_,expected["fine_resource_ancestry"]=h.lineage.resolve_receipt(fr,fine["pins_after"]);expected["fine_resource_receipt"]=pf["fine_resource_receipt"]
 if pf!=expected:raise ValueError("Prerequisite/occupancy proof changed")
 segments=resource_segments(name,pp,pm);resolved,_=h.lineage.resolve_receipt(rp,r["pins_after"])
 if segments!=resolved or segments!=r["numerical_segments"] or guard_projection(pp)!=r["guard_policy_projection"]:raise ValueError("Recursive full resource timeline or guard projection changed")
 for label,keys in (("conserved_metadata",h.CONSERVED),("guard_ledger",h.GUARD_KEYS)):
  wanted={k:old[k] for k in keys}
  if r[label+"_before"]!=wanted or r[label+"_after"]!=wanted:raise ValueError("Counter/ledger changed")
 if r["preserved_array_hashes"]!={k:h.arrayhash(v) for k,v in before.items() if k!="metadata_json"} or r["checkpoint_namespace"]!=sorted(f.checkpoint_namespace()) or r["checkpoint_total_keys"]!=28 or r["checkpoint_numeric_arrays"]!=27 or r["all_metadata_except_threads_preserved"] is not True:raise ValueError("All27numeric arrays plus metadata must be authenticated")
 for side in ("before","after"):
  if r["next_step_ps_"+side]!=old["next_step_ps"] or r["global_energy_error_"+side]!=g.global_error_at_anchor(old,before):raise ValueError("Adaptive next step or energy origin changed")
 path=f.owned(run["summary"]["path"]);live_payload=(path/"checkpoint_latest.npz").read_bytes();live=h.decode(live_payload);lm=f.live_prefix(initial,live);h.audit_arrays(live,run["manifest"])
 rr=run["manifest"]["resources"]
 if prepared:
  if run["manifest"]!=im or hashlib.sha256(live_payload).hexdigest()!=r["initial_checkpoint"]["sha256"] or rr!=r["resources_before"]:raise ValueError("Prepared child must still be the exact unlaunched immutable initial state")
 elif any(rr.get(k)!=v for k,v in dict(actual_threads=c["new"],requested_threads=c["new"],actual_affinity_cpus=c["cpus"],actual_processes=1,global_resource_ceiling_fraction=.9).items()):raise ValueError("Actual observed child resources differ")
 fa.geometry_relocation(original_manifest["geometry"],run["manifest"]["geometry"],S/c["parent"],path,prepared=prepared)
 for key in ("physicalparams","source","energy_reference","material","numerical"):
  if run["manifest"].get(key)!=parent["manifest"].get(key):raise ValueError("Child physical/solver identity changed")
 for prefix,table in (("trajectory__",run["history"]),("snapshots__",run["snap"])):
  for key,value in live.items():
   if key.startswith(prefix) and (key[len(prefix):] not in table or not guard.same_array(value,table[key[len(prefix):]])):raise ValueError("Audited tables differ from complete original checkpoint")
 guard.close(run["summary"]["final_time_ps"],lm["time_ps"],"Resource child original clock");fa.audit_energy_contract(run["summary"])
 current=g.validate_old_guard_history(run["manifest"],lm,live);h.audit_runtime(path/"runtime_sources",r["pins_after"])
 return dict(role=c["kind"],segments=r["guard_policy_projection"],anchor_time_ps=pp["anchor_time_ps"],guard_policy_projection_only=True,
  complete_resource_and_numerical_segments=segments,resource_transition_time_ps=pm["time_ps"],resource_transition_receipt=f.entry(rp),
  parent_guard_proof=pp,current_checkpoint_sha256=hashlib.sha256(live_payload).hexdigest(),current_guard_independently_verified=current,
  exact_saved_prefix=True,global_origin_preserved=True,future_cross_thread_bitwise_identity_claimed=False,admission=False)
def compare_with_resource_tail(kind,reference,candidate,protocol_path,audit,*,resource_receipt_path):
 f,_,_,guard,_=dependencies();r=f.read(f.owned(resource_receipt_path));name=r.get("case");c=case(name)
 if kind!=c["kind"]:raise ValueError("Explicit selector may not alias mesh, temporal or another axis")
 if f.digest(audit.__file__)!=f.BASELINE_SHA:raise ValueError("Frozen baseline changed")
 ctx=guard.load_protocol(protocol_path);ref=guard.verify_run(reference,ctx[0]["nominal"],"nominal",ctx);cand=verify_resource_run(name,candidate,resource_receipt_path,ctx,audit)
 result,curves,gapcurves=audit.compare_runs(kind,reference,candidate)
 if reference["manifest"]["source_pins_sha256"]!=candidate["manifest"]["source_pins_sha256"]:raise ValueError("Frozen13-source model differs")
 policy=guard.compare_guard_law(ref["segments"],cand["segments"],kind,ctx[3],result["common_trajectory_interval_ps"])
 baseline=copy.deepcopy(result);forbidden=copy.deepcopy(result["forbidden_parameter_changes"]);forbidden.pop("local_energy_tolerance",None)
 result["forbidden_parameter_changes"]=forbidden
 for key in ("checks","legacy_nodal_policy_checks"):result[key].update(same_nonrefined_parameters=not forbidden,verified_bounded_prospective_guard_policy=True,verified_explicit_resource_only_continuation=True)
 result["passed"]=all(result["checks"].values());result["legacy_nodal_policy_passed"]=all(result["legacy_nodal_policy_checks"].values());result["baseline_comparison_without_guard_protocol"]=baseline
 result["declared_numerical_control"]=dict(schema=guard.SCHEMA,protocol_path=str(guard.under(protocol_path)),protocol_sha256=guard.digest(protocol_path),
  prospective_decision_sha256=guard.DECLARATION_SHA,final_nominal_anchor_sha256=guard.ANCHOR_SHA,comparison_label=ctx[1]["comparison_label"],
  guard_law=policy,reference=ref,candidate=cand,observable_limits_unchanged=True,global_energy_origin_reset=False,
  original_completion_flags_preserved=True,stage5_admission_inferred=False,resource_only_continuation_explicit=True,
  resource_adapter_sha256=f.digest(__file__),resource_case=name)
 return result,curves,gapcurves
def execution_argv(name,resume_args):
 c=case(name);f,*_=dependencies();output=S/c["child"];runner=output/"runtime_sources/pipelines/05_effective_photon_run.py"
 if not isinstance(resume_args,list) or not all(isinstance(x,str) for x in resume_args) or resume_args[:5]!=["--output-dir",str(output),"--resume","--cpu-start",str(c["cpu_start"])]:raise ValueError("Reviewed resume_argv contains arguments only; exact child/affinity prefix required")
 if f.digest(runner)!=f.RUNNER_SHA:raise ValueError("Frozen child runner source changed")
 return [sys.executable,"-u",str(runner),*resume_args]
def launch_preflight(name,receipt_path,protocol_path,approval_path):
 c=case(name);f,fa,g,guard,audit=dependencies();rp=f.owned(receipt_path);r=f.read(rp)
 approval=f.read(f.owned(approval_path))
 if approval.get("status")!="ROOT_APPROVED_CONTROL_RESOURCE_TRANSFER" or approval.get("receipt_sha256")!=f.digest(rp) or approval.get("case")!=name or approval.get("global_host_occupancy_rechecked") is not True or approval.get("stage5") is not False:raise ValueError("Root-reviewed real transfer and full host occupancy approval required")
 ctx=guard.load_protocol(protocol_path);run=audit.audit_run("unlaunched_resource_child",S/c["child"])
 verify_resource_run(name,run,rp,ctx,audit,prepared=True);preflight(name,ctx,audit)
 meta=g.h.decode((S/c["child"]/"checkpoint_latest.npz").read_bytes());m=json.loads(str(meta["metadata_json"].item()))
 argv=g.h.resume_argv(m["parameter_identity"],S/c["child"],c["cpu_start"])
 if argv!=r["resume_argv"]:raise ValueError("Explicit resume command differs")
 return dict(case=name,resume_argv=argv,runner_argv=execution_argv(name,argv),receipt=f.entry(rp),root_approval=f.entry(approval_path),execution_started=False,stage5=False)
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument("--case",choices=sorted(CASES),required=True);sub=p.add_subparsers(dest="command",required=True)
 pre=sub.add_parser("preflight");pre.add_argument("--guard-protocol",required=True)
 launch=sub.add_parser("launch-preflight")
 for key in ("receipt","guard-protocol","root-approval"):launch.add_argument("--"+key,required=True)
 q=sub.add_parser("prepare")
 for key in ("anchor","output","guard-protocol"):q.add_argument("--"+key,required=True)
 q.add_argument("--prior-resource-receipt");q.add_argument("--num8-evidence");q.add_argument("--num8-approval")
 a=p.parse_args();f,fa,g,guard,audit=dependencies()
 if a.command=="preflight":r=preflight(a.case,guard.load_protocol(a.guard_protocol),audit)
 elif a.command=="launch-preflight":r=launch_preflight(a.case,a.receipt,a.guard_protocol,a.root_approval)
 else:r=prepare(a.case,a.anchor,a.output,a.guard_protocol,prior_receipt=a.prior_resource_receipt,num8=a.num8_evidence,num8_approval=a.num8_approval)
 print(json.dumps(r,indent=2,allow_nan=False))
if __name__=="__main__":main()
