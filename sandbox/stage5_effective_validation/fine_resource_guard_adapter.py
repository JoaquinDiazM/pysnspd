"""Strict mesh-only resource-tail adapter; original bounded guard stays unchanged."""
from __future__ import annotations
import copy,hashlib,importlib.util,json,sys
from pathlib import Path
import numpy as np
from sandbox.stage5_effective_validation import bounded_guard_protocol as guard
HELPER_SHA="51ad9818845fb89586e7d0fac2f00798fe7d3affa528ee32f3ce1c0f4477e938"
def helper():
 p=Path(__file__).with_name("fine_resource_transition.py")
 if guard.digest(p)!=HELPER_SHA:raise ValueError("Reviewed resource helper changed")
 spec=importlib.util.spec_from_file_location("_reviewed_fine_resource_transition",p)
 m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m
ARCHIVED_DC_MESH_ROOT=Path("runtime_sources/docs/implementation/stage5/prephoton_dc_20260924/meshes")
def mesh_path_without_symlinks(source):
 """Require a canonical absolute scratch path with no symlink components."""
 raw=Path(source)
 if not raw.is_absolute() or any(p.is_symlink() for p in (raw,*raw.parents)):
  raise ValueError("Mesh paths must be absolute and may not contain symlinks")
 resolved=guard.under(raw)
 if raw!=resolved:raise ValueError("Mesh paths must be canonical without parent traversal")
 return resolved
def geometry_relocation(original,current,parent_dir,child_dir,*,prepared=False):
 """Permit only exact generated or archived-DC mesh relocation on resume."""
 if not isinstance(original,dict) or not isinstance(current,dict):raise ValueError("Full original/current geometry required")
 old=copy.deepcopy(original);new=copy.deepcopy(current)
 old_source=old.pop("source_path",None);new_source=new.pop("source_path",None)
 if not isinstance(old_source,str) or not isinstance(new_source,str) or json.dumps(old,sort_keys=True,allow_nan=False)!=json.dumps(new,sort_keys=True,allow_nan=False):
  raise ValueError("Geometry fields/types may not change except source_path")
 parent_dir=mesh_path_without_symlinks(parent_dir);child_dir=mesh_path_without_symlinks(child_dir)
 old_path=mesh_path_without_symlinks(old_source)
 mesh_roots=(Path("meshes"),ARCHIVED_DC_MESH_ROOT)
 if not old_path.is_relative_to(parent_dir) or not any(old_path.is_relative_to(parent_dir/root) for root in mesh_roots):
  raise ValueError("Original mesh must remain within the declared generated or archived-DC parent mesh root")
 relative=old_path.relative_to(parent_dir);copied_path=child_dir/relative
 if mesh_path_without_symlinks(copied_path)!=copied_path or not any(copied_path.is_relative_to(child_dir/root) for root in mesh_roots):
  raise ValueError("Child mesh must be the exact copied relative mesh path")
 expected_source=old_source if prepared else str(copied_path)
 if new_source!=expected_source:raise ValueError("Only the exact parent-to-child mesh relocation is permitted")
 declared=original.get("mesh_sha256")
 if not isinstance(declared,str) or len(declared)!=64 or guard.digest(old_path)!=declared or guard.digest(copied_path)!=declared:
  raise ValueError("Original and child mesh bytes must both match the declared SHA256")
 return True
def audit_energy_contract(summary):
 energy=summary["energy"]
 if energy.get("relative_limit")!=.01 or not guard.finite_number(energy.get("maximum_relative_defect")) or energy["maximum_relative_defect"]>.01:
  raise ValueError("Independent global energy bound remains exactly1percent")
 if any(not (isinstance(summary["checks"].get(k),(bool,np.bool_)) and bool(summary["checks"][k])) for k in ("energy_relative_defect","independent_local_energy_guard")):
  raise ValueError("Independent global/local energy audit failed")
 return True

def verify_resource_run(run,receipt_path,parent_guard_entry,context,audit):
 m=helper();g,_,_=m.dependencies();h=g.h
 rp=m.owned(receipt_path);r=m.read(rp)
 if r.get("action")!=m.ACTION or r.get("helper_sha256")!=HELPER_SHA or r.get("admission") is not False:
  raise ValueError("Explicit reviewed resource-only receipt required")
 if r.get("output_dir")!=str(m.CHILD) or r["parent"]["output_dir"]!=str(m.PARENT):
  raise ValueError("Only the declared fine parent/child pair is supported")
 if r.get("parent_guard_entry")!=parent_guard_entry:raise ValueError("Original immutable guard entry changed")
 for k in ("frontend_changes","implementation_changes"):
  if r.get(k)!={}:raise ValueError("Resource transition hides source changes")
 for k in ("horizon_change","energy_guard_change"):
  if r.get(k) is not None:raise ValueError("Resource transition hides numerical/observation changes")
 for k in ("state_or_heat_correction","global_energy_origin_reset","guard_window_reset","future_cross_thread_bitwise_identity_claimed"):
  if r.get(k) is not False:raise ValueError("Forbidden state/ledger correction or bitwise claim")
 anchor=m.owned(r["parent"]["anchor"])
 if set(r.get("parent_original_audit_files",{}))!={"manifest.json","checkpoint_latest.npz","trajectory.npz","snapshots.npz"}:
  raise ValueError("Complete original immutable parent files required")
 for name,e in r["parent_original_audit_files"].items():
  if m.bound(e)!=anchor/name:raise ValueError("Parent audit evidence path differs")
 if m.bound(r["parent_runtime_archive_manifest"])!=anchor/"runtime_sources/archive_manifest.json":raise ValueError("Parent archive binding differs")
 original_manifest=m.read(anchor/"manifest_original.json")
 _,audit_view=m.audit_manifest_view(anchor,original_manifest)
 if not m.exact_json(audit_view,r.get("parent_audit_view_manifest")):raise ValueError("Authenticated derived parent audit view differs")
 parent=audit.audit_run("immutable_resource_parent",anchor)
 pp=guard.verify_run(parent,parent_guard_entry,"mesh",context)
 if pp!=r["parent_guard_verification"]:raise ValueError("Immutable parent guard proof changed")
 parent_payload=(anchor/"checkpoint_latest.npz").read_bytes();before=h.decode(parent_payload)
 if hashlib.sha256(parent_payload).hexdigest()!=r["parent"]["checkpoint_bytes_sha256"] or guard.digest(anchor/"manifest_original.json")!=r["parent"]["manifest_sha256"]:
  raise ValueError("Immutable full parent changed")
 pm=h.audit_arrays(before,parent["manifest"]);m.parent_contract(parent["manifest"],pm,before)
 if any(not m.true_check(v) for k,v in parent["summary"]["checks"].items() if k not in {"complete_declared_window","successful_saved_window_status"}):
  raise ValueError("Independent original parent physical/guard audit failed")
 initial=h.decode(m.bound(r["initial_checkpoint"]).read_bytes());im=m.read(m.bound(r["initial_manifest"]))
 old,new=m.transition_arrays(before,initial)
 if im!=m.prepared_manifest(original_manifest,new["parameter_identity"],r["resources_after"]):
  raise ValueError("Prepared child manifest altered physical, numerical or event identity")
 if r["new_checkpoint_bytes_sha256"]!=r["initial_checkpoint"]["sha256"]:
  raise ValueError("Initial prepared checkpoint SHA differs")
 if r["pins_before"]!=r["pins_after"] or r["pins_after"]!=pm["source_pins_sha256"] or len(r["pins_after"])!=13:
  raise ValueError("All13runtime sources must remain exactly equal")
 m.certificate(r["pins_after"])
 if r["existing_thread_evidence"]!=m.certificate(r["pins_after"]):raise ValueError("Existing certificate evidence changed")
 if r["resources_before"]!=parent["manifest"]["resources"]:raise ValueError("Previous observed resources changed")
 if r["resources_after"]!=dict(requested_threads=4,target_affinity_cpus=[0,1,2,3],cpu_start=0,
  preparation_only=True,execution_started=False,previous_observed_resources=parent["manifest"]["resources"]):
  raise ValueError("Only explicit Num4 affinity0..3 resource plan allowed")
 preflight=r["resource_preflight"];terminal=m.read(m.bound(preflight["nominal_terminal_manifest"]))
 expected_preflight=m.resource_preflight(terminal,preflight["active_pipeline_owners"],preflight["runtime_inventory"],preflight["scheduling_budget"])
 expected_preflight.update(nominal_terminal_manifest=preflight["nominal_terminal_manifest"],runtime_inventory=preflight["runtime_inventory"],scheduling_budget=preflight["scheduling_budget"])
 if expected_preflight!=preflight:raise ValueError("Recorded resource preflight/terminal nominal evidence differs")
 segments=m.resource_segments(pp["segments"],pm)
 resolved,_=context[5].resolve_receipt(rp,r["pins_after"])
 if resolved!=segments:raise ValueError("Recursive explicit resource lineage differs")
 if r["numerical_segments"]!=segments or r["guard_policy_projection"]!=m.guard_projection(segments,pm["time_ps"]):
  raise ValueError("Resource/numerical timeline or guard projection changed")
 for key in ("conserved_metadata","guard_ledger"):
  keys=h.CONSERVED if key=="conserved_metadata" else h.GUARD_KEYS
  expected={k:old[k] for k in keys}
  if r[key+"_before"]!=expected or r[key+"_after"]!=expected:raise ValueError("Conserved counters/ledger changed")
 expected_hashes={k:h.arrayhash(v) for k,v in before.items() if k!="metadata_json"}
 if r["preserved_array_hashes"]!=expected_hashes or r.get("checkpoint_total_keys")!=28 or r.get("checkpoint_numeric_arrays")!=27 or r.get("checkpoint_namespace")!=sorted(m.checkpoint_namespace()) or r["all_metadata_except_threads_preserved"] is not True:
  raise ValueError("All27 numeric arrays and metadata(28keys) must be explicitly checked")
 for k in ("before","after"):
  if r["next_step_ps_"+k]!=old["next_step_ps"] or r["global_energy_error_"+k]!=g.global_error_at_anchor(old,before):
   raise ValueError("Adaptive next step or global origin changed")
 live_path=m.owned(run["summary"]["path"])
 for rel in r["pins_after"]:m.owned(live_path/"runtime_sources"/rel)
 payload=(live_path/"checkpoint_latest.npz").read_bytes();live=h.decode(payload)
 lm=m.live_prefix(initial,live);h.audit_arrays(live,run["manifest"])
 if run["manifest"]["parameter_identity"]!=lm["parameter_identity"] or run["manifest"]["source_pins_sha256"]!=r["pins_after"]:
  raise ValueError("Live child numerical/source identity differs")
 geometry_relocation(original_manifest["geometry"],run["manifest"]["geometry"],m.PARENT,live_path)
 for key,value in original_manifest.items():
  if key in ("physicalparams","source","energy_reference","material","numerical") and run["manifest"].get(key)!=value:
   raise ValueError("Live child physical/solver identity changed: "+key)
 rr=run["manifest"]["resources"]
 if any(rr.get(k)!=v for k,v in dict(actual_threads=4,requested_threads=4,actual_affinity_cpus=[0,1,2,3],actual_processes=1,
  global_resource_ceiling_fraction=.9).items()):raise ValueError("Live child resources do not match the explicit Num4 plan")
 for prefix,table in (("trajectory__",run["history"]),("snapshots__",run["snap"])):
  for k,v in live.items():
   if k.startswith(prefix) and (k[len(prefix):] not in table or not guard.same_array(v,table[k[len(prefix):]])):
    raise ValueError("Audited child is not the original complete atomic checkpoint")
 guard.close(run["summary"]["final_time_ps"],lm["time_ps"],"resource child audited clock")
 audit_energy_contract(run["summary"])
 current=g.validate_old_guard_history(run["manifest"],lm,live)
 h.audit_runtime(live_path/"runtime_sources",r["pins_after"])
 return dict(role="mesh",segments=r["guard_policy_projection"],anchor_time_ps=pp["anchor_time_ps"],
  guard_policy_projection_only=True,complete_resource_and_numerical_segments=segments,
  resource_transition_time_ps=pm["time_ps"],resource_transition_receipt=m.entry(rp),
  parent_guard_proof=pp,current_checkpoint_sha256=hashlib.sha256(payload).hexdigest(),
  current_guard_independently_verified=current,exact_saved_prefix=True,global_origin_preserved=True,
  future_cross_thread_bitwise_identity_claimed=False,admission=False)
def compare_with_resource_tail(kind,reference,candidate,protocol_path,audit,*,resource_receipt_path):
 if kind!="mesh":raise ValueError("Resource-only continuation is limited to the declared fine mesh control")
 if guard.digest(audit.__file__)!=guard.AUDITOR_SHA:raise ValueError("Baseline auditor changed")
 context=guard.load_protocol(protocol_path);cfg,d,a,tstar,g,lineage=context
 ref=guard.verify_run(reference,cfg["nominal"],"nominal",context)
 cand=verify_resource_run(candidate,resource_receipt_path,cfg["controls"]["mesh"],context,audit)
 result,curves,gapcurves=audit.compare_runs(kind,reference,candidate)
 if reference["manifest"]["source_pins_sha256"]!=candidate["manifest"]["source_pins_sha256"]:raise ValueError("Frozen13-source model differs")
 policy=guard.compare_guard_law(ref["segments"],cand["segments"],kind,tstar,result["common_trajectory_interval_ps"])
 baseline=copy.deepcopy(result);forbidden=copy.deepcopy(result["forbidden_parameter_changes"]);forbidden.pop("local_energy_tolerance",None)
 result["forbidden_parameter_changes"]=forbidden
 for name in ("checks","legacy_nodal_policy_checks"):
  result[name].update(same_nonrefined_parameters=not forbidden,verified_bounded_prospective_guard_policy=True,
   verified_explicit_resource_only_continuation=True)
 result["passed"]=all(result["checks"].values());result["legacy_nodal_policy_passed"]=all(result["legacy_nodal_policy_checks"].values())
 result["baseline_comparison_without_guard_protocol"]=baseline
 result["declared_numerical_control"]=dict(schema=guard.SCHEMA,protocol_path=str(guard.under(protocol_path)),
  protocol_sha256=guard.digest(protocol_path),prospective_decision_sha256=guard.DECLARATION_SHA,
  final_nominal_anchor_sha256=guard.ANCHOR_SHA,comparison_label=d["comparison_label"],
  guard_law=policy,reference=ref,candidate=cand,observable_limits_unchanged=True,global_energy_origin_reset=False,
  original_completion_flags_preserved=True,stage5_admission_inferred=False,
  resource_only_continuation_explicit=True,resource_adapter_sha256=guard.digest(__file__))
 return result,curves,gapcurves
