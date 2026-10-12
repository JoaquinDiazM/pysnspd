"""Strict prospective local-energy guard tightening; no signals or evolution.

capture reads one coherent atomic checkpoint. tighten prepares a NEW output
directory with the same physical arrays, complete histories, source, controller,
13 source pins, resources and numerical parameters, except one strictly smaller
local_energy_tolerance. No old step is attributed to the new tolerance.
The global energy origin, signed/absolute work, residual sum and counters remain.
"""
from __future__ import annotations
import argparse,copy,importlib.util,json,math,shlex,shutil
from pathlib import Path
import numpy as np
BASE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location("_guard_window_verified",BASE/"extend_completed_control_v6.py")
w=importlib.util.module_from_spec(spec);spec.loader.exec_module(w)
h=w.h
EXPECTED_WINDOW_HELPER_SHA256="ee8a00b5972d2fcb8b93ee4214767316243ed1b50bda07709b7bec246ecd2596"
SCHEMA="pysnspd.explicit_checkpoint_transfer.v3"
WINDOW_KEYS={"guard_active_since_time_ps","max_guarded_accepted_local_energy_fraction","energy_guard_ledger_origin"}
IDENTITY_KEY="local_energy_tolerance"
def helper_integrity():
 w.helper_integrity()
 if h.filehash(BASE/"extend_completed_control_v6.py")!=EXPECTED_WINDOW_HELPER_SHA256:
  raise ValueError("Verified window helper changed")
def validate_old_guard_history(parent,m,arrays):
 """Independent persisted endpoint/work audit with adopted16/32eps bounds."""
 tol=float(m["parameter_identity"][IDENTITY_KEY]);ref=float(m["local_energy_reference_energy_J"])
 time=np.asarray(arrays["trajectory__time_s"],float);energy=np.asarray(arrays["trajectory__energy_total_J"],float)
 work=np.asarray(arrays["trajectory__external_energy_integral_J"],float)
 if time.ndim!=1 or len(time)<1 or energy.shape!=time.shape or work.shape!=time.shape or not all(np.all(np.isfinite(a)) for a in (time,energy,work)):
  raise ValueError("Finite coherent energy/work/time history required")
 if np.any(np.diff(time)<=0):raise ValueError("Strictly increasing accepted physical times required")
 if ref<=0 or not math.isfinite(ref):raise ValueError("Positive old guard reference required")
 h.close_scalar(float(time[-1])*1e12,float(m["time_ps"]),"guard endpoint clock")
 dw=np.diff(work);defect=np.diff(energy)-dw;denominator=ref+np.cumsum(abs(dw));fractions=abs(defect)/denominator
 anchor=float(m["guard_active_since_time_ps"])*1e-12;eps=np.finfo(float).eps
 time_roundoff=32*eps*max(float(abs(time).max()),1e-12)
 indices=np.flatnonzero(abs(time-anchor)<=time_roundoff);post=time[:-1]>=anchor-time_roundoff
 roundoff_J=16*eps*(abs(energy[:-1])+abs(energy[1:])+abs(work[:-1])+abs(work[1:]))
 fraction_roundoff=roundoff_J/denominator
 maximum=float(fractions.max()) if len(fractions) else 0.
 guarded_max=float(fractions[post].max()) if np.any(post) else 0.
 maximum_roundoff=float(fraction_roundoff.max()) if len(fractions) else 0.
 guarded_roundoff=float(fraction_roundoff[post].max()) if np.any(post) else 0.
 abs_work=math.fsum(abs(float(v)) for v in dw);sum_defect=math.fsum(abs(float(v)) for v in defect)
 work_sum_roundoff=16*eps*math.fsum(float(v) for v in abs(work[:-1])+abs(work[1:]))
 ledger_roundoff=max(work_sum_roundoff,32*eps*max(ref,abs_work))
 saved=parent.get("local_energy_control",{})
 count=float(m["local_energy_rejections"])
 checks=dict(strict_parameter_and_manifest_tolerance=parent.get("numerical",{}).get(IDENTITY_KEY)==tol,
  reference_matches_source=m["source"].get("retained_energy_J",0)<=0 or math.isclose(ref,m["source"]["retained_energy_J"],rel_tol=1e-13,abs_tol=0.),
  anchor_is_saved_accepted_time=bool(len(indices)==1 and anchor>=-time_roundoff and anchor<=time[-1]+time_roundoff),
  every_post_anchor_step_bounded=bool(np.all(fractions[post]<=tol+fraction_roundoff[post])),
  absolute_work_ledger=abs(float(m["absolute_external_energy_integral_J"])-abs_work)<=ledger_roundoff,
  absolute_defect_ledger=abs(float(m["sum_abs_local_energy_defect_J"])-sum_defect)<=ledger_roundoff,
  maximum_all_ledger=abs(float(m["max_accepted_local_energy_fraction"])-maximum)<=maximum_roundoff,
  maximum_guarded_ledger=abs(float(m["max_guarded_accepted_local_energy_fraction"])-guarded_max)<=guarded_roundoff,
  rejection_count_valid=count>=0 and not isinstance(m["local_energy_rejections"],bool) and count==int(count),
  checkpoint_manifest_guard_ledger_matches=all(k in saved and saved[k]==m[k] for k in h.GUARD_KEYS))
 checks={k:bool(v) for k,v in checks.items()}
 if not all(checks.values()):raise ValueError("Old guard history audit failed: "+",".join(k for k,v in checks.items() if not v))
 return dict(passed=True,checks=checks,tolerance=tol,guard_active_since_time_ps=m["guard_active_since_time_ps"],
  post_anchor_accepted_steps=int(post.sum()),legacy_accepted_steps=int((~post).sum()),
  absolute_work_reconstructed_J=abs_work,sum_abs_local_defect_reconstructed_J=sum_defect,
  maximum_all_reconstructed_fraction=maximum,maximum_guarded_reconstructed_fraction=guarded_max,
  ledger_sum_roundoff_bound_J=ledger_roundoff,maximum_fraction_roundoff=maximum_roundoff,
  maximum_guarded_fraction_roundoff=guarded_roundoff,
  roundoff_definition="Per step16eps*(abs(Eold)+abs(Enew)+abs(Wold)+abs(Wnew)); ledger16eps*sum(abs(Wold)+abs(Wnew)), lower arithmetic bound32eps*max(Eref,sumabsdW)",
  scope="The historical guard is reconstructed before its window closes; no retrospective certification under the new tolerance.")

def eligible(parent,m,arrays):
 if parent.get("run_complete") is not False or parent.get("status") not in ("RUNNING","INTERRUPTED","EXPLICIT_TRANSFER_PREPARED","EXPLICIT_GUARD_TIGHTENING_PREPARED"):
  raise ValueError("Only an unfinished pre-trigger checkpoint may tighten its guard")
 c=m.get("controller")
 if c is None or c["state"].get("status")!="ACTIVE":
  raise ValueError("Only an ACTIVE pre-trigger observation controller may continue")
 s=c["state"]
 if any(s.get(k) is not None for k in w.EVENT_KEYS) or s.get("failure_reason") is not None or int(s.get("rejected_candidate_count",0))!=0:
  raise ValueError("Pending, previous or terminal threshold events are excluded")
 w.no_history_crossing(arrays,m)
 old=m["parameter_identity"].get(IDENTITY_KEY)
 if isinstance(old,bool) or not isinstance(old,(int,float)) or not math.isfinite(old) or old<=0:
  raise ValueError("A positive finite existing explicit guard is required")
 w.validate_guard(m,parent,arrays)
 validate_old_guard_history(parent,m,arrays)
 if m["guard_active_since_time_ps"]>m["time_ps"]:raise ValueError("Old guard activates after the checkpoint")
 if len(m["source_pins_sha256"])!=13:raise ValueError("Complete13-pin namespace required")
 return w.resource_start(parent,m)
def guard_change(parent,m,arrays,tolerance):
 eligible(parent,m,arrays)
 old=m["parameter_identity"][IDENTITY_KEY]
 if isinstance(tolerance,bool) or not isinstance(tolerance,(int,float)) or not math.isfinite(tolerance) or not 0<tolerance<old:
  raise ValueError("New local-energy tolerance must be strictly positive and strictly smaller")
 ledger,change=h.energy_guard_transition(m,arrays,parent,float(tolerance))
 for k in h.GUARD_KEYS:
  if k not in WINDOW_KEYS and ledger[k]!=m[k]:
   raise ValueError("Global guard ledger changed: "+k)
 if ledger["guard_active_since_time_ps"]!=m["time_ps"] or ledger["max_guarded_accepted_local_energy_fraction"]!=0.:
  raise ValueError("New guard window must start exactly at the physical anchor")
 change.update(prospective_only=True,prior_steps_certified_under_new_tolerance=False,
  previous_guard_active_since_time_ps=m["guard_active_since_time_ps"],
  previous_guard_max_accepted_fraction=m["max_guarded_accepted_local_energy_fraction"],
  previous_guard_ledger={k:m[k] for k in h.GUARD_KEYS},
  previous_guard_history_validation=validate_old_guard_history(parent,m,arrays),
  global_error_origin_reset=False,global_denominator_reset=False,heat_or_state_correction_applied=False)
 return ledger,change
def _fresh_proof(run,arrays,parent,m,logpath,*,expected_output_dir=None):
 proof,payload=w.fresh_evidence(run,arrays,parent,m,logpath,expected_output_dir=expected_output_dir)
 proof["numerical_segments"][0]["duration_ps"]=m["parameter_identity"].get("duration_ps")
 return proof,payload
def capture(run,anchor,*,fresh_start_log=None):
 helper_integrity();run=h.scratch(run);anchor=h.scratch(anchor)
 payload=(run/"checkpoint_latest.npz").read_bytes();mb=(run/"manifest.json").read_bytes()
 arrays=h.decode(payload);parent=json.loads(mb);m=h.audit_arrays(arrays,parent)
 eligible(parent,m,arrays);archive=h.audit_runtime(run/"runtime_sources",m["source_pins_sha256"])
 proof=w.inherited_evidence(run,m,parent);logpayload=None
 if proof is None:proof,logpayload=_fresh_proof(run,arrays,parent,m,fresh_start_log)
 elif fresh_start_log is not None:raise ValueError("Inherited run may not be relabeled fresh")
 anchor=h.newdir(anchor);prov=anchor/"provenance";prov.mkdir()
 (anchor/"checkpoint_original.npz").write_bytes(payload);(anchor/"manifest_original.json").write_bytes(mb)
 shutil.copytree(run/"runtime_sources",anchor/"runtime_sources")
 for rel in ("dc","meshes"):
  if (run/rel).exists():shutil.copytree(run/rel,anchor/rel)
 for rel in ("material.npz","material.json","source_initial_state.npz"):shutil.copy2(run/rel,anchor/rel)
 dc=json.loads((anchor/"dc/effective_dc_receipt.json").read_text())
 if h.filehash(anchor/"dc/effective_dc.npz")!=dc["fields_sha256"]:raise ValueError("DC hash differs")
 if logpayload is not None:
  (prov/"fresh_original_console.log").write_bytes(logpayload);h.writejson(prov/"fresh_STARTED.json",proof["fresh_STARTED"])
 else:shutil.copy2(proof["parent_transfer_receipt_path"],prov/"parent_transfer_receipt.json")
 h.writejson(prov/"guard_parent_lineage.json",proof)
 for rel in ("parent_preparation_receipt.json",):
  if (run/rel).is_file():shutil.copy2(run/rel,anchor/rel)
 if (run/"diagnostic_provenance").is_dir():shutil.copytree(run/"diagnostic_provenance",anchor/"diagnostic_provenance")
 hashes={str(p.relative_to(anchor)):h.filehash(p) for p in anchor.rglob("*") if p.is_file() and "runtime_sources" not in p.parts}
 cr=dict(schema=SCHEMA,action="GUARDED_PRETRIGGER_CAPTURE_NO_SIGNAL_NO_EVOLUTION",parent_output_dir=str(run),
  checkpoint_bytes_sha256=h.digest(payload),manifest_sha256=h.digest(mb),runtime_archive=archive,
  physical_file_sha256=hashes,state_arrays_sha256={k:h.arrayhash(arrays[k]) for k in h.STATE},
  conserved_metadata={k:m[k] for k in h.CONSERVED},parameter_identity=m["parameter_identity"],
  source_pins_sha256=m["source_pins_sha256"],immutable_anchor=str(anchor),lineage_kind=proof["kind"],
  helper_sha256=h.filehash(Path(__file__)),live_source_checkpoint_may_advance=True)
 h.writejson(anchor/"capture_receipt.json",cr);return cr
def segments_for_change(proof,m,identity):
 segments=copy.deepcopy(proof["numerical_segments"]);h.lineage.segments_valid(segments)
 for k in h.lineage.KEYS:
  if segments[-1].get(k)!=m["parameter_identity"].get(k):raise ValueError("Ancestral numerical policy differs: "+k)
 now=float(m["time_ps"])
 if segments[-1]["start_time_ps"]>now:raise ValueError("Ancestral policy begins after anchor")
 nxt=copy.deepcopy(segments[-1]);segments[-1]["end_time_ps"]=now
 nxt.update(start_time_ps=now,end_time_ps=None,origin="prospective_local_energy_guard_tightening",
  local_energy_tolerance=identity[IDENTITY_KEY],energy_guard_enabled=True,duration_ps=identity.get("duration_ps"),
  guard_active_since_time_ps=now)
 segments.append(nxt);h.lineage.segments_valid(segments);return segments
def global_error_at_anchor(m,arrays):
 defect=float(arrays["trajectory__energy_total_J"][-1])-float(m["initial_total_energy_J"])-float(m["external_energy_integral_J"])
 denominator=float(m["local_energy_reference_energy_J"])+float(m["absolute_external_energy_integral_J"])
 if denominator<=0 or not math.isfinite(defect+denominator):raise ValueError("Invalid global error reference")
 return dict(defect_J=defect,denominator_J=denominator,absolute_relative_defect=abs(defect)/denominator,
  formula="abs(U(t)-U(0)-Wext(t))/(Eref + sum(abs(accepted delta Wext)))",global_origin_time_ps=0.)
def tighten(anchor,output,tolerance):
 helper_integrity();anchor=h.scratch(anchor);output=h.scratch(output)
 cr=json.loads((anchor/"capture_receipt.json").read_text());payload=(anchor/"checkpoint_original.npz").read_bytes()
 if h.digest(payload)!=cr["checkpoint_bytes_sha256"]:raise ValueError("Anchored checkpoint changed")
 for rel,want in cr["physical_file_sha256"].items():
  if h.filehash(anchor/rel)!=want:raise ValueError("Anchored physical evidence changed: "+rel)
 arrays=h.decode(payload);parent=json.loads((anchor/"manifest_original.json").read_text());old=h.audit_arrays(arrays,parent)
 cpustart=eligible(parent,old,arrays);ledger,change=guard_change(parent,old,arrays,tolerance)
 runtime_audit=h.audit_runtime(anchor/"runtime_sources",old["source_pins_sha256"])
 proof=json.loads((anchor/"provenance/guard_parent_lineage.json").read_text())
 if not proof["numerical_lineage_complete"]:raise ValueError("Incomplete source/numerical ancestry")
 if proof["kind"]=="fresh":
  _fresh_proof(anchor,arrays,parent,old,anchor/"provenance/fresh_original_console.log",expected_output_dir=cr["parent_output_dir"])
 elif proof["kind"]=="inherited":
  prior=anchor/"provenance/parent_transfer_receipt.json"
  if h.filehash(prior)!=proof["parent_transfer_receipt_sha256"]:raise ValueError("Captured ancestor receipt differs")
  h.lineage.verify_receipt(prior,old["source_pins_sha256"])
 else:raise ValueError("Unknown source ancestry kind")
 identity=copy.deepcopy(old["parameter_identity"]);identity[IDENTITY_KEY]=float(tolerance)
 if {k:v for k,v in identity.items() if k!=IDENTITY_KEY}!={k:v for k,v in old["parameter_identity"].items() if k!=IDENTITY_KEY}:
  raise ValueError("Only local-energy tolerance may change")
 segments=segments_for_change(proof,old,identity);global_before=global_error_at_anchor(old,arrays)
 output=h.newdir(output);prov=output/"provenance";prov.mkdir()
 for src,dst in [("checkpoint_original.npz","parent_checkpoint.npz"),("manifest_original.json","parent_manifest.json"),("capture_receipt.json","capture_receipt.json")]:
  shutil.copy2(anchor/src,prov/dst)
 for p in (anchor/"provenance").iterdir():
  if p.is_file():shutil.copy2(p,prov/p.name)
 shutil.copytree(anchor/"runtime_sources",output/"runtime_sources")
 for rel in ("dc","meshes"):
  if (anchor/rel).exists():shutil.copytree(anchor/rel,output/rel)
 for rel in ("material.npz","material.json","source_initial_state.npz"):shutil.copy2(anchor/rel,output/rel)
 new=copy.deepcopy(old);new.update(parameter_identity=identity,**ledger)
 arrays["metadata_json"]=np.array(json.dumps(new,allow_nan=False));np.savez_compressed(output/"checkpoint_latest.npz",**arrays)
 manifest=copy.deepcopy(parent);manifest.update(parameter_identity=identity,status="EXPLICIT_GUARD_TIGHTENING_PREPARED",
  run_complete=False,explicit_transfer_receipt="provenance/transfer_receipt.json",explicit_energy_guard_change=change,
  local_energy_control=copy.deepcopy(ledger))
 for k in ("execution","output_sha256","error"):manifest.pop(k,None)
 manifest["numerical"].update(local_energy_tolerance=float(tolerance))
 manifest["numerical"]["local_energy_guard"].update(enabled=True,reference_energy_J=ledger["local_energy_reference_energy_J"],
  acceptance=change["denominator"],pre_guard_history="Historical policy retained; new tolerance applies prospectively from explicit anchor only")
 h.writejson(output/"manifest.json",manifest)
 g,energy,argv=w.reconstruct(output/"runtime_sources",output,identity,cpustart,parent,arrays)
 np.savez_compressed(output/"trajectory.npz",**{k[12:]:v for k,v in arrays.items() if k.startswith("trajectory__")})
 snap={k[11:]:v for k,v in arrays.items() if k.startswith("snapshots__")}
 snap.update(coordinates_m=g.coordinates_m,triangles=g.triangles,edges=g.graph.edges);np.savez_compressed(output/"snapshots.npz",**snap)
 check=h.decode((output/"checkpoint_latest.npz").read_bytes());after=h.audit_arrays(check,manifest);before=h.decode(payload)
 for k,v in before.items():
  if k!="metadata_json" and h.arrayhash(v)!=h.arrayhash(check[k]):raise ValueError("Physical/history/snapshot array changed")
 for k,v in old.items():
  if k not in WINDOW_KEYS|{"parameter_identity"} and after[k]!=v:raise ValueError("Metadata changed beyond prospective guard window: "+k)
 if h.audit_runtime(output/"runtime_sources",old["source_pins_sha256"])!=runtime_audit:raise ValueError("Source archive changed")
 global_after=global_error_at_anchor(after,check)
 if global_after!=global_before:raise ValueError("Global energy origin, defect or denominator changed")
 receipt=dict(schema=SCHEMA,action="EXPLICIT_LOCAL_ENERGY_GUARD_TIGHTENING_PREPARED_NO_EVOLUTION_NO_SIGNAL",
  parent=dict(output_dir=cr["parent_output_dir"],anchor=str(anchor),checkpoint_bytes_sha256=h.digest(payload),
   manifest_sha256=cr["manifest_sha256"],runtime_archive_manifest_sha256=runtime_audit["manifest_sha256"]),
  output_dir=str(output),state_arrays_sha256_before={k:h.arrayhash(before[k]) for k in h.STATE},
  state_arrays_sha256_after={k:h.arrayhash(check[k]) for k in h.STATE},
  conserved_metadata_before={k:old[k] for k in h.CONSERVED},conserved_metadata_after={k:after[k] for k in h.CONSERVED},
  pins_before=old["source_pins_sha256"],pins_after=after["source_pins_sha256"],
  frontend_changes={},implementation_changes={},RHS_sources_identical=True,physical_RHS_bitwise_equivalence_certified=False,
  numerical_segments=segments,horizon_change=None,energy_guard_change=change,
  guard_ledger_before={k:old[k] for k in h.GUARD_KEYS},guard_ledger_after={k:after[k] for k in h.GUARD_KEYS},
  global_energy_error_before=global_before,global_energy_error_after=global_after,
  next_step_ps_before=old["next_step_ps"],next_step_ps_after=after["next_step_ps"],
  missing_old_adaptive_counters=[],adaptive_counter_policy="Every recorded counter and global ledger preserved exactly",
  preserved_history_sha256={k[12:]:h.arrayhash(v) for k,v in before.items() if k.startswith("trajectory__")},
  preserved_snapshots_sha256={k[11:]:h.arrayhash(v) for k,v in before.items() if k.startswith("snapshots__")},
  preserved_all_metadata_except_declared_guard_window=True,resources_before=parent["resources"],resources_after=manifest["resources"],
  lineage_kind=proof["kind"],parent_guard_lineage_sha256=h.filehash(anchor/"provenance/guard_parent_lineage.json"),
  fresh_initialization_proof=proof if proof["kind"]=="fresh" else None,parent_transfer_receipt_sha256=proof.get("parent_transfer_receipt_sha256"),
  physical_identity_verified=dict(parameters_equal=True,mesh_sha256=g.receipt["mesh_sha256"],
   DC_sha256=h.filehash(output/"dc/effective_dc.npz"),material_sha256=h.filehash(output/"material.npz"),
   reconstructed_total_energy_J=energy,photon_source_reapplied=False,no_collision_RHS_evaluated=True),
  new_checkpoint_bytes_sha256=h.filehash(output/"checkpoint_latest.npz"),
  new_runtime_archive_manifest_sha256=runtime_audit["manifest_sha256"],resume_argv=argv,
  parent_helper_sha256=w.EXPECTED_PARENT_HELPER_SHA256,window_helper_sha256=EXPECTED_WINDOW_HELPER_SHA256,
  helper_sha256=h.filehash(Path(__file__)),admission=False)
 h.writejson(prov/"transfer_receipt.json",receipt)
 cmd="env PYTHONDONTWRITEBYTECODE=1 PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/jdiaz/.conda/envs/snspd/bin/python "+shlex.quote(str(output/"runtime_sources/pipelines/05_effective_photon_run.py"))+" "+shlex.join(argv)
 (prov/"resume_command.txt").write_text(cmd+"\n");return receipt
def main():
 p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest="command",required=True)
 c=sub.add_parser("capture");c.add_argument("--run",required=True);c.add_argument("--anchor",required=True);c.add_argument("--fresh-start-log")
 t=sub.add_parser("tighten");t.add_argument("--anchor",required=True);t.add_argument("--output",required=True);t.add_argument("--local-energy-tolerance",required=True,type=float)
 a=p.parse_args()
 r=capture(a.run,a.anchor,fresh_start_log=a.fresh_start_log) if a.command=="capture" else tighten(a.anchor,a.output,a.local_energy_tolerance)
 print(json.dumps(dict(action=r["action"],admission=False),indent=2))
if __name__=="__main__":main()
