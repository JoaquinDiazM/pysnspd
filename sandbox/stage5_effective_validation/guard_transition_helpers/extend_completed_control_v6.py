"""Strict development-control window extension; never signals or evolves a job.

capture --run COMPLETED_CONTROL --anchor NEW_DIR [--fresh-start-log ORIGINAL_LOG]
extend --anchor IMMUTABLE_ANCHOR --output NEW_DIR --duration-ps 40|none

Only duration_ps may change, from exactly20ps to a larger value no greater than
the unchanged no-trigger ceiling (at most200ps), or None to use that ceiling.
All13 source pins, solver/guard/resources, physical arrays and metadata are fixed.
No terminal or pending detection event can be reopened. A fresh run requires
explicit STARTED-at-zero evidence; inherited controls require verified ancestry.
"""
from __future__ import annotations
import argparse,copy,importlib.util,json,math,shlex,shutil,sys
from pathlib import Path
import numpy as np
BASE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location("_window_v6_verified",BASE/"checkpoint_transfer_v5_v6.py")
h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h)
EXPECTED_PARENT_HELPER_SHA256="8f65a3ec53a2145d26ae39b206cac617839c3ec759da176c8f3a82b45874654c"
SCHEMA="pysnspd.explicit_checkpoint_transfer.v3"
EVENT_KEYS=("candidate_crossing_s","crossing_time_s","confirmation_time_s","stop_target_s","actual_stop_s")
RESOURCE_KEYS=("actual_processes","actual_threads","requested_threads","reserved_physical_core_count","affinity_cpus")
def helper_integrity():
 if h.filehash(BASE/"checkpoint_transfer_v5_v6.py")!=EXPECTED_PARENT_HELPER_SHA256:
  raise ValueError("Verified v6 parent helper changed")
def horizon(parent,m,requested):
 c=m.get("controller");ident=m["parameter_identity"]
 if parent.get("run_complete") is not True or parent.get("status")!="DEVELOPMENT_DURATION_REACHED":
  raise ValueError("Only completed DEVELOPMENT_DURATION_REACHED controls are eligible")
 if c is None or c["state"].get("status")!="ACTIVE":
  raise ValueError("Only ACTIVE observation controllers may continue")
 s=c["state"]
 if any(s.get(k) is not None for k in EVENT_KEYS) or s.get("failure_reason") is not None:
  raise ValueError("Pending/confirmed/terminal detection events may never be reopened")
 if int(s.get("rejected_candidate_count",0))!=0:
  raise ValueError("Control has a previous threshold crossing")
 prior=ident.get("duration_ps");ceiling=ident.get("ceiling_ps")
 if isinstance(prior,bool) or prior!=20.:
  raise ValueError("This workflow extends only a completed20ps control")
 if isinstance(ceiling,bool) or not isinstance(ceiling,(int,float)) or not math.isfinite(ceiling) or not 20.<ceiling<=200.:
  raise ValueError("Unchanged no-trigger ceiling must be finite, above20ps and at most200ps")
 h.close_scalar(float(m["time_ps"]),20.,"completed control duration")
 h.close_scalar(float(c["config"]["no_trigger_time_ceiling_s"]),float(ceiling)*1e-12,"controller ceiling")
 if requested is not None:
  if isinstance(requested,bool) or not isinstance(requested,(int,float)) or not math.isfinite(requested) or not 20.<requested<=ceiling:
   raise ValueError("New duration must exceed20ps and not exceed the unchanged ceiling")
  requested=float(requested)
 return dict(runtime_key="duration_ps",before_ps=20.,after_ps=requested,effective_time_ps=20.,
  parent_status=parent["status"],parent_run_complete=True,controller_status_preserved="ACTIVE",
  physical_clock_preserved=True,no_event_reopened=True,unchanged_no_trigger_ceiling_ps=ceiling,
  None_duration_semantics="No development cutoff; original trigger plus margin or finite censoring ceiling remains")
def validate_guard(m,parent,arrays):
 keys=(*h.GUARD_KEYS,"stable_accepted_steps","consecutive_rejections")
 if any(k not in m for k in keys):raise ValueError("Incomplete guarded/adaptive checkpoint")
 tol=m["parameter_identity"].get("local_energy_tolerance")
 if tol is None:raise ValueError("This workflow requires the existing explicit local-energy guard")
 ledger,_=h.energy_guard_transition(m,arrays,parent,tol)
 if ledger!={k:m[k] for k in h.GUARD_KEYS}:raise ValueError("Guard ledger changed during validation")
 return ledger
def recorded_resources(parent):
 r=parent.get("resources",{});result={k:r.get(k) for k in RESOURCE_KEYS}
 actual=r.get("actual_affinity_cpus");alias=r.get("affinity_cpus")
 if actual is not None and alias is not None and actual!=alias:
  raise ValueError("Ambiguous actual affinity metadata")
 result["affinity_cpus"]=actual if actual is not None else alias
 return result
def resource_start(parent,m):
 r=recorded_resources(parent);cpus=r.get("affinity_cpus")
 threads=m["parameter_identity"]["threads"]
 if not isinstance(cpus,list) or not cpus or len(cpus)!=threads or sorted(cpus)!=list(range(min(cpus),min(cpus)+threads)):
  raise ValueError("Parent resource affinity is incomplete or not a contiguous unchanged reservation")
 if r.get("requested_threads")!=threads or r.get("reserved_physical_core_count")!=threads:
  raise ValueError("Parent reserved resource metadata differs from identity")
 if r.get("actual_threads")!=threads or r.get("actual_processes")!=1:
  raise ValueError("Actual parent worker count must match the preserved implementation")
 return min(cpus)
def no_history_crossing(arrays,m):
 cfg=m["controller"]["config"];sgn=1 if cfg["polarity"]=="positive" else -1
 excursion=sgn*(arrays["trajectory__Vout_V"]-cfg["baseline_V"])
 if np.any(excursion>=cfg["threshold_V"]):raise ValueError("History contains a threshold crossing")
def _source_arrays(run,arrays,parent,m):
 initial=h.decode((run/"source_initial_state.npz").read_bytes())
 required=("delta_bar","p","n","temperature_K","node_deposition_J","coordinates_m")
 if not all(k in initial and np.all(np.isfinite(initial[k])) for k in required):
  raise ValueError("Fresh source state is incomplete/nonfinite")
 for k in ("delta_bar","p","n"):
  if initial[k].shape!=arrays[k].shape:raise ValueError("Fresh source quadrature differs")
 if np.any((initial["p"]<0)|(initial["p"]>1)) or np.any(initial["n"]<0):
  raise ValueError("Fresh source occupations violate support")
 Eu=parent["physicalparams"]["energy_scale_J"]
 if not np.array_equal(arrays["snapshots__delta_J"][0],initial["delta_bar"]*Eu):
  raise ValueError("Initial saved gap differs from source state")
 dc=h.decode((run/"dc/effective_dc.npz").read_bytes())
 if not np.array_equal(initial["delta_bar"],dc["delta_bar"]):
  raise ValueError("Fresh source did not preserve its fresh DC gap")
 if not np.array_equal(initial["coordinates_m"],dc["coordinates_m"]):
  raise ValueError("Fresh source and DC coordinates differ")
 h.close_scalar(float(np.sum(initial["node_deposition_J"])),m["source"]["retained_energy_J"],"fresh photon deposited energy")
 return initial
def fresh_evidence(run,arrays,parent,m,logpath,*,expected_output_dir=None):
 if logpath is None:raise ValueError("Fresh controls require explicit original STARTED log evidence")
 if m["guard_active_since_time_ps"]!=0 or m["energy_guard_ledger_origin"]!="accepted_midpoint_receipts_since_time_zero":
  raise ValueError("Fresh control guard must be anchored at time zero without inherited policy")
 if m["source"].get("old_prepared_photon_state_reused") is not False:
  raise ValueError("Fresh control source cannot be inherited")
 if parent.get("source")!=m["source"]:raise ValueError("Fresh manifest source differs")
 p=h.scratch(logpath);payload=p.read_bytes()
 lines=payload.decode("utf-8").splitlines()
 events=[]
 for line in lines:
  try:j=json.loads(line)
  except json.JSONDecodeError:continue
  if j.get("status")=="STARTED":events.append(j)
 if len(events)!=1:raise ValueError("Fresh proof requires exactly one STARTED record")
 e=events[0]
 if e.get("time_ps")!=0. or Path(e.get("output_dir","")).resolve()!=(run if expected_output_dir is None else Path(expected_output_dir).resolve()) or e.get("source")!=m["source"]:
  raise ValueError("STARTED must bind the same fresh source/output at physical time zero")
 for k in ("dt_ps",):
  if e.get(k)!=m["parameter_identity"].get(k):raise ValueError("STARTED numerical identity differs")
 for k in RESOURCE_KEYS:
  if e.get(k)!=recorded_resources(parent).get(k):raise ValueError("STARTED resources differ: "+k)
 if e.get("nodes")!=len(arrays["delta_bar"]) or e.get("electron_states")!=arrays["p"].shape[1] or e.get("phonon_modes")!=arrays["n"].shape[1]:
  raise ValueError("STARTED quadratures differ")
 _source_arrays(run,arrays,parent,m)
 seg=dict(start_time_ps=0.,end_time_ps=None,origin="verified_fresh_STARTED_zero",
  duration_ps=20.,energy_guard_enabled=True,**{k:m["parameter_identity"][k] for k in h.lineage.KEYS})
 return dict(kind="fresh",numerical_lineage_complete=True,numerical_segments=[seg],
  fresh_STARTED=e,original_log_path=str(p),original_log_sha256=h.digest(payload),
  source_initial_state_sha256=h.filehash(run/"source_initial_state.npz"),
  no_prior_transfer_invented=True,guard_anchored_at_zero=True),payload
def inherited_evidence(run,m,parent):
 prior=h.lineage.prior_receipt(run,parent)
 if prior is None:return None
 receipt=json.loads(prior.read_text())
 if receipt["pins_after"]!=m["source_pins_sha256"]:raise ValueError("Ancestral pins differ")
 proof=h.lineage.supplement(run)
 if not proof["numerical_lineage_complete"]:raise ValueError("Incomplete recovered ancestry: "+proof["reason"])
 for k in h.lineage.KEYS:
  if proof["numerical_segments"][-1].get(k)!=m["parameter_identity"].get(k):
   raise ValueError("Recovered current numerical policy differs: "+k)
 return dict(kind="inherited",numerical_lineage_complete=True,numerical_segments=proof["numerical_segments"],
  parent_transfer_receipt_path=str(prior),parent_transfer_receipt_sha256=h.filehash(prior),
  recovered_ancestry=proof)
def capture(run,anchor,*,fresh_start_log=None):
 helper_integrity();run=h.scratch(run);anchor=h.scratch(anchor)
 payload=(run/"checkpoint_latest.npz").read_bytes();mb=(run/"manifest.json").read_bytes()
 arrays=h.decode(payload);parent=json.loads(mb);m=h.audit_arrays(arrays,parent)
 horizon(parent,m,40. if m["parameter_identity"].get("ceiling_ps",0)>=40 else None)
 if len(m["source_pins_sha256"])!=13:raise ValueError("Complete13-pin namespace required")
 validate_guard(m,parent,arrays);resource_start(parent,m);no_history_crossing(arrays,m)
 archive=h.audit_runtime(run/"runtime_sources",m["source_pins_sha256"])
 proof=inherited_evidence(run,m,parent);logpayload=None
 if proof is None:proof,logpayload=fresh_evidence(run,arrays,parent,m,fresh_start_log)
 elif fresh_start_log is not None:raise ValueError("Inherited control may not be relabeled fresh")
 anchor=h.newdir(anchor);prov=anchor/"provenance";prov.mkdir()
 (anchor/"checkpoint_original.npz").write_bytes(payload);(anchor/"manifest_original.json").write_bytes(mb)
 shutil.copytree(run/"runtime_sources",anchor/"runtime_sources")
 for rel in ("dc","meshes"):
  if (run/rel).exists():shutil.copytree(run/rel,anchor/rel)
 for rel in ("material.npz","material.json","source_initial_state.npz"):shutil.copy2(run/rel,anchor/rel)
 dc=json.loads((anchor/"dc/effective_dc_receipt.json").read_text())
 if h.filehash(anchor/"dc/effective_dc.npz")!=dc["fields_sha256"]:raise ValueError("DC hash differs")
 if logpayload is not None:
  (prov/"fresh_original_console.log").write_bytes(logpayload)
  h.writejson(prov/"fresh_STARTED.json",proof["fresh_STARTED"])
 else:
  shutil.copy2(proof["parent_transfer_receipt_path"],prov/"parent_transfer_receipt.json")
 h.writejson(prov/"window_parent_lineage.json",proof)
 for rel in ("parent_preparation_receipt.json",):
  if (run/rel).is_file():shutil.copy2(run/rel,anchor/rel)
 if (run/"diagnostic_provenance").is_dir():shutil.copytree(run/"diagnostic_provenance",anchor/"diagnostic_provenance")
 hashes={str(p.relative_to(anchor)):h.filehash(p) for p in anchor.rglob("*") if p.is_file() and "runtime_sources" not in p.parts}
 cr=dict(schema=SCHEMA,action="COMPLETED_CONTROL_CAPTURE_NO_SIGNAL_NO_EVOLUTION",parent_output_dir=str(run),
  checkpoint_bytes_sha256=h.digest(payload),manifest_sha256=h.digest(mb),runtime_archive=archive,
  physical_file_sha256=hashes,state_arrays_sha256={k:h.arrayhash(arrays[k]) for k in h.STATE},
  conserved_metadata={k:m[k] for k in h.CONSERVED},parameter_identity=m["parameter_identity"],
  source_pins_sha256=m["source_pins_sha256"],immutable_anchor=str(anchor),lineage_kind=proof["kind"])
 h.writejson(anchor/"capture_receipt.json",cr);return cr
def reconstruction_energy(system,state):
 # Native functionals only. No collision RHS, no photon preparation, no time step.
 return float(system.functional.evaluate(state.delta_bar,state.p).energy_J+
  np.sum(system.kinetics.phonon_energy_weights_J*state.n)+system.circuit.energy(state.circuit_SI))
def reconstruct(runtime,output,identity,cpu_start,parent,arrays):
 runner=h.load_runner(runtime);argv=h.resume_argv(identity,output,cpu_start)
 args=runner.parser().parse_args(argv);args.threads=1
 system,dc,g,material,params=runner.build_system(args,output,resume=True)
 if params!=parent["physicalparams"] or material.manifest!=parent["material"]:
  raise ValueError("Reconstructed physical/material identity differs")
 if g.receipt["mesh_sha256"]!=parent["geometry"]["mesh_sha256"]:raise ValueError("Reconstructed mesh differs")
 runner.TriggerWindow.from_checkpoint(h.metadata(arrays)["controller"])
 state=runner.EffectiveState(*(arrays[k].copy() for k in h.STATE))
 energy=reconstruction_energy(system,state)
 if not math.isclose(float(system.kinetic_residual_energy_scale_J),h.metadata(arrays)["local_energy_reference_energy_J"],rel_tol=1e-13,abs_tol=0.):
  raise ValueError("Reconstructed energy-guard reference differs")
 initial=h.decode((output/"source_initial_state.npz").read_bytes())
 initial_cm=np.array([arrays["trajectory__Ib_A"][0],arrays["trajectory__Is_A"][0],arrays["trajectory__vc_V"][0]])
 initial_state=runner.EffectiveState(initial["delta_bar"],initial["p"],initial["n"],initial_cm)
 initial_energy=reconstruction_energy(system,initial_state)
 h.close_scalar(initial_energy,h.metadata(arrays)["initial_total_energy_J"],"original source energy from native quadratures")
 h.close_scalar(initial_energy,float(arrays["trajectory__energy_total_J"][0]),"initial history energy")
 h.close_scalar(energy,float(arrays["trajectory__energy_total_J"][-1]),"reconstructed total energy")
 return g,energy,argv
def extend(anchor,output,duration_ps):
 helper_integrity();anchor=h.scratch(anchor);output=h.scratch(output)
 cr=json.loads((anchor/"capture_receipt.json").read_text())
 payload=(anchor/"checkpoint_original.npz").read_bytes()
 if h.digest(payload)!=cr["checkpoint_bytes_sha256"]:raise ValueError("Anchored checkpoint changed")
 for rel,want in cr["physical_file_sha256"].items():
  if h.filehash(anchor/rel)!=want:raise ValueError("Anchored physical evidence changed: "+rel)
 arrays=h.decode(payload);parent=json.loads((anchor/"manifest_original.json").read_text())
 old=h.audit_arrays(arrays,parent);change=horizon(parent,old,duration_ps)
 if len(old["source_pins_sha256"])!=13:raise ValueError("Complete13-pin namespace required")
 validate_guard(old,parent,arrays);cpustart=resource_start(parent,old);no_history_crossing(arrays,old)
 runtime_audit=h.audit_runtime(anchor/"runtime_sources",old["source_pins_sha256"])
 proof=json.loads((anchor/"provenance/window_parent_lineage.json").read_text())
 if not proof["numerical_lineage_complete"]:raise ValueError("Incomplete source/numerical ancestry")
 h.lineage.segments_valid(proof["numerical_segments"])
 if proof["kind"]=="fresh":
  fresh_evidence(anchor,arrays,parent,old,anchor/"provenance/fresh_original_console.log",expected_output_dir=cr["parent_output_dir"])
 elif proof["kind"]=="inherited":
  prior=anchor/"provenance/parent_transfer_receipt.json"
  if h.filehash(prior)!=proof["parent_transfer_receipt_sha256"]:raise ValueError("Captured ancestor receipt differs")
  _,_,_,_=h.lineage.verify_receipt(prior,old["source_pins_sha256"])
 else:raise ValueError("Unknown source ancestry kind")
 identity=copy.deepcopy(old["parameter_identity"]);identity["duration_ps"]=change["after_ps"]
 for k,v in old["parameter_identity"].items():
  if k!="duration_ps" and identity[k]!=v:raise ValueError("Only duration may change")
 output=h.newdir(output);prov=output/"provenance";prov.mkdir()
 for src,dst in [("checkpoint_original.npz","parent_checkpoint.npz"),("manifest_original.json","parent_manifest.json"),("capture_receipt.json","capture_receipt.json")]:
  shutil.copy2(anchor/src,prov/dst)
 for p in (anchor/"provenance").iterdir():
  if p.is_file():shutil.copy2(p,prov/p.name)
 shutil.copytree(anchor/"runtime_sources",output/"runtime_sources")
 for rel in ("dc","meshes"):
  if (anchor/rel).exists():shutil.copytree(anchor/rel,output/rel)
 for rel in ("material.npz","material.json","source_initial_state.npz"):shutil.copy2(anchor/rel,output/rel)
 new=copy.deepcopy(old);new["parameter_identity"]=identity
 arrays["metadata_json"]=np.array(json.dumps(new,allow_nan=False))
 np.savez_compressed(output/"checkpoint_latest.npz",**arrays)
 manifest=copy.deepcopy(parent);manifest.update(parameter_identity=identity,status="EXPLICIT_WINDOW_EXTENSION_PREPARED",
  run_complete=False,explicit_transfer_receipt="provenance/transfer_receipt.json",explicit_horizon_change=change)
 for k in ("execution","output_sha256","error"):manifest.pop(k,None)
 manifest.setdefault("observation",{}).update(status="EXPLICIT_WINDOW_EXTENSION_PREPARED",
  development_duration_s=None if duration_ps is None else float(duration_ps)*1e-12)
 h.writejson(output/"manifest.json",manifest)
 g,energy,argv=reconstruct(output/"runtime_sources",output,identity,cpustart,parent,arrays)
 np.savez_compressed(output/"trajectory.npz",**{k[12:]:v for k,v in arrays.items() if k.startswith("trajectory__")})
 snap={k[11:]:v for k,v in arrays.items() if k.startswith("snapshots__")}
 snap.update(coordinates_m=g.coordinates_m,triangles=g.triangles,edges=g.graph.edges)
 np.savez_compressed(output/"snapshots.npz",**snap)
 check=h.decode((output/"checkpoint_latest.npz").read_bytes());after=h.audit_arrays(check,manifest)
 before=h.decode(payload)
 for k,v in before.items():
  if k!="metadata_json" and h.arrayhash(v)!=h.arrayhash(check[k]):raise ValueError("Physical/history/snapshot array changed")
 for k,v in old.items():
  if k!="parameter_identity" and after[k]!=v:raise ValueError("Metadata changed beyond duration: "+k)
 if h.audit_runtime(output/"runtime_sources",old["source_pins_sha256"])!=runtime_audit:
  raise ValueError("Source archive changed")
 segments=copy.deepcopy(proof["numerical_segments"]);segments[-1]["end_time_ps"]=20.
 nxt=copy.deepcopy(segments[-1]);nxt.update(start_time_ps=20.,end_time_ps=None,origin="explicit_control_window_only",duration_ps=duration_ps)
 segments.append(nxt);h.lineage.segments_valid(segments)
 receipt=dict(schema=SCHEMA,action="EXPLICIT_CONTROL_WINDOW_EXTENSION_PREPARED_NO_EVOLUTION_NO_SIGNAL",
  parent=dict(output_dir=cr["parent_output_dir"],anchor=str(anchor),checkpoint_bytes_sha256=h.digest(payload),
   manifest_sha256=cr["manifest_sha256"],runtime_archive_manifest_sha256=runtime_audit["manifest_sha256"]),
  output_dir=str(output),state_arrays_sha256_before={k:h.arrayhash(before[k]) for k in h.STATE},
  state_arrays_sha256_after={k:h.arrayhash(check[k]) for k in h.STATE},
  conserved_metadata_before={k:old[k] for k in h.CONSERVED},
  conserved_metadata_after={k:after[k] for k in h.CONSERVED},
  pins_before=old["source_pins_sha256"],pins_after=after["source_pins_sha256"],
  frontend_changes={},implementation_changes={},RHS_sources_identical=True,physical_RHS_bitwise_equivalence_certified=False,
  numerical_segments=segments,horizon_change=change,energy_guard_change=None,
  guard_ledger_before={k:old[k] for k in h.GUARD_KEYS},guard_ledger_after={k:after[k] for k in h.GUARD_KEYS},
  guard_ledger_sha256_before=h.digest(json.dumps({k:old[k] for k in h.GUARD_KEYS},sort_keys=True).encode()),
  guard_ledger_sha256_after=h.digest(json.dumps({k:after[k] for k in h.GUARD_KEYS},sort_keys=True).encode()),
  next_step_ps_before=old["next_step_ps"],next_step_ps_after=after["next_step_ps"],
  missing_old_adaptive_counters=[],adaptive_counter_policy="Every recorded counter preserved exactly",
  preserved_history_sha256={k[12:]:h.arrayhash(v) for k,v in before.items() if k.startswith("trajectory__")},
  preserved_snapshots_sha256={k[11:]:h.arrayhash(v) for k,v in before.items() if k.startswith("snapshots__")},
  preserved_all_metadata_except_duration=True,resources_before=parent["resources"],resources_after=manifest["resources"],
  lineage_kind=proof["kind"],parent_window_lineage_sha256=h.filehash(anchor/"provenance/window_parent_lineage.json"),
  fresh_initialization_proof=proof if proof["kind"]=="fresh" else None,
  parent_transfer_receipt_sha256=proof.get("parent_transfer_receipt_sha256"),
  physical_identity_verified=dict(parameters_equal=True,mesh_sha256=g.receipt["mesh_sha256"],
   DC_sha256=h.filehash(output/"dc/effective_dc.npz"),material_sha256=h.filehash(output/"material.npz"),
   reconstructed_total_energy_J=energy,photon_source_reapplied=False,no_collision_RHS_evaluated=True),
  new_checkpoint_bytes_sha256=h.filehash(output/"checkpoint_latest.npz"),
  new_runtime_archive_manifest_sha256=runtime_audit["manifest_sha256"],resume_argv=argv,
  parent_helper_sha256=EXPECTED_PARENT_HELPER_SHA256,helper_sha256=h.filehash(Path(__file__)),admission=False)
 h.writejson(prov/"transfer_receipt.json",receipt)
 cmd="env PYTHONDONTWRITEBYTECODE=1 PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/jdiaz/.conda/envs/snspd/bin/python "+shlex.quote(str(output/"runtime_sources/pipelines/05_effective_photon_run.py"))+" "+shlex.join(argv)
 (prov/"resume_command.txt").write_text(cmd+"\n")
 return receipt
def duration(value):
 return None if value.lower()=="none" else float(value)
def main():
 p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest="command",required=True)
 c=sub.add_parser("capture");c.add_argument("--run",required=True);c.add_argument("--anchor",required=True);c.add_argument("--fresh-start-log")
 e=sub.add_parser("extend");e.add_argument("--anchor",required=True);e.add_argument("--output",required=True);e.add_argument("--duration-ps",required=True,type=duration)
 a=p.parse_args()
 r=capture(a.run,a.anchor,fresh_start_log=a.fresh_start_log) if a.command=="capture" else extend(a.anchor,a.output,a.duration_ps)
 print(json.dumps(dict(action=r["action"],admission=False),indent=2))
if __name__=="__main__":main()
