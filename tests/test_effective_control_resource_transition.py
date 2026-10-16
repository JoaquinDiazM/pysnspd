"""Semantic contracts for three closed resource tails; no physical evolution."""
import copy,importlib.util,json,sys,unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
BASE=Path(__file__).resolve().parent
def module(name):
 p=BASE/(name+".py")
 if not p.exists():p=BASE.parent/"sandbox/stage5_effective_validation"/(name+".py")
 spec=importlib.util.spec_from_file_location("sandbox.stage5_effective_validation._resource_candidate_"+name,p)
 m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m
r=module("control_resource_transition");d=module("guard_protocol_dispatcher");n8=module("num8_resource_equivalence")
from sandbox.stage5_effective_validation import fine_resource_transition as f
def fixture(name):
 c=r.case(name)
 i=dict(threads=c["old"],local_energy_tolerance=2.5e-6,step_tolerance=2e-7,dt_ps=.02,
  wavelength_nm=775.,spacing_nm=5.,duration_ps=None,length_nm=480. if c["kind"]=="domain" else 320.,
  xi_cutoff=32. if c["kind"]=="domain" else 64.,xi_order=6 if c["kind"]=="domain" else 8,
  phonon_nodes=96 if c["kind"]=="domain" else 144,transport_order=48 if c["kind"]=="domain" else 96)
 m=dict(parameter_identity=i,source_pins_sha256={str(k):str(k)*64 for k in range(13)},time_ps=30.,
  accepted_steps=4,total_rejections=2,next_step_ps=.0025,initial_total_energy_J=1e-19,external_energy_integral_J=2e-20,
  source=dict(retained_energy_J=1e-19),controller=dict(state=dict(status="ACTIVE",rejected_candidate_count=0,**{k:None for k in f.EVENTS})),
  guard_active_since_time_ps=22.3,local_energy_reference_energy_J=1e-19,energy_guard_ledger_origin="original",
  local_energy_rejections=1,limited_event_energy_J=1e-22,totalabs_event_energy_J=1e-20,
  absolute_external_energy_integral_J=2e-20,sum_abs_local_energy_defect_J=1e-25,
  max_accepted_local_energy_fraction=1e-6,max_guarded_accepted_local_energy_fraction=1e-6)
 arrays={k:np.array([[1.,2.]]) if k.startswith("snapshots__") else np.array([1.,2.]) for k in f.checkpoint_namespace() if k!="metadata_json"}
 arrays["metadata_json"]=np.array(json.dumps(m))
 p=dict(status="INTERRUPTED",run_complete=False,parameter_identity=copy.deepcopy(i),source_pins_sha256=m["source_pins_sha256"],
  resources=dict(actual_threads=c["old"],requested_threads=c["old"],actual_affinity_cpus=c["old_cpus"]))
 return p,m,arrays
def child(name,before):
 out=copy.deepcopy(before);m=json.loads(str(out["metadata_json"].item()));m["parameter_identity"]["threads"]=r.case(name)["new"]
 out["metadata_json"]=np.array(json.dumps(m));return out
def inputs(name):
 c=r.case(name);mans={}
 for k in c["prerequisites"]:
  if k=="fine_moved":mans[k]=dict(status="RUNNING",run_complete=False,resources=dict(actual_threads=4,requested_threads=4,actual_affinity_cpus=[0,1,2,3]))
  else:mans[k]=dict(status="DEVELOPMENT_DURATION_REACHED" if k=="wide_terminal" else "TRIGGER_MARGIN_COMPLETE",run_complete=True)
 return mans,[],dict(physical_core_groups=[[k,k+16] for k in range(16)]),dict(selected_affinity_cpus=[x for k in range(14) for x in (k,k+16)])
class State(unittest.TestCase):
 def test_actual_producer_namespace_and_all_three_changes(self):
  for name in r.CASES:
   p,m,b=fixture(name)
   with self.subTest(case=name):
    self.assertEqual(len(b),28);self.assertEqual(len(b)-1,27);self.assertIn("trajectory__time_s",b);self.assertNotIn("trajectory__time_ps",b)
    self.assertTrue(r.parent_contract(name,p,m,b));old,new=r.transition_arrays(name,b,child(name,b))
    self.assertEqual(old["guard_active_since_time_ps"],new["guard_active_since_time_ps"])
 def test_all27_arrays_byte_exact_including_signed_zero(self):
  for name in r.CASES:
   _,_,b=fixture(name)
   for key in b:
    if key=="metadata_json":continue
    c=child(name,b);c[key].flat[0]+=1
    with self.subTest(case=name,key=key),self.assertRaises(ValueError):r.transition_arrays(name,b,c)
   for key in ("p","trajectory__time_s","snapshots__delta_J"):
    b[key].flat[0]=0.;c=child(name,b);c[key].flat[0]=-0.
    with self.subTest(case=name,signed_zero=key),self.assertRaises(ValueError):r.transition_arrays(name,b,c)
 def test_same_cardinality_wrong_names_and_dtype_fail(self):
  for name in r.CASES:
   _,_,b=fixture(name);c=child(name,b);c["trajectory__time_ps"]=c.pop("trajectory__time_s")
   with self.subTest(case=name),self.assertRaises(ValueError):r.transition_arrays(name,b,c)
   c=child(name,b);c["p"]=c["p"].astype(np.float32)
   with self.assertRaises(ValueError):r.transition_arrays(name,b,c)
 def test_all_metadata_and_physical_parameters_except_threads_preserved(self):
  keys=dict(guard_active_since_time_ps=30.,max_guarded_accepted_local_energy_fraction=0.,initial_total_energy_J=0.,
   accepted_steps=0,total_rejections=0,next_step_ps=.01,external_energy_integral_J=0.,
   absolute_external_energy_integral_J=0.,sum_abs_local_energy_defect_J=0.,source={},controller={},source_pins_sha256={})
  for name in r.CASES:
   _,_,b=fixture(name)
   for key,value in keys.items():
    c=child(name,b);m=json.loads(str(c["metadata_json"].item()));m[key]=value;c["metadata_json"]=np.array(json.dumps(m))
    with self.subTest(case=name,key=key),self.assertRaises(ValueError):r.transition_arrays(name,b,c)
   for key,value in dict(step_tolerance=2e-8,dt_ps=.1,local_energy_tolerance=1e-5,duration_ps=40.,wavelength_nm=1550.).items():
    c=child(name,b);m=json.loads(str(c["metadata_json"].item()));m["parameter_identity"][key]=value;c["metadata_json"]=np.array(json.dumps(m))
    with self.subTest(case=name,parameter=key),self.assertRaises(ValueError):r.transition_arrays(name,b,c)
 def test_types_and_signed_zero_metadata_fail(self):
  name="domain_2to4";_,_,b=fixture(name)
  for key,value in (("accepted_steps",4.),("external_energy_integral_J",-0.)):
   c=child(name,b);m=json.loads(str(c["metadata_json"].item()));m[key]=value;c["metadata_json"]=np.array(json.dumps(m))
   with self.assertRaises(ValueError):r.transition_arrays(name,b,c)
 def test_pretrigger_and_fixed_parent_resources_required(self):
  for name in r.CASES:
   p,m,b=fixture(name)
   for key in f.EVENTS:
    bad=copy.deepcopy(m);bad["controller"]["state"][key]=1.
    with self.subTest(case=name,event=key),self.assertRaises(ValueError):r.parent_contract(name,p,bad,b)
   for key,value in (("status","RUNNING"),("run_complete",True)):
    bad=copy.deepcopy(p);bad[key]=value
    with self.assertRaises(ValueError):r.parent_contract(name,bad,m,b)
   bad=copy.deepcopy(m);bad["parameter_identity"]["threads"]=True
   with self.assertRaises(ValueError):r.parent_contract(name,p,bad,b)
 def test_resource_chronology2_to4_to8_and_guard_anchor_unchanged(self):
  _,m,_=fixture("spectral_2to4")
  original=dict(segments=[dict(start_time_ps=0.,end_time_ps=22.3,threads=2,local_energy_tolerance=1e-5,guard_active_since_time_ps=12.4),
   dict(start_time_ps=22.3,end_time_ps=None,threads=2,local_energy_tolerance=2.5e-6,guard_active_since_time_ps=22.3,step_tolerance=2e-7,dt_ps=.02)])
  four=r.resource_segments("spectral_2to4",original,m)
  projection=r.guard_projection(original)
  proof=dict(segments=projection,complete_resource_and_numerical_segments=four)
  m["time_ps"]=35.;eight=r.resource_segments("spectral_4to8",proof,m)
  self.assertEqual([x["threads"] for x in eight],[2,2,4,8])
  self.assertEqual(eight[1]["end_time_ps"],30.);self.assertEqual(eight[2]["end_time_ps"],35.)
  self.assertEqual(eight[-1]["guard_active_since_time_ps"],22.3);self.assertEqual(r.guard_projection(proof),projection)
  self.assertEqual(original["segments"][-1]["end_time_ps"],None)
 def test_prepared_manifest_retains_observed_parent_resources(self):
  for name in r.CASES:
   p,m,b=fixture(name);identity=copy.deepcopy(m["parameter_identity"]);identity["threads"]=r.case(name)["new"]
   out=r.prepared_manifest(p,identity,dict(execution_started=False))
   self.assertEqual(out["resources"],p["resources"]);self.assertFalse(out["explicit_resource_transition"]["execution_started"])
class Budget(unittest.TestCase):
 def test_closed_cases_and_correct_targets(self):
  for name in r.CASES:
   mans,owners,inv,budget=inputs(name);proof=r.resource_preflight(name,mans,owners,inv,budget)
   self.assertEqual(proof["target_affinity_cpus"],r.case(name)["cpus"]);self.assertEqual(proof["reserved_physical_cores"],[14,15])
  for name in ("mesh","temporal","domain_4to8","spectral_2to8"):
   with self.assertRaises(ValueError):r.case(name)
 def test_prerequisites_must_really_finish_or_fine_really_move(self):
  for name in r.CASES:
   mans,o,i,b=inputs(name)
   for key in mans:
    bad=copy.deepcopy(mans)
    if key=="fine_moved":bad[key]["status"]="EXPLICIT_RESOURCE_TRANSITION_PREPARED"
    else:bad[key]["run_complete"]=False
    with self.subTest(case=name,key=key),self.assertRaises(ValueError):r.resource_preflight(name,bad,o,i,b)
   bad=copy.deepcopy(mans);key=next(iter(bad));bad.pop(key)
   with self.assertRaises(ValueError):r.resource_preflight(name,bad,o,i,b)
 def test_parent_child_target_prerequisite_owners_fail(self):
  for name in r.CASES:
   mans,o,i,b=inputs(name);c=r.case(name)
   for path,aff in ((str(r.S/c["parent"]),c["old_cpus"]),(str(r.S/c["child"]),c["cpus"]),("/other",[c["cpus"][0]+16]),("/reserved",[14])):
    with self.subTest(case=name,path=path),self.assertRaises(ValueError):r.resource_preflight(name,mans,[dict(output_dir=path,affinity=aff)],i,b)
   for key,path in c["prerequisites"].items():
    if key!="fine_moved":
     with self.assertRaises(ValueError):r.resource_preflight(name,mans,[dict(output_dir=str(r.S/path),affinity=[0])],i,b)
 def test_full_host_budget_no_outer_narrow_affinity(self):
  for name in r.CASES:
   mans,o,i,b=inputs(name)
   for cpus in ([4,5],list(range(32))):
    bad=dict(selected_affinity_cpus=cpus)
    with self.assertRaises(ValueError):r.resource_preflight(name,mans,o,i,bad)
 def test_num8_stays_blocked_and_num4_has_no_num8_alias(self):
  pins={}
  with patch.object(f,"certificate",return_value={"scope":"Num124 only"}):
   with self.assertRaises(ValueError):r.thread_evidence("spectral_4to8",pins)
   for name in ("domain_2to4","spectral_2to4"):
    self.assertEqual(r.thread_evidence(name,pins),{"scope":"Num124 only"})
    with self.assertRaises(ValueError):r.thread_evidence(name,pins,"fake","fake")
class Wiring(unittest.TestCase):
 def test_two_matching_num8_exceptions_must_still_match_reviewed_num4(self):
  expected=dict(error_class="ValueError",error="Supplied heating exhausts electron support")
  self.assertTrue(n8.require_reviewed_exception(expected,dict(expected),dict(expected)))
  for changed in (dict(error_class="RuntimeError",error=expected["error"]),dict(error_class="ValueError",error="different support")):
   with self.subTest(changed=changed),self.assertRaises(ValueError):n8.require_reviewed_exception(changed,dict(changed),expected)

 def test_full_python_runner_vector_preserves_args_only_contract(self):
  for name in r.CASES:
   c=r.case(name);output=r.S/c["child"];args=["--output-dir",str(output),"--resume","--cpu-start",str(c["cpu_start"]),"--threads",str(c["new"])]
   old=copy.deepcopy(args)
   with patch.object(r,"dependencies",return_value=r.dependencies()),patch.object(f,"digest",return_value=f.RUNNER_SHA):
    full=r.execution_argv(name,args)
    self.assertEqual(full,[sys.executable,"-u",str(output/"runtime_sources/pipelines/05_effective_photon_run.py"),*args])
   self.assertEqual(args,old)
   with patch.object(r,"dependencies",return_value=r.dependencies()),patch.object(f,"digest",return_value="tampered"),self.assertRaises(ValueError):r.execution_argv(name,args)
   with self.assertRaises(ValueError):r.execution_argv(name,["--threads",str(c["new"])])

 def test_default_and_fine_paths_unchanged(self):
  with patch.object(d.guard,"compare_with_bounded_guard",return_value=("baseline",{},{})) as call:
   self.assertEqual(d.compare_originals("domain",1,2,3),("baseline",{},{}));call.assert_called_once()
  class Adapter:
   def compare_with_resource_tail(self,*args,**kwargs):return args,kwargs
  with patch.object(d,"reviewed_analysis_module",return_value=Adapter()) as load:
   args,kw=d.compare_originals("mesh",1,2,3,resource_transition_receipt_path="fine")
   self.assertEqual(kw,dict(resource_receipt_path="fine"));self.assertEqual(load.call_args.args[0],"fine_resource_guard_adapter")
 def test_new_selector_is_explicit_domain_spectral_only(self):
  class Adapter:
   def compare_with_resource_tail(self,*args,**kwargs):return args,kwargs
  for kind in ("domain","spectral"):
   with patch.object(d,"reviewed_analysis_module",return_value=Adapter()) as load:
    args,kw=d.compare_originals(kind,1,2,3,control_resource_transition_receipt_path="control")
    self.assertEqual(args[0],kind);self.assertEqual(kw,dict(resource_receipt_path="control"));self.assertEqual(load.call_args.args[0],"control_resource_transition")
  for kind in ("mesh","temporal"):
   with self.assertRaises(ValueError):d.compare_originals(kind,1,2,3,control_resource_transition_receipt_path="wrong")
  with self.assertRaises(ValueError):d.compare_originals("spectral",1,2,3,resource_transition_receipt_path="fine",control_resource_transition_receipt_path="control")
 def test_finite_view_preserves_full_originals_failure_and_completion(self):
  ref=dict(manifest=dict(run_complete=True),history={},snap={});cand=dict(manifest=dict(run_complete=False),history={},snap={})
  control=dict(reference=dict(segments=[]),candidate=dict(segments=[]),guard_law=dict(nominal_anchor_time_ps=22.435),resource_only_continuation_explicit=True)
  full=dict(declared_numerical_control=control,passed=False)
  raw=dict(checks=dict(voltage_waveform=False),legacy_nodal_policy_checks=dict(voltage_waveform=False),forbidden_parameter_changes=dict(local_energy_tolerance=[1e-5,2.5e-6]),common_trajectory_interval_ps=[0.,30.])
  with patch.object(d,"compare_originals",return_value=(full,{},{})) as upstream,patch.object(d,"correct_physical_gap",side_effect=lambda x,y,z,*args:(x,y,z)),patch.object(d,"apply_verified_spectral_policy",side_effect=lambda kind,result,*args,**kwargs:result),patch.object(d.audit,"compare_runs",return_value=(raw,{},{})),patch.object(d.guard,"compare_guard_law",return_value=dict(passed=True)):
   compare,_=d.prepare_verified_comparator("domain",ref,cand,"protocol",control_resource_transition_receipt_path="control")
   upstream.assert_called_once_with("domain",ref,cand,"protocol",resource_transition_receipt_path=None,control_resource_transition_receipt_path="control")
   result,_,_=compare("domain",dict(ref),dict(cand))
   self.assertFalse(result["passed"]);self.assertFalse(result["checks"]["voltage_waveform"]);self.assertFalse(cand["manifest"]["run_complete"])
   self.assertEqual(result["declared_numerical_control"]["full_original_comparison"],full)
   bad=dict(cand);bad["history"]={}
   with self.assertRaises(ValueError):compare("domain",dict(ref),bad)
 def test_transfer_approval_required_before_any_audit(self):
  with patch.object(f,"owned",side_effect=lambda p:p),patch.object(f,"read",side_effect=[{},{}]),patch.object(f,"digest",return_value="sha"):
   with self.assertRaises(ValueError):r.launch_preflight("domain_2to4","receipt","protocol","approval")
 def test_num8_separate_approval_does_not_infer_admission(self):
  actual=dict(current_runtime_pins_sha256={"x":"y"});approve=dict(status="ROOT_APPROVED_NUM8_RESOURCE_EQUIVALENCE",certificate_sha256="sha",source_pins_sha256={"x":"y"},stage5=False)
  with patch.object(n8,"comparison",return_value=actual),patch.object(f,"owned",side_effect=lambda p:p),patch.object(f,"read",side_effect=[actual,dict(approve,status="OTHER")]),patch.object(f,"digest",return_value="sha"):
   with self.assertRaises(ValueError):n8.verify("certificate",{"x":"y"},"approval")

class CapturedStateTreeCopy(unittest.TestCase):
 def trees(self,tmp,archived=True):
  import hashlib
  base=Path(tmp).resolve();parent=base/"parent";anchor=base/"anchor";output=base/"output"
  relative=Path("runtime_sources/docs/implementation/stage5/prephoton_dc_20260924/meshes/l320_h5/mesh.npz") if archived else Path("meshes/l320_h5/mesh.npz")
  source=parent/relative;source.parent.mkdir(parents=True);source.write_bytes(b"exact captured DC mesh")
  captured=anchor/relative;captured.parent.mkdir(parents=True);captured.write_bytes(source.read_bytes())
  (anchor/"runtime_sources").mkdir(exist_ok=True);(anchor/"dc").mkdir()
  (anchor/"dc/bulk.npz").write_bytes(b"saved DC state");output.mkdir()
  sha=hashlib.sha256(source.read_bytes()).hexdigest()
  entries={str(relative.relative_to("runtime_sources")):dict(sha256=sha,bytes=source.stat().st_size)} if archived else {}
  (anchor/"runtime_sources/archive_manifest.json").write_text(json.dumps(dict(package_sources_and_small_geometry=entries)))
  geometry=dict(source_path=str(source),mesh_sha256=sha,n_nodes=7,new_mesh_generated=False)
  return base,parent,anchor,output,relative,geometry
 def test_missing_meshes_archived_copy_is_complete_and_authenticated(self):
  import tempfile
  from sandbox.stage5_effective_validation import fine_resource_guard_adapter as adapter
  with tempfile.TemporaryDirectory(prefix="resource_copy_") as tmp:
   base,parent,anchor,output,rel,geometry=self.trees(tmp)
   self.assertFalse((anchor/"meshes").exists())
   with patch.object(adapter.guard,"SCRATCH",base):r.copy_prepared_state_trees(anchor,output,geometry,parent)
   self.assertEqual((output/rel).read_bytes(),(anchor/rel).read_bytes())
   self.assertEqual((output/"dc/bulk.npz").read_bytes(),b"saved DC state");self.assertFalse((output/"meshes").exists())
 def test_generated_meshes_still_copy_and_missing_generated_tree_fails(self):
  import tempfile,shutil
  from sandbox.stage5_effective_validation import fine_resource_guard_adapter as adapter
  for missing in (False,True):
   with self.subTest(missing=missing),tempfile.TemporaryDirectory(prefix="resource_copy_") as tmp:
    base,parent,anchor,output,rel,geometry=self.trees(tmp,archived=False)
    if missing:shutil.rmtree(anchor/"meshes")
    with patch.object(adapter.guard,"SCRATCH",base):
     if missing:
      with self.assertRaises(ValueError):r.copy_prepared_state_trees(anchor,output,geometry,parent)
     else:
      r.copy_prepared_state_trees(anchor,output,geometry,parent)
      self.assertEqual((output/rel).read_bytes(),(parent/rel).read_bytes())
 def test_runtime_and_dc_remain_mandatory_with_archived_geometry(self):
  import tempfile,shutil
  from sandbox.stage5_effective_validation import fine_resource_guard_adapter as adapter
  for absent in ("runtime_sources","dc"):
   with self.subTest(absent=absent),tempfile.TemporaryDirectory(prefix="resource_copy_") as tmp:
    base,parent,anchor,output,rel,geometry=self.trees(tmp);shutil.rmtree(anchor/absent)
    with patch.object(adapter.guard,"SCRATCH",base),self.assertRaises(ValueError):r.copy_prepared_state_trees(anchor,output,geometry,parent)
 def test_archived_capture_or_archive_entry_corruption_fails(self):
  import tempfile
  from sandbox.stage5_effective_validation import fine_resource_guard_adapter as adapter
  for corrupt in ("captured_mesh","archive_sha","archive_missing"):
   with self.subTest(corrupt=corrupt),tempfile.TemporaryDirectory(prefix="resource_copy_") as tmp:
    base,parent,anchor,output,rel,geometry=self.trees(tmp)
    if corrupt=="captured_mesh":(anchor/rel).write_bytes(b"wrong captured mesh")
    else:
     ap=anchor/"runtime_sources/archive_manifest.json";ar=json.loads(ap.read_text())
     if corrupt=="archive_sha":ar["package_sources_and_small_geometry"][str(rel.relative_to("runtime_sources"))]["sha256"]="0"*64
     else:ar["package_sources_and_small_geometry"]={}
     ap.write_text(json.dumps(ar))
    with patch.object(adapter.guard,"SCRATCH",base),self.assertRaises(ValueError):r.copy_prepared_state_trees(anchor,output,geometry,parent)
 def test_unsupported_archive_root_or_symlink_cannot_hide_missing_generated_tree(self):
  import tempfile
  from sandbox.stage5_effective_validation import fine_resource_guard_adapter as adapter
  for bad in ("other_root","symlink"):
   with self.subTest(bad=bad),tempfile.TemporaryDirectory(prefix="resource_copy_") as tmp:
    base,parent,anchor,output,rel,geometry=self.trees(tmp)
    if bad=="other_root":geometry["source_path"]=str(parent/"runtime_sources/docs/other/mesh.npz")
    else:
     mesh=anchor/rel;payload=mesh.read_bytes();mesh.unlink();target=mesh.parent/"real.npz";target.write_bytes(payload);mesh.symlink_to(target)
    with patch.object(adapter.guard,"SCRATCH",base),self.assertRaises(ValueError):r.copy_prepared_state_trees(anchor,output,geometry,parent)
 def test_prepare_reaches_next_stage_with_missing_meshes_through_real_copy_helper(self):
  # Run the actual prepare control flow through copying; mock prior audited validators
  # and stop at the next material-copy stage before any numerical checkpoint is written.
  import tempfile,hashlib,shutil,types
  from sandbox.stage5_effective_validation import fine_resource_guard_adapter as adapter
  class AfterStateTreeCopy(Exception):pass
  with tempfile.TemporaryDirectory(prefix="resource_prepare_copy_") as tmp:
   base,parent0,anchor,output0,rel,geometry=self.trees(tmp)
   case=r.case("spectral_2to4");parent=base/case["parent"];parent0.rename(parent);geometry["source_path"]=str(parent/rel)
   output=base/case["child"];output0.rmdir()
   pm,meta,arrays=fixture("spectral_2to4");pm["geometry"]=geometry
   payload=b"synthetic audited immutable checkpoint";mb=json.dumps(pm).encode()
   (parent/"checkpoint_latest.npz").write_bytes(payload);(parent/"manifest.json").write_bytes(mb)
   (anchor/"checkpoint_original.npz").write_bytes(payload);(anchor/"manifest_original.json").write_bytes(mb)
   cr=dict(parent_output_dir=str(parent),checkpoint_bytes_sha256=hashlib.sha256(payload).hexdigest(),manifest_sha256=hashlib.sha256(mb).hexdigest(),physical_file_sha256={})
   (anchor/"capture_receipt.json").write_text(json.dumps(cr))
   for name in ("trajectory.npz","snapshots.npz"):(parent/name).write_bytes(name.encode())
   hh=types.SimpleNamespace(decode=lambda raw:arrays,audit_arrays=lambda arr,manifest:meta,audit_runtime=lambda *args:dict(manifest_sha256="fixture"),
    writejson=lambda p,v:p.write_text(json.dumps(v)),newdir=lambda p:(Path(p).mkdir(),Path(p))[1],lineage=types.SimpleNamespace(segments_valid=lambda x:True))
   gg=types.SimpleNamespace(h=hh,eligible=lambda *args:True)
   gd=types.SimpleNamespace(load_protocol=lambda *args:())
   audited=dict(summary=dict(checks=dict(physical=True)),manifest=copy.deepcopy(pm),history={},snap={},hotbelt={},traces={},spatial={})
   self.assertNotIn("geometry",audited)
   audit=types.SimpleNamespace(audit_run=lambda *args:audited)
   resources=dict(prerequisite_manifests={})
   copy2=shutil.copy2
   def stop_after_copy(src,dst,*args,**kwargs):
    if Path(src).name=="material.npz":raise AfterStateTreeCopy()
    return copy2(src,dst,*args,**kwargs)
   with patch.object(r,"S",base),patch.object(adapter.guard,"SCRATCH",base),patch.object(r,"dependencies",return_value=(f,adapter,gg,gd,audit)),patch.object(f,"owned",side_effect=lambda p:Path(p).resolve()),patch.object(f,"audit_manifest_view",return_value=({},{})),patch.object(r,"preflight",return_value=resources),patch.object(r,"thread_evidence",return_value={}),patch.object(r,"parent_proof",return_value={}),patch.object(r,"resource_segments",return_value=[]),patch.object(r.os,"sched_setaffinity"),patch.object(r.shutil,"copy2",side_effect=stop_after_copy),patch.object(r,"copy_prepared_state_trees",wraps=r.copy_prepared_state_trees) as actual_copy:
    with self.assertRaises(AfterStateTreeCopy):r.prepare("spectral_2to4",anchor,output,"fixture_protocol")
    actual_copy.assert_called_once_with(anchor,output,geometry,parent)
    audited["manifest"]["geometry"]["n_nodes"]+=1
    with self.assertRaisesRegex(ValueError,"Audited geometry differs"):
     r.prepare("spectral_2to4",anchor,output,"fixture_protocol")
    actual_copy.assert_called_once()
   self.assertEqual((output/rel).read_bytes(),(anchor/rel).read_bytes())
   self.assertTrue((output/"dc/bulk.npz").exists());self.assertFalse((output/"meshes").exists())
   self.assertFalse((output/"provenance/transfer_receipt.json").exists())

if __name__=="__main__":unittest.main()

# These contracts inspect the immutable historical checkpoint producer.
# The current 1.1.2 runner must not be substituted into an old namespace.
import pytest
@pytest.fixture(autouse=True)
def historical_producer_schema(monkeypatch):
    frozen_root=Path(__file__).resolve().parent/"fixtures/effective_v1_0_producer"
    monkeypatch.setattr(f,"ROOT",frozen_root)
