"""Portable semantic tests for the explicit fine resource-only bridge."""
import copy,hashlib,importlib.util,json,sys,unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
BASE=Path(__file__).resolve().parent
def module(name):
 p=BASE/(name+".py")
 if not p.exists():p=BASE.parent/"sandbox/stage5_effective_validation"/(name+".py")
 spec=importlib.util.spec_from_file_location("sandbox.stage5_effective_validation._candidate_"+name,p)
 m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m
r=module("fine_resource_transition");a=module("fine_resource_guard_adapter");d=module("guard_protocol_dispatcher");c=module("control_resource_transition")
def fixture():
 i=dict(threads=2,local_energy_tolerance=2.5e-6,step_tolerance=2e-7,dt_ps=.02,
  length_nm=320.,spacing_nm=3.5,wavelength_nm=775.,xi_cutoff=32.,xi_order=6,phonon_nodes=96,transport_order=48,duration_ps=None)
 m=dict(parameter_identity=i,source_pins_sha256={str(n):str(n)*64 for n in range(13)},
  time_ps=25.,accepted_steps=4,total_rejections=2,next_step_ps=.0025,
  initial_total_energy_J=1e-19,external_energy_integral_J=2e-20,source=dict(retained_energy_J=1e-19),
  controller=dict(state=dict(status="ACTIVE",rejected_candidate_count=0,**{k:None for k in r.EVENTS})),
  guard_active_since_time_ps=22.3,local_energy_reference_energy_J=1e-19,energy_guard_ledger_origin="original",
  local_energy_rejections=1,limited_event_energy_J=1e-22,totalabs_event_energy_J=1e-20,
  absolute_external_energy_integral_J=2e-20,sum_abs_local_energy_defect_J=1e-25,
  max_accepted_local_energy_fraction=1e-6,max_guarded_accepted_local_energy_fraction=1e-6)
 arrays={k:np.array([1.,2.]) for k in r.STATE_KEYS}
 arrays.update({name:np.array([1.,2.]) for name in r.checkpoint_namespace() if name.startswith("trajectory__")})
 arrays.update({name:np.array([[1.,2.]]) for name in r.checkpoint_namespace() if name.startswith("snapshots__")})
 arrays["metadata_json"]=np.array(json.dumps(m))
 manifest=dict(status="INTERRUPTED",run_complete=False,parameter_identity=copy.deepcopy(i),source_pins_sha256=m["source_pins_sha256"],
  resources=dict(actual_threads=2,requested_threads=2,actual_affinity_cpus=[12,13]),
  physicalparams=dict(Tb=.9),source=m["source"],energy_reference=dict(initial_total_energy_J=1e-19))
 return manifest,m,arrays
def child(before):
 out=copy.deepcopy(before);m=json.loads(str(out["metadata_json"].item()));m["parameter_identity"]["threads"]=4
 out["metadata_json"]=np.array(json.dumps(m));return out
class Contracts(unittest.TestCase):
 def test_exact28_keys_27_numeric_and_sole_threads_change(self):
  p,m,b=fixture();self.assertEqual(len(b),28);self.assertEqual(len(b)-1,27);self.assertTrue(r.parent_contract(p,m,b))
  old,new=r.transition_arrays(b,child(b));self.assertEqual(old["next_step_ps"],new["next_step_ps"])
 def test_complete_namespace_rejects_same_cardinality_wrong_field(self):
  _,_,b=fixture();c=child(b);c["trajectory__time_ps"]=c.pop("trajectory__time_s")
  self.assertEqual(len(c),28)
  with self.assertRaises(ValueError):r.transition_arrays(b,c)
 def test_all_metadata_must_be_preserved(self):
  _,m,b=fixture()
  for key,value in (("next_step_ps",.001),("initial_total_energy_J",0.),("total_rejections",0),
   ("accepted_steps",0),("guard_active_since_time_ps",25.),("max_guarded_accepted_local_energy_fraction",0.),
   ("absolute_external_energy_integral_J",0.),("source_pins_sha256",{}),("controller",{})):
   with self.subTest(key=key):
    c=child(b);cm=json.loads(str(c["metadata_json"].item()));cm[key]=value;c["metadata_json"]=np.array(json.dumps(cm))
    with self.assertRaises(ValueError):r.transition_arrays(b,c)
 def test_solver_and_physics_cannot_change(self):
  _,m,b=fixture()
  for key,value in (("local_energy_tolerance",1e-5),("step_tolerance",2e-8),("dt_ps",.04),("spacing_nm",5.),("duration_ps",20.)):
   with self.subTest(key=key):
    c=child(b);cm=json.loads(str(c["metadata_json"].item()));cm["parameter_identity"][key]=value;c["metadata_json"]=np.array(json.dumps(cm))
    with self.assertRaises(ValueError):r.transition_arrays(b,c)
 def test_every_physical_and_saved_array_checked(self):
  _,m,b=fixture()
  for key in b:
   if key=="metadata_json":continue
   with self.subTest(key=key):
    c=child(b);c[key]=c[key]+1
    with self.assertRaises(ValueError):r.transition_arrays(b,c)
 def test_signed_zero_is_not_byte_equal(self):
  for key in ("p","trajectory__time_s"):
   p,m,b=fixture();b[key][0]=0.;c=child(b);c[key][0]=-0.
   with self.subTest(key=key),self.assertRaises(ValueError):r.transition_arrays(b,c)
 def test_signed_zero_saved_prefix_rejected(self):
  _,_,b=fixture();b["snapshots__delta_J"][0,0]=0.;c=child(b);live=copy.deepcopy(c);live["snapshots__delta_J"][0,0]=-0.
  with self.assertRaises(ValueError):r.live_prefix(c,live)
 def test_metadata_numeric_types_and_signed_zero_preserved(self):
  _,_,b=fixture();bm=json.loads(str(b["metadata_json"].item()));bm["synthetic_zero"]=0.;b["metadata_json"]=np.array(json.dumps(bm))
  for key,value in (("synthetic_zero",-0.),("total_rejections",2.)):
   c=child(b);cm=json.loads(str(c["metadata_json"].item()));cm[key]=value;c["metadata_json"]=np.array(json.dumps(cm))
   with self.subTest(key=key),self.assertRaises(ValueError):r.transition_arrays(b,c)
 def test_dtype_and_namespace_changes_fail(self):
  _,_,b=fixture();c=child(b);c["p"]=c["p"].astype(np.float32)
  with self.assertRaises(ValueError):r.transition_arrays(b,c)
  c=child(b);c.pop("p")
  with self.assertRaises(ValueError):r.transition_arrays(b,c)
 def test_parent_strict_status_event_resources_and_anchor(self):
  p,m,b=fixture()
  for key,val in (("status","RUNNING"),("run_complete",True)):
   q=copy.deepcopy(p);q[key]=val
   with self.assertRaises(ValueError):r.parent_contract(q,m,b)
  for key in r.EVENTS:
   q=copy.deepcopy(m);q["controller"]["state"][key]=1.
   with self.subTest(key=key),self.assertRaises(ValueError):r.parent_contract(p,q,b)
  for key,val in (("threads",4),("threads",True),("local_energy_tolerance",1e-5)):
   q=copy.deepcopy(m);q["parameter_identity"][key]=val
   with self.assertRaises(ValueError):r.parent_contract(p,q,b)
  q=copy.deepcopy(m);q["guard_active_since_time_ps"]=12.4
  with self.assertRaises(ValueError):r.parent_contract(p,q,b)
 def test_advanced_prefix_preserves_guard_origin_and_history(self):
  _,_,b=fixture();c=child(b);self.assertEqual(r.live_prefix(c,c)["time_ps"],25.)
  for key in ("trajectory__time_s","snapshots__delta_J"):
   live=copy.deepcopy(c);live[key]=live[key]+1
   with self.assertRaises(ValueError):r.live_prefix(c,live)
 def test_same_clock_physical_state_cannot_change(self):
  _,_,b=fixture();c=child(b);live=copy.deepcopy(c);live["delta_bar"]*=.5
  with self.assertRaises(ValueError):r.live_prefix(c,live)
 def test_advanced_cumulative_counters_cannot_reset(self):
  _,_,b=fixture();c=child(b)
  for key in ("total_rejections","local_energy_rejections","sum_abs_local_energy_defect_J"):
   live=copy.deepcopy(c);m=json.loads(str(live["metadata_json"].item()));m["time_ps"]+=1;m["accepted_steps"]+=1;m[key]=0
   live["metadata_json"]=np.array(json.dumps(m))
   with self.subTest(key=key),self.assertRaises(ValueError):r.live_prefix(c,live)
 def test_full_resource_timeline_and_separate_guard_projection(self):
  _,m,_=fixture();rows=[dict(start_time_ps=0.,end_time_ps=22.3,threads=2,local_energy_tolerance=1e-5,
    guard_active_since_time_ps=12.4,origin="old",step_tolerance=2e-7,dt_ps=.02),
   dict(start_time_ps=22.3,end_time_ps=None,threads=2,local_energy_tolerance=2.5e-6,
    guard_active_since_time_ps=22.3,origin="prospective_local_energy_guard_tightening",step_tolerance=2e-7,dt_ps=.02)]
  full=r.resource_segments(rows,m);self.assertEqual(full[-1]["threads"],4)
  self.assertEqual(full[-1]["guard_active_since_time_ps"],22.3)
  projected=r.guard_projection(full,25.);self.assertNotIn("threads",projected[-1]);self.assertEqual(projected[-1]["end_time_ps"],None)
  for key,val in (("step_tolerance",2e-8),("local_energy_tolerance",1e-5),("guard_active_since_time_ps",25.)):
   bad=copy.deepcopy(full);bad[-1][key]=val
   with self.subTest(key=key),self.assertRaises(ValueError):r.guard_projection(bad,25.)
 def test_prepared_manifest_does_not_invent_observed_num4_execution(self):
  p,m,b=fixture();i=copy.deepcopy(m["parameter_identity"]);i["threads"]=4
  out=r.prepared_manifest(p,i,dict(requested_threads=4,execution_started=False))
  self.assertEqual(out["resources"]["actual_threads"],2);self.assertEqual(out["parameter_identity"]["threads"],4)
  self.assertFalse(out["explicit_resource_transition"]["execution_started"])
  self.assertEqual(out["physicalparams"],p["physicalparams"]);self.assertFalse(out["run_complete"])
class Resources(unittest.TestCase):
 def inputs(self):
  return dict(run_complete=True,status="TRIGGER_MARGIN_COMPLETE"),[],dict(physical_core_groups=[[n,n+16] for n in range(16)]),dict(selected_affinity_cpus=[n for c in range(14) for n in (c,c+16)])
 def test_nominal_completion_and_cores_really_free(self):
  n,o,i,b=self.inputs();self.assertEqual(r.resource_preflight(n,o,i,b)["target_affinity_cpus"],[0,1,2,3])
  n["run_complete"]=False
  with self.assertRaises(ValueError):r.resource_preflight(n,o,i,b)
 def test_parent_owner_or_target_occupancy_rejects(self):
  n,o,i,b=self.inputs()
  for owner in (dict(output_dir=str(r.PARENT),affinity=[12,13]),dict(output_dir="/other",affinity=[16])):
   with self.subTest(owner=owner),self.assertRaises(ValueError):r.resource_preflight(n,[owner],i,b)
 def test_reserve14_15_and_full_host_admission(self):
  n,o,i,b=self.inputs();b["selected_affinity_cpus"]=[0,1,2,3]
  with self.assertRaises(ValueError):r.resource_preflight(n,o,i,b)
  b["selected_affinity_cpus"]=list(range(32))
  with self.assertRaises(ValueError):r.resource_preflight(n,o,i,b)
 def test_global_ceiling_counts_nested_pipeline_affinity(self):
  n,o,i,b=self.inputs()
  for c in range(4,16):o.append(dict(output_dir="/other"+str(c),affinity=[c]))
  with self.assertRaises(ValueError):r.resource_preflight(n,o,i,b)
class Wiring(unittest.TestCase):
 def test_global_one_percent_and_independent_guard_unchanged(self):
  summary=dict(energy=dict(relative_limit=.01,maximum_relative_defect=.009),checks=dict(energy_relative_defect=True,independent_local_energy_guard=True))
  self.assertTrue(a.audit_energy_contract(summary))
  for key,value in (("relative_limit",.02),("maximum_relative_defect",.011),("maximum_relative_defect",float("nan"))):
   bad=copy.deepcopy(summary);bad["energy"][key]=value
   with self.subTest(key=key),self.assertRaises(ValueError):a.audit_energy_contract(bad)
  bad=copy.deepcopy(summary);bad["checks"]["independent_local_energy_guard"]=False
  with self.assertRaises(ValueError):a.audit_energy_contract(bad)
 def test_selector_default_unchanged(self):
  with patch.object(d.guard,"compare_with_bounded_guard",return_value=("original",{},{})) as fn:
   self.assertEqual(d.compare_originals("mesh",1,2,3),("original",{},{}));fn.assert_called_once()
 def test_explicit_selector_only_mesh(self):
  for kind in ("temporal","domain","spectral"):
   with self.subTest(kind=kind),self.assertRaises(ValueError):d.compare_originals(kind,1,2,3,resource_transition_receipt_path="receipt")
 def test_selector_calls_only_resource_adapter_when_explicit(self):
  class Adapter:
   def compare_with_resource_tail(self,*args,**kwargs):return args,kwargs
  with patch.object(d,"reviewed_analysis_module",return_value=Adapter()):
   args,kwargs=d.compare_originals("mesh",1,2,3,resource_transition_receipt_path="receipt")
   self.assertEqual(args[0],"mesh");self.assertEqual(kwargs,dict(resource_receipt_path="receipt"))
 def test_finite_view_keeps_resource_proof_and_physical_failure(self):
  ref=dict(manifest=dict(run_complete=True),history={},snap={});cand=dict(manifest=dict(run_complete=False),history={},snap={})
  control=dict(reference=dict(segments=[]),candidate=dict(segments=[]),guard_law=dict(nominal_anchor_time_ps=22.435),resource_only_continuation_explicit=True)
  full=dict(declared_numerical_control=control,passed=False)
  raw=dict(checks=dict(voltage_waveform=False),legacy_nodal_policy_checks=dict(voltage_waveform=False),
   forbidden_parameter_changes=dict(local_energy_tolerance=[1e-5,2.5e-6]),common_trajectory_interval_ps=[0.,20.])
  with patch.object(d,"compare_originals",return_value=(full,{},{})) as original,patch.object(d,"correct_physical_gap",side_effect=lambda x,y,z,*args:(x,y,z)),patch.object(d,"apply_verified_spectral_policy",side_effect=lambda kind,result,*args,**kwargs:result),patch.object(d.audit,"compare_runs",return_value=(raw,{},{})),patch.object(d.guard,"compare_guard_law",return_value=dict(passed=True)):
   compare,result=d.prepare_verified_comparator("mesh",ref,cand,"protocol",resource_transition_receipt_path="receipt")
   original.assert_called_once_with("mesh",ref,cand,"protocol",resource_transition_receipt_path="receipt")
   row,_,_=compare("mesh",dict(ref),dict(cand))
   self.assertTrue(row["checks"]["verified_explicit_resource_only_continuation"])
   self.assertFalse(row["checks"]["voltage_waveform"]);self.assertFalse(row["passed"])
   self.assertFalse(cand["manifest"]["run_complete"])
   self.assertTrue(row["declared_numerical_control"]["resource_only_continuation_explicit"])
   bad=dict(cand);bad["snap"]={}
   with self.assertRaises(ValueError):compare("mesh",dict(ref),bad)
 def test_adapter_rejects_axis_alias_and_other_controls(self):
  for kind in ("fine","temporal","domain","spectral"):
   with self.subTest(kind=kind),self.assertRaises(ValueError):a.compare_with_resource_tail(kind,None,None,None,None,resource_receipt_path="receipt")
class InterruptedAuditView(unittest.TestCase):
 def inputs(self):
  return dict(status="INTERRUPTED",run_complete=False,physicalparams=dict(Tb=.9)),{k:hashlib.sha256(k.encode()).hexdigest() for k in r.AUDIT_OUTPUTS}
 def test_only_true_python_or_numpy_booleans(self):
  for value in (True,np.bool_(True)):
   self.assertTrue(r.true_check(value))
  for value in (False,np.bool_(False),0,1,"True","False",None,[],{}):
   with self.subTest(value=value):self.assertFalse(r.true_check(value))
 def test_missing_interrupt_hashes_add_only_outputs_without_mutating_original(self):
  original,hashes=self.inputs();saved=copy.deepcopy(original)
  derived,reason=r.audit_manifest_projection(original,hashes)
  expected=copy.deepcopy(saved);expected["output_sha256"]=hashes
  self.assertTrue(r.exact_json(derived,expected));self.assertTrue(r.exact_json(original,saved));self.assertIn("INTERRUPTED",reason)
 def test_existing_valid_hashes_remain_exact(self):
  original,hashes=self.inputs();original["output_sha256"]=copy.deepcopy(hashes)
  derived,reason=r.audit_manifest_projection(original,hashes)
  self.assertTrue(r.exact_json(derived,original));self.assertIn("unchanged",reason)
 def test_existing_bad_or_incomplete_hashes_never_repaired(self):
  for bad in ({},None,False,dict(trajectory="a"*64),dict((k,"0"*64) for k in r.AUDIT_OUTPUTS)):
   original,hashes=self.inputs();original["output_sha256"]=bad;saved=copy.deepcopy(original)
   with self.subTest(bad=bad),self.assertRaises(ValueError):r.audit_manifest_projection(original,hashes)
   self.assertTrue(r.exact_json(original,saved))
 def test_missing_hashes_reject_other_status_or_complete(self):
  for status,complete in (("RUNNING",False),("TRIGGER_MARGIN_COMPLETE",True),("INTERRUPTED",True),("INTERRUPTED",0)):
   original,hashes=self.inputs();original.update(status=status,run_complete=complete)
   with self.subTest(status=status,complete=complete),self.assertRaises(ValueError):r.audit_manifest_projection(original,hashes)
 def test_io_projection_authenticates_original_derived_and_physical_files(self):
  import tempfile
  with tempfile.TemporaryDirectory(prefix="resource_audit_view_") as tmp,patch.object(r,"owned",side_effect=lambda p:Path(p).resolve()):
   anchor=Path(tmp);original,_=self.inputs();original_bytes=json.dumps(original).encode()
   (anchor/"manifest_original.json").write_bytes(original_bytes)
   physical={}
   for name in r.AUDIT_OUTPUTS:
    (anchor/name).write_bytes(name.encode());physical[name]=(anchor/name).read_bytes()
   writer=lambda path,data:path.write_text(json.dumps(data,sort_keys=True))
   derived,record=r.audit_manifest_view(anchor,original,write_json=writer)
   self.assertTrue(record["output_sha256_added"]);self.assertEqual((anchor/"manifest_original.json").read_bytes(),original_bytes)
   self.assertEqual(record["original_manifest_sha256"],hashlib.sha256(original_bytes).hexdigest())
   self.assertEqual(record["derived_manifest_sha256"],r.digest(anchor/"manifest.json"))
   self.assertEqual(r.audit_manifest_view(anchor,original)[1],record)
   for name,data in physical.items():self.assertEqual((anchor/name).read_bytes(),data)
   tampered=copy.deepcopy(derived);tampered["physicalparams"]["Tb"]=1.
   writer(anchor/"manifest.json",tampered)
   with self.assertRaises(ValueError):r.audit_manifest_view(anchor,original)
 def test_io_bad_original_output_hash_is_not_overwritten(self):
  import tempfile
  with tempfile.TemporaryDirectory(prefix="resource_bad_hash_") as tmp,patch.object(r,"owned",side_effect=lambda p:Path(p).resolve()):
   anchor=Path(tmp);original,hashes=self.inputs();original["output_sha256"]=hashes
   payload=json.dumps(original).encode();(anchor/"manifest_original.json").write_bytes(payload);(anchor/"manifest.json").write_bytes(payload)
   for name in r.AUDIT_OUTPUTS:(anchor/name).write_bytes(b"mismatching actual bytes")
   with self.assertRaises(ValueError):r.audit_manifest_view(anchor,original,write_json=lambda p,data:p.write_text(json.dumps(data)))
   self.assertEqual((anchor/"manifest.json").read_bytes(),payload);self.assertEqual((anchor/"manifest_original.json").read_bytes(),payload)
class ResourceGeometryRelocation(unittest.TestCase):
 def meshes(self,tmp):
  base=Path(tmp).resolve();parent=base/"parent";child=base/"child"
  old=parent/"meshes/l320/mesh.npz";new=child/"meshes/l320/mesh.npz"
  old.parent.mkdir(parents=True);new.parent.mkdir(parents=True);old.write_bytes(b"exact saved mesh");new.write_bytes(old.read_bytes())
  geometry=dict(source_path=str(old),mesh_sha256=hashlib.sha256(old.read_bytes()).hexdigest(),n_nodes=7,length_m=3.2e-7,new_mesh_generated=False)
  relocated=copy.deepcopy(geometry);relocated["source_path"]=str(new)
  return base,parent,child,old,new,geometry,relocated
 def test_actual_relocation_and_prepared_parent_path_are_distinct(self):
  import tempfile
  with tempfile.TemporaryDirectory(prefix="resource_geometry_") as tmp:
   base,parent,child,old,new,original,relocated=self.meshes(tmp)
   with patch.object(a.guard,"SCRATCH",base):
    self.assertTrue(a.geometry_relocation(original,relocated,parent,child))
    self.assertTrue(a.geometry_relocation(original,original,parent,child,prepared=True))
    for current,prepared in ((original,False),(relocated,True)):
     with self.subTest(prepared=prepared),self.assertRaises(ValueError):a.geometry_relocation(original,current,parent,child,prepared=prepared)
 def test_child_and_parent_bytes_must_both_match_declared_sha(self):
  import tempfile
  with tempfile.TemporaryDirectory(prefix="resource_geometry_") as tmp:
   base,parent,child,old,new,original,relocated=self.meshes(tmp)
   with patch.object(a.guard,"SCRATCH",base):
    new.write_bytes(b"different mesh")
    for prepared in (False,True):
     with self.subTest(prepared=prepared),self.assertRaises(ValueError):a.geometry_relocation(original,original if prepared else relocated,parent,child,prepared=prepared)
    new.write_bytes(b"exact saved mesh");old.write_bytes(b"wrong parent")
    with self.assertRaises(ValueError):a.geometry_relocation(original,relocated,parent,child)
 def test_scalar_types_values_and_geometry_namespace_remain_exact(self):
  import tempfile
  with tempfile.TemporaryDirectory(prefix="resource_geometry_") as tmp:
   base,parent,child,old,new,original,relocated=self.meshes(tmp)
   with patch.object(a.guard,"SCRATCH",base):
    for key,value in (("n_nodes",7.),("n_nodes",8),("length_m",3.3e-7),("new_mesh_generated",0),("mesh_sha256","0"*64),("extra",1)):
     changed=copy.deepcopy(relocated);changed[key]=value
     with self.subTest(key=key,value=value),self.assertRaises(ValueError):a.geometry_relocation(original,changed,parent,child)
 def test_same_hash_does_not_allow_different_paths_or_parent_escape(self):
  import tempfile
  with tempfile.TemporaryDirectory(prefix="resource_geometry_") as tmp:
   base,parent,child,old,new,original,relocated=self.meshes(tmp)
   alternate=child/"meshes/elsewhere.npz";alternate.write_bytes(new.read_bytes())
   outside=base/"outside.npz";outside.write_bytes(old.read_bytes())
   with patch.object(a.guard,"SCRATCH",base):
    for target in (str(old),str(alternate),str(outside)):
     changed=copy.deepcopy(relocated);changed["source_path"]=target
     with self.subTest(target=target),self.assertRaises(ValueError):a.geometry_relocation(original,changed,parent,child)
    changed=copy.deepcopy(original);changed["source_path"]=str(outside)
    with self.assertRaises(ValueError):a.geometry_relocation(changed,relocated,parent,child)
 def test_prepared_still_authenticates_the_child_mesh_copy(self):
  import tempfile
  with tempfile.TemporaryDirectory(prefix="resource_geometry_") as tmp:
   base,parent,child,old,new,original,relocated=self.meshes(tmp)
   new.unlink()
   with patch.object(a.guard,"SCRATCH",base),self.assertRaises(OSError):a.geometry_relocation(original,original,parent,child,prepared=True)

class ArchivedResourceGeometryRelocation(unittest.TestCase):
 def meshes(self,tmp,relative=None):
  base=Path(tmp).resolve();parent=base/"parent";child=base/"child"
  rel=Path(relative) if relative is not None else a.ARCHIVED_DC_MESH_ROOT/"l320_h5/mesh.npz"
  old=parent/rel;new=child/rel
  old.parent.mkdir(parents=True);new.parent.mkdir(parents=True)
  old.write_bytes(b"archived immutable DC mesh");new.write_bytes(old.read_bytes())
  original=dict(source_path=str(old),mesh_sha256=hashlib.sha256(old.read_bytes()).hexdigest(),
   n_nodes=3659,length_m=3.2e-7,new_mesh_generated=False)
  resumed=copy.deepcopy(original);resumed["source_path"]=str(new)
  return base,parent,child,old,new,original,resumed
 def test_archived_root_preserves_prepared_and_resumed_path_contract(self):
  import tempfile
  with tempfile.TemporaryDirectory(prefix="archived_geometry_") as tmp:
   base,parent,child,old,new,original,resumed=self.meshes(tmp)
   original_saved=copy.deepcopy(original);resumed_saved=copy.deepcopy(resumed)
   with patch.object(a.guard,"SCRATCH",base):
    self.assertTrue(a.geometry_relocation(original,original,parent,child,prepared=True))
    self.assertTrue(a.geometry_relocation(original,resumed,parent,child))
    for current,prepared in ((original,False),(resumed,True)):
     with self.subTest(prepared=prepared),self.assertRaises(ValueError):a.geometry_relocation(original,current,parent,child,prepared=prepared)
   self.assertEqual(original,original_saved);self.assertEqual(resumed,resumed_saved)
 def test_other_runtime_docs_roots_are_not_normalized(self):
  import tempfile
  forbidden=("runtime_sources/docs/implementation/stage5/other/meshes/l320_h5/mesh.npz",
   "runtime_sources/docs/implementation/stage5/prephoton_dc_20260924/other/l320_h5/mesh.npz",
   "runtime_sources/docs/implementation/stage5/prephoton_dc_20260925/meshes/l320_h5/mesh.npz",
   "runtime_sources/meshes/l320_h5/mesh.npz","docs/implementation/stage5/prephoton_dc_20260924/meshes/l320_h5/mesh.npz",
   "meshes_other/l320_h5/mesh.npz")
  for relative in forbidden:
   with self.subTest(relative=relative),tempfile.TemporaryDirectory(prefix="archived_geometry_") as tmp:
    base,parent,child,old,new,original,resumed=self.meshes(tmp,relative)
    with patch.object(a.guard,"SCRATCH",base),self.assertRaises(ValueError):a.geometry_relocation(original,resumed,parent,child)
 def test_archived_hashes_required_for_both_original_and_child_even_prepared(self):
  import tempfile
  for target in ("parent","child","declared"):
   for prepared in (False,True):
    with self.subTest(target=target,prepared=prepared),tempfile.TemporaryDirectory(prefix="archived_geometry_") as tmp:
     base,parent,child,old,new,original,resumed=self.meshes(tmp)
     if target=="parent":old.write_bytes(b"corrupt parent")
     elif target=="child":new.write_bytes(b"corrupt child")
     else:original["mesh_sha256"]=resumed["mesh_sha256"]="0"*64
     with patch.object(a.guard,"SCRATCH",base),self.assertRaises(ValueError):
      a.geometry_relocation(original,original if prepared else resumed,parent,child,prepared=prepared)
 def test_archived_same_hash_cannot_change_relative_child_path_or_root_class(self):
  import tempfile
  with tempfile.TemporaryDirectory(prefix="archived_geometry_") as tmp:
   base,parent,child,old,new,original,resumed=self.meshes(tmp)
   for relative in (Path("meshes/l320_h5/mesh.npz"),a.ARCHIVED_DC_MESH_ROOT/"other/mesh.npz"):
    target=child/relative;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(old.read_bytes())
    changed=copy.deepcopy(resumed);changed["source_path"]=str(target)
    with self.subTest(relative=relative),patch.object(a.guard,"SCRATCH",base),self.assertRaises(ValueError):
     a.geometry_relocation(original,changed,parent,child)
 def test_archived_fields_values_types_and_namespace_remain_exact(self):
  import tempfile
  with tempfile.TemporaryDirectory(prefix="archived_geometry_") as tmp:
   base,parent,child,old,new,original,resumed=self.meshes(tmp)
   for key,value in (("n_nodes",3659.),("n_nodes",3660),("length_m",3.3e-7),("new_mesh_generated",0),("extra",1)):
    changed=copy.deepcopy(resumed);changed[key]=value
    with self.subTest(key=key,value=value),patch.object(a.guard,"SCRATCH",base),self.assertRaises(ValueError):
     a.geometry_relocation(original,changed,parent,child)
 def test_symlink_file_and_directory_rejected_inside_and_outside_roots(self):
  import tempfile
  for side in ("parent","child"):
   for kind in ("file_inside","file_outside","directory_inside","directory_outside"):
    with self.subTest(side=side,kind=kind),tempfile.TemporaryDirectory(prefix="archived_geometry_") as tmp:
     base,parent,child,old,new,original,resumed=self.meshes(tmp)
     target_path=old if side=="parent" else new
     if kind.startswith("file"):
      content=target_path.read_bytes();target_path.unlink()
      target=(target_path.parent/"real.npz") if kind.endswith("inside") else (base/"outside.npz")
      target.write_bytes(content);target_path.symlink_to(target)
     else:
      folder=target_path.parent;content=target_path.read_bytes();target_path.unlink();folder.rmdir()
      target=(folder.parent/"real_folder") if kind.endswith("inside") else (base/"outside_folder")
      target.mkdir();(target/"mesh.npz").write_bytes(content);folder.symlink_to(target,target_is_directory=True)
     with patch.object(a.guard,"SCRATCH",base),self.assertRaises(ValueError):
      a.geometry_relocation(original,resumed,parent,child)
 def test_generated_mesh_symlinks_are_also_rejected(self):
  import tempfile
  with tempfile.TemporaryDirectory(prefix="archived_geometry_") as tmp:
   base,parent,child,old,new,original,resumed=self.meshes(tmp,"meshes/l320_h5/mesh.npz")
   saved=new.read_bytes();new.unlink();target=new.parent/"real.npz";target.write_bytes(saved);new.symlink_to(target)
   with patch.object(a.guard,"SCRATCH",base),self.assertRaises(ValueError):a.geometry_relocation(original,resumed,parent,child)
 def test_parent_traversal_and_symlink_roots_are_rejected(self):
  import tempfile
  with tempfile.TemporaryDirectory(prefix="archived_geometry_") as tmp:
   base,parent,child,old,new,original,resumed=self.meshes(tmp)
   changed=copy.deepcopy(original);changed["source_path"]=str(old.parent/".."/old.parent.name/old.name)
   with patch.object(a.guard,"SCRATCH",base),self.assertRaises(ValueError):a.geometry_relocation(changed,resumed,parent,child)
   alias=base/"parent_alias";alias.symlink_to(parent,target_is_directory=True)
   changed["source_path"]=str(alias/old.relative_to(parent))
   with patch.object(a.guard,"SCRATCH",base),self.assertRaises(ValueError):a.geometry_relocation(changed,resumed,alias,child)
 def test_existing_runtime_copytree_includes_archived_mesh_and_manifest_binding(self):
  import tempfile,shutil
  with tempfile.TemporaryDirectory(prefix="archived_geometry_") as tmp:
   base,parent,child,old,new,original,resumed=self.meshes(tmp)
   new.unlink();new.parent.rmdir();new.parent.parent.rmdir()
   shutil.rmtree(child)
   rel=old.relative_to(parent/"runtime_sources")
   archive=dict(package_sources_and_small_geometry={str(rel):dict(sha256=original["mesh_sha256"],bytes=old.stat().st_size)})
   (parent/"runtime_sources/archive_manifest.json").write_text(json.dumps(archive))
   child.mkdir();shutil.copytree(parent/"runtime_sources",child/"runtime_sources")
   self.assertEqual((child/"runtime_sources/archive_manifest.json").read_bytes(),(parent/"runtime_sources/archive_manifest.json").read_bytes())
   self.assertEqual(new.read_bytes(),old.read_bytes())
   with patch.object(a.guard,"SCRATCH",base):self.assertTrue(a.geometry_relocation(original,resumed,parent,child))


class HistoricalResourceReceiptOrigin(unittest.TestCase):
 def test_current_receipt_origin_still_matches_exact_current_source(self):
  receipt=dict(helper_sha256=hashlib.sha256(Path(c.__file__).read_bytes()).hexdigest())
  self.assertTrue(c.verify_receipt_helper_origin("spectral_2to4",Path("/not_a_historical_receipt"),receipt))
 def test_exact_real_domain_receipt_and_versionable_source_are_authenticated(self):
  self.assertEqual(hashlib.sha256(c.HISTORICAL_DOMAIN_HELPER_SOURCE.read_bytes()).hexdigest(),c.HISTORICAL_DOMAIN_HELPER_SHA)
  self.assertEqual(c.HISTORICAL_DOMAIN_HELPER_SOURCE.parent.name,"resource_transition_history")
  receipt=json.loads(c.HISTORICAL_DOMAIN_RECEIPT.read_text())
  self.assertEqual(hashlib.sha256(c.HISTORICAL_DOMAIN_RECEIPT.read_bytes()).hexdigest(),c.HISTORICAL_DOMAIN_RECEIPT_SHA)
  self.assertTrue(c.verify_receipt_helper_origin("domain_2to4",c.HISTORICAL_DOMAIN_RECEIPT,receipt))
 def test_historical_origin_is_not_a_generic_helper_allowlist(self):
  receipt=json.loads(c.HISTORICAL_DOMAIN_RECEIPT.read_text())
  for name,path in (("spectral_2to4",c.HISTORICAL_DOMAIN_RECEIPT),("spectral_4to8",c.HISTORICAL_DOMAIN_RECEIPT),
   ("domain_2to4",c.HISTORICAL_DOMAIN_RECEIPT.parent/"other_receipt.json")):
   with self.subTest(name=name,path=path),self.assertRaises(ValueError):c.verify_receipt_helper_origin(name,path,receipt)
  changed=copy.deepcopy(receipt);changed["helper_sha256"]="0"*64
  with self.assertRaises(ValueError):c.verify_receipt_helper_origin("domain_2to4",c.HISTORICAL_DOMAIN_RECEIPT,changed)
 def test_historical_receipt_or_source_tampering_is_rejected(self):
  import tempfile
  receipt=json.loads(c.HISTORICAL_DOMAIN_RECEIPT.read_text())
  with tempfile.TemporaryDirectory(prefix="resource_history_") as tmp:
   source=Path(tmp)/"historical_source.py";source.write_bytes(c.HISTORICAL_DOMAIN_HELPER_SOURCE.read_bytes())
   rp=Path(tmp)/"transfer_receipt.json";rp.write_bytes(c.HISTORICAL_DOMAIN_RECEIPT.read_bytes())
   with patch.object(c,"HISTORICAL_DOMAIN_RECEIPT",rp),patch.object(c,"HISTORICAL_DOMAIN_HELPER_SOURCE",source):
    self.assertTrue(c.verify_receipt_helper_origin("domain_2to4",rp,receipt))
    rp.write_bytes(rp.read_bytes()+b" ")
    with self.assertRaises(ValueError):c.verify_receipt_helper_origin("domain_2to4",rp,receipt)
   # Preserve original bytes independently because the patched path is the fixture.
   rp.write_bytes(json.dumps(receipt).encode())
   with patch.object(c,"HISTORICAL_DOMAIN_RECEIPT",rp),patch.object(c,"HISTORICAL_DOMAIN_RECEIPT_SHA",hashlib.sha256(rp.read_bytes()).hexdigest()),patch.object(c,"HISTORICAL_DOMAIN_HELPER_SOURCE",source):
    source.write_bytes(source.read_bytes()+b"\n")
    with self.assertRaises(ValueError):c.verify_receipt_helper_origin("domain_2to4",rp,receipt)
 def test_fine_origin_helper_stays_exact_and_outside_domain_exception(self):
  fine=Path("/home/jdiaz/scratch/stage5_effective_L320_h35_775_adaptive_v6_guard2b_resources4_20261011/provenance/transfer_receipt.json")
  receipt=json.loads(fine.read_text())
  self.assertEqual(hashlib.sha256(fine.read_bytes()).hexdigest(),"fd8289953ef1d6a8bdeacbd60c3ee8150b541d6677d43cfd4768752cadf6dc54")
  self.assertEqual(receipt["helper_sha256"],a.HELPER_SHA)
  self.assertEqual(hashlib.sha256(Path(r.__file__).read_bytes()).hexdigest(),a.HELPER_SHA)
  with self.assertRaises(ValueError):c.verify_receipt_helper_origin("domain_2to4",fine,receipt)
 def test_all_other_control_verification_statements_remain_ast_exact(self):
  import ast
  old=ast.parse(c.HISTORICAL_DOMAIN_HELPER_SOURCE.read_text())
  new=ast.parse(Path(c.__file__).read_text())
  old_fn=next(n for n in old.body if isinstance(n,ast.FunctionDef) and n.name=="verify_resource_run")
  new_fn=next(n for n in new.body if isinstance(n,ast.FunctionDef) and n.name=="verify_resource_run")
  # Original origin gate is replaced by the same case/action gate plus a closed provenance call.
  self.assertEqual(ast.dump(ast.Module(body=old_fn.body[6:],type_ignores=[])),
   ast.dump(ast.Module(body=new_fn.body[7:],type_ignores=[])))
  call=new_fn.body[6]
  self.assertIsInstance(call,ast.Expr);self.assertEqual(call.value.func.id,"verify_receipt_helper_origin")

if __name__=="__main__":unittest.main()
