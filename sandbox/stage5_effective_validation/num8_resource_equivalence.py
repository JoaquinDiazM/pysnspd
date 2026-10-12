"""Separate Num8 oracle/candidate certificate; never runs RHS or overwrites Num124."""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from sandbox.stage5_effective_validation import fine_resource_transition as f
DRIVER=Path("/home/jdiaz/scratch/stage5_node_fusion_candidate_20261011/verify_candidate.py")
DRIVER_SHA="b529544b287f5bd2c7c453aa8bc9a1d9d305668f0f5407f9597419182c82c76e"
BASE=DRIVER.parent
SCHEMA="pysnspd.num8_oracle_equivalence.v1"
RUNTIME=Path("/home/jdiaz/scratch/stage5_effective_runtime_v6_20261011")
LEGACY_FRONTENDS={"pysnspd/experimental/effective_bcs_preparation.py":("36a6d56d2e9fa520057c4c5511f8c890ef06d7e3ff8ad883c35f49c5b3d5a52b","a3c614baff8a6bdd80e89639e0b3ad53f963346373c6bc4da7a33e014c06f41c"),"pipelines/05_effective_photon_run.py":("234629c32b71bc39e96e20454e1c80cea565fe196e1067d4980f9e88f20bad73","70b7c5b5d39c74f19a35ad816f413ecc47451c46651817708549b8e9fc8e97aa")}
def require_reviewed_exception(oracle,candidate,reviewed):
 if oracle["error_class"]!=candidate["error_class"] or oracle["error"]!=candidate["error"] or oracle["error_class"]!=reviewed["error_class"] or oracle["error"]!=reviewed["error"]:raise ValueError("Num8 support exception differs from the original reviewedNum4 case")
 return True
def comparison():
 if f.digest(DRIVER)!=DRIVER_SHA or f.digest(f.CERTIFICATE)!=f.CERTIFICATE_SHA:raise ValueError("Existing oracle driver/certificate changed")
 original=f.read(f.CERTIFICATE)
 entries={role:f.entry(f.owned(BASE/"results"/(role+"_n8")/"receipt.json")) for role in ("oracle","candidate")}
 reports={role:f.read(f.bound(e)) for role,e in entries.items()}
 aa,bb=reports["oracle"],reports["candidate"]
 for role,r in reports.items():
  if r.get("status")!="COMPLETED" or r.get("mode")!=role or r.get("Numba_threads")!=8 or r.get("CPU_affinity")!=list(range(6,14)) or r.get("BLAS_threads")!=1 or r.get("driver_sha256")!=DRIVER_SHA or not 0<r.get("wall_s",301)<=300:raise ValueError("Each real Num8 oracle/candidate run must be complete, brief and use free6..13")
  source=r["source_sha256"];expected=original["sources_before" if role=="oracle" else "sources_after"]
  if any(source.get(k)!=v for k,v in expected.items()) or len(source)!=13:raise ValueError("Old/new frozen source identity differs")
 if {k:v for k,v in aa["source_sha256"].items() if k not in original["sources_before"]}!={k:v for k,v in bb["source_sha256"].items() if k not in original["sources_after"]}:raise ValueError("Other frozen physical sources differ")
 cases={r["case"]:r for r in aa["cases"]};other={r["case"]:r for r in bb["cases"]}
 reviewed={r["case"]:r for r in original["cases"] if r["thread_count"]==4}
 required=set(reviewed)
 if set(cases)!=set(other) or set(cases)!=required:raise ValueError("The complete previously reviewed case set is required")
 pins=f.read(RUNTIME/"archive_manifest.json")["runtime_pins_sha256"]
 if len(pins)!=13 or any(f.digest(RUNTIME/k)!=v for k,v in pins.items()):raise ValueError("Current runtime archive13sources changed")
 differences={k:(bb["source_sha256"].get(k),v) for k,v in pins.items() if bb["source_sha256"].get(k)!=v}
 if differences!=LEGACY_FRONTENDS:raise ValueError("Only the exact2historical frontend differences are allowed in the oracle-driver scope")
 evidence=[];cross=[]
 for name,a in cases.items():
  b=other[name]
  if a["status"]!=b["status"] or a["status"]!=reviewed[name]["status"] or a["input"]!=b["input"] or a["input"]!=reviewed[name]["input_identity"]:raise ValueError("Paired or original reviewed input/status differs")
  if a["status"]=="EXCEPTION":
   require_reviewed_exception(a,b,reviewed[name])
   evidence.append(dict(case=name,status="EXCEPTION",exception_class=a["error_class"],error=a["error"]));continue
  if a["status"]!="RESULT":raise ValueError("Unfinished comparison case")
  ap=f.owned(BASE/"results/oracle_n8"/(name+".npz"));bp=f.owned(BASE/"results/candidate_n8"/(name+".npz"))
  if f.digest(ap)!=a["output_sha256"] or f.digest(bp)!=b["output_sha256"]:raise ValueError("Recorded output SHA differs")
  with np.load(ap,allow_pickle=False) as x,np.load(bp,allow_pickle=False) as y:
   if set(x.files)!=set(y.files) or not x.files or any(not f.same_bytes(x[k],y[k]) for k in x.files):raise ValueError("WithinNum8 old/new full arrays or diagnostics differ")
   # Legacy reductions can differ across thread counts; report every such difference.
   cp=f.owned(BASE/"results/oracle_n4"/(name+".npz"))
   if str(cp)!=reviewed[name]["oracle_path"] or f.digest(cp)!=reviewed[name]["oracle_npz_sha256"]:raise ValueError("ExistingNum4 oracle bytes differ from original certificate")
   with np.load(cp,allow_pickle=False) as z:
    if set(z.files)!=set(x.files):raise ValueError("ExistingNum4 oracle output namespace differs")
    diffs=[dict(field=k,maximum_absolute_difference=float(np.max(abs(x[k].astype(np.complex128)-z[k].astype(np.complex128)))),num4_scale=float(np.max(abs(z[k])))) for k in x.files if not f.same_bytes(x[k],z[k])]
   cross.append(dict(case=name,num4_oracle=f.entry(cp),original_reduction_differences=diffs))
  evidence.append(dict(case=name,status="RESULT",oracle=f.entry(ap),candidate=f.entry(bp),all_fields_byte_equal=True))
 return dict(schema=SCHEMA,status="COMPARISON_PASSED_NOT_RESOURCE_ADMISSION",original_num124_certificate=f.entry(f.CERTIFICATE),
  existing_driver=f.entry(DRIVER),num8_receipts=entries,source_pins_sha256=bb["source_sha256"],within_thread_counts=[8],
  oracle_candidate_cases=evidence,original_num4_to8_diagnostics=cross,current_runtime_pins_sha256=pins,
  current_runtime_archive=f.entry(RUNTIME/"archive_manifest.json"),legacy_driver_frontend_mapping={k:dict(driver_sha256=v[0],current_sha256=v[1]) for k,v in LEGACY_FRONTENDS.items()},
  scope="Frozen implementation/support/midpoint cases only, including synthetic spectral80; not a current fullL480/80-state trajectory or future cross-thread bitwise identity",
  physical_evolution_of_production_jobs=False,original_num124_certificate_changed=False)
def verify(path,pins,approval_path):
 path=f.owned(path);saved=f.read(path);actual=comparison()
 if saved!=actual or saved["current_runtime_pins_sha256"]!=pins:raise ValueError("Num8 certificate changed or runtime pins differ")
 approval=f.read(f.owned(approval_path))
 if approval.get("status")!="ROOT_APPROVED_NUM8_RESOURCE_EQUIVALENCE" or approval.get("certificate_sha256")!=f.digest(path) or approval.get("source_pins_sha256")!=pins or approval.get("stage5") is not False:raise ValueError("Separate root approval of the real Num8 SHA is required")
 return dict(certificate=f.entry(path),root_approval=f.entry(approval_path),scope=saved["scope"])
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument("--output",type=Path,required=True);a=p.parse_args()
 output=a.output.resolve()
 if not output.is_relative_to(f.SCRATCH) or output.exists() or output==f.CERTIFICATE:raise ValueError("Fresh separate scratch certificate required")
 original=f.digest(f.CERTIFICATE);r=comparison();output.write_text(json.dumps(r,indent=2,allow_nan=False)+"\n")
 if f.digest(f.CERTIFICATE)!=original:raise ValueError("OriginalNum124 certificate changed")
 print(json.dumps(dict(path=str(output),sha256=f.digest(output),resource_admission=False)))
if __name__=="__main__":main()
