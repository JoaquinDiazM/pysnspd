"""Strict additive reconstruction of explicit-transfer numerical ancestry.
Reads immutable receipts/anchors only. Never signals jobs or rewrites history.
"""
from __future__ import annotations
import argparse,copy,hashlib,json,math
from pathlib import Path
import numpy as np
SCRATCH=Path("/home/jdiaz/scratch").resolve()
KEYS=("dt_ps","step_tolerance","threads","local_energy_tolerance")
STATE=("delta_bar","p","n","circuit_SI")

def digest(path):
 h=hashlib.sha256()
 with Path(path).open("rb") as f:
  for block in iter(lambda:f.read(8*1024*1024),b""):h.update(block)
 return h.hexdigest()
def read(path):return json.loads(Path(path).read_text())
def close(a,b,label):
 if a!=b and abs(float(a)-float(b))>64*max(math.ulp(float(a)),math.ulp(float(b))):
  raise ValueError("Lineage clock mismatch: "+label)
def under(path,root=SCRATCH):
 p=Path(path).resolve()
 if not p.is_relative_to(root) or p==root:raise ValueError("Lineage path escapes scratch")
 return p
def flags(argv):
 ans={}
 for i,x in enumerate(argv):
  if x.startswith("--") and i+1<len(argv) and not argv[i+1].startswith("--"):
   ans[x[2:].replace("-","_")]=argv[i+1]
 return ans
def segments_valid(segments):
 if not segments:raise ValueError("Empty numerical lineage")
 close(segments[0]["start_time_ps"],0.,"initial clock")
 for i,s in enumerate(segments):
  if s["start_time_ps"]<0:raise ValueError("Negative segment clock")
  if i:close(segments[i-1]["end_time_ps"],s["start_time_ps"],"segment continuity")
  if s["end_time_ps"] is None:
   if i!=len(segments)-1:raise ValueError("Only final segment may remain open")
  elif s["end_time_ps"]<s["start_time_ps"]:raise ValueError("Reversed numerical segment")
  for k in KEYS:
   if k=="local_energy_tolerance" and s.get(k) is None:continue
   v=s[k]
   if not math.isfinite(float(v)) or float(v)<=0:raise ValueError("Invalid segment parameter "+k)
 if segments[-1]["end_time_ps"] is not None:raise ValueError("Final continuation must remain open")

def prior_receipt(run,manifest,anchor=None):
 """A recomposed resume manifest may omit its descriptive receipt pointer."""
 run=under(run)
 candidates=[]
 declared=manifest.get("explicit_transfer_receipt")
 if declared:candidates.append(under(run/declared,run))
 if anchor is not None:
  saved=under(anchor)/"provenance/parent_transfer_receipt.json"
  if saved.exists():candidates.append(saved)
 fallback=run/"provenance/transfer_receipt.json"
 if fallback.exists():candidates.append(fallback)
 if not candidates:
  if declared:raise ValueError("Declared prior transfer receipt is missing")
  return None
 hashes=[]
 for p in candidates:
  if not p.is_file():raise ValueError("Declared prior transfer receipt is missing")
  hashes.append(digest(p))
 if len(set(hashes))!=1:raise ValueError("Preserved and original transfer receipts disagree")
 return candidates[0]

def verify_receipt(path,expected_after_pins=None):
 path=under(path);r=read(path)
 if r.get("schema") not in ("pysnspd.explicit_checkpoint_transfer.v1","pysnspd.explicit_checkpoint_transfer.v2","pysnspd.explicit_checkpoint_transfer.v3"):
  raise ValueError("Unrecognized transfer receipt schema")
 if r["state_arrays_sha256_before"]!=r["state_arrays_sha256_after"]:
  raise ValueError("Transfer changed physical state arrays")
 if r["conserved_metadata_before"]!=r["conserved_metadata_after"]:
  raise ValueError("Transfer changed conserved physical metadata")
 if expected_after_pins is not None and r["pins_after"]!=expected_after_pins:
  raise ValueError("Ancestor pins_after differ from descendant pins_before")
 if len(r["pins_before"])!=13 or set(r["pins_before"])!=set(r["pins_after"]):
  raise ValueError("Transfer must retain the 13-pin namespace")
 anchor=under(r["parent"]["anchor"])
 cp=anchor/"checkpoint_original.npz";manifest=anchor/"manifest_original.json";capture=anchor/"capture_receipt.json"
 if digest(cp)!=r["parent"]["checkpoint_bytes_sha256"]:raise ValueError("Ancestor checkpoint hash differs")
 if digest(manifest)!=r["parent"]["manifest_sha256"]:raise ValueError("Ancestor manifest hash differs")
 cr=read(capture);pm=read(manifest)
 if cr["checkpoint_bytes_sha256"]!=r["parent"]["checkpoint_bytes_sha256"] or cr["manifest_sha256"]!=r["parent"]["manifest_sha256"]:
  raise ValueError("Capture and transfer hashes differ")
 if cr["source_pins_sha256"]!=r["pins_before"] or pm["source_pins_sha256"]!=r["pins_before"]:
  raise ValueError("Captured source identity differs")
 if cr["conserved_metadata"]!=r["conserved_metadata_before"] or cr["state_arrays_sha256"]!=r["state_arrays_sha256_before"]:
  raise ValueError("Captured state/ledger differs")
 with np.load(cp,allow_pickle=False) as z:meta=json.loads(str(z["metadata_json"].item()))
 if meta["parameter_identity"]!=cr["parameter_identity"] or meta["source_pins_sha256"]!=r["pins_before"]:
  raise ValueError("Checkpoint numerical/source identity differs")
 for k,v in r["conserved_metadata_before"].items():
  if meta[k]!=v:raise ValueError("Checkpoint conserved metadata differs: "+k)
 archive=anchor/"runtime_sources/archive_manifest.json"
 if digest(archive)!=r["parent"]["runtime_archive_manifest_sha256"]:
  raise ValueError("Ancestor source archive manifest hash differs")
 ar=read(archive)
 if ar["runtime_pins_sha256"]!=r["pins_before"]:raise ValueError("Ancestor archive source pins differ")
 for rel,h in r["pins_before"].items():
  if digest(anchor/"runtime_sources"/rel)!=h:raise ValueError("Ancestor frozen source hash differs: "+rel)
 seg=copy.deepcopy(r["numerical_segments"]);segments_valid(seg)
 now=float(meta["time_ps"]);close(seg[-1]["start_time_ps"],now,"transfer activation")
 argv=flags(r["resume_argv"])
 for k in KEYS:
  actual=seg[-1].get(k)
  given=argv.get(k)
  if actual is None:
   if given is not None:raise ValueError("Undeclared continuation guard")
  elif given is None or float(actual)!=float(given):raise ValueError("Resume argv/segment differs: "+k)
 proof=dict(path=str(path),sha256=digest(path),parent_anchor=str(anchor),
  checkpoint_sha256=digest(cp),manifest_sha256=digest(manifest),
  capture_receipt_sha256=digest(capture),runtime_archive_sha256=digest(archive),
  pins_before=r["pins_before"],pins_after=r["pins_after"],
  arrays_and_conserved_metadata_equal=True,activation_time_ps=now,
  receipt_hash_references=[])
 # These immutable launch/stop intents can independently bind a receipt hash.
 origin=under(r["output_dir"])/"provenance"
 if origin.exists():
  for q in sorted(origin.glob("*intent*.json")):
   if proof["sha256"] in q.read_text():
    proof["receipt_hash_references"].append(dict(path=str(q),sha256=digest(q)))
 return r,meta,pm,proof

def resolve_receipt(path,expected_after_pins=None,seen=None):
 seen=set() if seen is None else set(seen)
 path=under(path)
 if path in seen:raise ValueError("Cycle in transfer ancestry")
 seen.add(path)
 r,meta,pm,proof=verify_receipt(path,expected_after_pins)
 anchor=under(r["parent"]["anchor"]);parent_run=under(r["parent"]["output_dir"])
 prior=prior_receipt(parent_run,pm,anchor)
 if prior is None:
  seg=copy.deepcopy(r["numerical_segments"]);proof["ancestry_kind"]="original run or explicitly documented diagnostic"
  # Diagnostic overrides are proven by their preserved source receipt, not by inference.
  for s in seg[:-1]:
   if s.get("origin")=="frozen_diagnostic_A":
    p=anchor/"parent_preparation_receipt.json"
    if digest(p)!=s["parent_preparation_receipt_sha256"]:raise ValueError("Diagnostic provenance hash differs")
    d=read(p)
    if d["checkpoint_sha256"]!=proof["checkpoint_sha256"]:raise ValueError("Diagnostic checkpoint binding differs")
  proofs=[proof]
 else:
  seg,proofs=resolve_receipt(prior,r["pins_before"],seen)
  for k in KEYS:
   if seg[-1].get(k)!=meta["parameter_identity"].get(k):
    raise ValueError("Ancestor final numerical segment differs from captured parent: "+k)
  close(seg[-1]["start_time_ps"],float(seg[-1]["start_time_ps"]),"ancestor segment")
  if seg[-1]["start_time_ps"]>meta["time_ps"]:raise ValueError("Ancestor starts after captured child")
  seg[-1]["end_time_ps"]=float(meta["time_ps"])
  seg.append(copy.deepcopy(r["numerical_segments"][-1]))
  proof["ancestry_kind"]="explicit ancestor receipt recovered without rewriting descendant"
  proof["ancestor_receipt_sha256"]=proofs[-1]["sha256"]
  proofs.append(proof)
 for s in seg:
  s.setdefault("local_energy_tolerance",None)
  s.setdefault("energy_guard_enabled",s["local_energy_tolerance"] is not None)
 segments_valid(seg)
 return seg,proofs

def supplement(run):
 run=under(run)
 path=run/"provenance/transfer_receipt.json"
 original_hash=digest(path)
 try:
  segments,proofs=resolve_receipt(path)
  return dict(schema="pysnspd.numerical_lineage_supplement.v1",run_dir=str(run),
   original_transfer_receipt_path=str(path),original_transfer_receipt_sha256=original_hash,
   numerical_lineage_complete=True,numerical_segments=segments,ancestor_proofs=proofs,
   originals_modified=False,physical_arrays_modified=False,runner_modified=False,
   explanation="Additive recovery from original receipts and SHA-verified frozen anchors; no retrospective numerical policy inferred.")
 except (ValueError,KeyError,FileNotFoundError,TypeError) as e:
  return dict(schema="pysnspd.numerical_lineage_supplement.v1",run_dir=str(run),
   original_transfer_receipt_path=str(path),original_transfer_receipt_sha256=original_hash,
   numerical_lineage_complete=False,numerical_segments=None,reason=str(e),
   originals_modified=False,physical_arrays_modified=False,runner_modified=False)
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument("--run",required=True);p.add_argument("--output",required=True)
 a=p.parse_args();out=under(a.output)
 if out.exists():raise FileExistsError("Never overwrite a lineage supplement")
 r=supplement(a.run);out.parent.mkdir(parents=True,exist_ok=True)
 out.write_text(json.dumps(r,indent=2,allow_nan=False)+"\n")
 print(json.dumps(dict(path=str(out),sha256=digest(out),complete=r["numerical_lineage_complete"])))
 if not r["numerical_lineage_complete"]:raise SystemExit(1)
if __name__=="__main__":main()
