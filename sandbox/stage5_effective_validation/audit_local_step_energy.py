"""Read-only local midpoint energy-error audit from coherent checkpoints."""
import argparse,hashlib,io,json
from pathlib import Path
import numpy as np
SCRATCH=Path("/home/jdiaz/scratch")
CASES={
 "A_strict_002":SCRATCH/"stage5_frozen_solver_probe_20261011/A_continuations/hmax_0p02/checkpoint_latest.npz",
 "A_practical_002":SCRATCH/"stage5_frozen_solver_probe_20261011/A_practical_norm_continuations/hmax_0p02/checkpoint_latest.npz",
 "live775_002":SCRATCH/"stage5_effective_L320_h5_775_dt002_v3_20261011/checkpoint_latest.npz",
 "live775_001":SCRATCH/"stage5_effective_L320_h5_775_dt001_v3_20261011/checkpoint_latest.npz",
 "live1550_002":SCRATCH/"stage5_effective_L320_h5_1550_dt002_v3_20261011/checkpoint_latest.npz"}
def audit(path,limit):
 raw=path.read_bytes()
 with np.load(io.BytesIO(raw),allow_pickle=False) as z:
  meta=json.loads(str(z["metadata_json"].item()))
  t=z["trajectory__time_s"].copy()*1e12
  energy=z["trajectory__energy_total_J"].copy()
  external=z["trajectory__external_energy_integral_J"].copy()
  defect=z["trajectory__energy_defect_J"].copy()
  hs=z["trajectory__accepted_step_ps"].copy()[1:]
 Eret=float(meta["source"]["retained_energy_J"]);work=np.diff(external)
 residual=np.diff(energy)-work;variation=np.r_[0.,np.cumsum(abs(work))]
 denominator=Eret+variation[1:];ratio=abs(residual)/denominator
 ix=int(np.argmax(ratio));windows={}
 for lo,hi in [(0,10),(10.4,12.4),(12.4,14),(14,15),(15,20),(20,200)]:
  select=(t[1:]>=lo-1e-8)&(t[1:]<=hi+1e-8)
  if not select.any():continue
  jj=np.flatnonzero(select);j=int(jj[np.argmax(ratio[select])])
  windows[f"{lo:g}_{hi:g}"]=dict(count=int(select.sum()),maximum_local_fraction=float(ratio[j]),
   maximum_time_ps=float(t[j+1]),maximum_residual_J=float(abs(residual[j])),
   steps_exceeding_proposal=int(np.sum(ratio[select]>limit)),global_defect_fraction_end=float(abs(defect[jj[-1]+1])/(Eret+variation[jj[-1]+1])))
 return dict(path=str(path),checkpoint_sha256=hashlib.sha256(raw).hexdigest(),final_time_ps=meta["time_ps"],
  retained_energy_J=Eret,total_absolute_external_work_J=float(variation[-1]),energy_denominator_final_J=float(Eret+variation[-1]),
  global_energy_defect_fraction=float(abs(defect[-1])/(Eret+variation[-1])),
  maximum_global_energy_defect_fraction=float(np.max(abs(defect)/(Eret+variation))),
  maximum_local_fraction=float(ratio[ix]),maximum_local_residual_J=float(abs(residual[ix])),
  maximum_local_time_ps=float(t[ix+1]),maximum_local_step_ps=float(hs[ix]),
  proposed_local_fraction_limit=limit,steps_exceeding_proposal=int(np.sum(ratio>limit)),
  local_fraction_quantiles={str(q):float(np.quantile(ratio,q)) for q in [.5,.9,.99,1.]},
  sum_absolute_local_residual_J=float(np.sum(abs(residual))),
  sum_absolute_local_residual_fraction=float(np.sum(abs(residual))/(Eret+variation[-1])),
  maximum_difference_from_incremental_global_defect_J=float(np.max(abs(residual-np.diff(defect)))),
  estimated_step_for_max_event_ps=float(hs[ix]*min(1.,(limit/ratio[ix])**(1/3))) if ratio[ix] else float(hs[ix]),
  step_estimate_is_asymptotic_not_admission=True,windows=windows)
def main():
 p=argparse.ArgumentParser();p.add_argument("--output",type=Path,required=True);p.add_argument("--fraction",type=float,default=1e-5)
 a=p.parse_args();result=dict(schema="pysnspd.local_midpoint_energy_error_audit.v1",proposal_fraction=a.fraction,
  proposal_equation="abs(Enew-Eold-midpoint_external_increment)<=fraction*(Eret+sum(abs(accepted_external_increments))+abs(trial_external_increment))",
  global_acceptance_unchanged_fraction=.01,read_only_no_evolution=True,cases={k:audit(v,a.fraction) for k,v in CASES.items()})
 a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+"\n")
 print(json.dumps(result,indent=2))
if __name__=="__main__":main()
