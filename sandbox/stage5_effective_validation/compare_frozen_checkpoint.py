"""Physical comparison of numerical-only continuation controls."""
from pathlib import Path
import argparse,hashlib,importlib.util,json
import numpy as np
def main():
 ap=argparse.ArgumentParser();ap.add_argument("root",type=Path);ap.add_argument("--partial",action="store_true");a=ap.parse_args();root=a.root
 spec=importlib.util.spec_from_file_location("checkpoint_probe",root/"probe_frozen_checkpoint_v3.py")
 probe=importlib.util.module_from_spec(spec);spec.loader.exec_module(probe)
 pipeline=probe.load_pipeline(root);system,_,initial,_,_,_,dc,geometry=probe.load_state(pipeline,root/"frozen_A")
 start=float(initial["time_ps"])*1e-12;Eret=system.kinetic_residual_energy_scale_J
 triangles=geometry.triangles
 xy=geometry.coordinates_m[triangles]
 orient=np.sign((xy[:,1,0]-xy[:,0,0])*(xy[:,2,1]-xy[:,0,1])-(xy[:,1,1]-xy[:,0,1])*(xy[:,2,0]-xy[:,0,0]))
 def vortex_counts(delta):
  values=delta[triangles];winding=np.rint(np.sum(np.angle(np.roll(values,-1,axis=1)*np.conj(values)),axis=1)/(2*np.pi))*orient
  return [int(np.sum(winding==1)),int(np.sum(winding==-1))]
 initial_sha=hashlib.sha256((root/"frozen_A/checkpoint_latest.npz").read_bytes()).hexdigest()
 def load(directory):
  with np.load(directory/"checkpoint_latest.npz",allow_pickle=False) as z:
   state=pipeline.EffectiveState(*(z[k].copy() for k in ("delta_bar","p","n","circuit_SI")))
   meta=json.loads(str(z["metadata_json"].item()))
   hist={k:z["trajectory__"+k].copy() for k in pipeline.HISTORY_KEYS}
   snaps={k:z["snapshots__"+k].copy() for k in pipeline.SNAPSHOT_KEYS}
  receipt=json.loads((directory/"receipt.json").read_text())
  if receipt["initial_checkpoint_sha256"]!=initial_sha:raise ValueError("Initial checkpoint differs")
  if receipt["source_pins_sha256"]!=initial["source_pins_sha256"]:raise ValueError("Physical source differs")
  if meta["source"]!=initial["source"]:raise ValueError("Photon source receipt differs")
  denominators=Eret+np.r_[0.,np.cumsum(abs(np.diff(hist["external_energy_integral_J"])))]
  denominator=float(denominators[-1])
  fraction=float(np.max(abs(hist["energy_defect_J"])/denominators))
  support=bool(state.p.min()>=0 and state.p.max()<=1 and state.n.min()>=0 and min(state.circuit_SI[:2])>0)
  return dict(directory=str(directory),state=state,metadata=meta,history=hist,snapshots=snaps,receipt=receipt,
   energy_denominator_J=denominator,maximum_energy_defect_fraction=fraction,endpoint_support_strict=support)
 def comparison(ref,test):
  stop=min(ref["history"]["time_s"][-1],test["history"]["time_s"][-1])
  times=np.unique(np.r_[ref["history"]["time_s"],test["history"]["time_s"]])
  times=times[(times>=start)&(times<=stop)]
  curves={}
  for key in ("Vout_V","Vdev_V","Is_A","domain_energy_J","energy_defect_J","limited_event_energy_J"):
   baseline=np.interp(times,ref["history"]["time_s"],ref["history"][key])
   candidate=np.interp(times,test["history"]["time_s"],test["history"][key])
   maximum=float(np.max(abs(candidate-baseline)));scale=float(np.max(abs(baseline)))
   curves[key]=dict(maximum_absolute_difference_SI=maximum,maximum_reference_signal_SI=scale,
    fraction_of_actual_reference_signal=maximum/scale if scale else None)
   if key in ("domain_energy_J","energy_defect_J","limited_event_energy_J"):
    curves[key]["fraction_of_retained_photon_energy"]=maximum/Eret
  rs,ts=ref["snapshots"],test["snapshots"]
  stimes=np.intersect1d(np.round(rs["time_s"]*1e12,7),np.round(ts["time_s"]*1e12,7))
  stimes=stimes[(stimes>=start*1e12-1e-7)&(stimes<=stop*1e12+1e-7)]
  gapdiff=gapfrac=gapmindiff=gapminfrac=map_l2dc=map_l2supp=Te=Tph=0.;vortices=[]
  roi=abs(geometry.coordinates_m[:,0]-.5*(geometry.coordinates_m[:,0].min()+geometry.coordinates_m[:,0].max()))<=60e-9
  volumes=system.functional.volumes_m3[roi]
  dcJ=abs(dc.delta_bar)*system.functional.Eu_J
  for t in stimes:
   i=np.argmin(abs(rs["time_s"]*1e12-t));j=np.argmin(abs(ts["time_s"]*1e12-t))
   diff=abs(abs(rs["delta_J"][i])-abs(ts["delta_J"][j]))
   vortices.append(dict(time_ps=float(t),reference_counts_plus_minus=vortex_counts(rs["delta_J"][i]),test_counts_plus_minus=vortex_counts(ts["delta_J"][j])))
   gapdiff=max(gapdiff,float(np.max(diff)))
   min_difference=float(abs(np.min(abs(rs["delta_J"][i]))-np.min(abs(ts["delta_J"][j]))))
   gapmindiff=max(gapmindiff,min_difference);gapminfrac=max(gapminfrac,min_difference/float(np.mean(dcJ)))
   numerator=float(np.sum(volumes*diff[roi]**2))
   map_l2dc=max(map_l2dc,float(np.sqrt(numerator/np.sum(volumes*dcJ[roi]**2))))
   suppression=float(np.sum(volumes*(dcJ[roi]-abs(rs["delta_J"][i][roi]))**2))
   if suppression>0:map_l2supp=max(map_l2supp,float(np.sqrt(numerator/suppression)))
   gapfrac=max(gapfrac,float(np.max(diff/(abs(dc.delta_bar)*system.functional.Eu_J))))
   Te=max(Te,float(np.max(abs(rs["Te_K"][i]-ts["Te_K"][j]))))
   Tph=max(Tph,float(np.max(abs(rs["Tph_K"][i]-ts["Tph_K"][j]))))
  common_endpoint=abs(float(ref["metadata"]["time_ps"])-float(test["metadata"]["time_ps"]))<1e-9
  endpoint=probe.compare_states(system,ref["state"],test["state"],.1) if common_endpoint else None
  signal=curves["Vout_V"]["fraction_of_actual_reference_signal"]
  passed=bool(common_endpoint and gapminfrac<=.02 and signal is not None and signal<=.02
    and all(curves[k]["fraction_of_actual_reference_signal"] is not None and curves[k]["fraction_of_actual_reference_signal"]<=.02 for k in ("Vdev_V","Is_A"))
    and ref["maximum_energy_defect_fraction"]<=.01 and test["maximum_energy_defect_fraction"]<=.01
    and ref["endpoint_support_strict"] and test["endpoint_support_strict"])
  return dict(reference=ref["directory"],test=test["directory"],common_window_ps=[start*1e12,stop*1e12],
   curves=curves,common_snapshot_count=len(stimes),maximum_map_gap_magnitude_difference_meV=gapdiff/1.602176634e-22,
   maximum_map_gap_difference_fraction_of_local_DC=gapfrac,
   maximum_gap_minimum_difference_meV=gapmindiff/1.602176634e-22,maximum_gap_minimum_difference_fraction_of_DC=gapminfrac,
   maximum_map_L2_difference_fraction_DC_ROI120nm=map_l2dc,maximum_map_L2_difference_fraction_actual_suppression_ROI120nm=map_l2supp,
   map_threshold_is_diagnostic_only=True,
   maximum_map_temperature_difference_K=Te,maximum_map_phonon_temperature_difference_K=Tph,
   common_endpoint=common_endpoint,endpoint=endpoint,population_force_response_step_ps=.1,
   vortex_counts=vortices,maximum_vortex_count_difference=int(max((max(abs(np.array(v["reference_counts_plus_minus"])-np.array(v["test_counts_plus_minus"]))) for v in vortices),default=0)),
   reference_energy_defect_fraction=ref["maximum_energy_defect_fraction"],test_energy_defect_fraction=test["maximum_energy_defect_fraction"],
   test_energy_denominator_J=test["energy_denominator_J"],reference_energy_denominator_J=ref["energy_denominator_J"],
   reference_elapsed_s=ref["receipt"]["elapsed_wall_s"],test_elapsed_s=test["receipt"]["elapsed_wall_s"],
   reference_new_steps=ref["receipt"]["new_steps"],test_new_steps=test["receipt"]["new_steps"],
   reference_new_rejections=ref["receipt"]["new_rejections"],test_new_rejections=test["receipt"]["new_rejections"],
   endpoint_support_strict=ref["endpoint_support_strict"] and test["endpoint_support_strict"],
   passes_predeclared_control_thresholds=passed,stage5_admission=False)
 strict=load(root/"A_continuations/hmax_0p02")
 paths=[root/"A_practical_norm_continuations"/("hmax_"+str(h).replace(".","p")) for h in (.02,.05,.1)]
 if not a.partial and any(not (p/"receipt.json").is_file() for p in paths):raise ValueError("Requested controls incomplete")
 practical=[load(p) for p in paths if (p/"receipt.json").is_file()]
 if not practical or Path(practical[0]["directory"]).name!="hmax_0p02":raise ValueError("Practical h.02 reference missing")
 result=dict(schema="pysnspd.frozen_numerical_control_comparison.v1",initial_checkpoint_sha256=initial_sha,
  source_pins_sha256=initial["source_pins_sha256"],retained_energy_J=Eret,complete_requested_caps=len(practical)==3,
  resource_segment_receipt=str(root/"affinity_segment_13_to_7.json"),cost_note="Strict baseline used contendedCPU13; practical case usedCPU13 thenCPU7; elapsed savings cannot be attributed entirely to tolerance",
  thresholds=dict(gap_minimum_fraction_DC=.02,signal_fraction_actual=.02,energy_fraction=.01,strict_support=True),
  tolerance_strict_vs_practical_2e7=comparison(strict,practical[0]),
  time_step_practical={Path(case["directory"]).name:comparison(practical[0],case) for case in practical[1:]},
  no_RHS_changes=True,no_photon_reinjection=True,no_live_run_mutations=True,no_admission=True)
 probe.writej(root/"A_practical_norm_continuations"/("comparison_partial.json" if a.partial else "comparison.json"),result)
 print(json.dumps(result,indent=2))
if __name__=="__main__":main()
