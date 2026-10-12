"""Warm cost profile of the frozen physical RHS and unchanged midpoint solver."""
import argparse,cProfile,hashlib,importlib.util,io,json,pstats,time,os
from pathlib import Path
import numpy as np
def profile_summary(prof):
 stats=pstats.Stats(prof);rows=[]
 for (filename,line,name),(primitive,total,selftime,cumtime,callers) in stats.stats.items():
  rows.append(dict(file=filename,line=line,function=name,primitive_calls=primitive,total_calls=total,
   self_time_s=selftime,cumulative_time_s=cumtime))
 return sorted(rows,key=lambda r:r["cumulative_time_s"],reverse=True)

def buffer_profile(node_counts,event_count):
 """Same five event-array constructors, no collision calculation or physical RHS."""
 from numba import njit
 @njit(cache=False)
 def packet(nn,nev):
  return (np.zeros((nn,nev)),np.zeros((nn,nev)),np.zeros((nn,nev)),
   np.zeros((nn,nev),np.int64),np.ones((nn,nev),np.int64))
 packet(1,1)
 rows=[]
 for nn in node_counts:
  samples=[]
  for repeat in range(3):
   start=time.monotonic();cpu=time.process_time();arrays=packet(nn,event_count)
   alloc_cpu=time.process_time()-cpu;alloc_wall=time.monotonic()-start
   start=time.monotonic();cpu=time.process_time()
   for index,arr in enumerate(arrays):arr.fill(1 if index==4 else 0)
   touch_cpu=time.process_time()-cpu;touch_wall=time.monotonic()-start
   samples.append(dict(constructor_cpu_s=alloc_cpu,constructor_wall_s=alloc_wall,
    supplemental_full_touch_cpu_s=touch_cpu,supplemental_full_touch_wall_s=touch_wall))
   del arrays
  rows.append(dict(nodes=nn,events=event_count,five_array_bytes=nn*event_count*5*8,samples=samples))
 return dict(rows=rows,
  interpretation='Isolated constructor packet is a measured buffer cost proxy, not an exact decomposition of the opaque acoustic_collisions dispatcher. Supplemental touch repeats writes and must not be added as actual RHS cost.')

def main():
 p=argparse.ArgumentParser();p.add_argument("root",type=Path);a=p.parse_args();root=a.root
 spec=importlib.util.spec_from_file_location("probe",root/"probe_frozen_checkpoint_v3.py")
 probe=importlib.util.module_from_spec(spec);spec.loader.exec_module(probe)
 pipeline=probe.load_pipeline(root);system,_,meta,_,_,_,_,_=probe.load_state(pipeline,root/"frozen_A")
 checkpoint=root/"A_practical_norm_continuations/hmax_0p1/checkpoint_latest.npz";raw=checkpoint.read_bytes()
 with np.load(io.BytesIO(raw),allow_pickle=False) as z:
  state=pipeline.EffectiveState(*(z[k].copy() for k in ("delta_bar","p","n","circuit_SI")))
  current=json.loads(str(z["metadata_json"].item()))
 out=root/"warm_profile";out.mkdir(exist_ok=True)
 system.evaluate(state,kinetic_budget=state,kinetic_dt_ps=.02)
 prof=cProfile.Profile(timer=time.process_time);start=time.monotonic();cpustart=time.process_time()
 try:
  candidate,receipt=prof.runcall(system.midpoint_step,state,.02,tolerance=2e-7,max_iterations=40)
  stepstatus="CONVERGED";error=None
 except (ArithmeticError,ValueError,RuntimeError) as exc:
  candidate=None;receipt=None;stepstatus="FAILED";error=repr(exc)
 stepwall=time.monotonic()-start;stepcpu=time.process_time()-cpustart;prof.dump_stats(str(out/"midpoint.prof"))
 with (out/"midpoint_top.txt").open("w") as s:pstats.Stats(prof,stream=s).sort_stats("cumtime").print_stats(60)
 if candidate is None:mid=state.copy()
 else:
  du=np.r_[np.column_stack(((candidate.delta_bar-state.delta_bar).real,(candidate.delta_bar-state.delta_bar).imag)).ravel()[system._free_cart],
    [np.angle(candidate.delta_bar[nodes[0]]/state.delta_bar[nodes[0]]) for nodes in system.cuts]]
  dp=candidate.p-state.p;dn=candidate.n-state.n
  dc=(candidate.circuit_SI-state.circuit_SI)/system.circuit_scales
  mid=system.shifted(state,du,dp,dn,dc,.5)
 system.evaluate(mid,kinetic_budget=state,kinetic_dt_ps=.02)
 evprof=cProfile.Profile(timer=time.process_time);times=[];cputimes=[]
 for _ in range(5):
  start=time.monotonic();cpustart=time.process_time();evprof.runcall(system.evaluate,mid,kinetic_budget=state,kinetic_dt_ps=.02);times.append(time.monotonic()-start);cputimes.append(time.process_time()-cpustart)
 evprof.dump_stats(str(out/"five_evaluate.prof"))
 with (out/"five_evaluate_top.txt").open("w") as s:pstats.Stats(evprof,stream=s).sort_stats("cumtime").print_stats(60)
 energies=np.hypot(system.kinetics.xi_bar[None,:],abs(state.delta_bar)[:,None])
 equilibrium_mask=system.kinetics._equilibrium_mask(energies,state.p,state.n,system.kinetics._temperature(energies,state.p))
 event_buffers=buffer_profile([len(state.delta_bar),7220],len(system.kinetics.event_i))
 result=dict(event_array_buffer_profile=event_buffers,event_count=len(system.kinetics.event_i),
  equilibrium_nodes_skipped_by_collision_kernel=int(np.count_nonzero(equilibrium_mask)),
  active_collision_nodes=int(np.count_nonzero(~equilibrium_mask)),status="PROFILE_ONLY",time_ps=current["time_ps"],checkpoint_sha256=hashlib.sha256(raw).hexdigest(),
  source_pins_sha256=meta["source_pins_sha256"],profile_helper_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
  threads=1,affinity_cpus=sorted(os.sched_getaffinity(0)),BLAS_threads=1,profile_step_ps=.02,profile_step_tolerance=2e-7,midpoint_status=stepstatus,midpoint_error=error,
  midpoint_wall_s=stepwall,midpoint_cpu_s=stepcpu,midpoint_receipt=receipt,evaluate_wall_s=times,evaluate_cpu_s=cputimes,
  evaluate_mean_wall_s=float(np.mean(times)),evaluate_mean_cpu_s=float(np.mean(cputimes)),profile_clock="process_time; excludes scheduler contention on the shared admitted CPU",evaluate_top=profile_summary(evprof)[:60],
  midpoint_top=profile_summary(prof)[:60],no_RHS_changes=True,no_admission=True,
  note="Cumulative timings overlap; do not sum nested entries. Numba kernels appear as opaque dispatcher calls.")
 probe.writej(out/"summary.json",result);print(json.dumps(result,indent=2))
if __name__=="__main__":main()
