"""User-run coupled first-step IMEX Rosenbrock candidate campaign.

Reuses prepared photon and saved spectral/Matsubara roots, never the old
ill-scaled truncated candidate. Each dt solves the physical linearized
gap/kinetic/neutrality/ports/CM equations. Finite states are re-solved on the
actual changed gap and links, when Pauli support and current direction permit.

This is the FIRST photon step only: C(y0)=0 at common local FD/BE is explicit.
It cannot evolve a general second state without the full material provider.
No full reset, report, detector admission or energy-residual heating.
"""
from __future__ import annotations
import argparse,hashlib,json,multiprocessing as mp,os,sys,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from sandbox.stage4_core.parallel_runtime import limit_thread_environment,linux_resources,resource_budget
limit_thread_environment()
import numpy as np
from scipy.constants import Boltzmann as kB
from pysnspd.experimental.joint_rosenbrock import JointRosenbrockSystem
from pysnspd.experimental.joint_state_tangent import ThermalBulkTangentBundle
from pysnspd.experimental.joint_linear_solve import solve_equilibrated_affine_candidate
from sandbox.stage5_photon.joint_snapshot_inputs import build_initial_system
from sandbox.stage5_photon.run_joint_snapshot import sha,write,atomic_array,init_worker
_SYSTEM=None;_BASE=None


def column(index):
    start=time.monotonic();v=np.zeros(_SYSTEM.size);v[index]=1.
    return index,_SYSTEM.evaluate(v).packed_rows-_BASE,time.monotonic()-start


def source_pins():
    paths=list((ROOT/"pysnspd/experimental").glob("*.py"))
    paths += [Path(__file__),ROOT/"sandbox/stage5_photon/joint_snapshot_inputs.py",
        ROOT/"sandbox/stage5_photon/joint_candidate_rebuild.py",
        ROOT/"sandbox/stage5_photon/run_joint_snapshot.py",
        ROOT/"sandbox/stage5_rhs/run_matrix_rhs_audit.py",
        ROOT/"sandbox/stage4_core/parallel_runtime.py"]
    return {str(p.relative_to(ROOT)):sha(p) for p in paths}


def build_thermal_bundle(snapshot):
    data=snapshot.initial_data;w=data["inputs"];p=w["summary"]["plan"]
    observer=Path(data["identity"]["plan"]["thermal_observer"])
    fields=[];eps=[];pins={}
    for i in range(256):
        file=observer/"spectral_receipts"/("mode_%05d.npz"%i)
        with np.load(file,allow_pickle=False) as z:
            fields.append(z["u"].copy());eps.append(float(z["epsilon_bar"]))
        pins[str(file)]=sha(file)
    bundle=ThermalBulkTangentBundle(snapshot.graph,snapshot.gap,p["T_K"]/p["Tc_K"],
        np.asarray(eps),np.asarray(fields),snapshot.alpha,fixed_nodes=snapshot.fixed,
        bulk_q_bar=snapshot.bulk.q_bar,energy_scale_J=w["observed"]["energy_scale_J"],
        current_unit_A=snapshot.rows.moments.I0)
    return bundle,pins


def initial_snapshot_check(snapshot,path):
    old=json.loads((path/"identity.json").read_text())
    if old["input_pins"]!=snapshot.provenance["pinned_files_sha256"]:
        raise ValueError("Completed snapshot and prepared numeric/reference inputs differ")
    base=np.load(path/"base.npy",allow_pickle=False)
    now=snapshot.evaluate(np.zeros(snapshot.size)).packed_rows
    difference=float(np.max(abs(now-base)))
    if difference>1e-7:
        raise ValueError("Current dt=0 physical equations disagree with the saved snapshot")
    # Two different physical columns guard the reused first-gradient convention.
    A=np.load(path/"matrix.npy",mmap_mode="r",allow_pickle=False)
    errors=[]
    for j in (0,2*snapshot.nf):
        v=np.zeros(snapshot.size);v[j]=1.
        error=float(np.max(abs(snapshot.evaluate(v).packed_rows-now-A[:,j])))
        if error>1e-7:raise ValueError("Current derivative convention changed from stored snapshot")
        errors.append(dict(column=j,maximum_difference=error))
    return dict(path=str(path),identity_sha256=sha(path/"identity.json"),
        summary_sha256=sha(path/"summary.json"),base_maximum_difference=difference,
        column_checks=errors,old_matrix_reassembled=False,
        old_candidate_not_used=True,old_producer_pins=old["source_pins"])


def finite_support(h,ht,dt_s):
    new=h+dt_s*ht
    u=np.stack((new[...,0]+new[...,1],new[...,0]-new[...,1]),axis=-1)
    excess=float(np.max(abs(u)-1.,initial=0.))
    return new,dict(maximum_h_plus_minus_excess=excess,
        exact_Pauli_support_ok=bool(excess<=32*np.finfo(float).eps),
        diagnostic_rebuild_support_ok=bool(excess<=1e-7),
        diagnostic_rebuild_excess_limit=1e-7,
        support_tolerance_not_relaxed_for_admission=True,no_clipping=True)


def memory_reservation(budget,size):
    rss=0
    for line in Path("/proc/self/status").read_text().splitlines():
        if line.startswith("VmRSS:"):rss=int(line.split()[1])*1024
    reservation=(budget["workers"]+1)*rss+5*8*size**2+2*1024**3
    if reservation>budget["memory_limit_bytes"]:
        raise MemoryError("Cache+thermal tangent+workers+QR exceed90% available memory; reduce workers")
    return reservation


def one_step(snapshot,bundle,output,dt_s,budget,args,identity):
    global _SYSTEM,_BASE
    start=time.monotonic()
    _SYSTEM=JointRosenbrockSystem(snapshot,dt_s,bundle,
        tau_escape_s=getattr(args,"tau_escape_ps",15.)*1e-12)
    _BASE=_SYSTEM.evaluate(np.zeros(_SYSTEM.size)).packed_rows
    record=dict(identity,dt_s=dt_s,source_pins=source_pins(),
        method="first-step IMEX linearly implicit Rosenbrock-Euler; not nonlinear backward Euler",
        original_residual_tolerance=1e-7)
    if output.exists():
        if not args.resume:raise FileExistsError("Fresh output or --resume required")
        old=json.loads((output/"identity.json").read_text())
        for key in ("dt_s","source_pins","input_pins","thermal_input_pins","snapshot_check"):
            if old[key]!=record[key]:raise ValueError("Resume identity changed: "+key)
        if (output/"summary.json").exists():
            return json.loads((output/"summary.json").read_text())
        if not np.array_equal(np.load(output/"base.npy"),_BASE):
            raise ValueError("Changed recomputed first-step base")
    else:
        output.mkdir(parents=True);(output/"columns").mkdir()
        write(output/"identity.json",record);atomic_array(output/"base.npy",_BASE)
    reservation=memory_reservation(budget,_SYSTEM.size)
    done=[];pending=[]
    for j in range(_SYSTEM.size):
        file=output/"columns"/("column_%05d.npy"%j);receipt=file.with_suffix(".json")
        if file.exists() and receipt.exists():
            r=json.loads(receipt.read_text());a=np.load(file,allow_pickle=False)
            if r["index"]!=j or r["sha256"]!=sha(file) or a.shape!=(_SYSTEM.size,) or np.any(~np.isfinite(a)):
                raise ValueError("Invalid changed column checkpoint")
            done.append(j)
        elif file.exists() or receipt.exists():
            raise ValueError("Incomplete checkpoint pair; preserve for inspection")
        else:pending.append(j)
    resources=linux_resources();usable=set(budget["selected_affinity_cpus"])
    ordered=[g[0] for g in resources["physical_core_groups"] if g[0] in usable]
    ordered += [c for g in resources["physical_core_groups"] for c in g[1:] if c in usable]
    cpus=[c for c in ordered if c!=budget["coordinator_cpu"]][:budget["workers"]]
    original=os.sched_getaffinity(0);os.sched_setaffinity(0,{budget["coordinator_cpu"]})
    assembly=time.monotonic();completed=0;counter=mp.get_context("fork").Value("i",0)
    try:
        with ProcessPoolExecutor(budget["workers"],mp_context=mp.get_context("fork"),
                initializer=init_worker,initargs=(cpus,counter)) as pool:
            futures={pool.submit(column,j):j for j in pending}
            for future in as_completed(futures):
                try:index,value,seconds=future.result()
                except BaseException:
                    for other in futures:other.cancel()
                    raise
                file=output/"columns"/("column_%05d.npy"%index)
                atomic_array(file,value);write(file.with_suffix(".json"),dict(index=index,sha256=sha(file),seconds=seconds))
                done.append(index);completed+=1
                if completed==1 or completed%25==0 or len(done)==_SYSTEM.size:
                    elapsed=time.monotonic()-assembly
                    print(json.dumps(dict(dt_s=dt_s,completed=len(done),total=_SYSTEM.size,
                        elapsed_s=elapsed,remaining_estimate_s=elapsed*(_SYSTEM.size-len(done))/completed)),flush=True)
        A=np.column_stack([np.load(output/"columns"/("column_%05d.npy"%j),allow_pickle=False) for j in range(_SYSTEM.size)])
        atomic_array(output/"matrix.npy",A)
        x,diag=solve_equilibrated_affine_candidate(A,_BASE,gauge_direction=_SYSTEM.gauge_direction,
            gauge_reference_index=_SYSTEM.size-1)
        value=_SYSTEM.evaluate(x,save_fields=True)
        h=np.asarray([sl.h for sl in _SYSTEM.slices]);hnew,support=finite_support(h,value.distribution_rate_per_s,dt_s)
        gapnew=_SYSTEM.gap+dt_s*value.fields["covariant_gap_rate_per_s"]
        alphanew=_SYSTEM.alpha+dt_s*value.fields["alpha_time_per_s"]
        cmnew=_SYSTEM.y+dt_s*value.circuit_rate_SI_per_s
        current_direction_ok=bool(cmnew[0]>0 and cmnew[1]>0)
        ratio=np.linalg.norm(value.fields["R1"].reshape(len(_SYSTEM.slices),-1),axis=1)/np.linalg.norm(value.fields["R0"].reshape(len(_SYSTEM.slices),-1),axis=1)
        kratio=np.linalg.norm(value.fields["K1"].reshape(len(_SYSTEM.slices),-1),axis=1)/np.maximum(np.linalg.norm(value.fields["K0"].reshape(len(_SYSTEM.slices),-1),axis=1),np.finfo(float).tiny)
        residual=float(np.max(abs(value.packed_rows)))
        np.savez_compressed(output/"candidate_fields.npz",x=x,gap_bar_new=gapnew,alpha_new=alphanew,
            h_new=hnew,circuit_state_new_SI=cmnew,h_initial=h,dt_s=dt_s,
            distribution_rate_per_s=value.distribution_rate_per_s,circuit_rate_SI_per_s=value.circuit_rate_SI_per_s,
            phi_V=value.phi_V,phonon_energy_J=_SYSTEM.phonon_energy_J,
            phonon_weights_J=_SYSTEM.phonon_weights_J,phonon_DOS_per_J_m3=_SYSTEM.phonon_DOS_per_J_m3,
            phonon_occupation_new=_SYSTEM.phonon_occupation,
            phonon_occupation_initial=_SYSTEM.phonon_occupation_initial,
            phonon_bath=_SYSTEM.phonon_bath,phonon_rate_per_s=_SYSTEM.phonon_rate,
            phonon_escape_power_W=_SYSTEM.phonon_escape_power_W,
            tau_escape_s=_SYSTEM.tau_escape_s,**value.fields)
        rebuild=dict(status="FINITE_STATE_NOT_REBUILT_OUTSIDE_DIAGNOSTIC_SUPPORT_LIMIT_OR_CURRENT_DIRECTION")
        if support["diagnostic_rebuild_support_ok"] and current_direction_ok:
            from sandbox.stage5_photon.joint_candidate_rebuild import rebuild_candidate
            os.sched_setaffinity(0,original)
            rebuild=rebuild_candidate(snapshot,output/"candidate_fields.npz",output/"rebuild",workers=budget["workers"])
        write(output/"spectral_receipts.json",dict(records=list(value.spectral_metrics)))
        summary=dict(status="FIRST_COUPLED_STEP_CANDIDATE_REQUIRES_REVIEW",
            dt_s=dt_s,algebraic=diag,recomputed_scaled_rows_max=residual,
            original_row_tolerance=1e-7,
            projected_rows_ok=bool(diag["algebraically_solved_full_rank"] and residual<=1e-7),
            kinetic_residual_max=value.kinetic_residual_max,full_nambu_residual_max_per_mass=value.full_nambu_residual_max_per_mass,
            maximum_R1_over_R0=float(np.max(ratio)),maximum_K1_over_K0=float(np.max(kratio)),
            support=support,current_direction_ok=current_direction_ok,
            circuit_state_new_SI=cmnew.tolist(),Vdev_V=value.moments.port_voltage_V,
            maximum_gap_increment_bar=float(np.max(abs(gapnew-_SYSTEM.gap))),
            maximum_population_increment=float(np.max(abs(hnew-h))),
            port_current_residual_A=value.moments.port_current_residual_A.tolist(),
            finite_state_rebuild=rebuild,source_applied_again=False,
            general_material_rhs_available=False,phonon_initial_transfer_rate_zero=True,
            phonon_escape_scope="Inherited effective interface input, freely configurable; infinity denotes closed control",
            tau_escape_s=None if np.isinf(_SYSTEM.tau_escape_s) else _SYSTEM.tau_escape_s,
            phonon_escape_power_W=_SYSTEM.phonon_escape_power_W,
            Born_spectral_ordering="Inherited weak-Born: no physical inelastic width resummed in R0; not generally admitted",
            stage51_complete=False,stage5_complete=False,energy_admission=False,time_admission=False,
            accepted_physical_time_steps=0,full_reset=False,elapsed_s=time.monotonic()-start,
            scratch_memory_reservation_bytes=reservation)
        write(output/"summary.json",summary);print(json.dumps(summary),flush=True)
        return summary
    except BaseException as error:
        write(output/("failure_%d.json"%time.time_ns()),dict(error=repr(error),completed_columns=len(done),
            resume_available=True,no_clipping_or_regularization=True))
        raise
    finally:os.sched_setaffinity(0,original)


def run(args):
    start=time.monotonic();output=args.output.resolve()
    if not output.is_relative_to(Path("/home/jdiaz/scratch")):
        raise ValueError("All new numeric data must remain in scratch")
    if any(not np.isfinite(t) or t<=0 for t in args.dt_ps) or len(set(args.dt_ps))!=len(args.dt_ps):
        raise ValueError("Distinct positive physical step sizes required")
    tau_ps=getattr(args,"tau_escape_ps",15.)
    if np.isnan(tau_ps) or tau_ps<=0:
        raise ValueError("Positive physical phonon escape time or inf required")
    budget=resource_budget(linux_resources(),.9,.9,args.workers)
    snapshot=build_initial_system(args.initial,resource_memory_limit_bytes=budget["memory_limit_bytes"]//(budget["workers"]+1))
    # This inexpensive old-state check reuses exactly two existing columns.
    check=initial_snapshot_check(snapshot,args.snapshot.resolve())
    bundle,pins=build_thermal_bundle(snapshot)
    identity=dict(schema="pysnspd.stage51.coupled_first_step.v1",initial=str(args.initial.resolve()),
        input_pins=snapshot.provenance["pinned_files_sha256"],thermal_input_pins=pins,
        snapshot_check=check,source_pins=source_pins(),resource_budget=budget,unknowns=snapshot.size,
        all_energies=len(snapshot.slices),all_nodes=snapshot.graph.n_nodes,dt_ps=args.dt_ps,
        tau_escape_ps=None if np.isinf(tau_ps) else tau_ps,
        escape_parameter_source="docs/modelo_v0_4/D_sintesis_plan_y_verificaciones_v0_4.md D.19:15ps effective inherited interface input, not Allmaras universal nominal",
        source_applied_again=False,full_reset=False,accepted_physical_time_steps=0,
        calibrators="Allmaras2020 nominal tau_ee(Tc)=5ps and tau_ep(Tc)=24.7ps remain free calibrators")
    if args.pilot:
        step=JointRosenbrockSystem(snapshot,args.dt_ps[0]*1e-12,bundle,tau_escape_s=tau_ps*1e-12)
        t=time.monotonic();base=step.evaluate(np.zeros(step.size)).packed_rows
        v=np.random.default_rng(511).normal(size=step.size)*1e-4
        value=step.evaluate(v);elapsed=time.monotonic()-t
        gauge=step.evaluate(v+.01*step.gauge_direction)
        if output.exists():raise FileExistsError("Fresh pilot output required")
        output.mkdir(parents=True);write(output/"identity.json",identity)
        write(output/"pilot.json",dict(status="FULL_COUPLED_STEP_OPERATOR_PILOT_ONLY",two_evaluations_s=elapsed,
            gauge_max_difference=float(np.max(abs(gauge.packed_rows-value.packed_rows))),
            size=step.size,energies=len(step.slices),matsubara_modes=len(bundle.modes),
            total_s=time.monotonic()-start,stage51_complete=False,
            memory_reservation_bytes=memory_reservation(budget,step.size)))
        print((output/"pilot.json").read_text(),flush=True);return
    if not args.execute:
        print(json.dumps(identity,indent=2));return
    if output.exists():
        if not args.resume:raise FileExistsError("Fresh campaign output or --resume required")
        old=json.loads((output/"identity.json").read_text())
        for key in ("schema","initial","input_pins","thermal_input_pins","source_pins","snapshot_check","dt_ps","tau_escape_ps"):
            if old[key]!=identity[key]:raise ValueError("Campaign identity changed: "+key)
        if (output/"summary.json").exists():raise ValueError("Campaign completed; preserve it")
    else:
        output.mkdir(parents=True);write(output/"identity.json",identity)
    summaries=[]
    for i,dt_ps in enumerate(args.dt_ps):
        summaries.append(one_step(snapshot,bundle,output/("step_%02d"%i),dt_ps*1e-12,budget,args,identity))
    summary=dict(status="COUPLED_FIRST_STEP_CAMPAIGN_COMPLETED_REVIEW_REQUIRED",
        steps=summaries,dt_ps=args.dt_ps,stage51_complete=False,stage5_complete=False,
        energy_admission=False,time_admission=False,accepted_physical_time_steps=0,full_reset=False,
        source_applied_again=False,elapsed_s=time.monotonic()-start,
        decision="Compare increment(dt) and finite-state residual; dt-independent jumps indicate an initial algebraic layer, not validated photon dynamics")
    write(output/"summary.json",summary);print(json.dumps(summary),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("initial","snapshot","output"):p.add_argument("--"+name,type=Path,required=True)
    p.add_argument("--dt-ps",type=float,nargs="+",default=[.001,.0005])
    p.add_argument("--tau-escape-ps",type=float,default=15.,
        help="Effective inherited D.19 interface input; free parameter. inf selects explicit closed control.")
    p.add_argument("--workers",type=int,default=8);p.add_argument("--execute",action="store_true")
    p.add_argument("--resume",action="store_true");p.add_argument("--pilot",action="store_true")
    a=p.parse_args()
    if a.execute and a.pilot:raise ValueError("Pilot and long campaign are distinct")
    run(a)


if __name__=="__main__":main()
