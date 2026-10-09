"""User-run full stage5.1 affine Schur/rank audit, no time steps or reset.

The already prepared photon is never deposited again. All 112 energy modes,
all interior gap/population fields, both superconducting cuts and CM states
participate. Long column assembly/SVD belongs to the user. Atomic column
checkpoints support an explicit --resume with unchanged inputs and producers.
No report or detector admission is generated.
"""
from __future__ import annotations
import argparse, hashlib, json, multiprocessing as mp, os, sys, time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from sandbox.stage4_core.parallel_runtime import limit_thread_environment,linux_resources,resource_budget
limit_thread_environment()
import numpy as np
from pysnspd.experimental.joint_evolution import solve_affine_candidate
from sandbox.stage5_photon.joint_snapshot_inputs import build_initial_system

_SYSTEM=None
_BASE=None


def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for part in iter(lambda:f.read(1024*1024),b""):h.update(part)
    return h.hexdigest()


def write(path,value):
    path=Path(path)
    with path.open("x") as stream:json.dump(value,stream,indent=2,allow_nan=False);stream.write("\n")


def atomic_array(path,array):
    path=Path(path);temporary=path.with_suffix(path.suffix+".partial")
    with temporary.open("wb") as stream:np.save(stream,array,allow_pickle=False)
    temporary.replace(path)


def init_worker(cpus,counter):
    with counter.get_lock():
        index=counter.value;counter.value+=1
    os.sched_setaffinity(0,{cpus[index%len(cpus)]})


def column(index):
    start=time.monotonic();v=np.zeros(_SYSTEM.size);v[index]=1
    value=_SYSTEM.evaluate(v)
    return index,value.packed_rows-_BASE,time.monotonic()-start


def support_metrics(h,ht):
    u=np.stack((h[...,0]+h[...,1],h[...,0]-h[...,1]),axis=-1)
    rate=np.stack((ht[...,0]+ht[...,1],ht[...,0]-ht[...,1]),axis=-1)
    bounds=np.full_like(rate,np.inf)
    np.divide(1-u,rate,out=bounds,where=rate>0)
    np.divide(-1-u,rate,out=bounds,where=rate<0)
    outward=((u>=1)&(rate>0))|((u<=-1)&(rate<0))
    finite=bounds[np.isfinite(bounds)]
    return dict(outward_exhausted_support_count=int(np.count_nonzero(outward)),
        outward_exhausted_support_max_per_s=float(np.max(abs(rate[outward]),initial=0)),
        maximum_formal_Euler_step_s=float(np.min(finite)) if finite.size else None,
        no_clip_applied=True)


def source_pins():
    paths=[Path(__file__),ROOT/"sandbox/stage5_photon/joint_snapshot_inputs.py",
        ROOT/"sandbox/stage5_rhs/run_matrix_rhs_audit.py"]
    paths += [ROOT/"pysnspd/experimental"/name for name in (
        "joint_evolution.py","joint_nambu_rows.py","first_spectral_response.py",
        "adiabatic_matrix_usadel.py","retarded_spatial_usadel.py",
        "bulk_current_reference.py","causal_kinetic_regulator.py",
        "matrix_kinetic_rhs.py","nambu_moments.py","electrical_ports.py",
        "effective_material.py","thermal_spatial_usadel.py")]
    return {str(p.relative_to(ROOT)):sha(p) for p in paths}


def run(args):
    global _SYSTEM,_BASE
    output=args.output.resolve()
    if not output.is_relative_to(Path("/home/jdiaz/scratch")):
        raise ValueError("All cache/matrix/log/output files must stay in scratch")
    budget=resource_budget(linux_resources(),.9,.9,args.workers)
    start=time.monotonic()
    _SYSTEM=build_initial_system(args.initial,resource_memory_limit_bytes=budget["memory_limit_bytes"]//(budget["workers"]+1))
    build_s=time.monotonic()-start
    zero=_SYSTEM.evaluate(np.zeros(_SYSTEM.size));_BASE=zero.packed_rows
    build_s=time.monotonic()-start
    identity=dict(schema="pysnspd.stage51.joint_snapshot.v1",initial=str(args.initial.resolve()),
        source_pins=source_pins(),input_provenance=_SYSTEM.provenance,
        input_pins=_SYSTEM.provenance["pinned_files_sha256"],
        unknowns=_SYSTEM.size,energies=len(_SYSTEM.slices),resource_budget=budget,
        full_reset=False,physical_time_steps=0,source_applied_again=False,
        coordinates="local electrochemical/Wilson temporal gauge; rates in SI",
        material_scope="Initial local common FD/BE, zero e-ph transfer collision only; not a general no-collision detector",
        build_and_zero_s=build_s)
    # A pilot measures the complete physical operator, without a long matrix.
    if args.pilot:
        v=np.random.default_rng(510).normal(size=_SYSTEM.size)*.001
        start=time.monotonic();value=_SYSTEM.evaluate(v);response_s=time.monotonic()-start
        global_shift=_SYSTEM.evaluate(v+.01*_SYSTEM.gauge_direction)
        gauge=float(np.max(abs(global_shift.packed_rows-value.packed_rows)))
        if output.exists():raise FileExistsError("Fresh scratch pilot output required")
        output.mkdir(parents=True)
        write(output/"identity.json",identity)
        write(output/"pilot.json",dict(status="FULL_OPERATOR_PILOT_NOT_RANK_ADMISSION",
            build_and_zero_s=build_s,response_s=response_s,unknowns=_SYSTEM.size,
            gauge_shift_scaled_row_difference=gauge,
            rough_serial_column_seconds=response_s*_SYSTEM.size,
            rough_parallel_column_seconds=response_s*_SYSTEM.size/budget["workers"],
            svd_time_not_included=True,stage51_complete=False,physical_time_steps=0))
        print(json.dumps(json.loads((output/"pilot.json").read_text())),flush=True)
        return
    if not args.execute:
        print(json.dumps(identity,indent=2));return
    if output.exists():
        if not args.resume:raise FileExistsError("Fresh output or explicit unchanged-input --resume required")
        saved=json.loads((output/"identity.json").read_text())
        for key in ("schema","initial","source_pins","input_pins","unknowns","energies"):
            if saved[key]!=identity[key]:raise ValueError("Resume input/producer identity changed: "+key)
        if not np.array_equal(np.load(output/"base.npy"),_BASE):
            raise ValueError("Recomputed base differs from saved affine state")
        if (output/"summary.json").exists():raise ValueError("Completed result already exists; do not overwrite")
    else:
        output.mkdir(parents=True);(output/"columns").mkdir()
        write(output/"identity.json",identity);atomic_array(output/"base.npy",_BASE)
    # Reserve actual cache RSS for every process conservatively, even though
    # fork shares most read-only spectral arrays. No 1-GiB guess for a112 cache.
    rss=0
    for line in Path("/proc/self/status").read_text().splitlines():
        if line.startswith("VmRSS:"):rss=int(line.split()[1])*1024
    dense_bytes=8*_SYSTEM.size**2
    reservation=(budget["workers"]+1)*rss+5*dense_bytes+2*1024**3
    if reservation>budget["memory_limit_bytes"]:
        raise MemoryError("Joint cache+workers+SVD reservation exceeds90% available RAM; reduce --workers")
    done=[];pending=[]
    for j in range(_SYSTEM.size):
        file=output/"columns"/("column_%05d.npy"%j);receipt=file.with_suffix(".json")
        if file.exists() and receipt.exists():
            record=json.loads(receipt.read_text())
            if record["index"]!=j or record["sha256"]!=sha(file):
                raise ValueError("Saved column identity changed")
            a=np.load(file,allow_pickle=False)
            if a.shape!=(_SYSTEM.size,) or np.any(~np.isfinite(a)):
                raise ValueError("Invalid saved joint column")
            done.append(j)
        elif file.exists() or receipt.exists():
            raise ValueError("Incomplete checkpoint pair; preserve and inspect it before resume")
        else:pending.append(j)
    counter=mp.get_context("fork").Value("i",0);original=os.sched_getaffinity(0)
    # Spread workers across physical cores before using SMT siblings.
    groups=linux_resources()["physical_core_groups"]
    usable=set(budget["selected_affinity_cpus"])
    ordered=[g[0] for g in groups if g[0] in usable]
    ordered += [cpu for g in groups for cpu in g[1:] if cpu in usable]
    cpus=[cpu for cpu in ordered if cpu!=budget["coordinator_cpu"]][:budget["workers"]]
    os.sched_setaffinity(0,{budget["coordinator_cpu"]})
    started=time.monotonic();completed_now=0
    try:
        with ProcessPoolExecutor(budget["workers"],mp_context=mp.get_context("fork"),
                initializer=init_worker,initargs=(cpus,counter)) as pool:
            futures={pool.submit(column,j):j for j in pending}
            for future in as_completed(futures):
                try:
                    index,value,seconds=future.result()
                except BaseException:
                    # Stop the unscheduled matrix work on the first failure;
                    # preserve completed columns for an explicit resume.
                    for other in futures:
                        other.cancel()
                    raise
                file=output/"columns"/("column_%05d.npy"%index)
                atomic_array(file,value);write(file.with_suffix(".json"),
                    dict(index=index,sha256=sha(file),seconds=seconds))
                done.append(index);completed_now+=1
                if completed_now==1 or completed_now%25==0 or len(done)==_SYSTEM.size:
                    elapsed=time.monotonic()-started
                    print(json.dumps(dict(completed=len(done),total=_SYSTEM.size,elapsed_s=elapsed,
                        remaining_estimate_s=elapsed*(_SYSTEM.size-len(done))/completed_now)),flush=True)
        A=np.column_stack([np.load(output/"columns"/("column_%05d.npy"%j),allow_pickle=False)
            for j in range(_SYSTEM.size)])
        atomic_array(output/"matrix.npy",A)
        print(json.dumps(dict(status="ASSEMBLY_COMPLETE_START_SVD",size=_SYSTEM.size)),flush=True)
        x,diagnostic=solve_affine_candidate(A,_BASE,gauge_direction=_SYSTEM.gauge_direction,
            gauge_reference_index=_SYSTEM.size-1)
        singular=diagnostic.pop("singular_values")
        atomic_array(output/"singular_values.npy",singular)
        value=_SYSTEM.evaluate(x,save_fields=True)
        h=np.asarray([sl.h for sl in _SYSTEM.slices])
        support=support_metrics(h,value.distribution_rate_per_s)
        f=value.fields
        ratio=np.linalg.norm(f["R1"].reshape(len(_SYSTEM.slices),-1),axis=1)/np.linalg.norm(f["R0"].reshape(len(_SYSTEM.slices),-1),axis=1)
        kratio=np.linalg.norm(f["K1"].reshape(len(_SYSTEM.slices),-1),axis=1)/np.linalg.norm(f["K0"].reshape(len(_SYSTEM.slices),-1),axis=1)
        np.savez_compressed(output/"joint_candidate_fields.npz",unknowns_scaled=x,
            gap_rate_per_s=value.gap_rate_per_s,phi_V=value.phi_V,
            distribution_rate_per_s=value.distribution_rate_per_s,
            circuit_rate_SI_per_s=value.circuit_rate_SI_per_s,
            neutrality_covariant_V=value.moments.neutrality_V,
            gap_force_J_per_gapbar=value.moments.gap_force_J_per_gapbar,
            edge_current_A=value.moments.edge_current_A,
            nodal_current_divergence_A=value.moments.nodal_current_divergence_A,
            energy_J=_SYSTEM.rows.moments.energy,weights_J=_SYSTEM.rows.moments.weights,
            coordinates_m=_SYSTEM.initial_data["coordinates_m"],h=h,**f)
        balance=_SYSTEM.circuit.power_balance(_SYSTEM.y,value.moments.port_voltage_V,
            derivative=value.circuit_rate_SI_per_s)
        recomputed_max=float(np.max(abs(value.packed_rows)))
        recomputed_ok=recomputed_max<=1e-7
        joint_ok=bool(diagnostic["algebraically_solved_full_rank"] and recomputed_ok)
        write(output/"spectral_receipts.json",dict(records=list(value.spectral_metrics)))
        summary=dict(status="JOINT_CANDIDATE_PROJECTED_ROWS_SOLVED_PENDING_PHYSICAL_ADMISSION" if joint_ok
            else "JOINT_CANDIDATE_SINGULAR_OR_INCONSISTENT_NO_EVOLUTION",
            joint_candidate_recomputed_ok=recomputed_ok,
            algebraic_and_recomputed_rows_ok=joint_ok,
            recomputed_scaled_row_tolerance=1e-7,
            unprojected_full_nambu_screen_ok=value.full_nambu_residual_max_per_mass<=1e-6,
            unprojected_screen_tolerance_per_mass=1e-6,
            algebraic=diagnostic,recomputed_scaled_rows_max=recomputed_max,
            kinetic_residual_max=value.kinetic_residual_max,
            full_nambu_residual_max_per_mass=value.full_nambu_residual_max_per_mass,
            cut_neutrality_covariant_max_V=value.cut_neutrality_max_V,
            port_current_residual_A=value.moments.port_current_residual_A.tolist(),
            current_divergence_interior_max_A=float(np.max(abs(value.moments.interior_current_divergence_A),initial=0)),
            Vdev_V=value.moments.port_voltage_V,CM_rate_SI_per_s=value.circuit_rate_SI_per_s.tolist(),
            circuit_power_residual_W=balance.residual_W,maximum_R1_over_R0=float(np.max(ratio)),
            maximum_K1_over_K0=float(np.max(kratio)),support=support,
            stage51_complete=False,stage5_complete=False,energy_admission=False,time_admission=False,
            gap_kinetic_neutrality_ports_CM_solved_together=True,
            phonon_rhs_computed=value.phonon_rate_per_s is not None,
            general_material_evolution_admitted=False,
            material_rhs_contract=value.material_rhs_contract,source_applied_again=False,
            physical_time_steps=0,full_reset=False,elapsed_s=time.monotonic()-start,
            scratch_memory_reservation_bytes=reservation)
        write(output/"summary.json",summary);print(json.dumps(summary),flush=True)
    except BaseException as error:
        failure=output/("failure_%d.json"%time.time_ns())
        write(failure,dict(error=repr(error),completed_columns=len(done),resume_available=True,
            no_clipping_or_regularization=True))
        raise
    finally:os.sched_setaffinity(0,original)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--initial",type=Path,required=True);p.add_argument("--output",type=Path,required=True)
    p.add_argument("--workers",type=int,default=8);p.add_argument("--execute",action="store_true")
    p.add_argument("--pilot",action="store_true");p.add_argument("--resume",action="store_true")
    a=p.parse_args()
    if a.pilot and a.execute:raise ValueError("Pilot and long execute are distinct")
    run(a)


if __name__=="__main__":main()
