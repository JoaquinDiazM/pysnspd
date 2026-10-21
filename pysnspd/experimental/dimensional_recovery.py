"""Authenticated, provenance-preserving continuation of the known v1.2.1 join defect.

Only the runner and dimensional projection may differ from the original
release. Every material/kinetic/circuit/geometry source and parameter stays
identical. Original run directories are never edited.
"""
import hashlib,json,os,shutil
from pathlib import Path
import h5py
import numpy as np
from .effective_coupled import EffectiveState
from .dimensional_transition import require_overlap

LEGACY={
    "pysnspd/experimental/dimensional_monitors.py":
        "9b52344b40ecd608170b0f96b5bc51f189f3d4a841df01bfb01bc67b7a026366",
    "pysnspd/experimental/dimensional_archive.py":
        "07bc6a0569cc0a853b2c3143c035e551f9ff16a9a59b22e6fd4c3196a5d08af5",
    "pipelines/13_dimensional_photon_run.py":
        "dd0c138cee8f191b540f5d77cb987279c736301b35162e0fb2c8fe8a51b80cb8",
    "pysnspd/experimental/dimensional_transition.py":
        "e99ee319d0828b93d68177926d9004ada3abfa6d55623f3d5daebfea77cb978b",
}
RECOVERY_SOURCE="pysnspd/experimental/dimensional_recovery.py"

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1024*1024),b""):h.update(chunk)
    return h.hexdigest()

def prepare_continuation(source,destination,pins,identity):
    source=Path(source).resolve();destination=Path(destination).resolve()
    scratch=Path("/home/jdiaz/scratch").resolve()
    if not source.is_relative_to(scratch) or not destination.is_relative_to(scratch):
        raise ValueError("Continuation data must stay in scratch")
    if destination.exists():raise FileExistsError("Continuation destination already exists")
    manifest=json.loads((source/"manifest.json").read_text())
    old=manifest["source_pins_sha256"]
    if set(old)!=set(pins)-{RECOVERY_SOURCE,"pysnspd/experimental/dimensional_1d.py","pysnspd/experimental/dimensional_cold.py"}:
        raise ValueError("Unsupported checkpoint source catalog")
    for key,value in old.items():
        if value!=(LEGACY[key] if key in LEGACY else pins[key]):
            raise ValueError("Unapproved source change: "+key)
    legacy_identity={k:v for k,v in identity.items() if k not in ("monitor_collar_energy_fraction","no_cold_recovery")}
    if identity.get("monitor_collar_energy_fraction",.001)!=.001:
        raise ValueError("Unsupported new collar energy margin")
    if manifest["parameter_identity"]!=legacy_identity:
        raise ValueError("Continuation physical/numerical parameters differ")
    if (manifest.get("model_version")!="1.2.1" or manifest.get("phase")!="2d" or
        manifest.get("status")!="POSSIBLE_DIMENSIONAL_ERROR" or
        not str(manifest.get("error","")).startswith("Projection energy tolerance exceeded:") or
        not manifest.get("observation",{}).get("margin_complete")):
        raise ValueError("Not the supported pre-handoff join failure")
    checkpoint=source/"checkpoint_latest.npz"
    if sha(checkpoint)!=manifest["checkpoint_sha256"] or sha(source/"dc1d.npz")!=manifest["dc1d_sha256"]:
        raise ValueError("Source checkpoint or DC hash mismatch")
    with np.load(checkpoint,allow_pickle=False) as z:
        arrays={k:z[k].copy() for k in z.files}
    meta=json.loads(str(arrays["metadata_json"]))
    if meta["source_pins_sha256"]!=old or meta["parameter_identity"]!=legacy_identity:
        raise ValueError("Manifest/checkpoint identity mismatch")
    if meta["phase"]!="2d" or not meta.get("shadow_active") or abs(meta["shadow_time_ps"]-meta["time_ps"])>1e-8:
        raise ValueError("A synchronized native shadow is required")
    # Refuse a live writer; no interruption or advisory-lock bypass.
    for proc in Path("/proc").iterdir():
        if not proc.name.isdigit() or int(proc.name)==os.getpid():continue
        try:
            words=[w.decode() for w in (proc/"cmdline").read_bytes().split(bytes([0])) if w]
            if any(Path(w).name=="13_dimensional_photon_run.py" for w in words):
                if "--output-dir" in words and Path(words[words.index("--output-dir")+1]).resolve()==source:
                    raise ValueError("Source case is still active")
        except (OSError,UnicodeError):pass
    shutil.copytree(source,destination)
    provenance=dict(schema="pysnspd.v121_join_continuation.v1",source_run=str(source),
        source_checkpoint_sha256=manifest["checkpoint_sha256"],source_dc1d_sha256=manifest["dc1d_sha256"],
        source_pins_sha256=old,reused_until_ps=meta["time_ps"],reused_phase="2d",
        source_archive_sha256=sha(source/"data.h5"),
        photon_redeposited=False,physical_kernel_changed=False,
        shadow_policy="replay native 1D overlap with current-following cold cuts; old rigid-cut shadow is invalid",
        boundary_condition_changed_1d=True,original_2d_physics_changed=False)
    if sha(destination/"data.h5")!=provenance["source_archive_sha256"]:
        raise ValueError("Archive copy hash mismatch")
    (destination/"pre_update_manifest.json").write_text(json.dumps(manifest,indent=2))
    (destination/"pre_update_checkpoint.npz").write_bytes(checkpoint.read_bytes())
    meta.update(source_pins_sha256=pins,parameter_identity=identity,terminal_diagnostic=False,reaudit_saved_overlap=False,
        replay_saved_overlap=True,
        continuation=provenance)
    arrays["metadata_json"]=np.array(json.dumps(meta,allow_nan=False))
    temp=destination/"checkpoint_recovery.tmp"
    with temp.open("wb") as f:np.savez_compressed(f,**arrays)
    os.replace(temp,destination/"checkpoint_latest.npz")
    manifest.update(source_pins_sha256=pins,parameter_identity=identity,checkpoint_sha256=sha(destination/"checkpoint_latest.npz"),
        status="RECOVERY_PREPARED",error=None,continuation=provenance)
    (destination/"manifest.json").write_text(json.dumps(manifest,indent=2))
    print(json.dumps(dict(status="RECOVERY_PREPARED",source_run=str(source),
        output_run=str(destination),reused_until_ps=meta["time_ps"],photon_redeposited=False)),flush=True)

def reaudit_overlap(directory,project,scale,tolerances):
    """Use exact 2D snapshots, bounding shadow-time mismatch by BOTH brackets.

    No interpolation of a rapidly rotating complex condensate or evolution
    history is invented. Require each bracketing saved shadow to pass every
    physical gate against the exact saved 2D state. Record their time offsets.
    """
    results=[]
    with np.load(Path(directory)/"checkpoint_latest.npz",allow_pickle=False) as z:
        saved_checks=json.loads(str(z["metadata_json"]))["shadow_checks"]
    saved_times=np.array([c["time_ps"] for c in saved_checks])
    with h5py.File(Path(directory)/"data.h5","r",swmr=True) as f:
        original=f["2d"];shadow=f["shadow1d"]
        ts=shadow["snapshot_time_ps"][:];to=original["snapshot_time_ps"][:]
        if len(ts)<2:raise ValueError("Insufficient saved native overlap")
        def state(group,index):
            t=float(group["snapshot_time_ps"][index]);trace=group["trace"][:]
            j=int(np.argmin(abs(trace["time_ps"]-t)))
            if abs(float(trace["time_ps"][j])-t)>1e-8:
                raise ValueError("Snapshot has no accepted circuit state")
            cm=np.array([trace["Ib_A"][j],trace["Is_A"][j],trace["vc_V"][j]])
            return EffectiveState(*(group[k][index] for k in ("delta_bar","p","n")),cm)
        for i,t in enumerate(to):
            if t<ts[0]-1e-9 or t>ts[-1]+1e-9:continue
            j=int(np.searchsorted(ts,t))
            indices=sorted({max(0,min(len(ts)-1,j-1)),min(len(ts)-1,j)})
            if min(abs(ts[indices]-t))<1e-9:
                indices=[indices[int(np.argmin(abs(ts[indices]-t)))]]
            elif ts[indices[-1]]-ts[indices[0]]>.11:
                raise ValueError("Saved shadow bracket too wide")
            st=state(original,i)
            checks=[]
            for k in indices:
                sh=state(shadow,k)
                c=project.errors(st,sh,scale)
                # Work/energy and circuit comparisons MUST use simultaneous
                # states. The old checkpoint already contains these exact-time
                # scalars. Taking energies from different times confuses real
                # Joule work with projection error.
                n=int(np.searchsorted(saved_times,t))
                same_time_indices=sorted({max(0,min(len(saved_checks)-1,n-1)),
                                         min(len(saved_checks)-1,n)})
                if min(abs(saved_times[same_time_indices]-t))>.11:
                    raise ValueError("Missing simultaneous scalar overlap checks")
                c.update(time_ps=float(t),shadow_sample_time_ps=float(ts[k]),
                    simultaneous_scalar_bracket_ps=[float(saved_times[q]) for q in same_time_indices])
                for key in ("film_energy_relative","voltage_error_V","current_relative"):
                    c[key]=max(saved_checks[q][key] for q in same_time_indices)
                require_overlap(c,tolerances)
                checks.append(c)
            worst={key:max(c[key] for c in checks) for key in checks[0]
                   if key not in ("time_ps","shadow_sample_time_ps","gauge_alignment_rad","simultaneous_scalar_bracket_ps")}
            worst.update(time_ps=float(t),shadow_bracket_ps=[float(ts[k]) for k in indices],
                simultaneous_scalar_bracket_ps=checks[0]["simultaneous_scalar_bracket_ps"],
                comparison="exact 2D, bracketing native shadow fields; stored simultaneous energy/circuit checks")
            results.append(worst)
    if len(results)<2:raise ArithmeticError("Insufficient re-audited overlap")
    return results

def replay_overlap(directory,project,reference,m,args):
    """Replay only the inexpensive 1D shadow; reuse every accepted 2D state.

    Start at the last stored 2D snapshot at/before crossing (<=0.1ps earlier),
    keeping the actual populations and circuit. Cold added tails are the
    explicit exterior hypothesis. The invalid old rigid-cut shadow is retained
    in the original archive, and replaced only in this continuation copy.
    """
    import time
    from .dimensional_transition import doubled_step
    from .dimensional_archive import PhaseArchive
    path=Path(directory)/"data.h5"
    with h5py.File(path,"r",swmr=True) as f:
        q=f["2d"];times=q["snapshot_time_ps"][:];trace=q["trace"][:]
        cross=m["controller"]["state"]["crossing_time_s"]*1e12
        begin=max(0,int(np.searchsorted(times,cross,side="right"))-1)
        indices=[i for i in range(begin,len(times)) if times[i]<=m["time_ps"]+1e-9]
        originals=[]
        for i in indices:
            t=float(times[i]);j=int(np.argmin(abs(trace["time_ps"]-t)))
            if abs(float(trace["time_ps"][j])-t)>1e-8:raise ValueError("Missing accepted 2D circuit snapshot")
            cm=np.array([trace["Ib_A"][j],trace["Is_A"][j],trace["vc_V"][j]])
            originals.append((t,EffectiveState(*(q[k][i] for k in ("delta_bar","p","n")),cm),
                              float(trace["device_absolute_J"][j]),float(trace["absolute_external_J"][j])))
    if len(originals)<3 or originals[0][0]>cross or cross-originals[0][0]>.11:
        raise ValueError("Insufficient saved 2D overlap for replay")
    from .dimensional_transition import TransitionTolerances
    tol=TransitionTolerances(**{k:getattr(args,"guard_"+k) for k in TransitionTolerances.__dataclass_fields__})
    shadow,initial=project.project(originals[0][1],reference+originals[0][2],tol)
    t=originals[0][0];hnext=args.dt_ps;s=project.s1;rhs=s.evaluate(shadow)
    counts=m["archive_counts"]
    counts["shadow1d"]={k:0 for k in counts["shadow1d"]}
    ar=PhaseArchive(path,{"2d":project.s2,"shadow1d":s,"1d":s},
        {"2d":project.g2,"shadow1d":project.g1,"1d":project.g1},
        {"2d":project.dc2,"shadow1d":project.dc1,"1d":project.dc1},resume=True,counts=counts)
    ar.file.attrs["model_version"]="1.2.2"
    checks=[];accepted=rejected=0;start=time.monotonic();last=start
    try:
        for target,original,device_work,absolute_work in originals:
            while t<target-1e-10:
                h=min(hnext,target-t)
                try:
                    new,rr,receipt=doubled_step(s,shadow,rhs,h,time_rtol=args.time_rtol,
                        time_atol=args.time_atol,tolerance=args.step_tolerance,
                        max_iterations=args.max_iterations,energy_tolerance=args.local_energy_tolerance,
                        reference=reference,absolute_work=absolute_work)
                except (ArithmeticError,ValueError):
                    rejected+=1;hnext=h/2
                    if hnext<args.minimum_adaptive_step_ps:raise ArithmeticError("Replayed shadow resolution exhausted")
                    continue
                shadow,rhs=new,rr;t+=h;accepted+=1
                hnext=min(args.dt_ps,h*min(2.,max(1.,.8/max(receipt["temporal_error"],1e-12)**(1/3))))
                if time.monotonic()-last>=30:
                    print(json.dumps(dict(status="REPLAYING_1D_OVERLAP",time_ps=t,target_ps=originals[-1][0],
                        accepted_steps=accepted,rejections=rejected,elapsed_wall_s=time.monotonic()-start)),flush=True)
                    last=time.monotonic()
            c=project.errors(original,shadow,reference+device_work)
            c.update(time_ps=t,voltage_error_V=abs(project.s2.circuit.R_load_ohm*(
                (original.circuit_SI[0]-original.circuit_SI[1])-(shadow.circuit_SI[0]-shadow.circuit_SI[1]))),
                current_relative=float(np.max(abs(original.circuit_SI[:2]-shadow.circuit_SI[:2]))/
                                       project.s2.current_scale_A))
            require_overlap(c,tol);checks.append(c)
            ar.append("shadow1d",dict(time_ps=t,Ib_A=shadow.circuit_SI[0],Is_A=shadow.circuit_SI[1],
                vc_V=shadow.circuit_SI[2],Vout_V=s.circuit.R_load_ohm*(shadow.circuit_SI[0]-shadow.circuit_SI[1]),
                Vdev_V=rhs.electrical.Vdev_V,energy_J=rhs.energy_J))
            ar.snapshot("shadow1d",s,shadow,rhs,t)
        m["archive_counts"]=ar.counts()
    finally:ar.close()
    if abs(t-m["time_ps"])>1e-8:raise ArithmeticError("Replay did not reach exact handoff")
    m["continuation"]["replayed_shadow"]=dict(start_ps=originals[0][0],end_ps=t,
        accepted_steps=accepted,rejections=rejected,wall_s=time.monotonic()-start,
        old_rigid_boundary_shadow_discarded=True)
    print(json.dumps(dict(status="1D_OVERLAP_REPLAY_COMPLETE",**m["continuation"]["replayed_shadow"])),flush=True)
    return shadow,checks
