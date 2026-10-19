"""v1.2: guarded full 2D -> full 1D evolution, native kinetic states retained."""
from __future__ import annotations
import argparse,importlib.util,json,math,os,signal,sys,time
from dataclasses import asdict
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="1"
os.environ.setdefault("NUMBA_CACHE_DIR","/home/jdiaz/scratch/tmp/pysnspd/numba_cache")
os.environ.setdefault("PYTHONPYCACHEPREFIX","/home/jdiaz/scratch/tmp/pysnspd/pycache")
early=argparse.ArgumentParser(add_help=False)
early.add_argument("--threads",type=int,default=int(os.environ.get("NUMBA_NUM_THREADS","4")))
os.environ["NUMBA_NUM_THREADS"]=str(early.parse_known_args()[0].threads)
spec=importlib.util.spec_from_file_location("_v12_physics",ROOT/"pipelines/09_thermal_stiffness_photon_run.py")
OLD=importlib.util.module_from_spec(spec);spec.loader.exec_module(OLD)
B=OLD.BASE
import numpy as np
from pysnspd.experimental.dimensional_transition import (
    TransitionTolerances,chain_geometry,DimensionalProjection,observed,
    checked_step,doubled_step,film_energy)
from pysnspd.experimental.dimensional_archive import PhaseArchive
from pysnspd.experimental.dimensional_geometry import verified_geometry
_LEGACY_GEOMETRY=OLD.calibrated_geometry
OLD.calibrated_geometry=lambda args,output,ell:verified_geometry(args,output,ell,_LEGACY_GEOMETRY)
from pysnspd.experimental.effective_coupled import EffectiveState
from pysnspd.experimental.resummed_thermal import ResummedThermalFunctional
from pysnspd.experimental.thermal_coupled import ThermalCoupledSystem

write_json=OLD._ORIGINAL_WRITE_JSON


def parser():
    p=OLD.parser()
    p.description=__doc__
    p.set_defaults(length_nm=300.,extra_time_ps=2.,confirmation_ps=.5,threads=4,
                   retained_fraction=.8,tau_escape_ps=20.,tau_ee_Tc_ps=4.)
    p.add_argument("--final-time-ns",type=float,default=15.)
    p.add_argument("--spacing-1d-nm",type=float,default=5.)
    p.add_argument("--maximum-1d-step-ps",type=float,default=10.)
    p.add_argument("--global-energy-tolerance",type=float,default=.01)
    p.add_argument("--time-rtol",type=float,default=2e-4)
    p.add_argument("--time-atol",type=float,default=2e-6)
    p.add_argument("--snapshot-1d-every-ps",type=float,default=2.,
        help="Target map cadence; 1D saves actual accepted times without limiting its integration step")
    p.add_argument("--checkpoint-1d-every-ps",type=float,default=100.)
    p.add_argument("--shadow-check-every-ps",type=float,default=.1)
    p.add_argument("--physical-cores",type=int,default=2)
    for key,value in asdict(TransitionTolerances()).items():
        p.add_argument("--guard-"+key.replace("_","-"),type=float,default=value)
    return p


def source_pins():
    pins=OLD.source_pins()
    for rel in ("pipelines/13_dimensional_photon_run.py",
        "pysnspd/experimental/dimensional_transition.py","pysnspd/experimental/dimensional_archive.py","pysnspd/experimental/dimensional_geometry.py"):
        pins[rel]=B.sha256(ROOT/rel)
    return pins


def identity(a):
    return B.safe_json({k:v for k,v in vars(a).items() if k not in
        ("output_dir","resume","cpu_start","duration_ps")})


def build_1d(a,s2,g2,output,resume):
    f2,k2=s2.functional,s2.kinetics
    g=chain_geometry(g2.length_m,g2.width_m,a.spacing_1d_nm*1e-9,f2.ell0_m)
    kwargs=dict(graph=g.graph,energy_scale_J=f2.Eu_J,N0_per_J_m3=f2.N0_per_J_m3,
        ell0_m=f2.ell0_m,thickness_m=f2.thickness_m,count_nodes=f2.xi_bar,
        count_weights=f2.weights_bar,gradient_law="bcs_radial",
        bath_temperature_K=a.bath_K,diffusion_m2_s=k2.D_m2_s)
    f=ResummedThermalFunctional(**kwargs)
    path=output/"dc1d.npz"
    if resume:
        with np.load(path,allow_pickle=False) as data:
            dc_d=data["delta_bar"].copy();dc_p=data["p"].copy()
    else:
        bath=B.EffectiveBCSFunctional(**kwargs)
        dc=B.prepare_effective_dc(bath,left_nodes=g.left_nodes,right_nodes=g.right_nodes,
            temperature_K=a.bath_K,current_A=a.current_uA*1e-6,width_m=g.width_m,
            sheet_resistance_ohm=a.sheet_resistance_ohm)
        dc_d,dc_p=dc.delta_bar,dc.p
        B.atomic_npz(path,dict(delta_bar=dc_d,p=dc_p))
    k=B.EffectiveBCSKinetics(xi_bar=f.xi_bar,weights_bar=f.weights_bar,
        phonons=k2.phonons,Eu_J=f.Eu_J,N0_per_J_m3=f.N0_per_J_m3,D_m2_s=k2.D_m2_s,
        volumes_m3=f.volumes_m3,edges=g.graph.edges,
        face_areas_m2=g.face_areas_m2(f.thickness_m),edge_lengths_m=g.edge_lengths_m,
        bath_temperature_K=a.bath_K,tau_ee_Tc_s=a.tau_ee_Tc_ps*1e-12,
        tau_escape_s=a.tau_escape_ps*1e-12,transport_order=a.transport_order,numba_threads=a.threads)
    ports=B.EffectiveNeutralPorts(g.graph,sheet_resistance_ohm=a.sheet_resistance_ohm,
        left_nodes=g.left_nodes,right_nodes=g.right_nodes)
    s=observed(ThermalCoupledSystem(f,k,ports,s2.circuit,tau_ee_Tc_ps=a.tau_ee_Tc_ps,
        tau_ep_Tc_ps=a.tau_ep_Tc_ps,current_scale_A=a.current_uA*1e-6,
        kinetic_residual_energy_scale_J=s2.kinetic_residual_energy_scale_J))
    state=EffectiveState(dc_d,dc_p,np.broadcast_to(k.bath_phonons,
        (g.graph.n_nodes,k.phonon_size)).copy(),np.array([a.current_uA*1e-6,a.current_uA*1e-6,0.]))
    s.evaluate(state)  # capacity, ellipticity, current and native support guards
    return s,g,state


def run(a):
    OLD.validate_args(a)
    for key in ("final_time_ns","spacing_1d_nm","maximum_1d_step_ps","time_rtol",
                "time_atol","snapshot_1d_every_ps","checkpoint_1d_every_ps","shadow_check_every_ps","global_energy_tolerance"):
        if not np.isfinite(getattr(a,key)) or getattr(a,key)<=0:raise ValueError(key+" must be positive")
    if a.extra_time_ps<=a.confirmation_ps:
        raise ValueError("The overlap requires extra time greater than confirmation time")
    if a.final_time_ns*1e3<=a.ceiling_ps+a.extra_time_ps:
        raise ValueError("Final time must follow the complete detection search")
    tol=TransitionTolerances(**{k:getattr(a,"guard_"+k) for k in asdict(TransitionTolerances())})
    if any(not np.isfinite(x) or x<=0 for x in asdict(tol).values()):raise ValueError("Positive finite guards required")
    from sandbox.stage4_core.parallel_runtime import linux_resources,resource_budget
    resources=linux_resources();budget=resource_budget(resources,.9,.9,max_workers=1)
    groups=[g for g in resources["physical_core_groups"] if set(g)<=set(budget["selected_affinity_cpus"])]
    selected=groups[a.cpu_start:a.cpu_start+a.physical_cores]
    cpus=[c for g in selected for c in g]
    if len(selected)!=a.physical_cores or a.threads>len(cpus):raise ValueError("Threads/whole cores exceed admitted budget")
    os.sched_setaffinity(0,set(cpus))
    out=B.scratch_path(a.output_dir);pins=source_pins();ident=identity(a)
    if a.resume:
        manifest=json.loads((out/"manifest.json").read_text())
        if manifest["run_complete"]:raise ValueError("Completed run is preserved")
        if manifest["source_pins_sha256"]!=pins or manifest["parameter_identity"]!=ident:
            raise ValueError("Exact sources/parameters changed; resume refused")
        if B.sha256(out/"dc1d.npz")!=manifest["dc1d_sha256"]:
            raise ValueError("Saved 1D DC reference changed")
        if B.sha256(out/"checkpoint_latest.npz")!=manifest["checkpoint_sha256"]:
            raise ValueError("Committed checkpoint changed")
        with np.load(out/"checkpoint_latest.npz",allow_pickle=False) as z:
            m=json.loads(str(z["metadata_json"].item()))
            state=EffectiveState(*(z[key].copy() for key in ("delta_bar","p","n","circuit_SI")))
            shadow=None if not m["shadow_active"] else EffectiveState(*(
                z["shadow_"+key].copy() for key in ("delta_bar","p","n","circuit_SI")))
        if m["source_pins_sha256"]!=pins or m["parameter_identity"]!=ident:
            raise ValueError("Checkpoint sources/parameters differ")
        if m.get("terminal_diagnostic"):raise ValueError("Diagnostic terminal state is preserved; use a fresh case")
    else:
        if out.exists():raise FileExistsError("Existing data preserved")
        out.mkdir(parents=True);m=None;shadow=None
    requested={"stop":False}
    def request(sig,frame):requested["stop"]=True
    signal.signal(signal.SIGTERM,request);signal.signal(signal.SIGINT,request)
    started=time.monotonic()
    manifest=dict(schema="pysnspd.dimensional_run.v1",model_version="1.2.0",
        physical_kernel_version="1.1.3",run_complete=False,status="PREPARATION",
        source_pins_sha256=pins,parameter_identity=ident,
        resources=dict(threads=a.threads,physical_cores=a.physical_cores,affinity_cpus=cpus,
                       runtime_inventory=resources),
        admission=dict(implemented=True,physical_validation=False),
        scope="Same finite strip, dy=0 after admitted overlap; no 1D tails or amplifier",
        snapshot_policy="2D fixed targets; 1D accepted states at target cadence without capping dt; actual timestamps stored")
    write_json(out/"manifest.json",manifest)
    setup_args=argparse.Namespace(**dict(vars(a),duration_ps=None))
    s2,dc,g2,material,params=OLD.build_system(setup_args,out/"setup2d",resume=a.resume);observed(s2)
    dc2=EffectiveState(dc.delta_bar.copy(),dc.p.copy(),
        np.broadcast_to(s2.kinetics.bath_phonons,(g2.graph.n_nodes,s2.kinetics.phonon_size)).copy(),
        np.array([a.current_uA*1e-6,a.current_uA*1e-6,0.]))
    s1,g1,dc1=build_1d(a,s2,g2,out,a.resume)
    project=DimensionalProjection(s2,g2,dc2,s1,g1,dc1)
    reference=float(s2.kinetic_residual_energy_scale_J)
    if not a.resume:
        state=dc2.copy()
        if not a.dc_only:
            state.p,state.n,temp,target,source=OLD.prepare_photon(s2.functional,s2.kinetics,
                state.delta_bar,dc2.p,g2,wavelength_nm=a.wavelength_nm,sigma_nm=a.sigma_nm,
                retained_fraction=a.retained_fraction,bath_temperature_K=a.bath_K)
            B.atomic_npz(out/"source_initial_state.npz",dict(delta_bar=state.delta_bar,p=state.p,
                n=state.n,temperature_K=temp,node_deposition_J=target,coordinates_m=g2.coordinates_m))
        else:source=dict(retained_energy_J=0.,mode="DC",source_must_not_be_applied_again=True)
        controller=B.TriggerWindow(B.TriggerWindowConfig(event_time_s=0.,baseline_V=0.,
            baseline_reference="Fresh current-carrying superconducting DC",threshold_V=a.threshold_mV*1e-3,
            polarity="positive",confirmation_s=a.confirmation_ps*1e-12,hysteresis_fraction=.05,
            post_trigger_margin_s=a.extra_time_ps*1e-12,no_trigger_time_ceiling_s=a.ceiling_ps*1e-12))
        controller.observe_accepted(0.,0.)
        m=dict(time_ps=0.,phase="2d",initial_energy_J=s2.evaluate(state).energy_J,
            external_J=0.,absolute_external_J=0.,cut_absolute_J=0.,device_absolute_J=0.,max_local_energy_fraction=0.,max_global_energy_fraction=0.,
            accepted_steps=0,rejections=0,consecutive_rejections=0,stable_steps=0,next_step_ps=a.dt_ps,next_snapshot_ps=a.snapshot_every_ps,next_checkpoint_ps=a.checkpoint_every_ps,
            next_shadow_check_ps=0.,shadow_time_ps=None,shadow_crossing_ps=None,
            shadow_next_step_ps=a.dt_ps,shadow_checks=[],shadow_error=None,projection=None,
            source=source,terminal_diagnostic=False)
    else:controller=B.TriggerWindow.from_checkpoint(m["controller"])
    systems={"2d":s2,"shadow1d":s1,"1d":s1}
    geoms={"2d":g2,"shadow1d":g1,"1d":g1};dcs={"2d":dc2,"shadow1d":dc1,"1d":dc1}
    archive=PhaseArchive(out/"data.h5",systems,geoms,dcs,resume=a.resume,counts=m.get("archive_counts"))
    system=systems[m["phase"]];rhs=system.evaluate(state)
    initial_total=m["initial_energy_J"]
    # Energy references change at handoff, never count a mesh reference shift as heat.
    offset=m.get("energy_reference_offset_J",0.)
    def domain_fraction():return m["cut_absolute_J"]/(reference+m["device_absolute_J"])
    def row(rr,st,t,h,receipt,shadow_row=False):
        Ib,Is,vc=st.circuit_SI
        return dict(time_ps=t,Vout_V=system.circuit.R_load_ohm*(Ib-Is),Vdev_V=rr.electrical.Vdev_V,
            Ib_A=Ib,Is_A=Is,vc_V=vc,energy_J=rr.energy_J,
            energy_defect_J=(np.nan if shadow_row else rr.energy_J+offset-initial_total-m["external_J"]),
            external_J=m["external_J"],absolute_external_J=m["absolute_external_J"],
            cut_absolute_J=m["cut_absolute_J"],device_absolute_J=m["device_absolute_J"],
            domain_fraction=domain_fraction(),escape_power_W=rr.escape_power_W,
            cuts_power_W=rr.reservoir_power_out_W,step_ps=h,
            iterations=receipt.get("iterations",0),temporal_error=receipt.get("temporal_error",0.))
    if not a.resume:
        archive.append("2d",row(rhs,state,0.,0.,{}));archive.snapshot("2d",s2,state,rhs,0.)
    status="RUNNING";error=None;last_checkpoint=time.monotonic();last_progress=last_checkpoint
    def save():
        nonlocal last_checkpoint
        archive.file.flush()
        m.update(controller=controller.checkpoint(),shadow_active=shadow is not None,
            archive_counts=archive.counts(),source_pins_sha256=pins,parameter_identity=ident,
            energy_reference_offset_J=offset)
        arrays={key:getattr(state,key) for key in ("delta_bar","p","n","circuit_SI")}
        if shadow is not None:
            arrays.update({"shadow_"+key:getattr(shadow,key) for key in ("delta_bar","p","n","circuit_SI")})
        arrays["metadata_json"]=np.array(json.dumps(B.safe_json(m),allow_nan=False))
        B.atomic_npz(out/"checkpoint_latest.npz",arrays)
        manifest.update(status=status,error=error,time_ps=m["time_ps"],phase=m["phase"],
            physicalparams=params,source=m["source"],observation=controller.summary(),
            domain_fraction=domain_fraction(),projection=m["projection"],
            shadow_checks=m["shadow_checks"],shadow_error=m["shadow_error"],
            archive_counts=m["archive_counts"],checkpoint_sha256=B.sha256(out/"checkpoint_latest.npz"),
            accepted_steps=m["accepted_steps"],rejections=m["rejections"],
            energy_defect_J=rhs.energy_J+offset-initial_total-m["external_J"],max_local_energy_fraction=m["max_local_energy_fraction"],max_global_energy_fraction=m["max_global_energy_fraction"],
            elapsed_session_wall_s=time.monotonic()-started,
            dc1d_sha256=B.sha256(out/"dc1d.npz"))
        write_json(out/"manifest.json",manifest);last_checkpoint=time.monotonic()
    save()
    try:
        while True:
            t=m["time_ps"]
            if requested["stop"]:status="INTERRUPTED_CHECKPOINTED";break
            if a.duration_ps is not None and t>=a.duration_ps-1e-10:
                status="DEVELOPMENT_DURATION_REACHED";break
            if m["max_global_energy_fraction"]>a.global_energy_tolerance:
                status="INDEPENDENT_ENERGY_AUDIT_FAILED";m["terminal_diagnostic"]=True;break
            if domain_fraction()>tol.cut_energy_relative:
                status="POSSIBLE_DOMAIN_ERROR";m["terminal_diagnostic"]=True;break
            if m["phase"]=="2d" and controller.terminal:
                if controller.status!="TRIGGER_MARGIN_COMPLETE":
                    status=controller.status;manifest["run_complete"]=(status=="TIME_CENSORED_NO_TRIGGER");break
                try:
                    if m["shadow_error"] is not None:raise ArithmeticError(m["shadow_error"])
                    if shadow is None or len(m["shadow_checks"])<2:
                        raise ArithmeticError("Insufficient independent overlap")
                    reduced,report=project.project(state,reference+m["device_absolute_J"],tol)
                    report["time_ps"]=t;m["projection"]=report
                    checks=m["shadow_checks"]+[report]
                    for c in checks:
                        if c["winding_count"] or c["gap_relative"]>tol.gap_relative or c["populations_energy_relative"]>tol.populations_energy_relative or c["current_pattern_relative"]>tol.current_pattern_relative:
                            raise ArithmeticError("Unresolved transverse/projection structure: "+str(c))
                        if c.get("voltage_error_V",0)>tol.shadow_voltage_V or c.get("current_relative",0)>tol.shadow_current_relative or c["film_energy_relative"]>tol.shadow_energy_relative:
                            raise ArithmeticError("Independent 1D overlap mismatch: "+str(c))
                    offset+=film_energy(s2,state)-film_energy(s1,reduced)
                    state=reduced;shadow=None;m["phase"]="1d";system=s1;rhs=s1.evaluate(state)
                    m["next_step_ps"]=a.dt_ps;m["next_snapshot_ps"]=t+a.snapshot_1d_every_ps;m["next_checkpoint_ps"]=t+a.checkpoint_1d_every_ps
                    archive.append("1d",row(rhs,state,t,0.,{}));archive.snapshot("1d",s1,state,rhs,t)
                    save();print(json.dumps(dict(status="HANDOFF_ADMITTED",time_ps=t,projection=report)),flush=True)
                except (ArithmeticError,ValueError) as exc:
                    status="POSSIBLE_DIMENSIONAL_ERROR";error=str(exc);m["terminal_diagnostic"]=True;break
            if m["phase"]=="1d" and t>=a.final_time_ns*1e3-1e-9:
                status="FINAL_TIME_REACHED";manifest["run_complete"]=True;break
            horizon=(controller.next_required_time_s()*1e12 if m["phase"]=="2d" else a.final_time_ns*1e3)
            if a.duration_ps is not None:horizon=min(horizon,a.duration_ps)
            h=min(m["next_step_ps"],horizon-t)
            if m["phase"]=="2d":h=min(h,m["next_snapshot_ps"]-t)
            if h<=0 or t+h==t:raise ArithmeticError("No representable positive physical step")
            kwargs=dict(tolerance=a.step_tolerance,max_iterations=a.max_iterations,
                energy_tolerance=a.local_energy_tolerance,reference=reference,
                absolute_work=m["absolute_external_J"])
            try:
                if m["phase"]=="2d":new,newrhs,receipt=checked_step(system,state,rhs,h,**kwargs)
                else:new,newrhs,receipt=doubled_step(system,state,rhs,h,
                    time_rtol=a.time_rtol,time_atol=a.time_atol,**kwargs)
            except (ArithmeticError,ValueError,FloatingPointError) as exc:
                m["rejections"]+=1;m["consecutive_rejections"]+=1;m["stable_steps"]=0;m["next_step_ps"]=h/2
                if m["consecutive_rejections"]>=a.max_rejections:
                    raise ArithmeticError("Consecutive rejection budget exhausted: "+str(exc)) from exc
                if h/2<a.minimum_adaptive_step_ps:
                    raise ArithmeticError("Declared adaptive resolution reached: "+str(exc)) from exc
                continue
            # Commit only complete, guarded physical increments.
            state=new;rhs=newrhs;m["time_ps"]=t+h;m["accepted_steps"]+=1;m["consecutive_rejections"]=0;m["stable_steps"]+=1
            m["external_J"]+=receipt["external_energy_increment_J"]
            m["absolute_external_J"]+=receipt["absolute_external_energy_increment_J"]
            m["max_local_energy_fraction"]=max(m["max_local_energy_fraction"],receipt["accepted_local_energy_fraction"])
            m["max_global_energy_fraction"]=max(m["max_global_energy_fraction"],
                abs(rhs.energy_J+offset-initial_total-m["external_J"])/(reference+m["absolute_external_J"]))
            m["cut_absolute_J"]+=receipt["absolute_cut_energy_increment_J"]
            m["device_absolute_J"]+=receipt["absolute_device_work_increment_J"]
            t=m["time_ps"]
            archive.append(m["phase"],row(rhs,state,t,h,receipt))
            cap=a.dt_ps if m["phase"]=="2d" else a.maximum_1d_step_ps
            factor=((2. if m["stable_steps"]>=8 else 1.) if m["phase"]=="2d" else min(2.,max(1.,.8/max(receipt["temporal_error"],1e-12)**(1/3))))
            if m["phase"]=="2d" and m["stable_steps"]>=8:m["stable_steps"]=0
            m["next_step_ps"]=min(cap,h*factor)
            if m["phase"]=="2d":
                controller.observe_accepted(t*1e-12,float(s2.circuit.R_load_ohm*(state.circuit_SI[0]-state.circuit_SI[1])))
                cross=(controller.crossing_time_s if controller.crossing_time_s is not None else controller.candidate_crossing_s)
                if cross is None:
                    shadow=None;m["shadow_crossing_ps"]=None;m["shadow_error"]=None;m["shadow_checks"]=[]
                elif m["shadow_error"] is None:
                    try:
                        crossing_ps=cross*1e12
                        scale=reference+m["device_absolute_J"]
                        if shadow is None or m["shadow_crossing_ps"]!=crossing_ps:
                            shadow,rep=project.project(state,scale,tol)
                            m["shadow_time_ps"]=t;m["shadow_crossing_ps"]=crossing_ps
                            m["shadow_start_lag_ps"]=t-crossing_ps;m["next_shadow_check_ps"]=t
                        if t>=m["next_shadow_check_ps"]-1e-10 or controller.terminal:
                            srhs=s1.evaluate(shadow)
                            while m["shadow_time_ps"]<t-1e-10:
                                sh=min(m["shadow_next_step_ps"],t-m["shadow_time_ps"])
                                try:
                                    candidate,cr,rec=doubled_step(s1,shadow,srhs,sh,
                                        time_rtol=a.time_rtol,time_atol=a.time_atol,**kwargs)
                                except (ArithmeticError,ValueError):
                                    m["shadow_next_step_ps"]=sh/2
                                    if sh/2<a.minimum_adaptive_step_ps:raise ArithmeticError("Shadow resolution exhausted")
                                    continue
                                shadow=candidate;srhs=cr;m["shadow_time_ps"]+=sh
                                m["shadow_next_step_ps"]=min(a.dt_ps,sh*1.4)
                            check=project.errors(state,shadow,scale)
                            check.update(time_ps=t,
                                voltage_error_V=abs(s2.circuit.R_load_ohm*((state.circuit_SI[0]-state.circuit_SI[1])-
                                    (shadow.circuit_SI[0]-shadow.circuit_SI[1]))),
                                current_relative=float(np.max(abs(state.circuit_SI[:2]-shadow.circuit_SI[:2]))/(a.current_uA*1e-6)))
                            m["shadow_checks"].append(check)
                            archive.append("shadow1d",row(srhs,shadow,t,0.,{},True))
                            archive.snapshot("shadow1d",s1,shadow,srhs,t)
                            m["next_shadow_check_ps"]=t+a.shadow_check_every_ps
                    except (ArithmeticError,ValueError) as exc:m["shadow_error"]=str(exc)
            if t>=m["next_snapshot_ps"]-1e-10:
                archive.snapshot(m["phase"],system,state,rhs,t)
                every=a.snapshot_every_ps if m["phase"]=="2d" else a.snapshot_1d_every_ps
                m["next_snapshot_ps"]=t+every
            if time.monotonic()-last_checkpoint>=a.checkpoint_wall_s or t>=m["next_checkpoint_ps"]-1e-10:
                m["next_checkpoint_ps"]=t+(a.checkpoint_every_ps if m["phase"]=="2d" else a.checkpoint_1d_every_ps)
                save()
            if time.monotonic()-last_progress>=30:
                print(json.dumps(dict(status="RUNNING",phase=m["phase"],time_ps=t,
                    accepted_steps=m["accepted_steps"],step_ps=h,
                    elapsed_session_wall_s=time.monotonic()-started,Vout_mV=1e3*system.circuit.R_load_ohm*(state.circuit_SI[0]-state.circuit_SI[1]),
                    domain_fraction=domain_fraction())),flush=True);last_progress=time.monotonic()
    except Exception as exc:
        error=type(exc).__name__+": "+str(exc)
        status="NUMERICAL_RESOLUTION_LIMIT" if "Declared adaptive resolution" in str(exc) else "INTEGRATION_FAILURE"
        m["terminal_diagnostic"]=True
    finally:
        archive.snapshot(m["phase"],system,state,rhs,m["time_ps"]);save();archive.close()
    print(json.dumps(dict(status=status,time_ps=m["time_ps"],phase=m["phase"],error=error)),flush=True)
    return 0 if status in ("FINAL_TIME_REACHED","TIME_CENSORED_NO_TRIGGER","DEVELOPMENT_DURATION_REACHED","INTERRUPTED_CHECKPOINTED") else 2


if __name__=="__main__":
    arguments=parser().parse_args()
    try:code=run(arguments)
    except Exception as exc:
        path=Path(arguments.output_dir)/"manifest.json"
        if path.is_file():
            prior=json.loads(path.read_text())
            if prior.get("status")=="PREPARATION":
                prior.update(status="PREPARATION_FAILURE",error=type(exc).__name__+": "+str(exc))
                write_json(path,prior)
        raise
    raise SystemExit(code)
