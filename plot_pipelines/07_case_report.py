"""One independent scientific PDF per effective SNSPD run, from saved data only.

No time evolution, kinetic RHS, fresh DC solve, fitting or inter-run comparison.
Completed inputs are authenticated. Optional live reports capture one atomic
checkpoint and remain explicitly provisional. All products stay in scratch.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
import math
import multiprocessing as mp
import os
from pathlib import Path
import re
import shutil
import sys
import textwrap
import time
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from sandbox.stage4_core.parallel_runtime import limit_thread_environment
limit_thread_environment()
os.environ.setdefault("MPLCONFIGDIR", "/home/jdiaz/scratch/tmp/pysnspd/matplotlib")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.tri import Triangulation
from scipy.constants import elementary_charge as e, hbar
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from pysnspd.plotting.style import apply_thesis_style, THESIS_DPI
from pysnspd.plotting.complex_gap import draw_complex_gap
from pysnspd.experimental.effective_bcs_functional import EffectiveBCSFunctional, bcs_count_quadrature
from pysnspd.experimental.effective_bcs_ports import EffectiveNeutralPorts
from pysnspd.experimental.thermal_spatial_usadel import ThermalGraph
from pysnspd.experimental.trigger_window import TriggerWindow

SCRATCH = Path("/home/jdiaz/scratch")
EV = 1.602176634e-19
MEV = EV*1e-3
SCHEMA = "pysnspd.individual_effective_case_report.v2"
RECON_SOURCES = ("pysnspd/experimental/effective_bcs_functional.py",
    "pysnspd/experimental/effective_bcs_ports.py","pysnspd/experimental/thermal_spatial_usadel.py")
THERMAL_CANDIDATE_VERSIONS = ("1.1.0-candidate","1.1.1-candidate","1.1.2-candidate","1.1.3-candidate")
DC_MATCHED_THERMAL_VERSION = "1.1.3-candidate"
TRACE = ("time_s","Vout_V","Vdev_V","Ib_A","Is_A","vc_V","energy_total_J",
    "energy_defect_J","escape_power_W","domain_energy_J","instant_power_residual_W",
    "reservoir_power_out_W","external_energy_integral_J","accepted_step_ps")
MAPS = ("time_s","delta_J","Te_K","Tph_K","u_excess_J_m3","electron_energy_J_m3","phonon_energy_J_m3")
INK, ACCENT = "#143448", "#16738b"


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda:stream.read(2**20),b""): h.update(block)
    return h.hexdigest()


def scratch_path(value):
    p=Path(value).expanduser().resolve()
    if not p.is_relative_to(SCRATCH.resolve()):
        raise ValueError("Products and captures must remain under /home/jdiaz/scratch")
    return p


def fmt(value,scale=1.,digits=6):
    if value is None: return "Unassigned"
    v=float(value)*scale
    return f"{v:.{digits}g}" if np.isfinite(v) else "infinity"


def case_id(path):
    p=Path(path).resolve()
    return re.sub(r"[^A-Za-z0-9_.-]+","_",p.name)+"_"+hashlib.sha256(str(p).encode()).hexdigest()[:8]


def memory_check(paths,limit):
    size=0
    for path,names in paths:
        with zipfile.ZipFile(path) as z:
            sizes={i.filename:i.file_size for i in z.infolist()}
            size+=sum(sizes[n+".npy"] for n in names)
    required=5*size+512*1024**2
    if required>limit: raise MemoryError("Decompressed report inputs exceed per-case RAM budget")
    return {"requested_uncompressed_bytes":size,"conservative_peak_estimate_bytes":required}


def existing_audit(run,manifest_sha):
    path=run.parent.parent/"campaign_execution.json"
    if not path.is_file(): return {"status":"No authenticated terminal audit"}
    record=json.loads(path.read_text()).get("cases",{}).get(run.name,{})
    audit=record.get("audit",{}); p=Path(audit.get("path",""))
    if record.get("manifest_sha256")!=manifest_sha or not p.is_file() or sha256(p)!=audit.get("sha256"):
        return {"status":"Missing or stale terminal audit"}
    content=json.loads(p.read_text())
    return dict(status="PASS" if content.get("all_requested_checks_passed") is True else "NO PASS",
        path=str(p),sha256=audit["sha256"],existing_audit_only=True)


def load_case(run,build,allow_incomplete,limit):
    run=Path(run).resolve(); raw=(run/"manifest.json").read_bytes(); m=json.loads(raw)
    if m.get("schema")!="pysnspd.effective_photon_run.v1" or m.get("closure",{}).get("microscopic_exact_Usadel") is not False:
        raise ValueError("Explicit effective pipeline05 case required")
    if m["parameter_identity"].get("dc_only"):
        raise ValueError("Photon case required: photon-free DC holds need a different clock and source report")
    inputs={"manifest.json":hashlib.sha256(raw).hexdigest()}
    dcfile=run/"dc/effective_dc.npz"; drfile=run/"dc/effective_dc_receipt.json"
    dr=json.loads(drfile.read_text())
    if sha256(dcfile)!=dr["fields_sha256"]: raise ValueError("DC checksum differs")
    inputs.update({"dc/effective_dc.npz":sha256(dcfile),"dc/effective_dc_receipt.json":sha256(drfile)})
    with np.load(dcfile,allow_pickle=False) as z:
        dc={k:z[k].copy() for k in ("delta_bar","coordinates_m","left_nodes","right_nodes")}
    complete=m.get("run_complete") is True
    if complete:
        if m["status"] not in {"TRIGGER_MARGIN_COMPLETE","TIME_CENSORED_NO_TRIGGER","DC_HOLD_COMPLETE","DEVELOPMENT_WINDOW_COMPLETE"}:
            raise ValueError("Successful finite-window terminal required")
        for name in ("trajectory.npz","snapshots.npz","dc/effective_dc.npz"):
            if sha256(run/name)!=m.get("output_sha256",{}).get(name): raise ValueError("Terminal checksum differs: "+name)
            inputs[name]=m["output_sha256"][name]
        memory=memory_check([(run/"trajectory.npz",TRACE),(run/"snapshots.npz",MAPS)],limit)
        with np.load(run/"trajectory.npz",allow_pickle=False) as z: trace={k:z[k].copy() for k in TRACE}
        with np.load(run/"snapshots.npz",allow_pickle=False) as z:
            maps={k:z[k].copy() for k in (*MAPS,"coordinates_m","triangles","edges")}
        obs=m["observation"].copy(); execution=m.get("execution",{})
    else:
        if not allow_incomplete: raise ValueError("Incomplete case: add --allow-incomplete for a provisional checkpoint report")
        capture=build/"checkpoint_capture.npz"
        # Read one open inode, coherent even if the running solver replaces it.
        with (run/"checkpoint_latest.npz").open("rb") as src,capture.open("xb") as dst:
            shutil.copyfileobj(src,dst,2**20)
        inputs["checkpoint_capture.npz"]=sha256(capture)
        memory=memory_check([(capture,tuple("trajectory__"+k for k in TRACE)+tuple("snapshots__"+k for k in MAPS))],limit)
        with np.load(capture,allow_pickle=False) as z:
            meta=json.loads(str(z["metadata_json"]))
            if meta["parameter_identity"]!=m["parameter_identity"] or meta["source_pins_sha256"]!=m["source_pins_sha256"]:
                raise ValueError("Checkpoint and manifest identities differ")
            trace={k:z["trajectory__"+k].copy() for k in TRACE}
            maps={k:z["snapshots__"+k].copy() for k in MAPS}
        if not math.isclose(trace["time_s"][-1]*1e12,meta["time_ps"],rel_tol=1e-12,abs_tol=1e-9):
            raise ValueError("Checkpoint time disagrees with captured trajectory")
        obs=TriggerWindow.from_checkpoint(meta["controller"]).summary() if meta.get("controller") else {}
        obs["trigger_confirmed"]=obs.get("detected",False)
        obs.setdefault("baseline_Vout_V",m["observation"].get("baseline_Vout_V",0.))
        m["source"]=meta["source"]; execution=meta
    candidates=[run/"meshes/mesh.npz",Path(m["geometry"]["source_path"])]
    meshfile=next((p for p in candidates if p.is_file() and sha256(p)==m["geometry"]["mesh_sha256"]),None)
    if meshfile is None: raise FileNotFoundError("Pinned original geometry is required; no remeshing")
    inputs["geometry_mesh"]=sha256(meshfile)
    with np.load(meshfile,allow_pickle=False) as z:
        mesh={k:z[k].copy() for k in ("coordinates_m","coordinates_bar","triangles","edges","area_weights",
            "conductance","dual_area_m2","edge_length_m","dual_face_length_m","fixed_nodes")}
    for k in ("coordinates_m","triangles","edges"):
        maps.setdefault(k,mesh[k])
        if not np.array_equal(maps[k],mesh[k]): raise ValueError("Saved geometry differs: "+k)
    if not np.array_equal(dc["coordinates_m"],maps["coordinates_m"]): raise ValueError("DC geometry differs")
    for family in (trace,maps):
        t=family["time_s"]
        if t.ndim!=1 or len(t)<1 or np.any(~np.isfinite(t)) or np.any(np.diff(t)<=0):
            raise ValueError("Finite increasing saved times required")
    for k in TRACE:
        if trace[k].shape!=trace["time_s"].shape or np.any(~np.isfinite(trace[k])): raise ValueError("Malformed trace: "+k)
    for k in MAPS[1:]:
        if maps[k].shape!=(len(maps["time_s"]),len(mesh["coordinates_m"])) or np.any(~np.isfinite(maps[k])):
            raise ValueError("Malformed snapshot: "+k)
    if maps["time_s"][0]<trace["time_s"][0]-1e-24 or maps["time_s"][-1]>trace["time_s"][-1]+1e-24:
        raise ValueError("Snapshots outside trajectory")
    if np.any(maps["Te_K"]<0) or np.any(maps["Tph_K"]<0): raise ValueError("Negative equivalent temperature")
    if complete:
        for k in ("trajectory.npz","snapshots.npz"):
            if sha256(run/k)!=inputs[k]: raise ValueError("Input changed during reading")
    for k in reconstruction_sources(m):
        if sha256(ROOT/k)!=m["source_pins_sha256"].get(k): raise ValueError("Saved reconstruction source differs: "+k)
    return dict(run=run,m=m,z=trace,s=maps,mesh=mesh,dc=dc,dc_receipt=dr,obs=obs,provisional=not complete,
        inputs=inputs,memory=memory,execution=execution,audit=existing_audit(run,inputs["manifest.json"]) if complete else {"status":"Provisional; no terminal audit inferred"})


def hotbelts(delta,dc,edges,xy,threshold=.2):
    """Saved-node graph diagnostic; no inference between snapshots."""
    if np.any(abs(dc)<=0): raise ValueError("Nonzero DC gap required")
    y=xy[:,1]; tol=1e-12*max(float(np.ptp(y)),1e-9)
    bottom=np.flatnonzero(abs(y-y.min())<=tol); top=np.flatnonzero(abs(y-y.max())<=tol)
    flags=[]
    for d in delta:
        low=abs(d)<=threshold*abs(dc)
        if not np.any(low[bottom]) or not np.any(low[top]): flags.append(False); continue
        face=edges[low[edges[:,0]] & low[edges[:,1]]]
        graph=coo_matrix((np.ones(2*len(face)),(np.r_[face[:,0],face[:,1]],np.r_[face[:,1],face[:,0]])),shape=(len(y),len(y))).tocsr()
        _,label=connected_components(graph,directed=False)
        flags.append(bool(set(label[bottom[low[bottom]]]) & set(label[top[low[top]]])))
    return np.asarray(flags,bool)



def snapshot_targets(times,obs,hotbelt,requested=None):
    """Exactly six ordered target slots; explicit requests never gain endpoints."""
    t=np.asarray(times,float)*1e12
    if t.ndim!=1 or not len(t) or np.any(~np.isfinite(t)) or np.any(np.diff(t)<=0):
        raise ValueError("Finite, increasing saved snapshot times required")
    if requested is not None:
        targets=np.asarray(list(requested)[:6],float)
        if not len(targets) or np.any(~np.isfinite(targets)):
            raise ValueError("At least one finite requested snapshot time is required")
    elif obs.get("trigger_confirmed",obs.get("detected",False)):
        cross=obs.get("first_crossing_time_s",obs.get("crossing_time_s"))
        conf=obs.get("confirmation_time_s")
        hb=t[np.flatnonzero(hotbelt)[0]] if np.any(hotbelt) else .5*t[-1]
        targets=np.array([t[0],min(1.,t[-1]),hb,cross*1e12,conf*1e12,t[-1]])
    else:
        targets=np.array([t[0],min(1.,t[-1]),min(10.,t[-1]),.25*t[-1],.5*t[-1],t[-1]])
    return np.pad(targets,(0,6-len(targets)),mode="edge")


def select_snapshots(times,obs,hotbelt,requested=None):
    targets=snapshot_targets(times,obs,hotbelt,requested)
    saved=np.asarray(times)*1e12
    # Preserve request order and repetitions; ties choose the earlier saved sample.
    return np.array([int(np.argmin(abs(saved-q))) for q in targets],int)


def snapshot_selection(times,targets,indices):
    t=np.asarray(times)*1e12;dt=np.diff(t)
    slots=[dict(slot=j+1,requested_ps=float(q),saved_ps=float(t[i]),
        offset_ps=float(t[i]-q),saved_index=int(i),repeated_saved_state=int(i) in indices[:j],
        outside_saved_window=bool(q<t[0] or q>t[-1])) for j,(q,i) in enumerate(zip(targets,indices))]
    queue=[]
    if len(t)<6:
        queue.append("Future acquisition: save at least six states with a stable physical-time cadence.")
    if len(np.unique(indices))<len(np.unique(np.clip(targets,t[0],t[-1]))):
        queue.append("Future acquisition: use a stable snapshot cadence finer than the requested time spacing.")
    return dict(policy="first six requests; pad with final request; nearest saved time; no extra endpoints",
        slots=slots,available_states=len(t),unique_selected_states=len(np.unique(indices)),
        saved_cadence_ps=dict(minimum=None if not len(dt) else float(dt.min()),
            median=None if not len(dt) else float(np.median(dt)),maximum=None if not len(dt) else float(dt.max())),
        future_visualization_queue=queue)


def phase_vorticity(delta,dc,xy,triangles,relative_floor=1e-6):
    """Signed CCW winding of complex P1 triangles with nonzero contour.

    Uses the complex order parameter, not curl(J). The contour floor is a
    diagnostic amplitude resolution relative to the local DC, not a model input.
    Cores are zeros of the same P1 interpolant. No temporal tracking is inferred.
    """
    d=np.asarray(delta,complex);dc=np.asarray(dc,complex);xy=np.asarray(xy,float)
    tri=np.asarray(triangles,int);z=d[tri];ref=np.min(abs(dc[tri]),axis=1)
    if not 0<relative_floor<1 or np.any(ref<=0) or np.any(~np.isfinite(z)):
        raise ValueError("Finite gap, nonzero DC reference and 0<floor<1 required")
    points=xy[tri];a=points[:,1]-points[:,0];b=points[:,2]-points[:,0]
    area2=a[:,0]*b[:,1]-a[:,1]*b[:,0]
    if np.any(area2==0):raise ValueError("Nondegenerate physical triangles required")
    # Normalize before phase products and zero finding to avoid SI underflow.
    scale=np.maximum(np.max(abs(z),axis=1),ref)
    v=z/scale[:,None];other=np.roll(v,-1,axis=1);edge=other-v
    length2=abs(edge)**2
    frac=np.divide(-np.real(np.conj(edge)*v),length2,out=np.zeros_like(length2),where=length2>0)
    closest=v+np.clip(frac,0,1)*edge
    minimum_contour_gap=np.min(abs(closest),axis=1)*scale
    valid=minimum_contour_gap>relative_floor*ref
    raw=np.angle(other*np.conj(v)).sum(axis=1)/(2*np.pi)*np.sign(area2)
    nearest=np.rint(raw).astype(int)
    valid &= (abs(raw-nearest)<1e-8)&(abs(nearest)<=1)
    charge=np.where(valid,nearest,0)
    candidates=np.flatnonzero(valid&(charge!=0));cores=[];signs=[];cells=[]
    for j in candidates:
        matrix=np.array([[np.real(v[j,1]-v[j,0]),np.real(v[j,2]-v[j,0])],
                         [np.imag(v[j,1]-v[j,0]),np.imag(v[j,2]-v[j,0])]])
        try:uv=np.linalg.solve(matrix,-np.array([v[j,0].real,v[j,0].imag]))
        except np.linalg.LinAlgError:
            valid[j]=False;charge[j]=0;continue
        weights=np.r_[1-uv.sum(),uv]
        if np.any(weights<-1e-9) or abs(weights@v[j])>1e-9:
            valid[j]=False;charge[j]=0;continue
        cores.append(weights@points[j]);signs.append(int(charge[j]));cells.append(int(j))
    return dict(winding=charge,valid=valid,core_points_m=np.asarray(cores).reshape(-1,2),
        core_charge=np.asarray(signs,int),core_cells=np.asarray(cells,int),
        contour_floor_relative_DC=relative_floor,positive=int(np.sum(charge[valid]==1)),
        negative=int(np.sum(charge[valid]==-1)),indeterminate_cells=int(np.sum(~valid)),
        resolved_cells=int(np.sum(valid)),minimum_contour_gap_J=minimum_contour_gap,
        convention="positive = counterclockwise phase winding in physical x,y; A=0",
        reference="https://www.mcs.anl.gov/papers/P5324-0315.pdf")


def metrics(c):
    z,m,o=c["z"],c["m"],c["obs"]; t=z["time_s"]*1e12
    baseline=float(o.get("baseline_Vout_V",o.get("config",{}).get("baseline_V",0.)))
    voltage=(z["Vout_V"]-baseline)*1e3
    load=m["physicalparams"]["circuit"]["R_load_ohm"]
    err=float(np.max(abs(z["Vout_V"]-load*(z["Ib_A"]-z["Is_A"]))))
    if err>1e-12+1e-10*max(float(abs(z["Vout_V"]).max()),1e-12): raise ValueError("Vout/readout-current identity failed")
    detected=o.get("trigger_confirmed",o.get("detected",False)) is True
    cross=o.get("first_crossing_time_s",o.get("crossing_time_s")) if detected else None
    if detected and (cross is None or cross>z["time_s"][-1]+1e-24): raise ValueError("Event outside captured trajectory")
    slope=None
    if cross is not None:
        window=abs(t-cross*1e12)<=.25
        if window.sum()>=3: slope=float(np.polyfit(t[window]-cross*1e12,voltage[window],1)[0])
    return dict(final_ps=float(t[-1]),maps_final_ps=float(c["s"]["time_s"][-1]*1e12),confirmed=detected,
        censored=bool(o.get("censored")) and not c["provisional"],provisional=c["provisional"],
        crossing_ps=None if cross is None else float(cross*1e12),
        confirmation_ps=None if not detected else float(o["confirmation_time_s"]*1e12),
        latency_ps=None if cross is None else float((cross-o.get("config",{}).get("event_time_s",0.))*1e12),
        max_observed_Vout_mV=float(voltage.max()),min_Is_uA=float(z["Is_A"].min()*1e6),
        local_slope_mV_ps=slope,readout_identity_error_V=err,
        positive_global_currents=bool(np.all(z["Ib_A"]>0) and np.all(z["Is_A"]>0)))


def thermal_candidate(manifest):
    mode=manifest.get("closure",{}).get("stiffness_temperature",
        manifest.get("parameter_identity",{}).get("stiffness_temperature"))
    version=manifest.get("model_version","")
    if mode=="electron" and version not in THERMAL_CANDIDATE_VERSIONS:
        raise ValueError("Explicit model_version from supported 1.1.x candidates required for Te reconstruction")
    return mode=="electron" and version in THERMAL_CANDIDATE_VERSIONS


def thermal_depairing_candidate(manifest):
    if not thermal_candidate(manifest) or manifest.get("model_version")!=DC_MATCHED_THERMAL_VERSION:
        return False
    # Explicit identity fields must agree. Early 1.1.3 manifests predate this
    # selectable control and identify the new closure by version alone.
    selections=[family["gradient_closure"] for family in
        (manifest.get("closure",{}),manifest.get("physicalparams",{}),
         manifest.get("parameter_identity",{})) if "gradient_closure" in family]
    if any(value not in ("thermal-depairing","quadratic") for value in selections):
        raise ValueError("Unsupported gradient_closure in saved model identity")
    if len(set(selections))>1:
        raise ValueError("Saved gradient_closure identity fields disagree")
    return not selections or selections[0]=="thermal-depairing"


def reconstruction_sources(manifest):
    if not thermal_candidate(manifest):return RECON_SOURCES
    sources=RECON_SOURCES+("pysnspd/experimental/thermal_stiffness.py",)
    if thermal_depairing_candidate(manifest):
        sources+=("pysnspd/experimental/general_thermal_gradient.py",
            "pysnspd/experimental/finite_depairing_sums.py","pysnspd/experimental/resummed_thermal.py")
    return sources


def model_contract(manifest):
    thermal=thermal_candidate(manifest);depairing=thermal_depairing_candidate(manifest)
    mode=manifest.get("closure",{}).get("stiffness_temperature",
        manifest.get("parameter_identity",{}).get("stiffness_temperature"))
    bath_control=(manifest.get("model_version") in THERMAL_CANDIDATE_VERSIONS and mode=="bath")
    return dict(model_label=("Effective neutral BCS/KWT with DC-matched thermal depairing correction" if depairing else
        "BCS/KWT Te candidate" if thermal else
        "BCS/KWT bath control, calibrated D" if bath_control else "Effective neutral BCS/KWT"),
        candidate=thermal,bath_control=bath_control,thermal_depairing=depairing,
        pauli_signed_heat=(manifest.get("model_version") in THERMAL_CANDIDATE_VERSIONS[1:]),
        stiffness_temperature=("saved Te_K" if thermal else "fixed bath"),
        stored_energy=("U_BCS+G-sum(Te*dG/dTe)+phonons+circuit" if thermal else
            "U_BCS+fixed-bath gradient+phonons+circuit"),
        force_current=("fixed-Te free-gradient derivatives; force=U_d-sum(Te*S_d)" if thermal else
            "fixed-bath internal-energy derivatives"),
        gradient_free_energy=("G_i=kappa(rho_i,Tb)*L_i+u_i*(phi(rho_i,Te_i/Tc,Q_i)-phi(rho_i,Tb/Tc,Q_i)); Q_i=L_i/(u_i*rho_i)"
            if depairing else "quadratic thermal gradient" if thermal else "fixed-bath gradient"),
        native_kinetics_limit=("Native quasiparticle spectra in collisions and transport use q=0; no dynamic pair-breaking DOS"
            if depairing else None),
        source=("condensate/Joule plus signed reversible Te*d(dG/dTe)/dt" if thermal else
            "condensate plus Joule heat"),
        admission=("candidate only; Ceff>0 guard, no capacity floor or global admission" if thermal else
            "individual saved case; no admission inferred"))


def reconstruct_functional(c):
    m,mesh=c["m"],c["mesh"];p=m["physicalparams"]
    graph=ThermalGraph(mesh["area_weights"],mesh["edges"],mesh["conductance"],
        mesh["coordinates_bar"],mesh["fixed_nodes"])
    xi,w=bcs_count_quadrature(m["parameter_identity"]["xi_cutoff"],m["parameter_identity"]["xi_order"])
    cls=EffectiveBCSFunctional
    if thermal_depairing_candidate(m):
        from pysnspd.experimental.resummed_thermal import ResummedThermalFunctional
        cls=ResummedThermalFunctional
    elif thermal_candidate(m):
        from pysnspd.experimental.thermal_stiffness import ThermalStiffnessFunctional
        cls=ThermalStiffnessFunctional
    return cls(graph,energy_scale_J=p["energy_scale_J"],N0_per_J_m3=p["N0_single_spin_per_J_m3"],
        ell0_m=m["geometry"]["ell0_m"],thickness_m=p["film_thickness_m"],count_nodes=xi,count_weights=w,
        gradient_law="bcs_radial",bath_temperature_K=p["bath_temperature_K"],diffusion_m2_s=p["diffusion_m2_s"])


def edge_supercurrent(f,d,temperature_K=None):
    """Exact gauge derivative of the identified closure at saved temperature."""
    closure=getattr(f,"closure_id",None)
    if closure=="TE_DC_matched_finite_depairing_thermal_gradient_v1":
        if temperature_K is None:raise ValueError("Saved Te_K is required for thermal-depairing current reconstruction")
        return f.gradient_evaluation(d,np.asarray(temperature_K,float)).edge_current_A
    tail,head=f.graph.edges.T
    if closure=="TE_free_gradient_internal_entropy_candidate_v1":
        if temperature_K is None:raise ValueError("Saved Te_K is required for thermal-stiffness current reconstruction")
        kap=f.rigidity(abs(d)**2,temperature_K)[0]
    else:kap=f.rigidity(abs(d)**2)[0]
    return 4*e/hbar*f.edge_energy_scale_J*(kap[tail]+kap[head])/2*np.imag(np.conj(d[tail])*d[head])


def triangle_vectors(mesh,flux,thickness):
    """Display-only weighted LS of conservative dual-face projections."""
    edges=mesh["edges"]; xy=mesh["coordinates_m"]
    unit=(xy[edges[:,1]]-xy[edges[:,0]])/mesh["edge_length_m"][:,None]
    areas=mesh["dual_face_length_m"]*thickness
    if np.any(areas<=0): raise ValueError("Positive dual face areas required")
    density=flux/areas
    lookup={tuple(sorted((int(a),int(b)))):i for i,(a,b) in enumerate(edges)}
    ids=np.asarray([[lookup[tuple(sorted((int(a),int(b))))] for a,b in ((t[0],t[1]),(t[1],t[2]),(t[2],t[0]))] for t in mesh["triangles"]])
    u=unit[ids]; weights=areas[ids]/areas[ids].max(axis=1)[:,None]
    matrix=np.einsum("tki,tk,tkj->tij",u,weights,u)
    rhs=np.einsum("tki,tk,tk->ti",u,weights,density[ids])
    vectors=np.linalg.solve(matrix,rhs[...,None])[...,0]
    residual=np.einsum("tki,ti->tk",u,vectors)-density[ids]
    rms=float(np.sqrt(np.sum(weights*residual**2)/weights.sum()))
    scale=float(np.sqrt(np.sum(weights*density[ids]**2)/weights.sum()))
    return vectors,dict(projection_RMS_A_m2=rms,relative_projection_RMS=None if scale==0 else rms/scale)


def currents(c,selected):
    m,mesh,z=c["m"],c["mesh"],c["z"]; p=m["physicalparams"]
    f=reconstruct_functional(c)
    ports=EffectiveNeutralPorts(f.graph,p["sheet_resistance_ohm"],c["dc"]["left_nodes"],c["dc"]["right_nodes"])
    results={}; checks=[]
    for i in selected:
        ts=c["s"]["time_s"][i]
        exact=np.flatnonzero(abs(z["time_s"]-ts)<=1e-23)
        if not len(exact): raise ValueError("No accepted electrical state at snapshot; temporal interpolation refused")
        j=int(exact[np.argmin(abs(z["time_s"][exact]-ts))]); Is=float(z["Is_A"][j])
        temperature=c["s"]["Te_K"][i] if thermal_candidate(m) else None
        el=ports.solve(edge_supercurrent(f,c["s"]["delta_J"][i]/f.Eu_J,temperature),Is)
        err=float(abs(el.Vdev_V-z["Vdev_V"][j]))
        if err>1e-11+1e-8*max(abs(el.Vdev_V),abs(z["Vdev_V"][j]),1e-8):
            raise ArithmeticError("Reconstructed Vdev differs from saved state")
        port=float(np.max(abs(el.port_current_residual_A)))
        if port>1e-12+1e-8*abs(Is): raise ArithmeticError("Cut current reconstruction failed")
        fields={"phi":el.phi_V}; fits={}
        for key,flux in (("super",el.I_super_A),("normal",el.I_normal_A),("total",el.I_total_A)):
            fields[key],fits[key]=triangle_vectors(mesh,flux,p["film_thickness_m"])
        results[int(i)]=fields
        checks.append(dict(snapshot_index=int(i),time_ps=float(ts*1e12),Vdev_error_V=err,
            cut_current_error_A=port,interior_current_error_A=float(np.max(abs(el.interior_current_residual_A),initial=0.)),
            display_vector_fits=fits,model_contract=model_contract(m)))
    return results,checks



class Report:
    def __init__(self,path,c,dpi):
        self.c=c;self.dpi=dpi;self.pages=[]
        self.pdf=PdfPages(path,metadata={"Title":"SNSPD - "+c["run"].name,"Author":"pysnspd",
            "Subject":"Individual effective saved case; no new evolution or inter-run comparison"})
    def page(self,title):
        fig=plt.figure(figsize=(11.69,8.27))
        fig.text(.055,.945,title,fontsize=20,weight="bold",color=INK)
        a=self.c["m"]["parameter_identity"]
        subtitle=f'{self.c["run"].name} | {a["wavelength_nm"]:g} nm | I0={a["current_uA"]:.7g} µA'
        if self.c["provisional"]:subtitle+=" | PROVISIONAL - "+self.c["m"]["status"]
        fig.text(.055,.912,subtitle,fontsize=11,color=ACCENT)
        return fig
    def save(self,fig,title,caption,source_url=None):
        caption_artist=fig.text(.055,.055,textwrap.fill(caption,145),fontsize=9.5,color=INK,va="bottom",linespacing=1.25)
        if source_url:caption_artist.set_url(source_url)
        fig.text(.055,.022,f'pysnspd | individual case | postcascade clock | window 0-{self.c["metrics"]["final_ps"]:.7g} ps',fontsize=8.5,color=INK)
        fig.text(.945,.022,str(len(self.pages)+1),fontsize=9,ha="right",color=INK)
        with plt.rc_context({"savefig.bbox":None}):self.pdf.savefig(fig,dpi=self.dpi)
        plt.close(fig);self.pages.append(dict(page=len(self.pages)+1,title=title,caption=caption))


def table(fig,bounds,title,rows,font=10):
    ax=fig.add_axes(bounds);ax.axis("off")
    ax.text(0,1.04,title,fontsize=13,color=ACCENT,weight="bold",transform=ax.transAxes)
    chars=[max(18,int(bounds[2]*11.69*72*part/(font*.5))) for part in (.58,.42)]
    data=[[textwrap.fill(str(x),chars[j]) for j,x in enumerate(row)] for row in rows]
    tab=ax.table(cellText=data,colWidths=[.58,.42],cellLoc="left",bbox=(0,0,1,1))
    tab.auto_set_font_size(False);tab.set_fontsize(font)
    for (i,j),cell in tab.get_celld().items():
        cell.set_edgecolor("#d6e1e5");cell.set_linewidth(.4);cell.PAD=.055
        if i%2==0:cell.set_facecolor("#f1f5f7")


def parameters(r):
    c=r.c;m=c["m"];p=m["physicalparams"];a=m["parameter_identity"];s=m["source"];g=m["geometry"]
    circ=p["circuit"];mat=m["material"];cfg=m["observation"].get("config",{})
    fig=r.page("1. Case parameters")
    table(fig,[.055,.525,.425,.34],"Geometry and material",[
        ("Width / resolved length (nm)",fmt(g["width_m"],1e9)+" / "+fmt(g["length_m"],1e9)),
        ("Film thickness / mesh step (nm)",fmt(p["film_thickness_m"],1e9)+" / "+fmt(a["spacing_nm"])),
        ("Tc / bath temperature (K)",fmt(p["Tc_K"])+" / "+fmt(p["bath_temperature_K"])),
        ("Diffusion D (m²/s)",fmt(p["diffusion_m2_s"])),("Sheet resistance (Ω)",fmt(p["sheet_resistance_ohm"])),
        ("Single-spin DOS N0 (1/(J m³))",fmt(p["N0_single_spin_per_J_m3"])),
        ("Sound speed / Debye (m/s; K)",fmt(mat.get("sound_speed_m_per_s"))+" / "+fmt(mat.get("Debye_temperature_K"))),
        ("Derived acoustic γ",fmt(mat.get("gamma_A20"))),("Nodes / edges",f'{g["n_nodes"]} / {g["n_edges"]}'),
        ("Model",model_contract(m)["model_label"])])
    table(fig,[.52,.525,.425,.34],"Photon and calibration parameters",[
        ("Wavelength (nm)",fmt(s.get("wavelength_nm",a["wavelength_nm"]))),
        ("Photon / retained energy (eV)",fmt(s.get("photon_energy_J"),1/EV)+" / "+fmt(s.get("retained_energy_J"),1/EV)),
        ("Retained energy fraction η",fmt(s.get("retained_fraction"))),("Deposition width σ (nm)",fmt(s.get("sigma_nm"))),
        ("Deposition (x,y) in nm",str([round(x*1e9,4) for x in s.get("center_m",[])])),
        ("Free τee(Tc) (ps)",fmt(p["tau_ee_Tc_s"],1e12)),("Free τep(Tc) (ps)",fmt(p["tau_ep_Tc_s"],1e12)),
        ("Derived e-ph τ0 (ns)",fmt(p["tau0_s"],1e9)),("Phonon escape time (ps)",fmt(p["tau_escape_s"],1e12)),
        ("Initial source","DC gap; common FD/BE")])
    table(fig,[.055,.16,.425,.295],"Circuit and ideal discriminator",[
        ("DC current (µA)",fmt(p["current_A"],1e6)),("Rb / Rload (Ω)",fmt(circ["R_bias_ohm"])+" / "+fmt(circ["R_load_ohm"])),
        ("Lb / fixed Lext (nH)",fmt(circ["L_bias_H"],1e9)+" / "+fmt(circ["Lk_ext_H"],1e9)),
        ("Capacitance / Vbias (pF; V)",fmt(circ["C_couple_F"],1e12)+" / "+fmt(circ["V_bias_V"])),
        ("Positive threshold (mV)",fmt(a["threshold_mV"])),("Confirmation / extra time (ps)",fmt(a["confirmation_ps"])+" / "+fmt(a["extra_time_ps"])),
        ("Hysteresis fraction / ceiling (ps)",fmt(cfg.get("hysteresis_fraction",.05))+" / "+fmt(a["ceiling_ps"])),
        ("Clock origin","Postcascade transfer")])
    table(fig,[.52,.16,.425,.295],"Numerical resolution",[
        ("Electronic states / phonon modes",f'{m["kinetics"]["electron_states"]} / {m["kinetics"]["phonon_modes"]}'),
        ("ξ cutoff / GL order / transport",f'{a["xi_cutoff"]:g} / {a["xi_order"]} / {a["transport_order"]}'),
        ("Step cap (ps)",fmt(a["dt_ps"])),("Nonlinear tolerance",fmt(a["step_tolerance"])),
        ("Local energy guard",fmt(a["local_energy_tolerance"])),
        ("Snapshots / checkpoints (ps)",fmt(a["snapshot_every_ps"])+" / "+fmt(a["checkpoint_every_ps"])),
        ("Simulation threads",str(a["threads"])),("Time integrator","Coupled midpoint")])
    r.save(fig,"Parameters","Mesh spacing is resolution, not film thickness. η is retained energy conditional on absorption, not optical efficiency/PCR. "
        "τee and τep are free calibration parameters, not measured constants or prescribed ranges. The postcascade source and acoustic material are exploratory assumptions.")


def observables(r):
    c=r.c;z=c["z"];m=c["m"];v=c["metrics"];t=z["time_s"]*1e12
    baseline=c["obs"].get("baseline_Vout_V",c["obs"].get("config",{}).get("baseline_V",0.))
    voltage=z["Vout_V"]-baseline;threshold=m["parameter_identity"]["threshold_mV"]*1e-3
    scale,unit=(1e6,"µV") if max(abs(voltage))<.02*threshold else (1e3,"mV")
    fig=r.page("2. Observables: readout and detection")
    ax=fig.add_axes([.085,.52,.38,.32]);bx=fig.add_axes([.575,.52,.36,.32])
    ax.plot(t,voltage*scale,color=ACCENT)
    if max(abs(voltage))>=.15*threshold:ax.axhline(threshold*scale,color="black",ls=":",label="Threshold")
    else:ax.text(.03,.91,f"Threshold = {threshold*1e3:g} mV (off scale)",transform=ax.transAxes,fontsize=11)
    ax.set(xlabel="Postcascade time (ps)",ylabel="Vout minus baseline ("+unit+")")
    bx.plot(t,(z["Ib_A"]-z["Is_A"])*1e6,color=ACCENT)
    bx.set(xlabel="Postcascade time (ps)",ylabel="Readout current Ib - Is (µA)")
    for a in (ax,bx):
        a.grid(alpha=.25)
        if v["confirmed"]:
            a.axvline(v["crossing_ps"],color="tab:green",ls=":",label="Crossing")
            a.axvline(v["confirmation_ps"],color="tab:green",ls="--",label="Confirmation")
            a.legend(fontsize=11)
    status=("Provisional; "+c["m"]["status"]+("; confirmed crossing" if v["confirmed"] else "")) if c["provisional"] else (
        "Confirmed detection" if v["confirmed"] else ("Censored without detection" if v["censored"] else "Development window"))
    table(fig,[.075,.16,.86,.245],"Saved events and observables",[
        ("Observation status",status),("Relative latency / crossing (ps)",fmt(v["latency_ps"],digits=9)),
        ("Accepted confirmation (ps)",fmt(v["confirmation_ps"],digits=9)),
        ("Trajectory / snapshot end (ps)",fmt(v["final_ps"],digits=9)+" / "+fmt(v["maps_final_ps"],digits=9)),
        ("Maximum observed Vout (mV)",fmt(v["max_observed_Vout_mV"],digits=8)),
        ("Local crossing slope (mV/ps)",fmt(v["local_slope_mV_ps"],digits=7)),("Minimum Is (µA)",fmt(v["min_Is_uA"],digits=8))],11)
    r.save(fig,"Readout and detection","Vout=Rload(Ib-Is), with DC baseline removed: no filter, amplifier, noise or measurement bandwidth is included. "
        "The slope uses accepted states within crossing ±0.25 ps; it does not predict jitter. The observed maximum is not the complete pulse peak. "
        "Without confirmed crossing, latency is unassigned; a partial checkpoint is not terminal censoring.")


def four_axes(r,title):
    fig=r.page(title)
    return fig,[fig.add_axes(b) for b in ((.085,.56,.37,.28),(.585,.56,.35,.28),(.085,.19,.37,.26),(.585,.19,.35,.26))]


def voltage_units(values):
    peak=float(np.max(abs(values)))
    return (1e3,"mV") if peak>=1e-3 else ((1e6,"µV") if peak>=1e-6 else (1e9,"nV"))


def circuit(r):
    c=r.c;z=c["z"];t=z["time_s"]*1e12
    fig,ax=four_axes(r,"2.1. Electrical circuit observables")
    ax[0].plot(t,z["Ib_A"]*1e6,label="Ib");ax[0].plot(t,z["Is_A"]*1e6,label="Is")
    ax[0].set_ylabel("Current (µA)");ax[0].legend();ax[0].ticklabel_format(axis="y",useOffset=False)
    scale,unit=voltage_units(z["Vdev_V"]);ax[1].plot(t,z["Vdev_V"]*scale);ax[1].set_ylabel("Strip voltage Vdev ("+unit+")")
    scale,unit=voltage_units(z["vc_V"]);ax[2].plot(t,z["vc_V"]*scale);ax[2].set_ylabel("Coupling voltage vc ("+unit+")")
    power=c["m"]["physicalparams"]["circuit"]["R_load_ohm"]*(z["Ib_A"]-z["Is_A"])**2
    scale,unit=(1e9,"nW") if power.max()>=1e-9 else ((1e12,"pW") if power.max()>=1e-12 else (1e15,"fW"))
    ax[3].plot(t,power*scale);ax[3].set_ylabel("Readout power ("+unit+")")
    for a in ax:a.set_xlabel("Postcascade time (ps)");a.grid(alpha=.25)
    r.save(fig,"Circuit","Full three-state circuit. Vdev is the voltage between the resolved strip cuts, distinct from Vout. "
        "Readout power is Rload(Ib-Is)². Only the saved observation window is shown; neither full recovery nor complete pulse energy is inferred.")


def internal(r):
    c=r.c;s=c["s"];t=s["time_s"]*1e12
    dc=abs(c["dc"]["delta_bar"])*c["m"]["physicalparams"]["energy_scale_J"]
    ratio=abs(s["delta_J"])/dc
    center=np.asarray(c["m"]["source"].get("center_m",np.mean(s["coordinates_m"],axis=0)))
    j=int(np.argmin(np.sum((s["coordinates_m"]-center)**2,axis=1)))
    fig,ax=four_axes(r,"3. Internal variables: state evolution")
    ax[0].plot(t,ratio[:,j],label="Node nearest deposition");ax[0].plot(t,ratio.min(axis=1),label="Nodal minimum")
    ax[0].set_ylabel("|Δ| / |ΔDC|");ax[0].legend(fontsize=11)
    ax[1].plot(t,s["Te_K"][:,j],label="Equivalent Te, center");ax[1].plot(t,s["Tph_K"][:,j],label="Equivalent Tph, center")
    ax[1].axhline(c["m"]["physicalparams"]["Tc_K"],color="black",ls=":",label="Tc")
    ax[1].set_ylabel("Equivalent temperature (K)");ax[1].legend(fontsize=11)
    volumes=c["mesh"]["dual_area_m2"]*c["m"]["physicalparams"]["film_thickness_m"]
    for key,label in (("electron_energy_J_m3","Quasiparticles"),("phonon_energy_J_m3","Phonons")):
        energy=s[key]@volumes;ax[2].plot(t,(energy-energy[0])/EV,label=label)
    ax[2].set_ylabel("Energy change from t=0 (eV)");ax[2].legend(fontsize=11)
    ax[3].step(t,c["hotbelt"].astype(int),where="post")
    ax[3].set(ylabel="Saved geometric hotbelt",yticks=[0,1],yticklabels=["No","Yes"],ylim=(-.1,1.1))
    for a in ax:a.set_xlabel("Postcascade time (ps)");a.grid(alpha=.25)
    xy=s["coordinates_m"][j]*1e9
    r.save(fig,"Internal variables",f"Nearest deposition node: ({xy[0]:.5g}, {xy[1]:.5g}) nm. Te/Tph are energy-equivalent observers, not imposed thermal states. "
        "Saved electronic energy is quasiparticle energy, excluding condensation/gradient. Hotbelt: nodes with |Δ|≤0.2|ΔDC| connect transverse edges; "
        "this is distinct from the electrical trigger and does not establish continuity between snapshots.")


def snapshot_axes(fig):
    return [fig.add_axes([x,y,.25,.21]) for y in (.64,.30) for x in (.075,.385,.695)]


def snapshot_title(c,slot,index):
    actual=c["s"]["time_s"][index]*1e12;target=c["snapshot_targets_ps"][slot]
    return f"t = {actual:.6g} ps\nrequested: {target:.6g} ps"


PALETTES=dict(gap="viridis",phase="twilight",electron_temperature="inferno",
    phonon_temperature="YlOrBr",energy="PuOr_r",signed_current="coolwarm",total_current="cividis",potential="GnBu")


def map_pages(r,indices,key,title,label,caption,kind="scalar",roi=False,flow=None,cmap=None,signed=False):
    if len(indices)!=6:raise ValueError("Every snapshot mosaic requires exactly six slots")
    c=r.c;s=c["s"];mesh=c["mesh"]
    center=np.asarray(c["m"]["source"].get("center_m",np.mean(s["coordinates_m"],axis=0)))
    xy=s["coordinates_m"].copy()*1e9;xy[:,0]-=center[0]*1e9
    tri=Triangulation(xy[:,0],xy[:,1],mesh["triangles"])
    bounds=(max(float(xy[:,0].min()),-60),min(float(xy[:,0].max()),60),float(xy[:,1].min()),float(xy[:,1].max())) if roi else (
        float(xy[:,0].min()),float(xy[:,0].max()),float(xy[:,1].min()),float(xy[:,1].max()))
    cmap=cmap or PALETTES["gap"]
    if kind=="gap":
        values={int(i):s["delta_J"][i] for i in indices};lo=0.;hi=max(float(abs(v).max()/MEV) for v in values.values())
    elif kind=="phase":
        dc=c["dc"]["delta_bar"];Eu=c["m"]["physicalparams"]["energy_scale_J"]
        values={int(i):np.ma.array(np.angle(s["delta_J"][i]*np.conj(dc)),mask=abs(s["delta_J"][i])<=.2*abs(dc)*Eu) for i in indices}
        lo,hi=-np.pi,np.pi;cmap=PALETTES["phase"]
    elif kind=="current":
        values={int(i):(np.linalg.norm(flow[int(i)][key],axis=1) if key=="total" else flow[int(i)][key][:,0])/1e10 for i in indices}
        hi=max(float(abs(v).max()) for v in values.values());lo=0. if key=="total" else -hi
        cmap=PALETTES["total_current" if key=="total" else "signed_current"]
        arrow_reference=max(float(np.linalg.norm(flow[int(i)][key],axis=1).max()) for i in indices)
    elif kind=="phi":
        values={int(i):flow[int(i)]["phi"]*1e3 for i in indices}
        lo=min(float(v.min()) for v in values.values());hi=max(float(v.max()) for v in values.values());cmap=PALETTES["potential"]
    else:
        factor=1e-3 if key.endswith("J_m3") else 1.
        values={int(i):s[key][i]*factor for i in indices}
        lo=min(float(v.min()) for v in values.values());hi=max(float(v.max()) for v in values.values())
        if signed:hi=max(abs(lo),abs(hi));lo=-hi;cmap=PALETTES["energy"]
    if hi<=lo:hi=lo+max(abs(lo)*1e-6,1e-12)
    fig=r.page(title);axes=snapshot_axes(fig);image=None
    for pos,(ax,i) in enumerate(zip(axes,indices)):
        v=values[int(i)]
        if kind=="gap":
            image=draw_complex_gap(ax,tri,v,unit_J=MEV,vmin=lo,vmax=hi,bounds=bounds,samples=401,cmap=cmap)
        elif kind=="phase":
            keep=~np.ma.getmaskarray(v)
            image=ax.scatter(xy[keep,0],xy[keep,1],c=np.asarray(v)[keep],s=6,edgecolors="none",cmap=cmap,vmin=lo,vmax=hi,rasterized=True)
        elif kind=="current":
            image=ax.tripcolor(tri,facecolors=v,shading="flat",cmap=cmap,vmin=lo,vmax=hi,rasterized=True)
            centers=xy[mesh["triangles"]].mean(axis=1);good=(centers[:,0]>=bounds[0])&(centers[:,0]<=bounds[1])
            ids=np.flatnonzero(good)[::max(1,int(good.sum()/85))];vectors=flow[int(i)][key][ids]
            norm=np.linalg.norm(vectors,axis=1);keep=norm>max(arrow_reference*1e-4,1e-30)
            direction=vectors[keep]/norm[keep,None]
            ax.quiver(centers[ids[keep],0],centers[ids[keep],1],direction[:,0]*6,direction[:,1]*6,angles="xy",scale_units="xy",scale=1,width=.003,color="black")
        else:image=ax.tripcolor(tri,v,shading="gouraud",cmap=cmap,vmin=lo,vmax=hi,rasterized=True)
        ax.set(xlim=bounds[:2],ylim=bounds[2:],xlabel="x - xdep (nm)",ylabel="y (nm)")
        ax.set_title(snapshot_title(c,pos,i),fontsize=12,pad=8);ax.set_aspect("equal")
        ax.tick_params(labelsize=12);ax.xaxis.label.set_size(13);ax.yaxis.label.set_size(13)
        ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5));ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3))
    cb=fig.colorbar(image,cax=fig.add_axes([.245,.215,.52,.02]),orientation="horizontal")
    cb.set_label(label,fontsize=14);cb.ax.tick_params(labelsize=12)
    r.save(fig,title,caption+" Physical coordinates; x centered on deposition. Exactly six nearest saved states, with repetitions retained; "
        "no temporal interpolation. One color range across the six slots for this quantity.")


def vorticity_page(r,indices):
    c=r.c;s=c["s"];mesh=c["mesh"];xy=s["coordinates_m"].copy()*1e9
    center=np.asarray(c["m"]["source"]["center_m"]);xy[:,0]-=center[0]*1e9
    tri=Triangulation(xy[:,0],xy[:,1],mesh["triangles"])
    dc=c["dc"]["delta_bar"]*c["m"]["physicalparams"]["energy_scale_J"]
    diagnostics={int(i):phase_vorticity(s["delta_J"][i],dc,s["coordinates_m"],mesh["triangles"]) for i in np.unique(indices)}
    fig=r.page("3.11. Vorticity: signed condensate phase winding")
    cmap=matplotlib.colors.ListedColormap(["#008176","#ffffff","#ba287f"]);cmap.set_bad("#b9b9b9")
    norm=matplotlib.colors.BoundaryNorm([-1.5,-.5,.5,1.5],3)
    for pos,(ax,i) in enumerate(zip(snapshot_axes(fig),indices)):
        d=diagnostics[int(i)];values=np.ma.array(d["winding"],mask=~d["valid"])
        # Explicit gray geometry remains visible where masked phase is undefined.
        ax.tripcolor(tri,facecolors=np.zeros(len(mesh["triangles"])),cmap=matplotlib.colors.ListedColormap(["#b9b9b9"]),vmin=0,vmax=1,rasterized=True)
        image=ax.tripcolor(tri,facecolors=values,cmap=cmap,norm=norm,shading="flat",rasterized=True)
        for charge,marker,color in [(1,"+","#ba287f"),(-1,"x","#008176")]:
            pts=d["core_points_m"][d["core_charge"]==charge]*1e9
            if len(pts):ax.scatter(pts[:,0]-center[0]*1e9,pts[:,1],marker=marker,s=65,c=color,linewidths=1.5)
        ax.set(xlabel="x - xdep (nm)",ylabel="y (nm)",xlim=(xy[:,0].min(),xy[:,0].max()),ylim=(xy[:,1].min(),xy[:,1].max()))
        ax.set_title(snapshot_title(c,pos,i),fontsize=12,pad=8);ax.set_aspect("equal")
        ax.tick_params(labelsize=12);ax.xaxis.label.set_size(13);ax.yaxis.label.set_size(13)
        ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5));ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3))
        ax.text(.03,.03,f'+: {d["positive"]}  -: {d["negative"]}  gray: {100*(~d["valid"]).mean():.1f}%',transform=ax.transAxes,
            fontsize=10,bbox=dict(facecolor="white",alpha=.85,edgecolor="none"))
    cb=fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm,cmap=cmap),
        cax=fig.add_axes([.245,.215,.52,.02]),orientation="horizontal",ticks=[-1,0,1])
    cb.set_label("Counterclockwise phase winding n (dimensionless)",fontsize=14);cb.ax.tick_params(labelsize=12)
    r.save(fig,"Phase-winding vorticity","n = sum of wrapped phase increments / (2π), around each counterclockwise triangle; +/x mark P1 vortex/antivortex cores. "
        "Gray: a contour gap ≤10^-6 of its local DC reference, or unresolved zero. This is phase topology, not curl(J) or magnetic flux. "
        "Saved A=0 states only; no trajectories or missed events are inferred. Method: Guo et al. (2015), Eq. (1), linked reference.",source_url="https://www.mcs.anl.gov/papers/P5324-0315.pdf")
    return [dict(snapshot_index=int(i),saved_time_ps=float(s["time_s"][i]*1e12),
        **{k:d[k] for k in ("positive","negative","indeterminate_cells","resolved_cells","contour_floor_relative_DC","convention","reference")},
        core_points_nm=(d["core_points_m"]*1e9).tolist(),core_charge=d["core_charge"].tolist())
        for i,d in diagnostics.items()]


def balance(r):
    c=r.c;z=c["z"];t=z["time_s"]*1e12;fig,ax=four_axes(r,"4. Independent energy balance")
    ax[0].plot(t,(z["energy_total_J"]-z["energy_total_J"][0])/EV,label="Total ΔE")
    ax[0].plot(t,(z["external_energy_integral_J"]-z["external_energy_integral_J"][0])/EV,"--",label="Net external work")
    ax[0].set_ylabel("Energy (eV)");ax[0].legend(fontsize=11)
    ref=c["m"]["source"].get("retained_energy_J",0.)
    if ref<=0:ref=c["m"]["numerical"]["kinetic_reference_energy_J"]
    work=np.r_[abs(z["external_energy_integral_J"][0]),abs(np.diff(z["external_energy_integral_J"]))].cumsum()
    defect=abs(z["energy_defect_J"])/(ref+work)*100
    ax[1].plot(t,defect,label="|Defect| / (Eref + Σ|ΔW|)");ax[1].axhline(1.,color="black",ls=":",label="Global criterion: 1 %")
    ax[1].set_ylabel("Global energy defect (%)");ax[1].legend(fontsize=11)
    ax[2].plot(t,z["escape_power_W"]*1e9,label="Escape")
    ax[2].plot(t,z["reservoir_power_out_W"]*1e9,label="Cuts, signed outward")
    ax[2].set_ylabel("Thermal power (nW)");ax[2].legend(fontsize=11)
    ax[3].plot(t,z["instant_power_residual_W"]);ax[3].set_ylabel("Power residual (W)")
    ax[3].ticklabel_format(axis="y",style="sci",scilimits=(0,0))
    for a in ax:a.set_xlabel("Postcascade time (ps)");a.grid(alpha=.25)
    c["balance"]=dict(maximum_global_energy_defect_percent=float(defect.max()),
        maximum_absolute_instant_power_residual_W=float(abs(z["instant_power_residual_W"]).max()))
    r.save(fig,"Energy balance",("Stored energy is U_BCS+G-sum(Te*dG/dTe), phonons and circuit. Force/current use free G at fixed Te. " if thermal_candidate(c["m"]) else "Total energy includes the BCS film, fixed-bath gradient, phonons and circuit. ")+
        "Net external work subtracts resistive loss, escape and cuts. "
        "Normalization: retained energy plus accumulated absolute external work. This report does not replace numerical admission or experimental calibration.")


def provenance(r,checks):
    c=r.c;z=c["z"];fig=r.page("5. Precision, provenance and scope")
    ax=fig.add_axes([.085,.63,.37,.20]);keep=z["accepted_step_ps"]>0
    ax.plot(z["time_s"][keep]*1e12,z["accepted_step_ps"][keep])
    ax.set(xlabel="Time (ps)",ylabel="Accepted step (ps)",yscale="log");ax.grid(alpha=.25)
    table(fig,[.55,.56,.39,.28],"Report checks",[
        ("Associated terminal audit",c["audit"]["status"]),("Accepted / rejected steps",str(len(z["time_s"])-1)+" / "+str(c["execution"].get("total_rejections","not saved"))),
        ("Session wall time (h)",fmt(c["execution"].get("elapsed_session_wall_s"),1/3600)),("Positive global currents",str(c["metrics"]["positive_global_currents"])),
        ("Readout identity error (V)",fmt(c["metrics"]["readout_identity_error_V"])),
        ("Maximum Vdev recovery error (V)",fmt(max(x["Vdev_error_V"] for x in checks))),
        ("Maximum global energy defect (%)",fmt(c["balance"]["maximum_global_energy_defect_percent"])),
        ("Trajectory / snapshot end (ps)",fmt(c["metrics"]["final_ps"])+" / "+fmt(c["metrics"]["maps_final_ps"]))],10)
    lines=["Data: "+str(c["run"]),"Manifest SHA256: "+c["inputs"]["manifest.json"],
        ("Captured checkpoint SHA256: "+c["inputs"]["checkpoint_capture.npz"]) if c["provisional"] else ("Trajectory SHA256: "+c["inputs"]["trajectory.npz"]),
        "Mesh SHA256: "+c["inputs"]["geometry_mesh"],"Pipeline SHA256: "+sha256(__file__),
        ("Currents: fixed-Te free-gradient gauge derivative using saved Te_K and algebraic neutrality from pinned sources. " if thermal_candidate(c["m"]) else "Currents: exact fixed-bath BCS gradient gauge derivative and algebraic neutrality from pinned sources. ") +
        "Display vectors use weighted dual-flux least squares; they are not additional model degrees of freedom.",
        "Complete p/n were not saved at every snapshot: collision-rate, nonclosure DOS and KWT-heating maps are unavailable.",
        ("Thermal depairing modifies G; native collision/transport spectra retain q=0. Acoustic material, KWT mobility and the postcascade source remain effective assumptions. "
         if thermal_depairing_candidate(c["m"]) else
         "Nonequilibrium radial rigidity/KWT, acoustic material and the postcascade source are effective assumptions. ")+
        "This window does not establish jitter, statistical efficiency, absolute optical latency or full reset."]
    y=.50
    for line in lines:
        text=textwrap.fill(line,128);fig.text(.075,y,text,fontsize=10,color=INK,va="top",linespacing=1.25)
        y-=.027*len(text.splitlines())+.012
    r.save(fig,"Provenance","Individual case; no inter-run comparison or new model admission. A provisional checkpoint is not a terminal result. "
        "Input hashes, requested/saved snapshot times, repetitions, winding diagnostics and reconstruction errors are retained in the JSON receipt.")


def equation_box(fig,x,y,equation,size=14):
    fig.text(x,y,equation,fontsize=size,ha="center",va="center",color=INK,
        bbox=dict(boxstyle="round,pad=.6",facecolor="#f1f5f7",edgecolor=ACCENT,linewidth=.8))


def model_summary(r):
    contract=model_contract(r.c["m"]);thermal=contract["candidate"];depairing=contract["thermal_depairing"]
    fig=r.page("6. Model summary: thermal-depairing candidate" if depairing else "6. Model summary: "+contract["model_label"])
    if depairing:fig.text(.075,.883,contract["model_label"],fontsize=12.5,color=ACCENT)
    fig.text(.075,.85 if depairing else .863,r"$X=\{d_i,p_{ik},n_{iq},I_b,I_s,v_c\},\quad d_i=\Delta_i/E_u,\quad E_u=k_BT_c,\quad E_{ik}=E_u\sqrt{\bar\xi_k^2+|d_i|^2}$",fontsize=14,color=INK)
    fig.text(.075,.807,"Coupled evolution",fontsize=14,weight="bold",color=ACCENT)
    equation_box(fig,.335,.744,r"$D_td_i=-\frac{z_i}{A_i\tau_K}\left[\frac{\mathrm{Re}\,f_i}{R_i}+iR_i\,\mathrm{Im}\,f_i\right]$",15)
    equation_box(fig,.335,.617,r"$\dot p=\mathcal{C}_{e-ph}+\mathcal{T}_E+\frac{f_{\rm FD}(E,T_E)-p}{\tau_{\rm BGK}}+\mathcal{H}$",14)
    equation_box(fig,.335,.497,r"$\dot n=\mathcal{C}_{ph-e}-\frac{n-n_{\rm BE}(\Omega,T_b)}{\tau_{\rm esc}}$",15)
    fig.text(.65,.807,"Algebraic neutrality",fontsize=14,weight="bold",color=ACCENT)
    equation_box(fig,.79,.746,r"$B_g[I^s+GB_g^{\mathsf{T}}\Phi]=s_g I_s$",13)
    fig.text(.65,.68,r"$D_t=\partial_t+2ie\phi/\hbar,\quad z=d/|d|$",fontsize=12.5)
    fig.text(.65,.625,r"$A_i=\sqrt{(1+\vartheta_i)/2},\quad \vartheta_i=\max(T_{E,i},T_b)/T_c$",fontsize=11.5)
    fig.text(.65,.57,r"$\tau_{{\rm sc},i}^{-1}=\vartheta_i/\tau_{ee}(T_c)+\vartheta_i^3/\tau_{ep}(T_c)$",fontsize=11.5)
    fig.text(.65,.515,r"$R_i=\sqrt{1+4E_u^2\tau_{{\rm sc},i}^2|d_i|^2/\hbar^2}$",fontsize=12.5)
    fig.text(.65,.46,r"$\tau_K=\pi\hbar/(8k_BT_c),\quad \tau_{\rm BGK}=\tau_{ee}(T_c)/\vartheta$",fontsize=11.5)
    if thermal:
        fig.text(.075,.405,"Free gradient, stored energy and thermodynamic force",fontsize=14,weight="bold",color=ACCENT)
        if depairing:
            fig.text(.075,.355,r"$G=\sum_i\{\kappa(\rho_i,T_b)L_i+u_i[\varphi(\rho_i,\theta_i,Q_i)-\varphi(\rho_i,\theta_b,Q_i)]\}$",fontsize=13)
            fig.text(.075,.315,r"$\rho_i=|d_i|^2,\quad \theta_i=T_i/T_c,\quad \theta_b=T_b/T_c,\quad Q_i=L_i/(u_i\rho_i),\quad U=U_{\rm BCS}+G-\sum_iT_iG_{T_i}$",fontsize=12)
        else:
            fig.text(.075,.35,r"$G=\sum_{\langle ij\rangle}C_{ij}\kappa_{ij}(T_E)|d_i-e^{-i\alpha_{ij}}d_j|^2,\quad U=U_{\rm BCS}+G-\sum_iT_iG_{T_i}$",fontsize=13)
        fig.text(.075,.267 if depairing else .29,r"$S=S_{\rm qp}-\sum_iG_{T_i},\quad {\cal F}_i=\partial_{d_i^*}U-\sum_jT_j\partial_{d_i^*}S_j=\left.\partial_{d_i^*}(U_{\rm BCS}+G)\right|_T$",fontsize=12.5)
        fig.text(.075,.225,r"$I^s_{ij}=-\frac{2e}{\hbar}\left.\partial_{\alpha_{ij}}G\right|_T,\quad C_{{\rm eff},i}=C_{{\rm FD},i}-T_iG_{T_iT_i}>0$",fontsize=13.5)
        if contract["pauli_signed_heat"]:
            fig.text(.075,.18,r"$\mathcal{H}_i=h_iP_i,\quad P_i=Q_{\Delta,i}+P_{J,i}+T_i\dot G_{T_i},\quad \epsilon_i\cdot h_i=1,\quad P_i(\sigma_i\cdot h_i-1/T_i)\geq0$",fontsize=11.5,color=INK)
            fig.text(.075,.145,"Pauli heat: f(1-p) for heating, p(1-f) for cooling; exact FD limit and nonnegative excess entropy.",fontsize=10.5,color=INK)
        else:
            fig.text(.075,.18,r"$\mathcal{H}_i=h_i[Q_{\Delta,i}+P_{J,i}+T_i\dot G_{T_i}],\quad \epsilon_i\cdot h_i=1,\quad \partial_{p_i}S_{\rm qp}\cdot h_i=1/T_i$",fontsize=11.5,color=INK)
            fig.text(.075,.145,"Historical 1.1.0 dual-moment heat carrier; superseded after a Pauli-viability failure.",fontsize=10.5,color=INK)
        caption=("The stored internal energy differs from free G; force/current use fixed-Te derivatives; native collision/transport spectra retain q=0. "
            "L is the nodal gradient measure; varphi is an equilibrium finite-depairing correction. "
            "Matching adds T-independent H=kappa(Tb)L-u*varphi(Tb,Q): it may dominate strong gradients and persist at fixed gap as Te increases. "
            "Ceff>0 is required without a floor; DC is matched at Tb. This does not establish microscopic accuracy or quantitative latency."
            if depairing else
            "Energy-equivalent Te is read directly from each saved state for current reconstruction. "
            "The stored internal energy differs from free G; force/current follow fixed-Te derivatives. "
            "The signed source closes energy; 1.1.1 adds Pauli-compatible excess entropy outside FD. Ceff>0 is required without a floor. "
            "This candidate has no global trajectory admission and is not exact microscopic Usadel.")
        r.save(fig,"Thermal-depairing candidate" if depairing else "Thermal-stiffness candidate",caption)
    else:
        fig.text(.075,.405,"Shared BCS energy and variational force",fontsize=14,weight="bold",color=ACCENT)
        fig.text(.075,.35,r"$U_{e+\nabla}=\sum_i u_i^*\left[\frac{|d_i|^2}{g}-2\sum_k w_k(e_{ik}-\bar\xi_k)+4\sum_k w_ke_{ik}p_{ik}\right]+U_\nabla$",fontsize=14)
        fig.text(.075,.29,r"$U_\nabla=\sum_{\langle ij\rangle}C_{ij}\kappa_{ij}|d_i-e^{-i\alpha_{ij}}d_j|^2,\quad \mathcal{F}_i=\partial_{d_i^*}U_{e+\nabla},\quad f_i=z_i^*\mathcal{F}_i/u_i^*$",fontsize=13.5)
        fig.text(.075,.225,r"$I^s_{ij}=-\frac{2e}{\hbar}\partial_{\alpha_{ij}}U_{e+\nabla},\quad u_i^*=N_0E_u^2V_i,\quad V_i=d_f A_i^{\rm dual},\quad \alpha_{ij}=0$",fontsize=13.5)
        fig.text(.075,.18,"p/n: quasiparticle/phonon occupations; C: shared conserving collisions; T: energy transport; H: condensate + Joule heat.",fontsize=10.5,color=INK)
        fig.text(.075,.145,"Bg: reduced incidence; G: normal conductance; Φ: reduced potential; sg: cut injection; κ: radial stiffness; g: BCS coupling.",fontsize=10.5,color=INK)
        r.save(fig,"Effective neutral BCS/KWT",
            "Neutral: quasistatic current balance determines φ, without a dynamical charge mode. Effective: radial stiffness, acoustic material and hot-state mobility are calibratable closures. "
            "BCS/KWT: a shared BCS energy drives Kramer-Watts-Tobin condensate relaxation and conservative kinetics. Te/Tph are energy-equivalent observers, not imposed thermal states. "
            "See Principal.pdf and MODELO_VIGENTE.md; this is not an exact microscopic Usadel evolution.")

def circuit_diagram(ax):
    from matplotlib.patches import Ellipse,Rectangle
    ax.set(xlim=(0,1),ylim=(0,1));ax.axis("off")
    def wire(x,y):ax.plot(x,y,color=INK,lw=1.7)
    def resistor(a,b,vertical=False):
        q=np.linspace(0,1,11);amp=np.array([0,0,1,-1,1,-1,1,-1,1,0,0])*.013
        if vertical:wire(np.full(11,a[0])+amp,a[1]+q*(b[1]-a[1]))
        else:wire(a[0]+q*(b[0]-a[0]),np.full(11,a[1])+amp*2.3)
    def inductor(a,b,vertical=False):
        q=np.linspace(0,1,160);coil=.025*np.sin(np.pi*(q*4%1))
        if vertical:wire(np.full_like(q,a[0])+coil*.5,a[1]+q*(b[1]-a[1]))
        else:wire(a[0]+q*(b[0]-a[0]),np.full_like(q,a[1])+coil)
    def arrow(a,b,label,offset=(0,0)):
        ax.annotate("",xy=b,xytext=a,arrowprops=dict(arrowstyle="->",color=ACCENT,lw=1.7))
        ax.text((a[0]+b[0])/2+offset[0],(a[1]+b[1])/2+offset[1],label,fontsize=13,ha="center",va="center",color=ACCENT)
    wire([.10,.86],[.14,.14]);wire([.10,.10],[.14,.40]);wire([.10,.10,.22],[.58,.83,.83])
    ax.add_patch(Ellipse((.10,.49),.064,.18,facecolor="white",edgecolor=INK,lw=1.7))
    ax.text(.10,.53,"+",ha="center",va="center",fontsize=14);ax.text(.10,.45,"-",ha="center",va="center",fontsize=14)
    ax.text(.045,.49,r"$V_b$",ha="right",va="center",fontsize=14)
    resistor((.22,.83),(.32,.83));wire([.32,.38],[.83,.83]);inductor((.38,.83),(.48,.83));wire([.48,.57],[.83,.83])
    ax.text(.27,.94,r"$R_b$",ha="center",fontsize=14);ax.text(.43,.94,r"$L_b$",ha="center",fontsize=14)
    arrow((.13,.72),(.21,.72),r"$I_b$",(0,-.07))
    wire([.57,.57],[.83,.70]);inductor((.57,.70),(.57,.54),True);wire([.57,.57],[.54,.43])
    ax.text(.49,.60,r"$L_{k,\rm ext}$",ha="right",va="center",fontsize=14)
    ax.add_patch(Rectangle((.51,.25),.12,.18,facecolor="#edf4f5",edgecolor=ACCENT,lw=1.7))
    ax.text(.57,.34,"2D strip",ha="center",va="center",fontsize=12)
    ax.text(.65,.34,r"$V_{\rm dev}$",ha="left",va="center",fontsize=14)
    wire([.57,.57],[.25,.14]);arrow((.64,.72),(.64,.53),r"$I_s$",(.035,0))
    wire([.57,.70],[.83,.83]);wire([.70,.70],[.75,.91]);wire([.74,.74],[.75,.91]);wire([.74,.86],[.83,.83])
    ax.text(.72,.98,r"$+\,v_c\,-$",ha="center",fontsize=14)
    ax.text(.72,.67,r"$C_c$",ha="center",fontsize=14)
    arrow((.77,.70),(.84,.70),r"$I_{\rm read}=I_b-I_s$",(0,-.12))
    wire([.86,.86],[.83,.62]);resistor((.86,.62),(.86,.35),True);wire([.86,.86],[.35,.14])
    ax.text(.90,.48,r"$R_l$",ha="left",va="center",fontsize=14)
    ax.text(.90,.80,r"$+\,V_{\rm out}$",ha="left",va="center",fontsize=14)
    ax.text(.90,.16,"-",ha="left",va="center",fontsize=14)
    ax.scatter([.57,.86],[.83,.83],s=20,c=INK)
    wire([.43,.43],[.14,.08])
    for yy,w in [(.08,.035),(.05,.023),(.02,.012)]:wire([.43-w,.43+w],[yy,yy])


def model_circuit_summary(r):
    fig=r.page("6.1. Circuit coupling and conservation")
    circuit_diagram(fig.add_axes([.055,.46,.88,.40]))
    fig.text(.075,.405,r"$R_l=R_{\rm load},\quad V_b=V_{\rm bias},\quad I_{\rm read}=I_b-I_s,\quad V_{\rm out}=R_l I_{\rm read}$",fontsize=14)
    equation_box(fig,.50,.325,r"$L_b\dot I_b=V_b-R_bI_b-v_c-R_l(I_b-I_s)$",15)
    equation_box(fig,.50,.245,r"$L_{k,\rm ext}\dot I_s=v_c+R_l(I_b-I_s)-V_{\rm dev},\qquad C_c\dot v_c=I_b-I_s$",15)
    fig.text(.075,.155,r"$\dot E_{\rm total}=V_b I_b-R_bI_b^2-R_l(I_b-I_s)^2-P_{\rm esc}-P_{\rm cuts}$",fontsize=15,color=INK)
    r.save(fig,"Circuit and conservation","The three circuit states evolve jointly with d, p and n; Vdev is obtained from the same algebraic neutrality solve. "
        "Exterior inductance is fixed from the DC reference. The 2D strip is resolved, not a time-varying lumped resistor/inductor. "+
        ("Condensate/Joule and signed reversible Te*dG_T/dt enter quasiparticles; stored U includes -sum(Te*G_T). " if thermal_candidate(r.c["m"]) else "Shared condensate/Joule heating enters quasiparticles; ") +
        "Phonon escape and superconducting cuts exchange energy explicitly. No filter/amplifier is included.")

def build_report(run,output,allow_incomplete=False,requested=None,dpi=160,memory_limit=2**33):
    start=time.monotonic();run=Path(run).resolve();output=scratch_path(output)
    rid=case_id(run);target=output/("case_report_"+rid+".pdf");receipt_path=target.with_suffix(".json")
    if target.exists() or receipt_path.exists():raise FileExistsError("Existing report preserved; choose a fresh destination")
    build=output/".build"/rid;build.mkdir(parents=True,exist_ok=False)
    c=load_case(run,build,allow_incomplete,memory_limit);c["metrics"]=metrics(c)
    dc=c["dc"]["delta_bar"]*c["m"]["physicalparams"]["energy_scale_J"]
    c["hotbelt"]=hotbelts(c["s"]["delta_J"],dc,c["mesh"]["edges"],c["mesh"]["coordinates_m"])
    indices=select_snapshots(c["s"]["time_s"],c["obs"],c["hotbelt"],requested)
    c["snapshot_targets_ps"]=snapshot_targets(c["s"]["time_s"],c["obs"],c["hotbelt"],requested)
    selection=snapshot_selection(c["s"]["time_s"],c["snapshot_targets_ps"],indices)
    print(json.dumps(dict(status="CASE_REPORT_STARTED",case=run.name,provisional=c["provisional"],snapshots=len(indices))),flush=True)
    flow,checks=currents(c,np.unique(indices))
    apply_thesis_style();plt.rcParams.update({"font.size":13,"axes.labelsize":14,"axes.titlesize":13,"legend.fontsize":11,"xtick.labelsize":12,"ytick.labelsize":12})
    temporary=build/"case_report.pdf";r=Report(temporary,c,dpi)
    try:
        parameters(r);observables(r);circuit(r);internal(r)
        map_pages(r,indices,"delta_J","3.1. Condensate: full domain","|Δ| (meV)","Magnitude of the spatial complex P1 interpolant; the nodal magnitude is not interpolated.",kind="gap")
        map_pages(r,indices,"delta_J","3.2. Condensate: deposition region","|Δ| (meV)","Central 120 nm view, full width; not a thermal boundary.",kind="gap",roi=True)
        map_pages(r,indices,"delta_J","3.3. Condensate: phase relative to DC","arg[Δ(t) ΔDC*] (rad)","Nodal phase without angle interpolation. Nodes with |Δ|≤0.2|ΔDC| are hidden. Wrapped ±π jumps alone do not establish a vortex.",kind="phase",roi=True)
        map_pages(r,indices,"Te_K","3.4. Equivalent electronic temperature","Equivalent Te (K)","Energy-equivalent observer of nonthermal electronic occupations; not an imposed temperature.",cmap=PALETTES["electron_temperature"])
        map_pages(r,indices,"Tph_K","3.5. Equivalent phonon temperature","Equivalent Tph (K)","Phonon energy-equivalent observer; local escape, without lateral phonon diffusion.",cmap=PALETTES["phonon_temperature"])
        map_pages(r,indices,"u_excess_J_m3","3.6. Energy density relative to DC","Excess energy density (kJ/m³)","Local BCS energy, apportioned gradient and phonons, minus the saved DC reference, divided by dual volume.",signed=True)
        for key,title,label in (("super","3.7. Reconstructed supercurrent density",r"$J_{s,x}$ ($10^{10}$ A/m$^2$)"),
            ("normal","3.8. Reconstructed normal current density",r"$J_{n,x}$ ($10^{10}$ A/m$^2$)"),("total","3.9. Reconstructed total current density",r"$|J_{\rm total}|$ ($10^{10}$ A/m$^2$)")):
            map_pages(r,indices,key,title,label,"Dual-flux least-squares display vectors. Arrows: direction, fixed 6 nm length; flows below 0.01 % of the selected maximum omitted. "
                "Js/Jn color: signed x; total: magnitude. Local sign does not imply global reversal.",kind="current",flow=flow)
        map_pages(r,indices,"phi","3.10. Algebraic electric potential","φ (mV)","Neutrality at the same saved state. Gauge: right cut φ=0; no new evolution.",kind="phi",flow=flow)
        vortex_checks=vorticity_page(r,indices)
        balance(r);provenance(r,checks);model_summary(r);model_circuit_summary(r)
    finally:r.pdf.close()
    receipt=dict(schema=SCHEMA,run_dir=str(run),case_name=run.name,report_path=str(target),report_sha256=sha256(temporary),
        generated_utc=datetime.now(timezone.utc).isoformat(),pages=r.pages,pipeline_sha256=sha256(__file__),
        inputs_sha256=c["inputs"],complete_input=not c["provisional"],provisional=c["provisional"],
        model_simulations_launched=0,kinetic_RHS_evaluations=0,new_DC_solves=0,cross_run_comparison=False,
        electrical_metrics=c["metrics"],energy_balance_metrics=c["balance"],existing_audit=c["audit"],
        selected_saved_times_ps=(c["s"]["time_s"][indices]*1e12).tolist(),snapshot_selection=selection,
        requested_times_ignored_after_six=0 if requested is None else max(0,len(requested)-6),
        palettes=PALETTES,vorticity_checks=vortex_checks,reconstruction_checks=checks,
        reconstruction_source_sha256={p:sha256(ROOT/p) for p in reconstruction_sources(c["m"])},
        model_contract=model_contract(c["m"]),
        hotbelt=dict(threshold_DC=.2,definition="saved-node graph connectivity",
            first_saved_time_ps=None if not np.any(c["hotbelt"]) else float(c["s"]["time_s"][np.flatnonzero(c["hotbelt"])[0]]*1e12)),
        memory=c["memory"],actual_cpu_affinity=sorted(os.sched_getaffinity(0)),dpi=dpi,
        absent_instrument_response="No filter, amplifier, measurement bandwidth or noise on Vout",
        unavailable_snapshot_fields=["complete distributions","collision-rate maps","KWT-heating maps"],
        scope="Single saved finite-window effective case; no calibration/convergence admission inferred",elapsed_s=time.monotonic()-start)
    temp=build/"receipt.json";temp.write_text(json.dumps(receipt,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    os.replace(temporary,target);os.replace(temp,receipt_path)
    print(json.dumps(dict(status="CASE_REPORT_COMPLETE",case=run.name,pages=len(r.pages),pdf=str(target),elapsed_s=receipt["elapsed_s"])),flush=True)
    return receipt


def initialize_worker(cpus,counter):
    with counter.get_lock():
        index=counter.value;counter.value+=1
    os.sched_setaffinity(0,{cpus[index]})


def worker(job):
    return build_report(**job)


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run-dir",type=Path,nargs="+",required=True,help="One independent PDF per directory")
    p.add_argument("--output-dir",type=Path,required=True,help="NEW directory under scratch")
    p.add_argument("--allow-incomplete",action="store_true",help="Capture one atomic checkpoint, explicitly provisional")
    p.add_argument("--snapshots-ps",type=float,nargs="+",default=None,help="Exactly six slots: first six requests, repeat last if fewer; nearest saved states")
    p.add_argument("--workers",type=int,default=1,help="Independent foreground reports; one physical core per worker")
    p.add_argument("--dpi",type=int,default=THESIS_DPI);a=p.parse_args(argv)
    if a.workers<1 or not 80<=a.dpi<=300:raise ValueError("Positive workers and80<=dpi<=300 required")
    runs=[x.resolve() for x in a.run_dir]
    if len(set(runs))!=len(runs):raise ValueError("Duplicate case input")
    output=scratch_path(a.output_dir)
    if output.exists():raise FileExistsError("Existing reports preserved; choose a new output directory")
    from sandbox.stage4_core.parallel_runtime import linux_resources,resource_budget
    resources=linux_resources();wanted=min(a.workers,len(runs));budget=resource_budget(resources,.9,.9,max_workers=wanted)
    allowed=set(budget["selected_affinity_cpus"]);coordinator=budget["coordinator_cpu"]
    cores=[g for g in resources["physical_core_groups"] if set(g)<=allowed and coordinator not in g]
    if wanted>budget["workers"] or wanted>len(cores):raise ValueError("Report concurrency exceeds90% budget")
    cpus=[g[0] for g in cores[-wanted:]][::-1];old=os.sched_getaffinity(0)
    output.mkdir(parents=True);os.sched_setaffinity(0,{coordinator})
    jobs=[dict(run=x,output=output,allow_incomplete=a.allow_incomplete,requested=a.snapshots_ps,dpi=a.dpi,
        memory_limit=budget["memory_limit_bytes"]//wanted) for x in runs]
    print(json.dumps(dict(status="CASE_REPORT_BATCH",cases=len(runs),workers=wanted,worker_cpus=cpus,
        coordinator_cpu=coordinator,output_dir=str(output),cross_run_comparison=False)),flush=True)
    try:
        if wanted==1:
            os.sched_setaffinity(0,{cpus[0]})
            results=[worker(job) for job in jobs]
        else:
            context=mp.get_context("spawn");counter=context.Value("i",0)
            with ProcessPoolExecutor(max_workers=wanted,mp_context=context,initializer=initialize_worker,initargs=(cpus,counter)) as pool:
                results=list(pool.map(worker,jobs))
        batch=dict(schema=SCHEMA+".batch",cross_run_comparison=False,
            reports=[{k:x[k] for k in ("run_dir","report_path","report_sha256","provisional")} for x in results],resources=budget,actual_workers=wanted,worker_cpus=cpus,coordinator_cpu=coordinator,model_simulations_launched=0)
        (output/"case_reports_manifest.json").write_text(json.dumps(batch,indent=2,ensure_ascii=False)+"\n")
        return batch
    finally:os.sched_setaffinity(0,old)


if __name__=="__main__":main()
