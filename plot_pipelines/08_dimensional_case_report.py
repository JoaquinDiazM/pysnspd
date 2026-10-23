"""Single-case v1.2 report from committed saved states; never evolves a model."""
from __future__ import annotations
import argparse, hashlib, importlib.util, json, os, re, textwrap, time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS","NUMBA_NUM_THREADS"):
    os.environ[key]="1"
os.environ.setdefault("MPLCONFIGDIR","/home/jdiaz/scratch/tmp/pysnspd/matplotlib")
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import h5py
import numpy as np
from scipy.constants import Boltzmann as kB, elementary_charge as e
from scipy.special import expit
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import TwoSlopeNorm, BoundaryNorm, ListedColormap
from matplotlib.tri import Triangulation
from matplotlib.patches import Rectangle, FancyBboxPatch
from cycler import cycler
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":12,"axes.labelsize":12,
    "axes.titlesize":12,"xtick.labelsize":10.5,"ytick.labelsize":10.5,"legend.fontsize":9,
    "figure.facecolor":"#f5f7fb","axes.facecolor":"white","axes.edgecolor":"#c0cbd8",
    "axes.spines.top":False,"axes.spines.right":False,"grid.color":"#b8c5d2",
    "axes.prop_cycle":cycler(color=["#146b83","#e58b36","#7562a2","#a6ad3d","#c34d65","#3c957c"])})
INK="#17354a"
ACCENT="#45c5b0"
def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda:stream.read(2**20),b""):h.update(block)
    return h.hexdigest()
def scratch(path):
    p=Path(path).expanduser().resolve()
    if not p.is_relative_to(Path("/home/jdiaz/scratch")):raise ValueError("Products must stay in scratch")
    return p
def slots(values,count=6):
    a=[float(x) for x in list(values)[:count]]
    if not a or any(not np.isfinite(x) or x<0 for x in a):raise ValueError("Finite nonnegative snapshot times required")
    return (a+[a[-1]]*count)[:count]
def nearest(saved,requested):
    t=np.asarray(saved)
    if len(t)==0:return []
    return [int(np.argmin(abs(t-v))) for v in requested]
def number(value,scale=1.):
    return "not saved" if value is None else f"{float(value)*scale:.6g}"
def legacy_module():
    spec=importlib.util.spec_from_file_location("individual_case_helpers",ROOT/"plot_pipelines/07_case_report.py")
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module

def first_voltage_crossing(trace,baseline_V=0.,threshold_mV=.1):
    """Timestamp a positive threshold crossing using accepted voltage samples only."""
    t=np.asarray(trace["time_ps"],float);v=np.asarray(trace["Vout_V"],float)-baseline_V
    if not len(t):return None
    if np.any(~np.isfinite(t)) or np.any(~np.isfinite(v)) or np.any(np.diff(t)<0):
        raise ValueError("Finite ordered accepted electrical samples required")
    found=np.flatnonzero(v>=threshold_mV*1e-3)
    if not len(found):return None
    j=int(found[0]);level=threshold_mV*1e-3
    if j==0:return float(t[j])
    return float(t[j-1]+(t[j]-t[j-1])*(level-v[j-1])/(v[j]-v[j-1]))

def event_snapshot_plan(q,m,requested,legacy):
    """Nine per-case 2D slots; scan only committed complex-gap snapshots.

    VS/HB are first SAVED detections, not exact continuous onset times.
    Missing events repeat the final saved 2D state and remain explicitly absent.
    Scalar trigger timing can interpolate accepted Vout; map fields never do.
    """
    count=m["archive_counts"]["2d"]["snapshot_time_ps"]
    saved=np.asarray(q["snapshot_time_ps"][:count],float)
    if not len(saved) or np.any(~np.isfinite(saved)) or np.any(np.diff(saved)<=0):
        raise ValueError("Finite increasing committed 2D snapshot times required")
    end=float(saved[-1]);xy=q["coordinates_m"][:];tri=q["triangles"][:]
    edges=q["edges"][:];dc=q["dc_delta_bar"][:]
    events={name:dict(time_ps=None,status="NOT_OBSERVED_IN_SAVED_2D",snapshot_index=None)
        for name in ("VS","HB")}
    indeterminate_before_vs=0
    for j,t in enumerate(saved):
        if all(events[key]["time_ps"] is not None for key in ("VS","HB")):break
        d=q["delta_bar"][j]
        if events["VS"]["time_ps"] is None:
            vortex=legacy.phase_vorticity(d,dc,xy,tri)
            if vortex["positive"]+vortex["negative"]:
                events["VS"].update(time_ps=float(t),status="FIRST_SAVED_DETECTION",snapshot_index=j,
                    previous_saved_time_ps=None if j==0 else float(saved[j-1]),
                    positive=vortex["positive"],negative=vortex["negative"],
                    indeterminate_cells_at_detection=vortex["indeterminate_cells"])
            else:indeterminate_before_vs+=vortex["indeterminate_cells"]
        if events["HB"]["time_ps"] is None and legacy.hotbelts(d[None,:],dc,edges,xy)[0]:
            events["HB"].update(time_ps=float(t),status="FIRST_SAVED_DETECTION",snapshot_index=j,
                previous_saved_time_ps=None if j==0 else float(saved[j-1]))
    events["VS"].update(definition="First saved resolved P1 core with nonzero signed CCW phase winding; not curl(J)",
        contour_floor_relative_DC=1e-6,indeterminate_cell_samples_before_detection=indeterminate_before_vs,
        onset_between_saved_states_unknown=True)
    events["HB"].update(definition="Connected saved-node path joining both transverse edges with abs(Delta)<=0.2*abs(DeltaDC)",
        threshold_relative_DC=.2,onset_between_saved_states_unknown=True)
    obs=m.get("observation",{});baseline=float(obs.get("config",{}).get("baseline_V",0.))
    crossing=first_voltage_crossing(traces(q,m,"2d"),baseline,.1)
    events["Trigger"]=dict(time_ps=crossing,
        status="ACCEPTED_TRACE_CROSSING" if crossing is not None else "NOT_REACHED_IN_SAVED_2D_TRACE",
        threshold_mV=.1,polarity="positive",baseline_V=baseline,
        timestamp_method="Linear interpolation between accepted Vout samples; no field interpolation",
        controller_crossing_time_ps=None if obs.get("crossing_time_s") is None else float(obs["crossing_time_s"])*1e12,
        controller_detected=bool(obs.get("detected",False)))
    events["Final 2D"]=dict(time_ps=end,status="LAST_COMMITTED_2D_SNAPSHOT")
    default=[("0 ps",0.,None),("1 ps",1.,None),
        ("VS",events["VS"]["time_ps"],"VS"),
        ("VS + 0.1 ps",None if events["VS"]["time_ps"] is None else events["VS"]["time_ps"]+.1,"VS"),
        ("VS + 0.2 ps",None if events["VS"]["time_ps"] is None else events["VS"]["time_ps"]+.2,"VS"),
        ("HB",events["HB"]["time_ps"],"HB"),
        ("HB + 2 ps",None if events["HB"]["time_ps"] is None else events["HB"]["time_ps"]+2.,"HB"),
        ("Trigger",crossing,"Trigger"),("Final 2D",end,"Final 2D")]
    entries=default if requested is None else [(f"User {i+1}",v,None) for i,v in enumerate(slots(requested,9))]
    plans=[];used=[]
    for label,want,event in entries:
        fallback=want is None;target=end if fallback else float(want)
        j=int(np.argmin(abs(saved-target)));actual=float(saved[j])
        plans.append(dict(slot=len(plans)+1,label=label,event=event,
            event_status=None if event is None else events[event]["status"],
            requested_ps=None if fallback else target,selection_target_ps=target,saved_ps=actual,snapshot_index=j,
            offset_ps=None if fallback else actual-target,
            outside_saved_window=bool(target<saved[0] or target>saved[-1]),
            repeated_saved_state=j in used,fallback="last committed 2D state; event not observed" if fallback else None))
        used.append(j)
    cadence=np.diff(saved);queue=[]
    if count<9:queue.append("Save at least nine 2D states for nine distinct panels when physically available.")
    if len(set(used))<len(set(p["selection_target_ps"] for p in plans)):
        queue.append("Stable finer 2D snapshot cadence is needed to separate requested event offsets; no fields were interpolated.")
    policy=dict(mode="per-case events" if requested is None else "explicit user times",
        slots=9,ignored_extra_requests=0 if requested is None else max(0,len(requested)-9),
        missing_event_policy="Repeat last committed 2D state; label event absent; never invent an onset",
        out_of_window_policy="Nearest endpoint; retain original request and flag truncation",
        tie_policy="Earlier saved state",repetition_policy="Keep repetitions and semantic/request order",
        timing_scope="VS/HB only sampled at saved states; transient or unresolved vortices can be missed",
        available_states=count,unique_selected_states=len(set(used)),
        saved_cadence_ps=dict(minimum=None if not len(cadence) else float(cadence.min()),
            median=None if not len(cadence) else float(np.median(cadence)),maximum=None if not len(cadence) else float(cadence.max())),
        future_visualization_queue=queue)
    return plans,events,policy

def snapshot_panel_title(plan):
    if plan["fallback"]:return f'{plan["label"]}: not observed\nsaved {plan["saved_ps"]:.5g} ps (fallback)'
    return f'{plan["label"]}\ntarget {plan["requested_ps"]:.5g} / saved {plan["saved_ps"]:.5g} ps'
class Report:
    def __init__(self,pdf,m,run):
        self.pdf,self.m,self.run=pdf,m,run
        self.pages=[]
    def page(self,title,rows=2,cols=2):
        title=re.sub(r"^\d+\.\s*","",title)
        fig,axes=plt.subplots(rows,cols,figsize=(11.69,8.27),squeeze=False)
        fig.patches.extend([Rectangle((0,.885),1,.115,transform=fig.transFigure,
            facecolor=INK,edgecolor="none",zorder=10),
            Rectangle((.035,.913),.007,.055,transform=fig.transFigure,facecolor=ACCENT,
                edgecolor="none",zorder=11)])
        fig.text(.055,.951,title,fontsize=17,weight="bold",color="white",zorder=12)
        p=self.m["parameter_identity"]
        fig.text(.055,.908,f"{self.run.name}    /    {p['wavelength_nm']:g} nm    /    {p['current_uA']:g} uA",
            fontsize=10,color="#cbdbe6",zorder=12)
        fig.text(.96,.91,"MODEL "+self.m["model_version"],fontsize=9,color=ACCENT,
            ha="right",weight="bold",zorder=12)
        return fig,axes.ravel()
    def save(self,fig,title,note="",layout=True):
        title=re.sub(r"^\d+\.\s*","",title)
        for a in fig.axes:
            label=a.get_ylabel()
            if len(label)>29 and "\n" not in label:a.set_ylabel("\n".join(textwrap.wrap(label,23)))
        if layout:fig.subplots_adjust(left=.095,right=.95,bottom=.22,top=.825,hspace=.51,wspace=.43)
        if note:
            fig.text(.055,.12,"\n".join(textwrap.wrap(note,162)),fontsize=8.5,va="top",color="#4e6072",linespacing=1.25)
        fig.lines.append(plt.Line2D([.055,.95],[.045,.045],transform=fig.transFigure,color="#d0dae5",lw=.8))
        fig.text(.055,.022,"PYSNSPD  /  INDIVIDUAL CASE  /  SAVED DATA ONLY",fontsize=8,color="#718397")
        fig.text(.95,.022,f"{len(self.pages)+1:02d}",ha="right",fontsize=10,color=INK,weight="bold")
        self.pdf.savefig(fig,dpi=150);plt.close(fig);self.pages.append(title)
    def events(self,ax):
        obs=self.m.get("observation",{})
        for v,color in [(obs.get("crossing_time_s"),"#c77c11"),(obs.get("stop_target_s"),"#4b639a")]:
            if v is not None:ax.axvline(v*1e12,color=color,ls=":",lw=1.)
        cold=self.m.get("cold_recovery") or {}
        if cold.get("active"):ax.axvline(cold["activation_time_ps"],color="#168b69",ls="--")
    def finish_axes(self,axes):
        for ax in axes:
            ax.set_xlabel("Time after postcascade transfer (ps)")
            ax.grid(alpha=.3);self.events(ax)
            handles,_=ax.get_legend_handles_labels()
            if handles:ax.legend(frameon=True,facecolor="white",edgecolor="#d9e1e9")

def card_table(fig,rect,title,rows):
    ax=fig.add_axes(rect);ax.set_axis_off()
    ax.add_patch(FancyBboxPatch((0,0),1,1,boxstyle="round,pad=.006,rounding_size=.018",
        transform=ax.transAxes,facecolor="white",edgecolor="#d8e2ec",lw=.8,clip_on=False))
    ax.add_patch(Rectangle((0,.89),1,.11,transform=ax.transAxes,facecolor=INK,edgecolor="none"))
    ax.text(.035,.945,title.upper(),color="white",weight="bold",fontsize=10.5,va="center")
    processed=[(textwrap.wrap(str(k),28) or [""],textwrap.wrap(str(v),25) or [""]) for k,v in rows]
    weights=[max(len(k),len(v))+.55 for k,v in processed];total=sum(weights);y=.875
    for i,((keys,values),weight) in enumerate(zip(processed,weights)):
        height=.85*weight/total
        if i%2==0:ax.add_patch(Rectangle((.012,y-height),.976,height,transform=ax.transAxes,
            facecolor="#f0f5f9",edgecolor="none"))
        ax.text(.035,y-height/2,"\n".join(keys),fontsize=9.7,color="#43576b",va="center",linespacing=1.05)
        ax.text(.575,y-height/2,"\n".join(values),fontsize=10.2,color=INK,va="center",linespacing=1.05)
        y-=height
    return ax

def traces(group,m,name):
    n=m["archive_counts"][name]["trace"]
    if len(group["trace"])<n:raise ValueError("Archive shorter than committed trace")
    return group["trace"][:n]
def phonon_data(run,p,Eu):
    from pysnspd.experimental.effective_material import EffectivePhononMaterial
    from pysnspd.experimental.effective_bcs_kinetics import acoustic_phonon_grid
    material=EffectivePhononMaterial.load(run/"setup2d/material.npz")
    grid=acoustic_phonon_grid(material,Eu,nodes=p["phonon_nodes"])
    return grid
def moment_temperature(target,omega,weights,Eu):
    lo=np.zeros_like(target);hi=np.full_like(target,100.)
    def moment(T):
        x=omega[None,:]*Eu/(kB*np.maximum(T[:,None],1e-300))
        with np.errstate(over="ignore",under="ignore"):
            return (np.exp(-x)/(-np.expm1(-x)))@weights
    for _ in range(12):
        need=moment(hi)<target
        if not np.any(need):break
        hi[need]*=2
    else:raise ArithmeticError("Phonon moment inversion lacks a bracket")
    for _ in range(48):
        mid=(lo+hi)/2;low=moment(mid)<target
        lo=np.where(low,mid,lo);hi=np.where(low,hi,mid)
    return np.where(target==0,0.,(lo+hi)/2)
def snapshot_data(q,indices,m,grid):
    ids=sorted(set(indices));p=m["physicalparams"];Eu=p["energy_scale_J"]
    if not np.allclose(q["phonon_energies_bar"][:],grid.energies,rtol=1e-13,atol=0):raise ValueError("Phonon quadrature mismatch")
    weights=p["N0_single_spin_per_J_m3"]*Eu**2*grid.capacities*grid.energies
    out={}
    for j in ids:
        d=q["delta_bar"][j];occ=q["p"][j];ph=q["n"][j]
        en=np.sqrt(q["xi_bar"][:][None,:]**2+abs(d[:,None])**2)
        ew=4*p["N0_single_spin_per_J_m3"]*Eu**2*q["xi_weights"][:]*en
        fd=expit(-en*Eu/(kB*p["bath_temperature_K"]))
        n0=q["dc_n"][:]
        qp_abs=np.sum(ew*abs(occ-fd),axis=1)
        ph_abs=abs(ph-n0)@weights
        out[j]=dict(d=d,Te=q["Te_K"][j],Tph=moment_temperature(ph@weights,grid.energies,weights,Eu),
            u=q["u_excess_J_m3"][j],qp_abs=qp_abs,ph_abs=ph_abs,
            qp=np.sum(ew*occ,axis=1),ph=ph@weights,p=occ,n=ph,en=en)
    return out

def total_display_vectors(mesh,flux,thickness):
    """Weighted display fit using only actual positive-area Voronoi faces."""
    edges=mesh["edges"];xy=mesh["coordinates_m"];tri=mesh["triangles"]
    areas=mesh["dual_face_length_m"]*thickness
    active=areas>0
    if np.max(abs(flux[~active]),initial=0.)>1e-18+1e-10*np.max(abs(flux),initial=0.):
        raise ArithmeticError("Finite total flux on a zero-area dual face")
    density=np.divide(flux,areas,out=np.zeros_like(flux),where=active)
    unit=(xy[edges[:,1]]-xy[edges[:,0]])/mesh["edge_length_m"][:,None]
    lookup={tuple(sorted((int(a),int(b)))):i for i,(a,b) in enumerate(edges)}
    ids=np.array([[lookup[tuple(sorted((int(a),int(b))))] for a,b in
        ((t[0],t[1]),(t[1],t[2]),(t[2],t[0]))] for t in tri])
    u=unit[ids];weights=areas[ids]/np.maximum(areas[ids].max(axis=1)[:,None],1e-300)
    matrix=np.einsum("tki,tk,tkj->tij",u,weights,u)
    rhs=np.einsum("tki,tk,tk->ti",u,weights,density[ids])
    good=np.linalg.det(matrix)>1e-14*np.trace(matrix,axis1=1,axis2=2)**2
    if not np.all(good):raise ValueError("Total-current display triangle lacks two independent positive-area faces")
    vectors=np.linalg.solve(matrix,rhs[...,None])[...,0]
    residual=np.einsum("tki,ti->tk",u,vectors)-density[ids]
    rms=float(np.sqrt(np.sum(weights*residual**2)/weights.sum()))
    scale=float(np.sqrt(np.sum(weights*density[ids]**2)/weights.sum()))
    return vectors,dict(projection_RMS_A_m2=rms,relative_projection_RMS=None if scale==0 else rms/scale,
        zero_area_faces_excluded=int(np.count_nonzero(~active)),unresolved_triangles=0)

def total_currents(run,q,selected,data,m,tr,legacy):
    # Static algebraic reconstruction at saved states; no DC solve or kinetic RHS.
    paths=[run/"setup2d/meshes/v12_boundary_repair/mesh.npz"]
    paths+=sorted((run/"setup2d/meshes/reference_geometry").glob("*/mesh.npz"))
    path=next((x for x in paths if x.exists()),None)
    if path is None:raise FileNotFoundError("Saved mesh needed for total-current maps")
    with np.load(path,allow_pickle=False) as z:mesh={k:z[k].copy() for k in z.files}
    if not np.allclose(mesh["coordinates_m"],q["coordinates_m"][:],rtol=0,atol=1e-20):raise ValueError("Mesh coordinates mismatch")
    from pysnspd.experimental.thermal_spatial_usadel import ThermalGraph
    from pysnspd.experimental.effective_bcs_functional import EffectiveBCSFunctional
    from pysnspd.experimental.thermal_stiffness import ThermalStiffnessFunctional
    from pysnspd.experimental.resummed_thermal import ResummedThermalFunctional
    from pysnspd.experimental.effective_bcs_ports import EffectiveNeutralPorts
    source_files=("effective_bcs_functional.py","effective_bcs_ports.py","thermal_spatial_usadel.py",
        "thermal_stiffness.py","general_thermal_gradient.py","finite_depairing_sums.py","resummed_thermal.py")
    for fn in source_files:
        rel="pysnspd/experimental/"+fn
        expected=m.get("source_pins_sha256",{}).get(rel)
        if expected is None or digest(ROOT/rel)!=expected:raise ValueError("Unidentified reconstruction source: "+rel)
    pars=m["physicalparams"];ident=m["parameter_identity"]
    graph=ThermalGraph(mesh["area_weights"],mesh["edges"],mesh["conductance"],mesh["coordinates_bar"],mesh["fixed_nodes"])
    cls=(ResummedThermalFunctional if ident.get("gradient_closure")=="thermal-depairing"
         else ThermalStiffnessFunctional if ident.get("stiffness_temperature")=="electron" else EffectiveBCSFunctional)
    f=cls(graph,energy_scale_J=pars["energy_scale_J"],N0_per_J_m3=pars["N0_single_spin_per_J_m3"],
        ell0_m=pars["ell0_m"],thickness_m=pars["film_thickness_m"],count_nodes=q["xi_bar"][:],
        count_weights=q["xi_weights"][:],gradient_law="bcs_radial",bath_temperature_K=pars["bath_temperature_K"],
        diffusion_m2_s=pars["diffusion_m2_s"])
    if not np.allclose(f.volumes_m3,q["volumes_m3"][:],rtol=1e-10,atol=0):raise ValueError("Mesh volume mismatch")
    ports=EffectiveNeutralPorts(graph,pars["sheet_resistance_ohm"],mesh["left_nodes"],mesh["right_nodes"])
    result={};checks=[]
    for j in sorted(set(selected)):
        t=float(q["snapshot_time_ps"][j]);k=int(np.argmin(abs(tr["time_ps"]-t)))
        if abs(tr["time_ps"][k]-t)>1e-8:raise ValueError("No accepted current state at saved snapshot")
        a=data[j];el=ports.solve(legacy.edge_supercurrent(f,a["d"],a["Te"]),float(tr["Is_A"][k]))
        err=abs(el.Vdev_V-tr["Vdev_V"][k])
        if err>1e-11+1e-8*max(abs(el.Vdev_V),abs(tr["Vdev_V"][k]),1e-8):raise ArithmeticError("Saved/reconstructed Vdev mismatch")
        vectors,fit=total_display_vectors(mesh,el.I_total_A,pars["film_thickness_m"])
        result[j]=vectors
        checks.append(dict(snapshot_index=j,time_ps=t,Vdev_error_V=float(err),display_fit=fit))
    return result,checks,path

def overview(r):
    m=r.m;p=m["parameter_identity"];v=m["physicalparams"];s=m.get("source",{});obs=m.get("observation",{})
    fig,axes=r.page("Run dossier - physical inputs")
    for ax in axes:ax.remove()
    status=m.get("status","not saved")
    status_text="TIME BUDGET / CHECKPOINT SAVED" if status=="WALL_BUDGET_CHECKPOINTED" else status
    fig.text(.06,.843,status_text,fontsize=11,weight="bold",color="#a46a18")
    fig.text(.95,.843,f"Stored endpoint  {number(m.get('time_ps'))} ps",ha="right",fontsize=11,color=INK)
    circ=v["circuit"];phase=m.get("circuit_by_phase",{});c1=phase.get("1d",circ)
    material=json.loads((r.run/"setup2d/material.json").read_text())
    sections=[
      ("Device and photon",[
        ("Bias / wavelength",f"{p['current_uA']:g} uA / {p['wavelength_nm']:g} nm"),
        ("Photon / retained energy",f"{number(s.get('photon_energy_J'),1/e)} / {number(s.get('retained_energy_J'),1/e)} eV"),
        ("Retained fraction eta",p["retained_fraction"]),("Deposition sigma",f"{p['sigma_nm']:g} nm"),
        ("Preparation profile","Reflecting Gaussian" if s.get("profile")=="reflecting_gaussian_images" else s.get("profile","Original Gaussian")),
        ("Width / thickness",f"80 / {v['film_thickness_m']*1e9:g} nm"),
        ("Domain 2D / 1D",f"{p['length_nm']:g} / {p.get('length_1d_nm',p['length_nm']):g} nm"),
        ("Deposition x / y",f"{number(s.get('center_m',[0,0])[0],1e9)} / {number(s.get('center_m',[0,0])[1],1e9)} nm")]),
      ("Material and kinetics",[
        ("Tc / bath",f"{v['Tc_K']:g} / {v['bath_temperature_K']:g} K"),
        ("Diffusion D",f"{v['diffusion_m2_s']*1e4:.6g} cm2/s"),("Sheet resistance",f"{p['sheet_resistance_ohm']:g} ohm/sq"),
        ("Single-spin N0",f"{v['N0_single_spin_per_J_m3']:.5g} /(J m3)"),
        ("Length scale ell0",f"{v['ell0_m']*1e9:.5g} nm"),
        ("tau_ee(Tc) / tau_ep(Tc)",f"{p['tau_ee_Tc_ps']:g} / {p['tau_ep_Tc_ps']:g} ps"),
        ("Escape / tau0",f"{p['tau_escape_ps']:g} ps / {v['tau0_s']*1e9:.5g} ns"),
        ("Target depairing current",f"{p.get('target_depairing_uA',0):.6g} uA"),
        ("Gradient / stiffness",f"{p.get('gradient_closure','not saved')} / {p.get('stiffness_temperature','not saved')}"),
        ("Sound speed / Debye T",f"{material.get('sound_speed_m_per_s',0):g} m/s / {material.get('Debye_temperature_K',0):g} K")]),
      ("Circuit",[
        ("Rb / RL",f"{circ['R_bias_ohm']:g} / {circ['R_load_ohm']:g} ohm"),
        ("Lb / Cc",f"{circ['L_bias_H']*1e6:g} uH / {circ['C_couple_F']*1e12:g} pF"),
        ("Source voltage",f"{circ['V_bias_V']:.6g} V"),
        ("Lext, 2D / 1D",f"{circ['Lk_ext_H']*1e9:.6g} / {c1['Lk_ext_H']*1e9:.6g} nH"),
        ("Readout","Vout = RL (Ib - Is)"),("Instrument response","No amplifier / filter")]),
      ("Observation and checkpoint",[
        ("Trigger / confirmation",f"+{p['threshold_mV']:g} mV / {p['confirmation_ps']:g} ps"),
        ("Crossing / handoff",f"{number(obs.get('crossing_time_s'),1e12)} / {number(obs.get('stop_target_s'),1e12)} ps"),
        ("Post-crossing margin",f"{p['extra_time_ps']:g} ps"),
        ("No-trigger / final target",f"{p['ceiling_ps']:g} ps / "+("detection + margin" if p.get("latency_only") else f"{p['final_time_ns']:g} ns")),
        ("Cold approximation active",str(bool((m.get("cold_recovery") or {}).get("active")))),
        ("1D session / budget",f"{number(m.get('runtime_1d_session_wall_s'),1/3600)} / {number(m.get('runtime_1d_wall_budget_hours'))} h")])]
    rects=[(.055,.495,.43,.315),(.515,.495,.43,.315),(.055,.155,.43,.315),(.515,.155,.43,.315)]
    for rect,(title,rows) in zip(rects,sections):card_table(fig,rect,title,rows)
    r.save(fig,"Physical input tables","Checkpointed by a wall budget means resumable incomplete data, not physical recovery. Time zero is the postcascade transfer. Diffusion and thermal relaxation inputs are effective calibrators, not independent measurements.",layout=False)
    fig,axes=r.page("Run dossier - numerical controls and guards")
    for ax in axes:ax.remove()
    cfg=m.get("monitor_config",{});latest=m.get("latest_monitor",{})
    sections=[
      ("Discretization and native state",[
        ("Spacing, 2D / 1D",f"{p['spacing_nm']:g} / {p['spacing_1d_nm']:g} nm"),
        ("xi cutoff / order",f"{p['xi_cutoff']:g} / {p['xi_order']:g}"),
        ("Electron / phonon bins",f"{v.get('diffusion_calibration',{}).get('quadrature',{}).get('electronic_states','see receipt')} / {p['phonon_nodes']:g}"),
        ("Transport quadrature",p["transport_order"]),("Initial / max 1D step",f"{p['dt_ps']:g} / {p['maximum_1d_step_ps']:g} ps"),
        ("Minimum adaptive step",f"{p['minimum_adaptive_step_ps']:g} ps"),
        ("Snapshots, 2D / 1D",f"{p['snapshot_every_ps']:g} / {p['snapshot_1d_every_ps']:g} ps")]),
      ("Solver and energy acceptance",[
        ("Temporal rtol / atol",f"{p['time_rtol']:g} / {p['time_atol']:g}"),
        ("Nonlinear tolerance",f"{p['step_tolerance']:g}"),("Maximum iterations",p["max_iterations"]),
        ("Local / global energy gate",f"{p['local_energy_tolerance']:g} / {p['global_energy_tolerance']:g}"),
        ("Max observed local / global",f"{number(m.get('max_local_energy_fraction'))} / {number(m.get('max_global_energy_fraction'))}"),
        ("Accepted / rejected steps",f"{m.get('accepted_steps')} / {m.get('rejections')}"),
        ("Final signed energy defect",f"{number(m.get('energy_defect_J'),1/e)} eV")]),
      ("Finite-domain diagnostics",[
        ("Outer free collar",f"{cfg.get('collar_nm',25):g} nm"),
        ("Thermal warning / persistence",f"{cfg.get('collar_temperature_rise_K',.5):g} K / {cfg.get('collar_confirmation_ps',2):g} ps"),
        ("Gap / collar energy margins",f"{cfg.get('collar_gap_fraction',.9):g} / {cfg.get('collar_energy_fraction',.001):g}"),
        ("Absolute cut / activity gate",f"{p['guard_cut_energy_relative']:g}"),
        ("Final cut / activity",f"{number(m.get('domain_fraction'),100)} %"),
        ("Energy basis","Photon + absolute device work" if m["model_version"] in ("1.2.3","1.2.4") else "Retained photon"),
        ("Gap reference","Current-carrying cold bulk" if m["model_version"] in ("1.2.3","1.2.4") else "Initial DC")]),
      ("Electrical and overlap guards",[
        ("Current / device voltage",f"{cfg.get('current_limit_uA',0):.5g} uA / {cfg.get('device_voltage_limit_mV',0):.5g} mV"),
        ("Readout / capacitor voltage",f"{cfg.get('readout_voltage_limit_mV',0):.5g} / {cfg.get('capacitor_voltage_limit_mV',0):.5g} mV"),
        ("Overlap amplitude / current",f"{p['guard_gap_relative']:g} / {p['guard_shadow_current_relative']:g}"),
        ("Overlap voltage / energy",f"{p['guard_shadow_voltage_V']*1e6:g} uV / {p['guard_shadow_energy_relative']:g}"),
        ("Projection / reconciliation",f"{p['guard_projection_energy_relative']:g} / {p['guard_reconciliation_energy_relative']:g}"),
        ("Workers per simulated case",f"{p['threads']} threads / {p['physical_cores']} cores touched"),
        ("Latest boundary alarm",str(latest.get("boundary_alarm")))])]
    rects=[(.055,.51,.43,.325),(.515,.51,.43,.325),(.055,.16,.43,.325),(.515,.16,.43,.325)]
    for rect,(title,rows) in zip(rects,sections):card_table(fig,rect,title,rows)
    r.save(fig,"Numerical control tables","All exact parameter identities, active monitor values and selected snapshot provenance are retained in the JSON receipt. Engineering guards do not establish experimental calibration or domain convergence. No saved state is changed by this report.",layout=False)

def electrical(r,official):
    fig,ax=r.page("2. Laboratory observables and circuit")
    for name,tr in official:
        t=tr["time_ps"]
        ax[0].plot(t,tr["Vout_V"]*1e3,label="Vout / "+name)
        ax[1].plot(t,tr["Vdev_V"]*1e3,label="Vdevice / "+name)
        for key,style in [("Ib_A","-"),("Is_A","--")]:ax[2].plot(t,tr[key]*1e6,style,label=key[:-2]+" / "+name)
        ax[3].plot(t,(tr["Ib_A"]-tr["Is_A"])*1e6,label="Iread / "+name)
    ax[0].axhline(r.m["parameter_identity"]["threshold_mV"],color=".5",ls="--",label="Trigger")
    for a,y in zip(ax,("Vout, unfiltered (mV)","Device voltage (mV)","Branch currents (uA)","Readout current (uA)")):a.set_ylabel(y)
    r.finish_axes(ax)
    r.save(fig,"Laboratory observables","Orange dotted: first voltage crossing. Blue dotted: official 2D-to-1D handoff. The post-handoff signal is also plotted in ps; no ns axis hides the terminal dynamics. Negative Vout alone can be capacitive undershoot.")
def balance(r,official):
    ref=r.m["source"]["retained_energy_J"];fig,ax=r.page("3. Thermal power and independent energy balance")
    for name,tr in official:
        t=tr["time_ps"]
        ax[0].plot(t,tr["Is_A"]*tr["Vdev_V"]*1e6,label="Device work / "+name)
        ax[0].plot(t,tr["escape_power_W"]*1e6,ls="--",label="Phonon escape / "+name)
        ax[0].plot(t,tr["cuts_power_W"]*1e6,ls=":",label="Cuts, signed outward / "+name)
        # The 1D representation has its own reference offset: compare changes
        # inside each phase; never silently join different absolute energies.
        ax[1].plot(t,(tr["energy_J"]-tr["energy_J"][0])/e,label="Stored delta / "+name)
        ax[1].plot(t,(tr["external_J"]-tr["external_J"][0])/e,ls="--",label="External delta / "+name)
        ax[2].plot(t,100*abs(tr["energy_defect_J"])/(ref+tr["absolute_external_J"]),label=name)
        ax[3].plot(t,100*tr["domain_fraction"],label=name)
    ax[2].axhline(100*r.m["parameter_identity"]["global_energy_tolerance"],color="r",ls="--")
    ax[3].axhline(100*r.m["parameter_identity"]["guard_cut_energy_relative"],color="r",ls="--")
    for a,y in zip(ax,("Power (uW)","Phase-relative energy change (eV)","Independent defect / activity (%)","Absolute cut energy / activity (%)")):a.set_ylabel(y)
    r.finish_axes(ax)
    r.save(fig,"Thermal power and independent balance","Native signed cut power includes reservoir/boundary work and is not a direct escaped-photon fraction. Cumulative absolute cut energy is independently normalized by retained photon energy plus absolute device work. Energy curves use each phase's first accepted state, avoiding a false representation jump.")
def solver(r,official):
    fig,ax=r.page("4. Numerical resolution and recorded solver effort")
    for name,tr in official:
        t=tr["time_ps"];positive=tr["step_ps"]>0
        ax[0].semilogy(t[positive],tr["step_ps"][positive],label=name)
        ax[1].plot(t,tr["iterations"],label=name)
        ax[2].plot(t,tr["temporal_error"],label=name)
        ax[3].plot(t,tr["energy_defect_J"]/e,label=name)
    ax[2].axhline(1.,ls="--",color="r",label="Temporal acceptance")
    for a,y in zip(ax,("Accepted step (ps)","Recorded iterations / stage count","Step-doubling normalized error","Signed independent energy defect (eV)")):a.set_ylabel(y)
    r.finish_axes(ax)
    r.save(fig,"Numerical resolution","2D trace stores no step-doubling error (zeros are not a measured convergence result). Cold DOP853 stores RHS evaluation count instead of native nonlinear iterations. Per-step wall time, rejected candidate history and capacity maps were not archived; no CPU speed curve is inferred.")

def maps_2d(r,q,data,selected,plan,flow,legacy):
    if len(selected)!=9 or len(plan)!=9:raise ValueError("Every 2D mosaic requires nine slots")
    xy=q["coordinates_m"][:]*1e9
    xy[:,0]-=r.m.get("source",{}).get("center_m",[np.mean(xy[:,0])*1e-9,0.])[0]*1e9
    tri=q["triangles"][:];mesh=Triangulation(xy[:,0],xy[:,1],tri)
    saved=q["snapshot_time_ps"][:]
    Eu=q.attrs["energy_scale_J"];meV=e*1e-3
    definitions=[("gap","5. 2D condensate snapshots","Gap magnitude (meV)","viridis"),
        ("Te","6. 2D electron energy-equivalent temperature","Te (K)","inferno"),
        ("Tph","7. 2D phonon energy-equivalent temperature","Tph (K)","YlOrBr"),
        ("u","8. 2D signed energy density relative to DC","Excess energy density (kJ/m3)","PuOr_r"),
        ("total","9. 2D total current density","Total current magnitude (1e10 A/m2)","cividis"),
        ("winding","10. 2D condensate phase-winding vorticity","Signed CCW phase winding","PiYG")]
    vortex={}
    for j in set(selected):vortex[j]=legacy.phase_vorticity(data[j]["d"],q["dc_delta_bar"][:],q["coordinates_m"][:],tri)
    for kind,title,label,cmap in definitions:
        fig,axes=r.page(title,3,3)
        values={}
        for j in set(selected):
            a=data[j]
            values[j]=(abs(a["d"])*Eu/meV if kind=="gap" else a[kind] if kind in ("Te","Tph")
                else a["u"]*1e-3 if kind=="u" else np.linalg.norm(flow[j],axis=1)/1e10 if kind=="total" else vortex[j]["winding"])
        lo=min(float(np.min(a)) for a in values.values());hi=max(float(np.max(a)) for a in values.values())
        if kind in ("gap","total"):lo=0.
        if kind=="u":hi=max(abs(lo),abs(hi),1e-15);lo=-hi
        if hi<=lo:hi=lo+max(abs(lo)*1e-6,1e-12)
        for ax,j,rec in zip(axes,selected,plan):
            value=values[j]
            if kind=="winding":
                d=vortex[j];palette=ListedColormap(["#008176","#ffffff","#ba287f"]);palette.set_bad("#b9b9b9")
                norm=BoundaryNorm([-1.5,-.5,.5,1.5],3)
                ax.tripcolor(mesh,facecolors=np.zeros(len(tri)),cmap=ListedColormap(["#b9b9b9"]),vmin=0,vmax=1,rasterized=True)
                image=ax.tripcolor(mesh,facecolors=np.ma.array(value,mask=~d["valid"]),cmap=palette,norm=norm,rasterized=True)
                ax.text(.02,.02,f"+: {d['positive']}  -: {d['negative']}",transform=ax.transAxes,fontsize=10)
            elif kind=="total":
                image=ax.tripcolor(mesh,facecolors=value,cmap=cmap,vmin=lo,vmax=hi,rasterized=True)
                centers=xy[tri].mean(axis=1);idx=np.arange(0,len(tri),max(1,len(tri)//65));v=flow[j][idx]
                norm=np.linalg.norm(v,axis=1);keep=norm>max(norm.max()*1e-4,1e-30)
                v=v[keep]/norm[keep,None]
                ax.quiver(centers[idx[keep],0],centers[idx[keep],1],6*v[:,0],6*v[:,1],
                    angles="xy",scale_units="xy",scale=1,width=.003,color="black")
            elif kind=="gap":
                # Interpolate the complex P1 order parameter before taking its
                # modulus: nodal |Delta| interpolation can hide vortex zeros.
                image=legacy.draw_complex_gap(ax,mesh,data[j]["d"]*Eu,unit_J=meV,vmin=lo,vmax=hi,
                    bounds=(xy[:,0].min(),xy[:,0].max(),xy[:,1].min(),xy[:,1].max()),samples=301,cmap=cmap)
            else:image=ax.tripcolor(mesh,value,shading="gouraud",cmap=cmap,vmin=lo,vmax=hi,rasterized=True)
            ax.set(title=snapshot_panel_title(rec),xlabel="x - xdep (nm)",ylabel="y (nm)",aspect="equal")
            ax.set_title(snapshot_panel_title(rec),fontsize=10.5,pad=7)
            ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5));ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3))
        fig.subplots_adjust(left=.075,right=.96,bottom=.27,top=.80,hspace=1.1,wspace=.45)
        cb=fig.colorbar(image,cax=fig.add_axes([.22,.18,.58,.018]),orientation="horizontal",ticks=[-1,0,1] if kind=="winding" else None)
        cb.set_label(label,fontsize=12)
        note=("Complex P1 gap magnitude, not interpolation of nodal magnitudes." if kind=="gap" else
            "Energy-equivalent moment of native nonthermal occupations; not an imposed thermodynamic temperature." if kind in ("Te","Tph") else
            "Native BCS internal energy, thermal gradient and phonons minus the phase DC reference." if kind=="u" else
            "Algebraic reconstruction at the exact saved state. Dual-flux least-squares display vectors; arrows show direction with fixed length. No separate normal/supercurrent snapshot pages." if kind=="total" else
            "Phase topology, not curl(J). Gray cells have an unresolved or vanishing-gap contour; no vortex trajectories are inferred between snapshots.")
        r.save(fig,title,note+" Nine nearest saved states; repeats retained. VS/HB are first saved detections; missing events use the labelled final-state fallback. No field interpolation. Shared color range.",layout=False)
    return [dict(snapshot_index=j,time_ps=float(saved[j]),positive=vortex[j]["positive"],negative=vortex[j]["negative"],
        indeterminate_cells=vortex[j]["indeterminate_cells"]) for j in sorted(set(selected))]

def profile_plan(f,m,requested):
    """Select the official representation at the requested physical time."""
    t2=f["2d"]["snapshot_time_ps"][:m["archive_counts"]["2d"]["snapshot_time_ps"]]
    t1=f["1d"]["snapshot_time_ps"][:m["archive_counts"]["1d"]["snapshot_time_ps"]]
    if not len(t2) and not len(t1):raise ValueError("No committed snapshots")
    hand=float(t1[0]) if len(t1) else float(t2[-1])
    end=float(t1[-1]) if len(t1) else float(t2[-1])
    wants=slots(requested if requested is not None else
        [0.,max(0.,hand-2),hand,hand+(end-hand)/3,hand+2*(end-hand)/3,end])
    plan=[]
    for want in wants:
        name="1d" if len(t1) and (want>=hand-1e-9 or not len(t2)) else "2d"
        saved=t1 if name=="1d" else t2
        j=int(np.argmin(abs(saved-want)))
        plan.append(dict(requested_ps=want,saved_ps=float(saved[j]),snapshot_index=j,phase=name,
            reconstruction="native official 1D" if name=="1d" else "2D-derived volume-slab average; display only",
            reconstructed_domain_nm=float(m["parameter_identity"]["length_1d_nm" if name=="1d" else "length_nm"])))
    return wants,plan

def reconstruct_profiles(f,m,plan,data1,data2):
    """Read-only scalar projection; never creates an evolved 1D state."""
    q2=f["2d"];xy=q2["coordinates_m"][:];xmin=float(xy[:,0].min());xmax=float(xy[:,0].max())
    q1=f["1d"];xx=q1["coordinates_m"][:,0];center1=.5*(xx.min()+xx.max())
    pars=m["physicalparams"];Eu=pars["energy_scale_J"];meV=e*1e-3
    remap=None;target=None;target_vol=None
    if any(v["phase"]=="2d" for v in plan):
        from pysnspd.experimental.dimensional_transition import volume_restriction
        source="pysnspd/experimental/dimensional_transition.py"
        if digest(ROOT/source)!=m["source_pins_sha256"][source]:raise ValueError("Projection source changed")
        count=max(3,int(round((xmax-xmin)/(m["parameter_identity"]["spacing_1d_nm"]*1e-9)))+1)
        target=np.linspace(xmin,xmax,count);edges=np.r_[xmin,.5*(target[:-1]+target[1:]),xmax]
        target_vol=80e-9*pars["film_thickness_m"]*np.diff(edges)
        remap=volume_restriction(xy[:,0],q2["volumes_m3"][:],target_vol)
        np.testing.assert_allclose(remap.sum(axis=1),1.,rtol=1e-10,atol=1e-12)
    profiles=[];receipts=[]
    for rec in plan:
        j=rec["snapshot_index"];phase=rec["phase"]
        if phase=="1d":
            a=data1[j];x=(xx-center1)*1e9
            values=dict(gap=abs(a["d"])*Eu/meV,Te=a["Te"],Tph=a["Tph"],u=a["u"]*1e-3)
            error=None
        else:
            a=data2[j];x=(target-.5*(xmin+xmax))*1e9
            fields=dict(gap=abs(a["d"])*Eu/meV,Te=a["Te"],Tph=a["Tph"],u=a["u"]*1e-3)
            values={key:remap@val for key,val in fields.items()}
            old_energy=float(q2["volumes_m3"][:]@a["u"])
            new_energy=float(target_vol@(values["u"]*1e3))
            error=abs(new_energy-old_energy)/max(m["source"]["retained_energy_J"],abs(old_energy))
            if error>1e-10:raise ArithmeticError("Display projection lost volume-integrated energy")
        profiles.append(dict(x_nm=x,values=values,**rec))
        receipts.append(dict(rec,energy_integral_relative_error=error))
    return profiles,receipts

def curve_page(r,q,profiles):
    fig,axes=r.page("Longitudinal snapshots - 2D-derived and native 1D")
    colors=["#166d89","#d68a31","#8b5f9e","#3c927f","#c34d65","#738d30"]
    x=q["coordinates_m"][:,0]*1e9;origin=.5*(x[0]+x[-1]);x-=origin
    cfg=r.m.get("monitor_config",{})
    for a,color in zip(profiles,colors):
        derived=a["phase"]=="2d";style="--" if derived else "-"
        label=f"{a['saved_ps']:.5g} ps / "+("2D avg" if derived else "1D")
        for ax,key in zip(axes,("gap","Te","Tph","u")):
            ax.plot(a["x_nm"],a["values"][key],c=color,ls=style,lw=1.7,label=label)
    for ax,y in zip(axes,("Gap magnitude (meV)","Te (K)","Tph (K)","Signed excess energy density (kJ/m3)")):
        ax.set(xlabel="x - absorption center (nm)",ylabel=y);ax.grid(alpha=.3)
        for attr in ("core_begin_m","core_end_m"):
            if attr in q.attrs:ax.axvline(q.attrs[attr]*1e9-origin,color="#98a8b8",ls=":",lw=1)
        ax.axvspan(x[0],x[0]+cfg.get("collar_nm",25.),color="#d86d6d",alpha=.07)
        ax.axvspan(x[-1]-cfg.get("collar_nm",25.),x[-1],color="#d86d6d",alpha=.07)
    fig.subplots_adjust(left=.095,right=.95,bottom=.285,top=.825,hspace=.48,wspace=.43)
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,ncol=3,loc="upper center",bbox_to_anchor=(.53,.205),frameon=False,fontsize=9)
    r.save(fig,"Longitudinal snapshots",
        "Dashed: volume-weighted scalar averages of saved 2D fields, restricted to the 300 nm lens. Solid: saved native 1D fields. No time interpolation, exterior filling or 1D evolution is inferred from a 2D-derived curve. Averaging |gap| and Te is a plotting diagnostic, not averaging the complex field or certifying a reduced thermodynamic state.",layout=False)

def longitudinal(r,q,data,profiles):
    p=r.m["parameter_identity"];cfg=r.m.get("monitor_config",{})
    count=r.m["archive_counts"]["1d"]["snapshot_time_ps"];tt=q["snapshot_time_ps"][:count]
    x=q["coordinates_m"][:,0]*1e9
    origin=.5*(x[0]+x[-1]);x-=origin
    gapDC=abs(q["dc_delta_bar"][:]);Eu=q.attrs["energy_scale_J"]
    fig,axes=r.page("11. Official 1D space-time evolution")
    for ax,key,label,cmap in zip(axes,("gap","Te","Tph","u"),("Gap magnitude (meV)","Te (K)","Tph (K)","Excess energy density (kJ/m3)"),("viridis","inferno","YlOrBr","PuOr_r")):
        values=np.array([abs(data[j]["d"])*Eu/(e*1e-3) if key=="gap" else data[j][key] for j in range(count)])
        if key=="u":values*=1e-3
        kwargs={}
        if key=="u":
            v=max(float(abs(values).max()),1e-15);kwargs=dict(vmin=-v,vmax=v)
        im=ax.pcolormesh(x,tt,values,shading="nearest",cmap=cmap,rasterized=True,**kwargs)
        for attr in ("core_begin_m","core_end_m"):
            if attr in q.attrs:ax.axvline(q.attrs[attr]*1e9-origin,color=".6",ls=":",lw=1)
        collar=cfg.get("collar_nm",25.)
        for bound in [x[0]+collar,x[-1]-collar]:ax.axvline(bound,color="#43bac0",ls="--",lw=1.)
        ax.set(xlabel="x - strip center (nm)",ylabel="Stored time (ps)",ylim=(tt[0],tt[-1]),title=label)
        cb=fig.colorbar(im,ax=ax);cb.ax.tick_params(labelsize=9)
        cold=r.m.get("cold_recovery") or {}
        if cold.get("active"):ax.axhline(cold["activation_time_ps"],color="white",ls="--")
    r.save(fig,"1D space-time evolution","Dotted vertical lines: former 300 nm lens interfaces. Cyan dashed: free outer-collar inner edges. Only stored snapshots are shown; block extent does not establish intermediate front motion. After a labelled cold activation, fields are equilibrium reconstructions.")
    curve_page(r,q,profiles)
    vol=q["volumes_m3"][:];dist=np.minimum(x-x[0],x[-1]-x)
    collar=(dist>1e-7)&(dist<=cfg.get("collar_nm",25.)+1e-7)
    masks=[collar&(x<0),collar&(x>0)]
    if any(not np.any(mask) for mask in masks):raise ValueError("No free left/right collar")
    ref=r.m["source"]["retained_energy_J"]
    bath=r.m["physicalparams"]["bath_temperature_K"]
    metrics=dict(time_ps=tt.tolist(),left_Te=[],right_Te=[],left_Tph=[],right_Tph=[],left_gap=[],right_gap=[],
        left_energy=[],right_energy=[],total_energy=[],qp_energy=[],ph_energy=[],hot_length_nm=[],
        left_gap_current=[],right_gap_current=[],activity_energy=[],left_activity_energy=[],right_activity_energy=[])
    for j in range(count):
        a=data[j];total=(a["qp_abs"]+a["ph_abs"])*vol
        for mask,side in zip(masks,("left","right")):
            metrics[side+"_Te"].append(float(max(a["Te"][mask])))
            metrics[side+"_Tph"].append(float(max(a["Tph"][mask])))
            metrics[side+"_gap"].append(float(min(abs(a["d"][mask])/gapDC[mask])))
            metrics[side+"_energy"].append(float(total[mask].sum()/ref))
        native_trace=traces(q,r.m,"1d")
        row=int(np.argmin(abs(native_trace["time_ps"]-tt[j])))
        if abs(native_trace["time_ps"][row]-tt[j])>1e-8:raise ValueError("No accepted work ledger at snapshot")
        activity=ref+native_trace["device_absolute_J"][row]
        # The cold continuation pins the terminal amplitude to the same
        # instantaneous stable bulk branch used by the simulator's guard.
        matching=float(np.mean(abs(a["d"][[0,-1]])/gapDC[[0,-1]]))
        for mask,side in zip(masks,("left","right")):
            metrics[side+"_gap_current"].append(float(min(abs(a["d"][mask])/(gapDC[mask]*matching))))
            metrics[side+"_activity_energy"].append(float(total[mask].sum()/activity))
        metrics["activity_energy"].append(float(total[collar].sum()/activity))
        metrics["total_energy"].append(float(total[collar].sum()/ref))
        metrics["qp_energy"].append(float((a["qp_abs"][collar]*vol[collar]).sum()/ref))
        metrics["ph_energy"].append(float((a["ph_abs"][collar]*vol[collar]).sum()/ref))
        # Exact saved-node FV hot length, matching the simulator's discrete guard.
        mask=abs(a["d"])<.2*gapDC;lengths=vol/(80e-9*r.m["physicalparams"]["film_thickness_m"])*1e9
        current=longest=0.
        for hot,dx in zip(mask,lengths):
            current=current+dx if hot else 0.;longest=max(longest,current)
        metrics["hot_length_nm"].append(float(longest))
    fig,axes=r.page("13. Boundary-collar diagnostics and stop criteria")
    for side,color in [("left","#4575b4"),("right","#d6604d")]:
        axes[0].plot(tt,metrics[side+"_Te"],c=color,label="Te / "+side)
        axes[0].plot(tt,metrics[side+"_Tph"],c=color,ls="--",label="Tph / "+side)
        current_basis=r.m.get("model_version")=="1.2.3"
        axes[1].plot(tt,metrics[side+("_gap_current" if current_basis else "_gap")],c=color,label=side)
        axes[2].plot(tt,100*np.array(metrics[side+("_activity_energy" if current_basis else "_energy")]),c=color,label=side)
    axes[0].axhline(bath+cfg.get("collar_temperature_rise_K",.5),color="r",ls=":",label="Thermal warning")
    axes[1].axhline(cfg.get("collar_gap_fraction",.9),color="r",ls="--",label="Suppression threshold")
    current_basis=r.m.get("model_version")=="1.2.3"
    axes[2].plot(tt,100*np.array(metrics["activity_energy" if current_basis else "total_energy"]),color="black",label="Both collars / active basis")
    if current_basis:
        axes[2].plot(tt,100*np.array(metrics["total_energy"]),color=".5",ls=":",label="Photon only / diagnostic")
    axes[2].axhline(100*cfg.get("collar_energy_fraction",.001),color="r",ls="--",label="Energy threshold")
    axes[3].plot(tt,metrics["hot_length_nm"],label="Longest |gap| < 0.2 DC segment")
    for a,y in zip(axes,("Maximum free-collar temperature (K)",
        "Free-collar |gap| / current bulk" if current_basis else "Minimum free-collar |gap| / DC",
        "Unsigned collar energy fraction (%)" if current_basis else "Unsigned collar energy / retained photon (%)",
        "Discrete suppressed length (nm)")):a.set_ylabel(y)
    r.finish_axes(axes)
    latest=r.m.get("latest_monitor",{})
    note=(f"Saved-snapshot reconstruction only; {cfg.get('collar_confirmation_ps',2):g} ps persistence is evaluated by the simulator at every accepted step. "
        "Thermal stop = (electron OR phonon thermal warning) AND significant unsigned native energy; gap suppression is an independent OR condition. "
        f"Final accepted monitor: Te={number(latest.get('collar_max_Te_K'))} K, gap/DC={number(latest.get('collar_min_gap_ratio'))}, collar energy={number(latest.get('collar_native_energy_fraction'),100)}%.")
    if current_basis:
        note+=" Active energy basis: retained photon plus cumulative absolute device work; photon-only curve is diagnostic. Gap reference: instantaneous cold bulk, recovered from the prescribed cut amplitude."
    r.save(fig,"Boundary diagnostics",note)
    return metrics
def spectra(r,q,data,grid):
    count=r.m["archive_counts"]["1d"]["snapshot_time_ps"];j=count-1;a=data[j];x=q["coordinates_m"][:,0]*1e9
    cfg=r.m.get("monitor_config",{});dist=np.minimum(x-x[0],x[-1]-x)
    collar=(dist>1e-7)&(dist<=cfg.get("collar_nm",25)+1e-7)
    center=int(np.argmin(abs(x-x.mean())));left=int(np.flatnonzero(collar&(x<x.mean()))[0]);right=int(np.flatnonzero(collar&(x>x.mean()))[-1])
    ids=[center,left,right];labels=["Center","Left free collar","Right free collar"]
    Eu=q.attrs["energy_scale_J"];meV=e*1e-3
    fig,ax=r.page("14. Native occupations at the last stored 1D state")
    for k,label in zip(ids,labels):
        en=a["en"][k]*Eu/meV
        ax[0].semilogy(en,np.maximum(a["p"][k],1e-300),label=label)
        fd=expit(-a["en"][k]*Eu/(kB*a["Te"][k]))
        ax[1].plot(en,a["p"][k]-fd,label=label)
        ax[2].semilogy(grid.energies*Eu/meV,np.maximum(a["n"][k],1e-300),label=label)
        T=a["Tph"][k];z=grid.energies*Eu/(kB*T)
        bose=np.exp(-z)/(-np.expm1(-z))
        weights=r.m["physicalparams"]["N0_single_spin_per_J_m3"]*Eu**2*grid.capacities*grid.energies
        ax[3].plot(grid.energies*Eu/meV,(a["n"][k]-bose)*weights,label=label)
    for a,xlab,ylab in zip(ax,("Quasiparticle energy (meV)","Quasiparticle energy (meV)","Phonon energy (meV)","Phonon energy (meV)"),
        ("Occupation p","p - FD(Te), energy-matched","Occupation n","Energy per bin: n - Bose(Tph) (J/m3)")):
        a.set(xlabel=xlab,ylabel=ylab);a.grid(alpha=.2);a.legend()
    for a in (ax[0],ax[2]):a.set_ylim(bottom=1e-16,top=max(1.,a.get_ylim()[1]))
    r.save(fig,"Last native occupations",f"Stored time {q['snapshot_time_ps'][j]:.6g} ps. Te/Tph match an energy moment, not the full spectrum. Log panels hide occupations below 1e-16 for display only; saved states are never clipped. Reconstructed fields after cold activation must not be interpreted as late quasiparticle kinetics.")
def overlap(r):
    checks=r.m.get("shadow_checks",[])
    if not checks:return
    fig,ax=r.page("15. Independent 2D-to-1D overlap",2,3)
    keys=("gap_amplitude_relative","supercurrent_relative","voltage_error_V","film_energy_relative","populations_energy_relative","current_relative")
    p=r.m["parameter_identity"]
    limits=[p["guard_gap_relative"],p["guard_gap_relative"],p["guard_shadow_voltage_V"],p["guard_shadow_energy_relative"],p["guard_populations_energy_relative"],p["guard_shadow_current_relative"]]
    labels=("Gap-amplitude discrepancy (%)","Supercurrent discrepancy (%)","Vout discrepancy (uV)","Film-energy discrepancy (%)","Unsigned population-energy discrepancy (%)","Branch-current discrepancy (%)")
    for a,key,limit,label in zip(ax,keys,limits,labels):
        fac=1e6 if key=="voltage_error_V" else 100
        a.plot([v["time_ps"] for v in checks],[v.get(key,np.nan)*fac for v in checks])
        a.axhline(limit*fac,color="r",ls="--");a.set(xlabel="Stored overlap time (ps)",ylabel=label);a.grid(alpha=.2)
    r.save(fig,"Independent overlap","Only independent overlap comparisons recorded by the runner are shown. Red: declared admission margins. Complex-gap phase discrepancy remains a gauge-sensitive diagnostic and does not replace the amplitude/current tests. Projection alone does not establish a valid long-domain approximation.")
def summary(r,legacy):
    fig,ax=r.page("16. Model, circuit and approximation scope",2,1)
    legacy.circuit_diagram(ax[0])
    for label in ax[0].texts:
        if label.get_text()=="2D strip":label.set_text("Resolved strip")
    ax[1].axis("off")
    circ=r.m["physicalparams"]["circuit"];phases=r.m.get("circuit_by_phase",{})
    ext1=phases.get("1d",circ)["Lk_ext_H"]
    fig.text(.095,.51,f"Rb={circ['R_bias_ohm']:g} ohm; RL={circ['R_load_ohm']:g} ohm; Lb={circ['L_bias_H']*1e6:g} uH; Cc={circ['C_couple_F']*1e12:g} pF; Lext(2D/1D)={circ['Lk_ext_H']*1e9:.5g}/{ext1*1e9:.5g} nH",fontsize=10)
    equations=[
        r"$\partial_t\Delta={\cal K}_{\rm KWT}[\Delta,p,T_e,\phi],\quad \partial_t p={\cal C}_e+{\cal T}_e+{\cal H}_e,\quad \partial_t n={\cal C}_{ph}-(n-n_B)/\tau_{\rm esc}$",
        r"$\nabla\cdot({\bf J}_s+{\bf J}_n)=0,\quad I_{s,e}=-\frac{2e}{\hbar}\,\partial_{\alpha_e}G,\quad {\bf J}_n=-\sigma_n\nabla\phi$",
        r"$V_{\rm out}=R_L(I_b-I_s),\quad U_{\rm film}=U_{\rm BCS}+G-\sum_i T_{e,i}\,\partial G/\partial T_{e,i}+U_{ph}$",
        r"$L_f(I_s)=\frac{\hbar L_1}{2e}\frac{dq_b}{dI_s},\quad \dot I_s=\frac{v_c+R_L(I_b-I_s)}{L_{\rm ext}+L_f(I_s)}$"]
    for i,eq in enumerate(equations):ax[1].text(.01,.95-i*.23,eq,fontsize=12,va="top")
    cold=r.m.get("cold_recovery") or {};policy=r.m.get("approximation_policy",{})
    r.save(fig,"Model and approximation scope","Neutral: quasistatic electrical continuity. Effective BCS/KWT: shared BCS energy, calibrated thermal-gradient closure and conserving native kinetic equations; not exact microscopic Usadel. 1D retains all native spectra until optional cold admission. Cold active: "+str(bool(cold.get("active")))+
        f"; earliest qualification {policy.get('minimum_time_ps',100):g} ps, independent overlap {policy.get('qualification_ps',10):g} ps. Instantaneous-current cold outer cuts are an exterior assumption; domain guards remain active.")

def report(job):
    run,out,requested2,requested1=job;started=time.monotonic()
    run=Path(run).resolve();out=scratch(out)
    manifest_path=run/"manifest.json";raw=manifest_path.read_bytes();m=json.loads(raw)
    if m.get("model_version") not in ("1.2.0","1.2.1","1.2.2","1.2.3","1.2.4"):raise ValueError("A supported v1.2 phase-aware archive is required")
    if m.get("status") in ("RUNNING","PREPARATION"):raise ValueError("Use a terminal committed case, not a live archive")
    if digest(run/"checkpoint_latest.npz")!=m.get("checkpoint_sha256"):raise ValueError("Committed checkpoint checksum mismatch")
    target=out/("case_report_"+run.name+".pdf");receipt=target.with_suffix(".json")
    if target.exists() or receipt.exists():raise FileExistsError("Existing report preserved; choose a fresh output directory")
    build=out/".build"/run.name;build.mkdir(parents=True,exist_ok=False)
    temporary=build/"case_report.pdf";pdf=PdfPages(temporary,metadata={"Title":run.name+" / v"+m["model_version"],"Author":"pysnspd"})
    result=dict(schema="pysnspd.dimensional_case_report.v4",run_dir=str(run),model_version=m["model_version"],
        status=m.get("status"),stop_reason=m.get("monitor_reason") or m.get("error"),model_simulations_launched=0,
        kinetic_rhs_evaluations=0,new_dc_solves=0,temporal_interpolation=False,
        field_temporal_interpolation=False,voltage_event_timestamp_interpolation=True,
        manifest_sha256=hashlib.sha256(raw).hexdigest(),checkpoint_sha256=m["checkpoint_sha256"],
        parameter_identity=m["parameter_identity"],monitor_config=m.get("monitor_config"),latest_monitor=m.get("latest_monitor"),
        pipeline_sha256=digest(__file__),excluded_snapshot_pages=["normal current","supercurrent"])
    r=Report(pdf,m,run);legacy=legacy_module()
    try:
        with h5py.File(run/"data.h5","r",swmr=True) as f:
            counts=m["archive_counts"]
            for name in ("2d","shadow1d","1d"):
                for key,n in counts[name].items():
                    if len(f[name][key])<n:raise ValueError("Archive shorter than committed checkpoint")
            official=[(name,traces(f[name],m,name)) for name in ("2d","1d") if counts[name]["trace"]]
            if not official:raise ValueError("No committed official trajectory")
            overview(r);electrical(r,official);balance(r,official);solver(r,official)
            grid=phonon_data(run,m["parameter_identity"],m["physicalparams"]["energy_scale_J"])
            wants1,plan1=profile_plan(f,m,requested1)
            data2={};data1={}
            if counts["2d"]["snapshot_time_ps"]:
                q=f["2d"];saved=q["snapshot_time_ps"][:counts["2d"]["snapshot_time_ps"]]
                plan2,events,policy=event_snapshot_plan(q,m,requested2,legacy)
                selected=[rec["snapshot_index"] for rec in plan2]
                extra=[rec["snapshot_index"] for rec in plan1 if rec["phase"]=="2d"]
                data2=snapshot_data(q,selected+extra,m,grid)
                flow,checks,mesh=total_currents(run,q,selected,data2,m,traces(q,m,"2d"),legacy)
                result["vorticity"]=maps_2d(r,q,data2,selected,plan2,flow,legacy)
                result["current_reconstruction_checks"]=checks;result["mesh_sha256"]=digest(mesh)
                result["snapshots_2d"]=plan2
                result["snapshot_events_2d"]=events;result["snapshot_policy_2d"]=policy
                result["event_detector_source_sha256"]=digest(ROOT/"plot_pipelines/07_case_report.py")
            if counts["1d"]["snapshot_time_ps"]:
                q=f["1d"];n=counts["1d"]["snapshot_time_ps"];saved=q["snapshot_time_ps"][:n]
                data1=snapshot_data(q,range(n),m,grid)
                profiles,selection=reconstruct_profiles(f,m,plan1,data1,data2)
                result["collar_snapshots"]=longitudinal(r,q,data1,profiles)
                spectra(r,q,data1,grid)
            else:
                profiles,selection=reconstruct_profiles(f,m,plan1,data1,data2)
                curve_page(r,f["1d"],profiles)
            result["snapshots_1d"]=selection
            result["longitudinal_profile_policy"]="official 1D after handoff; volume-slab scalar averages of actual 2D data before handoff"
            overlap(r);summary(r,legacy)
            result["archive_counts"]=counts;result["material_sha256"]=digest(run/"setup2d/material.npz")
        if manifest_path.read_bytes()!=raw:raise RuntimeError("Case changed during report; output not promoted")
    finally:pdf.close()
    result.update(pdf=str(target),pages=r.pages,pdf_sha256=digest(temporary),elapsed_wall_s=time.monotonic()-started,
        caveats=["Collar temporal curves are reconstructed only at saved snapshots; accepted-step persistence is not recreated",
            "VS/HB denote first saved detections, not exact physical onset; short-lived or unresolved vortices may be missed",
            "No instrument response; no full-reset or convergence admission inferred",
            "Saved 1D current continuity implies uniform longitudinal total current Is/(width*thickness)",
            "Post-cold-activation maps/spectra are equilibrium reconstructions, not native late kinetic evolution"])
    note=build/"receipt.json";note.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    os.replace(temporary,target);os.replace(note,receipt)
    print(json.dumps(dict(status="REPORT_COMPLETE",case=run.name,pages=len(r.pages),pdf=str(target),elapsed_s=result["elapsed_wall_s"])),flush=True)
    return str(target)
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run-dir",nargs="+",type=Path,required=True)
    p.add_argument("--output-dir",type=Path,required=True)
    p.add_argument("--snapshots-2d-ps","--snapshots-ps",dest="snapshots_2d_ps",nargs="+",type=float,
        help="Nine 2D slots: first nine requests, repeat last if fewer. Default per case: 0,1,VS,VS+.1,VS+.2,HB,HB+2,Trigger,final2D. Nearest saved fields; missing event repeats final2D. --snapshots-ps remains an alias")
    p.add_argument("--snapshots-1d-ps",nargs="+",type=float,help="Six curve slots; early times use saved 2D-derived scalar averages, later times use native 1D")
    p.add_argument("--workers",type=int,choices=(1,2),default=1)
    a=p.parse_args(argv);out=scratch(a.output_dir)
    for vals,count in [(a.snapshots_2d_ps,9),(a.snapshots_1d_ps,6)]:
        if vals is not None:slots(vals,count)
    if len({x.resolve().name for x in a.run_dir})!=len(a.run_dir):raise ValueError("Duplicate report names; use separate destinations")
    active=[str(x) for x in a.run_dir if json.loads((x/"manifest.json").read_text()).get("status") in ("RUNNING","PREPARATION")]
    if active:raise ValueError("No reports started: these cases are still active; rerun after they stop: "+", ".join(active))
    out.mkdir(parents=True,exist_ok=True)
    jobs=[(str(x),str(out),a.snapshots_2d_ps,a.snapshots_1d_ps) for x in a.run_dir]
    print(json.dumps(dict(status="SAVED_DATA_ONLY",cases=len(jobs),workers=min(a.workers,len(jobs)),
        numba_threads_per_worker=1,blas_threads_per_worker=1,simulations=0)),flush=True)
    if a.workers==1:
        for j in jobs:report(j)
    else:
        with ProcessPoolExecutor(max_workers=min(a.workers,len(jobs))) as pool:list(pool.map(report,jobs))
if __name__=="__main__":main()
