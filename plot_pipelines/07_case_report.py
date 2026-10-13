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
from pysnspd.plotting.style import apply_thesis_style
from pysnspd.plotting.complex_gap import draw_complex_gap
from pysnspd.experimental.effective_bcs_functional import EffectiveBCSFunctional, bcs_count_quadrature
from pysnspd.experimental.effective_bcs_ports import EffectiveNeutralPorts
from pysnspd.experimental.thermal_spatial_usadel import ThermalGraph
from pysnspd.experimental.trigger_window import TriggerWindow

SCRATCH = Path("/home/jdiaz/scratch")
EV = 1.602176634e-19
MEV = EV*1e-3
SCHEMA = "pysnspd.individual_effective_case_report.v1"
RECON_SOURCES = ("pysnspd/experimental/effective_bcs_functional.py",
    "pysnspd/experimental/effective_bcs_ports.py","pysnspd/experimental/thermal_spatial_usadel.py")
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
    if value is None: return "No asignado"
    v=float(value)*scale
    return f"{v:.{digits}g}" if np.isfinite(v) else "infinito"


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
    if not path.is_file(): return {"status":"Sin auditoría terminal autenticada"}
    record=json.loads(path.read_text()).get("cases",{}).get(run.name,{})
    audit=record.get("audit",{}); p=Path(audit.get("path",""))
    if record.get("manifest_sha256")!=manifest_sha or not p.is_file() or sha256(p)!=audit.get("sha256"):
        return {"status":"Auditoría ausente o de otra versión"}
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
    for k in RECON_SOURCES:
        if sha256(ROOT/k)!=m["source_pins_sha256"].get(k): raise ValueError("Saved reconstruction source differs: "+k)
    return dict(run=run,m=m,z=trace,s=maps,mesh=mesh,dc=dc,dc_receipt=dr,obs=obs,provisional=not complete,
        inputs=inputs,memory=memory,execution=execution,audit=existing_audit(run,inputs["manifest.json"]) if complete else {"status":"Provisional; sin auditoría terminal inferida"})


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


def select_snapshots(times,obs,hotbelt,requested=None):
    t=np.asarray(times)*1e12
    if requested is not None:
        q=np.asarray(requested,float)
        if not len(q) or len(q)>24 or np.any(~np.isfinite(q)) or np.any(q<t[0]) or np.any(q>t[-1]+1e-9):
            raise ValueError("At most24 requested times inside the saved snapshot window")
        targets=[t[0],*q,t[-1]]
    elif obs.get("trigger_confirmed",obs.get("detected",False)):
        cross=obs.get("first_crossing_time_s",obs.get("crossing_time_s"))
        conf=obs.get("confirmation_time_s")
        hb=t[np.flatnonzero(hotbelt)[0]] if np.any(hotbelt) else .5*t[-1]
        targets=[t[0],min(1.,t[-1]),hb,cross*1e12,conf*1e12,t[-1]]
    else:
        targets=[t[0],min(1.,t[-1]),min(10.,t[-1]),.25*t[-1],.5*t[-1],t[-1]]
    return np.asarray(sorted({0,len(t)-1,*[int(np.argmin(abs(t-q))) for q in targets]}),int)


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


def edge_supercurrent(f,d):
    """Exact gauge derivative of the saved gradient; independent of p and n."""
    tail,head=f.graph.edges.T
    kap=f.rigidity(abs(d)**2)[0]
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
    graph=ThermalGraph(mesh["area_weights"],mesh["edges"],mesh["conductance"],mesh["coordinates_bar"],mesh["fixed_nodes"])
    xi,w=bcs_count_quadrature(m["parameter_identity"]["xi_cutoff"],m["parameter_identity"]["xi_order"])
    f=EffectiveBCSFunctional(graph,energy_scale_J=p["energy_scale_J"],N0_per_J_m3=p["N0_single_spin_per_J_m3"],
        ell0_m=m["geometry"]["ell0_m"],thickness_m=p["film_thickness_m"],count_nodes=xi,count_weights=w,
        gradient_law="bcs_radial",bath_temperature_K=p["bath_temperature_K"],diffusion_m2_s=p["diffusion_m2_s"])
    ports=EffectiveNeutralPorts(graph,p["sheet_resistance_ohm"],c["dc"]["left_nodes"],c["dc"]["right_nodes"])
    results={}; checks=[]
    for i in selected:
        ts=c["s"]["time_s"][i]
        exact=np.flatnonzero(abs(z["time_s"]-ts)<=1e-23)
        if not len(exact): raise ValueError("No accepted electrical state at snapshot; temporal interpolation refused")
        j=int(exact[np.argmin(abs(z["time_s"][exact]-ts))]); Is=float(z["Is_A"][j])
        el=ports.solve(edge_supercurrent(f,c["s"]["delta_J"][i]/f.Eu_J),Is)
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
            display_vector_fits=fits))
    return results,checks


class Report:
    def __init__(self,path,c,dpi):
        self.c=c; self.dpi=dpi; self.pages=[]
        self.pdf=PdfPages(path,metadata={"Title":"SNSPD - "+c["run"].name,"Author":"pysnspd",
            "Subject":"Individual effective saved case; no new evolution or inter-run comparison"})
    def page(self,title):
        fig=plt.figure(figsize=(11.69,8.27))
        fig.text(.055,.945,title,fontsize=18,weight="bold",color=INK)
        a=self.c["m"]["parameter_identity"]
        subtitle=f'{self.c["run"].name} | {a["wavelength_nm"]:g} nm | I0={a["current_uA"]:.7g} µA'
        if self.c["provisional"]: subtitle+=" | PROVISIONAL - "+self.c["m"]["status"]
        fig.text(.055,.912,subtitle,fontsize=10,color=ACCENT)
        return fig
    def save(self,fig,title,caption):
        fig.text(.055,.065,textwrap.fill(caption,155),fontsize=8.5,color=INK,va="bottom",linespacing=1.35)
        fig.text(.055,.025,f'pysnspd | caso individual | reloj postcascada | ventana 0–{self.c["metrics"]["final_ps"]:.7g} ps',fontsize=8,color=INK)
        fig.text(.945,.025,str(len(self.pages)+1),fontsize=9,ha="right",color=INK)
        with plt.rc_context({"savefig.bbox":None}):
            self.pdf.savefig(fig,dpi=self.dpi)
        plt.close(fig)
        self.pages.append(dict(page=len(self.pages)+1,title=title,caption=caption))


def table(fig,bounds,title,rows,font=9):
    ax=fig.add_axes(bounds); ax.axis("off")
    ax.text(0,1.04,title,fontsize=12,color=ACCENT,weight="bold",transform=ax.transAxes)
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
    fig=r.page("1. Parámetros del caso")
    table(fig,[.055,.525,.425,.34],"Geometría y material",[
        ("Ancho / longitud resuelta (nm)",fmt(g["width_m"],1e9)+" / "+fmt(g["length_m"],1e9)),
        ("Espesor físico / malla (nm)",fmt(p["film_thickness_m"],1e9)+" / "+fmt(a["spacing_nm"])),
        ("Tc / baño (K)",fmt(p["Tc_K"])+" / "+fmt(p["bath_temperature_K"])),
        ("Difusión D (m²/s)",fmt(p["diffusion_m2_s"])),("Resistencia por cuadrado (Ω)",fmt(p["sheet_resistance_ohm"])),
        ("DOS N0 monoespín (1/(J m³))",fmt(p["N0_single_spin_per_J_m3"])),
        ("Velocidad acústica / Debye (m/s; K)",fmt(mat.get("sound_speed_m_per_s"))+" / "+fmt(mat.get("Debye_temperature_K"))),
        ("γ acústica derivada",fmt(mat.get("gamma_A20"))),("Nodos / aristas",f'{g["n_nodes"]} / {g["n_edges"]}'),
        ("Modelo","Neutral efectivo BCS/KWT")])
    table(fig,[.52,.525,.425,.34],"Fotón y calibradores",[
        ("Longitud de onda (nm)",fmt(s.get("wavelength_nm",a["wavelength_nm"]))),
        ("Energía fotónica / retenida (eV)",fmt(s.get("photon_energy_J"),1/EV)+" / "+fmt(s.get("retained_energy_J"),1/EV)),
        ("Fracción retenida η",fmt(s.get("retained_fraction"))),("Anchura σ (nm)",fmt(s.get("sigma_nm"))),
        ("Depósito (x,y) en nm",str([round(x*1e9,4) for x in s.get("center_m",[])])),
        ("τee(Tc), libre (ps)",fmt(p["tau_ee_Tc_s"],1e12)),("τep(Tc), libre (ps)",fmt(p["tau_ep_Tc_s"],1e12)),
        ("τ0 e-ph, derivado (ns)",fmt(p["tau0_s"],1e9)),("Escape fonónico (ps)",fmt(p["tau_escape_s"],1e12)),
        ("Fuente inicial","Gap DC; FD/BE comunes")])
    table(fig,[.055,.16,.425,.295],"Circuito y discriminador ideal",[
        ("Corriente DC (µA)",fmt(p["current_A"],1e6)),("Rb / Rlectura (Ω)",fmt(circ["R_bias_ohm"])+" / "+fmt(circ["R_load_ohm"])),
        ("Lb / Lext congelada (nH)",fmt(circ["L_bias_H"],1e9)+" / "+fmt(circ["Lk_ext_H"],1e9)),
        ("Capacidad / Vbias (pF; V)",fmt(circ["C_couple_F"],1e12)+" / "+fmt(circ["V_bias_V"])),
        ("Umbral positivo (mV)",fmt(a["threshold_mV"])),("Confirmación / extra (ps)",fmt(a["confirmation_ps"])+" / "+fmt(a["extra_time_ps"])),
        ("Histéresis (fracción) / techo (ps)",fmt(cfg.get("hysteresis_fraction",.05))+" / "+fmt(a["ceiling_ps"])),
        ("Origen del reloj","Transferencia postcascada")])
    table(fig,[.52,.16,.425,.295],"Resolución de la corrida",[
        ("Estados e / modos ph",f'{m["kinetics"]["electron_states"]} / {m["kinetics"]["phonon_modes"]}'),
        ("Corte ξ / orden GL / transporte",f'{a["xi_cutoff"]:g} / {a["xi_order"]} / {a["transport_order"]}'),
        ("Cap de paso (ps)",fmt(a["dt_ps"])),("Tolerancia no lineal",fmt(a["step_tolerance"])),
        ("Guarda energética local",fmt(a["local_energy_tolerance"])),
        ("Snapshots / checkpoint (ps)",fmt(a["snapshot_every_ps"])+" / "+fmt(a["checkpoint_every_ps"])),
        ("Hilos de la simulación",str(a["threads"])),("Método temporal","Punto medio acoplado")])
    r.save(fig,"Parámetros","h de malla es resolución, no espesor. η es energía retenida condicional a absorción, no eficiencia óptica/PCR. "
        "τee y τep son calibradores libres; no constantes medidas ni rangos. La fuente postcascada y el material acústico son hipótesis exploratorias.")


def observables(r):
    c=r.c;z=c["z"];m=c["m"];v=c["metrics"];t=z["time_s"]*1e12
    baseline=c["obs"].get("baseline_Vout_V",c["obs"].get("config",{}).get("baseline_V",0.))
    voltage=z["Vout_V"]-baseline;threshold=m["parameter_identity"]["threshold_mV"]*1e-3
    scale,unit=(1e6,"µV") if max(abs(voltage))<.02*threshold else (1e3,"mV")
    fig=r.page("2. Observables: lectura y detección")
    ax=fig.add_axes([.075,.51,.39,.33]);bx=fig.add_axes([.565,.51,.37,.33])
    ax.plot(t,voltage*scale,color=ACCENT)
    if max(abs(voltage))>=.15*threshold:ax.axhline(threshold*scale,color="black",ls=":",label="Umbral")
    else:ax.text(.03,.92,f"Umbral={threshold*1e3:g} mV (fuera de escala)",transform=ax.transAxes,fontsize=9)
    ax.set(xlabel="Tiempo postcascada (ps)",ylabel="Vout menos basal ("+unit+")")
    bx.plot(t,(z["Ib_A"]-z["Is_A"])*1e6,color=ACCENT)
    bx.set(xlabel="Tiempo postcascada (ps)",ylabel="Corriente de lectura Ib − Is (µA)")
    for a in (ax,bx):
        a.grid(alpha=.25)
        if v["confirmed"]:
            a.axvline(v["crossing_ps"],color="tab:green",ls=":",label="Cruce")
            a.axvline(v["confirmation_ps"],color="tab:green",ls="--",label="Confirmación")
            a.legend(fontsize=8)
    status=("Provisional; "+c["m"]["status"]+("; cruce confirmado" if v["confirmed"] else "")) if c["provisional"] else ("Detección confirmada" if v["confirmed"] else ("Censurado sin detección" if v["censored"] else "Ventana de desarrollo/DC"))
    table(fig,[.075,.16,.86,.245],"Eventos y magnitudes guardadas",[
        ("Estado de observación",status),("Latencia relativa / cruce (ps)",fmt(v["latency_ps"],digits=9)),
        ("Confirmación aceptada (ps)",fmt(v["confirmation_ps"],digits=9)),
        ("Final de trayectoria / mapas (ps)",fmt(v["final_ps"],digits=9)+" / "+fmt(v["maps_final_ps"],digits=9)),
        ("Máximo Vout observado (mV)",fmt(v["max_observed_Vout_mV"],digits=8)),
        ("Pendiente local al cruce (mV/ps)",fmt(v["local_slope_mV_ps"],digits=7)),("Corriente Is mínima (µA)",fmt(v["min_Is_uA"],digits=8))])
    r.save(fig,"Lectura y detección","Vout=Rlectura(Ib−Is), menos basal DC: NO incluye filtro, amplificador, ruido ni ancho de banda del setup. "
        "La pendiente usa estados aceptados en cruce ±0,25 ps, no predice jitter. Máximo observado no es pico completo. Sin cruce, latencia sin asignar; un checkpoint sin cruce no es censura a 200 ps.")


def four_axes(r,title):
    fig=r.page(title)
    return fig,[fig.add_axes(b) for b in ((.075,.55,.38,.29),(.565,.55,.37,.29),(.075,.18,.38,.27),(.565,.18,.37,.27))]


def voltage_units(values):
    peak=float(np.max(abs(values)))
    return (1e3,"mV") if peak>=1e-3 else ((1e6,"µV") if peak>=1e-6 else (1e9,"nV"))


def circuit(r):
    c=r.c;z=c["z"];t=z["time_s"]*1e12
    fig,ax=four_axes(r,"2.1. Observables eléctricos del circuito")
    ax[0].plot(t,z["Ib_A"]*1e6,label="Ib");ax[0].plot(t,z["Is_A"]*1e6,label="Is");ax[0].set_ylabel("Corriente (µA)");ax[0].legend();ax[0].ticklabel_format(axis="y",useOffset=False)
    scale,unit=voltage_units(z["Vdev_V"]);ax[1].plot(t,z["Vdev_V"]*scale);ax[1].set_ylabel("Voltaje sobre la cinta Vdev ("+unit+")")
    scale,unit=voltage_units(z["vc_V"]);ax[2].plot(t,z["vc_V"]*scale);ax[2].set_ylabel("Voltaje del acoplamiento vc ("+unit+")")
    power=c["m"]["physicalparams"]["circuit"]["R_load_ohm"]*(z["Ib_A"]-z["Is_A"])**2
    scale,unit=(1e9,"nW") if power.max()>=1e-9 else ((1e12,"pW") if power.max()>=1e-12 else (1e15,"fW"))
    ax[3].plot(t,power*scale);ax[3].set_ylabel("Potencia disipada en lectura ("+unit+")")
    for a in ax:a.set_xlabel("Tiempo postcascada (ps)");a.grid(alpha=.25)
    r.save(fig,"Circuito","Magnitudes eléctricas del circuito completo. Vdev es la caída entre los cortes de la cinta y no es Vout. "
        "La potencia de lectura es Rlectura(Ib−Is)². Sólo ventana guardada: no recuperación ni energía de un pulso completo.")


def internal(r):
    c=r.c;s=c["s"];t=s["time_s"]*1e12
    dc=abs(c["dc"]["delta_bar"])*c["m"]["physicalparams"]["energy_scale_J"]
    ratio=abs(s["delta_J"])/dc
    center=np.asarray(c["m"]["source"].get("center_m",np.mean(s["coordinates_m"],axis=0)))
    j=int(np.argmin(np.sum((s["coordinates_m"]-center)**2,axis=1)))
    fig,ax=four_axes(r,"3. Variables internas: evolución del estado")
    ax[0].plot(t,ratio[:,j],label="Nodo cercano al depósito");ax[0].plot(t,ratio.min(axis=1),label="Mínimo nodal")
    ax[0].set_ylabel("|Δ| / |ΔDC|");ax[0].legend(fontsize=8)
    ax[1].plot(t,s["Te_K"][:,j],label="Te equivalente, centro");ax[1].plot(t,s["Tph_K"][:,j],label="Tph equivalente, centro")
    ax[1].axhline(c["m"]["physicalparams"]["Tc_K"],color="black",ls=":",label="Tc");ax[1].set_ylabel("Temperatura equivalente (K)");ax[1].legend(fontsize=8)
    volumes=c["mesh"]["dual_area_m2"]*c["m"]["physicalparams"]["film_thickness_m"]
    for key,label in (("electron_energy_J_m3","Cuasipartículas"),("phonon_energy_J_m3","Fonones")):
        energy=s[key]@volumes;ax[2].plot(t,(energy-energy[0])/EV,label=label)
    ax[2].set_ylabel("Cambio de energía desde t=0 (eV)");ax[2].legend(fontsize=8)
    ax[3].step(t,c["hotbelt"].astype(int),where="post")
    ax[3].set(ylabel="Hotbelt geométrica en snapshot",yticks=[0,1],yticklabels=["No","Sí"],ylim=(-.1,1.1))
    for a in ax:a.set_xlabel("Tiempo postcascada (ps)");a.grid(alpha=.25)
    xy=s["coordinates_m"][j]*1e9
    r.save(fig,"Variables internas",f"Nodo próximo al depósito: ({xy[0]:.5g}, {xy[1]:.5g}) nm. Te/Tph son observadores de energía, no estados térmicos impuestos. "
        "La energía electrónica guardada es de cuasipartículas, sin condensación/gradiente. Hotbelt: nodos con |Δ|≤0,2|ΔDC| conectan ambos bordes; no equivale al gatillo ni prueba continuidad temporal.")


def map_pages(r,indices,key,title,label,caption,kind="scalar",roi=False,flow=None,cmap="viridis",signed=False):
    c=r.c;s=c["s"];mesh=c["mesh"]
    center=np.asarray(c["m"]["source"].get("center_m",np.mean(s["coordinates_m"],axis=0)))
    xy=s["coordinates_m"].copy()*1e9;xy[:,0]-=center[0]*1e9
    tri=Triangulation(xy[:,0],xy[:,1],mesh["triangles"])
    bounds=(max(float(xy[:,0].min()),-60),min(float(xy[:,0].max()),60),float(xy[:,1].min()),float(xy[:,1].max())) if roi else (
        float(xy[:,0].min()),float(xy[:,0].max()),float(xy[:,1].min()),float(xy[:,1].max()))
    if kind=="gap":
        values={int(i):s["delta_J"][i] for i in indices};lo=0.;hi=max(float(abs(v).max()/MEV) for v in values.values())
    elif kind=="phase":
        dc=c["dc"]["delta_bar"];Eu=c["m"]["physicalparams"]["energy_scale_J"]
        values={int(i):np.ma.array(np.angle(s["delta_J"][i]*np.conj(dc)),mask=abs(s["delta_J"][i])<=.2*abs(dc)*Eu) for i in indices}
        lo,hi=-np.pi,np.pi;cmap="twilight"
    elif kind=="current":
        values={int(i):(np.linalg.norm(flow[int(i)][key],axis=1) if key=="total" else flow[int(i)][key][:,0])/1e10 for i in indices}
        hi=max(float(abs(v).max()) for v in values.values());lo=0. if key=="total" else -hi
        cmap="magma" if key=="total" else "coolwarm"
        arrow_reference=max(float(np.linalg.norm(flow[int(i)][key],axis=1).max()) for i in indices)
    elif kind=="phi":
        values={int(i):flow[int(i)]["phi"]*1e3 for i in indices}
        lo=min(float(v.min()) for v in values.values());hi=max(float(v.max()) for v in values.values())
    else:
        factor=1e-3 if key.endswith("J_m3") else 1.
        values={int(i):s[key][i]*factor for i in indices}
        lo=min(float(v.min()) for v in values.values());hi=max(float(v.max()) for v in values.values())
        if signed:hi=max(abs(lo),abs(hi));lo=-hi;cmap="coolwarm"
    if hi<=lo:hi=lo+max(abs(lo)*1e-6,1e-12)
    for start in range(0,len(indices),6):
        chunk=indices[start:start+6];fig=r.page(title+("" if len(indices)<=6 else f" ({start//6+1})"))
        grid=fig.add_gridspec(2,3,left=.06,right=.93,bottom=.20,top=.85,wspace=.30,hspace=.43)
        axes=[];image=None
        for pos,i in enumerate(chunk):
            ax=fig.add_subplot(grid[pos//3,pos%3]);axes.append(ax);v=values[int(i)]
            if kind=="gap":
                image=draw_complex_gap(ax,tri,v,unit_J=MEV,vmin=lo,vmax=hi,bounds=bounds,samples=401)
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
            ax.set(xlim=bounds[:2],ylim=bounds[2:],xlabel="x − xdep (nm)",ylabel="y (nm)",title=f't={s["time_s"][i]*1e12:.6g} ps')
            ax.set_aspect("equal");ax.tick_params(labelsize=8)
        cb=fig.colorbar(image,ax=axes,orientation="horizontal",fraction=.045,pad=.17,shrink=.75);cb.set_label(label)
        r.save(fig,title,caption+" Coordenadas físicas, x centrado en el depósito. Tiempos guardados sin interpolación temporal. "
            "Escala de color común entre tiempos de esta magnitud, no necesariamente entre magnitudes distintas.")


def balance(r):
    c=r.c;z=c["z"];t=z["time_s"]*1e12;fig,ax=four_axes(r,"4. Balance energético independiente")
    ax[0].plot(t,(z["energy_total_J"]-z["energy_total_J"][0])/EV,label="ΔE total")
    ax[0].plot(t,(z["external_energy_integral_J"]-z["external_energy_integral_J"][0])/EV,"--",label="Trabajo neto externo")
    ax[0].set_ylabel("Energía (eV)");ax[0].legend(fontsize=8)
    ref=c["m"]["source"].get("retained_energy_J",0.)
    if ref<=0:ref=c["m"]["numerical"]["kinetic_reference_energy_J"]
    work=np.r_[abs(z["external_energy_integral_J"][0]),abs(np.diff(z["external_energy_integral_J"]))].cumsum()
    defect=abs(z["energy_defect_J"])/(ref+work)*100
    ax[1].plot(t,defect,label="|Defecto|/(Eref + Σ|ΔW|)");ax[1].axhline(1.,color="black",ls=":",label="Criterio global 1 %")
    ax[1].set_ylabel("Defecto global (%)");ax[1].legend(fontsize=8)
    ax[2].plot(t,z["escape_power_W"]*1e9,label="Escape")
    ax[2].plot(t,z["reservoir_power_out_W"]*1e9,label="Cortes, salida firmada")
    ax[2].set_ylabel("Potencia térmica (nW)");ax[2].legend(fontsize=8)
    ax[3].plot(t,z["instant_power_residual_W"]);ax[3].set_ylabel("Residuo instantáneo independiente (W)")
    ax[3].ticklabel_format(axis="y",style="sci",scilimits=(0,0))
    for a in ax:a.set_xlabel("Tiempo postcascada (ps)");a.grid(alpha=.25)
    c["balance"]=dict(maximum_global_energy_defect_percent=float(defect.max()),
        maximum_absolute_instant_power_residual_W=float(abs(z["instant_power_residual_W"]).max()))
    r.save(fig,"Balance energético","E total incluye película BCS+gradiente+fonones y circuito. Trabajo externo neto descuenta resistencias, escape y cortes. "
        "Normalización: energía retenida (DC: referencia cinética declarada) más trabajo externo absoluto acumulado. No sustituye auditoría ni ajuste experimental.")


def provenance(r,checks):
    c=r.c;z=c["z"];fig=r.page("5. Precisión, trazabilidad y alcance")
    ax=fig.add_axes([.075,.62,.38,.21]);keep=z["accepted_step_ps"]>0
    ax.plot(z["time_s"][keep]*1e12,z["accepted_step_ps"][keep])
    ax.set(xlabel="Tiempo (ps)",ylabel="Paso aceptado (ps)",yscale="log");ax.grid(alpha=.25)
    table(fig,[.55,.56,.39,.28],"Comprobaciones de este reporte",[
        ("Auditoría terminal asociada",c["audit"]["status"]),("Pasos aceptados / rechazos",str(len(z["time_s"])-1)+" / "+str(c["execution"].get("total_rejections","no guardado"))),("Tiempo de pared de sesión (h)",fmt(c["execution"].get("elapsed_session_wall_s"),1/3600)),("Corrientes globales positivas",str(c["metrics"]["positive_global_currents"])),
        ("Identidad de lectura, error (V)",fmt(c["metrics"]["readout_identity_error_V"])),
        ("Reconstrucción Vdev, error máx. (V)",fmt(max(x["Vdev_error_V"] for x in checks))),
        ("Defecto energético global máx. (%)",fmt(c["balance"]["maximum_global_energy_defect_percent"])),
        ("Final trayectoria / mapas (ps)",fmt(c["metrics"]["final_ps"])+" / "+fmt(c["metrics"]["maps_final_ps"]))],8.5)
    lines=["Datos: "+str(c["run"]),"Manifiesto SHA256: "+c["inputs"]["manifest.json"],
        ("Checkpoint capturado SHA256: "+c["inputs"]["checkpoint_capture.npz"]) if c["provisional"] else ("Trayectoria SHA256: "+c["inputs"]["trajectory.npz"]),
        "Malla SHA256: "+c["inputs"]["geometry_mesh"],"Pipeline SHA256: "+sha256(__file__),
        "Corrientes: derivada gauge exacta del gradiente BCS y neutralidad del mismo modelo fijado por hashes. "
        "Vectores J por mínimos cuadrados ponderados de flujos duales: aproximación de visualización, no nuevos grados de libertad.",
        "No se guardaron p/n completas en cada snapshot: no se inventan mapas de tasas de colisión, DOS fuera del cierre o calor KWT.",
        "Rigidez radial/KWT fuera del equilibrio, material acústico y fuente postcascada son hipótesis efectivas. "
        "No se infieren jitter, eficiencia estadística, latencia óptica absoluta ni reset completo de esta ventana."]
    y=.50
    for line in lines:
        text=textwrap.fill(line,145);fig.text(.075,y,text,fontsize=9,color=INK,va="top",linespacing=1.35)
        y-=.027*len(text.splitlines())+.014
    r.save(fig,"Trazabilidad","Reporte individual: no compara corridas ni promueve admisión. Un checkpoint provisional no es resultado terminal. "
        "Entradas, hashes, tiempos y errores de reconstrucción se conservan en el recibo JSON acompañante.")


def build_report(run,output,allow_incomplete=False,requested=None,dpi=160,memory_limit=2**33):
    start=time.monotonic();run=Path(run).resolve();output=scratch_path(output)
    rid=case_id(run);target=output/("case_report_"+rid+".pdf");receipt_path=target.with_suffix(".json")
    if target.exists() or receipt_path.exists():raise FileExistsError("Existing report preserved; choose a fresh destination")
    build=output/".build"/rid;build.mkdir(parents=True,exist_ok=False)
    c=load_case(run,build,allow_incomplete,memory_limit);c["metrics"]=metrics(c)
    dc=c["dc"]["delta_bar"]*c["m"]["physicalparams"]["energy_scale_J"]
    c["hotbelt"]=hotbelts(c["s"]["delta_J"],dc,c["mesh"]["edges"],c["mesh"]["coordinates_m"])
    indices=select_snapshots(c["s"]["time_s"],c["obs"],c["hotbelt"],requested)
    print(json.dumps(dict(status="CASE_REPORT_STARTED",case=run.name,provisional=c["provisional"],snapshots=len(indices))),flush=True)
    flow,checks=currents(c,indices)
    apply_thesis_style();plt.rcParams.update({"font.size":10,"axes.labelsize":10,"axes.titlesize":10,"legend.fontsize":9})
    temporary=build/"case_report.pdf";r=Report(temporary,c,dpi)
    try:
        parameters(r);observables(r);circuit(r);internal(r)
        map_pages(r,indices,"delta_J","3.1. Condensado: dominio completo","|Δ| (meV)","Magnitud después de interpolar espacialmente el gap complejo P1; no se interpola el módulo nodal.",kind="gap")
        map_pages(r,indices,"delta_J","3.2. Condensado: región del depósito","|Δ| (meV)","Recorte visual central de 120 nm, ancho completo; no frontera térmica.",kind="gap",roi=True)
        map_pages(r,indices,"delta_J","3.3. Condensado: fase relativa al DC","arg[Δ(t) ΔDC*] (rad)","Fase nodal sin interpolar ángulos. Se ocultan nodos con |Δ|≤0,2|ΔDC|. Saltos ±π de fase envuelta no prueban por sí solos un vórtice.",kind="phase",roi=True)
        map_pages(r,indices,"Te_K","3.4. Temperatura electrónica equivalente","Te equivalente (K)","Observador de energía de ocupaciones electrónicas no térmicas; no temperatura impuesta.",cmap="inferno")
        map_pages(r,indices,"Tph_K","3.5. Temperatura fonónica equivalente","Tph equivalente (K)","Observador de energía fonónica; el cierre es local con escape, sin difusión fonónica lateral.",cmap="inferno")
        map_pages(r,indices,"u_excess_J_m3","3.6. Densidad energética en exceso del DC","Energía en exceso (kJ/m³)","BCS local más gradiente repartido y fonones, menos la referencia DC guardada, dividido por volumen dual.",signed=True)
        for key,title,label in (("super","3.7. Corriente superconductora reconstruida",r"$J_{s,x}$ ($10^{10}$ A/m$^2$)"),
            ("normal","3.8. Corriente normal reconstruida",r"$J_{n,x}$ ($10^{10}$ A/m$^2$)"),("total","3.9. Corriente total reconstruida",r"$|J_{\rm total}|$ ($10^{10}$ A/m$^2$)")):
            map_pages(r,indices,key,title,label,"Vectores por mínimos cuadrados de flujos duales conservativos. Color físico; flechas de dirección con longitud visual fija de 6 nm; se omiten flujos menores a 0,01 % del máximo entre snapshots. "
                "Js/Jn: componente x firmada; total: magnitud. Signo local no es inversión global.",kind="current",flow=flow)
        map_pages(r,indices,"phi","3.10. Potencial eléctrico algebraico","φ (mV)","Neutralidad del mismo estado guardado. Gauge: corte derecho φ=0; sin nueva evolución.",kind="phi",flow=flow)
        balance(r);provenance(r,checks)
    finally:r.pdf.close()
    receipt=dict(schema=SCHEMA,run_dir=str(run),case_name=run.name,report_path=str(target),report_sha256=sha256(temporary),
        generated_utc=datetime.now(timezone.utc).isoformat(),pages=r.pages,pipeline_sha256=sha256(__file__),
        inputs_sha256=c["inputs"],complete_input=not c["provisional"],provisional=c["provisional"],
        model_simulations_launched=0,kinetic_RHS_evaluations=0,new_DC_solves=0,cross_run_comparison=False,
        electrical_metrics=c["metrics"],energy_balance_metrics=c["balance"],existing_audit=c["audit"],
        selected_saved_times_ps=(c["s"]["time_s"][indices]*1e12).tolist(),reconstruction_checks=checks,
        reconstruction_source_sha256={p:sha256(ROOT/p) for p in RECON_SOURCES},
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
    p.add_argument("--snapshots-ps",type=float,nargs="+",default=None,help="Nearest saved times; initial/final included; max24")
    p.add_argument("--workers",type=int,default=1,help="Independent foreground reports; one physical core per worker")
    p.add_argument("--dpi",type=int,default=160);a=p.parse_args(argv)
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
