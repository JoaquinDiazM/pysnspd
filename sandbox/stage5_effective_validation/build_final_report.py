"""Final physical stage-5 report, gated by an explicit scientific admission.

This program NEVER launches a simulation. It reads completed saved trajectories.
All figures, LaTeX, logs and PDF previews stay in scratch. No report is authored
unless the supplied JSON contains stage5: true. Publication additionally needs
a visual-QA receipt matching the compiled PDF SHA256.

Admission contract:
  {"stage5": true, "date": "2026-10-11",
   "primary_runs": {"775": "/home/jdiaz/scratch/...", "1550": "..."},
   "evidence_audits": ["/home/jdiaz/scratch/.../audit.json", ...]}
Entries may instead be {"path": "...", "sha256": "..."}.
Optional report: {"snapshot_times_ps": [...], "equivalent_times_ps": [...]}.
Optional dc_validation: {path:".../current_v6_DC_evidence.json",sha256:"..."};
current13-source basal hold only, never photon/global admission.
Optional primary_run_sha256: {"775": {"manifest.json": "...", ...}, ...}.
Optional physical_sensitivity_runs: {reference_wavelength_nm:775,
 common_interval_ps:[0,20],map_time_ps:15.3,cases:{narrow:{path,sha256:{...},
 audit:{path,sha256},original_STARTED_log:{path,sha256},original_output_dir},wide:{...}},
 bounded_guard_protocol:{path,sha256}}; the protocol is optional only when guards match.
 sha256 pins manifest,trajectory,snapshots,
 checkpoint and initial source. Physical sensitivity is never numerical admission.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import inspect
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import textwrap

REPO=Path(__file__).resolve().parents[2]
SCRATCH=Path("/home/jdiaz/scratch")
DEFAULT_BUILD=SCRATCH/"stage5_effective_final_report_20261011"
DEFAULT_PDF=REPO/"output/pdf/implementation/Informe_final_etapa5_modelo_efectivo_20261011.pdf"
sys.path.insert(0,str(REPO))
sys.pycache_prefix=str(SCRATCH/"pycache")
for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="1"
os.environ.setdefault("MPLCONFIGDIR",str(SCRATCH/"tmp/pysnspd/matplotlib"))
os.environ.setdefault("TMPDIR",str(SCRATCH/"tmp/jdiaz"))
EV_J=1.602176634e-19
MEV_J=EV_J*1e-3
COLORS={"775":"#1764a5","1550":"#bf4b16"}
from pysnspd.plotting.complex_gap import draw_complex_gap
COMPLEX_GAP_SOURCE = Path(inspect.getfile(draw_complex_gap)).resolve()

def sha(path):
    digest=hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda:stream.read(1024**2),b""):
            digest.update(block)
    return digest.hexdigest()

def inside(path,parent):
    path=Path(path).expanduser().resolve()
    if not path.is_relative_to(parent.resolve()):
        raise ValueError("Path must remain in "+str(parent)+": "+str(path))
    return path

def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf8"))

def pinned_entry(entry):
    if isinstance(entry,str):
        path,expected=Path(entry),None
    elif isinstance(entry,dict):
        path,expected=Path(entry["path"]),entry.get("sha256")
    else:
        raise ValueError("Evidence entry must be a path or a pinned path object")
    path=path.expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(path)
    actual=sha(path)
    if expected is not None and actual!=expected:
        raise ValueError("Evidence checksum changed: "+str(path))
    return path,actual

def admission_inputs(path):
    path=Path(path).expanduser().resolve()
    data=read_json(path)
    # An external audit explicitly says it cannot infer this scientific decision.
    # Never promote all_requested_checks_passed to stage5 admission implicitly.
    if data.get("stage5") is not True:
        raise ValueError("Final report blocked: admission JSON must explicitly contain stage5: true")
    primary=data.get("primary_runs",{})
    if set(primary)!={"775","1550"}:
        raise ValueError("Completed 775 and 1550 nm primary_runs are required")
    evidence=data.get("evidence_audits")
    if not isinstance(evidence,list) or not evidence:
        raise ValueError("Pinned evidence_audits are required")
    records=[];comparisons=[];pending_windows=[];input_hashes={str(path):sha(path)}
    for entry in evidence:
        item,digest=pinned_entry(entry)
        audit=read_json(item)
        if audit.get("schema")=="pysnspd.explicit_comparison_window.v1":
            pending_windows.append((item,digest,audit))
            continue
        if audit.get("all_requested_checks_passed") is False:
            raise ValueError("An explicitly failed saved-data audit cannot support closure: "+str(item))
        if audit.get("passed") is False or audit.get("errors"):
            raise ValueError("Failed evidence must be resolved before final reporting: "+str(item))
        for row in audit.get("comparisons",[]):
            if row.get("passed") is not True:
                raise ValueError("Unpassed convergence comparison: "+str(item))
            comparisons.append(row)
        input_hashes[str(item)]=digest
        records.append({"path":str(item),"sha256":digest,"schema":audit.get("schema")})
    from sandbox.stage5_effective_validation.audit_runs import audit_run
    runs={}
    for wavelength,entry in primary.items():
        runpath=Path(entry if isinstance(entry,str) else entry["path"]).expanduser().resolve()
        inside(runpath,SCRATCH)
        pins=data.get("primary_run_sha256",{}).get(wavelength,{})
        if isinstance(entry,dict):
            pins={**pins,**entry.get("sha256",{})}
        for filename,expected in pins.items():
            actual=sha(runpath/filename)
            if actual!=expected:
                raise ValueError("Primary run checksum changed: "+str(runpath/filename))
        run=audit_run(wavelength,runpath)
        if run["summary"]["passed"] is not True:
            raise ValueError("Primary run fails the independent saved-data audit: "+wavelength)
        declared=float(run["manifest"]["source"]["wavelength_nm"])
        if abs(declared-float(wavelength))>1e-10:
            raise ValueError("Wavelength label contradicts the primary run")
        if run["manifest"].get("run_complete") is not True:
            raise ValueError("Incomplete trajectory cannot appear as final evidence")
        for filename in ("manifest.json","trajectory.npz","snapshots.npz","checkpoint_latest.npz"):
            item=runpath/filename
            input_hashes[str(item)]=sha(item)
        runs[wavelength]=run
    if pending_windows:
        from sandbox.stage5_effective_validation.window_evidence import validate_window_evidence
        review_entry=data.get("finite_control_coverage_review")
        if not isinstance(review_entry,dict) or not isinstance(review_entry.get("sha256"),str) or len(review_entry["sha256"])!=64:
            raise ValueError("Finite windows require a separate SHA-pinned manual critical-coverage review")
        review_path,review_sha=pinned_entry(review_entry)
        reviews=read_json(review_path)
        input_hashes[str(review_path)]=review_sha
        def verify_window_input(entry):
            source,source_sha=pinned_entry(entry)
            input_hashes[str(source)]=source_sha
        for item,digest,window in pending_windows:
            wave=str(int(float(window["wavelength_nm"])))
            if wave not in runs:raise ValueError("Window wavelength absent from final primaries")
            row=validate_window_evidence(window,digest,reviews,
                verify_evidence=verify_window_input,primary_summary=runs[wave]["summary"])
            comparisons.append(row)
            input_hashes[str(item)]=digest
            records.append({"path":str(item),"sha256":digest,"schema":window["schema"],
                "scope":"Explicit finite comparison; original completion and events preserved; manual critical coverage mandatory"})
    kinds={row["kind"] for row in comparisons if "kind" in row}
    missing={"temporal","mesh","domain","spectral"}-kinds
    if missing:
        raise ValueError("Final report requires saved convergence evidence for "+", ".join(sorted(missing)))
    for row in comparisons:
        if row.get("robust_spatial") is not None:
            from sandbox.stage5_effective_validation import guard_protocol_dispatcher as dispatch
            correction=dispatch.physical_gap_module().validate_physical_gap_correction(row)
            for entry in correction["analysis_sources"].values():
                item,digest=pinned_entry(entry);input_hashes[str(item)]=digest
            input_hashes[str(Path(dispatch.__file__).resolve())]=sha(dispatch.__file__)
        if "spectral_observable_policy" in row:
            from sandbox.stage5_effective_validation import guard_protocol_dispatcher as dispatch
            dispatch.validate_spectral_policy_application(row)
            source=row.get("original_full_run_comparison",row)
            entries=[source["spectral_observable_policy"]["decision"],
                *source["verified_spectral_application_scope"]["original_manifests"].values()]
            if "exploratory_spectral_spatial_policy" in source:
                entries.append(source["exploratory_spectral_spatial_policy"]["decision"])
            if "exploratory_spectral_voltage_policy" in source:
                entries.append(source["exploratory_spectral_voltage_policy"]["decision"])
            for entry in entries:
                item,digest=pinned_entry(entry);input_hashes[str(item)]=digest
            modules=[dispatch,dispatch.spectral_policy_module()]
            if "exploratory_spectral_spatial_policy" in source:modules.append(dispatch.exploratory_policy_module())
            for module in modules:
                input_hashes[str(Path(module.__file__).resolve())]=sha(module.__file__)
    return data,runs,records,comparisons,input_hashes

def load_dc_validation(entry,reference,*,input_hashes,load_run=None):
    """Optional current-code basal evidence; no photon/global admission inferred."""
    from sandbox.stage5_effective_validation import audit_runs as audit
    import numpy as np
    if not isinstance(entry,dict) or len(entry.get("sha256",""))!=64:
        raise ValueError("dc_validation requires a SHA-pinned evidence receipt")
    receipt_path,receipt_sha=pinned_entry(entry);receipt=read_json(receipt_path)
    if receipt.get("schema")!="pysnspd.current_effective_DC_evidence.v1" or receipt.get("status")!="DC_HOLD_COMPLETED_AND_INDEPENDENTLY_AUDITED":
        raise ValueError("Current independently audited DC hold required")
    if receipt.get("no_photon") is not True or receipt.get("stage5_admission_inferred") is not False:
        raise ValueError("DC evidence cannot imply photon or global admission")
    root=inside(receipt["path"],SCRATCH)
    files=("manifest.json","trajectory.npz","snapshots.npz","checkpoint_latest.npz","dc/effective_dc.npz")
    for filename in files:
        path=root/filename;expected=receipt.get("input_sha256",{}).get(str(path))
        if not isinstance(expected,str) or len(expected)!=64:raise ValueError("DC evidence lacks required data SHA: "+filename)
        source,digest=pinned_entry({"path":str(path),"sha256":expected});input_hashes[str(source)]=digest
    ae=receipt.get("audit")
    if not isinstance(ae,dict) or len(ae.get("sha256",""))!=64:raise ValueError("DC audit must be SHA-pinned")
    ap,ad=pinned_entry(ae);document=read_json(ap)
    if document.get("schema")!="pysnspd.effective_saved_data_audit.v2" or document.get("all_requested_checks_passed") is not True or document.get("errors"):
        raise ValueError("DC evidence must pass the adopted saved-data audit")
    if document.get("auditor_source_sha256")!=sha(audit.__file__):raise ValueError("DC audit source differs from adopted auditor")
    found=[row for row in document.get("runs",{}).values() if Path(row["path"]).resolve()==root]
    if len(found)!=1 or found[0].get("passed") is not True:raise ValueError("DC audit must identify this exact completed run")
    refenergy=float(receipt["global_reference_energy_J"])
    if not math.isfinite(refenergy) or refenergy<=0:raise ValueError("Positive DC reference energy required")
    load_run=audit.audit_run if load_run is None else load_run
    run=load_run("dc_validation",root,dc_reference_energy_J=refenergy);sm=run["summary"];manifest=run["manifest"]
    if sm.get("passed") is not True or sm.get("declared_mode")!="DC" or sm.get("status")!="DC_HOLD_COMPLETED" or manifest.get("run_complete") is not True:
        raise ValueError("DC evidence must remain actually complete and photon-free")
    if found[0]["checks"]!=sm["checks"] or found[0]["final_time_ps"]!=sm["final_time_ps"]:
        raise ValueError("DC saved audit differs from independent recomputation")
    if not math.isclose(float(sm["final_time_ps"]),10.,rel_tol=0.,abs_tol=1e-10):raise ValueError("Declared DC validation requires10ps")
    pins=manifest.get("source_pins_sha256",{})
    if len(pins)!=13 or pins!=reference["manifest"].get("source_pins_sha256") or pins!=receipt.get("source_pins_sha256"):
        raise ValueError("DC validation requires the identical full13 current-source pins")
    if manifest["physicalparams"]!=reference["manifest"]["physicalparams"]:
        raise ValueError("DC material, operating point and full circuit must match the nominal primary")
    history=run["history"];snap=run["snap"]
    measured=dict(max_Vout_V=float(np.max(abs(history["Vout_V"]))),max_Vdev_V=float(np.max(abs(history["Vdev_V"]))),
        max_energy_defect_J=sm["energy"]["maximum_absolute_defect_J"],
        max_gap_drift_J=float(np.max(abs(abs(snap["delta_J"])-abs(snap["delta_J"][0])))),
        max_Te_drift_K=float(np.max(abs(snap["Te_K"]-snap["Te_K"][0]))),
        max_Tph_drift_K=float(np.max(abs(snap["Tph_K"]-snap["Tph_K"][0]))),
        max_Is_drift_A=float(np.max(abs(history["Is_A"]-history["Is_A"][0]))))
    if not all(math.isfinite(value) and value>=0 for value in measured.values()):raise ValueError("Nonfinite DC display quantity")
    if any(float(receipt.get(key,float("nan")))!=value for key,value in measured.items()):
        raise ValueError("DC receipt metrics differ from saved arrays")
    if int(receipt["accepted_steps"])!=len(history["time_s"])-1:raise ValueError("DC accepted-step count differs")
    for filename in files:pinned_entry({"path":str(root/filename),"sha256":receipt["input_sha256"][str(root/filename)]})
    input_hashes[str(receipt_path)]=receipt_sha;input_hashes[str(ap)]=ad
    ldc=float(manifest["physicalparams"]["exterior_inductance"]["total_H"])
    lref=float(reference["manifest"]["physicalparams"]["exterior_inductance"]["total_H"])
    display=dict(**measured,duration_ps=sm["final_time_ps"],accepted_steps=int(receipt["accepted_steps"]),
        same13pins=True,source_mode="DC",own_frozen_exterior_inductance_H=ldc,
        primary_exterior_inductance_H=lref,exterior_inductance_difference_fraction=(ldc-lref)/lref,
        stage5_admission_inferred=False,scope="Own current-code basal hold only; historical holds are distinct antecedents")
    record=dict(path=str(receipt_path),sha256=receipt_sha,schema=receipt["schema"],scope=display["scope"])
    return display,record

def dc_validation_latex(dc):
    if dc is None:return ""
    return (r"\textbf{Basal propio de esta versión.} "+number(dc["duration_ps"])+
        r" ps sin fotón ("+str(dc["accepted_steps"])+r" pasos): $\max|V_{\rm out}|="+
        number(dc["max_Vout_V"]*1e9)+r"$ nV, $\max|V_{\rm dev}|="+number(dc["max_Vdev_V"]*1e9)+
        r"$ nV; defecto energético máximo "+number(dc["max_energy_defect_J"])+
        r" J. Derivas de $|\Delta|$, $T_e$, $T_{\rm ph}$ e $I_s$: "+
        ("nulas en los datos guardados." if all(dc[key]==0 for key in ("max_gap_drift_J","max_Te_drift_K","max_Tph_drift_K","max_Is_drift_A")) else
         r"registradas en el recibo DC; no se sustituyen por cero.")+
        r" Este hold verifica el basal; su alcance se separa de la validación fotónica.")

def select_saved(time_s,count,requested_ps=None):
    times=np.asarray(time_s,float)*1e12
    if requested_ps is None:
        if len(times)<=count:
            return np.arange(len(times))
        positive=times[times>times[0]]
        first=max(float(positive[0]),min(.1,float(times[-1])/100))
        request=np.r_[times[0],np.geomspace(first,float(times[-1]),count-1)]
    else:
        request=np.asarray(requested_ps,float)
        if len(request)>count:
            raise ValueError("At most "+str(count)+" requested saved times fit this report page")
        if np.any(~np.isfinite(request)) or np.any(request<times[0]) or np.any(request>times[-1]):
            raise ValueError("Report snapshot requests must lie inside the saved interval")
        request=np.r_[times[0],request,times[-1]]
    indices=np.unique([int(np.argmin(abs(times-value))) for value in request])
    if len(indices)>count:
        indices=indices[np.unique(np.rint(np.linspace(0,len(indices)-1,count)).astype(int))]
    return indices

def save_figure(fig,build,name,records,*,units,reference,scope):
    target=build/"figures"/(name+".png")
    fig.savefig(target,dpi=190,facecolor="white")
    plt.close(fig)
    records.append({"path":str(target),"sha256":sha(target),"units":units,
                    "reference":reference,"scope":scope})
    return "figures/"+target.name


def add_voltage_overview_inset(ax,runs,threshold_mV):
    """Display-only inset for small signals; same saved times and DC subtraction."""
    threshold=float(threshold_mV)
    if not math.isfinite(threshold) or threshold<=0:
        raise ValueError("Positive saved electrical threshold required")
    peaks=[]
    for run in runs.values():
        baseline=float(run["manifest"].get("observation",{}).get("baseline_Vout_V",0.))
        peaks.append(float(abs(run["history"]["Vout_V"]-baseline).max())*1e3)
    peak=max(peaks,default=0.)
    if not 0<peak<.2*threshold:return False
    inset=ax.inset_axes([.15,.28,.81,.43])
    for wave,run in runs.items():
        baseline=float(run["manifest"].get("observation",{}).get("baseline_Vout_V",0.))
        inset.plot(run["history"]["time_s"]*1e12,
            (run["history"]["Vout_V"]-baseline)*1e6,color=COLORS[wave])
    inset.set(xlabel="Tiempo (ps)",ylabel=r"$V_{out}-V_{basal}$ (µV)",
        title="Misma señal, escala ampliada")
    inset.tick_params(labelsize=9)
    inset.xaxis.label.set_size(9);inset.yaxis.label.set_size(9)
    inset.title.set_size(10);inset.grid(alpha=.18)
    return True


def overview(runs,build,records):
    fig,axs=plt.subplots(2,2,figsize=(11,8),layout="constrained")
    for wave,run in runs.items():
        color=COLORS[wave];h=run["history"];s=run["snap"]
        t=h["time_s"]*1e12;ts=s["time_s"]*1e12
        baseline=float(run["manifest"].get("observation",{}).get("baseline_Vout_V",0))
        axs[0,0].plot(t,(h["Vout_V"]-baseline)*1e3,color=color,label=wave+" nm")
        axs[0,1].plot(t,h["Is_A"]*1e6,color=color,label=r"$I_s$, "+wave+" nm")
        axs[0,1].plot(t,h["Ib_A"]*1e6,color=color,ls="--",label=r"$I_b$, "+wave+" nm")
        axs[1,0].plot(ts,abs(s["delta_J"]).min(axis=1)/MEV_J,color=color,label=wave+" nm")
        axs[1,1].plot(ts,s["Te_K"].max(axis=1),color=color,label=r"$T_e^{eq}$, "+wave+" nm")
        axs[1,1].plot(ts,s["Tph_K"].max(axis=1),color=color,ls="--",label=r"$T_{ph}^{eq}$, "+wave+" nm")
    threshold=float(runs["775"]["manifest"].get("parameter_identity",{}).get("threshold_mV",.1))
    axs[0,0].axhline(threshold,color=".35",ls=":",label="Gatillo declarado")
    for ax,label,title in zip(axs.ravel(),
            [r"$V_{out}-V_{basal}$ (mV)","Corriente (µA)",r"$\min_{\rm nodos}|\Delta|$ (meV)","Temperatura equivalente (K)"],
            ["A. Señal de lectura intrínseca","B. Corrientes globales",
             "C. Pérdida local de condensado","D. Energía en electrones y fonones"]):
        ax.set(xlabel="Tiempo desde transferencia (ps)",ylabel=label,title=title)
        ax.grid(alpha=.18);ax.legend(fontsize=13)
    inset_added=add_voltage_overview_inset(axs[0,0],runs,threshold)
    return save_figure(fig,build,"01_electrical_overview",records,
        units="ps,mV,uA,meV,K"+(",uV" if inset_added else ""),reference="Vout minus the saved fresh-DC baseline; inset shows identical saved Vout in uV only if its peak is below20%of threshold; absolute gap and equivalent temperatures",
        scope="Completed deterministic post-cascade windows; electrical trigger comes only from accepted-step controller")

def gap_maps(wave,run,build,records,requested=None):
    from matplotlib.tri import Triangulation
    s=run["snap"];xy=s["coordinates_m"]*1e9;tri=Triangulation(xy[:,0],xy[:,1],s["triangles"])
    values=abs(s["delta_J"])/MEV_J
    indices=select_saved(s["time_s"],6,requested)
    rows=len(indices)
    fig,axs=plt.subplots(rows,2,figsize=(11,2.05*rows+1.2),
                         squeeze=False,layout="constrained")
    center=.5*(xy[:,0].min()+xy[:,0].max())
    vmax=float(values.max())
    for row,index in enumerate(indices):
        for column,ax in enumerate(axs[row]):
            bounds=None if column==0 else (center-60,center+60,xy[:,1].min(),xy[:,1].max())
            image=draw_complex_gap(ax,tri,s["delta_J"][index],unit_J=MEV_J,
                vmin=0,vmax=vmax,bounds=bounds)
            ax.set(aspect="equal",xlabel="x (nm)",ylabel="y (nm)")
            if column==1:
                ax.set_xlim(center-60,center+60)
            suffix="dominio completo" if column==0 else "ampliación central de 120 nm"
            ax.set_title(f'{s["time_s"][index]*1e12:.5g} ps | '+suffix,fontsize=14)
    fig.colorbar(image,ax=axs.ravel().tolist(),label=r"$|\Delta|$ (meV)",shrink=.75)
    fig.suptitle(wave+" nm: estados realmente guardados",fontsize=15)
    path=save_figure(fig,build,"gap_"+wave,records,units="nm,ps,meV",
        reference="Absolute |complex spatial P1 Delta| on nearest-cell raster (at most601 samples per axis); identical color scale over all saved snapshots; no time interpolation",
        scope="Full Delaunay geometry and central visual zoom, without restricting evolved thermal domain")
    return path,(s["time_s"][indices]*1e12).tolist()

def energy_maps(wave,run,build,records,requested=None):
    s=run["snap"];xy=s["coordinates_m"]*1e9;tri=s["triangles"]
    if "u_excess_J_m3" not in s:
        raise ValueError("Physical excess-energy density snapshots are required")
    if requested is None:
        times=s["time_s"]*1e12
        middle=run["summary"]["hotbelt"]["first_saved_width_spanning"]
        target=times[-1]/2 if middle is None else middle["first_saved_time_s"]*1e12
        requested=[float(times[0]),float(target),float(times[-1])]
    indices=select_saved(s["time_s"],3,requested)
    fig,axs=plt.subplots(3,len(indices),figsize=(11,8.7),squeeze=False,layout="constrained")
    center=.5*(xy[:,0].min()+xy[:,0].max())
    fields=(("Te_K",r"$T_e^{eq}$ (K)","inferno"),
            ("Tph_K",r"$T_{ph}^{eq}$ (K)","inferno"),
            ("u_excess_J_m3",r"$\delta u$ (J/m$^3$)","magma"))
    for row,(key,label,cmap) in enumerate(fields):
        values=s[key]
        lo,hi=float(values.min()),float(values.max())
        if key=="u_excess_J_m3" and lo<0:
            hi=max(abs(lo),abs(hi));lo=-hi;cmap="coolwarm"
        if hi<=lo:hi=lo+max(1e-25,abs(lo)*1e-12)
        for column,index in enumerate(indices):
            ax=axs[row,column]
            image=ax.tripcolor(xy[:,0],xy[:,1],tri,
                facecolors=values[index][tri].mean(axis=1),shading="flat",
                cmap=cmap,vmin=lo,vmax=hi,rasterized=True)
            ax.set(aspect="equal",xlim=(center-60,center+60),xlabel="x (nm)",ylabel="y (nm)")
            ax.set_title(f'{s["time_s"][index]*1e12:.5g} ps',fontsize=14)
        fig.colorbar(image,ax=axs[row].tolist(),label=label,shrink=.78)
    fig.suptitle(wave+" nm: energía del volumen alrededor de la absorción",fontsize=14)
    return save_figure(fig,build,"energy_maps_"+wave,records,units="nm,ps,K,J/m3",
        reference="Saved equivalent temperatures and saved energy excess above same fresh DC; absolute common color limits across all snapshots",
        scope="120nm visualization only; no truncated thermal evolution and no assumption of persistent FD/BE distributions")

def energy_balance(runs,build,records):
    fig,axs=plt.subplots(2,2,figsize=(11,7.6),layout="constrained")
    for wave,run in runs.items():
        h=run["history"];t=h["time_s"]*1e12;color=COLORS[wave]
        axs[0,0].plot(t,(h["domain_energy_J"]-h["domain_energy_J"][0])/EV_J,
                      color=color,label=wave+" nm, película")
        axs[0,1].plot(t,(h["energy_total_J"]-h["energy_total_J"][0])/EV_J,
                      color=color,label=wave+" nm, energía")
        axs[0,1].plot(t,h["external_energy_integral_J"]/EV_J,
                      color=color,ls="--",label=wave+" nm, trabajo")
        axs[1,0].plot(t,run["traces"]["energy_relative_defect"]*100,color=color,label=wave+" nm")
        axs[1,1].plot(t,h["escape_power_W"]*1e9,color=color,label=wave+" nm, sustrato")
        axs[1,1].plot(t,h["reservoir_power_out_W"]*1e9,color=color,ls="--",label=wave+" nm, cortes")
    axs[1,0].axhline(1,color=".35",ls=":",label="Límite 1%")
    labels=["Incremento de energía de película (eV)","Energía o trabajo neto (eV)",
            "|Defecto| / (Eret + Σ|ΔWext|) (%)","Potencia saliente firmada (nW)"]
    titles=["A. Redistribución tras el fotón","B. Película y circuito CM juntos",
            "C. Balance integrado","D. Salidas físicas de energía"]
    for ax,label,title in zip(axs.ravel(),labels,titles):
        ax.set(xlabel="Tiempo desde transferencia (ps)",ylabel=label,title=title)
        ax.grid(alpha=.18);ax.legend(fontsize=13)
    return save_figure(fig,build,"energy_balance",records,units="ps,eV,percent,nW",
        reference="Energy increments from first post-deposition saved state; defect normalized by Eret + integral(abs(net external power))",
        scope="Independent conserved-energy ledger; reservoir and substrate powers retain their signs")


def comparison_admission_metrics(row):
    """Read actual adopted margins per observable, retaining the original2% gate."""
    limit=float(row["relative_limit"])
    if not math.isfinite(limit) or limit<=0:raise ValueError("Positive finite comparison margin required")
    def ratio(value,margin=limit):
        value=float(value)
        if not math.isfinite(value) or value<0 or not math.isfinite(margin) or margin<=0:
            raise ValueError("Finite nonnegative comparison error and positive margin required")
        return value/margin
    spectral="spectral_observable_policy" in row
    nodal_limit=limit;field_limit=limit;deficit_limit=limit;voltage_limit=limit
    exploratory="exploratory_spectral_spatial_policy" in row
    if spectral:
        from sandbox.stage5_effective_validation.guard_protocol_dispatcher import validate_spectral_policy_application
        limits=validate_spectral_policy_application(row)
        nodal_limit=limits["nodal_minimum_over_DC"]
        field_limit=limits["central_field_L2_over_DC"]
        deficit_limit=limits["central_deficit_over_actual_reference_peak"]
        voltage_limit=limits.get("voltage_waveform_relative_limit",limit)
        if nodal_limit is None:nodal_limit=.03 # Historical3% display only; no active nodal gate.
    nodal=ratio(row["gap_min"]["relative_error"],nodal_limit)
    result=dict(Vout_margin_ratio=ratio(row["Vout"]["relative_error"],voltage_limit),
        voltage_gate_relative_limit=voltage_limit,latency_gate_relative_limit=limit,
        nodal_margin_ratio=nodal,nodal_diagnostic=row["kind"] in {"mesh","domain"} or exploratory,
        nodal_former_gate_passed=row["gap_min"].get("former_gate_passed",nodal<=1),
        spectral_multiple_observables=spectral,nodal_current_gate_relative_limit=None if exploratory else nodal_limit,
        nodal_display_is_historical=exploratory,field_gate_relative_limit=field_limit,deficit_gate_relative_limit=deficit_limit)
    if result["nodal_diagnostic"] or spectral:
        if row["kind"] in {"mesh","domain"} and (row.get("policy_version")!="spatial_observables_v2" or not row["gap_min"].get("diagnostic_only")):
            raise ValueError("Mesh/domain final report requires the adopted spatial audit policy")
        spatial=row["robust_spatial"]
        result.update(primary_gap_margin_ratio=ratio(spatial["L2_field"]["maximum_relative_error"],field_limit),
            deficit_margin_ratio=ratio(spatial["deficit"]["relative_error_actual_peak"],deficit_limit),
            primary_gap_name="L2 de magnitud en ROI / gap DC",
            field_normalization="Reference area-weighted initial DC gap in120x80nm ROI",
            deficit_normalization="Actual maximum reference nodal-deficit curve in common interval")
    else:
        result.update(primary_gap_margin_ratio=nodal,deficit_margin_ratio=None,
            primary_gap_name="Mínimo nodal / gap DC",
            field_normalization="Initial DC gap",deficit_normalization=None)
    return result

def convergence_caption_latex(comparisons):
    spectral=any("spectral_observable_policy" in row for row in comparisons)
    exploratory=any("exploratory_spectral_spatial_policy" in row for row in comparisons)
    same=(r"Paso: 2\% para señal y mínimo nodal. Espectro: señal y latencia 2\%, "
        r"mínimo nodal 3\% DC y, obligatoriamente, $L^2$/DC y déficit/pico 2\%. "
        r"Se conserva el fallo original del mínimo espectral con 2\%. " if spectral else
        r"Paso y espectro: 2\% para señal y mínimo nodal del gap. ")
    if exploratory:
        same=(r"Paso: señal y mínimo nodal 2\%. Espectro exploratorio: señal y latencia 2\%, "
            r"$L^2$/DC y déficit/pico real 5\%; mínimo nodal diagnóstico. "
            r"Se conservan íntegros los resultados anteriores de 2\% y 3/2/2. ")
    if any("exploratory_spectral_voltage_policy" in row for row in comparisons):
        same=(r"Paso: señal y mínimo nodal 2\%. Espectro exploratorio: señal 5\%, latencia 2\%, "
            r"$L^2$/DC y déficit/pico real 5\%; mínimo nodal diagnóstico. "
            r"La revisión de señal fue posterior a los datos y conserva el fallo previo de 2\%. "
            r"La resolución 54 es la opción exploratoria por defecto, validada en cuatro ejes a 775 nm "
            r"dentro de estos márgenes; 80 es control de sensibilidad, no verdad continua. ")
    return (r"\textbf{Figura 5.} Unidad = margen permitido. "+same+
        r"Malla y longitud: 5\% para señal, norma $L^2$ de $|\Delta|$ en la ROI física "
        r"de $120\times80$ nm dividida por el gap DC y error de déficit central "
        r"dividido por su pico real de referencia. La señal usa "
        r"$\max(\text{pico real},\text{umbral})$; se conserva también el error relativo "
        r"al pulso real en el recibo. La latencia procede sólo del cruce confirmado. "
        r"Sólo en el control temporal se permite un piso de error absoluto de 0.1 ps; "
        r"no se suma a la latencia. Para amplitudes: $P_1$ espacial complejo, módulo por snapshot "
        r"e interpolación temporal escalar; la métrica temporal compleja anterior queda preservada. "
        r"Las ventanas indicadas limitan la evidencia.")


def spectral_policy_history_latex(comparisons):
    """Current-window original policies first; immutable earlier evidence labelled."""
    rows=[r for r in comparisons if "exploratory_spectral_spatial_policy" in r]
    if not rows:return ""
    from sandbox.stage5_effective_validation import guard_protocol_dispatcher as dispatch
    table=[r"\begin{center}{\small\begin{tabular}{lrrrll}\hline ",
        r"Ventana común [ps] & Mín./DC [\%] & $L^2$/DC [\%] & Déficit/pico [\%] & 2\% & 3/2/2\\\hline "]
    completion={"reference_saved_window_passed","candidate_saved_window_passed"}
    for row in rows:
        raw=row.get("original_full_run_comparison",row)
        old3=raw["original_3_2_2_spectral_result"]
        old2=old3["original_2percent_spectral_result"]
        interval=old3["common_gap_snapshot_interval_ps"]
        values=[old3["gap_min"]["relative_error"],
            old3["robust_spatial"]["L2_field"]["maximum_relative_error"],
            old3["robust_spatial"]["deficit"]["relative_error_actual_peak"]]
        verdict=lambda r:"cumple" if all(v for k,v in r["checks"].items() if k not in completion) else "falla"
        table.append(number(interval[0])+"--"+number(interval[1])+" & "
            +" & ".join(number(v*100) for v in values)+" & "+verdict(old2)+" & "+verdict(old3)+r"\\")
    table.append(r"\hline\end{tabular}}\end{center}")
    entry=rows[0]["exploratory_spectral_spatial_policy"]["decision"]
    _,decision=dispatch.exploratory_policy_module().load_decision(entry["path"])
    old=decision["current_provisional_metrics"];interval=old["interval_gap_ps"]
    prior=[old["nodal_minimum_over_DC"]["value"],old["field_L2_over_DC"]["value"],
        old["deficit_over_partial_actual_peak"]["value"]]
    previous=decision["historical_failures_preserved"]["original_interpolation_metric"]["temporal_original_L2"]
    return ("\n".join(table)+r"Las políticas anteriores se recomputan arriba sobre la misma ventana final "
        r"con amplitud física: $P_1$ espacial complejo, módulo e interpolación temporal escalar. "
        r"El antecedente provisional "+number(interval[0])+"--"+number(interval[1])
        +r" ps conserva mínimo, $L^2$ y déficit "+"/".join(number(v*100) for v in prior)
        +r"\%, con FAIL 2\% y 3/2/2. No se recortan máximos ni se añade una cola hasta 200 ps para diluir el denominador. "
        r"El presupuesto exploratorio no acredita calibración ni convergencia asintótica. "
        +r"La antigua cuerda compleja temporal produjo "+number(previous*100)
        +r"\% DC: se conserva como artefacto del análisis anterior, no como fallo físico de amplitud.")

def compact_guard_segments(segments):
    if not isinstance(segments,list) or not segments:raise ValueError("Explicit numerical timeline required")
    merged=[];previous_end=0.
    for index,segment in enumerate(segments):
        start=segment.get("start_time_ps");end=segment.get("end_time_ps")
        if not isinstance(start,(int,float)) or isinstance(start,bool) or not math.isfinite(start) or start<0:
            raise ValueError("Finite nonnegative numerical segment clock required")
        if start!=previous_end:raise ValueError("Numerical timeline must be continuous from zero")
        if end is not None and (not isinstance(end,(int,float)) or isinstance(end,bool) or not math.isfinite(end) or end<start):
            raise ValueError("Invalid numerical segment endpoint")
        if end is None and index!=len(segments)-1:raise ValueError("Only final timeline segment may be open")
        tolerance=segment.get("local_energy_tolerance")
        if tolerance is not None and (not isinstance(tolerance,(int,float)) or isinstance(tolerance,bool) or not math.isfinite(tolerance) or tolerance<=0):
            raise ValueError("Positive finite guard tolerance required")
        row=dict(start_time_ps=start,end_time_ps=end,local_energy_tolerance=tolerance)
        if merged and merged[-1]["local_energy_tolerance"]==tolerance:merged[-1]["end_time_ps"]=end
        else:merged.append(row)
        previous_end=end
    return merged

def temporal_refinement_display_rows(comparisons):
    rows=[]
    for row in comparisons:
        if row.get("kind")!="temporal":continue
        control=row.get("declared_numerical_control")
        reference_segments=candidate_segments=None;transition_band=None
        if control is not None and control.get("schema")=="pysnspd.bounded_guard_comparison_protocol.v1":
            reference_segments=compact_guard_segments(control["reference"]["segments"])
            candidate_segments=compact_guard_segments(control["candidate"]["segments"])
            reference_guard=reference_segments[-1]["local_energy_tolerance"]
            candidate_guard=candidate_segments[-1]["local_energy_tolerance"]
            transition_band=control["guard_law"]["transition_band_ps"]
        else:
            reference_guard=None if control is None else control["reference_local_energy_tolerance"]
            candidate_guard=None if control is None else control["candidate_local_energy_tolerance"]
        stats=row.get("accepted_step_statistics_after_both_guard_anchors",{})
        rows.append(dict(reference=row["reference"],candidate=row["candidate"],
            reference_guard=reference_guard,candidate_guard=candidate_guard,
            reference_guard_segments=reference_segments,candidate_guard_segments=candidate_segments,
            transition_band_ps=transition_band,
            post_anchor_interval_ps=stats.get("interval_ps") if stats.get("available") else None,
            reference_accepted_steps=stats.get("reference",{}),
            candidate_accepted_steps=stats.get("candidate",{}),
            scope="Actual accepted-step statistics after both anchors; no global order or automatic2x refinement asserted"))
    return rows

def guard_tolerance_latex(value):
    text=format(float(value),".6g")
    if "e" not in text:return text
    mantissa,exponent=text.split("e")
    return mantissa+r"\times10^{"+str(int(exponent))+"}"

def guard_timeline_latex(comparisons):
    rows=temporal_refinement_display_rows(comparisons)
    rows=[r for r in rows if r["reference_guard_segments"] is not None]
    if not rows:return ""
    item=rows[0];lines=[r"\textbf{Política histórica, 775 nm.} "]
    for role,label in (("reference","nominal"),("candidate","control temporal")):
        parts=[]
        for segment in item[role+"_guard_segments"]:
            tolerance=segment["local_energy_tolerance"]
            value="sin guarda" if tolerance is None else r"$\varepsilon_E="+guard_tolerance_latex(tolerance)+"$"
            left=number(segment["start_time_ps"]);end=segment["end_time_ps"]
            interval="desde "+left+" ps" if end is None else left+"--"+number(end)+" ps"
            parts.append(value+" ("+interval+")")
        lines.append(label+": "+", ".join(parts)+". ")
    band=item["transition_band_ps"]
    lines.append("Las transiciones de los controles pertenecen a la banda "+
        number(band[0])+"--"+number(band[1])+" ps; se declara la transición numérica acotada, sin aislamiento exacto del eje dentro de ella.")
    return "".join(lines)


def finite_window_display_rows(comparisons):
    """Display validated finite-window metadata without altering original states."""
    rows=[]
    for row in comparisons:
        if row.get("comparison_scope")!="EXPLICIT_FINITE_COMPARISON_WINDOW":continue
        interval=row.get("declared_comparison_interval_ps")
        if not isinstance(interval,list) or len(interval)!=2 or not all(
            isinstance(v,(float,int)) and not isinstance(v,bool) and math.isfinite(v) for v in interval
        ) or interval[0]<0 or interval[1]<=interval[0]:
            raise ValueError("Finite comparison display requires its declared interval")
        review=row.get("critical_coverage_review",{})
        decision={"COVERED":"COVERED","QUANTIFIED_LATER_BOUND":"BOUND"}.get(review.get("decision"))
        if decision is None:raise ValueError("Finite comparison display requires admitted critical coverage")
        certified=row.get("latency_certified_in_this_window")
        if not isinstance(certified,bool):raise ValueError("Explicit boolean latency certification required")
        originals=row.get("original_full_run_states",{})
        if set(originals)!={"reference","candidate"} or any(
            not {"run_complete","terminal_status"}<=set(originals[side]) for side in originals
        ):raise ValueError("Original run completion/status must remain visible")
        rows.append(dict(axis=row["kind"],interval_ps=list(interval),decision=decision,
            latency_certified=certified,original_run_states=originals))
    return rows

def finite_window_scope_latex(comparisons):
    rows=finite_window_display_rows(comparisons)
    if not rows:return ""
    lines=[r"\begin{center}{\small\begin{tabular}{lrll}\hline ",
        r"Eje & Intervalo [ps] & Cobertura & Latencia cert.\\\hline "]
    for item in rows:
        interval=item["interval_ps"]
        lines.append(tex_escape(item["axis"])+" & "+number(interval[0])+"--"+number(interval[1])+
            " & "+item["decision"]+" & "+("sí" if item["latency_certified"] else "no")+r"\\")
    lines.append(r"\hline\end{tabular}}\end{center}")
    lines.append(r"\texttt{EXPLICIT\_FINITE\_COMPARISON\_WINDOW}: revisión crítica SHA "
        r"obligatoria, COVERED = cubierta, BOUND = cota cuantificada posterior. "
        r"Los estados originales \texttt{run\_complete} y terminales se conservan sin alteración.")
    critical=[row.get("critical_coverage_review",{}).get("critical_interval_end_ps") for row in comparisons]
    critical=[value for value in critical if isinstance(value,(int,float)) and not isinstance(value,bool) and math.isfinite(value)]
    if critical:
        lines.append(r"Formación y latencia: extremo crítico revisado hasta "+number(max(critical))+r" ps.")
    lines.append(r"\textbf{Todas las colas posteriores al intervalo declarado quedan fuera de la comparación de cuatro ejes.} "
        r"Cada original conserva su cruce+5 ps, con confirmación a cruce+0.5 ps, y sus auditorías individuales "
        r"de energía, soporte, eventos y signo de corriente. No se certifica el pico completo ni el reset.")
    return "\n".join(lines)

def expected_report_pages(config,has_sensitivity=False):
    return 11+2*bool(config.get("include_vortex_diagnostics",False))+int(bool(has_sensitivity))

def physical_width_maps(diagnostic,runs,build,records):
    """Display actual physical-width runs; no error gate or event interpolation."""
    import matplotlib.tri as mtri
    interval=diagnostic["common_interval_ps"]
    labels={"nominal":"σ nominal","narrow":"σ/√2","wide":"√2 σ"}
    colors={"nominal":"#1764a5","narrow":"#13868b","wide":"#bf7a16"}
    fig=plt.figure(figsize=(11,6.8),layout="constrained")
    grid=fig.add_gridspec(2,6,height_ratios=[1.,1.05])
    axv=fig.add_subplot(grid[0,:3]);axd=fig.add_subplot(grid[0,3:])
    nominal=runs["nominal"]["spatial"]
    nt=nominal["time_s"]*1e12;inside=(nt>=interval[0])&(nt<=interval[1])
    peak=float(max(nominal["deficit_J_m2"][inside]))
    if peak<=0:raise ValueError("Physical sensitivity requires a nonzero nominal central deficit to normalize")
    vmax=max(float(abs(run["snap"]["delta_J"][0]).max()/MEV_J) for run in runs.values())
    maps=[]
    for col,label in enumerate(("narrow","nominal","wide")):
        run=runs[label];t=run["traces"]["time_s"]*1e12;keep=(t>=interval[0])&(t<=interval[1])
        axv.plot(t[keep],run["traces"]["Vout_relative_V"][keep]*1e3,color=colors[label],label=labels[label])
        spatial=run["spatial"];st=spatial["time_s"]*1e12;keep=(st>=interval[0])&(st<=interval[1])
        axd.plot(st[keep],spatial["deficit_J_m2"][keep]/peak,color=colors[label],label=labels[label])
        ax=fig.add_subplot(grid[1,2*col:2*col+2]);maps.append(ax)
        xy=run["snap"]["coordinates_m"].copy()
        xy[:,0]-=run["manifest"]["source"]["center_m"][0];xy*=1e9
        tri=mtri.Triangulation(xy[:,0],xy[:,1],run["snap"]["triangles"])
        index=diagnostic["map_indices"][label]
        artist=draw_complex_gap(ax,tri,run["snap"]["delta_J"][index],unit_J=MEV_J,
            vmin=0,vmax=vmax,bounds=(-60,60,-40,40))
        ax.set(xlim=(-60,60),ylim=(-40,40),aspect="equal",xlabel="x−xγ [nm]",
            title=labels[label]+"; σ="+number(run["manifest"]["parameter_identity"]["sigma_nm"])+" nm")
        if col==0:ax.set_ylabel("y [nm]")
    for ax in (axv,axd):
        ax.set(xlim=interval,xlabel="Tiempo postcascada [ps]");ax.grid(alpha=.18)
    axv.set_ylabel("Vout−Vbasal [mV]");axv.legend(fontsize=10)
    axd.set_ylabel("Déficit central / pico nominal")
    fig.colorbar(artist,ax=maps,label="|P1 Δ complejo| [meV]",shrink=.8)
    return save_figure(fig,build,"physical_width_sensitivity",records,units="mV,ps,meV,nm",
        reference="Actual physical sigma/sqrt2,nominal,sqrt2*sigma; own saved voltage baselines; nodal central deficit/reference actual peak in declared common interval",
        scope="Physical sensitivity only. Maps sample |complex spatial P1 Delta| on a nearest-cell raster (at most601 samples per axis) at one actual common saved time, without temporal interpolation. Display sampling does not establish core resolution, numerical convergence or experimental calibration.")

def convergence(comparisons,build,records):
    names=[];v=[];nodal=[];spatial_y=[];field=[];deficit=[];same_mesh_y=[]
    for row in comparisons:
        if "Vout" not in row or "gap_min" not in row:continue
        metric=comparison_admission_metrics(row)
        interval=row.get("common_trajectory_interval_ps",[0,0])
        label=row["kind"]+": "+row["reference"]+" → "+row["candidate"]
        index=len(names)
        names.append(textwrap.fill(label,45)+f'\n[{interval[0]:g}, {interval[1]:g}] ps')
        v.append(metric["Vout_margin_ratio"])
        if metric["nodal_diagnostic"] or metric["spectral_multiple_observables"]:
            spatial_y.append(index);field.append(metric["primary_gap_margin_ratio"])
            deficit.append(metric["deficit_margin_ratio"])
        if not metric["nodal_diagnostic"]:
            same_mesh_y.append(index);nodal.append(metric["nodal_margin_ratio"])
    if not names:raise ValueError("No physical waveform convergence metrics supplied")
    y=np.arange(len(names))
    fig,ax=plt.subplots(figsize=(11,max(4.4,.75*len(names)+1.6)),layout="constrained")
    ax.barh(y-.30,v,.18,color=COLORS["775"],label="Vout / margen")
    if spatial_y:
        ax.barh(np.array(spatial_y)-.10,field,.18,color="#13868b",label="L2 de |Δ| / margen declarado")
        ax.barh(np.array(spatial_y)+.10,deficit,.18,color="#bf7a16",label="Déficit/pico / margen declarado")
    if same_mesh_y:
        ax.barh(np.array(same_mesh_y)+.30,nodal,.18,color="#a14651",label="Mínimo nodal vigente / margen")
    ax.axvline(1,color=".3",ls=":",label="Límite declarado")
    ax.set(yticks=y,yticklabels=names,xlabel="Error observado / margen (adimensional)",
           xlim=(0,max(1.15,1.08*max(v+field+deficit+nodal))),
           title="Comparaciones sobre la ventana común realmente guardada")
    ax.invert_yaxis();ax.grid(axis="x",alpha=.18);ax.legend(fontsize=10)
    return save_figure(fig,build,"convergence",records,units="dimensionless,ps",
        reference="Vout error/max(reference pulse,threshold); mesh/domain:field/DC5% and deficit/actual peak5%; temporal:nodal/DC2%; explicit spectral policy supplies its actual field/deficit margins, with nodal diagnostic under exploratory5%. Physical interpolation is spatial complex P1 then magnitude then scalar time interpolation; complete original2%,3/2/2 and complex-time metrics stay visible.",
        scope="Saved axis-specific comparisons only; no latency inferred from map connectivity, no false nodal failure relabeled")


def tex_escape(value):
    mapping={"\\":r"\textbackslash{}","&":r"\&","%":r"\%","$":r"\$",
             "#":r"\#","_":r"\_","{":r"\{","}":r"\}","~":r"\textasciitilde{}"}
    return "".join(mapping.get(char,char) for char in str(value))

def number(value,digits=5):
    if value is None:return "--"
    return f"{float(value):.{digits}g}"

def event_text(summary):
    event=summary["event"]
    if event["confirmed"]:
        return "Cruce confirmado: "+number(event["latency_s"]*1e12)+" ps"
    if event["censored"]:
        return "Censurado: sin cruce en "+number(summary["final_time_ps"])+" ps"
    return "Sin cruce confirmado en la ventana guardada"

def belt_text(run):
    item=run["summary"]["hotbelt"]["first_saved_width_spanning"]
    if item is None:
        return ("Sin conexión transversal en los snapshots guardados "
                "(criterio |Delta| <= 0.2 |Delta_DC|). No se asigna latencia eléctrica por este criterio.")
    a,b=item["bracket_s"]
    return ("Conexión transversal guardada a "+number(item["first_saved_time_s"]*1e12)
            +" ps (entre "+number(a*1e12)+" y "+number(b*1e12)
            +" ps). Criterio diagnóstico |Delta| <= 0.2 |Delta_DC|; no es el gatillo eléctrico.")

def make_tex(admission,runs,comparisons,figures,gap_times,thickness,failure,vortex=None,sensitivity=None,dc_validation=None):
    date=tex_escape(admission.get("date","2026-10-11"))
    p=runs["775"]["manifest"]["physicalparams"]
    identity=runs["775"]["manifest"]["parameter_identity"]
    geometry=runs["775"]["manifest"]["geometry"]
    source=runs["775"]["manifest"]["source"]
    bulk=thickness["reference_film"]
    effective_Idep_A=float(bulk["depairing_current_A"])
    effective_bias_fraction=float(p["current_A"])/effective_Idep_A
    if effective_Idep_A<=0 or not math.isfinite(effective_bias_fraction):
        raise ValueError("Positive effective depairing current required from the thickness audit")
    text=[r"""\documentclass[10pt,a4paper]{article}
\usepackage[utf8]{inputenc}\usepackage[T1]{fontenc}
\usepackage{lmodern,amsmath,amssymb,graphicx,geometry,color,array,fancyhdr,hyperref}
\geometry{margin=17mm,top=19mm,bottom=18mm}
\definecolor{ink}{RGB}{20,52,72}\definecolor{accent}{RGB}{22,115,139}
\hypersetup{colorlinks=true,urlcolor=accent}
\pagestyle{fancy}\fancyhf{}
\fancyhead[L]{\small\textcolor{accent}{pySNSPD / cierre de etapa 5}}
\fancyhead[R]{\small """+date+r"""}\fancyfoot[L]{\scriptsize Datos físicos preservados / geminga}
\fancyfoot[R]{\thepage}\setlength{\headheight}{13pt}
\setlength{\parindent}{0pt}\setlength{\parskip}{6pt}
\newcommand{\head}[1]{\section*{\textcolor{ink}{#1}}}
\newcommand{\plot}[2][.99]{\begin{center}\includegraphics[width=#1\linewidth]{#2}\end{center}}
\newcommand{\captiontext}[1]{{\small #1}\par}
\begin{document}
\head{Informe final de la etapa 5}
{\Large Modelo efectivo acoplado para un SNSPD}\par
\textbf{Decisión.} Un único sistema de evolución y su balance independiente quedan
admitidos dentro del alcance y márgenes del certificado científico adjunto.
Las trayectorias muestran la respuesta posterior a la transferencia fotónica
hasta el gatillo declarado más margen, o su censura al techo de observación.
El reloj empieza en esa transferencia; la cascada óptica anterior no se simula.

\textbf{Hipótesis adoptada.} Estados BCS nativos, energía común del condensado y
su gradiente gauge, movilidad KWT calibrable, cinética conservativa electrón--fonón,
continuidad cuasineutra y los tres estados circuitales de la memoria.
La rigidez radial a $T_b$ fijo es nuestra extensión energética efectiva de la
corriente dirty BCS; no se atribuye a una evolución microscópica Usadel/Moyal.
\[
 U_{\rm total}=U_{\rm BCS}+U_\nabla+U_{\rm ph}+U_{\rm CM},\qquad
 \dot U_{\rm total}=P_{\rm fuente}-P_{\rm bias}-P_{\rm lectura}
 -P_{\rm sustrato}-P_{\rm cortes}.
\]
\begin{center}\begin{tabular}{p{.44\linewidth}p{.48\linewidth}}\hline
Magnitud o calibrador & Valor de estas pruebas\\\hline
"""]
    add=text.append
    def row(a,b):add(a+" & "+b+r"\\")
    row(r"$W,\ d_f,\ L$",number(geometry["width_m"]*1e9)+", "+
        number(p["film_thickness_m"]*1e9)+", "+number(geometry["length_m"]*1e9)+" nm")
    row(r"$h_{\rm malla},\ T_b,\ T_c$",number(identity["spacing_nm"])+" nm; "+
        number(p["bath_temperature_K"])+" K; "+number(p["Tc_K"])+" K")
    row(r"$I_{\rm DC},\ R_\square,\ D$",number(p["current_A"]*1e6)+" "+r"$\mu$A; "+
        number(p["sheet_resistance_ohm"])+r" $\Omega$; "+number(p["diffusion_m2_s"])+r" m$^2$/s")
    row(r"$I_{\rm dep}^{\rm efectivo},\ I_{\rm DC}/I_{\rm dep}^{\rm efectivo}$",
        number(effective_Idep_A*1e6)+r" $\mu$A; "+number(effective_bias_fraction,6))
    row(r"$\tau_{ee}(T_c),\ \tau_{ep}(T_c)$",
        number(p["tau_ee_Tc_s"]*1e12)+", "+number(p["tau_ep_Tc_s"]*1e12)+
        " ps; calibradores libres")
    row(r"$\tau_{0,ep}$",number(p["tau0_s"]*1e9)+
        r" ns; ligado a $\tau_{ep}(T_c)$, no independiente")
    row(r"$\tau_{\rm esc},\ \eta_{\rm ret},\ \sigma_\gamma$",
        number(p["tau_escape_s"]*1e12)+" ps; "+number(source["retained_fraction"])+"; "+
        number(source["sigma_nm"])+" nm")
    row(r"$L_{\rm ext}$",number(p["exterior_inductance"]["total_H"]*1e9)+" nH; "+
        "fijo desde el nuevo equilibrio DC")
    add(r"""\hline\end{tabular}\end{center}
La corriente reducida usa el límite de depairing de este cierre efectivo,
no la fracción del cierre Usadel histórico. No se cambió la corriente física al
cambiar de cierre. Los tiempos nominales de Allmaras no fijan el centro ni el rango de calibración.
El material acústico, escape de interfaz, retención y anchura fotónica son
aproximaciones exploratorias declaradas, no parámetros ya ajustados al experimento.
\begin{center}\begin{tabular}{lrrrr}\hline
Fotón & Energía retenida [eV] & Señal máxima [mV] & Gap mínimo [meV] & Ventana [ps]\\\hline
""")
    for wave,run in runs.items():
        sm=run["summary"]
        add(" & ".join([wave+" nm",number(sm["energy"]["retained_energy_J"]/EV_J),
             number(sm["metrics"]["Vout_peak_mV"]),number(sm["metrics"]["gap_min_meV"]),
             number(sm["final_time_ps"])])+r"\\")
    add(r"\hline\end{tabular}\end{center}")
    for wave,run in runs.items():
        add(r"\textbf{"+wave+" nm:} "+tex_escape(event_text(run["summary"]))+r".\par")
    if dc_validation is not None:add(dc_validation_latex(dc_validation))
    add(r"\textbf{Alcance temporal.} La latencia es relativa al estado postcascada guardado: "
        r"$\tau_{\rm rel}=t_{\rm cruce}-t_0$. La cascada anterior y posibles retardos de "
        r"lectura experimental aportan offsets todavía desconocidos; no se asigna una "
        r"latencia óptica absoluta. La censura sin gatillo describe esta hipótesis y ventana, "
        r"no demuestra éxito ni fracaso experimental del dispositivo.")
    add(r"\textbf{Uso inmediato.} Barridos deterministas del modelo efectivo y comparación experimental posterior; jitter estadístico y reset total requieren sus campañas propias.")
    def page(title):add(r"\newpage\head{"+title+"}")
    def plot(name):add(r"\plot{"+figures[name]+"}")
    def caption(value):add(r"\captiontext{"+value+"}")
    page("1. Del condensado a la señal de lectura")
    plot("overview")
    caption(r"\textbf{Figura 1.} A: tensión intrínseca de lectura respecto del basal DC guardado. "
        r"B: corrientes reales del detector y bias, sin inversión global. "
        r"C: mínimo nodal del gap. D: temperaturas equivalentes por energía de electrones y fonones. "
        r"Cada curva usa sus propios tiempos aceptados o snapshots; no se alisan los datos. "
        r"El máximo de lectura de 1550 nm es $0.33138\,\mu$V; su detalle se muestra "
        r"ampliado en la pipeline visual individual. Aquí el eje conserva el gatillo en mV.")
    add(r"\textbf{Parte del SNSPD.} La pérdida local de condensado redistribuye corriente hacia "
        r"el canal disipativo y el circuito de lectura. El cruce eléctrico y su confirmación "
        r"se leen del controlador de pasos aceptados; una imagen de gap deprimido no sustituye esa detección.")
    add(r"\textbf{Puerto y lectura.} $V_{\rm dev}$ es la caída en el hilo resuelto; "
        r"$V_{\rm out}=R_L(I_b-I_s)$ es la tensión de lectura. Los estados del circuito "
        r"y la inductancia exterior conservada redistribuyen el impulso: una caída "
        r"intrínseca de mV durante el movimiento de vórtices puede dar una lectura de "
        r"$\mu$V. No se identifica el voltaje del hilo con el gatillo de lectura.")
    obs=runs["775"]["manifest"]["observation"]
    add(r"\textbf{Ventana.} Umbral positivo "+number(obs["threshold_voltage_V"]*1e3)+
        " mV, confirmación "+number(obs["confirmation_time_s_required"]*1e12)+
        " ps después del cruce y parada en el máximo entre confirmación y cruce + "+number(obs["extra_time_s"]*1e12)+" ps; techo sin gatillo "+
        number(obs["no_trigger_ceiling_s"]*1e12)+" ps. Circuito completo, sin acortar sus constantes.")
    for wave,run in runs.items():
        page("2. Hotbelt: gap físico a "+wave+" nm")
        plot("gap_"+wave)
        caption(r"\textbf{Estados guardados [ps]:} "+
            ", ".join(number(t) for t in gap_times[wave])+
            r". Izquierda: dominio completo. Derecha: zoom central de 120 nm. "
            r"$|\Delta|$ mide la integridad del condensado: escala absoluta común en meV, "
            r"módulo del interpolante espacial complejo P1, rasterizado para visualizar. Sin interpolación temporal ni ventana térmica impuesta.")
        add(tex_escape(belt_text(run)))
    for wave in ("775","1550"):
        page("3. Electrones, fonones y expansión energética: "+wave+" nm")
        plot("maps_"+wave)
        caption(r"Filas: $T_e^{eq}$, $T_{ph}^{eq}$ y exceso local de energía completa "
            r"BCS+gradiente+fonones respecto del nuevo DC, dividido por volumen dual [J/m$^3$]. "
            r"Cada fila conserva su escala absoluta en todos los snapshots; "
            r"la ampliación central es sólo visual.")
        add(r"\textbf{Parte del SNSPD.} El depósito primero eleva la energía local; "
            r"el transporte expande la perturbación y las colisiones intercambian energía con "
            r"los fonones, que escapan al sustrato. Las temperaturas equivalentes etiquetan "
            r"esa energía: no fuerzan a las distribuciones evolutivas $p_\xi,n_\Omega$ a permanecer FD/BE.")
    page("4. Balance independiente: película, cortes y circuito")
    plot("energy")
    caption(r"\textbf{Figura 4.} Incrementos respecto del primer estado postdepósito. "
        r"Trabajo neto = integral de la potencia externa del mismo paso aceptado. "
        r"Defecto = $|U(t)-U(0)-\int P_{\rm ext}\,dt|/"
        r"[E_{\gamma,\rm ret}+\int|P_{\rm ext}|\,dt]$. Las potencias salientes conservan su signo.")
    add(r"\begin{center}\begin{tabular}{lrr}\hline "
        r"Fotón & Máximo defecto integrado [\%] & Máximo residuo instantáneo [W]\\\hline ")
    for wave,run in runs.items():
        e=run["summary"]["energy"]
        add(wave+" nm & "+number(e["maximum_relative_defect"]*100)+" & "+
            number(e["maximum_instant_power_residual_W"])+r"\\")
    add(r"\hline\end{tabular}\end{center}")
    add(r"\textbf{Parte del SNSPD.} El calor KWT y Joule se entrega una sola vez a "
        r"los portadores. El sustrato y los cortes superconductores son salidas físicas "
        r"explícitas. No se añade una fuente igual al residuo para cerrar el balance.")
    page("5. Convergencia y admisión numérica")
    plot("convergence")
    caption(convergence_caption_latex(comparisons))
    add(r"Los cuatro refinamientos corresponden a 775 nm y a los intervalos indicados. "
        r"1550 nm se presenta como trayectoria nominal exploratoria auditada en energía y soporte; "
        r"su contraste de paso tardío cubre 29.62--31.62 ps, sin auditoría independiente "
        r"de los cuatro ejes durante 0--200 ps.")
    add(r"\begin{center}{\small\begin{tabular}{lrrrrll}\hline "
        r"Eje & Vout & $L^2$/mín. & Déficit & Mín.nod. & Previo & Lat.\\\hline ")
    for kind in ("temporal","mesh","domain","spectral"):
        rows=[r for r in comparisons if r.get("kind")==kind and "Vout" in r]
        if not rows:continue
        metrics=[comparison_admission_metrics(r) for r in rows]
        vr=max(m["Vout_margin_ratio"] for m in metrics)
        gr=max(m["primary_gap_margin_ratio"] for m in metrics)
        deficits=[m["deficit_margin_ratio"] for m in metrics if m["deficit_margin_ratio"] is not None]
        dr=None if not deficits else max(deficits)
        nr=max(m["nodal_margin_ratio"] for m in metrics)
        prior=("cumple" if all(m["nodal_former_gate_passed"] for m in metrics) else "falla") if kind in {"mesh","domain"} or any(m["spectral_multiple_observables"] for m in metrics) else "vigente"
        if any(m["nodal_display_is_historical"] for m in metrics):prior=r"2\%/3-2-2 previo"
        available=sum(bool(r.get("latency_comparison_complete")) for r in rows)
        add(tex_escape(kind)+" & "+number(vr)+" & "+number(gr)+" & "+number(dr)+
            " & "+number(nr)+" & "+prior+" & "+str(available)+"/"+str(len(rows))+r"\\")
    add(r"\hline\end{tabular}}\end{center}")
    add(r"Las cifras son error/margen; el mínimo espectral exploratorio usa el margen previo de 3\% "
        r"sólo como diagnóstico. En malla y longitud el mínimo nodal y "
        r"su criterio previo permanecen visibles como diagnóstico: un núcleo puede "
        r"moverse entre nodos aunque el interpolante complejo tenga cero. La decisión "
        r"espacial se fijó antes del control fino con vórtices; no se renombra un fallo "
        r"previo como un pase de la política antigua. En espectro se muestran los criterios previos "
        r"y el presupuesto vigente explícito; un mapa de hotbelt no sustituye el gatillo eléctrico.")
    page("5b. Política temporal y trazabilidad del balance")
    historical_spectral=spectral_policy_history_latex(comparisons)
    if historical_spectral:add(historical_spectral)
    scope_table=finite_window_scope_latex(comparisons)
    if scope_table:add(scope_table)
    temporal_rows=temporal_refinement_display_rows(comparisons)
    if temporal_rows:
        add(r"\begin{center}{\small\begin{tabular}{llrl}\hline "
            r"Corrida & $\varepsilon_E$ & $h_{\min}/h_{\rm med}/h_{\max}$ [ps] & Intervalo [ps]\\\hline ")
        for item in temporal_rows:
            interval=item["post_anchor_interval_ps"]
            window="--" if interval is None else number(interval[0])+"--"+number(interval[1])
            for label,guard,steps in ((item["reference"],item["reference_guard"],item["reference_accepted_steps"]),
                                     (item["candidate"],item["candidate_guard"],item["candidate_accepted_steps"])):
                effective="--" if not steps.get("available") else "/".join(
                    number(steps[key]) for key in ("minimum_ps","median_ps","maximum_ps"))
                add(tex_escape(label)+" & "+number(guard)+" & "+effective+" & "+window+r"\\")
        add(r"\hline\end{tabular}}\end{center}")
        add(r"Los pasos efectivos corresponden al intervalo posterior a ambas anclas "
            r"de guarda. El control temporal usa una guarda ocho veces más estricta "
            r"para buscar aproximadamente la mitad del paso crítico; no presupone una "
            r"razón dos en cada paso adaptativo ni establece orden dos global.")
    add(r"\textbf{Solución no lineal.} Punto medio conjunto, eventos compartidos "
        r"limitados conservativamente y soporte $0\le p\le1$, $n\ge0$. "
        r"Normas $L^1$ energéticas para $p,n$, divididas por $E_{\gamma,\rm ret}$ "
        r"(independiente de malla y dominio), más guarda independiente de su fuerza inducida "
        r"en el gap; tolerancia admitida $2\times10^{-7}$, contrastada con $2\times10^{-9}$. El limitador es continuo "
        r"$\min(1,0.95\,{\rm disponible}/{\rm pérdida})$, con factor 1 si la pérdida es cero.")
    timeline=guard_timeline_latex(comparisons)
    if timeline:add(timeline)
    local_tolerance=float(identity["local_energy_tolerance"])
    temporal_values={r["candidate_guard"] for r in temporal_rows if r["candidate_guard"] is not None}
    temporal_label=(r"$"+guard_tolerance_latex(next(iter(temporal_values)))+"$" if len(temporal_values)==1
        else r"$\varepsilon_{E,\mathrm{ref}}$ declarado por segmentos")
    add(r"\textbf{Guarda temporal de energía.} Cada nuevo paso aceptado satisface "
        r"$|E_{k+1}-E_k-\Delta W_{k+1}|\le\varepsilon_E"
        r"(E_{\gamma,\rm ret}+\sum_{j\le k}|\Delta W_j|+|\Delta W_{k+1}|)$, "
        r"con $\varepsilon_E="+guard_tolerance_latex(local_tolerance)+r"$ desde la última ancla nominal y "+temporal_label+" "
        r"en el control temporal declarado. Se reconstruye de los datos "
        r"guardados desde el ancla declarada; no certifica retrospectivamente pasos legacy. "
        r"El balance acumulado independiente conserva el margen de 1\%.")
    add(r"Los recibos conservan la corrección del salto previo de 5\% y el contraste "
        r"entre normas de ocupación y de energía. Cambia la resolución del sistema implícito, "
        r"sin añadir fuentes físicas ni recortar poblaciones.")
    page("6. Qué significa h y qué cambia al adelgazar la película")
    plot("thickness")
    caption(r"\textbf{Figura 6.} Comparación bulk continua del mismo cierre, manteniendo "
        r"$D,N_0,T_c,T_b,W$. Las ramas sólo se muestran mientras son estables. "
        r"$h=5/3.5$ nm es resolución de malla; ambas conservan $d_f=7$ nm y sección 560 nm$^2$.")
    a=thickness["reference_film"];b=thickness["thinner_film"]
    add(r"\begin{center}\begin{tabular}{lrr}\hline "
        r"Magnitud & $d_f=7$ nm & $d_f=4.9$ nm\\\hline ")
    rows=[
        (r"$R_\square$ [$\Omega$]",a["sheet_resistance_ohm"],b["sheet_resistance_ohm"]),
        (r"$J(I_{\rm ref})$ [A/m$^2$]",a["operating_density_A_m2"],b["requested_density_A_m2"]),
        (r"$I_{\rm dep}$ [$\mu$A]",a["depairing_current_A"]*1e6,b["depairing_current_A"]*1e6)]
    for label,v1,v2 in rows:row(label,number(v1,7)+" & "+number(v2,7))
    add(r"\hline\end{tabular}\end{center}")
    add(r"\textbf{Respuesta a la duda física.} Reducir el espesor real 30\% "
        r"a la misma corriente aumenta $J$ en 42.86\%. En este cierre la película de 4.9 nm "
        r"queda "+number(b["fraction_above_depairing"]*100)+r"\% por encima de su límite bulk: "
        r"no existe la referencia superconductora estable solicitada. No se le asigna artificialmente gap cero. "
        r"Este límite es una predicción del modelo efectivo, no una medida experimental. "
        r"A 7 nm, $I_{\rm DC}/I_{\rm dep}^{\rm efectivo}="+number(effective_bias_fraction,6)+
        r"$; no debe confundirse con la fracción reducida histórica de otro cierre.")
    add(r"""\[
        \sigma_n=2e^2N_0D,\qquad
        R_\square=\frac{1}{\sigma_n d_f},\qquad I=W d_f J .
        \]""")
    add(r"Si se cambia espesor manteniendo también $R_\square$, cambia la conductividad "
        r"del material y ya no se conserva simultáneamente $D,N_0$. Refinar la malla "
        r"no realiza ninguno de esos cambios físicos.")
    page("7. Decisión física, procedencia y siguientes usos")
    add(r"\textbf{Ruta descartada, sin borrar evidencia.} La evolución de primer gradiente "
        r"Moyal no se adoptó: el refinamiento del primer paso no conservó soporte ni produjo "
        r"incrementos convergentes en regiones frías. El candidato efectivo se validó con "
        r"su propio DC, energía, evolución y márgenes; no hereda automáticamente aquella admisión.")
    if failure:
        add(r"\begin{center}\begin{tabular}{lrrr}\hline "
            r"Paso histórico [fs] & Máx. $|\delta d|$ & Máx. exceso Pauli & Máx. $\|R_1\|/\|R_0\|$\\\hline ")
        for step in failure.get("steps",[]):
            domains=step.get("domains",{}).values()
            dg=max((r["gap_increment_bar_quantiles_50_99_100"][-1] for r in domains),default=0)
            domains=step.get("domains",{}).values()
            pauli=max((r["Pauli_excess_max"] for r in domains),default=0)
            add(" & ".join([number(step["dt_s"]*1e15),number(dg),number(pauli),
                number(step.get("spectral_R1_over_R0_max"))])+r"\\")
        add(r"\hline\end{tabular}\end{center}")
    add(r"\textbf{Límites explícitos.} No se evoluciona $h_T$ ni la DOS Usadel "
        r"espacial: $\Gamma_{\rm DOS}=0$ y el depairing entra en la energía de gradiente. "
        r"La transferencia inicial alcanza $T_{E,eq}\simeq20$ K $>T_c=8.65$ K: "
        r"movilidad KWT y rigidez $\kappa(T_b)$ se adoptan como cierre fenomenológico calibrable, "
        r"sin derivación microscópica válida para todo ese rango y sus distribuciones. "
        r"Fonones acústicos tridimensionales efectivos, "
        r"escape de interfaz y depósito postcascada requieren calibración del dispositivo. "
        r"Los primeros resultados deterministas no dan por sí solos jitter estadístico "
        r"ni tiempos totales de recuperación.")
    add(r"\textbf{Proveniencia del material.} La banda acústica y el anclaje temporal "
        r"siguen Allmaras; se usa $\tau_{0,ep}=720\zeta(5)\tau_{ep}(T_c)/\pi^2$ "
        r"(tesis, p.\ 81). Las actividades S/R conservativas siguen Vodolazov, ecs.\ 3--4. "
        r"El archivo externo \texttt{nbn-a2f-ph.dat} de QNN/MIT se conserva con su hash "
        r"y metadatos incompletos; no se calculó con DFPT/Quantum ESPRESSO ni se usa "
        r"su normalización absoluta desconocida para estas tasas acústicas.")
    add(r"\textbf{Fuentes primarias.}\par"
        r"\href{https://doi.org/10.1103/PhysRev.108.1175}{BCS (1957): teoría del estado superconductor.}\par"
        r"\href{https://doi.org/10.1103/PhysRevLett.40.1041}{Kramer--Watts-Tobin (1978): dinámica efectiva del condensado.}\par"
        r"\href{https://arxiv.org/pdf/1611.06060}{Vodolazov (2017): colisiones ecs.\ 3--4 y corriente dirty BCS ec.\ 33.}\par"
        r"\href{https://thesis.caltech.edu/13748/}{Allmaras (2020): material acústico, p.\ 81, KWT y calibradores ecs.\ 3.24 y pp.\ 94--98.}\par"
        r"\href{https://github.com/qnngroup/proj-KE-solver}{QNN/MIT: archivo externo conservado, procedencia independiente de nuestro código.}\par")
    add(r"\textbf{Documentos matemáticos.} Las ecuaciones principales 2D, los funcionales "
        r"evaluados durante la evolución y las cantidades precomputadas se entregan por separado "
        r"en los tres PDFs de \texttt{Matematica}. El recibo del presente informe fija "
        r"los hashes de todos los datos utilizados.")
    if vortex:
        page("Vórtices y señal de lectura: observables distintos")
        add(r"\begin{center}\begin{tabular}{lrrrr}\hline ")
        add(r"Fotón & Mín. nodal [meV] & Mín. complejo [meV] & Déficit [\%] & Impulso [$\Phi_0$]\\\hline ")
        for wave,diag in vortex.items():
            bracket=diag["topology"]["first_saved_vortex_bracket_ps"]
            first="--" if bracket is None else number(bracket[0])+"--"+number(bracket[1])
            impulse=diag["electrical"]["CM_impulse_at_last_snapshot_Phi0"]
            deficit=100*diag["rows"][-1]["central_deficit_fraction_DC"]
            nodal=diag["rows"][-1]["minimum_nodal_gap_meV"]
            linear=diag["rows"][-1]["minimum_complex_linear_gap_meV"]
            add(wave+" nm & "+number(nodal)+" & "+number(linear)+" & "+number(deficit)+" & "+number(impulse)+r"\\")
        add(r"\hline\end{tabular}\end{center}")
        for wave,diag in vortex.items():
            bracket=diag["topology"]["first_saved_vortex_bracket_ps"]
            first="sin winding guardado" if bracket is None else "primer winding entre "+number(bracket[0])+" y "+number(bracket[1])+" ps"
            certified=diag["electrical"]["phase_unwrap"]["certified"]
            phase_status="unwrap relativo certificado" if certified else "unwrap ambiguo; no se afirma la identidad desenvuelta"
            add(tex_escape(wave+" nm: "+first+"; "+phase_status)+r".\par")
        add(r"\[ q_T=\frac1{2\pi}\sum_{\partial T}\operatorname{Arg}(\Delta_j\Delta_i^*), \qquad \int V_{\rm dev}\,dt=\frac{\Phi_0}{2\pi}\Delta(\theta_R-\theta_L). \]")
        selected=str(admission.get("report",{}).get("vortex_wavelength","775"))
        if selected not in vortex:raise ValueError("Unknown optional vortex figure wavelength")
        plot("phase_diagnostic_"+selected)
        caption(r"Triángulos orientados CCW; $\mathbf A=0$. Impulso reconstruido "
            r"desde los tres estados CM guardados y contrastado con las fases Josephson "
            r"de los cortes. Los trapecios de voltajes de extremos se conservan por separado: "
            r"su diferencia es cuadratura de análisis. Si el incremento relativo entre snapshots "
            r"puede alcanzar $\pi$, se declara unwrap ambiguo y sólo se usa el impulso CM. "
            r"No se redondea el impulso a slips enteros. "
            r"Una pareja de vórtices no implica hotbelt ni gatillo eléctrico.")
        scale=vortex[selected]["electrical"]["one_Phi0_short_time_delta_Vout_V"]*1e6
        add(r"La escala inductiva de un impulso $\Phi_0$ es $R_L\Phi_0/L_{\rm ext}="
            +number(scale)+r"\,\mu{\rm V}$; el circuito completo conserva sus tres estados "
            r"y la inductancia exterior fija.")
        page("Movimiento de núcleos y supresión espacial del condensado")
        add(r"\begin{center}\includegraphics[width=.84\linewidth]{"+figures["vortex_maps_"+selected]+r"}\end{center}")
        caption(r"Estados realmente guardados de "+selected+r" nm. "
            r"Símbolos: ceros del interpolante complejo con winding $\pm1$; "
            r"no resuelven por sí solos el radio microscópico del núcleo. "
            r"El mínimo nodal puede subir al desplazarse un núcleo entre nodos, "
            r"sin recuperación del condensado. Se conserva ese mínimo y se contrasta "
            r"con el déficit integrado central y el campo complejo. El recorte "
            r"visual de $120\times80$ nm no restringe la evolución térmica.")
        add(r"\[ D_{h,\rm nod}(t)=\sum_i w_i^{\mathcal R}[|\Delta_{{\rm DC},i}|-|\Delta_i(t)|]_+, \qquad \mathcal R=[-60,60]\times[-40,40]\ {\rm nm}^2. \]")
        add(r"$w_i^{\mathcal R}$ son pesos geométricos del interpolante lineal del déficit "
            r"positivo nodal; su suma es $9600$ nm$^2$. Es una aproximación espacial al déficit "
            r"continuo, distinta de integrar la magnitud del interpolante complejo. Se integra "
            r"exactamente sólo ese diagnóstico nodal elegido "
            r"sobre triángulos recortados a la región física común; $x$ se centra en el "
            r"depósito y $y$ conserva su coordenada física. Se publican las normalizaciones "
            r"por gap DC y área, por pico real y el error $L^2$ del campo cuando exista "
            r"el control correspondiente; no sustituyen automáticamente la admisión.")
        add(r"\href{https://arxiv.org/html/1611.06060v3}{Vodolazov (2017), sección V}, "
            r"distingue el movimiento de vórtices de la expansión de un dominio normal "
            r"cuando la disipación vence al enfriamiento. Es contexto para esta clausura "
            r"efectiva, no su validación microscópica.")
    if sensitivity:
        page("Sensibilidad física del depósito: anchura exploratoria")
        plot("physical_width_sensitivity")
        caption(r"Una sola energía retenida y ubicación, variando sólo $\sigma$; "
            r"curvas en la ventana física común declarada y mapas nodales al tiempo "
            +number(sensitivity["map_time_ps"])+r" ps realmente guardado en los tres casos. "
            r"Déficit central: $D_{h,\rm nod}$ en $120\times80$ nm dividido por el pico "
            r"nominal de esa ventana. El color comparte escala; no se interpola el tiempo. "
            r"$T_{\max}(0)$ es la temperatura equivalente inicial guardada: una anchura menor concentra "
            r"la misma energía, sin cambiar los calibradores ni ajustar una temperatura a posteriori.")
        add(r"\begin{center}{\small\begin{tabular}{lrrrp{.22\linewidth}rr}\hline "
            r"Caso & $\sigma$ [nm] & $E_{ret}$ [aJ] & $T_{\max}(0)$ [K] & Estado (completa) & $t_f$ [ps] & Ventana [ps]\\\hline ")
        for row in sensitivity["cases"]:
            interval=row["common_interval_ps"]
            state=tex_escape(row["status"]).replace(r"\_",r"\_\allowbreak{}")
            state+=" ("+("sí" if row["run_complete"] else "no")+")"
            add(tex_escape(row["label"])+" & "+number(row["sigma_nm"])+" & "+
                number(row["retained_energy_J"]*1e18)+" & "+number(row["maximum_initial_temperature_K"])+" & "+state+" & "+
                number(row["saved_final_time_ps"])+" & "+number(interval[0])+"--"+number(interval[1])+r"\\")
        add(r"\hline\end{tabular}}\end{center}")
        add(r"$t_0=0$ es transferencia postcascada en los tres casos; "
            r"$\tau_{ee}(T_c),\tau_{ep}(T_c)$ y $\eta_{ret}$ permanecen iguales. "
            r"La preparación Fermi--Dirac/Bose--Einstein a temperatura común redistribuye "
            r"energía entre electrones y fonones al variar $\sigma$: no es una sensibilidad "
            r"geométrica pura. Un menor $T_{\max}(0)$ no ordena por sí solo $V_{\rm out}$ "
            r"ni la formación de hotbelt, pues cambian el reparto y la extensión del gap "
            r"deprimido. El circuito de polarización puede aportar energía y amplificar "
            r"después la perturbación. "
            r"Es sensibilidad física exploratoria, separada de los controles numéricos y "
            r"de la calibración experimental. Una corrida de desarrollo de 20 ps sin cruce "
            r"no se declara censurada a 200 ps ni recibe latencia ficticia.")
        for row in sensitivity["cases"]:
            value=row["latency_in_this_window_ps"]
            status="sin gatillo confirmado en esta ventana" if value is None else "latencia relativa guardada "+number(value)+" ps"
            add(tex_escape(row["label"]+": "+status)+r".\par")
    add(r"\end{document}")
    return "\n".join(text)+"\n"

def json_scalar(value):
    """Preserve NumPy scalar values when serializing saved-data receipts."""
    import numpy as np
    if isinstance(value,np.generic):
        return value.item()
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")

def compile_report(build,tex_path,receipt):
    compiler=Path("/usr/bin/pdflatex")
    poppler=Path("/home/jdiaz/.conda/envs/snspd/bin")
    if not compiler.is_file():raise FileNotFoundError(compiler)
    environment=os.environ.copy()
    environment["PATH"]=str(compiler.parent)+os.pathsep+environment.get("PATH","")
    for name,leaf in (("TEXMFVAR","texlive-var"),("TEXMFCONFIG","texlive-config"),("TEXMFHOME","texmf")):
        environment.setdefault(name,str(SCRATCH/"tmp/pysnspd"/leaf))
    for attempt in (1,2):
        proc=subprocess.run([str(compiler),"-interaction=nonstopmode","-halt-on-error",
            "-output-directory="+str(build),str(tex_path)],cwd=build,env=environment,
            stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=180)
        (build/f"compile_{attempt}.log").write_bytes(proc.stdout)
        if proc.returncode:
            raise RuntimeError("LaTeX failed; inspect "+str(build/f"compile_{attempt}.log"))
    pdf=tex_path.with_suffix(".pdf")
    info=subprocess.run([str(poppler/"pdfinfo"),str(pdf)],
        check=True,stdout=subprocess.PIPE,text=True).stdout
    pages=int(next(line.split(":",1)[1].strip() for line in info.splitlines() if line.startswith("Pages:")))
    expected=int(receipt["planned_pages"])
    if pages!=expected:
        raise ValueError(f"Report has {pages} pages; repair layout before delivery (planned {expected})")
    log=tex_path.with_suffix(".log").read_text(errors="replace")
    warnings=[line for line in log.splitlines() if "Overfull" in line or "Underfull" in line]
    preview=build/"preview";preview.mkdir()
    subprocess.run([str(poppler/"pdftoppm"),"-scale-to","1600","-png",
        str(pdf),str(preview/"report")],check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=180)
    receipt.update(compiled=True,pdf_path=str(pdf),pdf_sha256=sha(pdf),pdf_pages=pages,
        layout_warnings=warnings,preview_dir=str(preview),
        visual_QA_approved=False,
        artifact_marker="Caller must follow its PDF artifact workflow before --compile; no success asserted by this script")

def repository_receipt(receipt,full_receipt_path):
    """Keep detailed numerical comparison data in scratch, not in Git."""
    keys=("schema","created_UTC","stage5","admission_path","admission_sha256",
          "builder_path","builder_sha256","input_sha256",
          "selected_gap_snapshot_times_ps","figures","tex_path","tex_sha256",
          "compiled","published","simulations_launched","resources","planned_pages",
          "scope","pdf_path","pdf_sha256","pdf_pages","layout_warnings",
          "visual_QA_approved","published_pdf","visual_QA_receipt")
    result={key:receipt[key] for key in keys if key in receipt}
    result["full_numerical_receipt"]={"path":str(full_receipt_path),
                                    "sha256":sha(full_receipt_path),
                                    "storage":"scratch"}
    return result

def publish_report(admission_path,build,target,qa_path):
    receipt_path=build/"final_report_receipt.json"
    receipt=read_json(receipt_path)
    if receipt.get("compiled") is not True:
        raise ValueError("Compile and visually review the report before publication")
    if receipt.get("admission_sha256")!=sha(admission_path):
        raise ValueError("Admission changed after compilation")
    if receipt.get("builder_sha256")!=sha(__file__):
        raise ValueError("Builder source changed after compilation")
    pdf=Path(receipt["pdf_path"])
    if sha(pdf)!=receipt["pdf_sha256"]:
        raise ValueError("Compiled PDF changed after its receipt")
    for path,digest in receipt["input_sha256"].items():
        if sha(path)!=digest:raise ValueError("Evidence changed after compilation: "+path)
    qa=read_json(qa_path)
    if qa.get("approved") is not True or qa.get("pdf_sha256")!=receipt["pdf_sha256"]:
        raise ValueError("A matching approved visual-QA receipt is required")
    target=inside(target,REPO/"output/pdf/implementation")
    target.parent.mkdir(parents=True,exist_ok=True)
    if target.exists() and sha(target)!=receipt["pdf_sha256"]:
        history=build/"previous_final_pdfs";history.mkdir(exist_ok=True)
        oldsha=sha(target);old=history/(target.stem+"_"+oldsha[:16]+".pdf")
        if not old.exists():shutil.copy2(target,old)
        if sha(old)!=oldsha:raise RuntimeError("Previous PDF preservation failed")
    shutil.copy2(pdf,target)
    if sha(target)!=receipt["pdf_sha256"]:raise RuntimeError("Published PDF checksum mismatch")
    receipt.update(published=True,published_pdf=str(target),visual_QA_approved=True,
                   visual_QA_receipt=str(Path(qa_path).resolve()))
    receipt_path.write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+"\n")
    target.with_suffix(".receipt.json").write_text(json.dumps(repository_receipt(receipt,receipt_path),indent=2,ensure_ascii=False,allow_nan=False,default=json_scalar)+"\n")
    print(json.dumps({"published_pdf":str(target),"sha256":receipt["pdf_sha256"],
                      "pages":receipt["pdf_pages"]}))

def load_physical_sensitivity(config,reference,*,input_hashes):
    """Use full original SHA-verified guard policy only when explicitly declared."""
    from sandbox.stage5_effective_validation import physical_sensitivity as ps
    input_hashes[str(Path(ps.__file__).resolve())]=sha(ps.__file__)
    protocol=config.get("bounded_guard_protocol")
    if protocol is None:
        return ps.load_sensitivity(config,reference,input_hashes=input_hashes,verify_entry=pinned_entry)
    if not isinstance(protocol,dict) or set(protocol)!={"path","sha256"} or len(protocol.get("sha256",""))!=64:
        raise ValueError("Physical sensitivity guard protocol requires explicit path and SHA")
    from sandbox.stage5_effective_validation import width_guard_window_candidate as width
    from sandbox.stage5_effective_validation import guard_protocol_dispatcher as dispatch
    from sandbox.stage5_effective_validation import audit_runs as audit
    protocol_path,protocol_sha=pinned_entry(protocol)
    evidence=dispatch.evidence_inputs(protocol_path)
    for filename,digest in evidence.items():
        checked,actual=pinned_entry({"path":filename,"sha256":digest})
        input_hashes[str(checked)]=actual
    proofs={}
    def validate(reference_run,candidate_run,label):
        entry=config["cases"][label]
        if not isinstance(entry.get("original_STARTED_log"),dict) or len(entry["original_STARTED_log"].get("sha256",""))!=64:
            raise ValueError("Original width STARTED log must be SHA-pinned")
        logfile,logsha=pinned_entry(entry["original_STARTED_log"])
        input_hashes[str(logfile)]=logsha;evidence[str(logfile)]=logsha
        proofs[label]=width.verify_width_case(reference_run,candidate_run,label,protocol_path,entry,
            config["common_interval_ps"],audit)
    result,runs,records=ps.load_sensitivity(config,reference,input_hashes=input_hashes,
        verify_entry=pinned_entry,parameter_validator=validate)
    for filename,digest in evidence.items():pinned_entry({"path":filename,"sha256":digest})
    result["declared_guard_window_policies"]=proofs
    result["original_identities_preserved"]=True
    for module in (width,dispatch,audit):
        input_hashes[str(Path(module.__file__).resolve())]=sha(module.__file__)
    records.append(dict(path=str(protocol_path),sha256=protocol_sha,
        scope="Physical sensitivity original0--20ps guard proof only; no numerical or stage5 admission"))
    return result,runs,records

def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--admission",type=Path,required=True)
    parser.add_argument("--build-root",type=Path,default=DEFAULT_BUILD)
    action=parser.add_mutually_exclusive_group()
    action.add_argument("--prepare-only",action="store_true",
        help="After admission, create PNGs and TeX only (default)")
    action.add_argument("--compile",action="store_true",
        help="After admission, create PNGs/TeX and compile/render PDF; caller runs artifact marker first")
    action.add_argument("--publish",action="store_true",
        help="Publish an existing compiled, visually approved PDF; no recompilation")
    parser.add_argument("--visual-qa",type=Path,
        help='For publication: JSON {"approved":true,"pdf_sha256":"..."} after visual inspection')
    parser.add_argument("--output-file",type=Path,default=DEFAULT_PDF)
    parser.add_argument("--cpu",type=int,default=13)
    args=parser.parse_args(argv)
    # Gate runs before any output directory, plot or PDF is created.
    admission_path=args.admission.expanduser().resolve()
    declared=read_json(admission_path)
    if declared.get("stage5") is not True:
        raise ValueError("Final report blocked: admission JSON must explicitly contain stage5: true")
    build=inside(args.build_root,SCRATCH)
    if args.publish:
        if args.visual_qa is None:raise ValueError("--publish requires --visual-qa")
        publish_report(admission_path,build,args.output_file,args.visual_qa)
        return
    if build.exists():
        raise FileExistsError("Build evidence is preserved; choose a fresh --build-root")
    affinity=os.sched_getaffinity(0)
    if args.cpu not in affinity:
        raise ValueError("Requested one-core affinity is unavailable")
    os.sched_setaffinity(0,{args.cpu})
    try:
        global np,plt
        import numpy as np
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from pysnspd.plotting.style import apply_thesis_style
        apply_thesis_style()
        plt.rcParams.update({"font.size":15,"axes.titlesize":16,"axes.labelsize":15,
                             "legend.fontsize":13,"savefig.facecolor":"white"})
        admission,runs,evidence,comparisons,input_hashes=admission_inputs(admission_path)
        input_hashes[str(COMPLEX_GAP_SOURCE)]=sha(COMPLEX_GAP_SOURCE)
        config=admission.get("report",{})
        thickness_path=Path(config.get("thickness_audit",str(
            SCRATCH/"stage5_thickness_bulk_audit_20261011/thickness_audit.json")))
        thickness=read_json(thickness_path)
        thickness_png=Path(thickness["outputs"]["png"]["path"])
        if sha(thickness_png)!=thickness["outputs"]["png"]["sha256"]:
            raise ValueError("Thickness figure differs from its physical audit")
        input_hashes[str(thickness_path.resolve())]=sha(thickness_path)
        input_hashes[str(thickness_png.resolve())]=sha(thickness_png)
        failure_path=Path(config.get("failure_audit",str(
            SCRATCH/"stage51_failure_audit_20261011/audit.json")))
        failure=read_json(failure_path)
        input_hashes[str(failure_path.resolve())]=sha(failure_path)
        dc_validation=None
        if admission.get("dc_validation") is not None:
            dc_validation,dc_record=load_dc_validation(admission["dc_validation"],runs["775"],input_hashes=input_hashes)
            evidence.append(dc_record)
        sensitivity=None;sensitivity_runs=None
        if admission.get("physical_sensitivity_runs") is not None:
            sensitivity,sensitivity_runs,sensitivity_records=load_physical_sensitivity(
                admission["physical_sensitivity_runs"],runs["775"],input_hashes=input_hashes)
            evidence.extend(sensitivity_records)
        build.mkdir(parents=True);(build/"figures").mkdir()
        records=[];figures={};gap_times={}
        figures["overview"]=overview(runs,build,records)
        for wave,run in runs.items():
            requested=config.get("snapshot_times_ps")
            if isinstance(requested,dict):requested=requested.get(wave)
            figures["gap_"+wave],gap_times[wave]=gap_maps(wave,run,build,records,requested)
            eqtimes=config.get("equivalent_times_ps")
            if isinstance(eqtimes,dict):eqtimes=eqtimes.get(wave)
            figures["maps_"+wave]=energy_maps(wave,run,build,records,eqtimes)
        figures["energy"]=energy_balance(runs,build,records)
        figures["convergence"]=convergence(comparisons,build,records)
        if sensitivity:
            figures["physical_width_sensitivity"]=physical_width_maps(sensitivity,sensitivity_runs,build,records)
        target=build/"figures/thickness.png";shutil.copy2(thickness_png,target)
        figures["thickness"]="figures/thickness.png"
        records.append({"path":str(target),"sha256":sha(target),
            "units":"uA,A/m2,meV,nm","reference":"Same material at physical thickness7/4.9nm",
            "scope":"Continuum effective bulk; mesh h5/3.5 changes no physical film thickness"})
        vortex=None
        if config.get("include_vortex_diagnostics",False):
            # Optional diagnostics remain strictly after the scientific gate.
            from sandbox.stage5_effective_validation import vortex_diagnostics as vd
            input_hashes[str(Path(vd.__file__).resolve())]=sha(vd.__file__)
            vortex={}
            selected=str(config.get("vortex_wavelength","775"))
            root=build/"vortex_diagnostics";root.mkdir()
            for wave,run in runs.items():
                entry=admission["primary_runs"][wave]
                runpath=Path(entry if isinstance(entry,str) else entry["path"]).resolve()
                diag,vm,vh,vs,vt=vd.analyze_run(wave,runpath,root/wave)
                vortex[wave]=diag
                for record in diag["source_hashes"].values():
                    source=record["source"];digest=record["sha256"]
                    if source in input_hashes and input_hashes[source]!=digest:
                        raise ValueError("Optional diagnostics captured changed admitted input: "+source)
                    input_hashes[source]=digest
                if wave==selected:
                    products=vd.plots(diag,vm,vh,vs,vt,root/wave)
                    for record in products:
                        src=Path(record["path"]);target=build/"figures"/(wave+"_"+src.name)
                        shutil.copy2(src,target)
                        name="phase_diagnostic_"+wave if src.stem=="phase_voltage_energy" else "vortex_maps_"+wave
                        figures[name]="figures/"+target.name
                        records.append({"path":str(target),"sha256":sha(target),
                            "units":"nm,ps,meV,rad,J/m3,mV,Phi0",
                            "reference":"Saved A=0 phase; exact CCW winding; common physical central window",
                            "scope":record["caption"]})
            if selected not in vortex:raise ValueError("Optional vortex wavelength absent")
        tex=make_tex(admission,runs,comparisons,figures,gap_times,thickness,failure,vortex,sensitivity,dc_validation)
        tex_path=build/"Final_stage5_effective.tex";tex_path.write_text(tex,encoding="utf8")
        receipt={"schema":"pysnspd.stage5_final_report.v1",
            "created_UTC":datetime.now(timezone.utc).isoformat(),"stage5":True,
            "admission_path":str(admission_path),"admission_sha256":sha(admission_path),
            "builder_path":str(Path(__file__).resolve()),"builder_sha256":sha(__file__),
            "input_sha256":input_hashes,"evidence_audits":evidence,
            "primary_run_summaries":{wave:run["summary"] for wave,run in runs.items()},
            "comparisons":comparisons,"selected_gap_snapshot_times_ps":gap_times,
            "figures":records,"tex_path":str(tex_path),"tex_sha256":sha(tex_path),
            "compiled":False,"published":False,"simulations_launched":0,
            "resources":{"workers":1,"CPU":[args.cpu],"BLAS_threads":1},
            "planned_pages":expected_report_pages(config,sensitivity is not None),
            "physical_sensitivity":sensitivity,"dc_validation":dc_validation,
            "scope":"Admitted effective neutral post-cascade photon windows; no experimental calibration, jitter statistics or full reset claimed"}
        if args.compile:compile_report(build,tex_path,receipt)
        (build/"final_report_receipt.json").write_text(json.dumps(
            receipt,indent=2,ensure_ascii=False,allow_nan=False,default=json_scalar)+"\n")
        print(json.dumps({"build_root":str(build),"compiled":receipt["compiled"],
            "figures":len(records),"published":False,"simulations_launched":0}))
    finally:
        os.sched_setaffinity(0,affinity)

if __name__=="__main__":
    main()
