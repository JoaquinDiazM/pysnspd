"""Bounded independent review of Simon phonon shapes; no detector transient.

All temperatures/energies below are conditional on the source axis being THz.
No absolute DOS normalization or atom/cell basis is certified by this script.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import platform
import sys
import time
import urllib.request

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import pysnspd.experimental.material_preprocessing as preprocessing
from pysnspd.experimental.material_preprocessing import derive_phonon_shape

PIN = "5b6bd747f80016da5ccd51db73c110a8ecc6abf6"
REPO = "https://github.com/qnngroup/proj-KE-solver"
PUBLIC_SHA = "9aeea0948033d771deedae40da0cb4dc59fef80ac6ea41ac4d3a67b180efc610"
LOCAL_SHA = "e94f11273e6c42b7c8fad978b78d119166e241bec1614ee77e4176906fa8f5df"
H_MEV_THz = 4.135667696
KB_MEV_K = 0.08617333262145
SOURCE_H_EV_PS = 2*np.pi*658.2e-6
T_GRID = np.array([.5,1,2,4,8.65,10,20,40,80,160,300,600,1000])
TOLERANCE = 1e-3  # Frozen independently before execution: acceptance_criteria.json.


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)+"\n", encoding="utf-8")


def integrate(x, y):
    return float(np.trapezoid(y,x))


def quadrature(x, order=8):
    z,w = np.polynomial.legendre.leggauss(order)
    width = np.diff(x)
    q=(x[:-1,None]+x[1:,None])/2+width[:,None]*z/2
    weights=width[:,None]*w/2
    return q.ravel(), weights.ravel()


def thermal_weights(nu, temperature):
    energy=nu*H_MEV_THz
    z=energy/(KB_MEV_K*temperature)
    # exp(-z) prevents overflow. Quadrature nodes never contain zero.
    expm=np.exp(-z)
    den=-np.expm1(-z)
    n=expm/den
    return energy*n, KB_MEV_K*z*z*expm/(den*den)


def compare(value, reference, *, direct_change=None):
    absolute=float(value-reference if direct_change is None else direct_change)
    # Signed DOS is unphysical; relative changes are not labelled acceptance.
    if reference <= 0:
        return {"absolute_change":absolute,"relative_change":None,"small_effect":False,
                "reason":"nonpositive_reference; no relative PASS"}
    relative=absolute/reference
    return {"absolute_change":absolute,"relative_change":float(relative),
            "small_effect":bool(abs(relative)<=TOLERANCE)}


def source_audit(source_dir, online):
    filenames=["nbn-a2f-ph.dat","solver.m","f.m","g.m"]
    result={"checked_online":False,"revision":PIN,"sources":{}}
    if online:
        source_dir.mkdir(parents=True,exist_ok=True)
        for filename in filenames:
            url=f"https://raw.githubusercontent.com/qnngroup/proj-KE-solver/{PIN}/{filename}"
            req=urllib.request.Request(url,headers={"User-Agent":"pysnspd-material-closure"})
            (source_dir/filename).write_bytes(urllib.request.urlopen(req,timeout=15).read())
        url="https://api.github.com/repos/qnngroup/proj-KE-solver/commits/main"
        req=urllib.request.Request(url,headers={"User-Agent":"pysnspd-material-closure"})
        head=json.loads(urllib.request.urlopen(req,timeout=15).read())
        result.update(checked_online=True,current_main_sha=head["sha"],
                      current_commit_date=head["commit"]["committer"]["date"])
        req=urllib.request.Request(f"https://api.github.com/repos/qnngroup/proj-KE-solver/git/trees/{head['sha']}?recursive=1",
                                  headers={"User-Agent":"pysnspd-material-closure"})
        tree=json.loads(urllib.request.urlopen(req,timeout=15).read())
        result["published_nbn_paths"]=[p["path"] for p in tree["tree"] if "nbn" in p["path"].lower()]
    for filename in filenames:
        path=source_dir/filename
        if path.exists():
            result["sources"][filename]={"sha256":sha(path.read_bytes()),
                "url":f"{REPO}/blob/{PIN}/{filename}"}
    if (source_dir/"nbn-a2f-ph.dat").exists():
        assert sha((source_dir/"nbn-a2f-ph.dat").read_bytes())==PUBLIC_SHA
    result["code_evidence"]={
        "solver.m:92-112":"reads absent nbn-a2f-ph_2.dat, assumes THz and states/THz; stable first duplicates; zero F and alpha2 where F<1e-8 states/eV",
        "column_order_solver.m:95-97":"a2F=data.Var2; F=data.PhDOS_states_THz_/p.h. This supports column 2 alpha2F, column 3 phonon DOS under the expected header; inversion was not adopted.",
        "f.m:3-5":"electronic kernel reconstructs interpolated alpha2 times interpolated F",
        "g.m:29-54":"phonon kernel uses alpha2 directly",
        "scope":"code supports preprocessing, not the normalization of the differently named public numeric file",
    }
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nbn-path",type=Path,required=True)
    parser.add_argument("--output-root",type=Path,default=ROOT/"docs/implementation/stage1_closure")
    parser.add_argument("--source-dir",type=Path,default=ROOT/"tmp/stage1_closure/material_sources")
    parser.add_argument("--online",action="store_true")
    args=parser.parse_args()
    started=time.perf_counter()
    out=args.output_root; figdir=out/"figures"
    figdir.mkdir(parents=True,exist_ok=True)
    source=source_audit(args.source_dir,args.online)
    raw=args.nbn_path.read_bytes()
    assert sha(raw) in {LOCAL_SHA,PUBLIC_SHA}
    data=np.loadtxt(io.BytesIO(raw),comments="#")
    payload=raw.split(b"\n",1)[1] if sha(raw)==LOCAL_SHA else raw
    assert sha(payload)==PUBLIC_SHA
    derived=derive_phonon_shape(args.nbn_path,source_url=f"{REPO}/blob/{PIN}/nbn-a2f-ph.dat",source_revision=PIN)
    source_floor=derive_phonon_shape(args.nbn_path,source_url=f"{REPO}/blob/{PIN}/nbn-a2f-ph.dat",
                                    source_revision=PIN,dos_floor=1e-8*SOURCE_H_EV_PS)
    x=derived.axis; table=data[np.argsort(data[:,0],kind="stable")]
    original=table.copy()
    first=data[derived.source_rows].copy()
    alpha0=first[:,1]; dos0=first[:,2]
    # Keep first duplicate is the selected policy; last/mean are sensitivity only.
    unique,inverse=np.unique(table[:,0],return_inverse=True)
    last=table[np.r_[np.diff(table[:,0])!=0,True]]
    mean=np.column_stack([unique]+[np.bincount(inverse,weights=table[:,j])/np.bincount(inverse) for j in [1,2]])
    origin=first.copy(); origin[origin[:,0]==0,1]=0
    clip=origin.copy(); clip[:,2]=np.maximum(clip[:,2],0)
    combined=np.column_stack([x,derived.alpha2F,derived.phdos])
    variants={"first_regular_origin":origin,"last_regular_origin":last.copy(),
              "mean_regular_origin":mean.copy(),"clip_F_only":clip,"common_support_v1":combined,
              "source_floor":np.column_stack([x,source_floor.alpha2F,source_floor.phdos])}
    for value in variants.values(): value[value[:,0]==0,1]=0

    q,w=quadrature(x,8); q16,w16=quadrature(x,16)
    observed={}
    for name,table in variants.items():
        observed[name]=(np.interp(q,table[:,0],table[:,1]),np.interp(q,table[:,0],table[:,2]))
    alpha_q,dos_q=observed["common_support_v1"]
    # Exact published algorithm: interpolate nodal ratio first and multiply by F.
    source_ratio=np.divide(source_floor.alpha2F,source_floor.phdos,out=np.zeros_like(x),where=source_floor.phdos>0)
    source_alpha_eff=np.interp(q,x,source_ratio)*np.interp(q,x,source_floor.phdos)
    observed["source_product_of_interpolants"]=(source_alpha_eff,np.interp(q,x,source_floor.phdos))
    thermal=[]
    for temperature in T_GRID:
        uweight,cweight=thermal_weights(q,temperature)
        values={name:{"U":float(np.sum(w*f*uweight)),"C":float(np.sum(w*f*cweight))}
                for name,(a,f) in observed.items()}
        raw_v=values["first_regular_origin"]
        delta_f=observed["common_support_v1"][1]-observed["first_regular_origin"][1]
        thermal.append({"T_K_assuming_THz":float(temperature),"values_original_DOS_units":values,
                        "derived_vs_signed_reference":{key:compare(values["common_support_v1"][key],raw_v[key],direct_change=np.sum(w*delta_f*weight)) for key,weight in [("U",uweight),("C",cweight)]}})

    nonthermal=[]
    for center in [4.,15.5,18.5,None,7.45]:
        sigma=.04 if center==7.45 else .4
        n=np.ones_like(q) if center is None else np.exp(-.5*((q-center)/sigma)**2)
        values={name:float(np.sum(w*f*q*H_MEV_THz*n)) for name,(a,f) in observed.items()}
        coupling_values={name:float(np.sum(w*H_MEV_THz*a*q*H_MEV_THz*n)) for name,(a,f) in observed.items()}
        nonthermal.append({"profile":"flat n=1" if center is None else f"Gaussian center={center:g} THz sigma={sigma:g} THz amplitude=1",
            "center_THz":center,"values_original_DOS_units":values,
            "diagnostic_added_after_discovering_internal_gap":center==7.45,
            "coupling_energy_moment_meV2":coupling_values,
            "coupling_moment_derived_vs_reference":compare(coupling_values["common_support_v1"],coupling_values["first_regular_origin"]),
            "derived_vs_signed_reference":compare(values["common_support_v1"],values["first_regular_origin"]),
            "signed_reference_is_physical":bool(values["first_regular_origin"]>0),
            "removed_signed_mass_over_absolute_reference_mass":float(np.sum(w*(observed["common_support_v1"][1]-observed["first_regular_origin"][1])*q*H_MEV_THz*n)/np.sum(w*np.abs(observed["first_regular_origin"][1])*q*H_MEV_THz*n))})
    moments={}
    for name,(alpha,dos) in observed.items():
        moments[name]={"lambda_regular_origin":float(2*np.sum(w*alpha/q))}
        moments[name].update({f"moment_E_power_{power}":float(np.sum(w*H_MEV_THz*alpha*(q*H_MEV_THz)**power)) for power in [0,1,2]})
    moment_changes={name:{key:compare(value,moments["first_regular_origin"][key]) for key,value in values.items()}
                    for name,values in moments.items() if name!="first_regular_origin"}
    exchange=[]
    for te in [10.,40.,160.,600.]:
        def bose(t): return np.exp(-q*H_MEV_THz/(KB_MEV_K*t))/(-np.expm1(-q*H_MEV_THz/(KB_MEV_K*t)))
        weight=H_MEV_THz*(q*H_MEV_THz)**2*(bose(te)-bose(4.))
        values={name:float(np.sum(w*a*weight)) for name,(a,f) in observed.items()}
        exchange.append({"Te_K":te,"Tph_K":4.,"definition":"normal-state diagnostic integral E^2 alpha2F [nB(Te)-nB(Tph)] dE; no prefactors/no transient",
            "values_meV3":values,"derived_vs_reference":compare(values["common_support_v1"],values["first_regular_origin"])})

    zeros=np.flatnonzero((dos0==0)&(alpha0>0)&(x>0))
    zero_starts=zeros[np.r_[True,np.diff(zeros)>1]]
    zero_ends=zeros[np.r_[np.diff(zeros)>1,True]]
    first_cut=int(np.flatnonzero((dos0<=0)&(x>0))[0])
    ratio=np.divide(derived.alpha2F,derived.phdos,out=np.zeros_like(x),where=derived.phdos>0)
    cutoff=32./H_MEV_THz
    finite_grid_lambda=float(2*np.trapezoid(alpha0[x>0]/x[x>0],x[x>0]))
    regular_lambda=moments["first_regular_origin"]["lambda_regular_origin"]
    ratio_checks={
        "raw_has_nonintegrable_pole":bool(len(zeros)),
        "zero_DOS_nonzero_alpha_node_count":int(len(zeros)),
        "zero_DOS_nonzero_alpha_intervals":[{"frequency_range_THz_assumed":[float(x[i]),float(x[j])],"energy_range_meV_assumed":[float(x[i]*H_MEV_THz),float(x[j]*H_MEV_THz)],"alpha2F_at_start":float(alpha0[i])} for i,j in zip(zero_starts,zero_ends)],
        "raw_limit":"alpha2F/F diverges at the first positive-frequency zero of F; no finite global kernel admitted",
        "derived_ratio_bounded_by_max_positive_node":float(np.max(ratio)),
        "derived_ratio_max_below_32meV":float(np.max(ratio[x<=cutoff])),
        "derived_ratio_unit":"original alpha2F/DOS ordinate ratio; energy units not certified",
        "support_first_zero_THz_assumed":float(x[first_cut]),
        "support_first_zero_meV_assumed":float(x[first_cut]*H_MEV_THz),
        "minimum_signed_DOS":float(np.min(dos0)),
        "kernel_admission":"RESTRICTED_SHAPE_ONLY; finite interpolation after explicit common-support cutoff; no global physical rate certified",
        "ratio_identity_max_residual":float(np.max(np.abs(alpha_q-np.divide(alpha_q,dos_q,out=np.zeros_like(q),where=dos_q>0)*dos_q))),
    }
    low=q<=cutoff
    source_low=np.max(np.abs(source_alpha_eff[low]-alpha_q[low]))
    raw16=np.interp(q16,x,origin[:,1]); der16=np.interp(q16,x,derived.alpha2F)
    qcheck={"lambda_8_vs_16_relative_reference":float(abs(2*np.sum(w16*raw16/q16)-regular_lambda)/regular_lambda),
            "lambda_8_vs_16_relative_derived":float(abs(2*np.sum(w16*der16/q16)-moments["common_support_v1"]["lambda_regular_origin"])/moments["common_support_v1"]["lambda_regular_origin"])}
    operation_metrics={
        "signed_DOS_area_original_order":integrate(data[:,0],data[:,2]),
        "signed_DOS_area_sorted":integrate(original[:,0],original[:,2]),
        "signed_DOS_area_keep_first":integrate(x,dos0),
        "signed_DOS_area_keep_last":integrate(last[:,0],last[:,2]),
        "signed_DOS_area_mean_duplicates":integrate(mean[:,0],mean[:,2]),
        "derived_DOS_area":integrate(x,derived.phdos),
        "negative_DOS_area_removed":-integrate(x,np.minimum(dos0,0)),
        "negative_rows":int(np.sum(data[:,2]<0)),
        "origin_alpha_area_removed":float(alpha0[0]*x[1]/2),
        "alpha_area_removed_common_support_original_axis":integrate(x,origin[:,1]-derived.alpha2F),
        "alpha_area_removed_fraction":integrate(x,origin[:,1]-derived.alpha2F)/integrate(x,origin[:,1]),
        "finite_positive_grid_lambda_reference":finite_grid_lambda,
        "source_solver_lambda_parameter":1.2,
        "lambda_match_cannot_verify_DOS_units":True,
        "frequency_as_THz_DOS_as_states_meV_hypothesis_modes_raw":integrate(x,dos0)*H_MEV_THz,
        "frequency_as_THz_DOS_as_states_meV_hypothesis_modes_clipped":integrate(x,derived.phdos)*H_MEV_THz,
        "hypothesis_is_adopted":False,
        "first_negative_THz_assumed":float(x[np.flatnonzero(dos0<0)[0]]),
        "last_positive_THz_assumed":float(x[np.flatnonzero(dos0>0)[-1]]),
        "source_product_interpolants_vs_common_interpolant_max_absolute_alpha_below32meV":float(source_low),
    }
    # The negative tail is outside 32 meV, but a zero-DOS acoustic/optical gap
    # begins below 32 meV and DOES change coupling inside that reaction window.
    safe_cutoff=30./H_MEV_THz
    assert np.array_equal(origin[(x>0)&(x<=safe_cutoff),1],derived.alpha2F[(x>0)&(x<=safe_cutoff)])
    assert np.array_equal(dos0[x<=cutoff],derived.phdos[x<=cutoff])
    low_changes={}
    for limit in [30.,32.]:
        ql,wl=quadrature(np.r_[x[x<limit/H_MEV_THz],limit/H_MEV_THz])
        ar=np.interp(ql,x,origin[:,1]); ad=np.interp(ql,x,derived.alpha2F)
        low_changes[str(limit)]={f"M{power}_relative_change":float(np.sum(wl*(ad-ar)*ql**power)/np.sum(wl*ar*ql**power)) for power in [0,1,2]}
    assert np.max(np.abs(qcheck["lambda_8_vs_16_relative_derived"]))<1e-6
    assert args.nbn_path.read_bytes()==raw
    csv_path=out/"material_nbn_shape_v1.csv"
    np.savetxt(csv_path,np.column_stack([x,derived.alpha2F,derived.phdos,derived.source_rows]),delimiter=",",
               header="frequency_original_axis,alpha2F_derived,phonon_DOS_original_ordinate_units,source_row_zero_based",comments="",fmt=["%.17g","%.17g","%.17g","%d"])
    manifest=derived.manifest.copy()
    manifest.update(derived_csv_sha256=sha(csv_path.read_bytes()),derived_csv=csv_path.name,
                    source_numeric_payload_sha256=PUBLIC_SHA,
                    valid_experimental_use="shape sensitivity; common-support changes absent below 30 meV apart from origin regularity. A zero-DOS gap affects the 30.47-32 meV part of the R2 reaction window; not a globally admitted kinetic kernel. Conditional THz axis.",
                    absolute_SI_admission="BLOCKED: table-specific DOS units, atom/cell basis and matching number density unverified")
    dump(out/"material_nbn_shape_v1.manifest.json",manifest)
    result={
        "schema":"pysnspd.material-closure.v1","primary_sources":source,
        "status":"CONDITIONAL_SHAPE_ONLY; global physical phonon kinetics and absolute SI material remain unadmitted",
        "unit_scope":"Temperatures and meV are conditional on THz axis. U/C use the original DOS ordinate, without ion density; no absolute heat capacity prediction.",
        "threshold_small_relative_effect":TOLERANCE,
        "source_sha256":sha(raw),"source_numeric_payload_sha256":PUBLIC_SHA,
        "derived_csv_sha256":sha(csv_path.read_bytes()),
        "module_sha256":sha(Path(preprocessing.__file__).read_bytes()),"runner_sha256":sha(Path(__file__).read_bytes()),
        "operations":operation_metrics,"thermal":thermal,"nonthermal":nonthermal,
        "alpha_moments":moments,"alpha_moment_changes":moment_changes,"normal_exchange":exchange,
        "B39_phonon_injection":{
            "energy_fraction_retained_if_original_normalizer_reused":moments["common_support_v1"]["moment_E_power_1"]/moments["first_regular_origin"]["moment_E_power_1"],
            "energy_fraction_with_consistent_derived_normalizer":float(np.sum(w*H_MEV_THz*q*H_MEV_THz*alpha_q)/moments["common_support_v1"]["moment_E_power_1"]),
            "interpretation":"B.39 must recompute its alpha2F energy moment after preprocessing. This normalizes injected energy, not the DOS mode count; it does not certify absolute phonon occupations.",
        },
        "inverse_DOS_kernel":ratio_checks,"quadrature_check":qcheck,
        "electronic_reaction_scope":{"assumed_closed_phonon_window_meV":[0.,32.],
            "negative_DOS_tail_inside_window":False,"common_support_alpha_edit_inside_window":True,
            "origin_alpha_regularization_inside_window":True,"unmodified_common_support_below_meV":30.,
            "finite_window_alpha_moment_changes":low_changes,
            "warning":"Incoming phonons above the closed window can create electronic states outside the catalogue. No statement about the full photon cascade or unrestricted transient is made."},
        "platform":platform.platform(),"python":platform.python_version(),"numpy":np.__version__,
    }
    plot_shapes(x,alpha0,dos0,derived,ratio,figdir)
    plot_impacts(thermal,nonthermal,moment_changes,figdir)
    result["runtime_seconds"]=time.perf_counter()-started
    dump(out/"material_results.json",result)
    rows=[]
    for row in thermal:
        t=row["T_K_assuming_THz"]
        for name,vals in row["values_original_DOS_units"].items(): rows.append([t,name,vals["U"],vals["C"]])
    with (out/"material_thermal.csv").open("w",newline="",encoding="utf-8") as stream:
        writer=csv.writer(stream);writer.writerow(["T_K_conditional_THz","variant","U_original_DOS_units","C_original_DOS_units"]);writer.writerows(rows)
    print(json.dumps({"status":result["status"],"runtime_seconds":result["runtime_seconds"],"operations":operation_metrics,"inverse_DOS_kernel":ratio_checks},indent=2))


def plot_shapes(x,alpha,dos,derived,ratio,out):
    plt.rcParams.update({"font.size":11.5,"axes.titlesize":12,"axes.labelsize":11.5,"legend.fontsize":11,"xtick.labelsize":11,"ytick.labelsize":11})
    fig,axs=plt.subplots(2,2,figsize=(10,6.4),layout="constrained")
    a,b,c,d=axs.ravel()
    a.plot(x,dos,color="#293f54",lw=1.5,label="Simón: tabla publicada")
    a.plot(x,derived.phdos,color="#db802b",lw=1.8,ls="--",label="Derivado v1")
    a.axvspan(0,32/H_MEV_THz,color="#6bb5a0",alpha=.15,label="Ventana ≤32 meV")
    a.set(xlabel="Frecuencia [THz, supuesto]",ylabel="DOS fonónica [ordenada original]",title="a) Conservación de la forma y recorte")
    a.legend(loc="upper right",framealpha=.94)
    b.plot(x,dos,color="#293f54",lw=1.5,label="Publicada")
    b.plot(x,derived.phdos,color="#db802b",ls="--",lw=2,label="Derivada")
    b.axhline(0,color=".5",lw=.7); b.set(xlim=(16.1,20),ylim=(-.007,.011),xlabel="Frecuencia [THz, supuesto]",ylabel="DOS fonónica",title="b) Zoom de la cola negativa")
    b.legend(loc="upper right")
    c.plot(x,alpha,color="#293f54",lw=1.5,label="Bruta")
    c.plot(x,derived.alpha2F,color="#db802b",lw=1.8,ls="--",label="Derivada")
    c.fill_between(x,derived.alpha2F,alpha,color="#cc4455",alpha=.3,label="Área retirada")
    c.set(xlim=(7.25,7.6),ylim=(-.007,.18),xlabel="Frecuencia [THz, supuesto]",ylabel=r"$\alpha^2F$",title="c) También hay F=0 cerca de 30,47 meV")
    c.legend(loc="upper right")
    mask=(dos>0)&(x>0)
    rawratio=np.full_like(x,np.nan); rawratio[mask]=alpha[mask]/dos[mask]
    d.semilogy(x,rawratio,color="#293f54",lw=1.3,label="Bruto (F>0)")
    d.semilogy(x, np.where(ratio>0,ratio,np.nan),color="#db802b",lw=1.5,ls="--",label="Derivado: coinciden")
    d.set(xlim=(0,17),ylim=(1e-5,1e6),xlabel="Frecuencia [THz, supuesto]",ylabel=r"$\alpha^2F/F$ [unidades originales]",title="d) Polos del cociente bruto en los bordes")
    d.axvline(16.3917,color="#cc4455",lw=.9)
    d.axvline(7.36881,color="#cc4455",lw=.9)
    d.legend(loc="lower center")
    for ax in axs.ravel(): ax.grid(alpha=.2)
    fig.savefig(out/"material_support_review.png",dpi=180); fig.savefig(out/"material_support_review.pdf"); plt.close(fig)


def plot_impacts(thermal,nonthermal,moment_changes,out):
    fig,axs=plt.subplots(1,3,figsize=(10,4.5),layout="constrained",gridspec_kw={"width_ratios":[1.2,1.1,1.]})
    a,b,c=axs
    t=np.array([v["T_K_assuming_THz"] for v in thermal])
    for key,label,style in [("U","Energía U","-"),("C","Capacidad C","--")]:
        y=np.array([v["derived_vs_signed_reference"][key]["relative_change"] for v in thermal])*100
        a.semilogx(t,y,style,marker="o",ms=3,label=label)
    a.axhline(.1,color="#cc4455",lw=1,ls=":",label="Criterio 0,1 %")
    a.set(xlabel="Temperatura [K, supuesto THz]",ylabel="Cambio respecto a la tabla [%]",title="a) Sensibilidad térmica")
    a.legend(loc="upper left")
    labels=["4","15,5","18,5","Plano","7,45*"]
    yy=[100*v["removed_signed_mass_over_absolute_reference_mass"] for v in nonthermal]
    coupling=[100*abs(v["coupling_moment_derived_vs_reference"]["relative_change"]) for v in nonthermal]
    b.bar(np.arange(5)-.18,yy,width=.36,color="#318878",label="Energía U")
    b.bar(np.arange(5)+.18,coupling,width=.36,color="#db802b",label="Acoplamiento")
    b.set_xticks(range(5),labels,rotation=30,ha="right")
    b.set(ylabel="Sensibilidad [%]",xlabel="Centro [THz]; *ancho 0,04",title="b) Perfiles no térmicos",ylim=(0,140))
    b.legend(loc="upper left")
    changes=moment_changes["common_support_v1"]
    keys=["lambda_regular_origin","moment_E_power_0","moment_E_power_1","moment_E_power_2"]
    yy=[100*abs(changes[k]["relative_change"]) for k in keys]
    c.bar(range(4),yy,color="#db802b")
    c.axhline(.1,color="#cc4455",lw=1,ls=":")
    c.set_xticks(range(4),[r"$\lambda$",r"$M_0$",r"$M_1$",r"$M_2$"])
    c.set(ylabel="Magnitud del cambio [%]",title="c) Acoplamiento global")
    c.text(.03,.97,r"$M_k=\int E^k\alpha^2F\,dE$",transform=c.transAxes,va="top",fontsize=11)
    c.set_ylim(0,max(yy)*1.4)
    for ax in axs: ax.grid(axis="y",alpha=.2)
    fig.savefig(out/"material_impact_review.png",dpi=180); fig.savefig(out/"material_impact_review.pdf"); plt.close(fig)


if __name__=="__main__": main()
