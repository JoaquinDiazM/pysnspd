"""Thickness versus mesh-spacing audit of the admitted effective BCS bulk.

All computations use SI. The figures display physical meV, microampere and
A/m^2. The thin-film comparison changes only physical thickness and therefore
sheet resistance; D, N0, Tc, Tb and width remain unchanged.
"""
from __future__ import annotations
import argparse
from dataclasses import asdict
import csv
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.constants import Boltzmann as kB, elementary_charge as e
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pysnspd.experimental.effective_bcs_functional import (
    EffectiveBCSFunctional, bcs_count_quadrature)
from pysnspd.experimental.effective_bcs_preparation import get_effective_geometry

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=Path(
        "/home/jdiaz/scratch/stage5_thickness_bulk_audit_20261011"))
    args=parser.parse_args()
    root=args.output.resolve()
    if not str(root).startswith("/home/jdiaz/scratch/"):
        raise ValueError("Audit outputs must remain in /home/jdiaz/scratch")
    root.mkdir(parents=True,exist_ok=True)
    start=time.monotonic()
    Tc,Tb,D,N0=8.65,.9,5e-5,9.153297821093998e46
    width,film,thin=80e-9,7e-9,4.9e-9
    current=21.29720614489798e-6
    rs=608.
    xi,w=bcs_count_quadrature(cutoff_bar=32.,order=6)
    geom=get_effective_geometry(160e-9,5e-9)
    def make_model(thickness):
        return EffectiveBCSFunctional(
            geom.graph,energy_scale_J=kB*Tc,N0_per_J_m3=N0,
            ell0_m=geom.ell0_m,thickness_m=thickness,
            count_nodes=xi,count_weights=w,gradient_law="bcs_radial",
            bath_temperature_K=Tb,diffusion_m2_s=D)
    model=make_model(film)
    reference=model.bulk(temperature_K=Tb,current_A=current,width_m=width)
    qaxis=np.linspace(0.,reference.depairing_q_bar,201)
    sampled=np.array([model.bulk_current(Tb,q,width) for q in qaxis])
    I7,abar=sampled.T
    Ithin=I7*thin/film
    density=I7/(width*film)
    gap_J=abar*model.Eu_J
    gap_meV=gap_J/e*1e3
    idep7=reference.depairing_current_A
    idepthin=idep7*thin/film
    thin_model=make_model(thin)
    scaled_current,scaled_amplitude=thin_model.bulk_current(Tb,reference.q_bar,width)
    np.testing.assert_allclose(scaled_current,reference.current_A*thin/film,rtol=2e-14)
    np.testing.assert_allclose(scaled_amplitude,reference.amplitude_bar,rtol=0,atol=0)
    assert current>idepthin, "Update the admission explanation if the model changes."
    # Above the maximum stable branch, a zero gap would be an invented DC value:
    # report no admissible superconducting reference, not an assumed normal state.
    thin_admission={"admitted":False,"gap_J_at_requested_current":None,
                    "reason":"requested current exceeds this effective model depairing current"}
    mesh_rows=[]
    for h in (5e-9,3.5e-9):
        geometry=get_effective_geometry(160e-9,h)
        mesh_rows.append({
            "requested_mesh_spacing_m":h,"physical_film_thickness_m":film,
            "cross_section_m2":width*film,
            "nodes":len(geometry.coordinates_m),
            "geometry_path":geometry.source_path,
            "geometry_sha256":sha(geometry.source_path)})
    csvpath=root/"stable_bulk_branches.csv"
    with csvpath.open("w",newline="") as stream:
        writer=csv.writer(stream)
        writer.writerow(["q_m_inv","gap_J","gap_meV","current_7nm_A",
                         "current_4p9nm_A","current_density_A_m2"])
        writer.writerows(zip(qaxis/geom.ell0_m,gap_J,gap_meV,I7,Ithin,density))
    plt.rcParams.update({"font.size":11,"axes.titlesize":12,"axes.labelsize":12,
                         "savefig.facecolor":"white","font.family":"DejaVu Sans"})
    fig,axs=plt.subplots(1,2,figsize=(12.6,5.6),layout="constrained")
    color7,color49="#1764a5","#bf4b16"
    axs[0].plot(I7*1e6,gap_meV,color=color7,lw=2.4,label=r"$d_f=7.0$ nm")
    axs[0].plot(Ithin*1e6,gap_meV,color=color49,lw=2.4,label=r"$d_f=4.9$ nm")
    axs[0].axvline(current*1e6,color=".25",ls=":",lw=1.6,label=r"$I_{ref}=21.297206$ µA")
    axs[0].scatter([current*1e6],[reference.delta_J/e*1e3],
                   color=color7,s=40,zorder=4)
    axs[0].scatter([idep7*1e6,idepthin*1e6],[gap_meV[-1]]*2,
                   color=[color7,color49],s=38,zorder=4)
    axs[0].set(xlabel=r"Corriente $I$ (µA)",ylabel=r"Gap $|\Delta|$ (meV)",
               title="Mismo material; cambia el espesor físico")
    axs[0].legend(loc="upper right",fontsize=10,framealpha=.95)
    axs[0].text(.04,.06,
        "Límite de la rama estable:\n"
        f"7.0 nm: {idep7*1e6:.5f} µA\n"
        f"4.9 nm: {idepthin*1e6:.5f} µA < Iref\n"
        "4.9 nm: sin referencia DC estable a Iref",
        transform=axs[0].transAxes,fontsize=10,va="bottom",
        bbox={"boxstyle":"round,pad=0.5","facecolor":"white","edgecolor":".8"})
    axs[1].plot(density,gap_meV,color=color7,lw=3,label="7.0 nm")
    axs[1].plot(density,gap_meV,color=color49,lw=1.8,ls="--",label="4.9 nm: misma curva")
    J7=current/(width*film)
    Jthin=current/(width*thin)
    axs[1].axvline(J7,color=color7,ls=":",lw=1.6,
                   label=r"$J(I_{ref},7.0\,nm)$")
    axs[1].axvline(Jthin,color=color49,ls=":",lw=1.6,
                   label=r"$J(I_{ref},4.9\,nm)$")
    axs[1].scatter([J7],[reference.delta_J/e*1e3],color=color7,s=40,zorder=4)
    axs[1].set(xlabel=r"Densidad de corriente $J$ (A/m$^2$)",
               ylabel=r"Gap $|\Delta|$ (meV)",
               title="La depresión del gap sigue la densidad de corriente")
    axs[1].ticklabel_format(axis="x",style="sci",scilimits=(0,0),useMathText=True)
    axs[1].legend(loc="upper right",fontsize=9.6,framealpha=.95)
    for ax in axs:
        ax.grid(alpha=.18)
        ax.set_ylim(bottom=min(gap_meV)-.035,top=max(gap_meV)+.025)
    fig.suptitle(
        "BCS efectivo con rigidez radial: $T_b=0.9$ K, $T_c=8.65$ K, $W=80$ nm",
        fontsize=14)
    fig.supxlabel(
        "h = 5 o 3.5 nm es el espaciado de malla; ambos casos conservan d_f = 7 nm.\n"
        "Predicción del cierre efectivo; el límite de corriente no es un ajuste experimental.",
        fontsize=10)
    png=root/"gap_current_thickness.png"
    pdf=root/"gap_current_thickness.pdf"
    fig.savefig(png,dpi=170)
    fig.savefig(pdf)
    plt.close(fig)
    result={
        "schema":"pysnspd.effective_bcs.thickness_audit.v1",
        "status":"BULK_PHYSICAL_THICKNESS_AUDIT_NOT_TRANSIENT_ADMISSION",
        "model_contract":reference.contract,
        "theory":"BCS native thermal stationary branch with our shared-energy dirty-BCS radial stiffness",
        "stiffness_primary_reference":"https://arxiv.org/pdf/1611.06060 Eq. (33)",
        "radial_stiffness_is_our_effective_energy_hypothesis":True,
        "fixed_material_parameters":{
            "Tc_K":Tc,"Tb_K":Tb,"D_m2_s":D,"N0_per_J_m3":N0,
            "width_m":width,"normal_conductivity_S_m":1/(rs*film),
            "xi_cutoff_bar":32.,"quadrature_order":6,"xi_nodes":len(xi)},
        "reference_film":{
            "thickness_m":film,"sheet_resistance_ohm":rs,
            "cross_section_m2":width*film,
            "operating_current_A":current,"operating_density_A_m2":J7,
            "operating_gap_J":reference.delta_J,
            "operating_gap_meV":reference.delta_J/e*1e3,
            "depairing_current_A":idep7,
            "depairing_density_A_m2":idep7/(width*film),
            "bulk_reference":asdict(reference)},
        "thinner_film":{
            "thickness_m":thin,"sheet_resistance_ohm":rs*film/thin,
            "cross_section_m2":width*thin,
            "requested_current_A":current,"requested_density_A_m2":Jthin,
            "depairing_current_A":idepthin,
            "requested_current_over_depairing":current/idepthin,
            "fraction_above_depairing":current/idepthin-1,
            "DC_admission":thin_admission},
        "minimum_thickness_for_requested_current_m":film*current/idep7,
        "same_q_validation":{
            "gap_bar_absolute_difference":abs(scaled_amplitude-reference.amplitude_bar),
            "current_scaling_relative_error":abs(
                scaled_current/(reference.current_A*thin/film)-1)},
        "mesh_cases_same_physical_film":mesh_rows,
        "interpretation":[
            "h labels numerical in-plane mesh spacing, not physical film thickness.",
            "At fixed material D,N0,Tc and width, gap(q) and J(q) are thickness-independent.",
            "Physical film thickness changes I=W*d_f*J and Rs=1/(sigma_n*d_f).",
            "Holding Rs fixed while changing d_f would also change conductivity and inferred N0 or D.",
            "The 30 percent thinner film crosses this model's stable DC branch limit at the same operating current.",
            "No transient, quench dynamics or experimental critical-current calibration is inferred."],
        "outputs":{
            "source":str(Path(__file__).resolve()),
            "csv":{"path":str(csvpath),"sha256":sha(csvpath)},
            "png":{"path":str(png),"sha256":sha(png)},
            "pdf":{"path":str(pdf),"sha256":sha(pdf)}},
        "elapsed_seconds":time.monotonic()-start,
        "resources":{"physical_workers":1,"cpu_affinity":[8],"BLAS_threads":1}
    }
    (root/"thickness_audit.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({
        "output_root":str(root),"Idep_7nm_A":idep7,"Idep_4p9nm_A":idepthin,
        "reference_gap_meV":reference.delta_J/e*1e3,
        "thin_current_excess_fraction":current/idepthin-1,
        "elapsed_seconds":result["elapsed_seconds"]},indent=2))

if __name__=="__main__":
    main()
