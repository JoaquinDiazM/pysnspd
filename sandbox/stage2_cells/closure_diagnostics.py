"""Lightweight KWT/Debye checks against independent SI and integral references.

No full transient is launched. The inherited KWT clock and the synthetic Debye
parameters are explicit and independent of the kinetic relaxation time.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
import pysnspd.experimental.cell_closures as module
from pysnspd.experimental.cell_closures import CellScales,KWTMobility,DebyePhonons,ConditionalPhononShape
from pysnspd.experimental.cell_validation import ElectronicCell
from pysnspd.experimental.energy_catalog import OccupationEnergyCatalog

HBAR=1.054571817e-34
KB=1.380649e-23


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def scaled(a,b): return float(np.max(np.abs(np.asarray(a)-b))/max(1.,float(np.max(np.abs(b)))))


def dimensional_reference(z,force,theta,scales):
    delta=scales.delta0_J*np.asarray(z)
    tmob=max(theta*scales.delta0_J/KB,scales.Tb_K)
    ratio=tmob/scales.Tc_K
    tau=1e-12/(ratio/.50+ratio**3/2.47)
    tau0=np.pi*HBAR/(8*KB*scales.Tc_K)
    A0=scales.N0_per_J_m3*np.sqrt((1+ratio)/2)
    R=np.sqrt(1+4*tau*tau*np.dot(delta,delta)/HBAR**2)
    temporal=A0*tau0/R*(np.eye(2)+4*tau*tau*np.outer(delta,delta)/HBAR**2)
    velocity=np.linalg.solve(temporal,-scales.N0_per_J_m3*scales.delta0_J*np.asarray(force)/2)
    heat=2*A0*tau0/R*(np.dot(velocity,velocity)+tau*tau*(2*np.dot(delta,velocity))**2/HBAR**2)
    return velocity*scales.t_ref_s/scales.delta0_J,heat/scales.power_density_scale_W_m3


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog",type=Path,default=ROOT/"docs/implementation/stage1_r2/catalogs/occupation_catalog.npz")
    parser.add_argument("--output-root",type=Path,default=ROOT/"docs/implementation/stage2")
    args=parser.parse_args(); started=time.perf_counter()
    out=args.output_root; (out/"figures").mkdir(parents=True,exist_ok=True)
    catalog=OccupationEnergyCatalog.load(args.catalog)
    delta0=catalog.vacuum.delta0_J; N0=catalog.vacuum.N0_per_J_m3
    scales=CellScales(delta0,N0,8.65,.12*delta0/KB,1.)
    closure=KWTMobility(scales)
    algebra=[]
    for amp in [0.,.25,.6,1.,1.5]:
        for theta in [0.,.12,.3,.8]:
            for angle in [0.,.7]:
                z=amp*np.array([np.cos(angle),np.sin(angle)])
                force=np.array([.8,-.35]); response=closure.tensor_response(z,force,theta)
                vr,qr=dimensional_reference(z,force,theta,scales)
                algebra.append({"amplitude":amp,"temperature_bar":theta,"angle":angle,
                    "velocity_SI_comparison_scaled":scaled(response.velocity,vr),
                    "heat_SI_comparison_scaled":scaled(response.heat,qr),
                    "energy_identity_scaled":scaled(response.heat,-np.dot(force,response.velocity)),
                    "minimum_mobility_eigenvalue":float(np.linalg.eigvalsh(response.mobility)[0]),
                    "heat":response.heat})
    r2cases=[]
    for amp in [.35,.6,.9,1.2]:
        for gamma in [0.,.2]:
            cell=ElectronicCell(catalog,amp,gamma)
            for label,p in [("FD_0.12",cell.fermi_dirac(.12)),("packet",.2*np.exp(-.5*((cell.energies-2.5)/.45)**2))]:
                theta=cell.equivalent_temperature(p); energy,force,conjugate=cell.moments(p)
                response=closure.amplitude_response(amp,force,theta)
                r2cases.append({"amplitude":amp,"gamma":gamma,"population":label,
                    "equivalent_temperature_bar":theta,"Tmob_K":response.Tmob_K,
                    "force_from_R2":force,"amplitude_velocity_per_tbar":response.velocity,
                    "Q_per_tbar":response.heat,"R":response.R,"taupsi_ps":response.taupsi_ps,
                    "energy_balance_residual":abs(force*response.velocity+response.heat)})
    reference=[]; debye=[]
    for low in [.02,.01,.005]:
        model=DebyePhonons(scales,4.,10*N0*delta0,.1,"stage2 synthetic Debye; n_atom=10*N0*Delta0",low)
        coefficient=90./4**3
        thermal_refs={}
        for theta in [.03,.12,.6,2.]:
            def u(e):
                if e==0:return 0.
                q=e/theta
                return coefficient*e**3*np.exp(-q)/(-np.expm1(-q))
            def c(e):
                if e==0:return 0.
                q=e/theta
                return coefficient*e**2*q*q*np.exp(-q)/(-np.expm1(-q))**2
            thermal_refs[theta]=(quad(u,low,4,epsabs=1e-14,epsrel=1e-12)[0],quad(c,low,4,epsabs=1e-14,epsrel=1e-12)[0])
            full=quad(u,0,4,epsabs=1e-14,epsrel=1e-12)[0]
            reference.append({"infrared_cutoff_bar":low,"temperature_bar":theta,
                "retained_U_reference":thermal_refs[theta][0],"retained_C_reference":thermal_refs[theta][1],
                "full_Debye_U_reference":full,"thermal_U_fraction_omitted_IR":(full-thermal_refs[theta][0])/full})
        for nodes in [17,33,65,129]:
            rule=model.quadrature(nodes)
            thermal=[]
            for theta,(u_ref,c_ref) in thermal_refs.items():
                u,c=rule.thermal_energy(theta),rule.thermal_capacity(theta)
                thermal.append({"temperature_bar":theta,"U":u,"C":c,
                    "U_relative_error":abs(u-u_ref)/u_ref,"C_relative_error":abs(c-c_ref)/c_ref})
            target_modes=30*(1-(low/4)**3)
            target_lambda=.1*(1-(low/4)**2)
            debye.append({"nodes":nodes,"infrared_cutoff_bar":low,"metadata":model.metadata(),
                "minimum_capacity":float(np.min(rule.capacities)),
                "mode_integral":float(np.sum(rule.capacities)),
                "mode_relative_error":abs(float(np.sum(rule.capacities))-target_modes)/target_modes,
                "lambda_integral":float(2*np.dot(rule.weights,model.alpha2F(rule.energies)/rule.energies)),
                "lambda_relative_error":abs(float(2*np.dot(rule.weights,model.alpha2F(rule.energies)/rule.energies))-target_lambda)/target_lambda,
                "thermal":thermal})
    conditional=ConditionalPhononShape(
        ROOT/"docs/implementation/stage1_closure/material_nbn_shape_v1.csv",
        ROOT/"docs/implementation/stage1_closure/material_nbn_shape_v1.manifest.json",
        source_axis_to_energy_bar=1.,dos_ordinate_to_dos_bar=1.,
        synthetic_label="identity-coordinate loader check only; no NbN dimensional interpretation",infrared_cutoff_bar=.01)
    maxima={
        "KWT_velocity_SI_scaled":max(r["velocity_SI_comparison_scaled"] for r in algebra),
        "KWT_heat_SI_scaled":max(r["heat_SI_comparison_scaled"] for r in algebra),
        "KWT_energy_identity_scaled":max(r["energy_identity_scaled"] for r in algebra),
        "KWT_minimum_mobility_eigenvalue":min(r["minimum_mobility_eigenvalue"] for r in algebra),
        "R2_Q_plus_force_times_velocity_absolute":max(r["energy_balance_residual"] for r in r2cases),
        "Debye_mode_relative":max(r["mode_relative_error"] for r in debye),
        "Debye_lambda_relative":max(r["lambda_relative_error"] for r in debye),
        "Debye_129_nodes_thermal_relative":max(v[key] for r in debye if r["nodes"]==129 for v in r["thermal"] for key in ["U_relative_error","C_relative_error"]),
    }
    gates={name:(value>0 if "minimum_mobility" in name else value<=(1e-8 if "thermal" in name else 1e-11)) for name,value in maxima.items()}
    plot(out,closure,reference,debye)
    criteria=ROOT/"docs/implementation/stage2/acceptance_criteria.json"
    result={"schema":"pysnspd.stage2.cell_closures.v1","status":"PASS" if all(gates.values()) else "FAIL",
        "scope":"Inherited KWT D.11-13 and synthetic Debye algebra/quadrature only; no NbN phonon SI admission, calibrated latency, spatial stability or detector transient.",
        "scales":scales.metadata(),"kinetic_relaxation_time":"separate driver input; not taupsi",
        "catalog_sha256":sha(args.catalog),"module_sha256":sha(module.__file__),"runner_sha256":sha(__file__),
        "criteria_sha256":sha(criteria) if criteria.exists() else None,"maxima":maxima,"gates":gates,
        "KWT_algebra_cases":algebra,"KWT_with_actual_R2_forces":r2cases,
        "Debye_quadrature":debye,"Debye_IR_reference":reference,
        "conditional_NbN_loader":conditional.metadata(),
        "runtime_seconds":time.perf_counter()-started,"python":platform.python_version(),"platform":platform.platform(),"numpy":np.__version__}
    (out/"closure_results.json").write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({"status":result["status"],"runtime_seconds":result["runtime_seconds"],"maxima":maxima},indent=2))
    if not all(gates.values()): raise SystemExit(1)


def plot(out,closure,infrared,debye):
    plt.rcParams.update({"font.size":11,"axes.titlesize":12,"axes.labelsize":11,"legend.fontsize":11})
    fig,axes=plt.subplots(1,2,figsize=(10,4.6),layout="constrained")
    amps=np.linspace(0,1.5,160)
    for theta in [.12,.3,.8]:
        mobility=[closure.amplitude_response(a,1,theta).mobility for a in amps]
        axes[0].plot(amps,mobility,label=f"kBT/Δ₀ = {theta:g}")
    axes[0].set(xlabel="Amplitud |Δ|/Δ₀",ylabel="Movilidad radial [1/t̄]",title="a) Movilidad KWT heredada")
    axes[0].legend(loc="upper right")
    for theta in [.03,.12,.6]:
        rows=[r for r in infrared if r["temperature_bar"]==theta]
        axes[1].loglog([r["infrared_cutoff_bar"] for r in rows],[100*r["thermal_U_fraction_omitted_IR"] for r in rows],"o-",label=f"kBT/Δ₀ = {theta:g}")
    axes[1].set(xlabel="Corte infrarrojo Ωmín/Δ₀",ylabel="Energía térmica omitida [%]",title="b) Debye sintético: alcance del corte")
    axes[1].legend(loc="lower right")
    for ax in axes: ax.grid(alpha=.2)
    fig.savefig(out/"figures/closure_mobility_phonons.png",dpi=180)
    fig.savefig(out/"figures/closure_mobility_phonons.pdf")
    plt.close(fig)


if __name__=="__main__": main()
