"""Phase-aware individual-case report for v1.2 HDF5 trajectories, in English."""
from __future__ import annotations
import argparse,json,os
from pathlib import Path
for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):os.environ[key]="1"
os.environ.setdefault("MPLCONFIGDIR","/home/jdiaz/scratch/tmp/pysnspd/matplotlib")
import numpy as np,h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from concurrent.futures import ProcessPoolExecutor
plt.rcParams.update({"font.size":12,"axes.labelsize":13,"axes.titlesize":13,
                     "xtick.labelsize":11,"ytick.labelsize":11,"legend.fontsize":11})


def report(job):
    run,out,times=job
    run=Path(run);out=Path(out);m=json.loads((run/"manifest.json").read_text())
    if m.get("model_version") not in ("1.2.0","1.2.1","1.2.2"):raise ValueError("This report requires the v1.2 phase-aware archive")
    target=out/("case_report_"+run.name+".pdf")
    if target.exists():raise FileExistsError("Existing report preserved: "+str(target))
    counts=m["archive_counts"];p=m["parameter_identity"]
    with h5py.File(run/"data.h5","r",swmr=True) as f,PdfPages(target) as pdf:
        def done(fig):
            fig.tight_layout();pdf.savefig(fig);plt.close(fig)
        fig=plt.figure(figsize=(11.69,8.27))
        fig.text(.06,.92,"Individual case / neutral effective BCS-KWT / v"+m["model_version"],fontsize=18)
        lines=[f"Case: {run.name}",f"Status: {m['status']} / official phase: {m['phase']}",
            f"Bias {p['current_uA']:g} uA; wavelength {p['wavelength_nm']:g} nm; retained fraction {p['retained_fraction']:g}",
            f"Escape {p['tau_escape_ps']:g} ps; tau_ee(Tc) {p['tau_ee_Tc_ps']:g} ps; Rs {p['sheet_resistance_ohm']:g} ohm",
            f"2D: {p['length_nm']:g} x 80 nm; 1D: {p.get('length_1d_nm',p['length_nm']):g} nm; spacing {p['spacing_1d_nm']:g} nm",
            f"Threshold {p['threshold_mV']:g} mV; confirmation {p['confirmation_ps']:g} ps; extra {p['extra_time_ps']:g} ps",
            f"No-trigger ceiling {p['ceiling_ps']:g} ps; final target {p['final_time_ns']:g} ns",
            "Vout = Rload (Ib - Is): intrinsic circuit signal, without amplifier or filter.",
            f"Cut-energy / device-activity fraction: {m['domain_fraction']:.5g}",
            f"Independent energy defect: {m['energy_defect_J']:.5g} J",
            "A completed trajectory does not by itself establish physical or statistical validation."]
        cold=m.get("cold_recovery") or {}
        if cold.get("active"):
            lines.append(f"Late cold approximation from {cold['activation_time_ps']*1e-3:.5g} ns: spectra/maps reconstructed at the bath.")
            lines.append("The hot stage evolved native populations; late quasiparticle tails are outside this approximation.")
        for i,line in enumerate(lines):fig.text(.06,.83-.057*i,line,fontsize=12)
        pdf.savefig(fig);plt.close(fig)
        official=[]
        for name in ("2d","1d"):
            n=counts[name]["trace"]
            if n:official.append((name,f[name]["trace"][:n]))
        for full in (False,True):
            fig,axes=plt.subplots(2,2,figsize=(11.69,8.27))
            unit=1e-3 if full else 1.
            xlabel="Time after photon (ns)" if full else "Time after photon (ps)"
            for name,tr in official:
                if not full and name=="1d":continue
                x=tr["time_ps"]*unit
                axes[0,0].plot(x,tr["Vout_V"]*1e3,label=name)
                axes[0,1].plot(x,tr["Is_A"]*1e6,label="Is / "+name)
                axes[0,1].plot(x,(tr["Ib_A"]-tr["Is_A"])*1e6,ls="--",label="Iread / "+name)
                axes[1,0].plot(x,tr["energy_defect_J"]*1e18,label=name)
                axes[1,1].plot(x,tr["domain_fraction"],label=name)
            for ax,label in zip(axes.flat,("Vout (mV), unfiltered","Branch currents (uA)",
                "Independent energy defect (aJ)","Absolute cut energy / activity")):
                ax.set(xlabel=xlabel,ylabel=label);ax.grid(alpha=.25);ax.legend()
            axes[0,0].axhline(p["threshold_mV"],color=".5",ls=":")
            axes[1,1].axhline(p["guard_cut_energy_relative"],color="r",ls=":")
            done(fig)
        group=f["2d"];n=counts["2d"]["snapshot_time_ps"]
        if n:
            tt=group["snapshot_time_ps"][:n];xy=group["coordinates_m"][:]*1e9;tri=group["triangles"][:]
            dc=np.mean(abs(group["dc_delta_bar"][:]))
            fig,axes=plt.subplots(3,2,figsize=(11.69,8.27),sharex=True,sharey=True)
            for ax,want in zip(axes.flat,times):
                j=int(np.argmin(abs(tt-want)));z=abs(group["delta_bar"][j])/dc
                art=ax.tripcolor(xy[:,0],xy[:,1],tri,z,shading="gouraud",vmin=0,vmax=1,cmap="viridis")
                ax.set(title=f"Requested {want:g} ps / stored {tt[j]:.4g} ps",xlabel="x (nm)",ylabel="y (nm)",aspect="equal")
                fig.colorbar(art,ax=ax,label="|Delta| / DC gap")
            done(fig)
        group=f["1d"];n=counts["1d"]["snapshot_time_ps"]
        if n:
            tt=group["snapshot_time_ps"][:n]*1e-3;x=group["coordinates_m"][:,0]*1e9
            for key,label,cmap in (("delta_bar","|Delta| / DC gap","viridis"),
                                   ("Te_K","Equivalent electron temperature (K)","inferno"),
                                   ("u_excess_J_m3","Excess energy density (J/m3)","cividis")):
                values=group[key][:n]
                if key=="delta_bar":values=abs(values)/np.mean(abs(group["dc_delta_bar"][:]))
                fig,ax=plt.subplots(figsize=(11.69,8.27))
                im=ax.pcolormesh(x,tt,values,shading="nearest",cmap=cmap)
                ax.set(xlabel="x (nm)",ylabel="Time after photon (ns)",title="Official 1D evolution")
                if cold.get("active"):
                    ax.axhline(cold["activation_time_ps"]*1e-3,color="white",ls="--",lw=1.3)
                    ax.set_title("1D evolution; dashed line starts the late cold reconstruction")
                fig.colorbar(im,ax=ax,label=label);done(fig)
        checks=m.get("shadow_checks",[])
        if checks:
            fig,axes=plt.subplots(2,2,figsize=(11.69,8.27))
            for ax,key,label in zip(axes.flat,("gap_relative","populations_energy_relative","voltage_error_V","film_energy_relative"),
                ("Gap reconstruction RMS / DC","Unsigned population-energy discrepancy / activity",
                 "Independent Vout discrepancy (uV)","Film-energy discrepancy / activity")):
                ax.plot([c["time_ps"] for c in checks],[c[key]*(1e6 if key=="voltage_error_V" else 1) for c in checks])
                ax.set(xlabel="Time after photon (ps)",ylabel=label);ax.grid(alpha=.25)
            done(fig)
    return str(target)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run-dir",nargs="+",type=Path,required=True)
    p.add_argument("--output-dir",type=Path,required=True)
    p.add_argument("--snapshots-ps",nargs="+",type=float,default=[1,11,11.75,12.5,20,35])
    p.add_argument("--workers",type=int,choices=(1,2),default=1)
    a=p.parse_args();a.output_dir.mkdir(parents=True,exist_ok=True)
    times=(a.snapshots_ps+[a.snapshots_ps[-1]]*6)[:6]
    if any(not np.isfinite(t) or t<0 for t in times):raise ValueError("Nonnegative finite snapshot times required")
    jobs=[(str(run),str(a.output_dir),times) for run in a.run_dir]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for result in pool.map(report,jobs):print(result,flush=True)

if __name__=="__main__":main()
