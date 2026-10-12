"""Read-only physics diagnosis of coherent active-run checkpoints; no evolution."""
from __future__ import annotations
import argparse,csv,hashlib,io,json,os
from pathlib import Path
import sys
REPO=Path(__file__).resolve().parents[2]
SCRATCH=Path("/home/jdiaz/scratch")
sys.path.insert(0,str(REPO))
sys.pycache_prefix=str(SCRATCH/"pycache")
for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="1"
import numpy as np
from scipy.special import expit
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from pysnspd.experimental.effective_bcs_functional import bcs_count_quadrature
from pysnspd.experimental.effective_bcs_kinetics import acoustic_phonon_grid
from pysnspd.experimental.effective_material import EffectivePhononMaterial
EV=1.602176634e-19
CASES={
"775_dt002":"stage5_effective_L320_h5_775_dt002_v3_20261011",
"1550_dt002":"stage5_effective_L320_h5_1550_dt002_v3_20261011",
"775_dt001":"stage5_effective_L320_h5_775_dt001_v3_20261011",
"775_h35":"stage5_effective_L320_h35_775_dt002_v3_20261011"}

def cumulative_trapezoid(t,p):
    return np.r_[0,np.cumsum(.5*(p[1:]+p[:-1])*np.diff(t))]

def component_metric(mask,ratio,edges,xy,area):
    selected=edges[mask[edges].all(axis=1)]
    if not np.any(mask):
        return {"area_nm2":0.,"x_extent_nm":0.,"y_extent_nm":0.,"width_spanning":False}
    adj=coo_matrix((np.ones(2*len(selected)),
        (np.r_[selected[:,0],selected[:,1]],np.r_[selected[:,1],selected[:,0]])),
        shape=(len(xy),len(xy))).tocsr()
    _,labels=connected_components(adj,directed=False)
    dominant=int(np.argmin(np.where(mask,ratio,np.inf)))
    member=mask&(labels==labels[dominant])
    yy=xy[:,1];tol=max(np.ptp(yy)*1e-10,1e-9)
    low=mask&(abs(yy-yy.min())<=tol);high=mask&(abs(yy-yy.max())<=tol)
    spanning=bool(set(labels[low])&set(labels[high]))
    return {"area_nm2":float(area[mask].sum()*1e18),
        "x_extent_nm":float(np.ptp(xy[member,0])),
        "y_extent_nm":float(np.ptp(xy[member,1])),
        "width_spanning":spanning}

def capture(name,path):
    raw=(path/"checkpoint_latest.npz").read_bytes()
    with np.load(io.BytesIO(raw),allow_pickle=False) as z:
        meta=json.loads(str(z["metadata_json"].item()))
        h={k.split("__",1)[1]:z[k].copy() for k in z.files if k.startswith("trajectory__")}
        s={k.split("__",1)[1]:z[k].copy() for k in z.files if k.startswith("snapshots__")}
        p,n=z["p"].copy(),z["n"].copy()
    manifest=json.loads((path/"manifest.json").read_text())
    physical=manifest["physicalparams"];identity=meta["parameter_identity"]
    with np.load(manifest["geometry"]["source_path"],allow_pickle=False) as z:
        xy=z["coordinates_m"].copy()*1e9
        edges=z["edges"].copy();area=z["area_weights"].copy()*float(z["ell0_m"])**2
    volume=area*physical["film_thickness_m"]
    ratio=abs(s["delta_J"])/abs(s["delta_J"][0])
    times=s["time_s"]*1e12
    target=meta["source"]["retained_energy_J"]
    Eu,N0,Tb=physical["energy_scale_J"],physical["N0_single_spin_per_J_m3"],identity["bath_K"]
    xi,w=bcs_count_quadrature(cutoff_bar=identity["xi_cutoff"],order=identity["xi_order"])
    e0=np.sqrt(xi**2+(abs(s["delta_J"][0])/Eu)[:,None]**2)
    cold_e=4*N0*Eu**2*np.sum(w*e0*expit(-e0*Eu/(1.380649e-23*Tb)),axis=1)
    material=EffectivePhononMaterial.load(path/"material.npz")
    grid=acoustic_phonon_grid(material,Eu,nodes=identity["phonon_nodes"])
    thermal_x=grid.energies*Eu/(1.380649e-23*Tb)
    nb=np.zeros_like(thermal_x);keep=thermal_x<700;nb[keep]=1/np.expm1(thermal_x[keep])
    cold_ph=N0*Eu**2*np.sum(grid.capacities*grid.energies*nb)
    rows=[];center=np.array(meta["source"]["center_m"])*1e9
    icenter=int(np.argmin(np.sum((xy-center)**2,axis=1)))
    for index,t in enumerate(times):
        excess=s["u_excess_J_m3"][index]*volume
        positive=np.maximum(excess,0);total_positive=float(positive.sum())
        if total_positive>0:
            centroid=np.sum(positive[:,None]*xy,axis=0)/total_positive
            rms=np.sqrt(np.sum(positive[:,None]*(xy-centroid)**2,axis=0)/total_positive)
        else:centroid=np.array([np.nan,np.nan]);rms=np.array([np.nan,np.nan])
        electron=float(np.sum((s["electron_energy_J_m3"][index]-cold_e)*volume))
        phonon=float(np.sum((s["phonon_energy_J_m3"][index]-cold_ph)*volume))
        total=float(excess.sum())
        cutzone=(xy[:,0]<xy[:,0].min()+20)|(xy[:,0]>xy[:,0].max()-20)
        row={"time_ps":float(t),"gap_min_meV":float(abs(s["delta_J"][index]).min()/EV*1e3),
            "gap_center_meV":float(abs(s["delta_J"][index,icenter])/EV*1e3),
            "gap_min_fraction_DC":float(ratio[index].min()),
            "gap_center_fraction_DC":float(ratio[index,icenter]),
            "Te_max_K":float(s["Te_K"][index].max()),"Tph_max_K":float(s["Tph_K"][index].max()),
            "film_excess_eV":total/EV,"electron_excitation_excess_eV":electron/EV,
            "phonon_excess_eV":phonon/EV,"condensate_gradient_excess_eV":(total-electron-phonon)/EV,
            "film_remaining_fraction_retained":total/target,
            "positive_excess_rms_x_nm":float(rms[0]),"positive_excess_rms_y_nm":float(rms[1]),
            "negative_energy_inventory_eV":float(np.minimum(excess,0).sum()/EV),
            "positive_energy_fraction_near_20nm_cuts":float(positive[cutzone].sum()/total_positive)}
        for fraction in (.95,.8,.5,.2):
            metric=component_metric(ratio[index]<=fraction,ratio[index],edges,xy,area)
            row.update({"gap_"+str(fraction)+"_"+key:value for key,value in metric.items()})
        rows.append(row)
    ht=h["time_s"];esc=cumulative_trapezoid(ht,h["escape_power_W"])
    cuts=cumulative_trapezoid(ht,h["reservoir_power_out_W"])
    port=cumulative_trapezoid(ht,h["Is_A"]*h["Vdev_V"])
    totalwork=np.r_[abs(h["external_energy_integral_J"][0]),
        np.cumsum(abs(np.diff(h["external_energy_integral_J"])))+abs(h["external_energy_integral_J"][0])]
    imin=int(np.argmin([r["gap_min_fraction_DC"] for r in rows]))
    tail_start=max(0,len(rows)-6)
    def slope(key):
        if len(rows)-tail_start<2:return None
        return float(np.polyfit(times[tail_start:],np.array([r[key] for r in rows[tail_start:]]),1)[0])
    summary={"name":name,"run_path":str(path),"checkpoint_capture_sha256":hashlib.sha256(raw).hexdigest(),
        "captured_time_ps":float(meta["time_ps"]),"accepted_steps":meta["accepted_steps"],
        "manifest_status_at_read":manifest["status"],"scientific_admission":False,
        "source_retained_energy_eV":target/EV,
        "snapshot_count":len(rows),"snapshot_end_ps":float(times[-1]),
        "initial":rows[0],"latest":rows[-1],
        "maximum_gap_suppression":{"time_ps":rows[imin]["time_ps"],"gap_fraction_DC":rows[imin]["gap_min_fraction_DC"],
            "gap_min_meV":rows[imin]["gap_min_meV"]},
        "latest_trend_window_ps":[float(times[tail_start]),float(times[-1])],
        "latest_slopes_per_ps":{key:slope(key) for key in
            ("gap_min_fraction_DC","Te_max_K","Tph_max_K","film_excess_eV")},
        "Vout_peak_mV":float(h["Vout_V"].max()*1e3),"Vout_latest_mV":float(h["Vout_V"][-1]*1e3),
        "Vdev_latest_mV":float(h["Vdev_V"][-1]*1e3),"Is_latest_uA":float(h["Is_A"][-1]*1e6),
        "Ib_latest_uA":float(h["Ib_A"][-1]*1e6),
        "latest_powers_nW":{"escape":float(h["escape_power_W"][-1]*1e9),
             "reservoir_out":float(h["reservoir_power_out_W"][-1]*1e9),
             "device_input":float(h["Is_A"][-1]*h["Vdev_V"][-1]*1e9)},
        "integrated_power_diagnostic_eV":{"escape":float(esc[-1]/EV),"reservoir_out":float(cuts[-1]/EV),
             "device_input":float(port[-1]/EV),
             "quadrature":"Trapezoid of saved accepted endpoint powers; not the exact midpoint energy certificate"},
        "energy_defect_maximum_fraction":float(np.max(abs(h["energy_defect_J"])/(target+totalwork))),
        "event_limited_energy_fraction":float(meta["limited_event_energy_J"]/max(target,1e-300)),
        "support_final":{"p_min":float(p.min()),"p_max":float(p.max()),"n_min":float(n.min())},
        "parameter_identity":identity}
    return summary,rows,h

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args();out=args.output.resolve()
    if not out.is_relative_to(SCRATCH):raise ValueError("Write only under scratch")
    if out.exists():raise FileExistsError("Preserve previous captures")
    os.sched_setaffinity(0,{13})
    captured={name:capture(name,SCRATCH/path) for name,path in CASES.items()}
    out.mkdir(parents=True)
    report={"status":"PROVISIONAL_PARTIAL_SAVED_DATA_NOT_STAGE5_ADMISSION",
        "simulations_launched":0,"CPU":[13],"runs":{name:row[0] for name,row in captured.items()}}
    (out/"partial_physics.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    for name,(_,rows,_) in captured.items():
        with (out/(name+"_snapshots.csv")).open("w",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(2,2,figsize=(12,8),layout="constrained")
    for name,(summary,rows,h) in captured.items():
        t=[r["time_ps"] for r in rows]
        axs[0,0].plot(t,[r["gap_min_fraction_DC"] for r in rows],label=name)
        axs[0,1].plot(t,[r["Te_max_K"] for r in rows],label=name+" electrons")
        axs[0,1].plot(t,[r["Tph_max_K"] for r in rows],ls="--",label=name+" phonons")
        axs[1,0].plot(t,[r["film_remaining_fraction_retained"] for r in rows],label=name)
        axs[1,1].plot(h["time_s"]*1e12,h["Vout_V"]*1e3,label=name)
    for ax,ylabel,title in zip(axs.ravel(),
        ["Gap minimum / fresh DC","Equivalent temperature (K)","Film energy excess / retained photon","Vout (mV)"],
        ["Condensate suppression","Energy equivalent temperatures","Physical remaining energy","Intrinsic electrical signal"]):
        ax.set(xlabel="Post-cascade time (ps)",ylabel=ylabel,title=title)
        ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle("PARTIAL SAVED TRAJECTORIES - no stage5 admission",fontsize=15)
    fig.savefig(out/"partial_physics.png",dpi=160);plt.close(fig)
    print(json.dumps({name:{"time_ps":s["captured_time_ps"],"gap_fraction":s["latest"]["gap_min_fraction_DC"],
        "Te_K":s["latest"]["Te_max_K"],"Tph_K":s["latest"]["Tph_max_K"],
        "energy_remaining":s["latest"]["film_remaining_fraction_retained"],
        "Vout_mV":s["Vout_latest_mV"]} for name,(s,_,_) in captured.items()},indent=2))
if __name__=="__main__":main()
