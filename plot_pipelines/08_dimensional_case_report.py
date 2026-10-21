"""Single-case v1.2 report from committed saved states; never evolves a model."""
from __future__ import annotations
import argparse, hashlib, importlib.util, json, os, textwrap, time
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
plt.rcParams.update({"font.size":12,"axes.labelsize":13,"axes.titlesize":13,
    "xtick.labelsize":11,"ytick.labelsize":11,"legend.fontsize":10})
INK="#17354a"
def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda:stream.read(2**20),b""):h.update(block)
    return h.hexdigest()
def scratch(path):
    p=Path(path).expanduser().resolve()
    if not p.is_relative_to(Path("/home/jdiaz/scratch")):raise ValueError("Products must stay in scratch")
    return p
def slots(values):
    a=[float(x) for x in values][:6]
    if not a or any(not np.isfinite(x) or x<0 for x in a):raise ValueError("Finite nonnegative snapshot times required")
    return (a+[a[-1]]*6)[:6]
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
class Report:
    def __init__(self,pdf,m,run):
        self.pdf,self.m,self.run=pdf,m,run
        self.pages=[]
    def page(self,title,rows=2,cols=2):
        fig,axes=plt.subplots(rows,cols,figsize=(11.69,8.27),squeeze=False)
        fig.suptitle(title,fontsize=17,x=.06,ha="left",y=.965,color=INK)
        fig.text(.06,.925,f"{self.run.name} | {self.m['parameter_identity']['wavelength_nm']:g} nm | v{self.m['model_version']}",fontsize=10,color=INK)
        return fig,axes.ravel()
    def save(self,fig,title,note=""):
        for a in fig.axes:
            label=a.get_ylabel()
            if len(label)>29 and "\n" not in label:
                a.set_ylabel("\n".join(textwrap.wrap(label,23)))
        fig.subplots_adjust(left=.09,right=.95,bottom=.17,top=.86,hspace=.48,wspace=.42)
        if note:fig.text(.06,.075,"\n".join(textwrap.wrap(note,145)),fontsize=9,va="top")
        fig.text(.94,.025,str(len(self.pages)+1),ha="right",fontsize=10)
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
            ax.grid(alpha=.2);self.events(ax)
            handles,_=ax.get_legend_handles_labels()
            if handles:ax.legend()
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
    fig,ax=r.page("1. Case inputs and terminal status",1,2)
    for a in ax:a.axis("off")
    left=[("Status",m.get("status")),("Official phase / last time",f"{m.get('phase')} / {number(m.get('time_ps'))} ps"),
        ("Bias / wavelength",f"{p['current_uA']:g} uA / {p['wavelength_nm']:g} nm"),
        ("Retained fraction / energy",f"{p['retained_fraction']:g} / {number(s.get('retained_energy_J'),1/e)} eV"),
        ("Deposition sigma",f"{p['sigma_nm']:g} nm"),("Width / thickness",f"80 nm / {v['film_thickness_m']*1e9:g} nm"),
        ("Lengths 2D / 1D",f"{p['length_nm']:g} / {p.get('length_1d_nm',p['length_nm']):g} nm"),
        ("Mesh spacing 2D / 1D",f"{p['spacing_nm']:g} / {p['spacing_1d_nm']:g} nm"),
        ("Tc / bath",f"{v['Tc_K']:g} / {v['bath_temperature_K']:g} K"),
        ("D / sheet resistance",f"{v['diffusion_m2_s']*1e4:.5g} cm2/s / {p['sheet_resistance_ohm']:g} ohm/sq"),
        ("tau_ee(Tc) / tau_ep(Tc)",f"{p['tau_ee_Tc_ps']:g} / {p['tau_ep_Tc_ps']:g} ps"),
        ("Phonon escape",f"{p['tau_escape_ps']:g} ps")]
    right=[("Voltage trigger",f"+{p['threshold_mV']:g} mV"),("Crossing / confirmation",f"{number(obs.get('crossing_time_s'),1e12)} / {number(obs.get('confirmation_time_s'),1e12)} ps"),
        ("Official handoff",f"{number(obs.get('stop_target_s'),1e12)} ps"),
        ("No-trigger ceiling / final target",f"{p['ceiling_ps']:g} ps / {p['final_time_ns']:g} ns"),
        ("Independent energy error",f"{number(m.get('energy_defect_J'),1/e)} eV"),
        ("Maximum global / local defect",f"{number(m.get('max_global_energy_fraction'),100)} / {number(m.get('max_local_energy_fraction'),100)} %"),
        ("Latest cut/activity fraction",f"{number(m.get('domain_fraction'),100)} %"),
        ("Readout model","Intrinsic Rload (Ib - Is)"),
        ("Cold approximation active",str(bool((m.get('cold_recovery') or {}).get('active')))),
        ("1D session wall time",f"{number(m.get('runtime_1d_session_wall_s'),1/60)} min"),
        ("Accepted steps / rejections",f"{m.get('accepted_steps')} / {m.get('rejections')}")]
    for a,rows in zip(ax,(left,right)):
        for i,(key,value) in enumerate(rows):
            y=.99-i*.077;a.text(0,y,key,weight="bold",fontsize=10,va="top")
            a.text(0,y-.029,str(value),fontsize=11,va="top")
    r.save(fig,"Inputs and terminal status","Stop reason: "+str(m.get("monitor_reason") or m.get("error") or m.get("status"))+
        ". Failed/censored cases remain diagnostic evidence. No amplifier, filter, bandwidth or readout noise is applied.")
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

def maps_2d(r,q,data,selected,requested,flow,legacy):
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
        fig,axes=r.page(title,3,2)
        values={}
        for j in set(selected):
            a=data[j]
            values[j]=(abs(a["d"])*Eu/meV if kind=="gap" else a[kind] if kind in ("Te","Tph")
                else a["u"]*1e-3 if kind=="u" else np.linalg.norm(flow[j],axis=1)/1e10 if kind=="total" else vortex[j]["winding"])
        lo=min(float(np.min(a)) for a in values.values());hi=max(float(np.max(a)) for a in values.values())
        if kind in ("gap","total"):lo=0.
        if kind=="u":hi=max(abs(lo),abs(hi),1e-15);lo=-hi
        if hi<=lo:hi=lo+max(abs(lo)*1e-6,1e-12)
        for ax,j,want in zip(axes,selected,requested):
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
            ax.set(title=f"Requested {want:.4g} / saved {saved[j]:.5g} ps",xlabel="x - xdep (nm)",ylabel="y (nm)",aspect="equal")
            ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5));ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3))
        fig.subplots_adjust(left=.09,right=.93,bottom=.225,top=.86,hspace=.45,wspace=.3)
        cb=fig.colorbar(image,cax=fig.add_axes([.22,.13,.58,.018]),orientation="horizontal",ticks=[-1,0,1] if kind=="winding" else None)
        cb.set_label(label,fontsize=12)
        note=("Complex P1 gap magnitude, not interpolation of nodal magnitudes." if kind=="gap" else
            "Energy-equivalent moment of native nonthermal occupations; not an imposed thermodynamic temperature." if kind in ("Te","Tph") else
            "Native BCS internal energy, thermal gradient and phonons minus the phase DC reference." if kind=="u" else
            "Algebraic reconstruction at the exact saved state. Dual-flux least-squares display vectors; arrows show direction with fixed length. No separate normal/supercurrent snapshot pages." if kind=="total" else
            "Phase topology, not curl(J). Gray cells have an unresolved or vanishing-gap contour; no vortex trajectories are inferred between snapshots.")
        fig.text(.06,.067,"\n".join(textwrap.wrap(note+" Exactly six nearest saved states; repeated frames are retained. Shared color range per quantity.",145)),fontsize=9,va="top")
        fig.text(.94,.025,str(len(r.pages)+1),ha="right",fontsize=10)
        r.pdf.savefig(fig,dpi=150);plt.close(fig);r.pages.append(title)
    return [dict(snapshot_index=j,time_ps=float(saved[j]),positive=vortex[j]["positive"],negative=vortex[j]["negative"],
        indeterminate_cells=vortex[j]["indeterminate_cells"]) for j in sorted(set(selected))]
def longitudinal(r,q,data,selected,requested):
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
    fig,axes=r.page("12. Longitudinal profiles through the terminal dynamics")
    colors=plt.cm.viridis(np.linspace(.08,.9,6))
    for j,color in zip(selected,colors):
        a=data[j];label=f"{tt[j]:.5g} ps"
        axes[0].plot(x,abs(a["d"])*Eu/(e*1e-3),c=color,label=label)
        axes[1].plot(x,a["Te"],c=color,label=label)
        axes[2].plot(x,a["Tph"],c=color,label=label)
        axes[3].plot(x,a["u"]*1e-3,c=color,label=label)
    for ax,y in zip(axes,("Gap magnitude (meV)","Te (K)","Tph (K)","Signed excess energy density (kJ/m3)")):
        ax.set(xlabel="x - strip center (nm)",ylabel=y);ax.grid(alpha=.2)
        for attr in ("core_begin_m","core_end_m"):
            if attr in q.attrs:ax.axvline(q.attrs[attr]*1e9-origin,color=".5",ls=":")
        ax.axvspan(x[0],x[0]+cfg.get("collar_nm",25.),color="red",alpha=.07)
        ax.axvspan(x[-1]-cfg.get("collar_nm",25.),x[-1],color="red",alpha=.07)
    axes[0].legend(ncol=2,fontsize=9)
    r.save(fig,"Longitudinal profiles","Six nearest official 1D states, selected independently of the 2D mosaic. Shaded regions include the terminal cut, whereas the guard excludes the fixed terminal node. The added tails and original lens share one physical x axis; no 2D state is tiled into the tails.")
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
    fig.text(.09,.502,f"Rb={circ['R_bias_ohm']:g} ohm; RL={circ['R_load_ohm']:g} ohm; Lb={circ['L_bias_H']*1e6:g} uH; Cc={circ['C_couple_F']*1e12:g} pF; Lext(2D/1D)={circ['Lk_ext_H']*1e9:.5g}/{ext1*1e9:.5g} nH",fontsize=10)
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
    if m.get("model_version") not in ("1.2.0","1.2.1","1.2.2","1.2.3"):raise ValueError("A supported v1.2 phase-aware archive is required")
    if m.get("status") in ("RUNNING","PREPARATION"):raise ValueError("Use a terminal committed case, not a live archive")
    if digest(run/"checkpoint_latest.npz")!=m.get("checkpoint_sha256"):raise ValueError("Committed checkpoint checksum mismatch")
    target=out/("case_report_"+run.name+".pdf");receipt=target.with_suffix(".json")
    if target.exists() or receipt.exists():raise FileExistsError("Existing report preserved; choose a fresh output directory")
    build=out/".build"/run.name;build.mkdir(parents=True,exist_ok=False)
    temporary=build/"case_report.pdf";pdf=PdfPages(temporary,metadata={"Title":run.name+" / v"+m["model_version"],"Author":"pysnspd"})
    result=dict(schema="pysnspd.dimensional_case_report.v2",run_dir=str(run),model_version=m["model_version"],
        status=m.get("status"),stop_reason=m.get("monitor_reason") or m.get("error"),model_simulations_launched=0,
        kinetic_rhs_evaluations=0,new_dc_solves=0,temporal_interpolation=False,
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
            if counts["2d"]["snapshot_time_ps"]:
                q=f["2d"];saved=q["snapshot_time_ps"][:counts["2d"]["snapshot_time_ps"]]
                crossing=m.get("observation",{}).get("crossing_time_s")
                crossing=float(crossing)*1e12 if crossing is not None else float(saved[-1])
                wants=slots(requested2 if requested2 is not None else
                    [0.,max(0,crossing-10),max(0,crossing-4),max(0,crossing-1),crossing,float(saved[-1])])
                selected=nearest(saved,wants);data=snapshot_data(q,selected,m,grid)
                flow,checks,mesh=total_currents(run,q,selected,data,m,traces(q,m,"2d"),legacy)
                result["vorticity"]=maps_2d(r,q,data,selected,wants,flow,legacy)
                result["current_reconstruction_checks"]=checks;result["mesh_sha256"]=digest(mesh)
                result["snapshots_2d"]=[dict(requested_ps=float(w),saved_ps=float(saved[j]),snapshot_index=j) for w,j in zip(wants,selected)]
            if counts["1d"]["snapshot_time_ps"]:
                q=f["1d"];n=counts["1d"]["snapshot_time_ps"];saved=q["snapshot_time_ps"][:n]
                wants=slots(requested1 if requested1 is not None else np.linspace(saved[0],saved[-1],6))
                selected=nearest(saved,wants);data=snapshot_data(q,range(n),m,grid)
                result["collar_snapshots"]=longitudinal(r,q,data,selected,wants)
                spectra(r,q,data,grid)
                result["snapshots_1d"]=[dict(requested_ps=float(w),saved_ps=float(saved[j]),snapshot_index=j) for w,j in zip(wants,selected)]
            overlap(r);summary(r,legacy)
            result["archive_counts"]=counts;result["material_sha256"]=digest(run/"setup2d/material.npz")
        if manifest_path.read_bytes()!=raw:raise RuntimeError("Case changed during report; output not promoted")
    finally:pdf.close()
    result.update(pdf=str(target),pages=r.pages,pdf_sha256=digest(temporary),elapsed_wall_s=time.monotonic()-started,
        caveats=["Collar temporal curves are reconstructed only at saved snapshots; accepted-step persistence is not recreated",
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
    p.add_argument("--snapshots-ps",nargs="+",type=float,help="Six 2D slots; nearest saved states; default event-centered")
    p.add_argument("--snapshots-1d-ps",nargs="+",type=float,help="Six longitudinal-profile slots; default handoff through terminal state")
    p.add_argument("--workers",type=int,choices=(1,2),default=1)
    a=p.parse_args(argv);out=scratch(a.output_dir);out.mkdir(parents=True,exist_ok=True)
    for vals in (a.snapshots_ps,a.snapshots_1d_ps):
        if vals is not None:slots(vals)
    if len({x.resolve().name for x in a.run_dir})!=len(a.run_dir):raise ValueError("Duplicate report names; use separate destinations")
    jobs=[(str(x),str(out),a.snapshots_ps,a.snapshots_1d_ps) for x in a.run_dir]
    print(json.dumps(dict(status="SAVED_DATA_ONLY",cases=len(jobs),workers=min(a.workers,len(jobs)),
        numba_threads_per_worker=1,blas_threads_per_worker=1,simulations=0)),flush=True)
    if a.workers==1:
        for j in jobs:report(j)
    else:
        with ProcessPoolExecutor(max_workers=min(a.workers,len(jobs))) as pool:list(pool.map(report,jobs))
if __name__=="__main__":main()
