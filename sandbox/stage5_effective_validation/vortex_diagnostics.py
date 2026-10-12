"""Read saved effective-model arrays: phase topology, electrical impulse and heat.
No evolution, RHS, material reconstruction, admission or detection inference.
All numerical products remain in scratch. Winding assumes the saved gauge A=0.
"""
from __future__ import annotations
import argparse, hashlib, io, json, os
from collections import Counter
from pathlib import Path
for _key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[_key]="1"
import numpy as np

SCRATCH=Path("/home/jdiaz/scratch")
MEV_J=1.602176634e-22
EV_J=1.602176634e-19
PHI0=6.62607015e-34/(2*1.602176634e-19)

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def native(value):
    if isinstance(value, np.ndarray):return value.tolist()
    if isinstance(value, np.generic):return value.item()
    raise TypeError(type(value).__name__)

def increasing(t):
    t=np.asarray(t,float)
    if t.ndim!=1 or not len(t) or np.any(~np.isfinite(t)) or np.any(np.diff(t)<=0):
        raise ValueError("Finite increasing saved times required")
    return t

def ccw_triangles(xy,tri):
    tri=np.asarray(tri,int).copy()
    p=np.asarray(xy,float)[tri]
    u,v=p[:,1]-p[:,0],p[:,2]-p[:,0]
    area2=u[:,0]*v[:,1]-u[:,1]*v[:,0]
    if np.any(area2==0):raise ValueError("Degenerate spatial triangle")
    reverse=area2<0
    tri[reverse]=tri[reverse][:,[0,2,1]]
    return tri,int(reverse.sum())

def complex_triangle_minimum(delta,tri):
    """Exact minimum of |linear complex Delta| over every closed triangle."""
    z=np.asarray(delta,complex)[tri]
    a,b=z[:,1]-z[:,0],z[:,2]-z[:,0]
    cross=lambda u,v:u.real*v.imag-u.imag*v.real
    det=cross(a,b)
    scale=np.maximum(np.max(abs(z),axis=1)**2,np.finfo(float).tiny)
    valid=abs(det)>64*np.finfo(float).eps*scale
    l1=np.zeros(len(z));l2=np.zeros(len(z))
    l1[valid]=cross(-z[valid,0],b[valid])/det[valid]
    l2[valid]=cross(a[valid],-z[valid,0])/det[valid]
    bary=np.column_stack((1-l1-l2,l1,l2))
    inside=valid & np.all(bary>=-128*np.finfo(float).eps,axis=1)
    edge_min=[]
    for j,k in ((0,1),(1,2),(2,0)):
        d=z[:,k]-z[:,j];den=abs(d)**2
        t=np.zeros(len(z));ok=den>0
        t[ok]=np.clip(-np.real(np.conj(z[ok,j])*d[ok])/den[ok],0.,1.)
        edge_min.append(abs(z[:,j]+t*d))
    result=np.min(edge_min,axis=0)
    result[inside]=0.
    return result,inside,bary

def triangle_topology(delta,xy,tri):
    z=np.asarray(delta,complex)[tri]
    if np.any(abs(z)==0):raise ValueError("Nodal zero makes phase winding undefined")
    phases=np.angle(np.conj(z)*z[:,[1,2,0]])
    winding=phases.sum(axis=1)/(2*np.pi)
    charge=np.rint(winding).astype(int)
    if np.any(abs(winding-charge)>1e-9):raise ValueError("Noninteger discrete winding")
    minima,inside,bary=complex_triangle_minimum(delta,tri)
    rows=[]
    for index in np.flatnonzero(charge):
        polygon=np.asarray(xy)[tri[index]]
        distances=abs(polygon[:,None,:]-polygon[None,:,:])
        diameter=np.max(np.sqrt((distances**2).sum(axis=2)))
        center=polygon.mean(axis=0)
        linear_zero=(bary[index]@polygon) if inside[index] else None
        rows.append(dict(triangle_index=int(index),charge=int(charge[index]),
            centroid_nm=center*1e9,linear_complex_zero_nm=None if linear_zero is None else linear_zero*1e9,
            triangle_diameter_nm=float(diameter*1e9),
            minimum_nodal_gap_meV=float(abs(z[index]).min()/MEV_J),
            linear_complex_minimum_meV=float(minima[index]/MEV_J),
            nearest_principal_branch_margin_rad=float(np.pi-abs(phases[index]).max())))
    return charge,minima,rows,float(abs(winding-charge).max(initial=0.))

def central_linear_weights(xy,tri,center_m,halfwidth_x_m=60e-9):
    """Exact area weights of linear nodal fields after clipping to central window."""
    xy=np.asarray(xy,float)
    lo,hi=float(center_m[0]-halfwidth_x_m),float(center_m[0]+halfwidth_x_m)
    weights=np.zeros(len(xy))
    def clip(poly,bound,keep_above):
        result=[]
        for a,b in zip(poly,poly[1:]+poly[:1]):
            fa=a[0][0]-bound;fb=b[0][0]-bound
            ia=(fa>=0) if keep_above else (fa<=0)
            ib=(fb>=0) if keep_above else (fb<=0)
            if ia:result.append(a)
            if ia!=ib:
                q=fa/(fa-fb)
                result.append((a[0]+q*(b[0]-a[0]),a[1]+q*(b[1]-a[1])))
        return result
    for nodes in tri:
        poly=[(xy[nodes[k]],np.eye(3)[k]) for k in range(3)]
        poly=clip(clip(poly,lo,True),hi,False)
        if len(poly)<3:continue
        for j in range(1,len(poly)-1):
            vertices=(poly[0],poly[j],poly[j+1])
            u,v=vertices[1][0]-vertices[0][0],vertices[2][0]-vertices[0][0]
            area=abs(u[0]*v[1]-u[1]*v[0])/2
            bary=sum((v[1] for v in vertices))/3
            weights[nodes]+=area*bary
    return weights

def cumulative_trapezoid(y,t):
    return np.r_[0.,np.cumsum((np.asarray(y)[1:]+np.asarray(y)[:-1])/2*np.diff(t))]

def phase_unwrap_diagnostic(wrapped_relative_phase_rad, CM_impulse_V_s):
    """Only certify ordinary unwrapping when CM bounds each relative increment.
    CM remains valid if snapshots undersample a phase rotation. No extra 2pi
    turns are inferred or repaired from CM; ambiguity stays explicitly flagged.
    """
    wrapped=np.asarray(wrapped_relative_phase_rad,float)
    cm=np.asarray(CM_impulse_V_s,float)
    if wrapped.ndim!=1 or cm.shape!=wrapped.shape or np.any(~np.isfinite(wrapped)) or np.any(~np.isfinite(cm)):
        raise ValueError("Matching finite saved relative phase and CM impulse required")
    physical_increment=np.diff(cm)*2*np.pi/PHI0
    margin=1e-8
    ambiguous=np.flatnonzero(abs(physical_increment)>=np.pi-margin)
    unwrapped=np.unwrap(wrapped)
    info=dict(certified=not len(ambiguous),
        criterion="Absolute CM relative-phase increment strictly below pi-1e-8rad between consecutive saved snapshots",
        maximum_CM_phase_increment_rad=float(abs(physical_increment).max(initial=0.)),
        ambiguous_interval_indices=ambiguous.tolist(),
        comparison_scope="No integer branch correction inferred from CM; ordinary unwrap is untrusted if sampling is ambiguous",
        modulo_phase_increment_difference_rad=float(abs(np.angle(np.exp(1j*(np.diff(wrapped)-physical_increment)))).max(initial=0.)))
    return unwrapped,info

def boundary_edges_for_mesh(tri):
    counts=Counter(tuple(sorted(edge)) for row in tri for edge in zip(row,np.roll(row,-1)))
    edges=np.array([edge for row in tri for edge in zip(row,np.roll(row,-1))
        if counts[tuple(sorted(edge))]==1],int)
    return edges

def boundary_charge(delta,tri,edges=None):
    if edges is None:edges=boundary_edges_for_mesh(tri)
    phase=np.angle(np.conj(delta[edges[:,0]])*delta[edges[:,1]]).sum()/(2*np.pi)
    return int(round(float(phase))),float(phase)

def common_field_l2(snap_a,snap_b,index_a,index_b,center_a_m,center_b_m,dc_gap_J,
                    sample_spacing_nm=.5):
    """Analysis only: same central physical coordinates; complex interpolation first.
    Uniform midpoint spatial quadrature. Return error / initial DC gap; no alignment
    of vortex positions, field evolution, event interpolation or hidden smoothing.
    """
    import matplotlib.tri as mtri
    nx=int(np.ceil(120/sample_spacing_nm));ny=int(np.ceil(80/sample_spacing_nm))
    x=(-60+(np.arange(nx)+.5)*120/nx)*1e-9
    y=(-40+(np.arange(ny)+.5)*80/ny)*1e-9
    xx,yy=np.meshgrid(x,y);xq,yq=xx.ravel(),yy.ravel()
    fields=[]
    for s,index,c in ((snap_a,index_a,center_a_m),(snap_b,index_b,center_b_m)):
        xy=np.asarray(s["coordinates_m"]).copy()
        xy[:,0]-=float(c[0])
        mesh=mtri.Triangulation(xy[:,0],xy[:,1],s["triangles"])
        z=np.asarray(s["delta_J"])[index]
        real=mtri.LinearTriInterpolator(mesh,z.real)(xq,yq)
        imag=mtri.LinearTriInterpolator(mesh,z.imag)(xq,yq)
        if np.ma.getmaskarray(real).any() or np.ma.getmaskarray(imag).any():
            raise ValueError("Common analysis window extends outside a saved domain")
        fields.append(abs(np.asarray(real)+1j*np.asarray(imag)))
    return dict(relative_L2=float(np.sqrt(np.mean((fields[1]-fields[0])**2))/dc_gap_J),
        normalization="initial fresh DC gap (J)",
        region_nm=[-60,60,-40,40],sample_spacing_actual_nm=[120/nx,80/ny],
        interpolation="saved complex linear Delaunay interpolant, followed by magnitude")

def frozen_capture(runpath,out):
    """Read immutable byte images from atomically replaced live output files."""
    target=out/"inputs";target.mkdir(parents=True,exist_ok=False)
    hashes={}
    for name in ("manifest.json","snapshots.npz","checkpoint_latest.npz"):
        source=runpath/name;content=source.read_bytes();(target/name).write_bytes(content)
        hashes[name]=dict(source=str(source),sha256=hashlib.sha256(content).hexdigest(),
                          byte_count=len(content))
    manifest=json.loads((target/"manifest.json").read_text())
    with np.load(target/"checkpoint_latest.npz",allow_pickle=False) as z:
        metadata=json.loads(str(z["metadata_json"]))
        history={key.split("__",1)[1]:z[key].copy() for key in z.files if key.startswith("trajectory__")}
    with np.load(target/"snapshots.npz",allow_pickle=False) as z:
        snap={key:z[key].copy() for key in z.files}
    if manifest["source_pins_sha256"]!=metadata["source_pins_sha256"]:
        raise ValueError("Captured manifest/checkpoint source pins disagree")
    if manifest["parameter_identity"]!=metadata["parameter_identity"]:
        raise ValueError("Captured physical/numerical parameter identities disagree")
    if manifest["source"]!=metadata["source"]:
        raise ValueError("Captured photon source receipts disagree")
    meshpath=Path(manifest["geometry"]["source_path"])
    expected=manifest["geometry"]["mesh_sha256"]
    if sha(meshpath)!=expected:raise ValueError("Original pinned mesh changed")
    content=meshpath.read_bytes();(target/"mesh.npz").write_bytes(content)
    hashes["mesh.npz"]=dict(source=str(meshpath),sha256=expected,byte_count=len(content))
    with np.load(io.BytesIO(content),allow_pickle=False) as z:
        areas=z["dual_area_m2"].copy() if "dual_area_m2" in z else z["area_weights"]*float(z["ell0_m"])**2
        if not np.array_equal(z["coordinates_m"],snap["coordinates_m"]):
            raise ValueError("Saved snapshots and pinned mesh coordinates disagree")
    ht=increasing(history["time_s"]);st=increasing(snap["time_s"])
    keep=st<=ht[-1]+1e-24
    excluded=int((~keep).sum())
    for key,value in list(snap.items()):
        if key not in ("coordinates_m","triangles","edges") and value.shape[0]==len(st):
            snap[key]=value[keep]
    if snap["time_s"][0]<ht[0]-1e-24:raise ValueError("Snapshot predates captured history")
    return manifest,metadata,history,snap,areas,hashes,excluded

def analyze_run(label,runpath,out):
    out.mkdir(parents=True,exist_ok=False)
    manifest,metadata,h,s,areas,hashes,excluded=frozen_capture(runpath,out)
    xy=s["coordinates_m"];tri,reversed_count=ccw_triangles(xy,s["triangles"])
    t=increasing(h["time_s"]);ts=increasing(s["time_s"]);delta=s["delta_J"]
    center=np.asarray(manifest["source"]["center_m"])
    boundary_edges=boundary_edges_for_mesh(tri)
    weights=central_linear_weights(xy,tri,center)
    d0=abs(delta[0]);normalization=float(weights@d0)
    thickness=float(manifest["physicalparams"]["film_thickness_m"])
    volumes=areas*thickness;Eret=float(manifest["source"]["retained_energy_J"])
    left=np.flatnonzero(abs(xy[:,0]-xy[:,0].min())<1e-18)
    right=np.flatnonzero(abs(xy[:,0]-xy[:,0].max())<1e-18)
    phase=[];cut_coherence=[];cut_amplitude=[]
    for nodes in (left,right):
        ratio=delta[:,nodes]/delta[0,nodes]
        unit=ratio/abs(ratio)
        mean=unit.mean(axis=1)
        phase.append(np.unwrap(np.angle(mean)))
        cut_coherence.append(float(np.max(abs(unit-mean[:,None]))))
        cut_amplitude.append(float(np.max(abs(abs(ratio)-1))))
    relative_phase_wrapped=np.angle(np.exp(1j*(phase[1]-phase[0])))
    impulse=cumulative_trapezoid(h["Vdev_V"],t)
    impulse_at_snap=np.interp(ts,t,impulse)
    circ=manifest["physicalparams"]["circuit"];L=float(circ["Lk_ext_H"]);R=float(circ["R_load_ohm"])
    correction=cumulative_trapezoid(h["vc_V"]+h["Vout_V"],t)
    Is_delta=h["Is_A"]-h["Is_A"][0];Ib_delta=h["Ib_A"]-h["Ib_A"][0]
    Is_induction=(correction-impulse)/L
    Vout_induction=R*(Ib_delta-Is_induction)
    Vout_saved=h["Vout_V"]-h["Vout_V"][0]
    dt=np.diff(t)
    midpoint=lambda value:(value[1:]+value[:-1])/2
    Vdev_mid=midpoint(h["vc_V"])+R*(midpoint(h["Ib_A"])-midpoint(h["Is_A"]))-L*np.diff(h["Is_A"])/dt
    CM_impulse=np.r_[0.,np.cumsum(Vdev_mid*dt)]
    CM_work=np.r_[0.,np.cumsum(midpoint(h["Is_A"])*Vdev_mid*dt)]
    CM_impulse_at_snap=np.interp(ts,t,CM_impulse)
    phase_unwrapped,phase_info=phase_unwrap_diagnostic(relative_phase_wrapped,CM_impulse_at_snap)
    phase_flux=phase_unwrapped*PHI0/(2*np.pi)
    phase_info["ambiguous_saved_time_brackets_ps"]=[
        [float(ts[index]*1e12),float(ts[index+1]*1e12)]
        for index in phase_info["ambiguous_interval_indices"]]
    phase_certified=phase_info["certified"]
    Pb=float(circ["V_bias_V"])*h["Ib_A"]-float(circ["R_bias_ohm"])*h["Ib_A"]**2-R*(h["Ib_A"]-h["Is_A"])**2
    Pdevice=h["Is_A"]*h["Vdev_V"]
    Pescape=h["escape_power_W"];Pcuts=h["reservoir_power_out_W"]
    rows=[];largest_winding_error=0.;previous=[];next_id=1
    for index,time_s in enumerate(ts):
        charges,minima,vortices,error=triangle_topology(delta[index],xy,tri)
        largest_winding_error=max(largest_winding_error,error)
        total,boundary_raw=boundary_charge(delta[index],tri,boundary_edges)
        if total!=int(charges.sum()):raise ValueError("Boundary winding != sum of CCW triangle charges")
        unused=set(range(len(previous)))
        for v in vortices:
            pos=np.asarray(v["linear_complex_zero_nm"] if v["linear_complex_zero_nm"] is not None else v["centroid_nm"])
            options=[(np.linalg.norm(pos-np.asarray(previous[j]["position"])),j)
                     for j in unused if previous[j]["charge"]==v["charge"]]
            match=min(options,default=(np.inf,-1))
            if match[0]<=20.:
                v["track_id"]=previous[match[1]]["id"];unused.remove(match[1])
            else:v["track_id"]=next_id;next_id+=1
        previous=[dict(id=v["track_id"],charge=v["charge"],
            position=v["linear_complex_zero_nm"] if v["linear_complex_zero_nm"] is not None else v["centroid_nm"]) for v in vortices]
        deficit=float(weights@np.maximum(d0-abs(delta[index]),0))
        u=s["u_excess_J_m3"][index]
        positive=np.maximum(u,0)*volumes;mass=float(positive.sum())
        centroid=(positive@xy/mass) if mass>0 else np.zeros(2)
        sigma=np.sqrt(positive@((xy-centroid)**2)/mass) if mass>0 else np.zeros(2)
        near_core=[]
        for v in vortices:
            point=np.asarray(v["linear_complex_zero_nm"] if v["linear_complex_zero_nm"] is not None else v["centroid_nm"])*1e-9
            nodes=np.flatnonzero(np.linalg.norm(xy-point,axis=1)<=10e-9)
            near_core.append(dict(track_id=v["track_id"],radius_nm=10,
                excess_energy_J=float(volumes[nodes]@u[nodes]),
                Te_equivalent_max_K=float(s["Te_K"][index,nodes].max()),
                Tph_equivalent_max_K=float(s["Tph_K"][index,nodes].max()),
                note="Neighborhoods may overlap; their energies must not be summed"))
        rows.append(dict(time_s=float(time_s),time_ps=float(time_s*1e12),
            total_charge=int(charges.sum()),n_positive=int((charges>0).sum()),n_negative=int((charges<0).sum()),
            boundary_winding_raw=boundary_raw,vortices=vortices,
            minimum_nodal_gap_meV=float(abs(delta[index]).min()/MEV_J),
            minimum_complex_linear_gap_meV=float(minima.min()/MEV_J),
            central_deficit_meV_nm2=deficit/MEV_J*1e18,
            central_deficit_fraction_DC=deficit/normalization,
            electron_excitation_energy_J=float(volumes@s["electron_energy_J_m3"][index]),
            phonon_energy_J=float(volumes@s["phonon_energy_J_m3"][index]),
            excess_signed_energy_J=float(volumes@u),
            excess_positive_energy_J=mass,positive_excess_centroid_nm=centroid*1e9,
            positive_excess_sigma_nm=sigma*1e9,
            Te_equivalent_max_K=float(s["Te_K"][index].max()),
            Tph_equivalent_max_K=float(s["Tph_K"][index].max()),
            core_neighborhoods=near_core,
            phase_cut_impulse_Phi0=float(phase_flux[index]/PHI0) if phase_certified else None,
            saved_voltage_integral_Phi0=float(impulse_at_snap[index]/PHI0)))
    active=[i for i,r in enumerate(rows) if r["vortices"]]
    first=active[0] if active else None
    last=active[-1] if active else None
    first_bracket=None if first is None else [float(ts[max(0,first-1)]*1e12),float(ts[first]*1e12)]
    energy_integral=np.asarray([r["excess_signed_energy_J"] for r in rows])
    film_expected=Eret+np.interp(ts,t,h["domain_energy_J"]-h["domain_energy_J"][0])
    summary=dict(schema="pysnspd.saved_phase_vortex_audit.v1",label=label,run_directory=str(runpath),
        admission=False,stage5=False,scope="Saved-data diagnostics only; no evolution, gate or event inferred",
        gauge="A=0; CCW positive phase winding; no magnetic self-field or core resolution inferred",
        source_hashes=hashes,source_pins=metadata["source_pins_sha256"],
        captured_checkpoint_time_ps=float(metadata["time_ps"]),
        analyzed_snapshots_through_ps=float(ts[-1]*1e12),excluded_snapshots_later_than_checkpoint=excluded,
        geometry=dict(reoriented_CCW_triangles=reversed_count,central_region_relative_nm=[-60,60,-40,40],
            central_origin_m=center,central_area_nm2=float(weights.sum()*1e18),
            gap_deficit_rule="D_h,nodal=sum clipped-triangle weights * positive nodal gap deficit; exact integral of that linear nodal diagnostic, spatial approximation to continuum deficit",
            gap_deficit_normalization_J_m2=normalization),
        topology=dict(first_saved_vortex_bracket_ps=first_bracket,
            first_saved_vortex=None if first is None else rows[first],
            last_saved_vortex=None if last is None else rows[last],
            frames_with_vortices=len(active),largest_integer_winding_error=largest_winding_error,
            maximum_absolute_total_charge=max(abs(r["total_charge"]) for r in rows),
            limitations=["Appearance/exit only bracketed by saved snapshots",
                "Track IDs use nearest same-sign zero within 20nm between adjacent saved frames",
                "Linear complex zeros locate interpolation zeros, not resolved microscopic cores",
                "Nodal gap oscillations during core motion do not establish recovery",
                "Vortex pair is not a width-spanning hotbelt or confirmed electrical detection"]),
        electrical=dict(Phi0_V_s=PHI0,Lext_H=L,Rload_ohm=R,
            one_Phi0_short_time_delta_Is_A=-PHI0/L,
            one_Phi0_short_time_delta_Vout_V=R*PHI0/L,
            threshold_V=float(manifest["observation"]["config"]["threshold_V"]),
            threshold_in_one_Phi0_inductive_units=float(manifest["observation"]["config"]["threshold_V"]/(R*PHI0/L)),
            final_saved_voltage_integral_Phi0=float(impulse[-1]/PHI0),
            phase_unwrap=phase_info,
            final_cut_phase_impulse_Phi0=float(phase_flux[-1]/PHI0) if phase_certified else None,
            final_cut_phase_time_ps=float(ts[-1]*1e12),
            final_endpoint_integral_time_ps=float(t[-1]*1e12),
            endpoint_integral_at_last_snapshot_Phi0=float(impulse_at_snap[-1]/PHI0),
            CM_impulse_at_last_snapshot_Phi0=float(CM_impulse_at_snap[-1]/PHI0),
            max_phase_vs_CM_impulse_difference_Phi0=float(max(abs(phase_flux-CM_impulse_at_snap))/PHI0) if phase_certified else None,
            CM_reconstruction_rule="Vdev_mid=vc_mid+Rload*(Ib_mid-Is_mid)-Lext*deltaIs/dt, from saved linear CM states",
            maximum_cut_phase_vs_endpoint_trapezoid_difference_Phi0=float(max(abs(phase_flux-impulse_at_snap))/PHI0) if phase_certified else None,
            cut_phase_relation="Integral Vdev dt = Phi0*(change theta_right - change theta_left)/(2pi)",
            integration_rule="Trapezoid of saved endpoint voltages; actual midpoint voltage is not stored",
            cut_relative_amplitude_error=cut_amplitude,cut_phase_nonuniformity=cut_coherence,
            cut_max_adjacent_phase_increment_rad=[float(abs(np.diff(q)).max(initial=0)) for q in phase],
            max_abs_Vout_mV=float(max(abs(Vout_saved))*1e3),
            max_abs_Vdev_mV=float(max(abs(h["Vdev_V"]))*1e3),
            final_Vout_mV=float(Vout_saved[-1]*1e3),
            final_current_Is_A=float(h["Is_A"][-1]),
            final_current_Ib_A=float(h["Ib_A"][-1]),
            max_CM_endpoint_integral_current_defect_A=float(max(abs(Is_delta-Is_induction))),
            max_CM_endpoint_integral_readout_defect_V=float(max(abs(Vout_saved-Vout_induction))),
            CM_relation="Lext*delta Is = Integral(vc+Vout-Vdev)dt; Vout=Rload*(Ib-Is)",
            no_trigger_inferred=True),
        energy=dict(retained_energy_J=Eret,
            final_device_work_J=float(cumulative_trapezoid(Pdevice,t)[-1]),
            final_device_work_from_CM_midpoints_J=float(CM_work[-1]),
            final_escape_energy_J=float(cumulative_trapezoid(Pescape,t)[-1]),
            final_cuts_energy_out_J=float(cumulative_trapezoid(Pcuts,t)[-1]),
            final_net_circuit_power_integral_J=float(cumulative_trapezoid(Pb,t)[-1]),
            final_external_integral_J=float(h["external_energy_integral_J"][-1]),
            max_snapshot_vs_film_ledger_difference_J=float(max(abs(energy_integral-film_expected))),
            saved_Joule_heat_density_available=False,
            note="Electron population and phonon energies are saved absolute components, not full condensate energy; Te/Tph are energy-equivalent diagnostics. Signed excess retains negative regions."),
        rows=rows)
    np.savez_compressed(out/"analysis_traces.npz",time_s=t,snapshot_time_s=ts,
        voltage_integral_V_s=impulse,phase_cut_integral_V_s=phase_flux,
        phase_unwrap_certified=np.array(phase_certified),relative_phase_wrapped_rad=relative_phase_wrapped,
        CM_correction_integral_V_s=correction,CM_reconstructed_impulse_V_s=CM_impulse,CM_device_work_J=CM_work,
        Is_delta_A=Is_delta,Is_from_integrated_CM_A=Is_induction,
        Vout_relative_V=Vout_saved,Vout_from_integrated_CM_V=Vout_induction,
        Vout_short_time_V=R*impulse/L,
        device_work_J=cumulative_trapezoid(Pdevice,t),escape_energy_J=cumulative_trapezoid(Pescape,t),
        cut_energy_out_J=cumulative_trapezoid(Pcuts,t),net_circuit_energy_J=cumulative_trapezoid(Pb,t),
        central_linear_weights_m2=weights)
    (out/"analysis.json").write_text(json.dumps(summary,indent=2,default=native)+"\n")
    return summary,manifest,h,s,tri

def plots(summary,manifest,h,s,tri,out):
    os.environ.setdefault("MPLCONFIGDIR",str(SCRATCH/"tmp/pysnspd/matplotlib"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.tri import Triangulation
    import sys
    sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
    from pysnspd.plotting.complex_gap import draw_complex_gap
    rows=summary["rows"];times=np.array([r["time_ps"] for r in rows])
    with np.load(out/"analysis_traces.npz",allow_pickle=False) as z:a={k:z[k].copy() for k in z.files}
    t=a["time_s"]*1e12;ts=a["snapshot_time_s"]*1e12
    wave=str(float(manifest["source"]["wavelength_nm"]))
    title=wave+" nm | diagnóstico guardado A=0; sin admisión ni gatillo inferido"
    records=[]
    def save(fig,name,caption):
        path=out/(name+".png");fig.savefig(path,dpi=170,facecolor="white");plt.close(fig)
        records.append(dict(path=str(path),sha256=sha(path),caption=caption))
    fig,axes=plt.subplots(3,2,figsize=(11,9),layout="constrained")
    axes[0,0].plot(times,[r["minimum_nodal_gap_meV"] for r in rows],label="Mínimo nodal")
    axes[0,0].plot(times,[r["minimum_complex_linear_gap_meV"] for r in rows],label="Mínimo complejo lineal")
    axes[0,0].set_ylabel("Gap mínimo (meV)")
    axes[0,1].step(times,[r["n_positive"] for r in rows],where="post",label="Carga +1 (estados guardados)")
    axes[0,1].step(times,[r["n_negative"] for r in rows],where="post",label="Carga -1 (estados guardados)")
    axes[0,1].set_ylabel("Triángulos con winding")
    axes[1,0].plot(t,a["voltage_integral_V_s"]/PHI0,label="Integral Vdev, extremos")
    if summary["electrical"]["phase_unwrap"]["certified"]:
        axes[1,0].plot(ts,a["phase_cut_integral_V_s"]/PHI0,ls="--",label="Fases de cortes: unwrap certificado")
    else:
        axes[1,0].text(.02,.95,"FASE SUBMUESTREADA: unwrap ambiguo",transform=axes[1,0].transAxes,
            va="top",color="tab:red",fontsize=9)
    axes[1,0].plot(t,a["CM_reconstructed_impulse_V_s"]/PHI0,ls=":",label="Midpoints desde estados CM")
    axes[1,0].set_ylabel("Impulso eléctrico / Φ0")
    axes[1,1].plot(t,a["Vout_relative_V"]*1e3,label="Vout guardado")
    axes[1,1].plot(t,a["Vout_from_integrated_CM_V"]*1e3,ls="--",label="CM con trapecios de extremos")
    axes[1,1].plot(t,a["Vout_short_time_V"]*1e3,ls=":",label="Rload/Lext × integral Vdev")
    axes[1,1].set_ylabel("Señal de lectura (mV)")
    axes[2,0].plot(times,[r["central_deficit_fraction_DC"]*100 for r in rows])
    axes[2,0].set_ylabel("Déficit central / (gap DC × área) (%)")
    axes[2,1].plot(times,np.array([r["excess_signed_energy_J"] for r in rows])/EV_J,label="Exceso firmado")
    axes[2,1].plot(times,np.array([r["electron_excitation_energy_J"] for r in rows])/EV_J,label="Excitaciones electrónicas")
    axes[2,1].plot(times,np.array([r["phonon_energy_J"] for r in rows])/EV_J,label="Fonones")
    axes[2,1].set_ylabel("Energía del film (eV)")
    for ax in axes.ravel():
        ax.set_xlabel("Tiempo desde transferencia (ps)");ax.grid(alpha=.2)
        if ax.get_legend_handles_labels()[0]:ax.legend(fontsize=8)
    fig.suptitle(title,fontsize=12)
    save(fig,"phase_voltage_energy","Winding CCW con A=0. El mínimo complejo lineal puede ser cero entre nodos; no demuestra resolución del núcleo. La integral de Vdev incluye respuesta reactiva y no se redondea a un número de slips. CM conserva su inductancia exterior fija; un vórtice no es un gatillo.")
    active=[i for i,r in enumerate(rows) if r["vortices"]]
    if active:
        selected=sorted(set([max(0,active[0]-1),active[0],active[len(active)//2],active[-1],len(rows)-1]))
    else:selected=sorted(set([0,len(rows)//2,len(rows)-1]))
    xy=s["coordinates_m"]*1e9;center=np.asarray(manifest["source"]["center_m"])*1e9
    xrel=xy[:,0]-center[0];yrel=xy[:,1]
    gap_tri=Triangulation(xrel,yrel,tri)
    fig,axes=plt.subplots(len(selected),3,figsize=(11,2.6*len(selected)+.5),squeeze=False,layout="constrained")
    gapmax=float(abs(s["delta_J"]).max()/MEV_J)
    energymax=float(abs(s["u_excess_J_m3"]).max())
    for row,index in enumerate(selected):
        z=s["delta_J"][index]
        fields=(abs(z)/MEV_J,np.angle(z[tri].sum(axis=1)),s["u_excess_J_m3"][index])
        for col,ax in enumerate(axes[row]):
            if col==0:im=draw_complex_gap(ax,gap_tri,z,unit_J=MEV_J,
                vmin=0,vmax=gapmax,bounds=(-60,60,-40,40))
            elif col==1:im=ax.tripcolor(xrel,yrel,tri,facecolors=fields[col],shading="flat",cmap="twilight",vmin=-np.pi,vmax=np.pi)
            else:im=ax.tripcolor(xrel,yrel,tri,fields[col],shading="gouraud",cmap="coolwarm",vmin=-energymax,vmax=energymax)
            for v in rows[index]["vortices"]:
                pt=np.asarray(v["linear_complex_zero_nm"] if v["linear_complex_zero_nm"] is not None else v["centroid_nm"])-center
                ax.plot(pt[0],pt[1],"+" if v["charge"]>0 else "x",color="white" if col==0 else "black",ms=8,mew=1.5)
                ax.annotate(f'{v["charge"]:+d}',pt,xytext=(4,2),textcoords="offset points",color="white" if col==0 else "black",fontsize=8)
            ax.set(aspect="equal",xlim=(-60,60),ylim=(-40,40),xlabel="x − xγ (nm)",ylabel="y físico (nm)")
            ax.set_title(f'{times[index]:.5g} ps | '+("|P1 Δ complejo| (meV)" if col==0 else "arg Δ (rad), A=0" if col==1 else "δu (J/m³)"),fontsize=10)
            if row==len(selected)-1:fig.colorbar(im,ax=axes[:,col].tolist(),shrink=.5)
    fig.suptitle(title,fontsize=12)
    save(fig,"vortex_phase_gap_energy","Estados realmente guardados. Gap: módulo del interpolante espacial complejo P1 muestreado en raster de hasta601 puntos por eje, con celdas sin suavizado; esta resolución visual no prueba resolución del núcleo. Energía: interpolación visual nodal; fase: argumento de la media compleja de cada triángulo, sin interpolar ángulos a través de ±π. Símbolos: winding CCW y ceros del interpolante complejo, con incertidumbre de malla. Recorte de 120×80nm sólo visual. No se asimila una pareja a una hotbelt.")
    (out/"figure_manifest.json").write_text(json.dumps(dict(schema="pysnspd.saved_vortex_figures.v1",stage5=False,figures=records),indent=2)+"\n")
    return records

def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run",action="append",required=True,help="NAME=/home/jdiaz/scratch/run")
    parser.add_argument("--output-dir",required=True)
    parser.add_argument("--plots",action="store_true")
    parser.add_argument("--cpu",type=int,default=13)
    args=parser.parse_args(argv)
    if hasattr(os,"sched_setaffinity"):os.sched_setaffinity(0,{args.cpu})
    out=Path(args.output_dir).expanduser().resolve()
    if not out.is_relative_to(SCRATCH):raise ValueError("Analysis output must be in scratch")
    out.mkdir(parents=True,exist_ok=False);reports=[]
    for spec in args.run:
        name,path=spec.split("=",1);runpath=Path(path).resolve()
        if not runpath.is_relative_to(SCRATCH):raise ValueError("Saved inputs must be in scratch")
        summary,manifest,h,s,tri=analyze_run(name,runpath,out/name)
        if args.plots:plots(summary,manifest,h,s,tri,out/name)
        reports.append(dict(label=name,path=str(out/name/"analysis.json"),sha256=sha(out/name/"analysis.json"),
            first_saved_vortex_bracket_ps=summary["topology"]["first_saved_vortex_bracket_ps"],
            saved_through_ps=summary["analyzed_snapshots_through_ps"],
            endpoint_integral_Phi0=summary["electrical"]["final_saved_voltage_integral_Phi0"],
            Vout_max_mV=summary["electrical"]["max_abs_Vout_mV"]))
    receipt=dict(schema="pysnspd.saved_vortex_campaign.v1",stage5=False,
        helper_path=str(Path(__file__).resolve()),helper_sha256=sha(__file__),
        CPU_affinity=sorted(os.sched_getaffinity(0)),inputs=reports,
        independent_reference="https://arxiv.org/html/1611.06060v3",
        note="Vodolazov2017 motivates distinct hot-spot instability and normal-domain stages; these saved effective states are not a microscopic validation.")
    (out/"receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps(receipt,indent=2))
if __name__=="__main__":main()
