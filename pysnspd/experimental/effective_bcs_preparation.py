"""Fresh DC preparation on the inherited Delaunay--Voronoi strip.

The old mesh is reused exactly, but its microscopic DC state is not relabeled.
This module solves the effective BCS thermal gap force and the imposed net
current together. Only the two complete longitudinal cuts have prescribed
bulk amplitude; their phases are equipotential, with the left phase adjusted
to current and the right phase fixing the common gauge.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.constants import Boltzmann as kB, hbar
from scipy.sparse import bmat, coo_matrix, diags
from scipy.sparse.linalg import spsolve
from .thermal_spatial_usadel import ThermalGraph
from .effective_bcs_functional import EffectiveBCSFunctional, EffectiveBCSBulk

REPO=Path(__file__).resolve().parents[2]
SCRATCH=Path("/home/jdiaz/scratch")
BASE_MESH=REPO/"docs/implementation/stage4/practical_time_review_20260924/dual_mesh/resampled/mesh.npz"
BASE_MESH_SHA256="da1abb237cbe0fe6ac4260db9949144e5a2a7f4fecae6766a0def68fa5e2ccf4"
HISTORICAL_MESHES={
    (160.,5.):BASE_MESH,
    (160.,3.5):REPO/"docs/implementation/stage5/prephoton_dc_20260924/meshes/l160_h35/mesh.npz",
    (320.,5.):REPO/"docs/implementation/stage5/prephoton_dc_20260924/meshes/l320_h5/mesh.npz",
}


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


@dataclass(frozen=True)
class EffectiveBCSGeometry:
    graph: ThermalGraph
    coordinates_m: np.ndarray
    triangles: np.ndarray
    edge_lengths_m: np.ndarray
    face_areas_m2_per_thickness_m: np.ndarray
    left_nodes: np.ndarray
    right_nodes: np.ndarray
    ell0_m: float
    length_m: float
    width_m: float
    source_path: str
    receipt: dict

    @property
    def dual_face_lengths_m(self):
        return self.face_areas_m2_per_thickness_m

    def face_areas_m2(self, thickness_m):
        return self.dual_face_lengths_m*float(thickness_m)


@dataclass(frozen=True)
class EffectiveBCSDC:
    delta_bar: np.ndarray
    p: np.ndarray
    phi_V: np.ndarray
    bulk: EffectiveBCSBulk
    phase_left_rad: float
    receipt: dict

    @property
    def gap_bar(self):
        return self.delta_bar


def load_dual_geometry(path, *, expected_sha256=None, thickness_m=7e-9):
    """Read the original geometry arrays without rebuilding or remeshing."""
    path=Path(path)
    digest=_sha(path)
    if expected_sha256 is not None and digest!=expected_sha256:
        raise ValueError("Pinned Delaunay mesh changed: "+str(path))
    with np.load(path,allow_pickle=False) as data:
        arrays={name:data[name].copy() for name in data.files}
    for name in ("coordinates_m","coordinates_bar","triangles","edges",
                 "area_weights","conductance","edge_length_m","dual_face_length_m",
                 "fixed_nodes","left_nodes","right_nodes","ell0_m"):
        if name not in arrays:
            raise ValueError("Missing inherited dual geometry field: "+name)
    xy=arrays["coordinates_m"]
    ell=float(arrays["ell0_m"])
    left,right=arrays["left_nodes"],arrays["right_nodes"]
    graph=ThermalGraph(arrays["area_weights"],arrays["edges"],arrays["conductance"],
                       arrays["coordinates_bar"],arrays["fixed_nodes"])
    lengths=arrays["edge_length_m"];dual=arrays["dual_face_length_m"]
    if (xy.shape!=(graph.n_nodes,2) or np.any(~np.isfinite(xy))
            or lengths.shape!=(len(graph.edges),) or np.any(lengths<=0)
            or dual.shape!=lengths.shape or np.any(dual<0)
            or not np.allclose(graph.conductance,dual/lengths,rtol=1e-13,atol=1e-15)):
        raise ValueError("Geometry does not preserve the same dual/primal conductance")
    length=float(np.ptp(xy[:,0]));width=float(np.ptp(xy[:,1]))
    area_error=float(graph.area_weights.sum()*ell**2/(length*width)-1)
    if abs(area_error)>1e-10:
        raise ValueError("Dual areas do not cover the physical strip")
    receipt={"schema":"pysnspd.effective_bcs_geometry.v1","source_path":str(path),
             "mesh_sha256":digest,"n_nodes":graph.n_nodes,"n_edges":len(graph.edges),
             "length_m":length,"width_m":width,"ell0_m":ell,"area_relative_error":area_error,
             "triangulation":"inherited pyTDGL-family Delaunay--Voronoi",
             "new_mesh_generated":False}
    return EffectiveBCSGeometry(graph,xy,arrays["triangles"],lengths,dual,left,right,
                               ell,length,width,str(path),receipt)


def _boundary_sampling(length,width,spacing):
    perimeter=2*(length+width);half=spacing/2
    count=max(8,int(np.ceil(perimeter/half)))
    while True:
        nx=round(count*length/perimeter);ny=round(count*width/perimeter)
        if nx>1 and ny>1 and max(length/(nx-1),width/(ny-1))<=half*(1+1e-14):
            return count
        count+=1


def get_effective_geometry(length_m,spacing_m,*,width_m=80e-9,Tc_K=8.65,
                           diffusion_m2_s=5e-5,output_root=None):
    """Reuse existing cases, otherwise run exactly the production-family mesher.

    Every new numeric mesh is saved under scratch. Boundary sampling preserves
    the existing half-target-edge rule; no seed or outcome search is performed.
    """
    length,spacing,width=map(float,(length_m,spacing_m,width_m))
    if not np.all(np.isfinite([length,spacing,width,Tc_K,diffusion_m2_s])) or min(
            length,spacing,width,Tc_K,diffusion_m2_s)<=0:
        raise ValueError("positive physical geometry and material required")
    key=(round(length*1e9,8),round(spacing*1e9,8))
    historical=HISTORICAL_MESHES.get(key)
    ell=np.sqrt(hbar*diffusion_m2_s/(2*kB*Tc_K))
    if historical is not None and historical.is_file() and np.isclose(width,80e-9,rtol=1e-14,atol=0.):
        expected=BASE_MESH_SHA256 if historical==BASE_MESH else None
        if expected is None and (historical.parent/"morphology.json").is_file():
            expected=json.loads((historical.parent/"morphology.json").read_text())["mesh_sha256"]
        result=load_dual_geometry(historical,expected_sha256=expected)
        if not np.isclose(result.ell0_m,ell,rtol=1e-9,atol=0.):
            raise ValueError("Historical geometry uses a different physical length scale")
        return result
    output=(SCRATCH/"stage5_effective_bcs_development_20261011/meshes"
            if output_root is None else Path(output_root)).resolve()
    if not output.is_relative_to(SCRATCH.resolve()):
        raise ValueError("New heavy mesh arrays must live under /home/jdiaz/scratch")
    name=f"l{length*1e9:g}_w{width*1e9:g}_h{spacing*1e9:g}"
    directory=output/name
    target=directory/"mesh.npz"
    if target.is_file():
        metadata=json.loads((directory/"geometry.json").read_text())
        result=load_dual_geometry(target,expected_sha256=metadata["mesh_sha256"])
        if not np.isclose(result.ell0_m,ell,rtol=1e-9,atol=0.):
            raise ValueError("Cached geometry uses a different physical length scale")
        return result
    from pysnspd.mesh.pytdgl_like import (
        PyTDGLLikeMeshParameters,generate_rectangular_pytdgl_fvm_mesh_from_parameters)
    from .thermal_mesh_adapter import graph_from_mesh
    params=PyTDGLLikeMeshParameters(length_m=length,width_m=width,
        target_spacing_m=spacing,max_edge_length_m=spacing,seed=12345,
        boundary_points=_boundary_sampling(length,width,spacing))
    mesh=generate_rectangular_pytdgl_fvm_mesh_from_parameters(params)
    tolerance=64*np.finfo(float).eps*length
    left=np.flatnonzero(abs(mesh.sites[:,0])<=tolerance)
    right=np.flatnonzero(abs(mesh.sites[:,0]-length)<=tolerance)
    if not len(left) or not len(right):
        raise ValueError("Both full-width longitudinal cuts are required")
    fixed=np.unique(np.r_[left,right]);origin=np.array([length/2,0.])
    graph=graph_from_mesh(mesh,ell,fixed_nodes=fixed,origin_m=origin)
    directory.mkdir(parents=True,exist_ok=False)
    np.savez_compressed(target,coordinates_m=mesh.sites,
        coordinates_bar=graph.coordinates_bar,triangles=mesh.elements,edges=graph.edges,
        area_weights=graph.area_weights,conductance=graph.conductance,
        dual_area_m2=mesh.areas,edge_length_m=mesh.edge_mesh.edge_lengths,
        dual_face_length_m=mesh.edge_mesh.dual_edge_lengths,
        geometric_boundary_nodes=mesh.boundary_indices,fixed_nodes=fixed,left_nodes=left,
        right_nodes=right,origin_m=origin,ell0_m=np.array(ell))
    metadata={"schema":"pysnspd.effective_bcs_generated_dual_mesh.v1",
              "mesh_sha256":_sha(target),"parameters":asdict(params),
              "boundary_rule":"half target interior edge length",
              "seed_scope":"metadata retained; same meshpy construction, no seed search"}
    (directory/"geometry.json").write_text(json.dumps(metadata,indent=2)+"\n")
    result=load_dual_geometry(target,expected_sha256=metadata["mesh_sha256"])
    result.receipt.update(new_mesh_generated=True,parameters=asdict(params))
    return result


def prepare_effective_dc(model,*,left_nodes,right_nodes,temperature_K,current_A,
                         width_m,sheet_resistance_ohm=608.,bulk_reference=None,
                         tolerance=1e-10,max_iterations=60):
    """Solve thermal stationarity and net supercurrent together on this graph.

    The thermal Hessian includes the FD occupation response. It is the Jacobian
    of the same free-energy gap force, not a new stiffness or an old DC field.
    The result must pass independent quasineutral port reconstruction.
    """
    if not isinstance(model,EffectiveBCSFunctional):
        raise TypeError("EffectiveBCSFunctional required")
    graph=model.graph;n=graph.n_nodes
    left,right=np.asarray(left_nodes,int),np.asarray(right_nodes,int)
    if (not len(left) or not len(right) or len(np.intersect1d(left,right))
            or np.any(left<0) or np.any(left>=n) or np.any(right<0) or np.any(right>=n)):
        raise ValueError("Disjoint full-width cuts required")
    coordinate=graph.coordinates_bar[:,0]
    if not np.allclose(coordinate[left],coordinate.min(),rtol=0,atol=1e-12) or not np.allclose(
            coordinate[right],coordinate.max(),rtol=0,atol=1e-12):
        raise ValueError("Cuts must lie at the longitudinal graph extrema")
    fixed=np.unique(np.r_[left,right]);free=np.setdiff1d(np.arange(n),fixed)
    freecomp=np.column_stack((2*free,2*free+1)).ravel()
    bulk=(model.bulk(temperature_K=temperature_K,current_A=current_A,width_m=width_m)
          if bulk_reference is None else bulk_reference)
    d=bulk.amplitude_bar*np.exp(1j*bulk.q_bar*(coordinate-coordinate.max()))
    theta=float(bulk.q_bar*(coordinate.min()-coordinate.max()))
    scale=np.repeat(model.node_energy_scale_J[free],2)
    current_scale=max(abs(float(current_A)),1e-9)
    incidence_left=np.zeros(len(graph.edges))
    t,h=graph.edges.T
    incidence_left+=np.isin(t,left)
    incidence_left-=np.isin(h,left)
    thermal_t=float(temperature_K)/model.Tc_K
    history=[]
    def evaluated(gap,phase):
        gap=gap.copy()
        gap[left]=bulk.amplitude_bar*np.exp(1j*phase)
        gap[right]=bulk.amplitude_bar
        p=model.fermi_dirac(gap,temperature_K)
        value=model.evaluate(gap,p)
        force=value.gap_force_J_per_gapbar[free]/model.node_energy_scale_J[free]
        rows=np.r_[np.column_stack((force.real,force.imag)).ravel(),
                    (np.dot(incidence_left,value.edge_current_A)-current_A)/current_scale]
        return gap,p,value,rows
    for iteration in range(max_iterations+1):
        d,p,value,rows=evaluated(d,theta)
        norm=float(np.max(abs(rows)))
        history.append({"iteration":iteration,"scaled_residual_max":norm})
        if norm<=tolerance:
            break
        if iteration==max_iterations:
            raise RuntimeError("Effective DC Newton iteration limit; residual="+str(norm))
        H,Jcurrent=model.hessian_and_current_jacobian(d,p)
        en=value.electron_energies_J/model.Eu_J
        coefficient=-2*model.node_energy_scale_J*np.sum(
            model.weights_bar*p*(1-p)/(thermal_t*en**2),axis=1)
        z=np.column_stack((d.real,d.imag))
        correction=coefficient[:,None,None]*np.einsum("ni,nj->nij",z,z)
        component=2*np.arange(n)[:,None]+np.arange(2)
        rr=np.broadcast_to(component[:,:,None],correction.shape).ravel()
        cc=np.broadcast_to(component[:,None,:],correction.shape).ravel()
        H=H+coo_matrix((correction.ravel(),(rr,cc)),shape=H.shape).tocsr()
        theta_direction=np.zeros(2*n)
        theta_direction[2*left]=-d[left].imag
        theta_direction[2*left+1]=d[left].real
        selector=coo_matrix((np.ones(len(freecomp)),(freecomp,np.arange(len(freecomp)))),
                            shape=(2*n,len(freecomp))).tocsr()
        P=bmat([[selector,coo_matrix(theta_direction[:,None])]],format="csr")
        force_jac=diags(1/scale)@(H[freecomp]@P)
        current_row=coo_matrix(np.asarray(incidence_left@Jcurrent@P).ravel()[None,:]/current_scale)
        jac=bmat([[force_jac],[current_row]],format="csc")
        step=spsolve(jac,-rows)
        if np.any(~np.isfinite(step)):
            raise FloatingPointError("Nonfinite effective DC Newton step")
        dd=np.zeros(n,complex)
        dd[free]=step[:-1].reshape(-1,2)[:,0]+1j*step[:-1].reshape(-1,2)[:,1]
        amount=min(1.,.35/max(float(np.max(abs(dd))),.35),
                   .5/max(abs(float(step[-1])),.5))
        accepted=False
        for _ in range(24):
            trial=evaluated(d+amount*dd,theta+amount*step[-1])
            trialnorm=float(np.max(abs(trial[-1])))
            if trialnorm < norm*(1-1e-4*amount) or trialnorm <= tolerance:
                d,theta=trial[0],theta+amount*float(step[-1])
                accepted=True
                break
            amount/=2
        if not accepted:
            raise RuntimeError("Effective DC Newton line search failed; residual="+str(norm))
    from .effective_bcs_ports import EffectiveNeutralPorts
    ports=EffectiveNeutralPorts(graph,sheet_resistance_ohm,left,right)
    electrical=ports.solve(value.edge_current_A,current_A)
    receipt={"schema":"pysnspd.effective_bcs_dc_preparation.v1",
             "status":"FRESH_EFFECTIVE_BCS_DC_PREPARED",
             "current_A":float(current_A),"temperature_K":float(temperature_K),
             "n_nodes":n,"n_free":len(free),"iterations":iteration,
             "scaled_gap_force_max":float(np.max(abs(
                 value.gap_force_J_per_gapbar[free]/model.node_energy_scale_J[free]))),
             "net_supercurrent_left_A":float(np.dot(incidence_left,value.edge_current_A)),
             "current_relative_defect":float(abs(np.dot(incidence_left,value.edge_current_A)-current_A)/current_scale),
             "potential_max_abs_V":float(np.max(abs(electrical.phi_V))),
             "Vdev_V":electrical.Vdev_V,"normal_joule_W":electrical.normal_joule_W,
             "bulk":asdict(bulk),"phase_left_rad":theta,"tolerance":tolerance,
             "gap_amplitude_min_bar":float(abs(d).min()),"gap_amplitude_max_bar":float(abs(d).max()),
             "gap_amplitude_relative_span":float(np.ptp(abs(d))/np.mean(abs(d))),
             "old_DC_state_reused":False,"old_mesh_may_be_reused":True,
             "model":dict(energy_scale_J=model.Eu_J,N0_per_J_m3=model.N0_per_J_m3,
                          ell0_m=model.ell0_m,thickness_m=model.thickness_m,
                          gradient_law=model.gradient_law,bath_temperature_K=model.bath_temperature_K,
                          inverse_pairing_coupling=model.inverse_pairing_coupling,
                          xi_sha256=hashlib.sha256(model.xi_bar.tobytes()).hexdigest(),
                          weights_sha256=hashlib.sha256(model.weights_bar.tobytes()).hexdigest(),
                          functional_sha256=_sha(Path(__file__).with_name("effective_bcs_functional.py"))),
             "history":history}
    if (receipt["scaled_gap_force_max"]>2*tolerance
            or receipt["current_relative_defect"]>2*tolerance
            or receipt["potential_max_abs_V"]>1e-9):
        raise RuntimeError("Fresh DC fails force/current/nV baseline reconstruction: "+str(receipt))
    return EffectiveBCSDC(d,p,electrical.phi_V,bulk,theta,receipt)


def save_effective_dc(output,preparation,geometry):
    """Save numerical preparation only under scratch; preserve existing files."""
    output=Path(output).resolve()
    if not output.is_relative_to(SCRATCH.resolve()):
        raise ValueError("DC numerical data must be stored under scratch")
    output.mkdir(parents=True,exist_ok=True)
    target=output/"effective_dc.npz"
    if target.exists():
        raise FileExistsError("Existing preparation is preserved: "+str(target))
    np.savez_compressed(target,delta_bar=preparation.delta_bar,p=preparation.p,
                        phi_V=preparation.phi_V,coordinates_m=geometry.coordinates_m,
                        left_nodes=geometry.left_nodes,right_nodes=geometry.right_nodes)
    receipt=dict(preparation.receipt,geometry=geometry.receipt,fields_sha256=_sha(target))
    (output/"effective_dc_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
    return target
