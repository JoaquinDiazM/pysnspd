"""Prepare the two additional dual meshes for the pre-photon DC comparison.

The existing 160 x 80 nm, h=5 nm mesh is reused by the campaign, not copied.
This script only builds h=3.5 nm on that rectangle and h=5 nm on a 320 x 80 nm
rectangle. The pyTDGL-family production mesher and graph adapter are unchanged.
Boundary vertices are sampled no farther apart than half the target interior
edge length, following the previously admitted boundary-resampling rule.
No spectral, kinetic or time evolution is performed.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import json
import logging
import os
from pathlib import Path
import sys
import time

for _name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
              'NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[_name]='1'

import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from pysnspd.mesh.pytdgl_like import (
    PyTDGLLikeMeshParameters,generate_rectangular_pytdgl_fvm_mesh_from_parameters,
)
from pysnspd.gtdgl.tdgl_operators import build_laplacian
from pysnspd.experimental.thermal_mesh_adapter import graph_from_mesh
from pysnspd.experimental.bulk_current_reference import physical_scales

DEFAULT_OUTPUT=ROOT/'docs/implementation/stage5/prephoton_dc_20260924/meshes'
CASES=(('l160_h35',160e-9,80e-9,3.5e-9),
       ('l320_h5',320e-9,80e-9,5e-9))
BASE_MESH='docs/implementation/stage4/practical_time_review_20260924/dual_mesh/resampled/mesh.npz'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def boundary_sampling(length,width,target):
    """Smallest box point request meeting half-target spacing on both sides.

    ``geometry.box`` rounds its allocation to each side. Check those exact
    rounded counts, rather than assuming requested points equal final points.
    This is a deterministic geometry rule, not a search over mesh outcomes.
    """
    perimeter=2*(length+width)
    half=.5*target
    requested=max(8,int(np.ceil(perimeter/half)))
    while True:
        nx=round(requested*length/perimeter)
        ny=round(requested*width/perimeter)
        if nx>1 and ny>1 and max(length/(nx-1),width/(ny-1))<=half*(1+1e-14):
            return requested,nx,ny
        requested+=1


def create_mesh(output,name,length,width,spacing):
    target=output/name
    if target.exists():
        raise FileExistsError('Existing mesh output is preserved: '+str(target))
    requested,nx,ny=boundary_sampling(length,width,spacing)
    params=PyTDGLLikeMeshParameters(length_m=length,width_m=width,
        target_spacing_m=spacing,max_edge_length_m=spacing,seed=12345,
        boundary_points=requested,terminal_contact_mode='normal_left_right')
    start=time.perf_counter();warnings=[]
    class Capture(logging.Handler):
        def emit(self,record):
            if record.levelno>=logging.WARNING:
                warnings.append(record.getMessage())
    handler=Capture();logger=logging.getLogger('tdgl.finite_volume');logger.addHandler(handler)
    print(json.dumps(dict(event='MESH_START',case=name,parameters=asdict(params))),flush=True)
    try:
        mesh=generate_rectangular_pytdgl_fvm_mesh_from_parameters(params)
    finally:
        logger.removeHandler(handler)
    ell=physical_scales(Tc_K=8.65,diffusion_m2_s=5e-5,sheet_resistance_ohm=608.)['ell0_m']
    origin=np.array([length/2,0.])
    tolerance=64*np.finfo(float).eps*length
    left=np.flatnonzero(abs(mesh.sites[:,0])<=tolerance)
    right=np.flatnonzero(abs(mesh.sites[:,0]-length)<=tolerance)
    if not len(left) or not len(right):
        raise RuntimeError('Both intrinsic end reservoirs need mesh nodes')
    fixed=np.unique(np.r_[left,right])
    graph=graph_from_mesh(mesh,ell,fixed_nodes=fixed,origin_m=origin)
    xy=mesh.sites-origin
    profile=np.cos(np.pi*xy[:,0]/length)**2*np.cos(np.pi*xy[:,1]/width)**2
    profile[mesh.boundary_indices]=0.
    tail,head=graph.edges.T
    flux=graph.conductance*(profile[head]-profile[tail]);numerator=np.zeros(graph.n_nodes)
    np.add.at(numerator,tail,flux);np.add.at(numerator,head,-flux)
    laplacian,_=build_laplacian(mesh);direct=ell**2*(laplacian@profile)
    relative=float(np.linalg.norm(numerator/graph.area_weights-direct)/np.linalg.norm(direct))
    if relative>1e-12:
        raise RuntimeError('Graph does not reproduce the production-family Laplacian')
    centers=mesh.dual_sites
    outside=((centers[:,0]<-tolerance)|(centers[:,0]>length+tolerance)
             |(abs(centers[:,1])>width/2+tolerance))
    area_error=float(mesh.areas.sum()/(length*width)-1)
    if abs(area_error)>1e-12 or np.any(outside):
        raise RuntimeError('Dual area or circumcenter geometry requires review before use')
    vertices=mesh.sites[mesh.elements]
    side_lengths=np.stack([np.linalg.norm(vertices[:,(i+1)%3]-vertices[:,i],axis=1)
                           for i in range(3)],axis=1)
    cosines=np.empty_like(side_lengths)
    for i in range(3):
        a,b,c=side_lengths[:,i],side_lengths[:,(i+1)%3],side_lengths[:,(i+2)%3]
        cosines[:,i]=np.clip((b*b+c*c-a*a)/(2*b*c),-1.,1.)
    minimum_angle=float(np.degrees(np.arccos(cosines)).min())
    target.mkdir(parents=True)
    np.savez_compressed(target/'mesh.npz',coordinates_m=mesh.sites,
        coordinates_bar=graph.coordinates_bar,triangles=mesh.elements,edges=graph.edges,
        area_weights=graph.area_weights,conductance=graph.conductance,dual_area_m2=mesh.areas,
        edge_length_m=mesh.edge_mesh.edge_lengths,dual_face_length_m=mesh.edge_mesh.dual_edge_lengths,
        geometric_boundary_nodes=mesh.boundary_indices,fixed_nodes=fixed,left_nodes=left,
        right_nodes=right,origin_m=origin,ell0_m=np.array(ell),smooth_profile=profile)
    sources=['sandbox/stage5_prephoton/prepare_meshes.py','pysnspd/mesh/pytdgl_like.py',
        'pysnspd/mesh/finite_volume/mesh.py','pysnspd/mesh/finite_volume/edge_mesh.py',
        'pysnspd/experimental/thermal_mesh_adapter.py','pysnspd/gtdgl/tdgl_operators.py']
    receipt=dict(schema='pysnspd.stage5.prephoton_dual_mesh.v1',
        status='MESH_PREPARED_NO_SPECTRAL_OR_TIME_SOLVE',case=name,parameters=asdict(params),
        seed_scope='12345 retained as metadata; the meshpy construction path does not explicitly consume a random seed',
        material=dict(Tc_K=8.65,diffusion_m2_s=5e-5),ell0_m=ell,
        coordinate_origin_m=origin.tolist(),boundary_rule='Requested boundary spacing <= half the target interior edge length',
        requested_points_per_long_side=nx,requested_points_per_short_side=ny,
        requested_long_side_spacing_m=length/(nx-1),requested_short_side_spacing_m=width/(ny-1),
        terminal_label_scope='normal_left_right is inherited mesh metadata only; the DC solver imposes superconducting current-carrying bulk spectra at these cuts, not normal-metal contacts',
        fixed_node_policy='Only x=0 and x=L cuts are fixed; transverse side walls retain natural graph boundaries',
        n_nodes=graph.n_nodes,n_edges=len(graph.edges),n_triangles=len(mesh.elements),
        n_left_contact=len(left),n_right_contact=len(right),n_geometric_boundary=len(mesh.boundary_indices),
        bounds_m=dict(x=[float(mesh.sites[:,0].min()),float(mesh.sites[:,0].max())],
                      y=[float(mesh.sites[:,1].min()),float(mesh.sites[:,1].max())]),
        area_relative_error=area_error,minimum_dual_area_m2=float(mesh.areas.min()),
        minimum_conductance=float(graph.conductance.min()),
        maximum_edge_length_m=float(mesh.edge_mesh.edge_lengths.max()),
        minimum_triangle_angle_degrees=minimum_angle,
        n_circumcenters_outside_rectangle=int(outside.sum()),mesh_warnings=warnings,
        laplacian_relative_difference=relative,runtime_seconds=time.perf_counter()-start,
        profile='Saved smooth cosine profile is only an operator check; not a photon or an applied perturbation',
        base_mesh_reused_without_copy=BASE_MESH,base_mesh_sha256=sha(ROOT/BASE_MESH),
        mesh_sha256=sha(target/'mesh.npz'),source_sha256={name:sha(ROOT/name) for name in sources})
    (target/'morphology.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n',encoding='utf8')
    print(json.dumps(dict(event='MESH_COMPLETE',**receipt),allow_nan=False),flush=True)
    return receipt


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-root',type=Path,default=DEFAULT_OUTPUT)
    args=parser.parse_args()
    original=sorted(os.sched_getaffinity(0)) if hasattr(os,'sched_getaffinity') else None
    resource=dict(logical_cpus=os.cpu_count(),available_affinity=original,
                  mesh_workers=1,numerical_library_threads=1)
    if original:
        topology=[]
        for cpu in original:
            folder=Path(f'/sys/devices/system/cpu/cpu{cpu}/topology')
            topology.append((folder.joinpath('physical_package_id').read_text().strip(),
                             folder.joinpath('core_id').read_text().strip()))
        resource['physical_cores_in_affinity']=len(set(topology))
        quota=Path('/sys/fs/cgroup/cpu.max')
        if quota.exists():
            resource['cgroup_cpu_max']=quota.read_text().strip()
        else:
            for name in ('cpu.cfs_quota_us','cpu.cfs_period_us'):
                path=Path('/sys/fs/cgroup/cpu')/name
                resource[name]=path.read_text().strip() if path.exists() else None
        memory=Path('/proc/meminfo').read_text().splitlines()
        resource['memory_available']=[line for line in memory if line.startswith('MemAvailable:')]
        os.sched_setaffinity(0,{original[0]})
        resource['admitted_affinity']=[original[0]]
    print(json.dumps(dict(event='MESH_RESOURCES',**resource)),flush=True)
    try:
        for name,*geometry in CASES:
            create_mesh(args.output_root,name,*geometry)
    finally:
        if original:
            os.sched_setaffinity(0,original)


if __name__=='__main__':
    main()
