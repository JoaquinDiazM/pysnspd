"""Build the production-family dual mesh, with no spectral or time simulation."""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from pysnspd.mesh.pytdgl_like import PyTDGLLikeMeshParameters, generate_rectangular_pytdgl_fvm_mesh_from_parameters
from pysnspd.gtdgl.tdgl_operators import build_laplacian
from pysnspd.experimental.thermal_mesh_adapter import graph_from_mesh
from pysnspd.experimental.energy_catalog import HBAR_J_S, K_B_J_K


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-root', required=True, type=Path)
    args = parser.parse_args(); args.output_root.mkdir(parents=True, exist_ok=True)
    target = args.output_root/'mesh.npz'
    if target.exists():
        raise FileExistsError('Existing mesh is preserved: '+str(target))
    start = time.perf_counter()
    params = PyTDGLLikeMeshParameters(length_m=160e-9, width_m=80e-9,
        target_spacing_m=5e-9, max_edge_length_m=5e-9, seed=12345)
    print(json.dumps({'event':'MESH_START', 'parameters':asdict(params)}), flush=True)
    mesh = generate_rectangular_pytdgl_fvm_mesh_from_parameters(params)
    ell0 = float(np.sqrt(HBAR_J_S*5e-5/(2*K_B_J_K*8.65)))
    origin = np.array([params.length_m/2, 0.])
    # This threshold only identifies exact straight end-contact coordinates.
    tolerance = 64*np.finfo(float).eps*params.length_m
    left = np.flatnonzero(abs(mesh.sites[:, 0]) <= tolerance)
    right = np.flatnonzero(abs(mesh.sites[:, 0]-params.length_m) <= tolerance)
    if not len(left) or not len(right):
        raise RuntimeError('Both end-contact node sets are required')
    fixed = np.unique(np.r_[left, right])
    graph = graph_from_mesh(mesh, ell0, fixed_nodes=fixed, origin_m=origin)
    centered = mesh.sites-origin
    profile = np.cos(np.pi*centered[:, 0]/params.length_m)**2*np.cos(np.pi*centered[:, 1]/params.width_m)**2
    profile[mesh.boundary_indices] = 0.  # Exact prescribed boundary value of this cosine profile.
    tail, head = graph.edges.T
    flux = graph.conductance*(profile[head]-profile[tail])
    numerator = np.zeros(graph.n_nodes)
    np.add.at(numerator, tail, flux); np.add.at(numerator, head, -flux)
    laplacian, _ = build_laplacian(mesh)
    direct = ell0**2*(laplacian@profile)
    relative = float(np.linalg.norm(numerator/graph.area_weights-direct)/max(np.linalg.norm(direct), np.finfo(float).tiny))
    if relative > 1e-12:
        raise RuntimeError('Adapter does not reproduce the production-family Laplacian')
    np.savez_compressed(target, coordinates_m=mesh.sites, coordinates_bar=graph.coordinates_bar,
        triangles=mesh.elements, edges=graph.edges, area_weights=graph.area_weights,
        conductance=graph.conductance, dual_area_m2=mesh.areas,
        edge_length_m=mesh.edge_mesh.edge_lengths, dual_face_length_m=mesh.edge_mesh.dual_edge_lengths,
        geometric_boundary_nodes=mesh.boundary_indices, fixed_nodes=fixed, left_nodes=left,
        right_nodes=right, origin_m=origin, ell0_m=np.array(ell0), smooth_profile=profile)
    files = ['pysnspd/experimental/thermal_mesh_adapter.py', 'sandbox/stage4_core/prepare_dual_mesh.py',
        'pysnspd/mesh/pytdgl_like.py', 'pysnspd/mesh/finite_volume/mesh.py',
        'pysnspd/mesh/finite_volume/edge_mesh.py', 'pysnspd/gtdgl/tdgl_operators.py']
    receipt = {'schema':'pysnspd.stage4.production_dual_mesh_preparation.v1',
        'status':'MESH_PREPARED_NO_SPECTRAL_OR_TIME_SOLVE', 'parameters':asdict(params),
        'geometry_scope':'80 nm by 160 nm provisional control; 2W is not an admitted Korzh domain',
        'seed_scope':'12345 retained in parameters; this meshpy construction path does not explicitly consume a random seed',
        'coordinate_origin_m':origin.tolist(), 'ell0_m':ell0, 'D_m2_s':5e-5, 'Tc_K':8.65,
        'n_nodes':graph.n_nodes, 'n_triangles':len(mesh.elements), 'n_edges':len(graph.edges),
        'n_left_contact':len(left), 'n_right_contact':len(right),
        'n_geometric_boundary':len(mesh.boundary_indices),
        'fixed_node_policy':'Only x=0 and x=L end contacts; side walls remain natural graph boundaries',
        'profile':'cos^2(pi*(x-L/2)/L)*cos^2(pi*y/W), dimensionless amplitude 0..1; not a photon deposit',
        'profile_maximum':float(profile.max()), 'summed_dual_area_m2':float(mesh.areas.sum()),
        'minimum_dual_area_m2':float(mesh.areas.min()), 'minimum_conductance':float(graph.conductance.min()),
        'maximum_edge_length_m':float(mesh.edge_mesh.edge_lengths.max()),
        'laplacian_relative_difference':relative,
        'runtime_seconds':time.perf_counter()-start,
        'mesh_sha256':hashlib.sha256(target.read_bytes()).hexdigest(),
        'source_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in files}}
    (args.output_root/'receipt.json').write_text(json.dumps(receipt, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(receipt), flush=True)


if __name__ == '__main__':
    main()
