"""One deterministic boundary resampling of the same provisional rectangle.

The first mesh is preserved. Boundary facets are sampled at half the target
interior edge length; no random seed search or physical parameter change.
"""
import hashlib
import json
import logging
from pathlib import Path
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from pysnspd.mesh.pytdgl_like import PyTDGLLikeMeshParameters, generate_rectangular_pytdgl_fvm_mesh_from_parameters
from pysnspd.gtdgl.tdgl_operators import build_laplacian
from pysnspd.experimental.thermal_mesh_adapter import graph_from_mesh


def main():
    out = ROOT/'docs/implementation/stage4/practical_time_review_20260924/dual_mesh/resampled'
    out.mkdir(parents=True, exist_ok=True)
    if (out/'mesh.npz').exists():
        raise FileExistsError('Preserving previous resampling output')
    original = json.loads((out.parent/'receipt.json').read_text())
    p = original['parameters'].copy()
    # box allocates about 1/3 and 1/6 of points to each 2:1 rectangle side.
    # 196 gives 65 samples on x and 33 on y: 64+32 intervals of 2.5 nm.
    p['boundary_points'] = 196
    warnings = []
    class Capture(logging.Handler):
        def emit(self, record):
            if record.levelno >= logging.WARNING:
                warnings.append(record.getMessage())
    logger = logging.getLogger('tdgl.finite_volume'); handler = Capture(); logger.addHandler(handler)
    start = time.perf_counter()
    try:
        mesh = generate_rectangular_pytdgl_fvm_mesh_from_parameters(PyTDGLLikeMeshParameters(**p))
    finally:
        logger.removeHandler(handler)
    ell0 = original['ell0_m']; origin = np.asarray(original['coordinate_origin_m'])
    tolerance = 64*np.finfo(float).eps*p['length_m']
    left = np.flatnonzero(abs(mesh.sites[:, 0]) <= tolerance)
    right = np.flatnonzero(abs(mesh.sites[:, 0]-p['length_m']) <= tolerance)
    fixed = np.unique(np.r_[left, right])
    graph = graph_from_mesh(mesh, ell0, fixed_nodes=fixed, origin_m=origin)
    xy = mesh.sites-origin
    profile = np.cos(np.pi*xy[:, 0]/p['length_m'])**2*np.cos(np.pi*xy[:, 1]/p['width_m'])**2
    profile[mesh.boundary_indices] = 0.
    tail, head = graph.edges.T
    flux = graph.conductance*(profile[head]-profile[tail]); numerator = np.zeros(graph.n_nodes)
    np.add.at(numerator, tail, flux); np.add.at(numerator, head, -flux)
    laplacian, _ = build_laplacian(mesh); direct = ell0**2*(laplacian@profile)
    relative = float(np.linalg.norm(numerator/graph.area_weights-direct)/np.linalg.norm(direct))
    centers = mesh.dual_sites
    outside = ((centers[:, 0] < -tolerance) | (centers[:, 0] > p['length_m']+tolerance)
        | (abs(centers[:, 1]) > p['width_m']/2+tolerance))
    area_error = float(mesh.areas.sum()/(p['length_m']*p['width_m'])-1)
    np.savez_compressed(out/'mesh.npz', coordinates_m=mesh.sites,
        coordinates_bar=graph.coordinates_bar, triangles=mesh.elements, edges=graph.edges,
        area_weights=graph.area_weights, conductance=graph.conductance, dual_area_m2=mesh.areas,
        edge_length_m=mesh.edge_mesh.edge_lengths, dual_face_length_m=mesh.edge_mesh.dual_edge_lengths,
        geometric_boundary_nodes=mesh.boundary_indices, fixed_nodes=fixed, left_nodes=left,
        right_nodes=right, origin_m=origin, ell0_m=np.array(ell0), smooth_profile=profile)
    receipt = dict(schema='pysnspd.stage4.dual_mesh_boundary_resampling.v1',
        status='RESAMPLING_ONLY_REVIEW_GEOMETRY_METRICS', parameters=p,
        unchanged='Material, rectangular boundary, target interior edge length, equations; no spectra solved',
        boundary_rule='Half target interior length: 64 longitudinal and 32 transverse intervals before Triangle refinement',
        source_mesh_sha256=original['mesh_sha256'], n_nodes=graph.n_nodes,
        n_edges=len(graph.edges), n_triangles=len(mesh.elements),
        n_left_contact=len(left), n_right_contact=len(right),
        area_relative_error=area_error, n_circumcenters_outside_rectangle=int(outside.sum()),
        mesh_warnings=warnings, minimum_dual_area_m2=float(mesh.areas.min()),
        minimum_conductance=float(graph.conductance.min()), laplacian_relative_difference=relative,
        runtime_seconds=time.perf_counter()-start,
        mesh_sha256=hashlib.sha256((out/'mesh.npz').read_bytes()).hexdigest(),
        runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (out/'receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt), flush=True)


if __name__ == '__main__':
    main()
