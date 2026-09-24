"""Pass existing SI Delaunay/Voronoi geometry to the thermal graph unchanged.

No meshing, closure, boundary choice, or production operator is introduced here.
Fixed nodes are explicit: geometric side-wall nodes need not be contacts.
"""
import numpy as np
from .thermal_spatial_usadel import ThermalGraph


def _graph(areas, edges, lengths, dual_lengths, coordinates_m, ell0_m,
           fixed_nodes, origin_m):
    if not np.isscalar(ell0_m) or not np.isfinite(ell0_m) or ell0_m <= 0:
        raise ValueError('ell0_m must be one finite positive length in meters')
    lengths = np.asarray(lengths, float)
    dual_lengths = np.asarray(dual_lengths, float)
    if (lengths.ndim != 1 or dual_lengths.shape != lengths.shape
            or len(lengths) != len(edges) or np.any(~np.isfinite(lengths))
            or np.any(lengths <= 0) or np.any(~np.isfinite(dual_lengths))
            or np.any(dual_lengths < 0)):
        raise ValueError('Matching positive primal and nonnegative dual SI lengths required')
    origin = np.zeros(2) if origin_m is None else np.asarray(origin_m, float)
    if origin.shape != (2,) or np.any(~np.isfinite(origin)):
        raise ValueError('origin_m must contain two finite SI coordinates')
    return ThermalGraph(np.asarray(areas, float)/ell0_m**2, edges,
                        dual_lengths/lengths,
                        (np.asarray(coordinates_m, float)-origin)/ell0_m,
                        fixed_nodes)


def graph_from_mesh(mesh, ell0_m, *, fixed_nodes, origin_m=None):
    """Use a complete pyTDGL-compatible Mesh and its own dual measures.

    Mesh coordinates, areas and lengths must be in SI meters. The origin is a
    coordinate translation only; it does not change areas or link orientation.
    """
    if mesh.areas is None or mesh.edge_mesh is None:
        raise ValueError('Mesh must include its areas and EdgeMesh')
    edge_mesh = mesh.edge_mesh
    return _graph(mesh.areas, edge_mesh.edges, edge_mesh.edge_lengths,
                  edge_mesh.dual_edge_lengths, mesh.sites, ell0_m,
                  fixed_nodes, origin_m)


def graph_from_fv_operators(operators, coordinates_m, ell0_m, *, fixed_nodes,
                          origin_m=None):
    """Use one production FVOperators object, never mixing geometric families."""
    return _graph(operators.node_area_m2, operators.edges,
                  operators.edge_length_m, operators.dual_face_length_m,
                  coordinates_m, ell0_m, fixed_nodes, origin_m)
