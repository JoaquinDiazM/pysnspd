"""pyTDGL finite-volume Mesh, copied locally for pySNSPD.

Source compatibility target:
    loganbvh/py-tdgl, ``tdgl/finite_volume/mesh.py``
    MIT License, Copyright (c) 2022-2026 Logan Bishop-Van Horn.

Only the import paths differ. Coordinates remain in the units provided by the
caller; pySNSPD PRE meshing supplies SI meters.
"""

from __future__ import annotations

from typing import List, Sequence, Tuple, Union

import h5py
import matplotlib.pyplot as plt
import numpy as np

try:
    import cupy  # type: ignore
except ImportError:  # pragma: no cover - optional acceleration dependency
    cupy = None

from pysnspd.mesh.geometry import close_curve

from .edge_mesh import EdgeMesh
from .util import (
    compute_voronoi_polygon_areas,
    generate_voronoi_vertices,
    get_edges,
    get_voronoi_polygon_indices,
)


class Mesh:
    """A triangular mesh of a simply- or multiply-connected polygon."""

    def __init__(
        self,
        sites: Sequence[Tuple[float, float]],
        elements: Sequence[Tuple[int, int, int]],
        boundary_indices: Sequence[int],
        areas: Union[Sequence[float], None] = None,
        dual_sites: Union[Sequence[Tuple[float, float]], None] = None,
        edge_mesh: Union[EdgeMesh, None] = None,
        voronoi_polygons: Union[List[Sequence[Tuple[float, float]]], None] = None,
    ):
        self.sites = np.asarray(sites).squeeze()
        self.elements = np.asarray(elements, dtype=np.int64)
        self.boundary_indices = np.asarray(boundary_indices, dtype=np.int64)
        if areas is not None:
            areas = np.asarray(areas)
        if dual_sites is not None:
            dual_sites = np.asarray(dual_sites)
        self.areas = areas
        self.dual_sites = dual_sites
        self.edge_mesh = edge_mesh
        self.voronoi_polygons = voronoi_polygons
        self._center_of_mass: Union[Tuple[float, float], None] = None


    @staticmethod
    def from_triangulation(
        sites: Sequence[Tuple[float, float]],
        elements: Sequence[Tuple[int, int, int]],
        create_submesh: bool = True,
    ) -> "Mesh":
        """Create a triangular mesh from vertex coordinates and elements."""

        sites = np.asarray(sites).squeeze()
        elements = np.asarray(elements).squeeze()
        if sites.ndim != 2 or sites.shape[1] != 2:
            raise ValueError(
                f"The site coordinates must have shape (n, 2), got {sites.shape!r}"
            )
        if elements.ndim != 2 or elements.shape[1] != 3:
            raise ValueError(
                f"The elements must have shape (m, 3), got {elements.shape!r}."
            )
        boundary_indices = Mesh.find_boundary_indices(elements)
        dual_sites = edge_mesh = polygons = areas = None
        if create_submesh:
            dual_sites = generate_voronoi_vertices(sites, elements)
            edge_mesh = EdgeMesh.from_mesh(sites, elements, dual_sites)
            areas, polygons = Mesh.compute_voronoi_areas_polygons(
                sites, elements, dual_sites, edge_mesh, boundary_indices
            )
        return Mesh(
            sites=sites,
            elements=elements,
            boundary_indices=boundary_indices,
            edge_mesh=edge_mesh,
            voronoi_polygons=polygons,
            dual_sites=dual_sites,
            areas=areas,
        )

    @staticmethod
    def find_boundary_indices(elements: np.ndarray) -> np.ndarray:
        """Find the boundary vertices."""

        edges, is_boundary = get_edges(elements)
        boundary_edges = edges[is_boundary]
        return np.unique(boundary_edges.flatten())

    @staticmethod
    def compute_voronoi_areas_polygons(
        sites: np.ndarray,
        elements: np.ndarray,
        dual_sites: np.ndarray,
        edge_mesh: EdgeMesh,
        boundary_indices: np.ndarray,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Compute Voronoi cell areas and vertices for each site."""

        polygon_indices = get_voronoi_polygon_indices(elements, len(sites))
        areas, voronoi_polygons = compute_voronoi_polygon_areas(
            sites=sites,
            dual_sites=dual_sites,
            boundary=boundary_indices,
            edges=edge_mesh.edges,
            boundary_edge_indices=edge_mesh.boundary_edge_indices,
            polygons=polygon_indices,
        )
        return areas, voronoi_polygons


    def smooth(self, iterations: int, create_submesh: bool = True) -> "Mesh":
        """Perform Laplacian smoothing of the mesh."""

        mesh = self
        elements = mesh.elements
        edges, _ = get_edges(elements)
        n = len(mesh.sites)
        shape = (n, 2)
        boundary = mesh.boundary_indices
        for i in range(iterations):
            sites = mesh.sites
            num_neighbors = np.bincount(edges.ravel(), minlength=shape[0])
            new_sites = np.zeros(shape)
            vals = sites[edges[:, 1]].T
            new_sites += np.array(
                [np.bincount(edges[:, 0], val, minlength=n) for val in vals]
            ).T
            vals = sites[edges[:, 0]].T
            new_sites += np.array(
                [np.bincount(edges[:, 1], val, minlength=n) for val in vals]
            ).T
            new_sites /= num_neighbors[:, np.newaxis]
            new_sites[boundary] = sites[boundary]
            mesh = Mesh.from_triangulation(
                new_sites,
                elements,
                create_submesh=(create_submesh and (i == (iterations - 1))),
            )
        return mesh
