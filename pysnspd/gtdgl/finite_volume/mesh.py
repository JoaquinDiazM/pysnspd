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

from pysnspd.gtdgl.geometry import close_curve

from .edge_mesh import EdgeMesh
from .util import (
    compute_voronoi_polygon_areas,
    convex_polygon_centroid,
    generate_voronoi_vertices,
    get_edges,
    get_voronoi_polygon_indices,
    triangle_areas,
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

    @property
    def x(self) -> np.ndarray:
        """The x-coordinates of the mesh sites."""

        return self.sites[:, 0]

    @property
    def y(self) -> np.ndarray:
        """The y-coordinates of the mesh sites."""

        return self.sites[:, 1]

    @property
    def center_of_mass(self) -> Tuple[float, float]:
        """The ``(x, y)`` coordinates of the center of mass of the mesh."""

        if self._center_of_mass is None:
            sites = self.sites
            triangles = self.elements
            tri_areas = triangle_areas(sites, triangles)
            tri_centroids = sites[triangles].mean(axis=1)
            com = np.average(tri_centroids, axis=0, weights=tri_areas)
            self._center_of_mass = tuple(com)
        return self._center_of_mass

    def closest_site(self, xy: Tuple[float, float]) -> int:
        """Returns the index of the mesh site closest to ``(x, y)``."""

        return np.argmin(np.linalg.norm(self.sites - np.atleast_2d(xy), axis=1))

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

    def get_quantity_on_site(
        self,
        quantity_on_edge: np.ndarray,
        vector: bool = True,
        use_cupy: bool = False,
    ) -> np.ndarray:
        """Compute an edge quantity averaged over all edges connecting to each site."""

        normalized_directions = self.edge_mesh.normalized_directions
        edges = self.edge_mesh.edges
        if use_cupy:
            xp = cupy
            normalized_directions = xp.asarray(normalized_directions)
            edges = xp.asarray(edges)
        else:
            xp = np
        if vector:
            flux_x = quantity_on_edge * normalized_directions[:, 0]
            flux_y = quantity_on_edge * normalized_directions[:, 1]
        else:
            flux_x = flux_y = quantity_on_edge

        vertices = xp.concatenate([edges[:, 0], edges[:, 1]])
        x_values = xp.concatenate([flux_x, flux_x])
        y_values = xp.concatenate([flux_y, flux_y])
        counts = xp.bincount(vertices)
        x_group_values = xp.bincount(vertices, weights=x_values) / counts
        y_group_values = xp.bincount(vertices, weights=y_values) / counts
        vector_val = xp.array([x_group_values, y_group_values]).T / 2
        if vector:
            return vector_val
        return vector_val[:, 0]

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

    def plot(
        self,
        ax: Union[plt.Axes, None] = None,
        show_sites: bool = True,
        show_edges: bool = False,
        show_dual_edges: bool = True,
        show_voronoi_centroids: bool = False,
        site_color: Union[str, Sequence[float], None] = None,
        edge_color: Union[str, Sequence[float], None] = "k",
        centroid_color: Union[str, Sequence[float], None] = None,
        dual_edge_color: Union[str, Sequence[float], None] = "k",
        linewidth: float = 0.75,
        linestyle: str = "-",
        marker: str = ".",
    ) -> plt.Axes:
        """Plot the mesh."""

        if ax is None:
            _, ax = plt.subplots()
        ax.set_aspect("equal")
        x, y = self.sites.T
        tri = self.elements
        if show_edges:
            ax.triplot(x, y, tri, color=edge_color, ls=linestyle, lw=linewidth)
        if show_dual_edges:
            for poly in self.voronoi_polygons:
                ax.plot(
                    *close_curve(poly).T,
                    color=dual_edge_color,
                    ls=linestyle,
                    lw=linewidth,
                )
        if show_sites:
            ax.plot(x, y, marker=marker, ls="", color=site_color)
        if show_voronoi_centroids:
            centroids = [convex_polygon_centroid(p) for p in self.voronoi_polygons]
            ax.plot(*np.array(centroids).T, marker=marker, ls="", color=centroid_color)
        return ax

    def to_hdf5(self, h5group: h5py.Group, compress: bool = False) -> None:
        """Save the mesh to a :class:`h5py.Group`."""

        h5group["sites"] = self.sites
        h5group["elements"] = self.elements
        if not compress:
            h5group["boundary_indices"] = self.boundary_indices
            h5group["areas"] = self.areas
            self.edge_mesh.to_hdf5(h5group.create_group("edge_mesh"))
            if self.dual_sites is not None:
                h5group["dual_sites"] = self.dual_sites
            split_indices = np.cumsum(
                [len(polygon) for polygon in self.voronoi_polygons[:-1]]
            )
            polygons_flat = np.concatenate(self.voronoi_polygons, axis=0)
            h5group["voronoi_polygons_flat"] = polygons_flat
            h5group["voronoi_split_indices"] = split_indices

    @staticmethod
    def from_hdf5(h5group: h5py.Group) -> "Mesh":
        """Load a mesh from an HDF5 file."""

        if not ("sites" in h5group and "elements" in h5group):
            raise IOError("Could not load mesh due to missing data.")
        if Mesh.is_restorable(h5group):
            polygons_flat = np.array(h5group["voronoi_polygons_flat"])
            voronoi_indices = np.array(h5group["voronoi_split_indices"])
            voronoi_polygons = np.split(polygons_flat, voronoi_indices)
            return Mesh(
                sites=np.array(h5group["sites"]),
                elements=np.array(h5group["elements"], dtype=np.int64),
                boundary_indices=np.array(h5group["boundary_indices"], dtype=np.int64),
                areas=np.array(h5group["areas"]),
                dual_sites=np.array(h5group["dual_sites"]),
                voronoi_polygons=voronoi_polygons,
                edge_mesh=EdgeMesh.from_hdf5(h5group["edge_mesh"]),
            )
        return Mesh.from_triangulation(
            sites=np.array(h5group["sites"]).squeeze(),
            elements=np.array(h5group["elements"]),
        )

    @staticmethod
    def is_restorable(h5group: h5py.Group) -> bool:
        """Return whether a mesh can be restored from the given HDF5 group."""

        return (
            "sites" in h5group
            and "elements" in h5group
            and "boundary_indices" in h5group
            and "areas" in h5group
            and "edge_mesh" in h5group
            and "dual_sites" in h5group
            and "voronoi_polygons_flat" in h5group
            and "voronoi_split_indices" in h5group
        )
