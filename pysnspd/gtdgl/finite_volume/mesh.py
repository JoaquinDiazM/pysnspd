"""pyTDGL-like finite-volume Mesh in SI units."""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Sequence, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np

from .edge_mesh import EdgeMesh
from .util import (
    close_curve,
    compute_voronoi_polygon_areas,
    convex_polygon_centroid,
    generate_voronoi_vertices,
    get_edges,
    get_voronoi_polygon_indices,
    triangle_areas,
)


def _orient_elements_counterclockwise(sites: np.ndarray, elements: np.ndarray) -> np.ndarray:
    """Return triangle elements with counter-clockwise vertex ordering.

    The directed Voronoi-polygon construction used by pyTDGL assumes a
    consistent orientation of the Delaunay triangles.  ``meshpy.triangle`` may
    return valid triangles whose vertex order is not guaranteed after pySNSPD's
    wrapper/refinement path, so we normalize the ordering immediately when the
    SI Mesh is constructed.
    """
    sites = np.asarray(sites, dtype=float)
    out = np.asarray(elements, dtype=np.int64).copy()
    xy = sites[out]
    signed_area2 = (
        (xy[:, 1, 0] - xy[:, 0, 0]) * (xy[:, 2, 1] - xy[:, 0, 1])
        - (xy[:, 1, 1] - xy[:, 0, 1]) * (xy[:, 2, 0] - xy[:, 0, 0])
    )
    flip = signed_area2 < 0.0
    if np.any(flip):
        out[flip, 1], out[flip, 2] = out[flip, 2].copy(), out[flip, 1].copy()
    return out


@dataclass
class Mesh:
    """A triangular mesh of a simply- or multiply-connected polygon.

    This mirrors ``tdgl.finite_volume.Mesh``.  Coordinates are SI meters.
    """

    sites: np.ndarray
    elements: np.ndarray
    boundary_indices: np.ndarray
    areas: Union[np.ndarray, None] = None
    dual_sites: Union[np.ndarray, None] = None
    edge_mesh: Union[EdgeMesh, None] = None
    voronoi_polygons: Union[List[np.ndarray], None] = None

    def __post_init__(self) -> None:
        self.sites = np.asarray(self.sites, dtype=float).squeeze()
        self.elements = np.asarray(self.elements, dtype=np.int64)
        self.boundary_indices = np.asarray(self.boundary_indices, dtype=np.int64)
        if self.areas is not None:
            self.areas = np.asarray(self.areas, dtype=float)
        if self.dual_sites is not None:
            self.dual_sites = np.asarray(self.dual_sites, dtype=float)
        self._center_of_mass: Tuple[float, float] | None = None

    @property
    def x(self) -> np.ndarray:
        return self.sites[:, 0]

    @property
    def y(self) -> np.ndarray:
        return self.sites[:, 1]

    @property
    def center_of_mass(self) -> Tuple[float, float]:
        if self._center_of_mass is None:
            tri_areas = triangle_areas(self.sites, self.elements)
            tri_centroids = self.sites[self.elements].mean(axis=1)
            self._center_of_mass = tuple(np.average(tri_centroids, axis=0, weights=tri_areas))
        return self._center_of_mass

    def closest_site(self, xy: Tuple[float, float]) -> int:
        return int(np.argmin(np.linalg.norm(self.sites - np.atleast_2d(xy), axis=1)))

    @staticmethod
    def from_triangulation(
        sites: Sequence[Tuple[float, float]],
        elements: Sequence[Tuple[int, int, int]],
        create_submesh: bool = True,
    ) -> "Mesh":
        """Create a triangular mesh from sites and triangle elements."""
        sites = np.asarray(sites, dtype=float).squeeze()
        elements = np.asarray(elements, dtype=np.int64).squeeze()
        if sites.ndim != 2 or sites.shape[1] != 2:
            raise ValueError(f"The site coordinates must have shape (n, 2), got {sites.shape!r}")
        if elements.ndim != 2 or elements.shape[1] != 3:
            raise ValueError(f"The elements must have shape (m, 3), got {elements.shape!r}.")
        elements = _orient_elements_counterclockwise(sites, elements)
        boundary_indices = Mesh.find_boundary_indices(elements)
        dual_sites = edge_mesh = polygons = areas = None
        if create_submesh:
            dual_sites = generate_voronoi_vertices(sites, elements)
            edge_mesh = EdgeMesh.from_mesh(sites, elements, dual_sites)
            areas, polygons = Mesh.compute_voronoi_areas_polygons(
                sites,
                elements,
                dual_sites,
                edge_mesh,
                boundary_indices,
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
        edges, is_boundary = get_edges(elements)
        return np.unique(edges[is_boundary].flatten())

    @staticmethod
    def compute_voronoi_areas_polygons(
        sites: np.ndarray,
        elements: np.ndarray,
        dual_sites: np.ndarray,
        edge_mesh: EdgeMesh,
        boundary_indices: np.ndarray,
    ) -> Tuple[np.ndarray, list[np.ndarray]]:
        polygon_indices = get_voronoi_polygon_indices(elements, len(sites))
        return compute_voronoi_polygon_areas(
            sites=sites,
            dual_sites=dual_sites,
            boundary=boundary_indices,
            edges=edge_mesh.edges,
            boundary_edge_indices=edge_mesh.boundary_edge_indices,
            polygons=polygon_indices,
        )

    def get_quantity_on_site(
        self,
        quantity_on_edge: np.ndarray,
        vector: bool = True,
        use_cupy: bool = False,
    ) -> np.ndarray:
        if use_cupy:
            raise NotImplementedError("CuPy is not used in pySNSPD's SI backend yet.")
        if self.edge_mesh is None:
            raise ValueError("Mesh has no edge_mesh; create it with create_submesh=True.")
        normalized_directions = self.edge_mesh.normalized_directions
        edges = self.edge_mesh.edges
        if vector:
            flux_x = quantity_on_edge * normalized_directions[:, 0]
            flux_y = quantity_on_edge * normalized_directions[:, 1]
        else:
            flux_x = flux_y = quantity_on_edge
        vertices = np.concatenate([edges[:, 0], edges[:, 1]])
        x_values = np.concatenate([flux_x, flux_x])
        y_values = np.concatenate([flux_y, flux_y])
        counts = np.bincount(vertices)
        x_group_values = np.bincount(vertices, weights=x_values) / counts
        y_group_values = np.bincount(vertices, weights=y_values) / counts
        vector_val = np.array([x_group_values, y_group_values]).T / 2
        if vector:
            return vector_val
        return vector_val[:, 0]

    def smooth(self, iterations: int, create_submesh: bool = True) -> "Mesh":
        """Perform Laplacian smoothing while holding boundary sites fixed."""
        mesh = self
        elements = mesh.elements
        edges, _ = get_edges(elements)
        n = len(mesh.sites)
        shape = (n, 2)
        boundary = mesh.boundary_indices
        for i in range(int(iterations)):
            sites = mesh.sites
            num_neighbors = np.bincount(edges.ravel(), minlength=shape[0])
            new_sites = np.zeros(shape)
            vals = sites[edges[:, 1]].T
            new_sites += np.array([np.bincount(edges[:, 0], val, minlength=n) for val in vals]).T
            vals = sites[edges[:, 0]].T
            new_sites += np.array([np.bincount(edges[:, 1], val, minlength=n) for val in vals]).T
            new_sites /= num_neighbors[:, np.newaxis]
            new_sites[boundary] = sites[boundary]
            mesh = Mesh.from_triangulation(
                new_sites,
                elements,
                create_submesh=(create_submesh and (i == (int(iterations) - 1))),
            )
        return mesh

    def plot(
        self,
        ax: Union[plt.Axes, None] = None,
        show_sites: bool = True,
        show_edges: bool = False,
        show_dual_edges: bool = True,
        show_voronoi_centroids: bool = False,
        site_color=None,
        edge_color="k",
        centroid_color=None,
        dual_edge_color="k",
        linewidth: float = 0.75,
        linestyle: str = "-",
        marker: str = ".",
    ) -> plt.Axes:
        if ax is None:
            _, ax = plt.subplots()
        ax.set_aspect("equal")
        x, y = self.sites.T
        tri = self.elements
        if show_edges:
            ax.triplot(x, y, tri, color=edge_color, ls=linestyle, lw=linewidth)
        if show_dual_edges and self.voronoi_polygons is not None:
            for poly in self.voronoi_polygons:
                ax.plot(*close_curve(poly).T, color=dual_edge_color, ls=linestyle, lw=linewidth)
        if show_sites:
            ax.plot(x, y, marker=marker, ls="", color=site_color)
        if show_voronoi_centroids and self.voronoi_polygons is not None:
            centroids = [convex_polygon_centroid(p) for p in self.voronoi_polygons]
            ax.plot(*np.array(centroids).T, marker=marker, ls="", color=centroid_color)
        return ax
