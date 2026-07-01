"""Small pyTDGL-compatible geometry/device subset for pySNSPD PRE meshing.

The classes here implement only the pyTDGL API subset needed to generate a
rectangular PRE mesh through the same logical path as pyTDGL:

    geometry.box -> Polygon -> Device.make_mesh -> generate_mesh
    -> Mesh.from_triangulation

Coordinates are still SI meters in pySNSPD.  To keep the mesh sites in meters,
the compatibility Device uses a default coherence length of 1.0 when called by
``pysnspd.mesh.pytdgl_like``.

Source compatibility target:
    tdgl.device.Layer, tdgl.device.Polygon and tdgl.device.Device mesh helpers,
    pyTDGL, MIT License, Copyright (c) 2022-2026 Logan Bishop-Van Horn.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
from matplotlib import path
from shapely import affinity
from shapely import geometry as geo
from shapely.validation import explain_validity

from pysnspd.gtdgl.finite_volume.mesh import Mesh
from pysnspd.gtdgl.finite_volume.util import get_oriented_boundary
from pysnspd.gtdgl.geometry import close_curve
from pysnspd.gtdgl.finite_volume.meshing import generate_mesh

logger = logging.getLogger(__name__)

PolygonType = Union[
    "Polygon",
    np.ndarray,
    geo.linestring.LineString,
    geo.polygon.LinearRing,
    geo.polygon.Polygon,
]


@dataclass
class Layer:
    """A superconducting thin film, matching the pyTDGL layer fields."""

    london_lambda: float
    coherence_length: float
    thickness: float
    conductivity: Union[float, None] = None
    u: float = 5.79
    gamma: float = 10.0
    z0: float = 0.0

    @property
    def Lambda(self) -> float:
        """Effective magnetic penetration depth, lambda^2 / d."""

        return self.london_lambda**2 / self.thickness

    def copy(self) -> "Layer":
        """Create a deep copy of the layer."""

        return Layer(
            london_lambda=self.london_lambda,
            coherence_length=self.coherence_length,
            thickness=self.thickness,
            conductivity=self.conductivity,
            u=self.u,
            gamma=self.gamma,
            z0=self.z0,
        )


class Polygon:
    """A polygonal region located in a pyTDGL-like layer."""

    def __init__(
        self,
        name: Union[str, None] = None,
        *,
        points: PolygonType,
        mesh: bool = True,
    ):
        self.name = name
        self.points = points
        self.mesh = mesh

    @property
    def points(self) -> np.ndarray:
        """Counter-clockwise-oriented polygon vertices."""

        return self._points

    @points.setter
    def points(self, points) -> None:
        geom_types = (
            geo.linestring.LineString,
            geo.polygon.LinearRing,
            geo.polygon.Polygon,
        )
        if isinstance(points, Polygon):
            points = points.points
        if not isinstance(points, geom_types):
            points = np.asarray(points)
        points = geo.polygon.Polygon(points)
        points = geo.polygon.orient(points)
        if points.interiors:
            raise ValueError("Expected a simply-connected polygon.")
        if not points.is_valid:
            reason = explain_validity(points)
            raise ValueError(
                "The given points do not define a valid polygon for the following "
                f"reason: {reason}."
            )
        points = close_curve(np.array(points.exterior.coords))
        if points.ndim != 2 or points.shape[-1] != 2:
            raise ValueError(f"Expected shape (n, 2), but got {points.shape}.")
        self._points = points

    @property
    def is_valid(self) -> bool:
        """True if the polygon has a name and valid geometry."""

        polygon = self.polygon
        return self.name is not None and polygon.is_valid and not polygon.interiors

    @property
    def area(self) -> float:
        """The area of the polygon."""

        return self.polygon.area

    @property
    def bbox(self) -> Tuple[Tuple[float, float], Tuple[float, float]]:
        """Lower-left and upper-right corners of the bounding box."""

        minx, miny, maxx, maxy = self.polygon.bounds
        return (minx, miny), (maxx, maxy)

    @property
    def extents(self) -> Tuple[float, float]:
        """Total x and y extents of the polygon."""

        minx, miny, maxx, maxy = self.polygon.bounds
        return (maxx - minx), (maxy - miny)

    @property
    def polygon(self) -> geo.polygon.Polygon:
        """A shapely Polygon representation."""

        return geo.polygon.Polygon(self.points)

    @property
    def path(self) -> path.Path:
        """A matplotlib Path representation."""

        return path.Path(self.points, closed=True)

    def contains_points(
        self,
        points: np.ndarray,
        index: bool = False,
        radius: float = 0,
    ) -> Union[bool, np.ndarray]:
        """Determines whether points lie within the polygon."""

        bool_array = self.path.contains_points(np.atleast_2d(points), radius=radius)
        if index:
            return np.where(bool_array)[0]
        return bool_array

    def on_boundary(
        self,
        points: np.ndarray,
        radius: float = 1e-3,
        index: bool = False,
    ):
        """Determines whether points lie within a radius of the boundary."""

        points = np.atleast_2d(points)
        p = self.path
        in_outer = p.contains_points(points, radius=radius)
        in_inner = p.contains_points(points, radius=-radius)
        boundary = np.logical_and(in_outer, ~in_inner)
        if index:
            return np.where(boundary)[0]
        return boundary

    def make_mesh(
        self,
        min_points: Union[int, None] = None,
        smooth: int = 0,
        **meshpy_kwargs,
    ) -> Mesh:
        """Return a finite-volume mesh covering this polygon."""

        points, triangles = generate_mesh(
            self.points,
            min_points=min_points,
            convex_hull=False,
            **meshpy_kwargs,
        )
        if smooth:
            mesh = Mesh.from_triangulation(
                points,
                triangles,
                create_submesh=False,
            ).smooth(smooth)
        else:
            mesh = Mesh.from_triangulation(points, triangles)
        logger.debug(
            f"Finished generating mesh with {len(mesh.sites)} points and "
            f"{len(mesh.elements)} triangles."
        )
        return mesh

    def rotate(
        self,
        degrees: float,
        origin: Union[str, Tuple[float, float]] = (0.0, 0.0),
        inplace: bool = False,
    ) -> "Polygon":
        """Rotate the polygon counterclockwise by a given angle."""

        polygon = self if inplace else self.copy()
        polygon.points = affinity.rotate(
            self.polygon,
            degrees,
            origin=origin,
            use_radians=False,
        )
        return polygon

    def translate(
        self,
        dx: float = 0.0,
        dy: float = 0.0,
        inplace: bool = False,
    ) -> "Polygon":
        """Translate the polygon."""

        polygon = self if inplace else self.copy()
        polygon.points = affinity.translate(self.polygon, xoff=dx, yoff=dy)
        return polygon

    def copy(self) -> "Polygon":
        """Return a deep copy of the polygon."""

        return Polygon(self.name, points=self.points.copy(), mesh=self.mesh)


class Device:
    """Minimal pyTDGL-compatible device used for PRE meshing."""

    def __init__(
        self,
        name: str,
        *,
        layer: Layer,
        film: Polygon,
        holes: Union[List[Polygon], None] = None,
        terminals: Union[List[Polygon], None] = None,
        probe_points: Optional[Sequence[Tuple[float, float]]] = None,
        length_units: str = "m",
    ):
        self.name = name
        self.layer = layer
        self.film = film
        self.holes = holes or []
        self.terminals = tuple(terminals or [])

        terminal_names = set()
        for terminal in self.terminals:
            terminal.mesh = False
            if terminal.name is None or terminal.name in terminal_names:
                raise ValueError("All terminals must have a unique name")
            terminal_names.add(terminal.name)

        for polygon in [self.film] + self.holes:
            if not polygon.is_valid:
                raise ValueError(f"Invalid Polygon: {polygon!r}.")
        if len(self.holes) != len(set(hole.name for hole in self.holes)):
            raise ValueError("All holes must have a unique name.")

        if probe_points is not None:
            probe_points = np.asarray(probe_points).squeeze()
            if probe_points.ndim != 2 or probe_points.shape[1] != 2:
                raise ValueError(
                    f"Probe points must have shape (n, 2), got {probe_points.shape}."
                )
            if not self.contains_points(probe_points).all():
                raise ValueError("All probe points must lie within the film.")
        self.probe_points = probe_points
        self._length_units = length_units
        self.mesh: Optional[Mesh] = None

    @property
    def length_units(self) -> str:
        """Length units used for the device geometry."""

        return self._length_units

    @property
    def polygons(self) -> Tuple[Polygon, ...]:
        """Tuple of all polygons in the device."""

        return (self.film,) + tuple(self.holes) + self.terminals

    @property
    def points(self) -> Union[np.ndarray, None]:
        """The mesh vertex coordinates in length units."""

        if self.mesh is None:
            return None
        return self.mesh.sites * self.layer.coherence_length

    @property
    def triangles(self) -> Union[np.ndarray, None]:
        """Mesh triangle indices."""

        if self.mesh is None:
            return None
        return self.mesh.elements

    def boundary_sites(self) -> Union[Dict[str, np.ndarray], None]:
        """Return ``{polygon_name: oriented_boundary_indices}``."""

        if self.mesh is None:
            return None
        polygons = [self.film] + list(self.holes)
        points = self.points
        edge_mesh = self.mesh.edge_mesh
        boundary_edges = edge_mesh.edges[edge_mesh.boundary_edge_indices]
        boundary = {}
        for polygon in polygons:
            on_boundary = np.logical_and(
                polygon.on_boundary(points[boundary_edges[:, 0]], radius=1e-6),
                polygon.on_boundary(points[boundary_edges[:, 1]], radius=1e-6),
            )
            boundary_sites = get_oriented_boundary(points, boundary_edges[on_boundary])
            assert len(boundary_sites) == 1, len(boundary_sites)
            boundary[polygon.name] = boundary_sites[0]
        return boundary

    def contains_points(
        self,
        points: np.ndarray,
        index: bool = False,
        radius: float = 0,
    ) -> np.ndarray:
        """Determines whether points lie within the device."""

        if self.holes:
            holes_mask = np.logical_or.reduce(
                [hole.contains_points(points, radius=-radius) for hole in self.holes]
            )
        else:
            holes_mask = np.zeros(len(np.atleast_2d(points)), dtype=bool)
        mask = self.film.contains_points(points, radius=radius) & ~holes_mask
        if index:
            return np.where(mask)[0]
        return mask

    def copy(self, with_mesh: bool = True) -> "Device":
        """Copy this Device to create a new one."""

        holes = [hole.copy() for hole in self.holes]
        terminals = [term.copy() for term in self.terminals]
        if self.probe_points is None:
            probe_points = None
        else:
            probe_points = self.probe_points.copy()

        device = Device(
            self.name,
            layer=self.layer.copy(),
            film=self.film.copy(),
            holes=holes,
            terminals=terminals,
            probe_points=probe_points,
            length_units=self.length_units,
        )
        if with_mesh and self.mesh is not None:
            device.mesh = self.mesh
        return device

    def make_mesh(
        self,
        max_edge_length: Union[float, None] = None,
        min_points: Union[float, None] = None,
        smooth: int = 0,
        **meshpy_kwargs,
    ) -> None:
        """Generate and optimize the triangular mesh."""

        logger.info("Generating mesh...")
        if max_edge_length is None:
            max_edge_length = 1.0 * self.layer.coherence_length
        points, triangles = generate_mesh(
            self.film.points,
            hole_coords=[hole.points for hole in self.holes],
            min_points=min_points,
            max_edge_length=max_edge_length,
            boundary=self.film.points,
            **meshpy_kwargs,
        )
        if smooth:
            logger.info("Smoothing mesh.")
            mesh = Mesh.from_triangulation(
                points,
                triangles,
                create_submesh=False,
            ).smooth(smooth, create_submesh=False)
            points = mesh.sites
            triangles = mesh.elements
        logger.info("Creating Mesh object from triangulation.")
        self._create_dimensionless_mesh(points, triangles)
        logger.info(
            f"Finished generating mesh with {len(points)} points and "
            f"{len(triangles)} triangles."
        )

    def _create_dimensionless_mesh(
        self,
        points: np.ndarray,
        triangles: np.ndarray,
    ) -> Mesh:
        """Create the dimensionless mesh."""

        self.mesh = Mesh.from_triangulation(
            points / self.layer.coherence_length,
            triangles,
            create_submesh=True,
        )
        return self.mesh
