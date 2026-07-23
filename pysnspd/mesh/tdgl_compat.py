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

from pysnspd.mesh.finite_volume.mesh import Mesh
from pysnspd.mesh.geometry import close_curve
from pysnspd.mesh.finite_volume.meshing import generate_mesh

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
