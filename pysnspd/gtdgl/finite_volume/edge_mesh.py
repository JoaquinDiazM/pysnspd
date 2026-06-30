"""pyTDGL-like EdgeMesh container in SI units."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence, Tuple

import numpy as np

from .util import get_dual_edge_lengths, get_edges


@dataclass
class EdgeMesh:
    """A mesh composed of the edges in a triangular mesh.

    The constructor and attribute names mirror ``tdgl.finite_volume.EdgeMesh``.
    """

    centers: np.ndarray
    edges: np.ndarray
    boundary_edge_indices: np.ndarray
    directions: np.ndarray
    edge_lengths: np.ndarray
    dual_edge_lengths: np.ndarray

    def __post_init__(self) -> None:
        self.centers = np.asarray(self.centers, dtype=float)
        self.edges = np.asarray(self.edges, dtype=np.int64)
        self.boundary_edge_indices = np.asarray(self.boundary_edge_indices, dtype=np.int64)
        self.directions = np.asarray(self.directions, dtype=float)
        self.edge_lengths = np.asarray(self.edge_lengths, dtype=float)
        self.dual_edge_lengths = np.asarray(self.dual_edge_lengths, dtype=float)

    @property
    def x(self) -> np.ndarray:
        """The x-coordinates of the edge centers."""
        return self.centers[:, 0]

    @property
    def y(self) -> np.ndarray:
        """The y-coordinates of the edge centers."""
        return self.centers[:, 1]

    @property
    def normalized_directions(self) -> np.ndarray:
        """Unit vectors along each edge."""
        return self.directions / self.edge_lengths[:, None]

    @staticmethod
    def from_mesh(
        sites: Sequence[Tuple[float, float]],
        elements: Sequence[Tuple[int, int, int]],
        dual_sites: Sequence[Tuple[float, float]],
    ) -> "EdgeMesh":
        """Create edge mesh from mesh, following pyTDGL's data model."""
        sites = np.asarray(sites, dtype=float)
        elements = np.asarray(elements, dtype=np.int64)
        dual_sites = np.asarray(dual_sites, dtype=float)
        edges, is_boundary = get_edges(elements)
        centers = sites[edges].mean(axis=1)
        directions = sites[edges[:, 1]] - sites[edges[:, 0]]
        edge_lengths = np.linalg.norm(directions, axis=1)
        boundary_edge_indices = np.where(is_boundary)[0]
        dual_edge_lengths = get_dual_edge_lengths(
            centers,
            elements,
            dual_sites,
            edges,
            len(sites),
        )
        return EdgeMesh(
            centers=centers,
            edges=edges,
            boundary_edge_indices=boundary_edge_indices,
            directions=directions,
            edge_lengths=edge_lengths,
            dual_edge_lengths=dual_edge_lengths,
        )
