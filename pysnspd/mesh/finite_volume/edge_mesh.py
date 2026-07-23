"""pyTDGL finite-volume EdgeMesh, copied locally for pySNSPD.

Source compatibility target:
    loganbvh/py-tdgl, ``tdgl/finite_volume/edge_mesh.py``
    MIT License, Copyright (c) 2022-2026 Logan Bishop-Van Horn.
"""

from __future__ import annotations

from typing import Sequence, Tuple

import h5py
import numpy as np

from .util import get_dual_edge_lengths, get_edges


class EdgeMesh:
    """A mesh composed of the edges in a triangular mesh."""

    def __init__(
        self,
        centers: Sequence[Tuple[float, float]],
        edges: Sequence[Tuple[int, int]],
        boundary_edge_indices: Sequence[int],
        directions: Sequence[Tuple[float, float]],
        edge_lengths: Sequence[float],
        dual_edge_lengths: Sequence[float],
    ):
        self.centers = np.asarray(centers)
        self.edges = np.asarray(edges)
        self.boundary_edge_indices = np.asarray(boundary_edge_indices, dtype=np.int64)
        self.directions = np.asarray(directions)
        self.normalized_directions = (
            self.directions / np.linalg.norm(self.directions, axis=1)[:, np.newaxis]
        )
        self.edge_lengths = np.asarray(edge_lengths)
        self.dual_edge_lengths = np.asarray(dual_edge_lengths)


    @staticmethod
    def from_mesh(
        sites: np.ndarray,
        elements: np.ndarray,
        dual_sites: np.ndarray,
    ) -> "EdgeMesh":
        """Create edge mesh from mesh."""

        edges, is_boundary = get_edges(elements)
        boundary_edge_indices = np.where(is_boundary)[0]
        edge_coords = sites[edges]
        edge_centers = edge_coords.mean(axis=1)
        directions = np.diff(edge_coords, axis=1).squeeze()
        edge_lengths = np.linalg.norm(directions, axis=1)
        dual_edge_lengths = get_dual_edge_lengths(
            edge_centers,
            elements,
            dual_sites,
            edges,
            len(sites),
        )
        return EdgeMesh(
            edge_centers,
            edges,
            boundary_edge_indices,
            directions,
            edge_lengths,
            dual_edge_lengths,
        )
