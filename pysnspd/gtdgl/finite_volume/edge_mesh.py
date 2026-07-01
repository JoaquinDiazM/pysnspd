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

    @property
    def x(self) -> np.ndarray:
        """The x-coordinates of the edge centers."""

        return self.centers[:, 0]

    @property
    def y(self) -> np.ndarray:
        """The y-coordinates of the edge centers."""

        return self.centers[:, 1]

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

    def to_hdf5(self, h5group: h5py.Group) -> None:
        """Save the data to a HDF5 file."""

        h5group["centers"] = self.centers
        h5group["edges"] = self.edges
        h5group["boundary_edge_indices"] = self.boundary_edge_indices
        h5group["directions"] = self.directions
        h5group["edge_lengths"] = self.edge_lengths
        h5group["dual_edge_lengths"] = self.dual_edge_lengths

    @classmethod
    def from_hdf5(cls, h5group: h5py.Group) -> "EdgeMesh":
        """Load edge mesh from file."""

        if not (
            "centers" in h5group
            and "edges" in h5group
            and "boundary_edge_indices" in h5group
            and "directions" in h5group
            and "edge_lengths" in h5group
            and "dual_edge_lengths" in h5group
        ):
            raise IOError("Could not load edge mesh due to missing data.")
        return EdgeMesh(
            centers=np.array(h5group["centers"]),
            edges=np.array(h5group["edges"], dtype=np.int64),
            boundary_edge_indices=np.array(h5group["boundary_edge_indices"], np.int64),
            directions=np.array(h5group["directions"]),
            edge_lengths=np.array(h5group["edge_lengths"]),
            dual_edge_lengths=np.array(h5group["dual_edge_lengths"]),
        )
