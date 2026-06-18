"""
Edge extraction and boundary tagging for pySNSPD meshes.

OE2 scope:
- Extract unique triangle edges.
- Identify boundary edges.
- Tag left/right contacts and top/bottom insulating boundaries.

No finite-volume operator is assembled here yet.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np


@dataclass(frozen=True)
class EdgeData:
    """
    Container for mesh edge connectivity and boundary tags.

    Attributes
    ----------
    edges:
        Array ``(n_edges, 2)`` with sorted node indices.
    edge_triangles:
        Array ``(n_edges, 2)`` containing neighboring triangle indices.
        Boundary edges have one entry equal to ``-1``.
    is_boundary:
        Boolean array marking boundary edges.
    tags:
        String array with values ``interior``, ``left``, ``right``, ``top``,
        ``bottom`` or ``boundary_unknown``.
    midpoints:
        Array ``(n_edges, 2)`` with edge midpoints.
    lengths:
        Array ``(n_edges,)`` with edge lengths in meters.
    """
    edges: np.ndarray
    edge_triangles: np.ndarray
    is_boundary: np.ndarray
    tags: np.ndarray
    midpoints: np.ndarray
    lengths: np.ndarray

    @property
    def n_edges(self) -> int:
        """Return number of unique edges."""
        return int(self.edges.shape[0])

    @property
    def n_boundary_edges(self) -> int:
        """Return number of boundary edges."""
        return int(np.count_nonzero(self.is_boundary))


def extract_unique_edges(triangles: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Extract unique edges and neighboring triangles from a triangle table.

    Parameters
    ----------
    triangles:
        Array of shape ``(n_triangles, 3)``.

    Returns
    -------
    edges, edge_triangles:
        ``edges`` has shape ``(n_edges, 2)``.
        ``edge_triangles`` has shape ``(n_edges, 2)``.
    """
    edge_map: dict[tuple[int, int], list[int]] = {}

    for t_idx, tri in enumerate(np.asarray(triangles, dtype=np.int64)):
        local_edges = [
            (int(tri[0]), int(tri[1])),
            (int(tri[1]), int(tri[2])),
            (int(tri[2]), int(tri[0])),
        ]

        for a, b in local_edges:
            key = (a, b) if a < b else (b, a)
            edge_map.setdefault(key, []).append(int(t_idx))

    edges = np.array(sorted(edge_map.keys()), dtype=np.int64)
    edge_triangles = -np.ones((edges.shape[0], 2), dtype=np.int64)

    for i, key in enumerate(edges):
        neighbors = edge_map[(int(key[0]), int(key[1]))]
        if len(neighbors) > 2:
            raise ValueError(f"Non-manifold edge {tuple(key)} belongs to {len(neighbors)} triangles.")
        edge_triangles[i, : len(neighbors)] = neighbors

    return edges, edge_triangles


def build_edge_data(
    nodes: np.ndarray,
    triangles: np.ndarray,
    *,
    length_m: float,
    width_m: float,
    boundary_tol_m: float | None = None,
) -> EdgeData:
    """
    Build edge connectivity and assign boundary tags.

    Boundary tags are based on the midpoint position and require the edge to be
    a true boundary edge.
    """
    nodes = np.asarray(nodes, dtype=float)
    triangles = np.asarray(triangles, dtype=np.int64)

    edges, edge_triangles = extract_unique_edges(triangles)

    p0 = nodes[edges[:, 0]]
    p1 = nodes[edges[:, 1]]
    midpoints = 0.5 * (p0 + p1)
    lengths = np.linalg.norm(p1 - p0, axis=1)

    is_boundary = np.any(edge_triangles < 0, axis=1)

    if boundary_tol_m is None:
        positive_lengths = lengths[lengths > 0.0]
        scale = float(np.min(positive_lengths)) if positive_lengths.size else 1.0
        boundary_tol_m = max(1.0e-15, 1.0e-6 * scale)

    tags = np.full(edges.shape[0], "interior", dtype="<U32")

    xmin = 0.0
    xmax = float(length_m)
    ymin = -0.5 * float(width_m)
    ymax = 0.5 * float(width_m)

    left = is_boundary & np.isclose(midpoints[:, 0], xmin, atol=boundary_tol_m, rtol=0.0)
    right = is_boundary & np.isclose(midpoints[:, 0], xmax, atol=boundary_tol_m, rtol=0.0)
    bottom = is_boundary & np.isclose(midpoints[:, 1], ymin, atol=boundary_tol_m, rtol=0.0)
    top = is_boundary & np.isclose(midpoints[:, 1], ymax, atol=boundary_tol_m, rtol=0.0)

    tags[left] = "left"
    tags[right] = "right"
    tags[bottom] = "bottom"
    tags[top] = "top"

    unknown = is_boundary & (tags == "interior")
    tags[unknown] = "boundary_unknown"

    return EdgeData(
        edges=edges,
        edge_triangles=edge_triangles,
        is_boundary=is_boundary,
        tags=tags,
        midpoints=midpoints,
        lengths=lengths,
    )


def edge_summary(edge_data: EdgeData) -> dict[str, Any]:
    """
    Build a compact summary dictionary for edge diagnostics.
    """
    tags, counts = np.unique(edge_data.tags, return_counts=True)
    tag_counts = {str(tag): int(count) for tag, count in zip(tags, counts)}

    return {
        "n_edges": edge_data.n_edges,
        "n_boundary_edges": edge_data.n_boundary_edges,
        "n_interior_edges": int(edge_data.n_edges - edge_data.n_boundary_edges),
        "edge_length_min_m": float(np.min(edge_data.lengths)),
        "edge_length_max_m": float(np.max(edge_data.lengths)),
        "edge_length_mean_m": float(np.mean(edge_data.lengths)),
        "tag_counts": tag_counts,
    }


def save_edges_npz(edge_data: EdgeData, path: str | Path) -> Path:
    """
    Save edge arrays to compressed ``.npz``.
    """
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)

    np.savez_compressed(
        output,
        edges=edge_data.edges,
        edge_triangles=edge_data.edge_triangles,
        is_boundary=edge_data.is_boundary,
        tags=edge_data.tags,
        midpoints=edge_data.midpoints,
        lengths=edge_data.lengths,
    )

    return output


def load_edges_npz(path: str | Path) -> EdgeData:
    """
    Load edge arrays saved by :func:`save_edges_npz`.
    """
    source = Path(path)

    with np.load(source) as data:
        return EdgeData(
            edges=np.asarray(data["edges"], dtype=np.int64),
            edge_triangles=np.asarray(data["edge_triangles"], dtype=np.int64),
            is_boundary=np.asarray(data["is_boundary"], dtype=bool),
            tags=np.asarray(data["tags"]).astype("<U32"),
            midpoints=np.asarray(data["midpoints"], dtype=float),
            lengths=np.asarray(data["lengths"], dtype=float),
        )


def assert_edge_data_consistent(edge_data: EdgeData) -> None:
    """
    Raise ValueError if the edge table has obvious consistency problems.
    """
    if edge_data.edges.ndim != 2 or edge_data.edges.shape[1] != 2:
        raise ValueError("edges must have shape (n_edges, 2).")

    if edge_data.edge_triangles.shape != (edge_data.edges.shape[0], 2):
        raise ValueError("edge_triangles must have shape (n_edges, 2).")

    if edge_data.is_boundary.shape[0] != edge_data.edges.shape[0]:
        raise ValueError("is_boundary length must match number of edges.")

    if edge_data.tags.shape[0] != edge_data.edges.shape[0]:
        raise ValueError("tags length must match number of edges.")

    if np.any(edge_data.lengths <= 0.0):
        raise ValueError("All edge lengths must be positive.")

    expected_boundary = np.any(edge_data.edge_triangles < 0, axis=1)
    if not np.array_equal(edge_data.is_boundary, expected_boundary):
        raise ValueError("is_boundary does not match edge_triangles.")