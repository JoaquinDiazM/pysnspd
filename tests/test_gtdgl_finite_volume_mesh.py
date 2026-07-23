"""Finite-volume mesh infrastructure consolidated under ``pysnspd.mesh``.

These tests follow the pyTDGL-style finite-volume split:

    geometry.close_curve
    finite_volume.util.get_edges
    finite_volume.Mesh.from_triangulation

The local pySNSPD package keeps coordinates in the caller units, but the
finite-volume topology is pyTDGL-like.
"""

from __future__ import annotations

import numpy as np

from pysnspd.mesh.delaunay import triangle_areas
from pysnspd.mesh.finite_volume.mesh import Mesh
from pysnspd.mesh.finite_volume.util import get_edges
from pysnspd.mesh.geometry import close_curve


def _square_with_center_mesh() -> tuple[np.ndarray, np.ndarray]:
    sites = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
            [0.5, 0.5],
        ],
        dtype=float,
    )

    elements = np.array(
        [
            [0, 1, 4],
            [1, 2, 4],
            [2, 3, 4],
            [3, 0, 4],
        ],
        dtype=np.int64,
    )

    return sites, elements


def test_finite_volume_mesh_from_triangulation_has_edges_and_areas():
    sites, elements = _square_with_center_mesh()

    mesh = Mesh.from_triangulation(sites, elements)

    assert mesh.sites.shape == (5, 2)
    assert mesh.elements.shape == (4, 3)
    assert mesh.edge_mesh is not None
    assert mesh.edge_mesh.edges.shape[1] == 2
    assert mesh.areas.shape == (5,)
    assert np.all(np.isfinite(mesh.areas))
    assert np.all(mesh.areas >= 0.0)
    assert len(mesh.boundary_indices) == 4


def test_finite_volume_utilities_are_consistent():
    sites, elements = _square_with_center_mesh()

    edges, is_boundary = get_edges(elements)

    assert edges.shape[1] == 2
    assert np.count_nonzero(is_boundary) == 4

    areas = triangle_areas(sites, elements)

    assert areas.shape == (4,)
    assert np.allclose(np.abs(areas), 0.25)

    closed = close_curve(sites)

    assert closed.shape == (sites.shape[0] + 1, 2)
    assert np.allclose(closed[0], closed[-1])


def test_mesh_boundary_indices_are_square_corners():
    sites, elements = _square_with_center_mesh()

    mesh = Mesh.from_triangulation(sites, elements)

    assert set(mesh.boundary_indices.tolist()) == {0, 1, 2, 3}
    assert 4 not in set(mesh.boundary_indices.tolist())
