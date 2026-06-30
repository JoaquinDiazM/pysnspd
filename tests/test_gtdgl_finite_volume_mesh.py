"""Finite-volume mesh infrastructure promoted to ``pysnspd.gtdgl``."""
from __future__ import annotations

import numpy as np

from pysnspd.gtdgl.finite_volume import Mesh
from pysnspd.gtdgl.finite_volume.util import close_curve, get_edges, triangle_areas


def test_finite_volume_mesh_from_triangulation_has_edges_and_areas():
    sites = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
        ],
        dtype=float,
    )
    elements = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.int64)
    mesh = Mesh.from_triangulation(sites, elements)
    assert mesh.sites.shape == (4, 2)
    assert mesh.edge_mesh is not None
    assert mesh.edge_mesh.edges.shape[1] == 2
    assert mesh.areas.shape == (4,)
    assert np.all(mesh.areas >= 0.0)
    assert len(mesh.boundary_indices) == 4


def test_finite_volume_utilities_are_consistent():
    elements = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.int64)
    edges, is_boundary = get_edges(elements)
    assert edges.shape[1] == 2
    assert np.count_nonzero(is_boundary) == 4
    sites = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
    areas = triangle_areas(sites, elements)
    assert np.allclose(np.abs(areas), 0.5)
    closed = close_curve(sites)
    assert np.allclose(closed[0], closed[-1])
