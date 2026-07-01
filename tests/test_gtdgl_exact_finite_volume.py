"""Tests for the pyTDGL-exact finite-volume mesh data model."""

from __future__ import annotations

import numpy as np

from pysnspd.gtdgl.finite_volume import Mesh
from pysnspd.gtdgl.finite_volume.edge_mesh import EdgeMesh
from pysnspd.gtdgl.finite_volume.util import (
    generate_voronoi_vertices,
    get_edges,
    get_voronoi_polygon_indices,
)


def test_exact_mesh_from_single_square_triangulation_has_edge_mesh():
    sites = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
        ]
    )
    elements = np.array(
        [
            [0, 1, 2],
            [0, 2, 3],
        ],
        dtype=np.int64,
    )

    mesh = Mesh.from_triangulation(sites, elements, create_submesh=True)

    assert mesh.sites.shape == (4, 2)
    assert mesh.elements.shape == (2, 3)
    assert isinstance(mesh.edge_mesh, EdgeMesh)
    assert mesh.areas.shape == (4,)
    assert len(mesh.voronoi_polygons) == 4
    assert np.all(mesh.areas >= 0.0)


def test_get_edges_boundary_flags_match_pytdgl_convention():
    elements = np.array(
        [
            [0, 1, 2],
            [0, 2, 3],
        ],
        dtype=np.int64,
    )

    edges, is_boundary = get_edges(elements)

    assert edges.shape == (5, 2)
    assert is_boundary.sum() == 4
    assert (~is_boundary).sum() == 1


def test_voronoi_polygon_indices_use_directed_triangle_indices():
    elements = np.array(
        [
            [0, 1, 2],
            [0, 2, 3],
        ],
        dtype=np.int64,
    )

    dual = generate_voronoi_vertices(
        np.array(
            [
                [0.0, 0.0],
                [1.0, 0.0],
                [1.0, 1.0],
                [0.0, 1.0],
            ]
        ),
        elements,
    )
    polygons = get_voronoi_polygon_indices(elements, 4)

    assert dual.shape == (2, 2)
    assert len(polygons) == 4
    assert all(np.all(poly >= 0) for poly in polygons if len(poly))
