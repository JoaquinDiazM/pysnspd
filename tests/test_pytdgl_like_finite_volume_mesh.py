from __future__ import annotations

import importlib.util

import numpy as np
import pytest

from pysnspd.gtdgl.pytdgl_like.finite_volume.mesh import Mesh
from pysnspd.mesh.pytdgl_like import rectangular_boundary_points


def test_mesh_from_triangulation_builds_edge_and_voronoi_data():
    sites = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
            [0.5, 0.5],
        ]
    )
    elements = np.array([[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]])
    mesh = Mesh.from_triangulation(sites, elements)
    assert mesh.edge_mesh is not None
    assert mesh.dual_sites.shape == (4, 2)
    assert mesh.areas.shape == (5,)
    assert len(mesh.voronoi_polygons) == 5
    assert np.all(mesh.edge_mesh.edge_lengths > 0)
    assert np.all(mesh.edge_mesh.dual_edge_lengths >= 0)


@pytest.mark.skipif(importlib.util.find_spec("meshpy") is None, reason="meshpy is required for pyTDGL-faithful mesh generation")
def test_generate_meshpy_rectangular_mesh_has_voronoi_submesh():
    from pysnspd.mesh.pytdgl_like import (
        PyTDGLLikeMeshParameters,
        generate_rectangular_pytdgl_fvm_mesh_from_parameters,
    )

    params = PyTDGLLikeMeshParameters(
        length_m=2.4e-7,
        width_m=1.2e-7,
        target_spacing_m=1.2e-8,
        max_edge_length_m=1.2e-8,
        seed=123,
        smooth=2,
        min_angle_deg=28.0,
    )
    mesh = generate_rectangular_pytdgl_fvm_mesh_from_parameters(params)
    assert mesh.sites.shape[1] == 2
    assert mesh.elements.shape[1] == 3
    assert mesh.edge_mesh is not None
    assert np.all(mesh.areas >= 0)
    assert np.max(mesh.edge_mesh.edge_lengths) <= 1.7 * params.max_edge_length_m


def test_rectangular_boundary_points_are_open_curve():
    pts = rectangular_boundary_points(2.4e-7, 1.2e-7, 1.2e-8)
    assert pts.shape[1] == 2
    assert not np.allclose(pts[0], pts[-1])
    assert np.isclose(pts[:, 0].min(), 0.0)
    assert np.isclose(pts[:, 0].max(), 2.4e-7)
