from __future__ import annotations

import numpy as np

from pysnspd.mesh.delaunay import generate_rectangular_delaunay_mesh, triangle_areas
from pysnspd.mesh.edges import build_edge_data, assert_edge_data_consistent


def _cfg():
    return {
        "project": {"default_run_name": "mesh_test"},
        "paths": {"big_data_root": "/tmp"},
        "material": {
            "name": "NbN",
            "Tc_K": 8.65,
            "D_m2_s": 1.58e-4,
            "sigma_n_S_m": 4.2e5,
            "lambda_L_m": 5.4e-7,
            "thickness_m": 7e-9,
            "width_m": 1.2e-7,
        },
        "bias": {"T_bias_K": 0.9, "I_bias_A": 35e-6},
        "mesh": {"target_spacing_m": 1.2e-8, "seed": 123, "length_m": 2.4e-7},
    }


def test_pytdgl_like_mesh_is_unstructured_and_boundary_tagged():
    mesh = generate_rectangular_delaunay_mesh(
        _cfg(),
        jitter_fraction=0.20,
        boundary_guard_layers=2,
    )

    assert mesh.triangulation_method == "pytdgl_generate_mesh_meshpy_triangle_v1"
    assert mesh.n_nodes > 20
    assert mesh.n_triangles > 20
    assert np.all(triangle_areas(mesh.nodes, mesh.triangles) > 0.0)

    edge_data = build_edge_data(
        mesh.nodes,
        mesh.triangles,
        length_m=mesh.length_m,
        width_m=mesh.width_m,
    )
    assert_edge_data_consistent(edge_data)

    tags = set(np.asarray(edge_data.tags).astype(str))
    assert {"left", "right", "top", "bottom"}.issubset(tags)


def test_pytdgl_like_mesh_is_reproducible():
    a = generate_rectangular_delaunay_mesh(_cfg(), jitter_fraction=0.17, boundary_guard_layers=1)
    b = generate_rectangular_delaunay_mesh(_cfg(), jitter_fraction=0.17, boundary_guard_layers=1)
    assert np.allclose(a.nodes, b.nodes)
    assert np.array_equal(a.triangles, b.triangles)
