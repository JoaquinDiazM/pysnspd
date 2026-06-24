from __future__ import annotations

import numpy as np

from pysnspd.gtdgl.operators import (
    FVOperators,
    edge_scalar_to_node_vector_least_squares,
)


def test_least_squares_reconstructs_uniform_vector_from_edge_projections():
    nodes = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
        ],
        dtype=float,
    )

    edges = np.array(
        [
            [0, 1],
            [0, 2],
            [1, 3],
            [2, 3],
            [0, 3],
            [1, 2],
        ],
        dtype=np.int64,
    )

    edge_i = edges[:, 0]
    edge_j = edges[:, 1]
    edge_vec = nodes[edge_j] - nodes[edge_i]
    edge_len = np.linalg.norm(edge_vec, axis=1)
    edge_unit = edge_vec / edge_len[:, None]

    ops = FVOperators(
        edges=edges,
        edge_i=edge_i,
        edge_j=edge_j,
        edge_vec_m=edge_vec,
        edge_length_m=edge_len,
        edge_unit=edge_unit,
        dual_face_length_m=np.ones(edges.shape[0]),
        node_area_m2=np.ones(nodes.shape[0]),
    )

    j_true = np.array([3.2, -1.7])
    edge_projection = edge_unit @ j_true

    jx, jy = edge_scalar_to_node_vector_least_squares(edge_projection, ops)

    assert np.allclose(jx, j_true[0], rtol=1.0e-12, atol=1.0e-12)
    assert np.allclose(jy, j_true[1], rtol=1.0e-12, atol=1.0e-12)