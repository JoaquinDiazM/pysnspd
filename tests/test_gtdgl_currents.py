"""Current adapter tests for flat gTDGL."""
from __future__ import annotations

import numpy as np

from pysnspd.gtdgl.currents import edge_to_node_weighted_average, pairbreaking_ratio_edges


def test_pairbreaking_ratio_edges_is_finite(small_strip_mesh_bundle, gtdgl_material):
    _, _, ops = small_strip_mesh_bundle
    Te = np.full(ops.n_nodes, 0.9)
    ratio = pairbreaking_ratio_edges(
        edge_current_A_m2=np.ones(ops.n_edges),
        material=gtdgl_material,
        Te_K=Te,
        ops=ops,
    )
    assert ratio.shape == (ops.n_edges,)
    assert np.all(np.isfinite(ratio))
    assert np.all(ratio >= 0.0)


def test_edge_to_node_weighted_average_shape(small_strip_mesh_bundle):
    _, _, ops = small_strip_mesh_bundle
    node_values = edge_to_node_weighted_average(np.arange(ops.n_edges, dtype=float), ops)
    assert node_values.shape == (ops.n_nodes,)
    assert np.all(np.isfinite(node_values))
