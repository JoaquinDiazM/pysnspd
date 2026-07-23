"""NPZ persistence tests for gTDGL states and histories."""
from __future__ import annotations

import numpy as np

from pysnspd.gtdgl.state import CurrentFields, GTDGLStationaryState
from pysnspd.solver.state_io import save_relaxation_history_npz, save_stationary_state_npz


def _zero_currents(n_nodes: int, n_edges: int) -> CurrentFields:
    edge = np.zeros(n_edges)
    node = np.zeros(n_nodes)
    return CurrentFields(
        edge_Q_m_inv=edge.copy(),
        edge_js_us_A_m2=edge.copy(),
        edge_js_gl_A_m2=edge.copy(),
        edge_jn_A_m2=edge.copy(),
        edge_jtot_A_m2=edge.copy(),
        node_div_js_us_A_m3=node.copy(),
        node_div_js_gl_A_m3=node.copy(),
        node_div_jtot_A_m3=node.copy(),
        node_js_us_x_A_m2=node.copy(),
        node_js_us_y_A_m2=node.copy(),
        node_jn_x_A_m2=node.copy(),
        node_jn_y_A_m2=node.copy(),
        node_jtot_x_A_m2=node.copy(),
        node_jtot_y_A_m2=node.copy(),
        edge_pairbreaking_ratio=edge.copy(),
        node_pairbreaking_ratio=node.copy(),
    )


def test_save_stationary_state_npz(tmp_path, small_strip_mesh_bundle, gtdgl_material):
    mesh, _, ops = small_strip_mesh_bundle
    state = GTDGLStationaryState(
        psi_J=np.full(mesh.n_nodes, gtdgl_material.delta0_J, dtype=np.complex128),
        phi_V=np.zeros(mesh.n_nodes),
        Te_K=np.full(mesh.n_nodes, 0.9),
        Tph_K=np.full(mesh.n_nodes, 0.9),
        currents=_zero_currents(mesh.n_nodes, ops.n_edges),
        metadata={"backend": "flat-gtdgl-test"},
    )
    out = save_stationary_state_npz(state, tmp_path / "state.npz")
    assert out.exists()
    with np.load(out) as data:
        assert "psi_real_J" in data.files
        assert "metadata_json" in data.files


def test_save_relaxation_history_npz(tmp_path):
    out = save_relaxation_history_npz({"t_s": np.array([0.0, 1.0])}, tmp_path / "history.npz")
    assert out.exists()
    with np.load(out) as data:
        assert np.array_equal(data["t_s"], np.array([0.0, 1.0]))
