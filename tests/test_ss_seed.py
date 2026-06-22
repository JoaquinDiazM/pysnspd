from types import SimpleNamespace

import numpy as np

from pysnspd.gtdgl.seed import (
    build_stationary_seed,
    compute_boundary_currents,
    select_bias_state_from_usadel,
)


def _toy_mesh():
    length = 2.0
    width = 1.0
    nodes = np.array(
        [
            [0.0, -0.5],
            [1.0, -0.5],
            [2.0, -0.5],
            [0.0, 0.0],
            [1.0, 0.0],
            [2.0, 0.0],
            [0.0, 0.5],
            [1.0, 0.5],
            [2.0, 0.5],
        ],
        dtype=float,
    )
    triangles = np.array(
        [
            [0, 1, 4],
            [0, 4, 3],
            [1, 2, 5],
            [1, 5, 4],
            [3, 4, 7],
            [3, 7, 6],
            [4, 5, 8],
            [4, 8, 7],
        ],
        dtype=np.int64,
    )
    return SimpleNamespace(
        nodes=nodes,
        triangles=triangles,
        length_m=length,
        width_m=width,
    )


def _toy_edges():
    # Boundary edges only are enough for this OE6 seed test.
    edges = np.array(
        [
            [0, 3],
            [3, 6],
            [2, 5],
            [5, 8],
            [0, 1],
            [1, 2],
            [6, 7],
            [7, 8],
        ],
        dtype=np.int64,
    )
    tags = np.array(
        ["left", "left", "right", "right", "bottom", "bottom", "top", "top"]
    )
    lengths = np.ones(edges.shape[0], dtype=float) * 0.5
    lengths[4:] = 1.0

    return SimpleNamespace(
        edges=edges,
        tags=tags,
        lengths=lengths,
    )


def _toy_usadel_catalog():
    q = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    current = np.array([0.0, 1.0, 2.0, 2.5, 2.2]) * 1.0e-6
    width = 1.0
    thickness = 0.1
    current_density = current / (width * thickness)
    delta = np.array([1.3, 1.25, 1.1, 0.9, 0.5]) * 1.0e-22
    gamma = np.array([0.0, 0.1, 0.4, 0.9, 1.6]) * 1.0e-22

    return SimpleNamespace(
        calibration_q_values_m_inv=q,
        calibration_current_values_A=current,
        calibration_current_density_values_A_m2=current_density,
        calibration_delta_eq_values_J=delta,
        calibration_gamma_values_J=gamma,
        delta_values_J=np.linspace(0.0, 1.3e-22, 4),
        metadata={
            "I_bias_A": 1.5e-6,
            "T_bias_K": 0.9,
            "thickness_m": thickness,
            "width_m": width,
            "delta0_J": 1.3e-22,
        },
    )


def test_select_bias_state_uses_stable_branch():
    catalog = _toy_usadel_catalog()
    state = select_bias_state_from_usadel(catalog)

    assert state.I_bias_A == 1.5e-6
    assert np.isclose(state.Ic_A, 2.5e-6)
    assert state.q_bias_m_inv < state.q_critical_m_inv
    assert state.branch_policy == "stable_usadel_branch_q_le_qc"
    assert state.delta_bias_J > 0.0


def test_boundary_currents_for_uniform_longitudinal_current():
    edges = _toy_edges()
    out = compute_boundary_currents(
        edge_data=edges,
        jx_A_m2=10.0,
        jy_A_m2=0.0,
        thickness_m=0.1,
    )

    assert np.isclose(out["left_A"], -1.0)
    assert np.isclose(out["right_A"], 1.0)
    assert np.isclose(out["top_A"], 0.0)
    assert np.isclose(out["bottom_A"], 0.0)
    assert np.isclose(out["net_boundary_current_A"], 0.0)


def test_stationary_seed_is_uniform_and_zero_voltage():
    mesh = _toy_mesh()
    edges = _toy_edges()
    catalog = _toy_usadel_catalog()

    seed = build_stationary_seed(
        mesh=mesh,
        edge_data=edges,
        usadel_catalog=catalog,
    )

    assert np.allclose(seed.node_Te_K, 0.9)
    assert np.allclose(seed.node_Tph_K, 0.9)
    assert np.all(seed.node_delta_J > 0.0)
    assert np.allclose(seed.node_phi_electric_V, 0.0)
    assert np.allclose(seed.node_jn_x_A_m2, 0.0)
    assert np.allclose(seed.node_jn_y_A_m2, 0.0)
    assert np.allclose(seed.node_div_j_A_m3, 0.0)

    meta = seed.metadata
    assert np.isclose(meta["terminal_voltage_V"], 0.0)
    assert meta["right_current_error_rel"] < 1.0e-12
    assert meta["left_current_error_rel"] < 1.0e-12