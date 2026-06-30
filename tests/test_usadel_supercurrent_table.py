import numpy as np

from pysnspd.usadel.supercurrent_table import compute_dirty_usadel_supercurrent_table


def test_dirty_usadel_supercurrent_table_has_expected_sign_and_scale():
    delta0 = 1.3148844131376842 * 1.602176634e-22
    q = np.array([0.0, 5.772939103e7])
    delta = np.array([0.0, 0.891 * delta0])
    table = compute_dirty_usadel_supercurrent_table(
        delta_axis_J=delta,
        q_axis_m_inv=q,
        T_K=0.9,
        D_m2_s=1.5813313145524903e-4,
        sigma_n_S_m=4.2e5,
        n_matsubara=500,
    )
    assert table.shape == (2, 2)
    assert np.allclose(table[:, 0], 0.0)
    assert np.allclose(table[0, :], 0.0)
    assert 3.0e10 < table[1, 1] < 5.5e10


def test_dirty_usadel_supercurrent_table_is_odd_in_q():
    delta0 = 1.3148844131376842 * 1.602176634e-22
    q = np.array([-5.0e7, 0.0, 5.0e7])
    delta = np.array([delta0])
    table = compute_dirty_usadel_supercurrent_table(
        delta_axis_J=delta,
        q_axis_m_inv=q,
        T_K=0.9,
        D_m2_s=1.58e-4,
        sigma_n_S_m=4.2e5,
        n_matsubara=200,
    )
    assert table[0, 0] < 0.0
    assert table[0, 2] > 0.0
    assert np.isclose(table[0, 0], -table[0, 2], rtol=1e-10, atol=0.0)
