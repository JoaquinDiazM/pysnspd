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

from pathlib import Path

from pysnspd.usadel.supercurrent_table import append_usadel_supercurrent_table_to_npz


def test_dirty_usadel_supercurrent_table_parallel_matches_serial():
    delta0 = 1.3148844131376842 * 1.602176634e-22
    q = np.array([-5.0e7, 0.0, 5.0e7])
    delta = np.array([0.0, 0.5 * delta0, delta0])
    serial = compute_dirty_usadel_supercurrent_table(
        delta_axis_J=delta,
        q_axis_m_inv=q,
        T_K=0.9,
        D_m2_s=1.58e-4,
        sigma_n_S_m=4.2e5,
        n_matsubara=80,
        workers=1,
    )
    parallel = compute_dirty_usadel_supercurrent_table(
        delta_axis_J=delta,
        q_axis_m_inv=q,
        T_K=0.9,
        D_m2_s=1.58e-4,
        sigma_n_S_m=4.2e5,
        n_matsubara=80,
        workers=2,
    )
    assert np.allclose(serial, parallel)


def test_append_usadel_supercurrent_table_preserves_object_arrays(tmp_path: Path):
    npz = tmp_path / "usadel_dos_catalog.npz"
    summary = tmp_path / "usadel_dos_summary.yaml"
    out = tmp_path / "usadel_supercurrent_table_summary.yaml"
    np.savez_compressed(
        npz,
        rho=np.ones((2, 3, 4)),
        metadata=np.array([{"source": "test"}], dtype=object),
    )
    summary.write_text(
        "\n".join(
            [
                "n_delta: 2",
                "n_q: 3",
                "delta_min_J: 0.0",
                "delta_max_J: 2.10667708314e-22",
                "q_min_m_inv: 0.0",
                "q_max_m_inv: 5.0e7",
                "D_m2_s: 1.58e-4",
            ]
        ),
        encoding="utf-8",
    )
    cfg = {
        "bias": {"T_bias_K": 0.9},
        "material": {"D_m2_s": 1.58e-4, "sigma_n_S_m": 4.2e5},
    }
    result = append_usadel_supercurrent_table_to_npz(
        npz,
        config=cfg,
        summary_yaml=summary,
        output_summary_yaml=out,
        n_matsubara=40,
        workers=2,
        progress=False,
    )
    assert result.table_shape == (2, 3)
    with np.load(npz, allow_pickle=True) as data:
        assert "metadata" in data.files
        assert "js_A_m2" in data.files
        assert data["js_A_m2"].shape == (2, 3)
