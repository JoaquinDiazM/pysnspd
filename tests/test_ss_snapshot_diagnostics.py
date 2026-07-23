from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np

from pysnspd.analysis.snapshots import (
    compute_ss_snapshot_power_diagnostics,
    save_ss_snapshot_bundle_npz,
)


def test_ss_snapshot_power_diagnostics_shapes(tmp_path: Path) -> None:
    Te = np.array([0.9, 4.0])
    Tph = np.array([0.9, 4.0])
    delta = np.array([0.0, 2.0e-22])
    q = np.array([0.0, 5.0e7])
    shape = (Te.size, Tph.size, delta.size, q.size)
    P = np.ones(shape) * 3.0
    cat = tmp_path / "power_table_catalog.npz"
    np.savez_compressed(
        cat,
        Te_values_K=Te,
        Tph_values_K=Tph,
        delta_values_J=delta,
        q_values_m_inv=q,
        P_S_W_m3=P,
        P_R_W_m3=2.0 * P,
        P_total_W_m3=3.0 * P,
        u_e_J_m3=np.ones((Te.size, delta.size, q.size)),
        C_e_J_m3_K=2.0 * np.ones((Te.size, delta.size, q.size)),
        kappa_s_W_m_K=np.ones((Te.size, delta.size)),
        u_ph_J_m3=np.array([0.1, 10.0]),
        C_ph_J_m3_K=np.array([0.2, 20.0]),
        P_esc_W_m3=np.array([0.0, 30.0]),
        metadata=np.array({"backend": "test"}, dtype=object),
    )

    history = {
        "snapshot_t_s": np.array([0.0, 1.0e-12, 2.0e-12]),
        "delta_snapshot_meV": np.array([[0.0, 1.0], [1.0, 1.0], [0.0, 0.0]]),
        "edge_Q_snapshot_m_inv": np.array([[0.0], [5.0e7], [1.0e7]]),
        "edge_jn_snapshot_A_m2": np.array([[0.0], [1.0], [2.0]]),
        "edge_jtot_snapshot_A_m2": np.array([[0.0], [3.0], [4.0]]),
        "edge_i": np.array([0]),
        "edge_j": np.array([1]),
        "edge_length_m": np.array([1.0]),
        "dual_face_length_m": np.array([1.0]),
    }
    state = SimpleNamespace(Te_K=np.array([0.9, 4.0]), Tph_K=np.array([0.9, 4.0]))

    out = compute_ss_snapshot_power_diagnostics(
        history=history,
        state=state,
        power_table_npz=cat,
        sigma_n_S_m=2.0,
    )

    assert out["P_total_snapshot_W_m3"].shape == (3, 2)
    assert out["u_e_snapshot_J_m3"].shape == (3, 2)
    assert out["kappa_s_snapshot_W_m_K"].shape == (3, 2)
    assert out["P_esc_snapshot_W_m3"].shape == (3, 2)
    assert out["joule_snapshot_W_m3"].shape == (3, 2)
    assert np.all(np.isfinite(out["P_total_snapshot_W_m3"]))

    snapshots = save_ss_snapshot_bundle_npz(history, tmp_path / "stationary_snapshots.npz")
    assert snapshots.exists()
    with np.load(snapshots) as data:
        assert "delta_snapshot_meV" in data.files
        assert "edge_i" in data.files
