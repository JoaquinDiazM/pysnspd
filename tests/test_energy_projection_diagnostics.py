from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np

from pysnspd.analysis.energy_projection import (
    CompressedNpzRowStream,
    extract_energy_projection_diagnostics,
)
from pysnspd.plotting.energy_projection import write_energy_projection_figures


def test_compressed_npz_row_stream_reconstructs_members(tmp_path: Path) -> None:
    path = tmp_path / "rows.npz"
    first = np.arange(30, dtype=np.float64).reshape(6, 5)
    second = first + 100.0
    np.savez_compressed(path, first=first, second=second)

    chunks = list(CompressedNpzRowStream(path, ("first", "second")).iter_chunks(chunk_rows=4))

    assert [offset for offset, _ in chunks] == [0, 4]
    np.testing.assert_array_equal(np.concatenate([chunk["first"] for _, chunk in chunks]), first)
    np.testing.assert_array_equal(np.concatenate([chunk["second"] for _, chunk in chunks]), second)


def test_energy_projection_closes_storage_decomposition_and_writes_three_pdfs(
    tmp_path: Path,
) -> None:
    snapshots_path = tmp_path / "snapshots.npz"
    table_path = tmp_path / "power_table.npz"
    current_table_path = tmp_path / "usadel_current.npz"
    nodes, triangles, ops = _small_mesh()
    n_snapshots = 6
    time_ps = np.asarray([0.0, 1.0, 2.0, 3.0, 4.0, 4.0])
    x = nodes[:, 0]
    amplitudes = 1.0 - 0.04 * np.arange(n_snapshots, dtype=float)[:, None]
    phases = 0.10 * np.arange(n_snapshots, dtype=float)[:, None] * x[None, :]
    psi = amplitudes * np.exp(1j * phases)
    Te = 1.2 + 0.05 * np.arange(n_snapshots, dtype=float)[:, None] + np.zeros((n_snapshots, nodes.shape[0]))
    Tph = 1.0 + 0.01 * np.arange(n_snapshots, dtype=float)[:, None] + np.zeros_like(Te)
    phi = 0.02 * np.arange(n_snapshots, dtype=float)[:, None] * x[None, :]
    np.savez_compressed(
        snapshots_path,
        snapshot_t_ps=time_ps,
        psi_real_snapshot_J=np.real(psi),
        psi_imag_snapshot_J=np.imag(psi),
        phi_snapshot_V=phi,
        Te_snapshot_K=Te,
        Tph_snapshot_K=Tph,
    )
    _write_small_power_table(table_path)
    _write_small_current_table(current_table_path)
    history = {
        "t_ps": time_ps,
        "V_tdgl_center_V": np.asarray([0.0, 0.2, 1.0, 0.4, 0.1, 0.1]),
        "photon_applied": np.asarray([False, True, False, False, False, False]),
    }

    result = extract_energy_projection_diagnostics(
        snapshots_npz=snapshots_path,
        power_table_npz=table_path,
        usadel_current_npz=current_table_path,
        history=history,
        nodes_m=nodes,
        triangles=triangles,
        ops=ops,
        sigma_n_S_m=2.0,
        thickness_m=0.5,
        Tc_K=3.0,
        xi_m=0.2,
        window_m=2.0,
        chunk_rows=2,
    )

    np.testing.assert_allclose(
        result["integrated_P_spec_W"],
        result["integrated_P_delta_W"] + result["integrated_P_q_W"],
        rtol=1.0e-11,
        atol=1.0e-11,
    )
    assert result["time_ps"].shape == (n_snapshots - 2,)
    assert result["selected_times_ps"].shape == (3,)
    assert result["selected_P_spec_W_m3"].shape == (3, nodes.shape[0])
    assert int(result["dropped_duplicate_time_count"][0]) == 1
    np.testing.assert_allclose(
        result["selected_delta_over_delta0"][-1],
        np.abs(psi[-1]) / 2.0,
    )
    assert not bool(result["truncated"][0])
    assert float(result["strict_q_max_m_inv"][0]) == 0.25
    assert np.nanmax(result["strict_q_clipped_fraction"]) > 0.0
    finite_ratios = result["strict_q_clipped_js_p95_over_catalog_max"]
    assert np.all((finite_ratios[np.isfinite(finite_ratios)] >= 0.0))

    saved = write_energy_projection_figures(result, tmp_path / "figures", dpi=72)
    assert set(saved) == {"colormaps", "temporal", "profiles"}
    assert all(path.exists() and path.stat().st_size > 1000 for path in saved.values())


def _small_mesh() -> tuple[np.ndarray, np.ndarray, SimpleNamespace]:
    nodes = np.asarray(
        [[-0.5, -0.5], [0.5, -0.5], [0.5, 0.5], [-0.5, 0.5]],
        dtype=float,
    )
    triangles = np.asarray([[0, 1, 2], [0, 2, 3]], dtype=np.int64)
    edge_i = np.asarray([0, 1, 2, 3, 0], dtype=np.int64)
    edge_j = np.asarray([1, 2, 3, 0, 2], dtype=np.int64)
    vectors = nodes[edge_j] - nodes[edge_i]
    lengths = np.linalg.norm(vectors, axis=1)
    return nodes, triangles, SimpleNamespace(
        edge_i=edge_i,
        edge_j=edge_j,
        edge_length_m=lengths,
        dual_face_length_m=np.ones_like(lengths),
        node_area_m2=np.full(nodes.shape[0], 0.25),
        edge_unit=vectors / lengths[:, None],
        n_nodes=nodes.shape[0],
    )


def _write_small_power_table(path: Path) -> None:
    Te = np.asarray([1.0, 2.0, 3.0])
    Tph = np.asarray([1.0, 2.0, 3.0])
    delta = np.asarray([0.0, 1.0, 2.0])
    q = np.asarray([0.0, 1.0, 2.0])
    ue = (
        Te[:, None, None]
        + 2.0 * delta[None, :, None]
        + 3.0 * q[None, None, :]
    )
    P_total = (
        Te[:, None, None, None]
        - Tph[None, :, None, None]
        + np.zeros((Te.size, Tph.size, delta.size, q.size))
    )
    np.savez_compressed(
        path,
        Te_values_K=Te,
        Tph_values_K=Tph,
        delta_values_J=delta,
        q_values_m_inv=q,
        u_e_J_m3=ue,
        C_e_J_m3_K=np.ones_like(ue),
        P_total_W_m3=P_total,
        kappa_s_W_m_K=np.ones((Te.size, delta.size)),
        u_ph_J_m3=Tph,
        P_esc_W_m3=Tph - 1.0,
        metadata=np.asarray(
            {
                "N0_J_m3": 1.0,
                "delta0_J": 2.0,
                "T_bath_K": 1.0,
            },
            dtype=object,
        ),
    )


def _write_small_current_table(path: Path) -> None:
    Te = np.asarray([1.0, 2.0, 3.0])
    delta = np.asarray([0.0, 1.0, 2.0])
    q = np.asarray([0.0, 0.125, 0.25])
    js = (
        (1.0 + Te[:, None, None])
        * (1.0 + delta[None, :, None])
        * q[None, None, :]
    )
    np.savez_compressed(
        path,
        Te_axis_K=Te,
        delta_axis_J=delta,
        q_axis_m_inv=q,
        js_A_m2=js,
    )
