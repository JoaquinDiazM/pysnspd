from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib.tri as mtri

from pysnspd.analysis.photon_snapshots import (
    compute_photon_snapshot_plot_diagnostics,
)
from pysnspd.plotting.photon_figures import phase_gradient_q_abs_m_inv


def test_photon_snapshot_diagnostics_use_all_snapshots_for_limits(
    tmp_path: Path,
    small_strip_mesh_bundle,
) -> None:
    mesh, _, ops = small_strip_mesh_bundle
    n_nodes = mesh.n_nodes
    x = np.asarray(mesh.nodes, dtype=float)[:, 0]
    x_unit = (x - np.min(x)) / max(float(np.ptp(x)), 1.0e-300)
    phase_scale = np.array([0.1, 0.4, 0.2])
    voltage_scale = np.array([1.0, 4.0, 2.0]) * 1.0e-4
    phase = phase_scale[:, None] * x_unit[None, :]
    amplitude_J = np.full((3, n_nodes), 1.0e-22)
    Te = np.vstack(
        [
            np.full(n_nodes, 0.9),
            0.9 + 0.5 * x_unit,
            0.9 + 0.2 * x_unit,
        ]
    )
    Tph = np.vstack(
        [
            np.full(n_nodes, 0.9),
            0.9 + 0.2 * x_unit,
            0.9 + 0.1 * x_unit,
        ]
    )
    snapshots = {
        "snapshot_t_ps": np.array([0.0, 1.0, 2.0]),
        "psi_real_snapshot_J": amplitude_J * np.cos(phase),
        "psi_imag_snapshot_J": amplitude_J * np.sin(phase),
        "delta_snapshot_meV": amplitude_J / 1.602176634e-22,
        "phi_snapshot_V": voltage_scale[:, None] * x_unit[None, :],
        "Te_snapshot_K": Te,
        "Tph_snapshot_K": Tph,
    }

    Te_axis = np.array([0.9, 4.0])
    Tph_axis = np.array([0.9, 4.0])
    delta_axis = np.array([0.0, 2.0e-22])
    q_axis = np.array([0.0, 1.0e9])
    shape4 = (2, 2, 2, 2)
    table = np.ones(shape4)
    catalog = tmp_path / "power_table_catalog.npz"
    np.savez_compressed(
        catalog,
        Te_values_K=Te_axis,
        Tph_values_K=Tph_axis,
        delta_values_J=delta_axis,
        q_values_m_inv=q_axis,
        P_S_W_m3=table,
        P_R_W_m3=2.0 * table,
        P_total_W_m3=3.0 * table,
        u_e_J_m3=np.ones((2, 2, 2)),
        C_e_J_m3_K=2.0 * np.ones((2, 2, 2)),
        kappa_s_W_m_K=np.ones((2, 2)),
        u_ph_J_m3=np.array([0.1, 1.0]),
        C_ph_J_m3_K=np.array([0.2, 2.0]),
        P_esc_W_m3=np.array([0.0, 3.0]),
    )

    result = compute_photon_snapshot_plot_diagnostics(
        snapshots=snapshots,
        selected_indices=[0, 2],
        power_table_npz=str(catalog),
        ops=ops,
        sigma_n_S_m=4.2e5,
        thermal_active_mask=np.ones(n_nodes, dtype=bool),
        thermal_bath_K=0.9,
        chunk_size=2,
    )

    assert result["joule_snapshot_W_m3"].shape == (2, n_nodes)
    assert result["P_diff_snapshot_W_m3"].shape == (2, n_nodes)
    assert result["u_e_snapshot_J_m3"].shape == (2, n_nodes)
    limits = result["snapshot_global_limits"]
    assert limits["joule_snapshot_W_m3"][1] > np.max(
        result["joule_snapshot_W_m3"]
    )
    assert limits["q_abs_snapshot_m_inv"][1] > np.max(
        result["q_abs_snapshot_m_inv"]
    )
    triangulation = mtri.Triangulation(
        1.0e9 * mesh.nodes[:, 0],
        1.0e9 * mesh.nodes[:, 1],
        mesh.triangles,
    )
    q_reference = phase_gradient_q_abs_m_inv(
        triangulation,
        snapshots["psi_real_snapshot_J"][2]
        + 1j * snapshots["psi_imag_snapshot_J"][2],
        x_nm=triangulation.x,
        y_nm=triangulation.y,
    )
    np.testing.assert_allclose(result["q_abs_snapshot_m_inv"][1], q_reference)
    weighted_diffusion = np.sum(
        result["P_diff_snapshot_W_m3"] * ops.node_area_m2[None, :], axis=1
    )
    np.testing.assert_allclose(weighted_diffusion, 0.0, atol=1.0e-15)

    reused = compute_photon_snapshot_plot_diagnostics(
        snapshots=snapshots,
        selected_indices=[0, 2],
        power_table_npz=str(catalog),
        ops=ops,
        sigma_n_S_m=4.2e5,
        thermal_active_mask=np.ones(n_nodes, dtype=bool),
        thermal_bath_K=0.9,
        global_limits_override=limits,
        chunk_size=2,
    )
    assert reused["global_limits_reused"] is True
    np.testing.assert_allclose(
        reused["joule_snapshot_W_m3"], result["joule_snapshot_W_m3"]
    )
    for key in limits:
        np.testing.assert_allclose(reused["snapshot_global_limits"][key], limits[key])
