"""Functional smoke for photon transient closure and output serialization."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import yaml

from pysnspd.circuit.readout import CircuitParams
from pysnspd.excitation.photon import PhotonBubbleParams
from pysnspd.solver.transient import CoupledTransientConfig, run_coupled_transient


def test_coupled_transient_reaches_final_serialization(
    tmp_path,
    small_strip_mesh_bundle,
    gtdgl_material,
    stationary_seed_factory,
):
    mesh, edge_data, ops = small_strip_mesh_bundle
    seed = stationary_seed_factory(mesh, gtdgl_material)
    initial_state = tmp_path / "stationary_state.npz"
    np.savez_compressed(
        initial_state,
        psi_real_J=seed.node_psi_real_J,
        psi_imag_J=seed.node_psi_imag_J,
        phi_V=seed.node_phi_electric_V,
        Te_K=seed.node_Te_K,
        Tph_K=seed.node_Tph_K,
    )

    usadel_catalog = SimpleNamespace(
        js_A_m2=np.zeros((2, 2, 2), dtype=float),
        js_stiffness_A_per_m_J2=np.ones((2, 2, 2), dtype=float) * 1.0e48,
        Te_axis_K=np.array([0.9, 1.0]),
        delta_axis_J=np.array([0.0, gtdgl_material.delta0_J]),
        delta2_axis_J2=np.array([0.0, gtdgl_material.delta0_J**2]),
        q_axis_m_inv=np.array([0.0, 1.0e8]),
    )
    output_dir = tmp_path / "photon"
    summary = run_coupled_transient(
        mesh=mesh,
        edge_data=edge_data,
        ops=ops,
        material=gtdgl_material,
        initial_state_npz=initial_state,
        initial_current_A=0.0,
        usadel_catalog=usadel_catalog,
        power_table_npz=None,
        output_dir=output_dir,
        config=CoupledTransientConfig(
            total_time_s=2.0e-18,
            mesoscopic_dt_s=1.0e-18,
            chunk_time_s=1.0e-18,
            n_snapshots=2,
            thermal_enabled=False,
            terminal_psi=0.0,
            early_stop_mode="none",
            progress=False,
        ),
        circuit_params=CircuitParams(),
        photon_params=PhotonBubbleParams(enabled=False),
    )

    assert summary["stop_reason"] == "requested_time_reached"
    assert summary["circuit"]["params"]["R_load_ohm"] == 50.0
    assert summary["circuit"]["params"]["V_bias_V"] == 0.0
    assert summary["circuit"]["final_state"]["I_s_A"] == 0.0
    assert (output_dir / "final_state.npz").exists()
    assert (output_dir / "transient_history.npz").exists()
    assert (output_dir / "transient_snapshots.npz").exists()
    assert (output_dir / "timing_summary.yaml").exists()
    summary_path = output_dir / "photon_summary.yaml"
    assert summary_path.exists()
    with summary_path.open("r", encoding="utf-8") as stream:
        persisted = yaml.safe_load(stream)
    assert persisted["circuit"]["params"] == summary["circuit"]["params"]
