"""Stationary adapter tests for the strict 3D Usadel-Poisson supercurrent law."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from pysnspd.gtdgl.adapter import solve_stationary_pytdgl_like


def _strict_3d_catalog(material) -> SimpleNamespace:
    q_axis = np.array([0.0, 1.0e7, 2.0e7, 4.0e7], dtype=float)
    delta_axis = np.array([0.0, 0.5, 0.9, 1.0], dtype=float) * material.delta0_J
    Te_axis = np.array([0.5, 0.9, 1.2], dtype=float)
    js = np.empty((Te_axis.size, delta_axis.size, q_axis.size), dtype=float)
    for iT, T in enumerate(Te_axis):
        temp_factor = 1.0 - 0.05 * (T - 0.9)
        for iD, delta in enumerate(delta_axis):
            js[iT, iD, :] = temp_factor * (delta / material.delta0_J) * 1.0e3 * q_axis
    return SimpleNamespace(
        js_A_m2=js,
        Te_axis_K=Te_axis,
        delta_axis_J=delta_axis,
        q_axis_m_inv=q_axis,
    )


def test_usadel_poisson_law_accepts_strict_3d_catalog_table(
    small_strip_mesh_bundle,
    gtdgl_material,
    stationary_seed_factory,
):
    mesh, edge_data, ops = small_strip_mesh_bundle
    seed = stationary_seed_factory(mesh, gtdgl_material, q0_m_inv=1.0e7)
    catalog = _strict_3d_catalog(gtdgl_material)

    result = solve_stationary_pytdgl_like(
        mesh=mesh,
        edge_data=edge_data,
        seed=seed,
        material=gtdgl_material,
        ops=ops,
        steps=2,
        dt_s=1.0e-18,
        target_current_A=0.0,
        terminal_psi=0.0,
        adaptive=False,
        n_snapshots=2,
        usadel_catalog=catalog,
        supercurrent_law="usadel_poisson",
    )

    assert result.summary["supercurrent_law"] == "usadel_poisson"
    assert result.summary["usadel_current_available"] is True
    assert "Te" in result.summary["usadel_current_backend"]
    assert "delta" in result.summary["usadel_current_backend"]
    assert "q" in result.summary["usadel_current_backend"]
    assert "edge_js_usadel_snapshot_A_m2" in result.history
    assert result.history["edge_js_usadel_snapshot_A_m2"].shape == (2, ops.n_edges)


def test_usadel_poisson_law_rejects_legacy_delta_q_table(
    small_strip_mesh_bundle,
    gtdgl_material,
    stationary_seed_factory,
):
    mesh, edge_data, ops = small_strip_mesh_bundle
    seed = stationary_seed_factory(mesh, gtdgl_material, q0_m_inv=1.0e7)
    q_axis = np.array([0.0, 1.0e7, 2.0e7, 4.0e7])
    delta_axis = np.array([0.0, 0.5, 0.9, 1.0]) * gtdgl_material.delta0_J
    js = np.outer(delta_axis / gtdgl_material.delta0_J, 1.0e3 * q_axis)
    catalog = SimpleNamespace(js_A_m2=js, q_axis_m_inv=q_axis, delta_axis_J=delta_axis)

    with pytest.raises(RuntimeError, match="3D|Te_axis|supercurrent table"):
        solve_stationary_pytdgl_like(
            mesh=mesh,
            edge_data=edge_data,
            seed=seed,
            material=gtdgl_material,
            ops=ops,
            steps=2,
            dt_s=1.0e-18,
            target_current_A=0.0,
            terminal_psi=0.0,
            adaptive=False,
            n_snapshots=2,
            usadel_catalog=catalog,
            supercurrent_law="usadel_poisson",
        )
