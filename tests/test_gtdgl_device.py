"""Device adapter tests for the flattened gTDGL backend."""
from __future__ import annotations

import numpy as np

from pysnspd.mesh.device import build_pytdgl_like_device
from pysnspd.solver.options import SolverOptions
from pysnspd.solver.core import TDGLSolver


def test_device_adapter_uses_si_current_conversion(small_strip_mesh_bundle, gtdgl_material):
    mesh, edge_data, ops = small_strip_mesh_bundle
    Te = np.full(mesh.n_nodes, 0.9)
    device = build_pytdgl_like_device(
        mesh=mesh,
        edge_data=edge_data,
        material=gtdgl_material,
        ops=ops,
        Te_K=Te,
        target_current_A=3.5e-6,
    )
    assert device.length_scale_m > 0.0
    assert device.voltage_scale_V > 0.0
    assert device.current_scale_A > 0.0
    assert device.terminal_neumann_current_unit_A > 0.0
    terms = {t.name: t for t in device.terminal_info()}
    assert set(terms) == {"left", "right"}
    assert terms["left"].length > 0.0
    assert terms["right"].length > 0.0

    left_mu = device.terminal_mu_boundary_value(
        terminal=terms["left"],
        terminal_currents_A={"left": -3.5e-6, "right": 3.5e-6},
    )
    right_mu = device.terminal_mu_boundary_value(
        terminal=terms["right"],
        terminal_currents_A={"left": -3.5e-6, "right": 3.5e-6},
    )
    assert np.isfinite(left_mu)
    assert np.isfinite(right_mu)
    assert np.sign(left_mu) == -np.sign(right_mu)


def test_tdgl_solver_zero_current_smoke(small_strip_mesh_bundle, gtdgl_material):
    mesh, edge_data, ops = small_strip_mesh_bundle
    Te = np.full(mesh.n_nodes, 0.9)
    device = build_pytdgl_like_device(
        mesh=mesh,
        edge_data=edge_data,
        material=gtdgl_material,
        ops=ops,
        Te_K=Te,
        target_current_A=0.0,
    )
    options = SolverOptions(
        solve_time=2.0e-4,
        dt_init=1.0e-4,
        dt_max=1.0e-4,
        adaptive=False,
        terminal_psi=None,
    )
    solver = TDGLSolver(device, options, terminal_currents={"left": 0.0, "right": 0.0})
    sol = solver.solve()
    assert sol is not None
    assert np.all(np.isfinite(sol.tdgl_data.psi))
    assert np.all(np.isfinite(sol.tdgl_data.mu))
