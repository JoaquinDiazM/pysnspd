"""SI current adapters for the pyTDGL-like OE7 backend.

The solver core intentionally follows pyTDGL's CPU/no-screening structure, whose
edge currents live in the same native operator space as ``mu`` and the mesh
coordinates.  This module is the only place where those native edge currents are
converted back to pySNSPD SI diagnostics.

No user-facing bias current is made dimensionless here: ``target_current_A`` is
always interpreted in amperes.  The conversion below follows directly from
``phi = V0 * mu`` and ``x = L0 * x'``:

    j_n = -sigma_n grad(phi) = -(sigma_n V0 / L0) grad'(mu).

The native supercurrent must use the same current-density scale to be consistent
with the pyTDGL-like Poisson projection.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from pysnspd.gtdgl.material import GTDGLMaterial
from pysnspd.gtdgl.operators import (
    FVOperators,
    boundary_currents_from_edge_scalar_least_squares,
    divergence_from_edge_scalar,
    edge_average,
    edge_phase_gradient_from_psi,
    edge_scalar_to_node_vector_least_squares,
    terminal_boundary_accum_A_m,
)
from pysnspd.gtdgl.state import CurrentFields


@dataclass(frozen=True)
class NativeCurrentDiagnostics:
    """Extra diagnostics for native-current SI conversion."""

    current_scale_A_m2: float
    boundary_accum_A_m: np.ndarray
    residual_no_boundary_rms_A_m3: float
    residual_plus_boundary_rms_A_m3: float
    residual_minus_boundary_rms_A_m3: float
    selected_boundary_sign: float
    boundary_currents_from_total_A: dict[str, float]


def native_current_scale_A_m2(device) -> float:
    """Return the SI current-density scale for native pyTDGL-like currents."""

    sigma = float(device.material.sigma_n_S_m)
    V0 = float(device.voltage_scale_V)
    L0 = float(device.length_scale_m)
    return sigma * V0 / max(L0, 1.0e-300)


def native_edge_currents_to_current_fields(
    *,
    psi_dimensionless: np.ndarray,
    native_supercurrent: np.ndarray,
    native_normal_current: np.ndarray,
    device,
    mesh,
    edge_data,
    ops: FVOperators,
    material: GTDGLMaterial,
    Te_K: np.ndarray,
    target_current_A: float,
    boundary_sign: float | None = None,
) -> tuple[CurrentFields, NativeCurrentDiagnostics]:
    """Convert native pyTDGL-like edge currents into a ``CurrentFields`` object.

    The returned currents are SI current densities in A/m^2.  The node-vector
    reconstruction and FV divergence use the same pySNSPD operators used by the
    rest of OE7, but the edge currents themselves come directly from the
    pyTDGL-like solver instead of being re-derived by the older notebook/Usadel
    formula diagnostics.
    """

    psi = np.asarray(psi_dimensionless, dtype=np.complex128)
    js_native = np.asarray(native_supercurrent, dtype=float)
    jn_native = np.asarray(native_normal_current, dtype=float)
    if js_native.shape != (ops.n_edges,):
        raise ValueError(f"native_supercurrent must have shape ({ops.n_edges},), got {js_native.shape}.")
    if jn_native.shape != (ops.n_edges,):
        raise ValueError(f"native_normal_current must have shape ({ops.n_edges},), got {jn_native.shape}.")

    scale = native_current_scale_A_m2(device)
    edge_js = scale * js_native
    edge_jn = scale * jn_native
    edge_jtot = edge_js + edge_jn

    boundary = terminal_boundary_accum_A_m(
        edge_data,
        n_nodes=ops.n_nodes,
        target_current_A=float(target_current_A),
        thickness_m=material.thickness_m,
    )
    div_no = divergence_from_edge_scalar(edge_jtot, ops)
    div_plus = divergence_from_edge_scalar(edge_jtot, ops, boundary_accum_A_m=boundary)
    div_minus = divergence_from_edge_scalar(edge_jtot, ops, boundary_accum_A_m=-boundary)

    rms_no = _rms(div_no)
    rms_plus = _rms(div_plus)
    rms_minus = _rms(div_minus)
    if boundary_sign is None:
        # The native pyTDGL-like Poisson equation is solved in the convention
        #
        #     L_mu mu = div(J_s) + boundary_rhs,
        #
        # while pySNSPD's physical divergence diagnostic expects zero for the
        # closed total-current balance.  Therefore the boundary contribution
        # must be applied with the sign that actually minimizes div(J_s+J_n)
        # in the SI diagnostic.  Keeping this automatic is intentional while we
        # migrate contact assembly: it exposes sign-convention mistakes instead
        # of hiding them in the plotting layer.
        boundary_sign = -1.0 if rms_minus <= rms_plus else +1.0
    boundary_sign = float(boundary_sign)
    div_total = divergence_from_edge_scalar(
        edge_jtot,
        ops,
        boundary_accum_A_m=boundary_sign * boundary,
    )

    node_js_x, node_js_y = edge_scalar_to_node_vector_least_squares(edge_js, ops)
    node_jn_x, node_jn_y = edge_scalar_to_node_vector_least_squares(edge_jn, ops)
    node_jtot_x, node_jtot_y = edge_scalar_to_node_vector_least_squares(edge_jtot, ops)

    edge_Q = edge_phase_gradient_from_psi(psi, ops)
    edge_pb = pairbreaking_ratio_edges(
        Q_edge_m_inv=edge_Q,
        Te_edge_K=edge_average(Te_K, ops),
        material=material,
    )
    node_pb = edge_to_node_weighted_average(edge_pb, ops)

    div_js = divergence_from_edge_scalar(edge_js, ops)
    div_gl = div_js.copy()

    currents = CurrentFields(
        edge_Q_m_inv=edge_Q,
        edge_js_us_A_m2=edge_js,
        edge_js_gl_A_m2=edge_js.copy(),
        edge_jn_A_m2=edge_jn,
        edge_jtot_A_m2=edge_jtot,
        node_div_js_us_A_m3=div_js,
        node_div_js_gl_A_m3=div_gl,
        node_div_jtot_A_m3=div_total,
        node_js_us_x_A_m2=node_js_x,
        node_js_us_y_A_m2=node_js_y,
        node_jn_x_A_m2=node_jn_x,
        node_jn_y_A_m2=node_jn_y,
        node_jtot_x_A_m2=node_jtot_x,
        node_jtot_y_A_m2=node_jtot_y,
        edge_pairbreaking_ratio=edge_pb,
        node_pairbreaking_ratio=node_pb,
    )

    bc_from_total = boundary_currents_from_edge_scalar_least_squares(
        mesh=mesh,
        edge_data=edge_data,
        ops=ops,
        edge_current_i_to_j=edge_jtot,
        thickness_m=material.thickness_m,
    )
    diag = NativeCurrentDiagnostics(
        current_scale_A_m2=float(scale),
        boundary_accum_A_m=boundary,
        residual_no_boundary_rms_A_m3=float(rms_no),
        residual_plus_boundary_rms_A_m3=float(rms_plus),
        residual_minus_boundary_rms_A_m3=float(rms_minus),
        selected_boundary_sign=float(boundary_sign),
        boundary_currents_from_total_A=bc_from_total,
    )
    return currents, diag


def _rms(values: np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        return 0.0
    return float(np.sqrt(np.nanmean(arr * arr)))

# These two helpers were previously imported from the legacy top-level
# ``pysnspd.gtdgl.fields`` module.  They are kept here because they are used by
# the pyTDGL-like current adapter and do not require the obsolete legacy solver.
def pairbreaking_ratio_edges(
    *,
    Q_edge_m_inv: np.ndarray,
    Te_edge_K: np.ndarray,
    material: GTDGLMaterial,
) -> np.ndarray:
    """Return xi^2 Q^2/(1 - T/Tc) on edges."""

    Q = np.asarray(Q_edge_m_inv, dtype=float)
    Te = np.asarray(Te_edge_K, dtype=float)
    a = np.maximum(1.0 - Te / material.Tc_K, 1.0e-30)
    xi2 = material.xi_mod_squared_m2(Te)
    return np.asarray(xi2 * Q * Q / a, dtype=float)


def edge_to_node_weighted_average(edge_values: np.ndarray, ops: FVOperators) -> np.ndarray:
    """Average edge quantities to nodes using dual/length finite-volume weights."""

    values = np.asarray(edge_values, dtype=float)
    if values.shape != (ops.n_edges,):
        raise ValueError(f"edge_values must have shape ({ops.n_edges},), got {values.shape}.")
    weights = ops.dual_face_length_m / np.maximum(ops.edge_length_m, 1.0e-300)
    weights = np.maximum(weights, 1.0e-300)
    out = np.zeros(ops.n_nodes, dtype=float)
    wsum = np.zeros(ops.n_nodes, dtype=float)
    np.add.at(out, ops.edge_i, weights * values)
    np.add.at(out, ops.edge_j, weights * values)
    np.add.at(wsum, ops.edge_i, weights)
    np.add.at(wsum, ops.edge_j, weights)
    return out / np.maximum(wsum, 1.0e-300)

