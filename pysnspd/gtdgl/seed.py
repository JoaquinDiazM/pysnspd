"""OE6 stationary analytic seed for the gTDGL/Poisson solver.

OE6 does not evolve the gTDGL equation yet. It builds a physically consistent
initial condition for the later stationary relaxation stage:

    Te(r)      = T_bias,
    Tph(r)     = T_bias,
    Delta(r)   = Delta_eq(T_bias, q_bias),
    theta(r)   = q_bias * (x - x0),
    Q(r)       = grad(theta) = q_bias xhat,
    phi_el(r)  = 0,
    j_n(r)     = 0,
    j_s(r)     = j_Usadel(q_bias,T_bias) xhat.

The bias q is selected from the stable branch of the OE3 Usadel calibration
sweep:

    I_s(q_bias,T_bias) = I_bias,
    q_bias <= q_c,

where q_c is the q value at the maximum Usadel current. This avoids using the
post-critical branch, where the same current can correspond to a different
superfluid momentum.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np


HBAR_J_S = 1.054571817e-34
MEV_J = 1.602176634e-22


@dataclass(frozen=True)
class BiasUsadelState:
    """Selected stable-branch Usadel state for the imposed bias current."""

    I_bias_A: float
    Ic_A: float
    I_bias_over_Ic: float
    q_bias_m_inv: float
    q_critical_m_inv: float
    gamma_bias_J: float
    gamma_bias_meV: float
    delta_bias_J: float
    delta_bias_meV: float
    current_density_bias_A_m2: float
    branch_policy: str
    interpolation_policy: str


@dataclass(frozen=True)
class StationarySeed:
    """Node-based OE6 stationary seed fields."""

    node_Te_K: np.ndarray
    node_Tph_K: np.ndarray
    node_delta_J: np.ndarray
    node_R_J: np.ndarray
    node_R_normalized: np.ndarray
    node_theta_rad: np.ndarray
    node_psi_real_J: np.ndarray
    node_psi_imag_J: np.ndarray
    node_phi_electric_V: np.ndarray
    node_qx_m_inv: np.ndarray
    node_qy_m_inv: np.ndarray
    node_js_x_A_m2: np.ndarray
    node_js_y_A_m2: np.ndarray
    node_jn_x_A_m2: np.ndarray
    node_jn_y_A_m2: np.ndarray
    node_jtot_x_A_m2: np.ndarray
    node_jtot_y_A_m2: np.ndarray
    node_div_j_A_m3: np.ndarray
    metadata: dict[str, Any]


def select_bias_state_from_usadel(
    usadel_catalog,
    *,
    I_bias_A: float | None = None,
    relative_overbias_tol: float = 1.0e-9,
) -> BiasUsadelState:
    """Select q_bias from the stable branch of the Usadel current sweep.

    The calibration sweep is treated as the authoritative map

        q -> I_s(q,T_bias).

    Since I_s(q) is not globally one-to-one, only points up to q_c are used.
    """
    q = np.asarray(usadel_catalog.calibration_q_values_m_inv, dtype=float)
    current = np.asarray(usadel_catalog.calibration_current_values_A, dtype=float)
    current_density = np.asarray(
        usadel_catalog.calibration_current_density_values_A_m2,
        dtype=float,
    )
    delta = np.asarray(usadel_catalog.calibration_delta_eq_values_J, dtype=float)
    gamma = np.asarray(usadel_catalog.calibration_gamma_values_J, dtype=float)

    finite = (
        np.isfinite(q)
        & np.isfinite(current)
        & np.isfinite(current_density)
        & np.isfinite(delta)
        & np.isfinite(gamma)
    )
    finite &= q >= 0.0
    finite &= current >= 0.0
    finite &= delta >= 0.0
    finite &= gamma >= 0.0

    if np.count_nonzero(finite) < 3:
        raise ValueError("Usadel calibration sweep has fewer than 3 valid points.")

    q = q[finite]
    current = current[finite]
    current_density = current_density[finite]
    delta = delta[finite]
    gamma = gamma[finite]

    order_q = np.argsort(q)
    q = q[order_q]
    current = current[order_q]
    current_density = current_density[order_q]
    delta = delta[order_q]
    gamma = gamma[order_q]

    idx_ic = int(np.argmax(current))
    Ic_A = float(current[idx_ic])
    q_c = float(q[idx_ic])

    if Ic_A <= 0.0:
        raise ValueError("Usadel calibration critical current is not positive.")

    if I_bias_A is None:
        I_bias_A = float(usadel_catalog.metadata.get("I_bias_A", np.nan))
    I_bias_A = float(I_bias_A)

    if not np.isfinite(I_bias_A):
        raise ValueError("I_bias_A is not finite.")
    if I_bias_A < 0.0:
        raise ValueError("I_bias_A must be non-negative.")
    if I_bias_A > Ic_A * (1.0 + relative_overbias_tol):
        raise ValueError(
            f"I_bias_A={I_bias_A:.6e} A exceeds Usadel Ic={Ic_A:.6e} A. "
            "OE6-v1 only builds superconducting stationary seeds below Ic."
        )

    target_I = min(I_bias_A, Ic_A)

    q_stable = q[: idx_ic + 1]
    current_stable = current[: idx_ic + 1]
    current_density_stable = current_density[: idx_ic + 1]
    delta_stable = delta[: idx_ic + 1]
    gamma_stable = gamma[: idx_ic + 1]

    # Interpolate using current as coordinate. Numerical sweeps should be
    # monotone on the stable branch; unique sorting makes this robust to tiny
    # repeated values.
    order_I = np.argsort(current_stable)
    I_sorted = current_stable[order_I]
    q_sorted_by_I = q_stable[order_I]

    I_unique, unique_idx = np.unique(I_sorted, return_index=True)
    q_unique_by_I = q_sorted_by_I[unique_idx]

    if I_unique.size < 2:
        raise ValueError("Stable Usadel branch is not usable for interpolation.")

    q_bias = float(np.interp(target_I, I_unique, q_unique_by_I))

    delta_bias = float(np.interp(q_bias, q_stable, delta_stable))
    gamma_bias = float(np.interp(q_bias, q_stable, gamma_stable))
    j_bias = float(np.interp(q_bias, q_stable, current_density_stable))

    return BiasUsadelState(
        I_bias_A=I_bias_A,
        Ic_A=Ic_A,
        I_bias_over_Ic=float(I_bias_A / Ic_A),
        q_bias_m_inv=q_bias,
        q_critical_m_inv=q_c,
        gamma_bias_J=gamma_bias,
        gamma_bias_meV=float(gamma_bias / MEV_J),
        delta_bias_J=delta_bias,
        delta_bias_meV=float(delta_bias / MEV_J),
        current_density_bias_A_m2=j_bias,
        branch_policy="stable_usadel_branch_q_le_qc",
        interpolation_policy="linear_interpolation_on_stable_I_s_q_branch",
    )


def build_stationary_seed(
    *,
    mesh,
    edge_data,
    usadel_catalog,
    I_bias_A: float | None = None,
    T_bias_K: float | None = None,
    phase_origin: str = "center",
) -> StationarySeed:
    """Build the OE6 analytic stationary seed on mesh nodes."""
    nodes = np.asarray(mesh.nodes, dtype=float)
    n_nodes = int(nodes.shape[0])

    if T_bias_K is None:
        T_bias_K = float(usadel_catalog.metadata["T_bias_K"])
    T_bias_K = float(T_bias_K)

    if T_bias_K <= 0.0:
        raise ValueError("T_bias_K must be positive.")

    state = select_bias_state_from_usadel(usadel_catalog, I_bias_A=I_bias_A)

    delta0_J = float(usadel_catalog.metadata.get("delta0_J", np.nan))
    if not np.isfinite(delta0_J) or delta0_J <= 0.0:
        delta0_J = float(np.max(usadel_catalog.delta_values_J))

    x = nodes[:, 0]
    if phase_origin == "center":
        x0 = 0.5 * (float(np.min(x)) + float(np.max(x)))
    elif phase_origin == "left":
        x0 = float(np.min(x))
    else:
        raise ValueError("phase_origin must be 'center' or 'left'.")

    theta = state.q_bias_m_inv * (x - x0)

    delta = np.full(n_nodes, state.delta_bias_J, dtype=float)
    R = delta.copy()
    R_norm = R / delta0_J if delta0_J > 0.0 else np.zeros_like(R)

    psi_real = R * np.cos(theta)
    psi_imag = R * np.sin(theta)

    qx = np.full(n_nodes, state.q_bias_m_inv, dtype=float)
    qy = np.zeros(n_nodes, dtype=float)

    js_x = np.full(n_nodes, state.current_density_bias_A_m2, dtype=float)
    js_y = np.zeros(n_nodes, dtype=float)

    jn_x = np.zeros(n_nodes, dtype=float)
    jn_y = np.zeros(n_nodes, dtype=float)

    jtot_x = js_x + jn_x
    jtot_y = js_y + jn_y

    phi = np.zeros(n_nodes, dtype=float)

    # For this first analytic seed the current is exactly uniform, so the
    # continuum divergence is zero. A real finite-volume divergence diagnostic
    # will belong to OE7, where the current becomes spatially nonuniform.
    div_j = np.zeros(n_nodes, dtype=float)

    boundary = compute_boundary_currents(
        edge_data=edge_data,
        jx_A_m2=state.current_density_bias_A_m2,
        jy_A_m2=0.0,
        thickness_m=float(usadel_catalog.metadata["thickness_m"]),
    )

    voltage = compute_terminal_voltage(
        nodes=nodes,
        phi_electric_V=phi,
        length_m=float(mesh.length_m),
    )

    I_right = boundary["right_A"]
    I_left = boundary["left_A"]
    I_target = state.I_bias_A

    metadata = {
        "backend": "oe6_stationary_analytic_seed_v1",
        "description": (
            "Analytic stationary seed for later gTDGL relaxation. No gTDGL time "
            "evolution is performed in OE6-v1."
        ),
        "T_bias_K": T_bias_K,
        "Tph_bias_K": T_bias_K,
        "phase_origin": phase_origin,
        "amplitude_convention": "dimensional_delta_J",
        "normal_current_seed": "zero",
        "electric_potential_seed": "zero",
        "phase_seed": "linear_longitudinal_theta=q_bias*(x-x0)",
        "gauge_policy": "A=0, Q=grad(theta)",
        "continuity_policy": (
            "Uniform analytic seed has continuum div(j)=0. Discrete Poisson "
            "continuity enforcement begins in OE7."
        ),
        "bias_state": state.__dict__,
        "boundary_currents_A": boundary,
        "terminal_voltage_V": voltage,
        "right_current_error_rel": _relative_error(I_right, I_target),
        "left_current_error_rel": _relative_error(-I_left, I_target),
        "target_current_A": I_target,
        "integrated_right_current_A": I_right,
        "integrated_left_current_A": I_left,
        "phase_span_rad": float(np.max(theta) - np.min(theta)),
        "delta0_J": delta0_J,
        "delta0_meV": float(delta0_J / MEV_J),
        "seed_checks": {
            "delta_is_uniform": True,
            "phi_is_zero": True,
            "normal_current_is_zero": True,
            "divergence_is_zero_analytic": True,
        },
    }

    return StationarySeed(
        node_Te_K=np.full(n_nodes, T_bias_K, dtype=float),
        node_Tph_K=np.full(n_nodes, T_bias_K, dtype=float),
        node_delta_J=delta,
        node_R_J=R,
        node_R_normalized=R_norm,
        node_theta_rad=theta,
        node_psi_real_J=psi_real,
        node_psi_imag_J=psi_imag,
        node_phi_electric_V=phi,
        node_qx_m_inv=qx,
        node_qy_m_inv=qy,
        node_js_x_A_m2=js_x,
        node_js_y_A_m2=js_y,
        node_jn_x_A_m2=jn_x,
        node_jn_y_A_m2=jn_y,
        node_jtot_x_A_m2=jtot_x,
        node_jtot_y_A_m2=jtot_y,
        node_div_j_A_m3=div_j,
        metadata=metadata,
    )


def compute_boundary_currents(
    *,
    edge_data,
    jx_A_m2: float,
    jy_A_m2: float,
    thickness_m: float,
) -> dict[str, float]:
    """Integrate a uniform current density through tagged boundaries."""
    tags = np.asarray(edge_data.tags).astype(str)
    lengths = np.asarray(edge_data.lengths, dtype=float)

    normals = {
        "left": (-1.0, 0.0),
        "right": (1.0, 0.0),
        "bottom": (0.0, -1.0),
        "top": (0.0, 1.0),
    }

    out: dict[str, float] = {}
    for tag, normal in normals.items():
        mask = tags == tag
        nx, ny = normal
        flux_density = float(jx_A_m2) * nx + float(jy_A_m2) * ny
        out[f"{tag}_A"] = float(thickness_m * np.sum(lengths[mask] * flux_density))

    out["net_boundary_current_A"] = float(sum(out.values()))
    return out


def compute_terminal_voltage(
    *,
    nodes: np.ndarray,
    phi_electric_V: np.ndarray,
    length_m: float,
) -> float:
    """Return <phi>_right - <phi>_left."""
    nodes = np.asarray(nodes, dtype=float)
    phi = np.asarray(phi_electric_V, dtype=float)

    x = nodes[:, 0]
    tol = max(1.0e-15, 1.0e-9 * float(length_m))

    left = np.abs(x - np.min(x)) <= tol
    right = np.abs(x - np.max(x)) <= tol

    if not np.any(left) or not np.any(right):
        return float("nan")

    return float(np.mean(phi[right]) - np.mean(phi[left]))


def save_stationary_seed_npz(seed: StationarySeed, path: str | Path) -> Path:
    """Save a StationarySeed to NPZ."""
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)

    np.savez_compressed(
        output,
        node_Te_K=seed.node_Te_K,
        node_Tph_K=seed.node_Tph_K,
        node_delta_J=seed.node_delta_J,
        node_R_J=seed.node_R_J,
        node_R_normalized=seed.node_R_normalized,
        node_theta_rad=seed.node_theta_rad,
        node_psi_real_J=seed.node_psi_real_J,
        node_psi_imag_J=seed.node_psi_imag_J,
        node_phi_electric_V=seed.node_phi_electric_V,
        node_qx_m_inv=seed.node_qx_m_inv,
        node_qy_m_inv=seed.node_qy_m_inv,
        node_js_x_A_m2=seed.node_js_x_A_m2,
        node_js_y_A_m2=seed.node_js_y_A_m2,
        node_jn_x_A_m2=seed.node_jn_x_A_m2,
        node_jn_y_A_m2=seed.node_jn_y_A_m2,
        node_jtot_x_A_m2=seed.node_jtot_x_A_m2,
        node_jtot_y_A_m2=seed.node_jtot_y_A_m2,
        node_div_j_A_m3=seed.node_div_j_A_m3,
        metadata=np.array(seed.metadata, dtype=object),
    )
    return output


def seed_summary(seed: StationarySeed) -> dict[str, Any]:
    """Build a compact YAML/console summary for the seed."""
    meta = seed.metadata
    bias = meta["bias_state"]
    boundary = meta["boundary_currents_A"]

    return {
        "backend": meta["backend"],
        "n_nodes": int(seed.node_Te_K.size),
        "T_bias_K": float(meta["T_bias_K"]),
        "I_bias_A": float(bias["I_bias_A"]),
        "Ic_A": float(bias["Ic_A"]),
        "I_bias_over_Ic": float(bias["I_bias_over_Ic"]),
        "q_bias_m_inv": float(bias["q_bias_m_inv"]),
        "q_critical_m_inv": float(bias["q_critical_m_inv"]),
        "gamma_bias_meV": float(bias["gamma_bias_meV"]),
        "delta_bias_meV": float(bias["delta_bias_meV"]),
        "current_density_bias_A_m2": float(bias["current_density_bias_A_m2"]),
        "phase_span_rad": float(meta["phase_span_rad"]),
        "terminal_voltage_V": float(meta["terminal_voltage_V"]),
        "integrated_left_current_A": float(boundary["left_A"]),
        "integrated_right_current_A": float(boundary["right_A"]),
        "net_boundary_current_A": float(boundary["net_boundary_current_A"]),
        "left_current_error_rel": float(meta["left_current_error_rel"]),
        "right_current_error_rel": float(meta["right_current_error_rel"]),
        "branch_policy": str(bias["branch_policy"]),
        "phase_seed": str(meta["phase_seed"]),
        "electric_potential_seed": str(meta["electric_potential_seed"]),
        "normal_current_seed": str(meta["normal_current_seed"]),
        "divergence_rms_A_m3": float(np.sqrt(np.mean(seed.node_div_j_A_m3**2))),
        "delta_min_meV": float(np.min(seed.node_delta_J) / MEV_J),
        "delta_max_meV": float(np.max(seed.node_delta_J) / MEV_J),
        "theta_min_rad": float(np.min(seed.node_theta_rad)),
        "theta_max_rad": float(np.max(seed.node_theta_rad)),
        "js_x_min_A_m2": float(np.min(seed.node_js_x_A_m2)),
        "js_x_max_A_m2": float(np.max(seed.node_js_x_A_m2)),
    }


def _relative_error(value: float, reference: float) -> float:
    if abs(reference) <= 0.0:
        return 0.0 if abs(value) <= 0.0 else float("inf")
    return float(abs(value - reference) / abs(reference))