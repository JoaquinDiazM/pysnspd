"""Shared state containers for OE7 stationary gTDGL relaxation."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

MEV_J = 1.602176634e-22


@dataclass
class CurrentFields:
    """Edge and node current diagnostics."""

    edge_Q_m_inv: np.ndarray
    edge_js_us_A_m2: np.ndarray
    edge_js_gl_A_m2: np.ndarray
    edge_jn_A_m2: np.ndarray
    edge_jtot_A_m2: np.ndarray
    node_div_js_us_A_m3: np.ndarray
    node_div_js_gl_A_m3: np.ndarray
    node_div_jtot_A_m3: np.ndarray
    node_js_us_x_A_m2: np.ndarray
    node_js_us_y_A_m2: np.ndarray
    node_jn_x_A_m2: np.ndarray
    node_jn_y_A_m2: np.ndarray
    node_jtot_x_A_m2: np.ndarray
    node_jtot_y_A_m2: np.ndarray
    edge_pairbreaking_ratio: np.ndarray
    node_pairbreaking_ratio: np.ndarray


@dataclass
class FormulaFields:
    """Notebook formula-field bundle evaluated from Delta, varphi and Te."""

    Te_K: np.ndarray
    rho: np.ndarray
    tau_sc_s: np.ndarray
    alpha_kwt_J_inv2: np.ndarray
    xi_mod_m: np.ndarray
    delta_mod_J: np.ndarray
    delta_abs_J: np.ndarray
    edge_Q_m_inv: np.ndarray
    edge_js_us_A_m2: np.ndarray
    edge_js_gl_A_m2: np.ndarray
    node_div_js_us_A_m3: np.ndarray
    node_div_js_gl_A_m3: np.ndarray
    node_div_correction_A_m3: np.ndarray
    node_lap_delta_J_m2: np.ndarray
    forcing_J: np.ndarray


@dataclass
class PoissonResult:
    """Poisson projection result, following the notebook naming."""

    phi_V: np.ndarray
    edge_jn_A_m2: np.ndarray
    edge_jtot_A_m2: np.ndarray
    node_div_jtot_A_m3: np.ndarray
    node_div_js_A_m3: np.ndarray
    node_div_jn_A_m3: np.ndarray
    lambda_mean: float


@dataclass
class StepInfo:
    """One accepted KWT/Poisson step diagnostic."""

    dt_eff_s: float
    retries: int
    discr_min: float
    max_amp2_change_rel: float
    max_Q_m_inv: float
    p95_Q_m_inv: float
    median_Q_m_inv: float
    max_js_A_m2: float
    max_j_A_m2: float
    rms_divj_rel: float
    max_divj_rel: float
    delta_min_norm: float
    delta_max_norm: float
    phi_ptp_V: float


@dataclass
class GTDGLStationaryState:
    """Node-based stationary gTDGL state."""

    psi_J: np.ndarray
    phi_V: np.ndarray
    Te_K: np.ndarray
    Tph_K: np.ndarray
    currents: CurrentFields
    metadata: dict[str, Any]


@dataclass
class RelaxationResult:
    """Final state and compact history for one stationary relaxation run."""

    state: GTDGLStationaryState
    history: dict[str, np.ndarray]
    summary: dict[str, Any]


@dataclass
class _PoissonOperator:
    A_aug: Any
    solver: Any


@dataclass
class PhiBoundaryConditions:
    """Discrete electric boundary constraints for varphi.

    Each row enforces, on the first inward normal edge b -> k,

        phi_b - phi_k = ell/sigma_n * (j_target - j_s,bk).

    Here ``edge_sign_b_to_inner`` converts the solver edge orientation into the
    boundary-to-inner orientation used by the constraint.
    """

    boundary_nodes: np.ndarray
    inner_nodes: np.ndarray
    edge_index: np.ndarray
    edge_sign_b_to_inner: np.ndarray
    target_edge_A_m2: np.ndarray
    edge_length_m: np.ndarray
    boundary_mask: np.ndarray

