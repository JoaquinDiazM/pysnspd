"""Adapters between pySNSPD OE7 data and the pyTDGL-like solver core."""
from __future__ import annotations

from types import SimpleNamespace
from typing import Any

import numpy as np

from pysnspd.gtdgl.material import GTDGLMaterial
from pysnspd.mesh.operators import FVOperators, terminal_voltage, edge_scalar_to_node_vector_least_squares
from pysnspd.gtdgl.state import GTDGLStationaryState, RelaxationResult
from pysnspd.solver.diagnostics import (
    current_residual,
    current_density_maxima_A_m2,
    seed_target_current_A,
    target_current_density_A_m2,
)
from pysnspd.gtdgl.currents import native_edge_currents_to_current_fields, native_current_scale_A_m2
from pysnspd.gtdgl.usadel_current import compute_usadel_supercurrent_diagnostic
from pysnspd.gtdgl.allmaras import (
    PhaseDriveContinuationSolver,
    allmaras_coefficients,
    compute_allmaras_appendix_b_diagnostic,
    compute_allmaras_forcing_dimensionless,
    rms as _allmaras_rms,
    max_abs as _allmaras_max_abs,
)
from pysnspd.mesh.device import build_pytdgl_like_device
from .options import SolverOptions, SparseSolver
from .core import TDGLSolver
from pysnspd.thermal.evolution import ThermalRuntimeConfig, ThermalRuntimeController, thermal_stationarity_diagnostics
from .targets import (
    apply_terminal_proximity_seed,
    contact_recovery_diagnostics,
    continuity_diagnostics,
    dynamic_stationarity_diagnostics,
    stationarity_diagnostics,
)

MEV_J = 1.602176634e-22

def _terminal_site_mask_from_device(device, n_nodes: int) -> np.ndarray:
    """Return a boolean mask for metallic normal-terminal sites."""
    mask = np.zeros(int(n_nodes), dtype=bool)
    try:
        terminal_info = device.terminal_info()
    except Exception:
        terminal_info = []
    for terminal in terminal_info:
        idx = np.asarray(getattr(terminal, "site_indices", []), dtype=np.int64)
        idx = idx[(idx >= 0) & (idx < mask.size)]
        if idx.size:
            mask[idx] = True
    return mask


def _terminal_edge_mask_from_device(device, ops: FVOperators) -> np.ndarray:
    """Return edges incident on normal-terminal sites.

    The GL current is automatically zero on these edges when terminal psi is
    clamped to zero.  Usadel-Poisson uses an external constitutive table, so we
    explicitly block the same contact edges to keep the metallic terminal
    condition consistent.
    """
    node_mask = _terminal_site_mask_from_device(device, ops.n_nodes)
    return node_mask[np.asarray(ops.edge_i, dtype=np.int64)] | node_mask[np.asarray(ops.edge_j, dtype=np.int64)]



def _normalize_supercurrent_law(value: str) -> str:
    law = str(value).strip().lower().replace("-", "_")
    aliases = {
        "gl": "gl",
        "pytdgl": "gl",
        "native_gl": "gl",
        "usadel": "usadel_poisson",
        "usadel_poisson": "usadel_poisson",
        "poisson_usadel": "usadel_poisson",
    }
    if law not in aliases:
        raise ValueError(
            "supercurrent_law must be one of gl or usadel_poisson "
            f"(got {value!r})."
        )
    return aliases[law]


def _build_usadel_poisson_supercurrent_override(
    *,
    usadel_catalog: Any | None,
    device,
    material: GTDGLMaterial,
    Te_K: np.ndarray,
    ops: FVOperators,
):
    scale = native_current_scale_A_m2(device)
    blocked_edge_mask = _terminal_edge_mask_from_device(device, ops)

    def usadel_poisson_supercurrent(psi_dimensionless: np.ndarray, gl_supercurrent_native: np.ndarray) -> np.ndarray:
        del gl_supercurrent_native
        diag = compute_usadel_supercurrent_diagnostic(
            usadel_catalog=usadel_catalog,
            psi_dimensionless=psi_dimensionless,
            material=material,
            Te_K=Te_K,
            ops=ops,
            blocked_edge_mask=blocked_edge_mask,
        )
        if not diag.available:
            raise RuntimeError(
                "--ss-supercurrent-law usadel-poisson requires a PRE Usadel "
                f"supercurrent table. Diagnostic reason: {diag.reason}"
            )
        return np.asarray(diag.edge_js_usadel_A_m2, dtype=float) / max(scale, 1.0e-300)

    return usadel_poisson_supercurrent


def _build_allmaras_forcing_callback(
    *,
    usadel_catalog: Any | None,
    device,
    material: GTDGLMaterial,
    Te_K: np.ndarray,
    ops: FVOperators,
    blocked_edge_mask: np.ndarray,
    require_usadel: bool,
    phase_drive_continuation: PhaseDriveContinuationSolver,
):
    """Build the Appendix-B forcing callback used by ``TDGLSolver``.

    The callback is explicit in the current order parameter.  For the official
    ``usadel_poisson`` path it uses the PRE Matsubara/Usadel supercurrent table
    for both Poisson and the Allmaras current-divergence correction.
    """

    Te = np.asarray(Te_K, dtype=float)
    L0 = float(device.length_scale_m)

    def callback(psi_dimensionless: np.ndarray, psi_laplacian) -> np.ndarray:
        psi = np.asarray(psi_dimensionless, dtype=np.complex128)
        edge_js = None
        if require_usadel:
            diag = compute_usadel_supercurrent_diagnostic(
                usadel_catalog=usadel_catalog,
                psi_dimensionless=psi,
                material=material,
                Te_K=Te,
                ops=ops,
                blocked_edge_mask=blocked_edge_mask,
            )
            if not diag.available:
                raise RuntimeError(
                    "Appendix-B Allmaras update with usadel_poisson requires a PRE "
                    f"Matsubara supercurrent table. Diagnostic reason: {diag.reason}"
                )
            edge_js = diag.edge_js_usadel_A_m2

        forcing = compute_allmaras_forcing_dimensionless(
            psi_dimensionless=psi,
            psi_laplacian_dimensionless=psi_laplacian @ psi,
            material=material,
            Te_K=Te,
            ops=ops,
            length_scale_m=L0,
            edge_js_usadel_A_m2=edge_js,
            blocked_edge_mask=blocked_edge_mask,
            phase_drive_continuation=phase_drive_continuation,
        )
        info = forcing.phase_drive_convergence
        callback.last_convergence_diagnostics = {
            "converged": bool(info.converged),
            "iterations": int(info.iterations),
            "residual_rel": float(info.residual_rel),
            "direct_node_count": int(info.direct_node_count),
            "continued_node_count": int(info.continued_node_count),
            "zero_amplitude_node_count": int(info.zero_amplitude_node_count),
        }
        return np.asarray(forcing.forcing_dimensionless, dtype=np.complex128)

    callback.last_convergence_diagnostics = {}
    return callback
