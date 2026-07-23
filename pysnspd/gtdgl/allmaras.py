"""Appendix-B Allmaras diagnostics for the pyTDGL-like backend.

This module implements Appendix-B coefficient fields, the
current-divergence correction, and the dimensionless forcing injected into the
promoted flat gTDGL solver.  The physical Usadel/Matsubara current catalogue is
owned by ``usadel_current.py``; when available, its edge current is supplied to
this module so the Allmaras correction uses the same supercurrent closure as the
Poisson projection.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from pysnspd.gtdgl.material import GTDGLMaterial, K_B_J_K
from pysnspd.mesh.operators import (
    FVOperators,
    divergence_from_edge_scalar,
    edge_average,
    edge_phase_gradient_from_psi,
    edge_scalar_to_node_vector_least_squares,
)

HBAR_J_S = 1.054571817e-34
E_ABS_C = 1.602176634e-19
MEV_J = 1.602176634e-22


@dataclass(frozen=True)
class AllmarasCoefficientFields:
    """Node-local Appendix-B material coefficients."""

    tau_sc_s: np.ndarray
    rho_kwt: np.ndarray
    gamma_kwt_dimensionless: np.ndarray
    alpha_kwt_inv_J2: np.ndarray
    xi_mod2_m2: np.ndarray
    delta_mod2_J2: np.ndarray
    delta_mod_over_delta0: np.ndarray
    solver_epsilon: np.ndarray
    correction_C_J2_m3_A: np.ndarray


@dataclass(frozen=True)
class PhaseDriveConvergenceInfo:
    """Convergence metadata for the low-amplitude phase-drive continuation."""

    converged: bool
    iterations: int
    residual_rel: float
    direct_node_count: int
    continued_node_count: int
    zero_amplitude_node_count: int
    direct_amplitude_fraction: float
    tolerance: float
    max_iterations: int


@dataclass(frozen=True)
class AllmarasDiagnosticFields:
    """Snapshot diagnostic fields for the Appendix-B correction."""

    coefficients: AllmarasCoefficientFields
    edge_Q_m_inv: np.ndarray
    edge_R_J: np.ndarray
    edge_Te_K: np.ndarray
    edge_js_us_allmaras_A_m2: np.ndarray
    edge_js_gl_allmaras_A_m2: np.ndarray
    node_js_us_allmaras_x_A_m2: np.ndarray
    node_js_us_allmaras_y_A_m2: np.ndarray
    node_js_gl_allmaras_x_A_m2: np.ndarray
    node_js_gl_allmaras_y_A_m2: np.ndarray
    node_div_js_us_allmaras_A_m3: np.ndarray
    node_div_js_gl_allmaras_A_m3: np.ndarray
    node_mismatch_divergence_A_m3: np.ndarray
    node_phase_drive_complex_J: np.ndarray
    node_phase_drive_abs_J: np.ndarray
    node_phase_drive_abs_over_delta0: np.ndarray
    bulk_node_mask: np.ndarray
    phase_drive_convergence: PhaseDriveConvergenceInfo



@dataclass(frozen=True)
class AllmarasForcingFields:
    """Appendix-B forcing fields used by the local KWT update.

    ``forcing_dimensionless`` is the complete bracket in Eq. (139), divided by
    ``Delta0`` and written in the dimensionless coordinates of the promoted
    pyTDGL-like solver.  The local solver then multiplies it by
    ``rho_KWT * dt/tau0`` through the existing ``w_i, z_i`` algebra.
    """

    coefficients: AllmarasCoefficientFields
    forcing_dimensionless: np.ndarray
    diffusion_dimensionless: np.ndarray
    reaction_dimensionless: np.ndarray
    phase_drive_dimensionless: np.ndarray
    edge_Q_m_inv: np.ndarray
    edge_R_J: np.ndarray
    edge_Te_K: np.ndarray
    edge_js_us_A_m2: np.ndarray
    edge_js_gl_A_m2: np.ndarray
    node_div_js_us_A_m3: np.ndarray
    node_div_js_gl_A_m3: np.ndarray
    node_mismatch_divergence_A_m3: np.ndarray
    phase_drive_convergence: PhaseDriveConvergenceInfo


def allmaras_coefficients(
    *,
    psi_dimensionless: np.ndarray,
    material: GTDGLMaterial,
    Te_K: np.ndarray,
) -> AllmarasCoefficientFields:
    """Return Appendix-B local coefficient fields in SI units.

    The solver currently advances ``psi = Delta / Delta0``.  Therefore the KWT
    gamma used by the pyTDGL-like local algebra is

        gamma_i = 2 Delta0 tau_sc(Te_i) / hbar,

    so that ``sqrt(1 + gamma_i^2 |psi_i|^2)`` equals
    ``rho_KWT(Te_i, |Delta_i|)``.
    """

    psi = np.asarray(psi_dimensionless, dtype=np.complex128)
    Te = np.asarray(Te_K, dtype=float).reshape(-1)
    if psi.shape != Te.shape:
        raise ValueError(f"psi and Te must have the same node shape, got {psi.shape} and {Te.shape}.")
    Te_safe = np.maximum(Te, 1.0e-12)
    R_J = np.abs(psi) * float(material.delta0_J)

    tau_sc = np.asarray(material.tau_sc_s(Te_safe), dtype=float)
    rho = np.asarray(material.rho_kwt(Te_safe, np.maximum(R_J, 0.0)), dtype=float)
    gamma = 2.0 * float(material.delta0_J) * tau_sc / HBAR_J_S
    alpha = 2.0 * tau_sc * tau_sc / (HBAR_J_S * HBAR_J_S)

    delta_mod2 = np.asarray(material.delta_mod_squared_J2(Te_safe), dtype=float)
    delta0 = float(material.delta0_J)
    solver_epsilon = np.clip(delta_mod2 / max(delta0 * delta0, 1.0e-300), 0.0, 1.0)
    delta_mod_over_delta0 = np.sqrt(np.maximum(delta_mod2, 0.0)) / max(delta0, 1.0e-300)

    xi_mod2 = (
        np.pi
        * HBAR_J_S
        * float(material.D_m2_s)
        / (
            4.0
            * np.sqrt(2.0)
            * K_B_J_K
            * float(material.Tc_K)
            * np.sqrt(np.maximum(1.0 + Te_safe / float(material.Tc_K), 1.0e-300))
        )
    )

    correction_C = (
        HBAR_J_S
        * E_ABS_C
        * float(material.D_m2_s)
        / (
            float(material.sigma_n_S_m)
            * np.sqrt(2.0)
            * np.sqrt(np.maximum(1.0 + Te_safe / float(material.Tc_K), 1.0e-300))
        )
    )

    return AllmarasCoefficientFields(
        tau_sc_s=tau_sc,
        rho_kwt=rho,
        gamma_kwt_dimensionless=gamma,
        alpha_kwt_inv_J2=alpha,
        xi_mod2_m2=np.asarray(xi_mod2, dtype=float),
        delta_mod2_J2=delta_mod2,
        delta_mod_over_delta0=delta_mod_over_delta0,
        solver_epsilon=solver_epsilon,
        correction_C_J2_m3_A=np.asarray(correction_C, dtype=float),
    )


def compute_allmaras_appendix_b_diagnostic(
    *,
    psi_dimensionless: np.ndarray,
    material: GTDGLMaterial,
    Te_K: np.ndarray,
    ops: FVOperators,
    terminal_node_mask: np.ndarray | None = None,
    blocked_edge_mask: np.ndarray | None = None,
    edge_js_usadel_A_m2: np.ndarray | None = None,
    phase_drive_continuation: PhaseDriveContinuationSolver | None = None,
) -> AllmarasDiagnosticFields:
    """Compute the Appendix-B current-divergence correction as a diagnostic.

    This evaluates the Appendix-B current correction on the current snapshot
    using the same normalized harmonic-continuation rule as the ``w_i, z_i``
    update.
    """

    psi = np.asarray(psi_dimensionless, dtype=np.complex128)
    Te = np.asarray(Te_K, dtype=float).reshape(-1)
    coeff = allmaras_coefficients(psi_dimensionless=psi, material=material, Te_K=Te)

    Q = edge_phase_gradient_from_psi(psi, ops)
    R_node_J = np.abs(psi) * float(material.delta0_J)
    R_edge_J = edge_average(R_node_J, ops)
    Te_edge = np.maximum(edge_average(Te, ops), 1.0e-12)

    if edge_js_usadel_A_m2 is None:
        # Analytic Allmaras fallback retained for diagnostics and tests.  The
        # production usadel_poisson path supplies the Matsubara-table current
        # through ``edge_js_usadel_A_m2``.
        js_us = analytic_allmaras_usadel_current_edges(
            edge_Q_m_inv=Q,
            edge_R_J=R_edge_J,
            edge_Te_K=Te_edge,
            material=material,
        )
    else:
        js_us = np.asarray(edge_js_usadel_A_m2, dtype=float).reshape(-1)
        if js_us.shape != (ops.n_edges,):
            raise ValueError(f"edge_js_usadel_A_m2 must have shape ({ops.n_edges},), got {js_us.shape}.")
    js_gl = gl_supercurrent_edges(
        edge_Q_m_inv=Q,
        edge_R_J=R_edge_J,
        material=material,
    )

    if blocked_edge_mask is not None:
        mask = np.asarray(blocked_edge_mask, dtype=bool).reshape(-1)
        if mask.size != ops.n_edges:
            raise ValueError(f"blocked_edge_mask has length {mask.size}, expected {ops.n_edges}.")
        if np.any(mask):
            js_us = js_us.copy()
            js_gl = js_gl.copy()
            js_us[mask] = 0.0
            js_gl[mask] = 0.0

    js_us_x, js_us_y = edge_scalar_to_node_vector_least_squares(js_us, ops)
    js_gl_x, js_gl_y = edge_scalar_to_node_vector_least_squares(js_gl, ops)
    div_us = divergence_from_edge_scalar(js_us, ops)
    div_gl = divergence_from_edge_scalar(js_gl, ops)
    mismatch = div_us - div_gl

    terminal_mask = (
        np.zeros(ops.n_nodes, dtype=bool)
        if terminal_node_mask is None
        else np.asarray(terminal_node_mask, dtype=bool).reshape(-1).copy()
    )
    if terminal_mask.size != ops.n_nodes:
        raise ValueError(f"terminal_node_mask has length {terminal_mask.size}, expected {ops.n_nodes}.")
    bulk_mask = nonterminal_node_mask(terminal_mask, n_nodes=ops.n_nodes)

    continuation = phase_drive_continuation or PhaseDriveContinuationSolver.from_operators(ops)
    drive_dimensionless, convergence = continuation.solve(
        psi_dimensionless=psi,
        mismatch_divergence_A_m3=mismatch,
        correction_C_J2_m3_A=coeff.correction_C_J2_m3_A,
        delta0_J=float(material.delta0_J),
    )
    drive_complex = drive_dimensionless * float(material.delta0_J)
    drive_abs = np.abs(drive_complex)

    return AllmarasDiagnosticFields(
        coefficients=coeff,
        edge_Q_m_inv=Q,
        edge_R_J=R_edge_J,
        edge_Te_K=Te_edge,
        edge_js_us_allmaras_A_m2=np.asarray(js_us, dtype=float),
        edge_js_gl_allmaras_A_m2=np.asarray(js_gl, dtype=float),
        node_js_us_allmaras_x_A_m2=np.asarray(js_us_x, dtype=float),
        node_js_us_allmaras_y_A_m2=np.asarray(js_us_y, dtype=float),
        node_js_gl_allmaras_x_A_m2=np.asarray(js_gl_x, dtype=float),
        node_js_gl_allmaras_y_A_m2=np.asarray(js_gl_y, dtype=float),
        node_div_js_us_allmaras_A_m3=np.asarray(div_us, dtype=float),
        node_div_js_gl_allmaras_A_m3=np.asarray(div_gl, dtype=float),
        node_mismatch_divergence_A_m3=np.asarray(mismatch, dtype=float),
        node_phase_drive_complex_J=np.asarray(drive_complex, dtype=np.complex128),
        node_phase_drive_abs_J=np.asarray(drive_abs, dtype=float),
        node_phase_drive_abs_over_delta0=np.asarray(np.abs(drive_dimensionless), dtype=float),
        bulk_node_mask=bulk_mask,
        phase_drive_convergence=convergence,
    )



def analytic_allmaras_usadel_current_edges(
    *,
    edge_Q_m_inv: np.ndarray,
    edge_R_J: np.ndarray,
    edge_Te_K: np.ndarray,
    material: GTDGLMaterial,
) -> np.ndarray:
    """Analytic Allmaras current closure on edges.

    This is kept as a fallback diagnostic.  The production path for the present
    thesis uses the PRE Matsubara/Usadel table instead.
    """

    Q = np.asarray(edge_Q_m_inv, dtype=float)
    R = np.asarray(edge_R_J, dtype=float)
    Te = np.maximum(np.asarray(edge_Te_K, dtype=float), 1.0e-12)
    return (
        np.pi
        * float(material.sigma_n_S_m)
        / (2.0 * E_ABS_C)
        * R
        * np.tanh(R / (2.0 * K_B_J_K * Te))
        * Q
    )


def gl_supercurrent_edges(
    *,
    edge_Q_m_inv: np.ndarray,
    edge_R_J: np.ndarray,
    material: GTDGLMaterial,
) -> np.ndarray:
    """Auxiliary GL supercurrent used only in the Allmaras mismatch."""

    Q = np.asarray(edge_Q_m_inv, dtype=float)
    R = np.asarray(edge_R_J, dtype=float)
    return (
        np.pi
        * float(material.sigma_n_S_m)
        * R
        * R
        / (4.0 * E_ABS_C * K_B_J_K * float(material.Tc_K))
        * Q
    )


@dataclass(frozen=True)
class PhaseDriveContinuationSolver:
    """Continue the Allmaras phase quotient through low-amplitude regions.

    The exact dimensionless correction contains

        (div j_s^Us - div j_s^GL) * psi / |psi|^2.

    That quotient is evaluated directly where ``|psi|`` is resolved.  Where the
    condensate amplitude is small, its finite continuous limit is obtained as a
    discrete harmonic continuation on the same finite-volume graph.  A
    Jacobi-preconditioned conjugate-gradient solve acts only on the unresolved
    nodes and reports its actual residual and iteration count.
    """

    graph_laplacian: sp.csr_matrix
    direct_amplitude_fraction: float = 1.0e-2
    tolerance: float = 1.0e-3
    max_iterations: int = 64
    zero_amplitude_fraction: float = 1.0e-12

    @classmethod
    def from_operators(
        cls,
        ops: FVOperators,
        *,
        direct_amplitude_fraction: float = 1.0e-2,
        tolerance: float = 1.0e-3,
        max_iterations: int = 64,
    ) -> "PhaseDriveContinuationSolver":
        direct = float(direct_amplitude_fraction)
        tol = float(tolerance)
        iterations = int(max_iterations)
        if not np.isfinite(direct) or not (0.0 < direct < 1.0):
            raise ValueError("direct_amplitude_fraction must lie in (0, 1).")
        if not np.isfinite(tol) or tol <= 0.0:
            raise ValueError("phase-drive convergence tolerance must be positive and finite.")
        if iterations < 1:
            raise ValueError("phase-drive max_iterations must be at least one.")

        edge_i = np.asarray(ops.edge_i, dtype=np.int64)
        edge_j = np.asarray(ops.edge_j, dtype=np.int64)
        edge_length = np.maximum(np.asarray(ops.edge_length_m, dtype=float), 1.0e-300)
        dual_length = np.asarray(ops.dual_face_length_m, dtype=float)
        weights = dual_length / edge_length
        if np.any(~np.isfinite(weights)) or np.any(weights <= 0.0):
            raise ValueError("Finite-volume continuation weights must be positive and finite.")

        rows = np.concatenate((edge_i, edge_j))
        cols = np.concatenate((edge_j, edge_i))
        data = np.concatenate((weights, weights))
        degree = np.bincount(rows, weights=data, minlength=ops.n_nodes)
        if np.any(degree <= 0.0):
            raise ValueError("Every node must have at least one positive continuation weight.")
        adjacency = sp.csr_matrix(
            (data, (rows, cols)),
            shape=(ops.n_nodes, ops.n_nodes),
            dtype=float,
        )
        graph_laplacian = (sp.diags(degree, format="csr") - adjacency).tocsr()
        return cls(
            graph_laplacian=graph_laplacian,
            direct_amplitude_fraction=direct,
            tolerance=tol,
            max_iterations=iterations,
        )

    def solve(
        self,
        *,
        psi_dimensionless: np.ndarray,
        mismatch_divergence_A_m3: np.ndarray,
        correction_C_J2_m3_A: np.ndarray,
        delta0_J: float,
    ) -> tuple[np.ndarray, PhaseDriveConvergenceInfo]:
        psi = np.asarray(psi_dimensionless, dtype=np.complex128).reshape(-1)
        mismatch = np.asarray(mismatch_divergence_A_m3, dtype=float).reshape(-1)
        correction = np.asarray(correction_C_J2_m3_A, dtype=float)
        if psi.size != self.graph_laplacian.shape[0] or mismatch.shape != psi.shape:
            raise ValueError("Phase-drive continuation arrays must match the FV node count.")
        if correction.ndim == 0:
            correction = np.full(psi.shape, float(correction), dtype=float)
        else:
            correction = correction.reshape(-1)
        if correction.shape != psi.shape:
            raise ValueError("correction_C_J2_m3_A must be scalar or node shaped.")
        if np.any(~np.isfinite(psi)) or np.any(~np.isfinite(mismatch)) or np.any(~np.isfinite(correction)):
            raise ValueError("Phase-drive continuation inputs must be finite.")

        delta0 = float(delta0_J)
        if not np.isfinite(delta0) or delta0 <= 0.0:
            raise ValueError("delta0_J must be positive and finite.")

        amplitude = np.abs(psi)
        direct_mask = amplitude >= float(self.direct_amplitude_fraction)
        zero_mask = amplitude <= float(self.zero_amplitude_fraction)
        continuation_mask = ~direct_mask

        quotient = np.zeros(psi.shape, dtype=np.complex128)
        if np.any(direct_mask):
            quotient[direct_mask] = (
                mismatch[direct_mask]
                * psi[direct_mask]
                / np.maximum(amplitude[direct_mask] ** 2, 1.0e-300)
            )

        # A bounded Tikhonov estimate gives CG a local, finite initial guess.
        # The converged harmonic extension is independent of this initialization.
        if np.any(continuation_mask):
            threshold2 = float(self.direct_amplitude_fraction) ** 2
            quotient[continuation_mask] = (
                mismatch[continuation_mask]
                * psi[continuation_mask]
                / (amplitude[continuation_mask] ** 2 + threshold2)
            )

        converged = True
        residual_rel = 0.0
        iterations = 0
        if np.any(continuation_mask):
            if not np.any(direct_mask):
                # With no resolved condensate anywhere, phase has no Dirichlet
                # reference and the finite correction is set to its null limit.
                quotient[continuation_mask] = 0.0
            else:
                continued = np.flatnonzero(continuation_mask)
                direct = np.flatnonzero(direct_mask)
                operator = self.graph_laplacian[continued][:, continued].tocsr()
                rhs = -self.graph_laplacian[continued][:, direct] @ quotient[direct]
                diagonal = np.asarray(operator.diagonal(), dtype=float)
                if np.any(~np.isfinite(diagonal)) or np.any(diagonal <= 0.0):
                    raise RuntimeError("Harmonic continuation operator has an invalid diagonal.")
                preconditioner = spla.LinearOperator(
                    operator.shape,
                    matvec=lambda value: np.asarray(value, dtype=float) / diagonal,
                    dtype=float,
                )

                def solve_component(rhs_component: np.ndarray, initial: np.ndarray):
                    rhs_component = np.asarray(rhs_component, dtype=float)
                    rhs_norm = float(np.linalg.norm(rhs_component))
                    if rhs_norm <= 1.0e-300:
                        return np.zeros(rhs_component.shape, dtype=float), True, 0, 0.0
                    counter = [0]

                    def count_iteration(_value) -> None:
                        counter[0] += 1

                    solution, status = spla.cg(
                        operator,
                        rhs_component,
                        x0=np.asarray(initial, dtype=float),
                        rtol=float(self.tolerance),
                        atol=0.0,
                        maxiter=int(self.max_iterations),
                        M=preconditioner,
                        callback=count_iteration,
                    )
                    residual = float(
                        np.linalg.norm(operator @ solution - rhs_component)
                        / max(rhs_norm, 1.0e-300)
                    )
                    component_converged = bool(
                        status == 0
                        and np.isfinite(residual)
                        and residual <= float(self.tolerance)
                    )
                    return np.asarray(solution, dtype=float), component_converged, counter[0], residual

                real, real_ok, real_iterations, real_residual = solve_component(
                    np.real(rhs),
                    np.real(quotient[continued]),
                )
                imag, imag_ok, imag_iterations, imag_residual = solve_component(
                    np.imag(rhs),
                    np.imag(quotient[continued]),
                )
                quotient[continued] = real + 1j * imag
                converged = bool(real_ok and imag_ok)
                iterations = max(int(real_iterations), int(imag_iterations))
                residual_rel = max(float(real_residual), float(imag_residual))

        if not converged:
            raise RuntimeError(
                "Allmaras phase-drive harmonic continuation did not converge: "
                f"residual={residual_rel:.6e}, tolerance={float(self.tolerance):.6e}, "
                f"iterations={iterations}/{int(self.max_iterations)}."
            )

        phase_drive = 1j * correction * quotient / (delta0 * delta0)
        if np.any(~np.isfinite(phase_drive)):
            raise RuntimeError("Controlled Allmaras phase-drive continuation produced non-finite values.")

        info = PhaseDriveConvergenceInfo(
            converged=bool(converged),
            iterations=int(iterations),
            residual_rel=float(residual_rel),
            direct_node_count=int(np.count_nonzero(direct_mask)),
            continued_node_count=int(np.count_nonzero(continuation_mask)),
            zero_amplitude_node_count=int(np.count_nonzero(zero_mask)),
            direct_amplitude_fraction=float(self.direct_amplitude_fraction),
            tolerance=float(self.tolerance),
            max_iterations=int(self.max_iterations),
        )
        return np.asarray(phase_drive, dtype=np.complex128), info


def compute_allmaras_forcing_dimensionless(
    *,
    psi_dimensionless: np.ndarray,
    psi_laplacian_dimensionless: np.ndarray,
    material: GTDGLMaterial,
    Te_K: np.ndarray,
    ops: FVOperators,
    length_scale_m: float,
    edge_js_usadel_A_m2: np.ndarray | None = None,
    blocked_edge_mask: np.ndarray | None = None,
    phase_drive_continuation: PhaseDriveContinuationSolver | None = None,
    r_epsilon_fraction: float = 1.0e-6,
) -> AllmarasForcingFields:
    """Return the Appendix-B forcing injected into the local KWT update.

    The continuous equation used here is

        tau0 [D_t Delta + alpha_KWT d_t |Delta|^2 Delta]
        = rho_KWT F[Delta],

    with

        F/Delta0 = (xi_mod^2/L0^2) L' psi
                 + [1 - Te/Tc - |Delta|^2/Delta_mod^2] psi
                 + i (C/Delta0^2) H[(div j_s^Us-div j_s^GL) psi/|psi|^2].

    ``L'`` is the promoted pyTDGL-like Laplacian in dimensionless coordinates
    x' = x/L0 and ``H`` denotes the controlled harmonic continuation through
    low-amplitude nodes.  The returned array is dimensionless; the solver
    multiplies it by ``rho_KWT * dt`` because its time variable is already
    scaled by ``tau0``.
    """

    psi = np.asarray(psi_dimensionless, dtype=np.complex128).reshape(-1)
    lap = np.asarray(psi_laplacian_dimensionless, dtype=np.complex128).reshape(-1)
    Te = np.asarray(Te_K, dtype=float).reshape(-1)
    if psi.shape != (ops.n_nodes,):
        raise ValueError(f"psi_dimensionless must have shape ({ops.n_nodes},), got {psi.shape}.")
    if lap.shape != (ops.n_nodes,):
        raise ValueError(f"psi_laplacian_dimensionless must have shape ({ops.n_nodes},), got {lap.shape}.")
    if Te.shape != (ops.n_nodes,):
        raise ValueError(f"Te_K must have shape ({ops.n_nodes},), got {Te.shape}.")

    L0 = float(length_scale_m)
    if not np.isfinite(L0) or L0 <= 0.0:
        raise ValueError("length_scale_m must be positive and finite.")

    coeff = allmaras_coefficients(psi_dimensionless=psi, material=material, Te_K=Te)
    Q = edge_phase_gradient_from_psi(psi, ops)
    Delta_J = psi * float(material.delta0_J)
    R_node_J = np.abs(Delta_J)
    R_edge_J = edge_average(R_node_J, ops)
    Te_edge = np.maximum(edge_average(Te, ops), 1.0e-12)

    if edge_js_usadel_A_m2 is None:
        js_us = analytic_allmaras_usadel_current_edges(
            edge_Q_m_inv=Q,
            edge_R_J=R_edge_J,
            edge_Te_K=Te_edge,
            material=material,
        )
    else:
        js_us = np.asarray(edge_js_usadel_A_m2, dtype=float).reshape(-1)
        if js_us.shape != (ops.n_edges,):
            raise ValueError(f"edge_js_usadel_A_m2 must have shape ({ops.n_edges},), got {js_us.shape}.")
    js_gl = gl_supercurrent_edges(edge_Q_m_inv=Q, edge_R_J=R_edge_J, material=material)

    if blocked_edge_mask is not None:
        mask = np.asarray(blocked_edge_mask, dtype=bool).reshape(-1)
        if mask.size != ops.n_edges:
            raise ValueError(f"blocked_edge_mask has length {mask.size}, expected {ops.n_edges}.")
        if np.any(mask):
            js_us = js_us.copy()
            js_gl = js_gl.copy()
            js_us[mask] = 0.0
            js_gl[mask] = 0.0

    div_us = divergence_from_edge_scalar(js_us, ops)
    div_gl = divergence_from_edge_scalar(js_gl, ops)
    mismatch = div_us - div_gl

    diffusion = (coeff.xi_mod2_m2 / (L0 * L0)) * lap

    delta_mod2 = np.asarray(coeff.delta_mod2_J2, dtype=float)
    delta0 = float(material.delta0_J)
    delta_mod_floor2 = (float(r_epsilon_fraction) * delta0) ** 2
    denom = np.maximum(delta_mod2, delta_mod_floor2)
    reaction = (1.0 - Te / float(material.Tc_K) - (np.abs(Delta_J) ** 2) / denom) * psi

    continuation = phase_drive_continuation or PhaseDriveContinuationSolver.from_operators(ops)
    phase_drive, phase_drive_convergence = continuation.solve(
        psi_dimensionless=psi,
        mismatch_divergence_A_m3=mismatch,
        correction_C_J2_m3_A=coeff.correction_C_J2_m3_A,
        delta0_J=delta0,
    )

    forcing = diffusion + reaction + phase_drive
    forcing = np.asarray(forcing, dtype=np.complex128)
    forcing[~np.isfinite(forcing)] = 0.0

    return AllmarasForcingFields(
        coefficients=coeff,
        forcing_dimensionless=forcing,
        diffusion_dimensionless=np.asarray(diffusion, dtype=np.complex128),
        reaction_dimensionless=np.asarray(reaction, dtype=np.complex128),
        phase_drive_dimensionless=np.asarray(phase_drive, dtype=np.complex128),
        edge_Q_m_inv=Q,
        edge_R_J=R_edge_J,
        edge_Te_K=Te_edge,
        edge_js_us_A_m2=np.asarray(js_us, dtype=float),
        edge_js_gl_A_m2=np.asarray(js_gl, dtype=float),
        node_div_js_us_A_m3=np.asarray(div_us, dtype=float),
        node_div_js_gl_A_m3=np.asarray(div_gl, dtype=float),
        node_mismatch_divergence_A_m3=np.asarray(mismatch, dtype=float),
        phase_drive_convergence=phase_drive_convergence,
    )


def nonterminal_node_mask(terminal_node_mask: np.ndarray, *, n_nodes: int) -> np.ndarray:
    """Return the diagnostic domain with only exact terminal nodes excluded."""

    terminal = np.asarray(terminal_node_mask, dtype=bool).reshape(-1)
    if terminal.size != int(n_nodes):
        raise ValueError(f"terminal_node_mask has length {terminal.size}, expected {n_nodes}.")
    bulk = ~terminal
    return bulk if np.any(bulk) else np.ones(int(n_nodes), dtype=bool)


def rms(values: np.ndarray, mask: np.ndarray | None = None) -> float:
    """Finite-safe RMS helper for summaries."""

    arr = np.asarray(values, dtype=float)
    if mask is not None:
        m = np.asarray(mask, dtype=bool)
        arr = arr[m]
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return float("nan")
    return float(np.sqrt(np.mean(arr * arr)))


def max_abs(values: np.ndarray, mask: np.ndarray | None = None) -> float:
    """Finite-safe max-absolute helper for summaries."""

    arr = np.asarray(values, dtype=float)
    if mask is not None:
        m = np.asarray(mask, dtype=bool)
        arr = arr[m]
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return float("nan")
    return float(np.max(np.abs(arr)))
