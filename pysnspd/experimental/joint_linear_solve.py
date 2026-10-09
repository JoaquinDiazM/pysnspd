"""Diagonal equilibration and rank-revealing QR for a stored joint RHS.

Inputs already carry the declared physical SI row/column units. Numerical
equilibration only changes the coordinates of the linear algebra; the
acceptance residual is ALWAYS recomputed in the original physical row scales.
No physical equation, kinetic field or force is removed. Exactly one checked
gauge column may be eliminated. This solver makes no temporal, energy,
gradient-expansion, or material admission.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import lstsq
from .joint_nambu_rows import checked_gauge_reduction


def _finite_real(value, name):
    raw = np.asarray(value)
    if np.iscomplexobj(raw) or raw.dtype.kind == "b":
        raise ValueError(name+" must be finite real")
    result = np.asarray(raw, dtype=float)
    if np.any(~np.isfinite(result)):
        raise ValueError(name+" must be finite real")
    return result


def equilibrate_matrix(matrix, *, iterations=8):
    """Return B=A/(row[:,None]*column[None,:]) without mutating A.

    Alternating square-root infinity norms are a Ruiz-style diagonal
    equilibration. Zero rows/columns retain unit scales and their deficiency
    remains visible to QR. No floors, regularization or source terms are added.
    """
    A = _finite_real(matrix, "Matrix")
    if A.ndim != 2 or not all(A.shape):
        raise ValueError("A nonempty two-dimensional matrix is required")
    if (isinstance(iterations, (bool, np.bool_))
            or not isinstance(iterations, (int, np.integer))
            or iterations < 0):
        raise ValueError("A nonnegative integer iteration count is required")
    B = np.array(A, copy=True)
    row = np.ones(A.shape[0])
    column = np.ones(A.shape[1])
    for _ in range(iterations):
        maximum = np.max(abs(B), axis=1)
        scale = np.ones_like(maximum)
        np.sqrt(maximum, out=scale, where=maximum > 0)
        B /= scale[:, None]
        row *= scale
        maximum = np.max(abs(B), axis=0)
        scale = np.ones_like(maximum)
        np.sqrt(maximum, out=scale, where=maximum > 0)
        B /= scale[None, :]
        column *= scale
        if (np.any(~np.isfinite(B)) or np.any(~np.isfinite(row))
                or np.any(~np.isfinite(column))
                or np.any(row <= 0) or np.any(column <= 0)):
            raise FloatingPointError("Numerical equilibration overflow; no hidden scale replacement")
    return B, row, column


def solve_equilibrated_affine_candidate(matrix, base, *,
        relative_rank_tolerance=1e-12, residual_tolerance=1e-7,
        gauge_direction=None, gauge_reference_index=None,
        equilibration_iterations=8):
    """Solve A*x+base=0 by pivoted complete orthogonal factorization.

    The reported rank is that of the STORED matrix at a declared QR threshold,
    not an assertion that its producer is physically accurate to that threshold.
    In a deficient case GELSY minimizes the norm of the equilibrated unknowns;
    that candidate is retained for diagnostics and cannot be marked admitted.
    Original row residuals, rather than the balanced least-squares norm, decide
    the algebraic residual check. A large physical rate remains the caller's
    responsibility even when this purely algebraic check passes.
    """
    A = _finite_real(matrix, "Matrix")
    b = _finite_real(base, "Base")
    if A.ndim != 2 or not all(A.shape) or b.shape != (A.shape[0],):
        raise ValueError("A nonempty matrix and one base value per physical row are required")
    if (not np.isfinite(relative_rank_tolerance)
            or not 0 < relative_rank_tolerance < 1
            or not np.isfinite(residual_tolerance) or residual_tolerance <= 0):
        raise ValueError("Positive explicit QR rank and original-row residual limits required")
    reduced = A
    reduction = None
    if gauge_direction is not None:
        reduction = checked_gauge_reduction(A, gauge_direction,
            reference_index=gauge_reference_index, row_scales=np.ones(A.shape[0]),
            column_scales=np.ones(A.shape[1]), tolerance=1e-10)
        reduced = np.asarray(reduction.matrix)
    elif gauge_reference_index is not None:
        raise ValueError("A gauge reference requires the complete checked gauge direction")
    B, row, column = equilibrate_matrix(reduced, iterations=equilibration_iterations)
    z, _, rank, _ = lstsq(np.array(B, order="F"), -b/row,
        cond=relative_rank_tolerance, overwrite_a=True,
        check_finite=False, lapack_driver="gelsy")
    reduced_x = z/column
    x = reduced_x if reduction is None else reduction.expand(reduced_x)
    if np.any(~np.isfinite(x)):
        raise FloatingPointError("Nonfinite QR candidate; no clipping")
    residual = A@x+b
    maximum = float(np.max(abs(residual), initial=0.))
    balanced_maximum = float(np.max(abs(residual/row), initial=0.))
    norm_base = float(np.linalg.norm(b))
    relative_residual = float(np.linalg.norm(residual)/norm_base) if norm_base else None
    denominator = float(np.linalg.norm(A, ord=np.inf)*np.max(abs(x), initial=0.)
        + np.max(abs(b), initial=0.))
    backward_error = maximum/denominator if denominator else 0.
    full_rank = int(rank) == reduced.shape[1]
    original_rows_ok = maximum <= residual_tolerance
    diagnostic = dict(method="diagonal Ruiz-style equilibration + pivoted QR GELSY",
        rank=int(rank), size=A.shape[1], physical_columns=reduced.shape[1],
        physical_rows=A.shape[0], scaled_residual_max=maximum,
        original_scaled_residual_max=maximum,
        balanced_residual_max=balanced_maximum,
        relative_original_residual_l2=relative_residual,
        normwise_backward_error_inf=backward_error,
        relative_rank_tolerance=float(relative_rank_tolerance),
        original_row_residual_tolerance=float(residual_tolerance),
        full_column_rank=bool(full_rank), original_rows_ok=bool(original_rows_ok),
        algebraically_solved_full_rank=bool(full_rank and original_rows_ok),
        temporal_admission=False, energy_admission=False,
        no_regularization=True, no_physical_rows_removed=True,
        gauge_reduced=reduction is not None,
        gauge_ward_max=None if reduction is None else reduction.scaled_ward_residual,
        gauge_reference_index=gauge_reference_index,
        equilibration_iterations=int(equilibration_iterations),
        row_equilibration_min=float(np.min(row)), row_equilibration_max=float(np.max(row)),
        column_equilibration_min=float(np.min(column)), column_equilibration_max=float(np.max(column)),
        maximum_unknown_scaled=float(np.max(abs(x), initial=0.)),
        condition_number=None,
        rank_scope="Stored-matrix QR rank only; producer accuracy and physical rate admissibility remain independent",
        deficient_candidate_norm="equilibrated unknowns, not physical fields",
        acceptance_scope="original declared physical row scales only")
    return x, diagnostic
