"""Reusable first-gradient spectral response with the full 2x2 algebra.

The rectangular spectral/normalization operator is assembled once at a fixed
R0. A sparse normal-equation factorization is a performance option only:
there is no diagonal shift, damping, eta bath, row deletion beyond prescribed
cuts, or rank completion. Its answers are checked against the ORIGINAL
rectangular system, and LSMR on that original system is the explicit fallback.
Every answer also re-evaluates the independent full nonlinear Jet residual.
Solver status and small residuals do not constitute physical admission.
"""
from __future__ import annotations

from dataclasses import replace
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import lsmr, splu
from . import adiabatic_matrix_usadel as ad


def _positive(value, name):
    raw = np.asarray(value)
    if (raw.ndim or np.iscomplexobj(raw) or raw.dtype.kind == 'b'
            or not np.isfinite(raw) or raw <= 0):
        raise ValueError(name+' must be a positive finite real scalar')
    return float(raw)


def _nodes(value, n):
    raw = np.asarray(value)
    if raw.size == 0:
        return np.empty(0, int)
    if (raw.ndim != 1 or raw.dtype.kind not in 'iu' or np.any(raw < 0)
            or np.any(raw >= n) or len(np.unique(raw)) != len(raw)):
        raise ValueError('Unique integer prescribed nodes required')
    return raw.astype(int, copy=True)


def _assemble(graph, R, B, alpha):
    n = graph.n_nodes
    basis = np.eye(4, dtype=complex).reshape(4, 2, 2)
    rows, columns, values = [], [], []
    def block(i, j, entries, offset=0):
        entries = np.asarray(entries).reshape(4, 4).T
        rr, cc = np.indices((4, 4))
        rows.extend((offset+4*i+rr).ravel())
        columns.extend((4*j+cc).ravel())
        values.extend(entries.ravel())
    for i in range(n):
        block(i, i, -.5*(B[i]@basis-basis@B[i]))
        block(i, i, R[i]@basis+basis@R[i], offset=4*n)
    links = ad.links(alpha).value
    for edge, (i, j) in enumerate(graph.edges):
        P = links[edge]
        Pd = P.conj().T
        Rj, Ri = P@R[j]@Pd, Pd@R[i]@P
        transported, backwards = P@basis@Pd, Pd@basis@P
        c = graph.conductance[edge]/2
        block(i, i, c/graph.area_weights[i]*(basis@Rj-Rj@basis))
        block(i, j, c/graph.area_weights[i]*(R[i]@transported-transported@R[i]))
        block(j, j, c/graph.area_weights[j]*(basis@Ri-Ri@basis))
        block(j, i, c/graph.area_weights[j]*(R[j]@backwards-backwards@R[j]))
    result = coo_matrix((np.asarray(values), (np.asarray(rows), np.asarray(columns))),
        shape=(8*n, 4*n)).tocsr()
    result.eliminate_zeros()
    return result


class FirstSpectralResponse:
    """Repeated R1 solves for the same graph, R0, link phases and gap state.

    Time derivatives are per second and kappa_s=hbar/(2*Eu) is in seconds.
    R0 and B0 are Jets; only their zeroth values enter the cached operator.
    Later energy/time/first derivatives enter each right-hand side. A scalar
    identity change of B0 leaves its commutator invariant and is permitted.
    All four complex matrix entries are independent unknowns, including the
    identity part required by star normalization under a temporal gauge.
    """
    def __init__(self, graph, R0, B0, alpha=None, fixed_nodes=(), *, kappa_s,
                 factor_normal_equations=True, minimum_pivot_ratio=1e-12):
        self.graph = graph
        self.kappa_s = _positive(kappa_s, 'kappa_s')
        self.minimum_pivot_ratio = _positive(minimum_pivot_ratio, 'minimum_pivot_ratio')
        n = graph.n_nodes
        if R0.value.shape != (n, 2, 2) or B0.value.shape != (n, 2, 2):
            raise ValueError('One full R0 and B0 matrix per node required')
        self.R0_value = np.array(R0.value, complex, copy=True)
        self.B0_value = np.array(B0.value, complex, copy=True)
        self.alpha = np.zeros(len(graph.edges)) if alpha is None else np.asarray(alpha)
        if (np.iscomplexobj(self.alpha) or self.alpha.shape != (len(graph.edges),)
                or np.any(~np.isfinite(self.alpha))):
            raise ValueError('One finite real fixed link phase per edge required')
        self.alpha = np.array(self.alpha, float, copy=True)
        self.fixed_nodes = _nodes(fixed_nodes, n)
        self.free = np.ones(n, bool)
        self.free[self.fixed_nodes] = False
        self.columns = np.repeat(self.free, 4)
        self.rows = np.tile(self.columns, 2)
        self.full_matrix = _assemble(graph, self.R0_value, self.B0_value, self.alpha)
        self.system = self.full_matrix[self.rows][:, self.columns].tocsr()
        self.adjoint = self.system.conj().T.tocsr()
        self.factor = None
        self.factorization_reason = 'normal equations disabled explicitly'
        self.normal_pivot_ratio = None
        self.apply_calls = 0
        self.fallback_calls = 0
        if factor_normal_equations and self.system.shape[1]:
            try:
                factor = splu((self.adjoint@self.system).tocsc())
                pivots = abs(factor.U.diagonal())
                ratio = float(np.min(pivots)/np.max(pivots)) if len(pivots) else 0.
                self.normal_pivot_ratio = ratio
                if not np.isfinite(ratio) or ratio < self.minimum_pivot_ratio:
                    self.factorization_reason = 'small LU pivot; use original rectangular LSMR'
                else:
                    self.factor = factor
                    self.factorization_reason = 'unshifted sparse factor; every solve checks original system'
            except RuntimeError:
                self.factorization_reason = 'singular unshifted normal matrix; use original rectangular LSMR'
        elif not self.system.shape[1]:
            self.factorization_reason = 'all nodes prescribed; no unknown spectral correction'

    @property
    def estimated_storage_bytes(self):
        """Conservative array ledger for cache budgeting; excludes graph/caller Jets."""
        arrays = [self.R0_value, self.B0_value, self.alpha, self.fixed_nodes,
                  self.free, self.columns, self.rows]
        matrices = [self.full_matrix, self.system, self.adjoint]
        if self.factor is not None:
            matrices.extend((self.factor.L, self.factor.U))
        return int(sum(a.nbytes for a in arrays)+sum(m.data.nbytes+m.indices.nbytes
            +m.indptr.nbytes for m in matrices))

    def apply(self, retarded, B, *, alpha_time_per_s=None, fixed_first=None,
              tolerance=1e-10, maxiter=None):
        """Return (full corrected Jet, independently evaluated residual metrics).

        Explicit fixed_first is required when cuts are prescribed; it must be
        the correction of the same evolving superconducting continuation.
        A failed least-squares residual is returned visibly, not admitted.
        """
        tol = _positive(tolerance, 'tolerance')
        n = self.graph.n_nodes
        if retarded.value.shape != self.R0_value.shape or not np.array_equal(retarded.value, self.R0_value):
            raise ValueError('Cached response requires the identical instantaneous R0')
        if B.value.shape != self.B0_value.shape:
            raise ValueError('B must have the cached matrix shape')
        difference = B.value-self.B0_value
        traceless = difference-np.trace(difference, axis1=-2, axis2=-1)[..., None, None]*ad.I2/2
        if np.max(abs(traceless), initial=0.) > 32*np.finfo(float).eps*max(1., float(np.max(abs(B.value)))):
            raise ValueError('B commutator changed; rebuild the spectral response')
        if maxiter is not None and (type(maxiter) is not int or maxiter < 1):
            raise ValueError('Positive integer LSMR iteration limit required')
        zero = replace(retarded, first=np.zeros_like(retarded.value))
        base, norm = ad.spectral_residual(self.graph, zero, B, alpha=self.alpha,
            alpha_time=alpha_time_per_s, kappa=self.kappa_s)
        rhs = -np.concatenate(((base.first/self.graph.area_weights[:, None, None]).ravel(),
            norm.first.ravel()))
        known = np.zeros((n, 2, 2), complex)
        if len(self.fixed_nodes):
            if fixed_first is None:
                raise ValueError('Explicit first correction of the evolving cuts required')
            value = np.asarray(fixed_first, complex)
            if value.shape != (len(self.fixed_nodes), 2, 2) or np.any(~np.isfinite(value)):
                raise ValueError('One finite full correction per prescribed cut node required')
            known[self.fixed_nodes] = value
        elif fixed_first is not None and np.asarray(fixed_first).size:
            raise ValueError('Cut correction supplied without prescribed nodes')
        rhs = (rhs-self.full_matrix@known.ravel())[self.rows]
        self.apply_calls += 1
        solver, lsmr_stop, iterations = 'prescribed', None, 0
        answer = np.empty(0, complex)
        if self.system.shape[1]:
            answer = None
            if self.factor is not None:
                candidate = self.factor.solve(self.adjoint@rhs)
                mismatch = self.system@candidate-rhs
                bound = 10*tol*max(1., float(np.max(abs(rhs), initial=0.)))
                if np.all(np.isfinite(candidate)) and np.max(abs(mismatch), initial=0.) <= bound:
                    answer = candidate
                    solver = 'unshifted_normal_factor_checked_on_original_system'
            if answer is None:
                self.fallback_calls += 1
                result = lsmr(self.system, rhs, atol=tol, btol=tol, maxiter=maxiter)
                answer, lsmr_stop, iterations = result[0], int(result[1]), int(result[2])
                solver = 'lsmr_original_rectangular_system'
        solved = known.ravel()
        solved[self.columns] = answer
        solved = solved.reshape(n, 2, 2)
        if np.any(~np.isfinite(solved)):
            raise FloatingPointError('Nonfinite spectral response; no repair applied')
        result = replace(retarded, first=solved)
        residual, normalization = ad.spectral_residual(self.graph, result, B,
            alpha=self.alpha, alpha_time=alpha_time_per_s, kappa=self.kappa_s)
        linear_residual = self.system@answer-rhs
        metrics = dict(solver=solver, lsmr_stop=lsmr_stop, iterations=iterations,
            linear_original_max=float(np.max(abs(linear_residual), initial=0.)),
            spectral_first_max=float(np.max(abs(residual.first[self.free]
                /self.graph.area_weights[self.free, None, None]), initial=0.)),
            normalization_first_max=float(np.max(abs(normalization.first[self.free]), initial=0.)),
            normalization_first_all_max=float(np.max(abs(normalization.first), initial=0.)),
            spectral_zero_max=float(np.max(abs(residual.value[self.free]
                /self.graph.area_weights[self.free, None, None]), initial=0.)),
            normalization_zero_max=float(np.max(abs(normalization.value), initial=0.)),
            identity_first_max=float(np.max(abs(np.trace(solved, axis1=-2, axis2=-1)/2), initial=0.)),
            normal_pivot_ratio=self.normal_pivot_ratio,
            factorization_reason=self.factorization_reason,
            apply_calls=self.apply_calls, fallback_calls=self.fallback_calls,
            admission=False)
        return result, metrics


class ChartFirstSpectralResponse:
    """Full R1 through its exact normalization part and two chart tangents.

    R1=P+DR0[da,db], P=R0*(-i*kappa*[RE,Rt])/2. For normalized traceless
    R0 and its consistent tangents P is a scalar identity field, whose
    zeroth spectral commutator is zero even if the scalar varies by node.
    The two remaining components use the original retarded Riccati Jacobian:

      Fa=-(1+a*b)*(L01-a*L00), Fb=(1+a*b)*(L10-b*L00).

    This identity is exact at a stationary root. Reuse of a numerically
    stationary root retains its measured root error; full residuals remain
    independent diagnostics and are never replaced by the two chart rows.
    No normal-equation matrix, extra factor or hidden least-squares solve is
    built. Inconsistent full sources/cut normalization stay visibly failed.
    """
    def __init__(self, graph, R0, B0, alpha=None, fixed_nodes=(), *, kappa_s,
                 retarded_tangent=None, spectral_tolerance=1e-7,
                 practical_spectral_tolerance=1e-7):
        from types import SimpleNamespace
        from . import retarded_spatial_usadel as ret
        self.graph = graph
        self.kappa_s = _positive(kappa_s, 'kappa_s')
        self.practical_spectral_tolerance = _positive(practical_spectral_tolerance, 'practical_spectral_tolerance')
        self.R0_value = np.array(R0.value, complex, copy=True)
        self.B0_value = np.array(B0.value, complex, copy=True)
        n = graph.n_nodes
        if self.R0_value.shape != (n, 2, 2) or self.B0_value.shape != (n, 2, 2):
            raise ValueError('One full R0/B0 matrix per graph node required')
        self.alpha = np.zeros(len(graph.edges)) if alpha is None else np.asarray(alpha)
        if (np.iscomplexobj(self.alpha) or self.alpha.shape != (len(graph.edges),)
                or np.any(~np.isfinite(self.alpha))):
            raise ValueError('Finite real fixed link phase per edge required')
        self.alpha = np.array(self.alpha, float, copy=True)
        self.fixed_nodes = _nodes(fixed_nodes, n)
        denominator = 1+self.R0_value[:, 0, 0]
        if np.any(abs(denominator) < 1e-14):
            raise FloatingPointError('Chart pole; no clipping or rank repair')
        self.a, self.b = self.R0_value[:, 0, 1]/denominator, self.R0_value[:, 1, 0]/denominator
        self.g = self.R0_value[:, 0, 0]
        # This class is a response at a supplied physical BCS/Usadel state;
        # the scalar identity potential does not enter the stationary root.
        if retarded_tangent is None:
            diagonal = (self.B0_value[:, 0, 0]-self.B0_value[:, 1, 1])/2
            if not np.allclose(diagonal, diagonal[0], rtol=1e-13, atol=1e-14):
                raise ValueError('One fixed causal energy per slice is required')
            gap = self.B0_value[:, 0, 1]
            if not np.allclose(self.B0_value[:, 1, 0], gap.conj(), rtol=1e-13, atol=1e-14):
                raise ValueError('Physical stationary gap generator required')
            solution = SimpleNamespace(z=complex(diagonal[0]), a=self.a, b=self.b)
            retarded_tangent = ret.RetardedTangent(graph, gap, solution, self.alpha,
                fixed_nodes=self.fixed_nodes, spectral_tolerance=spectral_tolerance)
            self.reuses_tangent = False
        else:
            self.reuses_tangent = True
        self.tangent = retarded_tangent
        if (self.tangent.graph.signature != graph.signature
                or not np.array_equal(self.tangent.fixed_nodes, self.fixed_nodes)
                or not np.array_equal(self.tangent.alpha, self.alpha)
                or not np.allclose(self.tangent.matrix, self.R0_value, rtol=1e-12, atol=1e-12)):
            raise ValueError('Reused retarded tangent does not match the same graph, R0, links and cuts')
        base_traceless = self.B0_value-np.trace(self.B0_value, axis1=-2, axis2=-1)[..., None, None]*ad.I2/2
        expected = ad.generator(self.tangent.d, -self.tangent.z.imag, eta=self.tangent.z.real).value
        if not np.allclose(base_traceless, expected, rtol=1e-12, atol=1e-12):
            raise ValueError('Reused retarded tangent has a different stationary generator')
        self.free = self.tangent.free
        self.components = self.tangent.components
        self.jacobian = self.tangent.jacobian
        self.factor = self.tangent.factor
        self.apply_calls = 0
        self.fallback_calls = 0

    @property
    def estimated_storage_bytes(self):
        """Owned array ledger; caller accounts a shared tangent only once."""
        arrays = (self.R0_value, self.B0_value, self.alpha, self.fixed_nodes,
                  self.a, self.b, self.g)
        # Shared tangent storage belongs to the caller and is accounted once
        # there. Standalone construction owns and counts its Jacobian/factor.
        matrices = [] if self.reuses_tangent else [self.jacobian]
        if not self.reuses_tangent and self.factor is not None:
            matrices.extend((self.factor.L, self.factor.U))
        return int(sum(value.nbytes for value in arrays)+sum(matrix.data.nbytes
            +matrix.indices.nbytes+matrix.indptr.nbytes for matrix in matrices))

    def apply(self, retarded, B, *, alpha_time_per_s=None, fixed_first=None,
              tolerance=1e-10, maxiter=None):
        from . import retarded_spatial_usadel as ret
        tol = _positive(tolerance, 'tolerance')
        if not np.array_equal(retarded.value, self.R0_value):
            raise ValueError('Cached chart response requires the identical instantaneous R0')
        difference = B.value-self.B0_value
        traceless = difference-np.trace(difference, axis1=-2, axis2=-1)[..., None, None]*ad.I2/2
        if np.max(abs(traceless), initial=0.) > 32*np.finfo(float).eps*max(1., float(np.max(abs(B.value)))):
            raise ValueError('B commutator changed; rebuild the chart response')
        n = self.graph.n_nodes
        particular = ad.normalization_particular(retarded, self.kappa_s)
        base = replace(retarded, first=particular)
        spectral, _ = ad.spectral_residual(self.graph, base, B, alpha=self.alpha,
            alpha_time=alpha_time_per_s, kappa=self.kappa_s)
        L = spectral.first
        mapped = np.column_stack((-(1+self.a*self.b)*(L[:, 0, 1]-self.a*L[:, 0, 0]),
            (1+self.a*self.b)*(L[:, 1, 0]-self.b*L[:, 0, 0])))
        rhs = -mapped.ravel()
        coordinates = np.zeros((n, 2), complex)
        cut_error = 0.
        if len(self.fixed_nodes):
            if fixed_first is None:
                raise ValueError('Explicit first correction of the evolving cuts required')
            fixed = np.asarray(fixed_first, complex)
            if fixed.shape != (len(self.fixed_nodes), 2, 2) or np.any(~np.isfinite(fixed)):
                raise ValueError('One finite full correction per prescribed cut node required')
            delta = fixed-particular[self.fixed_nodes]
            dg = (delta[:, 0, 0]-delta[:, 1, 1])/2
            a, b, g = self.a[self.fixed_nodes], self.b[self.fixed_nodes], self.g[self.fixed_nodes]
            da = (delta[:, 0, 1]-a*dg)/(1+g)
            db = (delta[:, 1, 0]-b*dg)/(1+g)
            coordinates[self.fixed_nodes] = np.column_stack((da, db))
            reconstructed = ret.field_direction(a, b, da, db)[3]
            cut_error = float(np.max(abs(reconstructed-delta), initial=0.))
        elif fixed_first is not None and np.asarray(fixed_first).size:
            raise ValueError('Cut correction supplied without prescribed nodes')
        raw = coordinates.ravel()
        if self.factor is not None:
            raw[self.components] = self.factor.solve((rhs-self.jacobian@raw)[self.components])
        if np.any(~np.isfinite(raw)):
            raise FloatingPointError('Nonfinite chart spectral response; no repair applied')
        result = replace(retarded, first=particular+ret.field_direction(
            self.a, self.b, coordinates[:, 0], coordinates[:, 1])[3])
        # Retain supplied cuts exactly, even when their normalization is wrong.
        # Their failure is then seen by the independent all-node norm below.
        if len(self.fixed_nodes):
            result.first[self.fixed_nodes] = fixed
        residual, normalization = ad.spectral_residual(self.graph, result, B,
            alpha=self.alpha, alpha_time=alpha_time_per_s, kappa=self.kappa_s)
        chart_residual = (self.jacobian@raw-rhs).reshape(n, 2)
        full_error = float(np.max(abs(residual.first[self.free]
            /self.graph.area_weights[self.free, None, None]), initial=0.))
        norm_error = float(np.max(abs(normalization.first[self.free]), initial=0.))
        all_norm = float(np.max(abs(normalization.first), initial=0.))
        source_scale = max(1., float(np.max(abs(spectral.first[self.free]
            /self.graph.area_weights[self.free, None, None]), initial=0.)),
            float(np.max(abs(result.first), initial=0.)))
        zeroth_error = float(np.max(abs(residual.value[self.free]
            /self.graph.area_weights[self.free, None, None]), initial=0.))
        first_scale = max(1., float(np.max(abs(result.first), initial=0.)))
        tangent_product = self.kappa_s*float(np.max(abs(retarded.energy), initial=0.))*float(
            np.max(abs(retarded.time), initial=0.))
        link_product = 0. if alpha_time_per_s is None else self.kappa_s*float(
            np.max(abs(alpha_time_per_s), initial=0.))*float(np.max(abs(retarded.energy), initial=0.))
        # Measured scaling estimate, not a rigorous error bound: nonlinear
        # chart/graph coefficients and temporal Ward cancellations matter.
        root_error_scale = zeroth_error*(first_scale+tangent_product+link_product)
        self.apply_calls += 1
        metrics = dict(solver='exact_stereographic_tangent_factor', lsmr_stop=None, iterations=0,
            linear_original_max=max(full_error, norm_error),
            linear_chart_max=float(np.max(abs(chart_residual[self.free]), initial=0.)),
            linear_chart_per_mass_max=float(np.max(abs(chart_residual[self.free]
                /self.graph.area_weights[self.free, None]), initial=0.)),
            spectral_first_max=full_error, normalization_first_max=norm_error,
            normalization_first_all_max=all_norm,
            spectral_zero_max=float(np.max(abs(residual.value[self.free]
                /self.graph.area_weights[self.free, None, None]), initial=0.)),
            normalization_zero_max=float(np.max(abs(normalization.value), initial=0.)),
            identity_first_max=float(np.max(abs(np.trace(result.first, axis1=-2, axis2=-1)/2), initial=0.)),
            particular_traceless_max=float(np.max(abs(particular-
                np.trace(particular, axis1=-2, axis2=-1)[..., None, None]*ad.I2/2), initial=0.)),
            cut_tangent_reconstruction_max=cut_error,
            full_equations_within_requested_scaled_tolerance=bool(max(full_error, all_norm, cut_error) <= 10*tol*source_scale),
            full_spectral_source_scale=source_scale,
            full_spectral_relative_residual=full_error/source_scale,
            normalization_first_relative_residual=all_norm/first_scale,
            cached_zeroth_spectral_per_mass=zeroth_error,
            root_accuracy_scale_estimate=root_error_scale,
            root_accuracy_scale_is_certified_bound=False,
            practical_spectral_relative_tolerance=self.practical_spectral_tolerance,
            practical_tolerance_reason='retained cached R0 precision; full residual and algebraic normalization independently reported',
            consistent_with_cached_root_accuracy=bool(full_error/source_scale <= self.practical_spectral_tolerance
                and max(all_norm, cut_error)/first_scale <= 10*tol),
            normal_pivot_ratio=None, factorization_reason='same unshifted original Riccati Jacobian/factor; no normal matrix',
            reuses_retarded_tangent=self.reuses_tangent, apply_calls=self.apply_calls,
            fallback_calls=0, admission=False)
        return result, metrics
