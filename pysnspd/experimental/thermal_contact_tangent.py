"""Exact real Cartesian thermal Usadel tangent with moving BCS contacts.

This is separate from the existing fixed-contact SpectralTangent. It keeps the
declared boundary law u_contact=d_contact/epsilon and its nonzero derivative.
Contact reactions are returned, not silently dropped from a reduced energy.
A complex vector encodes two REAL Cartesian gap variations; a Fourier phasor
must be assembled by complex-linear combination of separate physical columns.
"""
from __future__ import annotations
from dataclasses import dataclass
from types import SimpleNamespace
import numpy as np
from scipy.sparse.linalg import splu
from . import thermal_spatial_usadel as thermal


@dataclass(frozen=True)
class MovingContactDirection:
    du: np.ndarray
    df: np.ndarray
    dg: np.ndarray
    current_derivative: np.ndarray
    residual_derivative: np.ndarray
    spectral_gradient_derivative: np.ndarray
    fixed_nodes: np.ndarray
    boundary_residual_derivative: np.ndarray
    boundary_gradient_derivative: np.ndarray
    free_equation_residual_max: float
    fixed_law_residual_max: float

    @property
    def equation_residual(self):
        return self.free_equation_residual_max

    @property
    def boundary_reaction(self):
        return self.boundary_gradient_derivative


class ThermalMovingContactTangent:
    """Reusable exact sparse Jacobian with prescribed moving contact values."""
    def __init__(self, graph, d, solution, alpha=None, *, contact_tolerance=1e-10):
        self.graph = graph
        self.d = thermal._complex_vector(d, graph.n_nodes, 'd').copy()
        self.alpha = thermal._alpha(graph, alpha).copy()
        self.solution = solution
        if (not solution.converged or
                solution.source_signature != thermal._signature(graph, self.d, self.alpha)):
            raise ValueError('A stationary solution of the same gap, graph and links is required')
        if not np.isfinite(contact_tolerance) or contact_tolerance <= 0:
            raise ValueError('Positive finite contact tolerance required')
        self.fixed_nodes = thermal._fixed_indices(solution.fixed_nodes, graph.n_nodes)
        if not len(self.fixed_nodes):
            raise ValueError('Explicit BCS contact nodes required')
        self.free = np.ones(graph.n_nodes, bool)
        self.free[self.fixed_nodes] = False
        self.components = np.flatnonzero(np.repeat(self.free, 2))
        self.fixed_components = np.flatnonzero(np.repeat(~self.free, 2))
        self.epsilon = float(solution.epsilon)
        contact_scale = np.maximum(1., abs(self.d[self.fixed_nodes]/self.epsilon))
        contact_error = abs(solution.u[self.fixed_nodes]-self.d[self.fixed_nodes]/self.epsilon)
        if np.any(contact_error > contact_tolerance*contact_scale):
            raise ValueError('Reference contacts do not obey the declared BCS law u=d/epsilon')
        value, jacobian = thermal.spectral_residual_jacobian(
            graph, self.d, self.epsilon, solution.u, self.alpha)
        self.value, self.jacobian = value, jacobian.tocsr()
        self.f, self.g, self.u = value.f, value.g, solution.u
        self.factor = splu(self.jacobian[self.components][:, self.components].tocsc()) if len(self.components) else None

    def apply(self, gap_direction):
        """Return the full-node derivative, including spectral contact work.

        Free rows solve Jff du_f=m*delta_d-Jfb du_b, with
        du_b=delta_d_b/epsilon. The full residual derivative need not vanish
        at a prescribed reservoir: those are its explicit reactions.
        """
        dd = thermal._complex_vector(gap_direction, self.graph.n_nodes, 'gap_direction')
        source = self.graph.area_weights*dd
        real_source = np.column_stack((source.real, source.imag)).ravel()
        derivative = np.zeros(2*self.graph.n_nodes)
        fixed_du = dd[self.fixed_nodes]/self.epsilon
        derivative[self.fixed_components] = np.column_stack((fixed_du.real, fixed_du.imag)).ravel()
        if self.factor is not None:
            rhs = real_source[self.components]-self.jacobian[self.components][:, self.fixed_components]@derivative[self.fixed_components]
            derivative[self.components] = self.factor.solve(rhs)
        du = derivative[0::2]+1j*derivative[1::2]
        projection = np.real(np.conj(self.u)*du)
        dg = -self.g**3*projection
        df = self.g*du+self.u*dg
        tail, head = self.graph.edges.T
        transport = np.exp(-1j*self.alpha)
        current = 2*self.graph.conductance*np.imag(
            np.conj(df[tail])*transport*self.f[head]+
            np.conj(self.f[tail])*transport*df[head])
        residual_real = self.jacobian@derivative-real_source
        dr = residual_real[0::2]+1j*residual_real[1::2]
        # Exact derivative of gradient_u=2g(R-g^2*u*Re(conj(u)R)).
        r = self.value.residual
        p = np.real(np.conj(self.u)*r)
        dp = np.real(np.conj(du)*r+np.conj(self.u)*dr)
        inner = r-self.g**2*self.u*p
        gradient_derivative = 2*dg*inner+2*self.g*(
            dr-2*self.g*dg*self.u*p-self.g**2*du*p-self.g**2*self.u*dp)
        if not all(np.all(np.isfinite(value)) for value in (du, df, dg, current, dr, gradient_derivative)):
            raise FloatingPointError('Nonfinite moving-contact spectral tangent')
        return MovingContactDirection(du, df, dg, current, dr, gradient_derivative,
            self.fixed_nodes.copy(), dr[self.fixed_nodes], gradient_derivative[self.fixed_nodes],
            float(np.max(abs(dr[self.free]), initial=0.)),
            float(np.max(abs(du[self.fixed_nodes]-fixed_du), initial=0.)))


class ThermalContactTangent(ThermalMovingContactTangent):
    """Field-based entry point for resident spectral workers.

    Validate the supplied stationary fields directly; no saved solution object
    or invented convergence receipt is needed from the caller. One sparse
    factorization is retained for all apply(direction) calls at this frequency.
    """
    def __init__(self, graph, d, epsilon, u, alpha=None, *, fixed_nodes=None,
                 spectral_tolerance=1e-7, contact_tolerance=1e-10):
        gap = thermal._complex_vector(d, graph.n_nodes, 'd')
        field = thermal._complex_vector(u, graph.n_nodes, 'u')
        links = thermal._alpha(graph, alpha)
        if fixed_nodes is None:
            fixed_nodes = graph.boundary_nodes
        if fixed_nodes is None:
            raise ValueError('Explicit contact nodes required')
        fixed = thermal._fixed_indices(fixed_nodes, graph.n_nodes)
        free = np.ones(graph.n_nodes, bool); free[fixed] = False
        if not np.isfinite(spectral_tolerance) or spectral_tolerance <= 0:
            raise ValueError('Positive finite spectral tolerance required')
        value = thermal.spectral_energy_gradient(graph, gap, epsilon, field, links)
        residual = float(np.max(abs(value.residual[free])/
            (graph.area_weights[free]*np.maximum(1., abs(gap[free]))), initial=0.))
        if residual > spectral_tolerance:
            raise ValueError('Supplied spectral fields are not stationary on the free nodes')
        verified = SimpleNamespace(converged=True, epsilon=float(epsilon), u=field.copy(),
            f=value.f, g=value.g, fixed_nodes=fixed,
            source_signature=thermal._signature(graph, gap, links))
        self.reference_spectral_residual = residual
        super().__init__(graph, gap, verified, links, contact_tolerance=contact_tolerance)
