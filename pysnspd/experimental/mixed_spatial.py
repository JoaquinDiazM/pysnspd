"""One static GLL energy for a rectangle and two uniform 1D continuations.

The prolongation R identifies every field value on each rectangular end section
with the corresponding 1D endpoint. This admits ONLY the uniform transverse
trace mode, not higher transverse modes or their radiation into the leads.
Restriction of the variational force is R.T: stationary free interface DOFs
balance the integrated sector reactions of this same energy. This is the field
part of D.25 under its transverse-reduction assumption, not its full kinetic
condition. Populations remain sector-quadrature data; p(x) is never silently
identified across different spectra. Equal-energy f(E) remapping, phonon state
ownership and time evolution are outside this module.

Both derivative directions use compatible positive GLL stiffnesses. Gauge
currents are derivatives of the same energy, computed line by line, without
inverting a divergence or assigning an arbitrary circulating current. The
identified interface has no independent transverse gauge edge. Gauge changes
must respect its uniform trace, as any vector of independent DOF phases does.

Dense directional matrices intentionally target small manufactured/static
experiments. The positive electrical face graph is provided as a separate
geometric operator; its Laplacian is not the high-order GLL stiffness.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.sparse import coo_matrix

from .energy_catalog import E_CHARGE_C, HBAR_J_S
from .spatial_functional import PeriodicSpatialFunctional, SpatialAdmissibilityError
from .spatial_open import lobatto_operators


@dataclass(frozen=True)
class MixedSampledFields:
    delta_quadrature_bar: np.ndarray
    derivative_quadrature_bar: np.ndarray
    amplitude_quadrature_bar: np.ndarray
    gamma_quadrature_bar: np.ndarray
    q_delta_quadrature_bar: np.ndarray


@dataclass(frozen=True)
class MixedPrincipalSymbol:
    matrix: np.ndarray
    eigenvalues: np.ndarray
    uncertainty: float
    stable: bool
    amplitude_bar: float
    gamma_bar: float
    spatial_dimensions: int
    fd_steps: tuple
    gamma_curvatures: tuple


@dataclass(frozen=True)
class MixedSpatialEvaluation:
    energy_bar: float
    energy_J: float
    free_energy_bar: float | None
    gradient_cartesian_bar: np.ndarray
    quadrature_gradient_cartesian_bar: np.ndarray
    link_current_bar: np.ndarray
    current_A: np.ndarray
    sampled_fields: MixedSampledFields
    electronic_derivatives: np.ndarray
    p_quadrature: np.ndarray
    principal_symbols: tuple | None
    gradient_remainder_bar: float
    noether_residual_bar: np.ndarray
    global_phase_residual_bar: float
    line_phase_residuals_bar: np.ndarray
    interface_reactions_bar: dict


class MixedSpatialFunctional(PeriodicSpatialFunctional):
    """Rectangle (x in [0,L], y in [-W/2,W/2]) with two same-section leads.

    Energy units are N0*Delta0**2*(width*thickness)*ell0, matching the open
    longitudinal functional. ``mass_bar`` is the assembled DOF measure;
    ``quadrature_mass_bar`` is its unassembled sector measure. Their sums both
    equal total longitudinal length/ell0. Endpoint quadrature contributions
    account for neighbouring subvolumes, not extra interface material.
    """
    scheme = "mixed_GLL_rectangle_uniform_trace_1D_v1"

    def __init__(self, catalog, rectangle_length_m, width_m, thickness_m,
                 left_length_m, right_length_m, *, elements_x, elements_y,
                 left_elements, right_elements, degree=4, Tc_K=None):
        for value in (rectangle_length_m, width_m, thickness_m, left_length_m, right_length_m):
            if not np.isfinite(value) or value <= 0:
                raise ValueError("all geometric lengths must be finite and positive")
        super().__init__(catalog, rectangle_length_m+left_length_m+right_length_m,
                         width_m*thickness_m, 3, Tc_K=Tc_K)
        self.width_m = float(width_m)
        self.thickness_m = float(thickness_m)
        self.rectangle_length_m = float(rectangle_length_m)
        self.left_length_m, self.right_length_m = float(left_length_m), float(right_length_m)
        self.degree = degree
        ell = self.ell0_m
        xl, ml, dl, kl = lobatto_operators(left_elements, left_length_m/ell, degree)
        xr, mr, dr, kr = lobatto_operators(right_elements, right_length_m/ell, degree)
        x, mx, dx, kx = lobatto_operators(elements_x, rectangle_length_m/ell, degree)
        y, my, dy, ky = lobatto_operators(elements_y, width_m/ell, degree)
        xl -= left_length_m/ell
        xr += rectangle_length_m/ell
        y -= width_m/(2*ell)
        nl, nr, nx, ny = len(xl), len(xr), len(x), len(y)
        self.rectangle_shape = (nx, ny)
        self.quadrature_size = nl+nx*ny+nr
        self.sector_slices = dict(left=slice(0, nl), rectangle=slice(nl, nl+nx*ny),
                                  right=slice(nl+nx*ny, self.quadrature_size))
        self.rectangle_quadrature = np.arange(nl, nl+nx*ny).reshape(nx, ny)
        rectangle_map = np.empty((nx, ny), dtype=int)
        rectangle_map[0] = nl-1
        rectangle_map[1:-1] = np.arange(nl, nl+(nx-2)*ny).reshape(nx-2, ny)
        right_interface = nl+(nx-2)*ny
        rectangle_map[-1] = right_interface
        self.quadrature_to_dof = np.r_[np.arange(nl), rectangle_map.ravel(),
                                       np.arange(right_interface, right_interface+nr)]
        self.cells = right_interface+nr
        self.h_bar = None  # GLL/mixed geometry has no single uniform spacing.
        self.interface_dofs = dict(left=nl-1, right=right_interface)
        self.terminal_dofs = (0, self.cells-1)
        self.prolongation = coo_matrix((np.ones(self.quadrature_size),
            (np.arange(self.quadrature_size), self.quadrature_to_dof)),
            shape=(self.quadrature_size, self.cells)).tocsr()
        self.quadrature_coordinates_bar = np.vstack((np.column_stack((xl, np.zeros(nl))),
            np.column_stack((np.repeat(x, ny), np.tile(y, nx))), np.column_stack((xr, np.zeros(nr)))))
        self.dof_coordinates_bar = np.empty((self.cells, 2))
        self.dof_coordinates_bar[self.quadrature_to_dof] = self.quadrature_coordinates_bar
        self.dof_coordinates_bar[nl-1] = [0., 0.]
        self.dof_coordinates_bar[right_interface] = [rectangle_length_m/ell, 0.]
        width_bar = width_m/ell
        self.quadrature_mass_bar = np.r_[ml, np.kron(mx, my/width_bar), mr]
        self.mass_bar = np.asarray(self.prolongation.T@self.quadrature_mass_bar)
        volume_scale = self.cross_section_m2*ell
        self.quadrature_volumes_m3 = volume_scale*self.quadrature_mass_bar
        self.node_volumes_m3 = volume_scale*self.mass_bar
        self.spatial_dimensions = np.r_[np.ones(nl, dtype=int), np.full(nx*ny, 2), np.ones(nr, dtype=int)]
        self._derivatives = np.zeros((2, self.quadrature_size, self.quadrature_size))
        self._stiffnesses = np.zeros_like(self._derivatives)
        sl, sc, sr = (self.sector_slices[name] for name in ("left", "rectangle", "right"))
        self._derivatives[0, sl, sl] = dl
        self._derivatives[0, sc, sc] = np.kron(dx, np.eye(ny))
        self._derivatives[0, sr, sr] = dr
        self._derivatives[1, sc, sc] = np.kron(np.eye(nx), dy)
        self._stiffnesses[0, sl, sl] = kl
        self._stiffnesses[0, sc, sc] = np.kron(kx, np.diag(my/width_bar))
        self._stiffnesses[0, sr, sr] = kr
        self._stiffnesses[1, sc, sc] = np.kron(np.diag(mx/width_bar), ky)
        self._lines = []
        edges, lengths, areas = [], [], []

        def line(indices, axis, area):
            indices = np.asarray(indices, dtype=int)
            mapped = self.quadrature_to_dof[indices]
            if np.all(mapped == mapped[0]):
                edge_ids = np.array([], dtype=int)  # Uniform interface trace.
            else:
                if np.any(mapped[:-1] == mapped[1:]):
                    raise ValueError("partly collapsed line is unsupported")
                edge_ids = np.arange(len(edges), len(edges)+len(indices)-1)
                edges.extend(zip(mapped[:-1], mapped[1:]))
                lengths.extend(np.diff(self.quadrature_coordinates_bar[indices, axis])*ell)
                areas.extend([area]*(len(indices)-1))
            self._lines.append((indices, axis, edge_ids))

        line(np.arange(nl), 0, self.cross_section_m2)
        for j in range(ny):
            line(self.rectangle_quadrature[:, j], 0, thickness_m*ell*my[j])
        line(np.arange(sr.start, sr.stop), 0, self.cross_section_m2)
        for i in range(nx):
            line(self.rectangle_quadrature[i], 1, thickness_m*ell*mx[i])
        self.graph_edges = np.asarray(edges, dtype=int)
        self.edge_lengths_m = np.asarray(lengths)
        self.face_areas_m2 = np.asarray(areas)
        nedge = len(edges)
        self.incidence = coo_matrix((np.r_[np.ones(nedge), -np.ones(nedge)],
            (np.r_[self.graph_edges[:, 0], self.graph_edges[:, 1]], np.tile(np.arange(nedge), 2))),
            shape=(self.cells, nedge)).tocsr()

    def conductance(self, conductivity_S_m):
        """Positive SI graph conductances sigma*S/length; not a SEM Laplacian."""
        raw = np.asarray(conductivity_S_m)
        if np.iscomplexobj(raw):
            raise ValueError("conductivity must be real")
        sigma = np.asarray(raw, dtype=float)
        if sigma.ndim == 0:
            sigma = np.full(len(self.graph_edges), float(sigma))
        if sigma.shape != (len(self.graph_edges),) or np.any(~np.isfinite(sigma)) or np.any(sigma <= 0):
            raise ValueError("one finite positive conductivity per edge, or a scalar, required")
        return sigma*self.face_areas_m2/self.edge_lengths_m

    def _operators(self, delta_bar, link_phases):
        z = np.asarray(delta_bar, dtype=complex)
        if z.shape != (self.cells,) or np.any(~np.isfinite(z)):
            raise ValueError("one finite complex field per independent DOF required")
        raw = np.zeros(len(self.graph_edges)) if link_phases is None else np.asarray(link_phases)
        if np.iscomplexobj(raw) or raw.shape != (len(self.graph_edges),) or np.any(~np.isfinite(raw)):
            raise ValueError("one finite real gauge phase per graph edge required")
        phases = np.asarray(raw, dtype=float)
        derivatives = self._derivatives.astype(complex)
        stiffnesses = self._stiffnesses.astype(complex)
        for indices, axis, edge_ids in self._lines:
            if not len(edge_ids):
                continue
            transport = np.exp(1j*np.r_[0., np.cumsum(phases[edge_ids])])
            conjugation = np.conj(transport)[:, None]*transport[None, :]
            block = np.ix_(indices, indices)
            derivatives[axis][block] *= conjugation
            stiffnesses[axis][block] *= conjugation
        return z, np.asarray(self.prolongation@z), derivatives, stiffnesses

    def _sample(self, zq, derivatives):
        d = np.column_stack([operator@zq for operator in derivatives])
        amplitude = abs(zq)
        q = np.imag(np.conj(zq)[:, None]*d)/(amplitude[:, None]**2+self.delta_regularizer_bar**2)
        gamma = np.sum(q*q, axis=1)/self.gap_ratio
        return MixedSampledFields(zq, d, amplitude, gamma, q)

    def sample_fields(self, delta_bar, *, link_phases=None):
        """Sector samples: rows retain their local spectrum, including interfaces."""
        _, zq, derivatives, _ = self._operators(delta_bar, link_phases)
        return self._sample(zq, derivatives)

    def _mixed_symbol(self, z, q, gamma, p, cache, dimensions, gamma_step=2e-4):
        a, rho = abs(z), abs(z)**2
        s = rho+self.delta_regularizer_bar**2
        fg = self._potential(a, gamma, p, cache)[2]
        q = np.asarray(q[:dimensions])
        steps, estimates, fgg, error = (), (), 0., 0.
        if np.any(q) and rho != 0:
            lo, hi = self.catalog.vacuum.gamma_axis[[0, -1]]
            step = min(gamma_step*max(1., gamma), float((hi-lo)/8))
            if gamma-lo >= step and hi-gamma >= step:
                def stencil(h):
                    v = [self._potential(a, gamma+j*h, p, cache)[2] for j in (-1, 1)]
                    return (v[1]-v[0])/(2*h), sum(abs(t) for t in v)/h
            else:
                sign = 1. if hi-gamma >= gamma-lo else -1.
                step = min(step, float((hi-gamma if sign > 0 else gamma-lo)/3))
                def stencil(h):
                    v = [self._potential(a, gamma+sign*j*h, p, cache)[2] for j in (0, 1, 2)]
                    return sign*(-3*v[0]+4*v[1]-v[2])/(2*h), sum(abs(c*t) for c, t in zip((3, 4, 1), v))/h
            if step <= 0 or gamma+step == gamma:
                raise ValueError("Gamma support cannot resolve the mixed principal symbol")
            steps = (step, step/2, step/4)
            computed = [stencil(h) for h in steps]
            estimates = tuple(float(v[0]) for v in computed)
            fgg = estimates[-1]
            error = 8*max(abs(estimates[1]-estimates[0]), abs(estimates[2]-estimates[1]))
            error += 32*np.finfo(float).eps*max(v[1] for v in computed)
        jz = np.array([-z.imag, z.real])
        spatial = ((2*fg/self.gap_ratio-2*self.kappa*rho)*np.eye(dimensions)
                   +4*fgg/self.gap_ratio**2*np.outer(q, q))
        matrix = 2*self.kappa*np.eye(2*dimensions)+np.kron(spatial, np.outer(jz, jz)/s**2)
        eigen = np.linalg.eigvalsh(matrix)
        uncertainty = rho/s**2*4*np.dot(q, q)/self.gap_ratio**2*error
        uncertainty += 64*np.finfo(float).eps*max(1., np.linalg.norm(matrix, 2))
        stable = bool(np.all(np.isfinite(matrix)) and np.isfinite(uncertainty) and eigen[0] > uncertainty)
        return MixedPrincipalSymbol(matrix, eigen, float(uncertainty), stable, float(a), float(gamma),
                                    dimensions, steps, estimates)

    def principal_symbol(self, delta_bar, gradient_bar, occupation, *, spatial_dimensions=2):
        """Local gradient Hessian in order (dx_Re,dx_Im,dy_Re,dy_Im)."""
        if spatial_dimensions not in (1, 2):
            raise ValueError("one or two spatial dimensions required")
        z = self._complex(delta_bar, "delta_bar")
        d = np.asarray(gradient_bar, dtype=complex)
        if d.shape != (spatial_dimensions,) or np.any(~np.isfinite(d)):
            raise ValueError("one finite complex derivative per spatial dimension required")
        raw = np.asarray(occupation)
        if np.iscomplexobj(raw):
            raise ValueError("occupations must be real")
        p = self._population(raw)
        q = np.imag(np.conj(z)*d)/(abs(z)**2+self.delta_regularizer_bar**2)
        return self._mixed_symbol(z, q, float(np.dot(q, q)/self.gap_ratio), p, {}, spatial_dimensions)

    def evaluate(self, delta_bar, p_quadrature, *, link_phases=None, require_stability=True, on_quadrature=None):
        return self._evaluate(delta_bar, p_quadrature, None, link_phases, require_stability, on_quadrature)

    def evaluate_thermal(self, delta_bar, bath_theta, *, link_phases=None,
                         require_stability=True, on_quadrature=None):
        """Thermal FREE energy gradient; occupations are prepared at each sector spectrum."""
        if not np.isfinite(bath_theta) or bath_theta < 0:
            raise ValueError("finite nonnegative bath temperature required")
        return self._evaluate(delta_bar, None, float(bath_theta), link_phases, require_stability, on_quadrature)

    def _evaluate(self, delta_bar, p_quadrature, theta, link_phases, require_stability, callback):
        z, zq, derivatives, stiffnesses = self._operators(delta_bar, link_phases)
        fields = self._sample(zq, derivatives)
        shape = (self.quadrature_size, len(self.catalog.count_nodes))
        if theta is None:
            raw = np.asarray(p_quadrature)
            if np.iscomplexobj(raw):
                raise ValueError("quadrature populations must be real")
            p = np.array(raw, dtype=float, copy=True)
            if p.shape != shape or np.any(~np.isfinite(p)) or np.any((p < 0)|(p > 1)):
                raise ValueError("physical sector-quadrature populations required; no implicit DOF remapping")
        else:
            p = np.zeros(shape)
        cache, electronic, free = {}, [], []
        symbols = [] if require_stability else None
        for i in range(self.quadrature_size):
            a, gamma = fields.amplitude_quadrature_bar[i], fields.gamma_quadrature_bar[i]
            if theta is not None and theta > 0:
                key = (float(a), float(gamma))
                if key not in cache:
                    cache[key] = [np.asarray(self.catalog.vacuum.evaluate(a, gamma)), None]
                if cache[key][1] is None:
                    cache[key][1] = tuple(np.asarray(k) for k in self.catalog.energy_kernel(a, gamma))
                vacuum, kernels = cache[key]
                exponent = np.exp(-kernels[0]/theta)
                p[i] = exponent/(1+exponent)
                free.append(float(vacuum[0]-4*theta*np.dot(self.catalog.count_weights,
                                                          np.logaddexp(0., -kernels[0]/theta))))
            row = self._potential(a, gamma, p[i], cache)
            electronic.append(row)
            if theta == 0:
                free.append(row[0])
            if require_stability:
                symbol = self._mixed_symbol(zq[i], fields.q_delta_quadrature_bar[i], gamma, p[i], cache,
                                            int(self.spatial_dimensions[i]))
                symbols.append(symbol)
                if not symbol.stable:
                    raise SpatialAdmissibilityError(i, symbol)
            if callback is not None:
                callback(i, self.quadrature_size)
        electronic = np.asarray(electronic)
        mass = self.quadrature_mass_bar
        rho = abs(zq)**2
        q = fields.q_delta_quadrature_bar
        v = (2*q*electronic[:, 2, None]/self.gap_ratio-2*self.kappa*rho[:, None]*q)/(rho[:, None]+.01)
        radial = np.zeros_like(zq)
        np.divide(electronic[:, 1]*zq, abs(zq), out=radial, where=abs(zq) != 0)
        directional = []
        for axis in range(2):
            qa, va, d = q[:, axis], v[:, axis], fields.derivative_quadrature_bar[:, axis]
            gz = -2*self.kappa*qa*qa*zq-1j*va*d-2*va*qa*zq
            gd = 1j*va*zq
            directional.append(mass*gz+derivatives[axis].conj().T@(mass*gd)+2*self.kappa*(stiffnesses[axis]@zq))
        raw_gradient = mass*radial+directional[0]+directional[1]
        gradient = np.asarray(self.prolongation.T@raw_gradient)
        current = np.zeros(len(self.graph_edges))
        line_residuals = []
        for indices, axis, edge_ids in self._lines:
            gphase = np.imag(np.conj(zq[indices])*directional[axis][indices])
            line_residuals.append(float(np.sum(gphase)))
            if len(edge_ids):
                current[edge_ids] = -np.cumsum(gphase)[:-1]
        phase_gradient = np.imag(np.conj(z)*gradient)
        noether = phase_gradient+self.incidence@current
        stiffness_energy = sum(np.real(np.vdot(zq, k@zq)) for k in stiffnesses)
        remainder = float(self.kappa*(stiffness_energy-np.dot(mass, rho*np.sum(q*q, axis=1))))
        energy = float(np.dot(mass, electronic[:, 0])+remainder)
        free_energy = None if theta is None else float(np.dot(mass, free)+remainder)
        left = self.sector_slices["left"]
        right = self.sector_slices["right"]
        reactions = dict(
            left_lead=raw_gradient[left.stop-1], left_rectangle=np.sum(raw_gradient[self.rectangle_quadrature[0]]),
            right_rectangle=np.sum(raw_gradient[self.rectangle_quadrature[-1]]), right_lead=raw_gradient[right.start],
            interpretation="Restricted integrated force contributions, including endpoint quadrature local terms; not isolated pointwise continuum tractions. Their sums are the shared DOF forces.")
        return MixedSpatialEvaluation(energy, self.energy_scale_J*energy, free_energy,
            np.column_stack((gradient.real, gradient.imag)), np.column_stack((raw_gradient.real, raw_gradient.imag)),
            current, (2*E_CHARGE_C/HBAR_J_S)*self.energy_scale_J*current, fields, electronic, p,
            tuple(symbols) if symbols is not None else None, remainder, np.asarray(noether),
            float(np.sum(phase_gradient)), np.asarray(line_residuals), reactions)
