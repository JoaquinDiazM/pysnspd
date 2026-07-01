"""pyTDGL-like TDGL solver core for OE7 stationary comparisons.

The public names and method signatures in this module intentionally mirror
``tdgl.solver.solver``.  The implementation is reduced to the no-screening,
CPU/SuperLU path required for a stationary SNSPD comparison backend, while the
local nonlinear update is where pySNSPD can substitute its modified
``w_i^n``/``z_i^n`` physics.
"""
from __future__ import annotations

import inspect
import itertools
import logging
import numbers
import sys
from dataclasses import dataclass, field
from datetime import datetime
from typing import Callable, Dict, List, NamedTuple, Optional, Sequence, Tuple, Union

import numpy as np
import scipy.sparse as sp

from .device import PySNSPDTDGLDevice as Device, TerminalInfo
from .tdgl_operators import MeshOperators
from .options import SolverOptions, SparseSolver

logger = logging.getLogger("pytdgl_like_solver")


def validate_terminal_currents(
    terminal_currents: Union[Callable, Dict[str, float]],
    terminal_info: Sequence[TerminalInfo],
    solver_options: SolverOptions,
    num_evals: int = 100,
) -> None:
    """Ensure that the terminal currents satisfy current conservation."""

    def check_total_current(currents: Dict[str, float]):
        names = set([t.name for t in terminal_info])
        unknown = set(currents).difference(names)
        if unknown:
            raise ValueError(f"Unknown terminal(s) in terminal currents: {list(unknown)}.")
        total_current = sum(currents.values())
        if abs(total_current) > 1.0e-13:
            raise ValueError(
                f"The sum of all terminal currents must be 0 (got {total_current:.2e})."
            )

    if callable(terminal_currents):
        times = np.random.default_rng(12345).random(num_evals) * solver_options.solve_time
        for t in times:
            check_total_current(terminal_currents(float(t)))
    else:
        check_total_current(terminal_currents)


class SolverResult(NamedTuple):
    """A container for the results of a single solve step."""

    dt: float
    psi: np.ndarray
    mu: np.ndarray
    supercurrent: np.ndarray
    normal_current: np.ndarray
    A_induced: np.ndarray
    A_applied: Optional[np.ndarray] = None
    epsilon: Optional[np.ndarray] = None


@dataclass
class RunningState:
    data: dict[str, list[np.ndarray | float]] = field(default_factory=dict)

    def append(self, name: str, value) -> None:
        self.data.setdefault(name, []).append(np.array(value, copy=True) if np.ndim(value) else float(value))


@dataclass
class PyTDGLLikeSolution:
    device: Device
    options: SolverOptions
    tdgl_data: SolverResult
    history: dict[str, np.ndarray]
    total_seconds: float


class TDGLSolver:
    """Solver for a TDGL model."""

    def __init__(
        self,
        device: Device,
        options: SolverOptions,
        applied_vector_potential: Union[Callable, float] = 0.0,
        terminal_currents: Union[Callable, Dict[str, float], None] = None,
        disorder_epsilon: Union[Callable, float] = 1.0,
        seed_solution: Optional[object] = None,
        progress: bool = False,
        supercurrent_override: Optional[Callable[[np.ndarray, np.ndarray], np.ndarray]] = None,
        supercurrent_law: str = "gl",
        allmaras_forcing_callback: Optional[Callable[[np.ndarray, sp.spmatrix], np.ndarray]] = None,
        stop_eta: Optional[float] = None,
        stop_min_steps: int = 0,
        stop_on_convergence: bool = False,
    ):
        self.device = device
        self.options = options
        self.options.validate()
        self.terminal_currents = terminal_currents
        self.seed_solution = seed_solution
        self.progress = bool(progress)
        self.supercurrent_override = supercurrent_override
        self.supercurrent_law = str(supercurrent_law)
        self.allmaras_forcing_callback = allmaras_forcing_callback
        self.last_allmaras_forcing_dimensionless = None
        self.stop_eta = None if stop_eta is None else float(stop_eta)
        self.stop_min_steps = max(0, int(stop_min_steps))
        self.stop_on_convergence = bool(stop_on_convergence)
        self.converged = False
        self.convergence_reason = "not_evaluated"
        self.eta_converged = False
        self.eta_convergence_step = -1
        self.eta_convergence_time = float("nan")
        self.stop_reason = "not_started"
        self.xp = np
        self.use_cupy = False

        mesh = self.device.mesh
        edges = mesh.edge_mesh.edges
        self.num_edges = len(edges)
        self.normalized_directions = mesh.edge_mesh.normalized_directions
        self.sites = mesh.sites
        self.edge_centers = mesh.edge_mesh.centers
        self.z0 = np.zeros(len(self.edge_centers), dtype=float)
        self.applied_vector_potential = applied_vector_potential
        self.dynamic_vector_potential = False

        current_A_applied = np.zeros_like(mesh.edge_mesh.directions, dtype=float)

        if callable(disorder_epsilon):
            argspec = inspect.getfullargspec(disorder_epsilon)
            self.dynamic_epsilon = "t" in argspec.kwonlyargs
            self.vectorized_epsilon = (
                argspec.kwonlydefaults is not None
                and argspec.kwonlydefaults.get("vectorized", False)
            )
        else:
            _disorder_epsilon = float(disorder_epsilon)

            def disorder_epsilon(r):
                return _disorder_epsilon * np.ones(len(r), dtype=float)

            self.vectorized_epsilon = True
            self.dynamic_epsilon = False
        self.disorder_epsilon = disorder_epsilon
        epsilon = disorder_epsilon(self.sites) if self.vectorized_epsilon else np.array([float(disorder_epsilon(r)) for r in self.sites])
        if np.any(epsilon > 1.5):
            logger.warning("epsilon contains values > 1; continuing for pySNSPD diagnostics.")

        self.terminal_info = device.terminal_info()
        self.terminal_names = [term.name for term in self.terminal_info]
        for term_info in self.terminal_info:
            if term_info.length == 0:
                raise ValueError(
                    f"Terminal {term_info.name!r} does not contain any points on the boundary of the mesh."
                )

        if terminal_currents is None:
            terminal_currents = {name: 0.0 for name in self.terminal_names}
        if callable(terminal_currents):
            current_func = terminal_currents
        else:
            terminal_currents = {name: float(terminal_currents.get(name, 0.0)) for name in self.terminal_names}

            def current_func(t):
                return terminal_currents

        self.current_func = current_func
        validate_terminal_currents(self.current_func, self.terminal_info, self.options)
        terminal_indices = [t.site_indices for t in self.terminal_info]
        if terminal_indices:
            normal_boundary_index = np.concatenate(terminal_indices).astype(np.int64, copy=False)
        else:
            normal_boundary_index = np.array([], dtype=np.int64)
        normal_boundary_index = np.unique(normal_boundary_index).astype(np.int64, copy=False)
        self.normal_boundary_index = normal_boundary_index
        self.terminal_current_densities = {name: 0.0 for name in self.terminal_names}

        terminal_psi = options.terminal_psi
        self.terminal_psi_value = None if terminal_psi is None else complex(terminal_psi)
        operators = MeshOperators(
            mesh,
            options.sparse_solver,
            use_cupy=False,
            fixed_sites=normal_boundary_index,
            fix_psi=(terminal_psi is not None),
        )
        operators.build_operators()
        operators.set_link_exponents(current_A_applied)
        self.operators = operators

        psi_init = np.ones(len(mesh.sites), dtype=np.complex128)
        psi_init = self.apply_terminal_psi(psi_init)
        mu_init = np.zeros(len(mesh.sites), dtype=float)
        mu_boundary = np.zeros_like(mesh.edge_mesh.boundary_edge_indices, dtype=float)

        self.psi_init = psi_init
        self.mu_init = mu_init
        self.epsilon = np.asarray(epsilon, dtype=float)
        self.mu_boundary = mu_boundary
        self.current_A_applied = current_A_applied
        self.new_A_induced = None
        self.areas = None
        self.d_psi_sq_vals: list[float] = []
        self.tentative_dt = options.dt_init
        self.dt_max = options.dt_max if options.adaptive else options.dt_init
        # Adaptive-Euler diagnostics. pyTDGL shrinks the tentative step when
        # the local algebraic |psi|^2 solve fails; these fields expose that
        # retry logic without changing the public adaptive_euler_step return
        # signature.
        self.last_adaptive_dt_attempt = float(options.dt_init)
        self.last_adaptive_dt_accepted = float(options.dt_init)
        self.last_adaptive_retries = 0
        self.last_adaptive_rejected_attempts = 0
        self.last_adaptive_target_dt = float("nan")
        self.last_adaptive_next_dt = float(options.dt_init)
        self.last_adaptive_window_mean_d_abs_sq = float("nan")

    def apply_terminal_psi(self, psi: np.ndarray) -> np.ndarray:
        """Apply the metallic-normal-terminal Dirichlet condition to psi.

        For pySNSPD SNSPD windows the left/right terminals represent normal
        metallic contacts.  The scalar-potential problem still injects the
        imposed current through a Neumann condition, but the superconducting
        order parameter must be clamped on the terminal sites,

            psi_i = terminal_psi,  i in Gamma_N.

        With the default ``terminal_psi=0`` this enforces Delta=0 at the
        metallic terminals.  The operation is intentionally repeated on the
        initial seed and after every local nonlinear update because externally
        supplied seeds are not guaranteed to satisfy the terminal constraint.
        """

        arr = np.asarray(psi, dtype=np.complex128)
        if self.terminal_psi_value is not None and self.normal_boundary_index.size:
            arr[self.normal_boundary_index] = self.terminal_psi_value
        return arr

    def update_mu_boundary(self, time: float) -> None:
        """Compute terminal Neumann values for the scalar-potential solve.

        ``self.current_func`` returns pySNSPD terminal currents in SI amperes.
        The device adapter converts those physical currents to the internal
        dimensionless boundary derivative used by the pyTDGL-like Poisson
        operator.  This keeps the solver API comparable to pyTDGL while avoiding
        the previous artificial ``I_norm = 1`` current normalization.
        """

        currents_A = {name: float(value) for name, value in self.current_func(time).items()}
        terminal_current_densities = self.terminal_current_densities
        for terminal in self.terminal_info:
            mu_boundary_value = self.device.terminal_mu_boundary_value(
                terminal=terminal,
                terminal_currents_A=currents_A,
            )
            if mu_boundary_value != terminal_current_densities[terminal.name]:
                terminal_current_densities[terminal.name] = mu_boundary_value
                # ``mu_boundary`` is indexed by the compact boundary-edge vector.
                bmap = self.device.mesh.edge_mesh.boundary_edge_indices
                pos = {int(edge): k for k, edge in enumerate(np.asarray(bmap, dtype=int))}
                for edge in terminal.boundary_edge_indices:
                    k = pos.get(int(edge))
                    if k is not None:
                        self.mu_boundary[k] = mu_boundary_value

    def update_applied_vector_potential(self, time: float) -> np.ndarray:
        """Evaluates the time-dependent vector potential."""

        return np.zeros_like(self.current_A_applied)

    def update_epsilon(self, time: float) -> np.ndarray:
        """Evaluates the time-dependent disorder parameter epsilon."""

        if self.vectorized_epsilon:
            epsilon = self.disorder_epsilon(self.sites, t=time)
        else:
            epsilon = np.array([float(self.disorder_epsilon(r, t=time)) for r in self.sites])
        return np.asarray(epsilon, dtype=float)

    @staticmethod
    def solve_for_psi_squared(
        *,
        psi: np.ndarray,
        abs_sq_psi: np.ndarray,
        mu: np.ndarray,
        epsilon: np.ndarray,
        gamma: float,
        u: float,
        dt: float,
        psi_laplacian: sp.spmatrix,
        forcing_dimensionless: Optional[np.ndarray] = None,
    ) -> Union[Tuple[np.ndarray, np.ndarray], None]:
        """Solves for psi^{n+1} and |psi^{n+1}|^2.

        This is the pyTDGL local algebraic solve with the same public arguments.
        In this backend, ``epsilon``, ``gamma`` and ``u`` are supplied by the
        pySNSPD adapter.  If ``forcing_dimensionless`` is supplied, it replaces
        the native GL bracket and injects the Appendix-B Allmaras functional
        directly into the same local ``w_i, z_i`` solve.
        """

        xp = np
        U = xp.exp(-1j * mu * dt)
        z = U * gamma**2 / 2 * psi
        with np.errstate(all="raise"):
            try:
                if forcing_dimensionless is None:
                    forcing = (epsilon - abs_sq_psi) * psi + psi_laplacian @ psi
                else:
                    forcing = xp.asarray(forcing_dimensionless, dtype=xp.complex128)
                    if forcing.shape != psi.shape:
                        raise ValueError(
                            "forcing_dimensionless must have the same shape as psi, "
                            f"got {forcing.shape} and {psi.shape}."
                        )
                w = z * abs_sq_psi + U * (
                    psi
                    + (dt / u)
                    * xp.sqrt(1 + gamma**2 * abs_sq_psi)
                    * forcing
                )
                c = w.real * z.real + w.imag * z.imag
                two_c_1 = 2 * c + 1
                w2 = xp.absolute(w) ** 2
                discriminant = two_c_1**2 - 4 * xp.absolute(z) ** 2 * w2
            except Exception:
                logger.warning("Unable to solve for |psi|^2.", exc_info=True)
                return None
        if xp.any(discriminant < -1.0e-13):
            return None
        discriminant = xp.maximum(discriminant, 0.0)
        new_sq_psi = (2 * w2) / (two_c_1 + xp.sqrt(discriminant))
        psi = w - z * new_sq_psi
        if not np.all(np.isfinite(psi)) or not np.all(np.isfinite(new_sq_psi)):
            return None
        return psi, new_sq_psi

    def adaptive_euler_step(
        self,
        step: int,
        psi: np.ndarray,
        abs_sq_psi: np.ndarray,
        mu: np.ndarray,
        epsilon: np.ndarray,
        dt: float,
    ) -> Tuple[np.ndarray, np.ndarray, float]:
        """Updates the order parameter and time step in an adaptive Euler step."""

        options = self.options
        forcing_dimensionless = None
        if self.allmaras_forcing_callback is not None:
            forcing_dimensionless = self.allmaras_forcing_callback(
                np.asarray(psi, dtype=np.complex128),
                self.operators.psi_laplacian,
            )
            forcing_dimensionless = np.asarray(forcing_dimensionless, dtype=np.complex128)
            if forcing_dimensionless.shape != psi.shape:
                raise ValueError(
                    "allmaras_forcing_callback returned shape "
                    f"{forcing_dimensionless.shape}, expected {psi.shape}."
                )
            self.last_allmaras_forcing_dimensionless = forcing_dimensionless.copy()

        kwargs = dict(
            psi=psi,
            abs_sq_psi=abs_sq_psi,
            mu=mu,
            epsilon=epsilon,
            gamma=self.device.layer.gamma,
            u=self.device.layer.u,
            dt=dt,
            psi_laplacian=self.operators.psi_laplacian,
            forcing_dimensionless=forcing_dimensionless,
        )
        self.last_adaptive_dt_attempt = float(dt)
        self.last_adaptive_dt_accepted = float(dt)
        self.last_adaptive_retries = 0
        self.last_adaptive_rejected_attempts = 0

        result = self.solve_for_psi_squared(**kwargs)
        retries_used = 0
        for retries in itertools.count():
            if result is not None:
                retries_used = int(retries)
                break
            if not options.adaptive or retries > options.max_solve_retries:
                raise RuntimeError(
                    f"Solver failed to converge in {options.max_solve_retries} retries "
                    f"at step {step} with dt = {dt:.2e}. Try using a smaller dt_init."
                )
            kwargs["dt"] = dt = dt * options.adaptive_time_step_multiplier
            result = self.solve_for_psi_squared(**kwargs)

        self.last_adaptive_dt_accepted = float(dt)
        self.last_adaptive_retries = int(retries_used)
        self.last_adaptive_rejected_attempts = int(retries_used)
        psi, new_sq_psi = result
        psi = self.apply_terminal_psi(psi)
        if self.terminal_psi_value is not None and self.normal_boundary_index.size:
            new_sq_psi = np.asarray(new_sq_psi, dtype=float)
            new_sq_psi[self.normal_boundary_index] = np.abs(self.terminal_psi_value) ** 2
        return psi, new_sq_psi, dt

    def solve_for_observables(
        self,
        psi: np.ndarray,
        dA_dt: Union[float, np.ndarray],
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Solves for mu, supercurrent and normal current.

        This follows the pyTDGL no-screening CPU order and also stores the
        native Poisson terms used internally by the solver.  These diagnostics
        are deliberately kept in the same dimensionless operator space as the
        pyTDGL-like sparse system so we can distinguish a true solver residual
        from a pySNSPD adapter/plotting convention mismatch.
        """

        operators = self.operators
        psi = self.apply_terminal_psi(np.asarray(psi, dtype=np.complex128).copy())
        gl_supercurrent = operators.get_supercurrent(psi)
        supercurrent = gl_supercurrent
        if self.supercurrent_override is not None:
            supercurrent = self.supercurrent_override(psi, gl_supercurrent)
            supercurrent = np.asarray(supercurrent, dtype=float)
            if supercurrent.shape != gl_supercurrent.shape:
                raise ValueError(
                    "supercurrent_override returned shape "
                    f"{supercurrent.shape}, expected {gl_supercurrent.shape}."
                )
            if not np.all(np.isfinite(supercurrent)):
                raise ValueError("supercurrent_override returned non-finite values.")

        div_supercurrent = operators.divergence @ (supercurrent - dA_dt)
        boundary_rhs = operators.mu_boundary_laplacian @ self.mu_boundary
        rhs = div_supercurrent - boundary_rhs
        mu = operators.mu_laplacian_lu(rhs)
        lhs = operators.mu_laplacian @ mu
        poisson_residual = lhs - rhs
        normal_current = -(operators.mu_gradient @ mu) - dA_dt

        # Native pyTDGL-like diagnostics.  These arrays are not converted to SI;
        # they live in the exact linear system solved above.  ``last_supercurrent``
        # is the current actually used in Poisson.  ``last_gl_supercurrent`` is
        # always the native pyTDGL/GL current reconstructed from psi.
        self.last_gl_supercurrent = np.asarray(np.real_if_close(gl_supercurrent, tol=1000), dtype=float)
        self.last_supercurrent = np.asarray(np.real_if_close(supercurrent, tol=1000), dtype=float)
        self.last_normal_current = np.asarray(np.real_if_close(normal_current, tol=1000), dtype=float)
        self.last_div_supercurrent = np.asarray(np.real_if_close(div_supercurrent, tol=1000), dtype=float)
        self.last_boundary_rhs = np.asarray(np.real_if_close(boundary_rhs, tol=1000), dtype=float)
        self.last_poisson_rhs = np.asarray(np.real_if_close(rhs, tol=1000), dtype=float)
        self.last_poisson_lhs = np.asarray(np.real_if_close(lhs, tol=1000), dtype=float)
        self.last_poisson_residual = np.asarray(np.real_if_close(poisson_residual, tol=1000), dtype=float)
        self.last_mu_boundary = np.asarray(self.mu_boundary, dtype=float).copy()

        mu = np.real_if_close(mu, tol=1000)
        supercurrent = np.real_if_close(supercurrent, tol=1000)
        normal_current = np.real_if_close(normal_current, tol=1000)
        return np.asarray(mu, dtype=float), np.asarray(supercurrent, dtype=float), np.asarray(normal_current, dtype=float)

    def get_induced_vector_potential(
        self,
        current_density: np.ndarray,
        A_induced_vals: List[np.ndarray],
        velocity: List[np.ndarray],
    ) -> Tuple[np.ndarray, float]:
        """Placeholder matching pyTDGL's method name; screening is disabled."""

        return A_induced_vals[-1], 0.0

    def update(
        self,
        state: Dict[str, numbers.Real],
        running_state: RunningState,
        dt: float,
        *,
        psi: np.ndarray,
        mu: np.ndarray,
        supercurrent: np.ndarray,
        normal_current: np.ndarray,
        induced_vector_potential: np.ndarray,
        applied_vector_potential: Optional[np.ndarray] = None,
        epsilon: Optional[np.ndarray] = None,
    ) -> SolverResult:
        """Called at each time step to update the state of the system."""

        options = self.options
        operators = self.operators
        step = int(state["step"])
        time = float(state["time"])

        self.update_mu_boundary(time)
        dA_dt = 0.0
        current_A_applied = self.current_A_applied
        operators.set_link_exponents(current_A_applied)

        if self.dynamic_epsilon:
            self.epsilon = self.update_epsilon(time)
        epsilon = self.epsilon

        old_sq_psi = np.absolute(psi) ** 2
        dt = self.tentative_dt if options.adaptive else dt
        psi, abs_sq_psi, dt = self.adaptive_euler_step(
            step, psi, old_sq_psi, mu, epsilon, dt
        )
        mu, supercurrent, normal_current = self.solve_for_observables(psi, dA_dt)

        running_state.append("dt", dt)
        running_state.append("max_abs_psi", float(np.max(np.abs(psi))))
        running_state.append("min_abs_psi", float(np.min(np.abs(psi))))
        running_state.append("max_d_abs_sq_psi", float(np.max(np.abs(abs_sq_psi - old_sq_psi))))
        running_state.append("mu_ptp", float(np.ptp(mu)))
        running_state.append("max_supercurrent", float(np.max(np.abs(supercurrent))) if supercurrent.size else 0.0)
        running_state.append("max_normal_current", float(np.max(np.abs(normal_current))) if normal_current.size else 0.0)
        forcing = getattr(self, "last_allmaras_forcing_dimensionless", None)
        running_state.append(
            "allmaras_update_forcing_max_abs",
            float(np.max(np.abs(forcing))) if forcing is not None and np.size(forcing) else 0.0,
        )
        rhs = getattr(self, "last_poisson_rhs", np.array([], dtype=float))
        residual = getattr(self, "last_poisson_residual", np.array([], dtype=float))
        div_s = getattr(self, "last_div_supercurrent", np.array([], dtype=float))
        b_rhs = getattr(self, "last_boundary_rhs", np.array([], dtype=float))
        mu_b = getattr(self, "last_mu_boundary", np.array([], dtype=float))
        rhs_norm = float(np.linalg.norm(rhs)) if rhs.size else 0.0
        res_norm = float(np.linalg.norm(residual)) if residual.size else 0.0
        running_state.append("poisson_rhs_norm", rhs_norm)
        running_state.append("poisson_residual_norm", res_norm)
        running_state.append("poisson_residual_rel", res_norm / max(rhs_norm, 1.0e-300))
        running_state.append("poisson_residual_max_abs", float(np.max(np.abs(residual))) if residual.size else 0.0)
        running_state.append("div_supercurrent_norm", float(np.linalg.norm(div_s)) if div_s.size else 0.0)
        running_state.append("boundary_rhs_norm", float(np.linalg.norm(b_rhs)) if b_rhs.size else 0.0)
        running_state.append("mu_boundary_max_abs", float(np.max(np.abs(mu_b))) if mu_b.size else 0.0)

        if options.adaptive:
            d_abs_sq = float(np.absolute(abs_sq_psi - old_sq_psi).max())
            self.d_psi_sq_vals.append(d_abs_sq)
            window = max(1, int(options.adaptive_window))
            mean_d_abs_sq = float(np.mean(self.d_psi_sq_vals[-window:]))

            if step > window:
                target_dt = float(options.dt_init / max(1e-10, mean_d_abs_sq))
            else:
                target_dt = float("nan")

            # The first adaptive diagnostic version let a tiny mean change push
            # the tentative step directly to dt_max after every accepted step.
            # For this SS problem that wasted time: the local cubic solve then
            # rejected the same oversized step on essentially every iteration.
            # Grow toward the pyTDGL window target, but cap per-step growth so
            # the next tentative dt stays close to recently accepted values.
            growth_factor = float(getattr(options, "adaptive_growth_factor", 1.5))
            growth_factor = max(1.0, growth_factor)
            desired_dt = self.dt_max if not np.isfinite(target_dt) else min(float(target_dt), self.dt_max)
            next_dt = float(min(desired_dt, max(float(dt), float(dt) * growth_factor)))
            next_dt = float(np.clip(next_dt, options.dt_init, self.dt_max))

            self.last_adaptive_window_mean_d_abs_sq = mean_d_abs_sq
            self.last_adaptive_target_dt = target_dt
            self.last_adaptive_next_dt = next_dt
            self.tentative_dt = next_dt
        else:
            self.last_adaptive_window_mean_d_abs_sq = float("nan")
            self.last_adaptive_target_dt = float("nan")
            self.last_adaptive_next_dt = float(dt)

        running_state.append("dt_attempt", float(self.last_adaptive_dt_attempt))
        running_state.append("dt_accepted", float(self.last_adaptive_dt_accepted))
        running_state.append("dt_next", float(self.last_adaptive_next_dt))
        running_state.append("adaptive_retries", int(self.last_adaptive_retries))
        running_state.append("adaptive_rejected_attempts", int(self.last_adaptive_rejected_attempts))
        running_state.append("adaptive_window_mean_d_abs_sq", float(self.last_adaptive_window_mean_d_abs_sq))
        running_state.append("adaptive_target_dt", float(self.last_adaptive_target_dt))

        return SolverResult(dt, psi, mu, supercurrent, normal_current, induced_vector_potential)

    def solve(self) -> Optional[PyTDGLLikeSolution]:
        """Runs the solver."""

        start_time = datetime.now()
        options = self.options
        options.validate()
        num_edges = self.num_edges

        if self.seed_solution is None:
            parameters = {
                "psi": self.psi_init.copy(),
                "mu": self.mu_init.copy(),
                "supercurrent": np.zeros(num_edges, dtype=float),
                "normal_current": np.zeros(num_edges, dtype=float),
                "induced_vector_potential": np.zeros((num_edges, 2), dtype=float),
            }
        else:
            seed_data = getattr(self.seed_solution, "tdgl_data", self.seed_solution)

            def _seed_get(name, default):
                if isinstance(seed_data, dict):
                    return seed_data.get(name, default)
                return getattr(seed_data, name, default)

            parameters = {
                "psi": np.asarray(_seed_get("psi", self.psi_init), dtype=np.complex128).copy(),
                "mu": np.asarray(_seed_get("mu", self.mu_init), dtype=float).copy(),
                "supercurrent": np.asarray(_seed_get("supercurrent", np.zeros(num_edges)), dtype=float).copy(),
                "normal_current": np.asarray(_seed_get("normal_current", np.zeros(num_edges)), dtype=float).copy(),
                "induced_vector_potential": np.zeros((num_edges, 2), dtype=float),
            }

        parameters["psi"] = self.apply_terminal_psi(parameters["psi"])

        running_state = RunningState()
        state: dict[str, numbers.Real] = {"step": 0, "time": 0.0}
        dt = float(options.dt_init)
        self.update_mu_boundary(0.0)
        mu0, supercurrent0, normal_current0 = self.solve_for_observables(parameters["psi"], 0.0)
        result = SolverResult(
            dt,
            parameters["psi"],
            mu0,
            supercurrent0,
            normal_current0,
            parameters["induced_vector_potential"],
        )

        # Minimal trajectory snapshots for pySNSPD diagnostics.  pyTDGL writes
        # full HDF5 frames; this lightweight backend only stores a small fixed
        # number of in-memory frames so the existing OE7 plotting utilities can
        # inspect the actual evolution instead of repeated final states.
        snapshot_count = max(2, int(getattr(self, "snapshot_count", 2)))
        snapshot_times = np.linspace(0.0, float(options.solve_time), snapshot_count)
        next_snapshot = 0

        def append_snapshot(time_value: float, frame: SolverResult) -> None:
            running_state.append("snapshot_t", float(time_value))
            running_state.append("psi_snapshot", frame.psi)
            running_state.append("mu_snapshot", frame.mu)
            running_state.append("supercurrent_snapshot", frame.supercurrent)
            running_state.append("gl_supercurrent_snapshot", getattr(self, "last_gl_supercurrent", frame.supercurrent))
            running_state.append("normal_current_snapshot", frame.normal_current)
            running_state.append("poisson_rhs_snapshot", getattr(self, "last_poisson_rhs", np.zeros(len(self.sites))))
            running_state.append("poisson_lhs_snapshot", getattr(self, "last_poisson_lhs", np.zeros(len(self.sites))))
            running_state.append("poisson_residual_snapshot", getattr(self, "last_poisson_residual", np.zeros(len(self.sites))))
            running_state.append("div_supercurrent_snapshot", getattr(self, "last_div_supercurrent", np.zeros(len(self.sites))))
            running_state.append("boundary_rhs_snapshot", getattr(self, "last_boundary_rhs", np.zeros(len(self.sites))))
            running_state.append("mu_boundary_snapshot", getattr(self, "last_mu_boundary", np.zeros(len(self.device.mesh.edge_mesh.boundary_edge_indices))))

        append_snapshot(0.0, result)
        next_snapshot = 1

        progress_next_fraction = 0.0

        def emit_progress(final: bool = False) -> None:
            nonlocal progress_next_fraction
            if not self.progress:
                return
            total = max(float(options.solve_time), 1.0e-300)
            frac = min(1.0, max(0.0, float(state["time"]) / total))
            if (not final) and frac < progress_next_fraction:
                return
            width = 32
            filled = int(round(width * frac))
            bar = "#" * filled + "-" * (width - filled)
            try:
                tau0 = float(getattr(self.device.material, "tau0_GL_s", 1.0))
            except Exception:
                tau0 = 1.0
            t_ps = float(state["time"]) * tau0 / 1.0e-12
            sys.stderr.write(
                f"\rSS pyTDGL-like [{bar}] {100.0 * frac:6.2f}% "
                f"step={int(state['step'])} t={t_ps:.4g} ps"
            )
            sys.stderr.flush()
            progress_next_fraction = min(1.0, frac + 0.01)
            if final:
                sys.stderr.write("\n")
                sys.stderr.flush()

        emit_progress()

        while float(state["time"]) < options.solve_time:
            result = self.update(
                state,
                running_state,
                dt,
                psi=result.psi,
                mu=result.mu,
                supercurrent=result.supercurrent,
                normal_current=result.normal_current,
                induced_vector_potential=result.A_induced,
            )
            dt = float(result.dt)
            state["time"] = float(state["time"]) + dt
            state["step"] = int(state["step"]) + 1

            while next_snapshot < snapshot_times.size and float(state["time"]) >= float(snapshot_times[next_snapshot]):
                append_snapshot(float(state["time"]), result)
                next_snapshot += 1

            emit_progress()

            if self.stop_eta is not None and int(state["step"]) >= self.stop_min_steps:
                eta_step = float(running_state.data.get("max_d_abs_sq_psi", [float("inf")])[-1])
                if np.isfinite(eta_step) and eta_step <= float(self.stop_eta):
                    if not self.eta_converged:
                        self.eta_converged = True
                        self.eta_convergence_step = int(state["step"])
                        self.eta_convergence_time = float(state["time"])
                    self.converged = True
                    self.convergence_reason = f"max_d_abs_sq_psi<={float(self.stop_eta):.3e}"
                    if self.stop_on_convergence:
                        self.stop_reason = "eta_stop"
                        break

            if int(state["step"]) > 10_000_000:
                self.stop_reason = "safety_step_limit"
                raise RuntimeError("pytdgl_like solve exceeded 10,000,000 steps.")

        if self.stop_reason == "not_started":
            if float(state["time"]) >= float(options.solve_time):
                self.stop_reason = "requested_time_reached"
            else:
                self.stop_reason = "loop_exited_before_requested_time"

        emit_progress(final=True)

        if next_snapshot < snapshot_times.size:
            append_snapshot(float(state["time"]), result)

        end_time = datetime.now()
        hist = {key: np.asarray(vals) for key, vals in running_state.data.items()}
        hist["final_step"] = np.array([int(state["step"])], dtype=int)
        hist["final_time"] = np.array([float(state["time"])], dtype=float)
        hist["supercurrent_law"] = np.array([self.supercurrent_law], dtype=object)
        hist["converged"] = np.array([bool(self.converged)], dtype=bool)
        hist["convergence_reason"] = np.array([self.convergence_reason], dtype=object)
        hist["eta_converged"] = np.array([bool(self.eta_converged)], dtype=bool)
        hist["eta_convergence_step"] = np.array([int(self.eta_convergence_step)], dtype=int)
        hist["eta_convergence_time"] = np.array([float(self.eta_convergence_time)], dtype=float)
        hist["stop_reason"] = np.array([self.stop_reason], dtype=object)
        hist["stop_on_convergence"] = np.array([bool(self.stop_on_convergence)], dtype=bool)
        hist["stop_eta"] = np.array([float(self.stop_eta) if self.stop_eta is not None else np.nan], dtype=float)
        hist["stop_min_steps"] = np.array([int(self.stop_min_steps)], dtype=int)
        hist["adaptive_enabled"] = np.array([bool(options.adaptive)], dtype=bool)
        hist["adaptive_window"] = np.array([int(options.adaptive_window)], dtype=int)
        hist["adaptive_time_step_multiplier"] = np.array([float(options.adaptive_time_step_multiplier)], dtype=float)
        hist["adaptive_growth_factor"] = np.array([float(getattr(options, "adaptive_growth_factor", np.nan))], dtype=float)
        hist["dt_init"] = np.array([float(options.dt_init)], dtype=float)
        hist["dt_max"] = np.array([float(options.dt_max)], dtype=float)
        rejected = np.asarray(hist.get("adaptive_rejected_attempts", []), dtype=float)
        hist["total_rejected_attempts"] = np.array([int(np.nansum(rejected))], dtype=int)
        retries = np.asarray(hist.get("adaptive_retries", []), dtype=float)
        hist["max_adaptive_retries_per_step"] = np.array([int(np.nanmax(retries)) if retries.size else 0], dtype=int)
        return PyTDGLLikeSolution(
            device=self.device,
            options=options,
            tdgl_data=result,
            history=hist,
            total_seconds=(end_time - start_time).total_seconds(),
        )
