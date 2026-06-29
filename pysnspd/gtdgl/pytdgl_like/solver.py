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
from dataclasses import dataclass, field
from datetime import datetime
from typing import Callable, Dict, List, NamedTuple, Optional, Sequence, Tuple, Union

import numpy as np
import scipy.sparse as sp

from .device import PySNSPDTDGLDevice as Device, TerminalInfo
from .operators import MeshOperators
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
    ):
        self.device = device
        self.options = options
        self.options.validate()
        self.terminal_currents = terminal_currents
        self.seed_solution = seed_solution
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
        self.terminal_current_densities = {name: 0.0 for name in self.terminal_names}

        terminal_psi = options.terminal_psi
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
        if terminal_psi is not None and normal_boundary_index.size:
            psi_init[normal_boundary_index] = terminal_psi
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

    def update_mu_boundary(self, time: float) -> None:
        """Computes terminal current density and updates scalar-potential BCs."""

        currents = self.current_func(time)
        terminal_current_densities = self.terminal_current_densities
        for terminal in self.terminal_info:
            current_density = (-1 / terminal.length) * sum(
                currents.get(name, 0.0)
                for name in self.terminal_names
                if name != terminal.name
            )
            if current_density != terminal_current_densities[terminal.name]:
                terminal_current_densities[terminal.name] = current_density
                # ``mu_boundary`` is indexed by the compact boundary-edge vector.
                bmap = self.device.mesh.edge_mesh.boundary_edge_indices
                pos = {int(edge): k for k, edge in enumerate(np.asarray(bmap, dtype=int))}
                for edge in terminal.boundary_edge_indices:
                    k = pos.get(int(edge))
                    if k is not None:
                        self.mu_boundary[k] = current_density

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
    ) -> Union[Tuple[np.ndarray, np.ndarray], None]:
        """Solves for psi^{n+1} and |psi^{n+1}|^2.

        This is the pyTDGL local algebraic solve with the same public arguments.
        In this comparison backend, ``epsilon``, ``gamma`` and ``u`` are supplied
        by the pySNSPD adapter and encode the modified KWT/Allmaras coefficients.
        """

        xp = np
        U = xp.exp(-1j * mu * dt)
        z = U * gamma**2 / 2 * psi
        with np.errstate(all="raise"):
            try:
                w = z * abs_sq_psi + U * (
                    psi
                    + (dt / u)
                    * xp.sqrt(1 + gamma**2 * abs_sq_psi)
                    * ((epsilon - abs_sq_psi) * psi + psi_laplacian @ psi)
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
        kwargs = dict(
            psi=psi,
            abs_sq_psi=abs_sq_psi,
            mu=mu,
            epsilon=epsilon,
            gamma=self.device.layer.gamma,
            u=self.device.layer.u,
            dt=dt,
            psi_laplacian=self.operators.psi_laplacian,
        )
        result = self.solve_for_psi_squared(**kwargs)
        for retries in itertools.count():
            if result is not None:
                break
            if not options.adaptive or retries > options.max_solve_retries:
                raise RuntimeError(
                    f"Solver failed to converge in {options.max_solve_retries} retries "
                    f"at step {step} with dt = {dt:.2e}. Try using a smaller dt_init."
                )
            kwargs["dt"] = dt = dt * options.adaptive_time_step_multiplier
            result = self.solve_for_psi_squared(**kwargs)
        psi, new_sq_psi = result
        return psi, new_sq_psi, dt

    def solve_for_observables(
        self,
        psi: np.ndarray,
        dA_dt: Union[float, np.ndarray],
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Solves for mu, supercurrent and normal current."""

        operators = self.operators
        supercurrent = operators.get_supercurrent(psi)
        rhs = (operators.divergence @ (supercurrent - dA_dt)) - (
            operators.mu_boundary_laplacian @ self.mu_boundary
        )
        mu = operators.mu_laplacian_lu(rhs)
        normal_current = -(operators.mu_gradient @ mu) - dA_dt
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

        if options.adaptive:
            self.d_psi_sq_vals.append(float(np.absolute(abs_sq_psi - old_sq_psi).max()))
            window = max(1, int(options.adaptive_window))
            if step > window:
                new_dt = options.dt_init / max(1e-10, np.mean(self.d_psi_sq_vals[-window:]))
                self.tentative_dt = float(np.clip(0.5 * (new_dt + dt), 0, self.dt_max))

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

        running_state = RunningState()
        state: dict[str, numbers.Real] = {"step": 0, "time": 0.0}
        dt = float(options.dt_init)
        result = SolverResult(
            dt,
            parameters["psi"],
            parameters["mu"],
            parameters["supercurrent"],
            parameters["normal_current"],
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
            running_state.append("normal_current_snapshot", frame.normal_current)

        append_snapshot(0.0, result)
        next_snapshot = 1

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

            if int(state["step"]) > 10_000_000:
                raise RuntimeError("pytdgl_like solve exceeded 10,000,000 steps.")

        if next_snapshot < snapshot_times.size:
            append_snapshot(float(state["time"]), result)

        end_time = datetime.now()
        hist = {key: np.asarray(vals) for key, vals in running_state.data.items()}
        hist["final_step"] = np.array([int(state["step"])], dtype=int)
        hist["final_time"] = np.array([float(state["time"])], dtype=float)
        return PyTDGLLikeSolution(
            device=self.device,
            options=options,
            tdgl_data=result,
            history=hist,
            total_seconds=(end_time - start_time).total_seconds(),
        )
