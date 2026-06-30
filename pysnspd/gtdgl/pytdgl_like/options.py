"""pyTDGL-compatible solver options for the pySNSPD OE7 comparison backend.

This module intentionally mirrors the public names and constructor arguments of
``tdgl.solver.options`` so that later source-to-source comparisons against
pyTDGL are straightforward.  The implementation is local to pySNSPD and keeps
only the CPU/SuperLU path needed for stationary SNSPD tests.

pyTDGL attribution: function/class names and the option layout follow the MIT
licensed project ``loganbvh/py-tdgl``.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Union


class SolverOptionsError(ValueError):
    pass


class SparseSolver(Enum):
    """Supported sparse linear solvers."""

    SUPERLU: str = "superlu"
    UMFPACK: str = "umfpack"
    PARDISO: str = "pardiso"
    CUPY: str = "cupy"


@dataclass
class SolverOptions:
    """Options for the TDGL solver.

    The fields and defaults mirror pyTDGL's ``SolverOptions``.  Only the CPU
    SuperLU branch is implemented in this comparison backend; unsupported
    backends are rejected by ``validate`` instead of silently changing meaning.
    """

    solve_time: float
    skip_time: float = 0.0
    dt_init: float = 1e-6
    dt_max: float = 1e-1
    adaptive: bool = True
    adaptive_window: int = 10
    max_solve_retries: int = 10
    adaptive_time_step_multiplier: float = 0.25
    output_file: Union[str, None] = None
    terminal_psi: Union[float, complex, None] = 0.0
    gpu: bool = False
    sparse_solver: Union[SparseSolver, str] = SparseSolver.SUPERLU
    pause_on_interrupt: bool = True
    save_every: int = 100
    progress_interval: int = 0
    monitor: bool = False
    monitor_update_interval: float = 1.0
    field_units: str = "mT"
    current_units: str = "uA"
    include_screening: bool = False
    max_iterations_per_step: int = 1000
    screening_tolerance: float = 1e-3
    screening_step_size: float = 0.1
    screening_step_drag: float = 0.5

    def validate(self) -> None:
        if self.dt_init <= 0:
            raise SolverOptionsError("dt_init must be positive.")
        if self.dt_max <= 0:
            raise SolverOptionsError("dt_max must be positive.")
        if self.solve_time <= 0:
            raise SolverOptionsError("solve_time must be positive.")
        if self.dt_init > self.dt_max:
            raise SolverOptionsError("dt_init must be less than or equal to dt_max.")
        if self.terminal_psi is not None and not (0 <= abs(self.terminal_psi) <= 1):
            raise SolverOptionsError(
                "terminal_psi must be None or have absolute value in [0, 1]"
                f" (got {self.terminal_psi})."
            )
        if not (0 < self.adaptive_time_step_multiplier < 1):
            raise SolverOptionsError(
                "adaptive_time_step_multiplier must be in (0, 1)"
                f" (got {self.adaptive_time_step_multiplier})."
            )
        if not (0 < self.screening_step_drag <= 1):
            raise SolverOptionsError(
                "screening_step_drag must be in (0, 1]"
                f" (got {self.screening_step_drag})."
            )
        if self.screening_step_size <= 0:
            raise SolverOptionsError(
                "screening_step_size must be in > 0"
                f" (got {self.screening_step_size})."
            )
        if self.screening_tolerance <= 0:
            raise SolverOptionsError(
                "screening_tolerance must be in > 0"
                f" (got {self.screening_tolerance})."
            )

        solver = self.sparse_solver
        if isinstance(solver, str):
            try:
                solver = SparseSolver[solver.upper()]
            except KeyError as exc:
                valid_solvers = list(SparseSolver.__members__.keys())
                raise SolverOptionsError(
                    f"sparse solver must be one of {valid_solvers!r}, got {solver}."
                ) from exc
        self.sparse_solver = solver

        if self.gpu or self.sparse_solver is SparseSolver.CUPY:
            raise SolverOptionsError(
                "The pySNSPD pytdgl_like backend implements the CPU/SuperLU path only."
            )
        if self.sparse_solver is not SparseSolver.SUPERLU:
            raise SolverOptionsError(
                "Only SparseSolver.SUPERLU is enabled in this comparison backend."
            )
        if self.include_screening:
            raise SolverOptionsError(
                "include_screening=True is intentionally not implemented in the "
                "first pySNSPD pytdgl_like stationary backend."
            )
