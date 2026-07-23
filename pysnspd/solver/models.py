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

from pysnspd.mesh.device import PySNSPDTDGLDevice as Device, TerminalInfo
from pysnspd.gtdgl.tdgl_operators import MeshOperators
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
