"""Minimal readout-circuit model for pySNSPD photon transients.

The circuit is deliberately small and SI-only.  It is used by
``pipelines/03_photon_run_template.py`` and does not alter the existing SS
pipeline.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Mapping

import numpy as np


@dataclass(frozen=True)
class CircuitParams:
    """Lumped bias/readout circuit parameters in SI units."""

    R_load_ohm: float = 50.0
    R_bias_ohm: float = 1.0e4
    L_bias_H: float = 1.0e-6
    L_k_H: float = 10.0e-9
    C_couple_F: float = 100.0e-12
    V_bias_V: float | None = None

    def validated(self) -> "CircuitParams":
        for name in ("R_load_ohm", "R_bias_ohm", "L_bias_H", "L_k_H", "C_couple_F"):
            value = float(getattr(self, name))
            if not np.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be positive and finite; got {value!r}.")
        if self.V_bias_V is not None and not np.isfinite(float(self.V_bias_V)):
            raise ValueError("V_bias_V must be finite when provided.")
        return self


@dataclass(frozen=True)
class CircuitState:
    """State of the lumped circuit.

    ``I_s_A`` is the current imposed on the SNSPD terminals during the next
    mesoscopic chunk.
    """

    I_b_A: float
    I_s_A: float
    v_c_V: float

    def as_dict(self) -> dict[str, float]:
        return asdict(self)


@dataclass(frozen=True)
class CircuitRuntimeConfig:
    """Runtime activation and stationarity policy for an SS circuit."""

    start_time_s: float = 5.0e-12
    hold_time_s: float = 5.0e-12
    current_relative_tolerance: float = 1.0e-2
    current_absolute_tolerance_A: float = 0.05e-6
    voltage_relative_tolerance: float = 1.0e-2
    voltage_absolute_tolerance_V: float = 10.0e-6

    def validated(self) -> "CircuitRuntimeConfig":
        if not np.isfinite(self.start_time_s) or self.start_time_s < 0.0:
            raise ValueError("Circuit start time must be finite and non-negative.")
        if not np.isfinite(self.hold_time_s) or self.hold_time_s <= 0.0:
            raise ValueError("Circuit hold time must be positive and finite.")
        for name in (
            "current_relative_tolerance",
            "current_absolute_tolerance_A",
            "voltage_relative_tolerance",
            "voltage_absolute_tolerance_V",
        ):
            value = float(getattr(self, name))
            if not np.isfinite(value) or value < 0.0:
                raise ValueError(f"{name} must be finite and non-negative.")
        return self


class CircuitRuntimeController:
    """Evolve the lumped circuit inside a monolithic SS solver."""

    def __init__(
        self,
        *,
        I_ss_A: float,
        V_tdgl_ss_V: float,
        params: CircuitParams,
        config: CircuitRuntimeConfig,
    ) -> None:
        self.config = config.validated()
        self.state, self.params = initialize_circuit_from_ss(
            I_ss_A=float(I_ss_A),
            V_tdgl_ss_V=float(V_tdgl_ss_V),
            params=params,
        )
        self.initial_state = self.state
        self.last_V_tdgl_V = float(V_tdgl_ss_V)
        self.last_diagnostics = self._diagnostics(
            time_s=0.0,
            active=False,
            V_tdgl_V=float(V_tdgl_ss_V),
        )

    def terminal_current_A(self, time_s: float) -> float:
        if float(time_s) < float(self.config.start_time_s):
            return float(self.initial_state.I_s_A)
        return float(self.state.I_s_A)

    def step(
        self,
        *,
        time_s: float,
        dt_s: float,
        V_tdgl_V: float,
    ) -> dict[str, float | bool]:
        end_s = float(time_s)
        start_s = end_s - float(dt_s)
        active_dt = max(0.0, end_s - max(start_s, float(self.config.start_time_s)))
        active = active_dt > 0.0
        if active:
            self.state = step_circuit_rk2(
                self.state,
                V_tdgl_V=float(V_tdgl_V),
                dt_s=active_dt,
                params=self.params,
            )
        self.last_V_tdgl_V = float(V_tdgl_V)
        self.last_diagnostics = self._diagnostics(
            time_s=end_s,
            active=active,
            V_tdgl_V=float(V_tdgl_V),
        )
        return dict(self.last_diagnostics)

    def snapshot_payload(self) -> dict[str, np.ndarray]:
        obs = circuit_observables(
            self.state,
            params=self.params,
            V_tdgl_V=self.last_V_tdgl_V,
        )
        return {
            "circuit_I_b_A_snapshot": np.asarray(float(obs["I_b_A"])),
            "circuit_I_s_A_snapshot": np.asarray(float(obs["I_s_A"])),
            "circuit_I_rf_A_snapshot": np.asarray(float(obs["I_rf_A"])),
            "circuit_V_out_V_snapshot": np.asarray(float(obs["V_out_V"])),
            "circuit_v_c_V_snapshot": np.asarray(float(obs["v_c_V"])),
            "circuit_V_tdgl_center_V_snapshot": np.asarray(
                float(obs["V_tdgl_center_V"])
            ),
        }

    def _diagnostics(
        self,
        *,
        time_s: float,
        active: bool,
        V_tdgl_V: float,
    ) -> dict[str, float | bool]:
        obs = circuit_observables(
            self.state,
            params=self.params,
            V_tdgl_V=V_tdgl_V,
        )
        rhs = circuit_rhs(self.state, V_tdgl_V=V_tdgl_V, params=self.params)
        return {
            "circuit_time_s": float(time_s),
            "circuit_enabled": True,
            "circuit_active": bool(
                active or time_s >= float(self.config.start_time_s)
            ),
            **{f"circuit_{key}": float(value) for key, value in obs.items()},
            "circuit_dI_b_A_s": float(rhs.I_b_A),
            "circuit_dI_s_A_s": float(rhs.I_s_A),
            "circuit_dv_c_V_s": float(rhs.v_c_V),
        }


def circuit_stationarity_diagnostics(
    history: Mapping[str, Any],
    *,
    config: CircuitRuntimeConfig,
) -> dict[str, Any]:
    """Classify a stationary circuit tail using values and circuit RHS."""

    cfg = config.validated()
    time = np.asarray(history.get("circuit_time_s", []), dtype=float).reshape(-1)
    active = np.asarray(history.get("circuit_active", []), dtype=bool).reshape(-1)
    if time.size == 0 or active.size == 0:
        return {
            "passes": False,
            "reason": "circuit runtime history is unavailable",
            "sufficient_history": False,
        }
    n = min(time.size, active.size)
    time = time[:n]
    active = active[:n]
    end_s = float(time[-1])
    time_slack_s = max(1.0e-24, float(cfg.hold_time_s) * 1.0e-12)
    tail = active & (
        time >= end_s - float(cfg.hold_time_s) - time_slack_s
    )
    duration_s = (
        float(time[tail][-1] - time[tail][0]) if np.count_nonzero(tail) >= 2 else 0.0
    )
    sufficient = bool(duration_s >= 0.999 * float(cfg.hold_time_s))

    initial_I = abs(
        float(np.asarray(history.get("circuit_I_s_A", [0.0]), dtype=float).reshape(-1)[0])
    )
    current_tol = max(
        float(cfg.current_absolute_tolerance_A),
        float(cfg.current_relative_tolerance) * initial_I,
    )
    voltage_scale = abs(
        float(
            np.asarray(
                history.get("circuit_V_tdgl_center_V", [0.0]),
                dtype=float,
            ).reshape(-1)[0]
        )
    )
    voltage_tol = max(
        float(cfg.voltage_absolute_tolerance_V),
        float(cfg.voltage_relative_tolerance) * voltage_scale,
    )
    component_tolerances = {
        "circuit_I_b_A": current_tol,
        "circuit_I_s_A": current_tol,
        "circuit_I_rf_A": current_tol,
        "circuit_V_out_V": voltage_tol,
        "circuit_v_c_V": voltage_tol,
        "circuit_V_tdgl_center_V": voltage_tol,
    }
    spans: dict[str, float] = {}
    values_pass = sufficient
    for key, tolerance in component_tolerances.items():
        values = np.asarray(history.get(key, []), dtype=float).reshape(-1)[:n]
        if values.size != n or not np.any(tail):
            span = float("inf")
        else:
            span = float(np.nanmax(values[tail]) - np.nanmin(values[tail]))
        spans[key] = span
        values_pass = bool(values_pass and np.isfinite(span) and span <= tolerance)

    rhs_tolerances = {
        "circuit_dI_b_A_s": current_tol / float(cfg.hold_time_s),
        "circuit_dI_s_A_s": current_tol / float(cfg.hold_time_s),
        "circuit_dv_c_V_s": voltage_tol / float(cfg.hold_time_s),
    }
    final_rhs: dict[str, float] = {}
    rhs_pass = True
    for key, tolerance in rhs_tolerances.items():
        values = np.asarray(history.get(key, []), dtype=float).reshape(-1)
        value = abs(float(values[min(values.size, n) - 1])) if values.size else float("inf")
        final_rhs[key] = value
        rhs_pass = bool(rhs_pass and np.isfinite(value) and value <= tolerance)

    passes = bool(sufficient and values_pass and rhs_pass)
    return {
        "diagnostic": "lumped_circuit_stationarity_v1",
        "passes": passes,
        "reason": (
            "circuit values and RHS are stationary over the hold window"
            if passes
            else "circuit stationarity tolerances are not all satisfied"
        ),
        "sufficient_history": sufficient,
        "tail_duration_ps": duration_s / 1.0e-12,
        "hold_time_ps": float(cfg.hold_time_s / 1.0e-12),
        "component_spans": spans,
        "component_tolerances": component_tolerances,
        "final_rhs_absolute": final_rhs,
        "rhs_tolerances": rhs_tolerances,
    }


def initialize_circuit_from_ss(
    *,
    I_ss_A: float,
    V_tdgl_ss_V: float,
    params: CircuitParams,
) -> tuple[CircuitState, CircuitParams]:
    """Return an exactly stationary circuit state for the SS solution.

    The stationary condition is

    ``I_b = I_s = I_ss``, ``v_c = V_tdgl_ss`` and
    ``V_bias = R_bias I_ss + V_tdgl_ss``.
    """

    p = params.validated()
    I_ss_A = float(I_ss_A)
    V_tdgl_ss_V = float(V_tdgl_ss_V)
    if not np.isfinite(I_ss_A):
        raise ValueError("I_ss_A must be finite.")
    if not np.isfinite(V_tdgl_ss_V):
        raise ValueError("V_tdgl_ss_V must be finite.")

    if p.V_bias_V is None:
        p = CircuitParams(
            R_load_ohm=float(p.R_load_ohm),
            R_bias_ohm=float(p.R_bias_ohm),
            L_bias_H=float(p.L_bias_H),
            L_k_H=float(p.L_k_H),
            C_couple_F=float(p.C_couple_F),
            V_bias_V=float(p.R_bias_ohm) * I_ss_A + V_tdgl_ss_V,
        )

    return CircuitState(I_b_A=I_ss_A, I_s_A=I_ss_A, v_c_V=V_tdgl_ss_V), p.validated()


def circuit_rhs(
    state: CircuitState,
    *,
    V_tdgl_V: float,
    params: CircuitParams,
) -> CircuitState:
    p = params.validated()
    I_b = float(state.I_b_A)
    I_s = float(state.I_s_A)
    v_c = float(state.v_c_V)
    V_tdgl = float(V_tdgl_V)
    V_bias = float(p.V_bias_V if p.V_bias_V is not None else 0.0)

    I_rf = I_b - I_s
    V_load = float(p.R_load_ohm) * I_rf

    dI_b = (V_bias - float(p.R_bias_ohm) * I_b - v_c - V_load) / float(p.L_bias_H)
    dI_s = (v_c + V_load - V_tdgl) / float(p.L_k_H)
    dv_c = I_rf / float(p.C_couple_F)

    return CircuitState(I_b_A=dI_b, I_s_A=dI_s, v_c_V=dv_c)


def step_circuit_rk2(
    state: CircuitState,
    *,
    V_tdgl_V: float,
    dt_s: float,
    params: CircuitParams,
) -> CircuitState:
    """Second-order explicit circuit step.

    The mesoscopic solver supplies ``V_tdgl_V`` at the end of the chunk.  This
    RK2 update is cheap and stable for the small chunks used in pipeline 03.
    """

    dt = float(dt_s)
    if not np.isfinite(dt) or dt <= 0.0:
        raise ValueError("dt_s must be positive and finite.")

    k1 = circuit_rhs(state, V_tdgl_V=V_tdgl_V, params=params)
    mid = CircuitState(
        I_b_A=float(state.I_b_A) + 0.5 * dt * float(k1.I_b_A),
        I_s_A=float(state.I_s_A) + 0.5 * dt * float(k1.I_s_A),
        v_c_V=float(state.v_c_V) + 0.5 * dt * float(k1.v_c_V),
    )
    k2 = circuit_rhs(mid, V_tdgl_V=V_tdgl_V, params=params)

    return CircuitState(
        I_b_A=float(state.I_b_A) + dt * float(k2.I_b_A),
        I_s_A=float(state.I_s_A) + dt * float(k2.I_s_A),
        v_c_V=float(state.v_c_V) + dt * float(k2.v_c_V),
    )


def central_tdgl_voltage_V(
    *,
    nodes_m: np.ndarray,
    phi_V: np.ndarray,
    center_width_m: float = 100.0e-9,
    probe_band_m: float | None = None,
) -> float:
    """Voltage seen by the circuit: right-minus-left across the center window.

    This is intentionally *not* the full terminal voltage.  It uses the two
    x-probe cuts at the ends of the central window, averaged over all y nodes
    close to each cut.
    """

    nodes = np.asarray(nodes_m, dtype=float)
    if nodes.ndim != 2 or nodes.shape[1] < 2:
        raise ValueError("nodes_m must have shape (n_nodes, >=2).")
    phi = np.asarray(phi_V, dtype=float).reshape(-1)
    if phi.size != nodes.shape[0]:
        raise ValueError(f"phi_V has length {phi.size}, expected {nodes.shape[0]}.")

    x = nodes[:, 0]
    x_mid = 0.5 * (float(np.nanmin(x)) + float(np.nanmax(x)))
    half = 0.5 * float(center_width_m)
    x_left = x_mid - half
    x_right = x_mid + half

    if probe_band_m is None:
        unique_x = np.unique(np.round(x[np.isfinite(x)], decimals=15))
        if unique_x.size > 1:
            dx = float(np.nanmedian(np.diff(np.sort(unique_x))))
            probe_band_m = max(0.55 * abs(dx), 1.0e-12)
        else:
            probe_band_m = max(0.02 * float(center_width_m), 1.0e-12)

    left_mask = np.abs(x - x_left) <= float(probe_band_m)
    right_mask = np.abs(x - x_right) <= float(probe_band_m)

    if not np.any(left_mask):
        left_mask = np.abs(x - x_left) == np.nanmin(np.abs(x - x_left))
    if not np.any(right_mask):
        right_mask = np.abs(x - x_right) == np.nanmin(np.abs(x - x_right))

    return float(np.nanmean(phi[right_mask]) - np.nanmean(phi[left_mask]))
    #return float(np.nanmean(phi[left_mask]) - np.nanmean(phi[right_mask]))

def circuit_observables(
    state: CircuitState,
    *,
    params: CircuitParams,
    V_tdgl_V: float,
) -> dict[str, float]:
    p = params.validated()
    I_rf = float(state.I_b_A) - float(state.I_s_A)
    return {
        "I_b_A": float(state.I_b_A),
        "I_s_A": float(state.I_s_A),
        "I_rf_A": float(I_rf),
        "v_c_V": float(state.v_c_V),
        "V_tdgl_center_V": float(V_tdgl_V),
        "V_out_V": float(p.R_load_ohm) * I_rf,
        "V_bias_V": float(p.V_bias_V if p.V_bias_V is not None else np.nan),
    }


__all__ = [
    "CircuitParams",
    "CircuitRuntimeConfig",
    "CircuitRuntimeController",
    "CircuitState",
    "central_tdgl_voltage_V",
    "circuit_observables",
    "circuit_stationarity_diagnostics",
    "initialize_circuit_from_ss",
    "step_circuit_rk2",
]
