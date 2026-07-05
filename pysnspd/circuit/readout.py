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

    def as_dict(self) -> dict[str, float | None]:
        return asdict(self)


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

    #return float(np.nanmean(phi[right_mask]) - np.nanmean(phi[left_mask]))
    return float(np.nanmean(phi[left_mask]) - np.nanmean(phi[right_mask]))

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
    "CircuitState",
    "central_tdgl_voltage_V",
    "circuit_observables",
    "initialize_circuit_from_ss",
    "step_circuit_rk2",
]
