"""SS circuit activation and stationarity policy."""

from __future__ import annotations

import numpy as np

from pysnspd.circuit.readout import (
    CircuitParams,
    CircuitRuntimeConfig,
    CircuitRuntimeController,
    circuit_stationarity_diagnostics,
)


def test_runtime_circuit_stays_fixed_before_shared_start():
    controller = CircuitRuntimeController(
        I_ss_A=30.0e-6,
        V_tdgl_ss_V=2.0e-6,
        params=CircuitParams(),
        config=CircuitRuntimeConfig(start_time_s=5.0e-12),
    )

    diagnostics = controller.step(
        time_s=4.0e-12,
        dt_s=1.0e-12,
        V_tdgl_V=100.0e-6,
    )

    assert diagnostics["circuit_active"] is False
    assert controller.state == controller.initial_state


def test_runtime_circuit_uses_only_post_start_fraction_of_crossing_step():
    controller = CircuitRuntimeController(
        I_ss_A=30.0e-6,
        V_tdgl_ss_V=0.0,
        params=CircuitParams(),
        config=CircuitRuntimeConfig(start_time_s=5.0e-12),
    )

    controller.step(
        time_s=6.0e-12,
        dt_s=2.0e-12,
        V_tdgl_V=1.0e-3,
    )

    assert controller.state != controller.initial_state
    assert controller.last_diagnostics["circuit_active"] is True


def test_stationary_circuit_tail_passes_value_and_rhs_criteria():
    config = CircuitRuntimeConfig(start_time_s=5.0e-12, hold_time_s=5.0e-12)
    controller = CircuitRuntimeController(
        I_ss_A=30.0e-6,
        V_tdgl_ss_V=2.0e-6,
        params=CircuitParams(),
        config=config,
    )
    rows = [
        controller.step(
            time_s=float(time_ps) * 1.0e-12,
            dt_s=1.0e-12,
            V_tdgl_V=2.0e-6,
        )
        for time_ps in range(1, 13)
    ]
    history = {
        key: np.asarray([row[key] for row in rows])
        for key in rows[0]
    }

    result = circuit_stationarity_diagnostics(history, config=config)

    assert result["sufficient_history"] is True
    assert result["passes"] is True
