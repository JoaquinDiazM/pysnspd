"""Configuration tests for the coupled photon transient."""

from __future__ import annotations

import pytest

from pysnspd.gtdgl.photon_transient import CoupledTransientConfig


def _config(**kwargs) -> CoupledTransientConfig:
    values = {
        "total_time_s": 1.0e-12,
        "mesoscopic_dt_s": 1.0e-15,
        "chunk_time_s": 1.0e-13,
    }
    values.update(kwargs)
    return CoupledTransientConfig(**values)


def test_photon_transient_can_freeze_thermal_dynamics():
    cfg = _config(thermal_enabled=False).validated()
    assert cfg.thermal_enabled is False


@pytest.mark.parametrize(
    ("key", "value"),
    [
        ("allmaras_phase_direct_amplitude_fraction", 0.0),
        ("allmaras_phase_convergence_tol", 0.0),
        ("allmaras_phase_convergence_max_iterations", 0),
    ],
)
def test_photon_transient_rejects_invalid_phase_continuation_controls(key, value):
    with pytest.raises(ValueError):
        _config(**{key: value}).validated()
