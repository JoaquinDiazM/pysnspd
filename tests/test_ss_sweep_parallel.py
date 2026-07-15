from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


def _pipeline_module():
    path = Path(__file__).resolve().parents[1] / "pipelines" / "02_ss_run_template.py"
    spec = importlib.util.spec_from_file_location("ss_run_template_for_test", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_unsigned_extra_current_arguments_are_absolute_currents():
    pipeline = _pipeline_module()
    offsets = pipeline._resolve_extra_current_specs(
        base_current_uA=20.0,
        specs=["35", "38", "39"],
    )
    assert offsets == [15.0, 18.0, 19.0]


def test_explicit_plus_arguments_are_relative_offsets():
    pipeline = _pipeline_module()
    offsets = pipeline._resolve_extra_current_specs(
        base_current_uA=20.0,
        specs=["+12", "+45", "+90"],
    )
    assert offsets == [12.0, 45.0, 90.0]


def test_sweep_offsets_reject_duplicate_parallel_output_directory():
    pipeline = _pipeline_module()
    with pytest.raises(ValueError, match="Duplicate sweep current 39"):
        pipeline._resolve_extra_current_specs(
            base_current_uA=20.0,
            specs=["39", "39"],
        )
