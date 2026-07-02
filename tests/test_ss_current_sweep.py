import importlib.util
from pathlib import Path


def _load_pipeline_module():
    path = Path(__file__).resolve().parents[1] / "pipelines" / "02_ss_run_template.py"
    spec = importlib.util.spec_from_file_location("ss_run_template_for_test", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_current_sweep_run_name_is_stable_and_informative():
    module = _load_pipeline_module()

    name = module._current_sweep_run_name("base", 2.5, 30.5)

    assert name == "base_dI_plus2p5uA_I30p5uA"


def test_resolve_sweep_workers_uses_config_but_caps_to_extra_count():
    module = _load_pipeline_module()
    cfg = {"parallel": {"enabled": True, "workers": 16}}

    assert module._resolve_sweep_workers(cfg, None, n_extra=3) == 3
    assert module._resolve_sweep_workers(cfg, 2, n_extra=3) == 2
    assert module._resolve_sweep_workers(cfg, None, n_extra=0) == 1


def test_resolve_seed_current_clamps_overcritical_target():
    module = _load_pipeline_module()

    class Catalog:
        calibration_current_values_A = [0.0, 1.0e-6, 2.0e-6]

    seed_current_A, policy = module._resolve_seed_current_for_target(
        usadel_catalog=Catalog(),
        target_current_A=3.0e-6,
        overcritical_policy="clamp-to-ic",
        overcritical_seed_fraction=0.95,
    )

    assert seed_current_A == 0.95 * 2.0e-6
    assert policy["seed_is_overcritical_clamped"] is True
    assert policy["target_over_pre_usadel_Ic"] == 1.5
    assert policy["analytic_seed_over_pre_usadel_Ic"] == 0.95


def test_resolve_seed_current_uses_target_below_ic():
    module = _load_pipeline_module()

    class Catalog:
        calibration_current_values_A = [0.0, 1.0e-6, 2.0e-6]

    seed_current_A, policy = module._resolve_seed_current_for_target(
        usadel_catalog=Catalog(),
        target_current_A=1.5e-6,
        overcritical_policy="clamp-to-ic",
        overcritical_seed_fraction=0.95,
    )

    assert seed_current_A == 1.5e-6
    assert policy["seed_is_overcritical_clamped"] is False
    assert policy["analytic_seed_over_pre_usadel_Ic"] == 0.75


def test_resolve_seed_current_error_policy_rejects_overcritical_target():
    module = _load_pipeline_module()

    class Catalog:
        calibration_current_values_A = [0.0, 1.0e-6, 2.0e-6]

    try:
        module._resolve_seed_current_for_target(
            usadel_catalog=Catalog(),
            target_current_A=3.0e-6,
            overcritical_policy="error",
            overcritical_seed_fraction=0.95,
        )
    except ValueError as exc:
        assert "exceeds PRE Usadel Ic" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("Expected overcritical target to be rejected with policy='error'.")
