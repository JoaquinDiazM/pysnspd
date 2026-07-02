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
