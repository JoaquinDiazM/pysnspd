"""Synthetic staged-campaign admission, affinity, ownership and resume tests."""
import importlib.util
import io
import json
from pathlib import Path
import threading
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("thermal_candidate_campaign_test", ROOT/"pipelines/10_thermal_stiffness_basecases.py")
campaign = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(campaign)


def resources(memory=128*1024**3, quota=None):
    return dict(logical_cpus=list(range(32)), physical_core_groups=[[i, i+16] for i in range(16)],
        available_memory_bytes=memory, cpu_quota_units=quota)


def test_twelve_cases_and_commands_have_exact_physical_identities():
    cases = campaign.catalog()
    assert len(cases) == 12
    assert sum(case["stage"] == "anchors" for case in cases) == 2
    assert sum(case["stage"] == "controls" for case in cases) == 10
    by_pair = {}
    for case in cases:
        by_pair.setdefault(case["pair"], []).append(case)
        assert case["parameters"]["length_nm"] == 320.
        assert case["parameters"]["spacing_nm"] == 5.
        assert case["parameters"]["current_uA"] == 21.5
        assert case["parameters"]["tau_escape_ps"] == 20.
        command = campaign.command(case, Path("/home/jdiaz/scratch/fixture"), threads=4 if case["stage"] == "anchors" else 1, cpu_start=3)
        args = campaign.RUNNER.parser().parse_args(command[3:])
        identity = campaign.RUNNER.BASE.parameter_identity(args)
        expected = dict(case["parameters"], threads=args.threads)
        assert identity == expected
    assert len(by_pair) == 6
    for pair in by_pair.values():
        first, second = (dict(case["parameters"]) for case in pair)
        assert {first.pop("wavelength_nm"), second.pop("wavelength_nm")} == {775., 1550.}
        assert first == second


def test_global_allocation_reserves_complete_cores_and_distinct_worker_groups():
    plan = campaign.allocation("anchors", resources=resources(), active=[])
    assert plan["worker_physical_cores"] == 12
    assert [slot["threads"] for slot in plan["slots"]] == [6, 6]
    sets = [{cpu for group in slot["core_groups"] for cpu in group} for slot in plan["slots"]]
    assert not sets[0]&sets[1]
    assert not (sets[0]|sets[1])&set(plan["coordinator_core_group"])
    assert plan["budget"]["physical_cores_reserved"] == 2
    # Foreign activity on a normally reserved core still consumes global budget.
    outside = [dict(pid=91, affinity_cpus=[14, 30])]
    shared = campaign.allocation("anchors", resources=resources(), active=outside)
    assert shared["worker_physical_cores"]+1+1 <= 14


def test_active_baseline_jobs_preserved_and_low_utilization_not_forced():
    active = [dict(pid=100+i, affinity_cpus=[i, i+16]) for i in (2, 4, 6, 8, 10, 13)]
    anchors = campaign.allocation("anchors", resources=resources(), active=active)
    assert [slot["threads"] for slot in anchors["slots"]] == [2, 2]
    assert anchors["worker_physical_cores"] == 4
    allocated = campaign.allocation("controls", resources=resources(), active=active)
    occupied = {cpu for item in active for cpu in item["affinity_cpus"]}
    assert allocated["worker_physical_cores"] == 7
    assert not set(allocated["coordinator_core_group"])&occupied
    for slot in allocated["slots"]:
        assert not set(slot["core_groups"][0])&occupied
    small = campaign.allocation("controls", resources=resources(memory=16*1024**3), active=[])
    assert len(small["slots"]) == 2
    quota = campaign.allocation("controls", resources=resources(quota=16.), active=[])
    assert quota["worker_physical_cores"]+1 <= 7


def test_process_ownership_does_not_count_timeout_taskset_or_shell_twice():
    python = [b"/home/jdiaz/.conda/envs/snspd/bin/python", b"-u", b"/home/jdiaz/pysnspd/pipelines/09_thermal_stiffness_photon_run.py"]
    assert campaign.known_runner(python) == "09_thermal_stiffness_photon_run.py"
    for wrapper in (b"timeout", b"taskset", b"bash"):
        assert campaign.known_runner([wrapper, b"240s", *python]) is None
    assert campaign.known_runner([b"python", b"-m", b"pytest"]) is None


def fake_terminal(tmp_path):
    case = campaign.catalog()[0]
    threads = 4
    mode = case["parameters"]["stiffness_temperature"]
    definition = dict(source_pins_sha256={mode: {"f.py": "identity"}})
    run = tmp_path/"runs"/case["name"]
    (run/"dc").mkdir(parents=True)
    files = ["snapshots.npz", "checkpoint_latest.npz", "dc/effective_dc.npz", "material.npz", "material.json", "source_initial_state.npz"]
    for name in files:
        (run/name).write_bytes(b"authenticated fixture")
    np.savez(run/"trajectory.npz", time_s=np.array([0., 200e-12]), energy_total_J=np.array([2., 2.]), external_energy_integral_J=np.array([0., 0.]))
    calibration = {"diffusion_m2_s": 1.4103757147e-4}
    identity = dict(case["parameters"], threads=threads)
    physics = dict(calibration=calibration, parameter_identity=identity, source_pins_sha256=definition["source_pins_sha256"][mode])
    campaign.write_json(run/"calibration.json", calibration)
    campaign.write_json(run/"candidate_physics_identity.json", physics)
    manifest = dict(parameter_identity=identity, source_pins_sha256=physics["source_pins_sha256"],
        run_complete=True, status="TIME_CENSORED_NO_TRIGGER", model_version="1.1.0-candidate",
        closure=dict(stiffness_temperature=mode), physicalparams=dict(diffusion_calibration=calibration),
        observation=dict(threshold_voltage_V=.1e-3, confirmation_time_s_required=.5e-12,
            extra_time_s=5e-12, no_trigger_ceiling_s=200e-12, full_reset_requested=False,
            trigger_confirmed=False, first_crossing_time_s=None),
        output_sha256={name: campaign.sha(run/name) for name in files+["trajectory.npz"]},
        source=dict(retained_energy_J=1., relative_total_energy_defect=0., maximum_nodal_energy_defect_J=0.,
            positive_capacity_checked_throughout_preparation=True, minimum_preparation_effective_capacity_ratio=.9),
        thermal_candidate_diagnostics=dict(min_successful_rhs_evaluations=.8),
        local_energy_control=dict(max_guarded_accepted_local_energy_fraction=1e-7),
        numerical=dict(local_energy_tolerance=2.5e-6))
    campaign.write_json(run/"manifest.json", manifest)
    return case, definition, run, manifest, threads


def test_terminal_gate_checks_hash_capacity_source_and_tolerance(tmp_path):
    case, definition, run, manifest, threads = fake_terminal(tmp_path)
    result = campaign.audit_case(case, tmp_path, definition, expected_threads=threads)
    assert Path(result["path"]).is_file()
    manifest["thermal_candidate_diagnostics"]["min_successful_rhs_evaluations"] = -1.
    campaign.write_json(run/"manifest.json", manifest)
    with pytest.raises(ValueError, match="capacity"):
        campaign.audit_case(case, tmp_path, definition, expected_threads=threads)
    manifest["thermal_candidate_diagnostics"]["min_successful_rhs_evaluations"] = .8
    manifest["local_energy_control"]["max_guarded_accepted_local_energy_fraction"] = 1e-3
    campaign.write_json(run/"manifest.json", manifest)
    with pytest.raises(ValueError, match="local-energy"):
        campaign.audit_case(case, tmp_path, definition, expected_threads=threads)
    manifest["local_energy_control"]["max_guarded_accepted_local_energy_fraction"] = 1e-7
    campaign.write_json(run/"manifest.json", manifest)
    (run/"trajectory.npz").write_bytes(b"tampered")
    with pytest.raises(ValueError, match="hash mismatch"):
        campaign.audit_case(case, tmp_path, definition, expected_threads=threads)


def test_controls_refuse_missing_anchor_receipts(tmp_path):
    with pytest.raises(ValueError, match="both completed"):
        campaign.require_anchors(tmp_path, dict(cases=campaign.catalog()), dict(cases={}))


def test_foreground_workers_use_distinct_slots_and_completed_resume_is_exact(tmp_path, monkeypatch):
    # Fake child processes only: no physics, subprocess, CPU pin or long run.
    monkeypatch.setattr(campaign, "SCRATCH", tmp_path)
    monkeypatch.setattr(campaign, "linux_resources", resources)
    monkeypatch.setattr(campaign, "active_known_processes", lambda exclude=(): [])
    monkeypatch.setattr(campaign.os, "sched_getaffinity", lambda pid: set(range(32)))
    monkeypatch.setattr(campaign.os, "sched_setaffinity", lambda pid, cpus: None)
    monkeypatch.setattr(campaign, "audit_case", lambda case, *args, **kwargs:
                        dict(path=str(tmp_path/"fake_audit.json"), sha256="fixture", manifest_sha256="fixture"))
    mutex = threading.Lock()
    barrier = threading.Barrier(2)
    commands = []
    active_sets = []
    class FakeProcess:
        def __init__(self, cmd, **kwargs):
            with mutex:
                self.pid = 99000+len(commands)
                commands.append(cmd)
            self.cmd = cmd
            self.stdout = self.lines()
        def lines(self):
            start = int(self.cmd[self.cmd.index("--cpu-start")+1])
            count = int(self.cmd[self.cmd.index("--threads")+1])
            mine = set(range(start, start+count))
            with mutex:
                assert all(not mine&other for other in active_sets)
                active_sets.append(mine)
            barrier.wait(timeout=5)
            yield '{"status":"synthetic fixture"}\n'
            with mutex:
                active_sets.remove(mine)
        def wait(self):
            return 0
        def send_signal(self, sig):
            raise AssertionError("No synthetic failure expected")
    monkeypatch.setattr(campaign.subprocess, "Popen", FakeProcess)
    args = ["--campaign-dir", str(tmp_path), "--stage", "anchors"]
    assert campaign.main(args) == 0
    assert len(commands) == 2
    assert campaign.main(args+["--resume"]) == 0
    assert len(commands) == 2
    receipt = json.loads((tmp_path/"campaign_execution.json").read_text())
    assert all(item["status"] == "COMPLETE_AUDITED" for item in receipt["cases"].values())
    assert all(item["threads"] == 6 for item in receipt["cases"].values())