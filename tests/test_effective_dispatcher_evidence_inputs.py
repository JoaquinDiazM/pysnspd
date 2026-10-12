"""Repository baseline evidence is opt-in for the voltage decision only; no evolution."""
import hashlib,importlib.util,json,os,sys
from pathlib import Path
from types import SimpleNamespace
import pytest

REPO=Path(__file__).resolve().parents[1]
if not (REPO/"pysnspd").is_dir():REPO=Path.cwd().resolve()
if str(REPO) not in sys.path:sys.path.insert(0,str(REPO))
SOURCE=Path(os.environ.get("PYSNSPD_TEST_DISPATCHER_SOURCE",REPO/"sandbox/stage5_effective_validation/guard_protocol_dispatcher.py"))
spec=importlib.util.spec_from_file_location("sandbox.stage5_effective_validation.repository_evidence_dispatcher_test",SOURCE)
d=importlib.util.module_from_spec(spec);spec.loader.exec_module(d)
# A candidate is reviewed in scratch before being installed at this repository path.
d.ROOT=REPO
DECISION_SHA="bda09a6aaf095a175440dc77ac76faef8eaab98072fce2f5b5aa8984eb35efa6"
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write_protocol(path,entry):
    path.write_text(json.dumps(dict(nested=[dict(reference=entry)])));return path

@pytest.fixture
def isolated(tmp_path,monkeypatch):
    repo=tmp_path/"repo";repo.mkdir()
    scratch=tmp_path/"scratch";scratch.mkdir()
    baseline=repo/"audit_runs.py";baseline.write_text("# reviewed fixture baseline\n")
    monkeypatch.setattr(d,"ROOT",repo)
    monkeypatch.setattr(d.audit,"__file__",str(baseline))
    monkeypatch.setattr(d.guard,"AUDITOR_SHA",sha(baseline))
    monkeypatch.setattr(d.guard,"SCRATCH",scratch.resolve())
    return repo,scratch,baseline

def test_repository_source_is_rejected_by_default_and_allowed_only_with_explicit_opt_in(isolated):
    repo,scratch,baseline=isolated
    protocol=write_protocol(scratch/"protocol.json",dict(path=str(baseline),sha256=sha(baseline)))
    with pytest.raises(ValueError,match="scratch"):d.evidence_inputs(protocol)
    assert d.evidence_inputs(protocol,allow_repository_sources=True)=={str(protocol):sha(protocol),str(baseline):sha(baseline)}

def test_scratch_reference_keeps_existing_bound_validation(isolated):
    _,scratch,_=isolated
    evidence=scratch/"receipt.json";evidence.write_text("{}")
    protocol=write_protocol(scratch/"protocol.json",dict(path=str(evidence),sha256=sha(evidence)))
    for opt_in in [False,True]:
        assert d.evidence_inputs(protocol,allow_repository_sources=opt_in)[str(evidence)]==sha(evidence)
        evidence.write_text("{ }")
        with pytest.raises(ValueError,match="SHA"):d.evidence_inputs(protocol,allow_repository_sources=opt_in)
        evidence.write_text("{}")

@pytest.mark.parametrize("fault",["changed_hash","changed_source","other_repository_source","outside","repository_symlink_escape","scratch_symlink_escape","outside_alias","inside_alias"])
def test_bad_repository_and_scratch_references_fail_closed(isolated,tmp_path,fault):
    repo,scratch,baseline=isolated
    path=baseline;digest=sha(path)
    if fault=="changed_hash":digest="0"*64
    elif fault=="changed_source":baseline.write_text("# changed\n")
    elif fault=="other_repository_source":
        path=repo/"unreviewed.py";path.write_text("# other source\n");digest=sha(path)
    elif fault=="outside":
        path=tmp_path/"outside.py";path.write_text("# external\n");digest=sha(path)
    elif fault=="repository_symlink_escape":
        outside=tmp_path/"outside.py";outside.write_bytes(baseline.read_bytes())
        baseline.unlink();baseline.symlink_to(outside)
    elif fault=="scratch_symlink_escape":
        outside=tmp_path/"outside.json";outside.write_text("{}")
        path=scratch/"escape.json";path.symlink_to(outside);digest=sha(path)
    elif fault=="outside_alias":
        path=tmp_path/"alias.py";path.symlink_to(baseline);digest=sha(path)
    elif fault=="inside_alias":
        path=repo/"alias.py";path.symlink_to(baseline);digest=sha(path)
    protocol=write_protocol(scratch/"protocol.json",dict(path=str(path),sha256=digest))
    with pytest.raises(ValueError):d.evidence_inputs(protocol,allow_repository_sources=True)

def real_decision():
    value=os.environ.get("PYSNSPD_TEST_REAL_VOLTAGE_DECISION")
    if not value:pytest.skip("Optional actual SHA-bound decision traversal; no trajectories or RHS")
    path=Path(value);assert sha(path)==DECISION_SHA
    return path

def test_actual_adopted_voltage_decision_and_all_preserved_dependency_hashes():
    path=real_decision();value=json.loads(path.read_text())
    got=d.evidence_inputs(path,allow_repository_sources=True)
    assert got[str(path.resolve())]==DECISION_SHA
    assert len(got)==1+len(value["preserved_dependencies"])
    for entry in value["preserved_dependencies"].values():
        assert got[str(Path(entry["path"]).resolve())]==entry["sha256"]==sha(entry["path"])
    with pytest.raises(ValueError,match="scratch"):d.evidence_inputs(path)

def test_cli_only_opts_in_voltage_decision_before_independent_audits(tmp_path,monkeypatch):
    path=real_decision();value=json.loads(path.read_text())
    protocol=tmp_path/"guard.json";protocol.write_text("{}")
    resource=tmp_path/"resource.json";resource.write_text("{}")
    window=tmp_path/"window.json";window.write_text(json.dumps(dict(start_ps=0.,end_ps=36.6)))
    for name in ["reference","candidate"]:
        run=tmp_path/name;run.mkdir()
        for fn in ["manifest.json","trajectory.npz","snapshots.npz","checkpoint_latest.npz"]:(run/fn).write_bytes(b"fixture only, no physical arrays")
    helper=SimpleNamespace(__file__=str(SOURCE))
    physical=SimpleNamespace(__file__=str(SOURCE),source_integrity=lambda audit:(helper,{}))
    policy=SimpleNamespace(__file__=str(SOURCE),VOLTAGE_DECISION_SHA256=DECISION_SHA,load_decision=lambda p:(Path(p),{}))
    monkeypatch.setattr(d,"__file__",str(REPO/"sandbox/stage5_effective_validation/guard_protocol_dispatcher.py"))
    monkeypatch.setattr(d,"physical_gap_module",lambda:physical)
    monkeypatch.setattr(d,"spectral_policy_module",lambda:policy)
    monkeypatch.setattr(d,"exploratory_policy_module",lambda:policy)
    calls=[];original=d.evidence_inputs
    def recorded(p,**kw):
        calls.append((Path(p).resolve(),kw.get("allow_repository_sources",False)))
        return original(p,**kw)
    monkeypatch.setattr(d,"evidence_inputs",recorded)
    class ReachedIndependentAudit(Exception):pass
    def independent_audit(*args):raise ReachedIndependentAudit
    monkeypatch.setattr(d.audit,"audit_run",independent_audit)
    output=tmp_path/"not_created"
    with pytest.raises(ReachedIndependentAudit):
        d.main(["--kind","spectral","--reference",str(tmp_path/"reference"),"--candidate",str(tmp_path/"candidate"),
            "--guard-protocol",str(protocol),"--window-declaration",str(window),
            "--spectral-policy-decision",value["preserved_dependencies"]["prior_observable_decision"]["path"],
            "--exploratory-spectral-decision",value["preserved_dependencies"]["prior_spatial_decision"]["path"],
            "--exploratory-voltage-decision",str(path),"--control-resource-transition-receipt",str(resource),
            "--output-dir",str(output)])
    assert calls==[(protocol.resolve(),False),(path.resolve(),True),(resource.resolve(),False)]
    assert not output.exists()
