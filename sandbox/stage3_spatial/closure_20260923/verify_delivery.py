"""Read-only verification of scoped stage-3 closure and the stage-3.5 handoff.

No --write option, physical calculation or implicit regeneration. The previous
delivery is checked against its preserved front-door copies, not today's text.
"""
from pathlib import Path, PurePosixPath
import argparse
import hashlib
import json
import runpy

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT/'docs/implementation/stage3/closure_20260923'
PREVIOUS = DATA/'previous_delivery_exact'
MANIFEST = DATA/'delivery_manifest.json'
COUPLED_MANIFEST = 'docs/implementation/stage3/coupled_20260923/delivery_manifest.json'


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def project_path(name):
    if not isinstance(name, str) or not name or '\\' in name:
        raise ValueError('Expected a nonempty repository-relative POSIX path')
    relative = PurePosixPath(name)
    if relative.is_absolute() or '..' in relative.parts or ':' in relative.parts[0]:
        raise ValueError('Manifest path must remain inside the repository: '+name)
    path = (ROOT/name).resolve()
    path.relative_to(ROOT)
    return path


def entries(value):
    """Canonical files list; entries is accepted as an unambiguous alias."""
    if 'files' in value and 'entries' in value and value['files'] != value['entries']:
        raise ValueError('Manifest files and entries disagree')
    result = value.get('files', value.get('entries'))
    if not isinstance(result, list) or not result:
        raise ValueError('Manifest requires a nonempty files (or entries) list')
    return result


def verify(records, overrides=None, require_size=False):
    if not isinstance(records, list):
        raise ValueError('Manifest records must be a list')
    seen = set()
    for record in records:
        name = record['path']
        expected = project_path(name)
        if name in seen:
            raise ValueError('Duplicate manifest path: '+name)
        seen.add(name)
        path = (overrides or {}).get(name, expected)
        digest = record['sha256']
        if (not isinstance(digest, str) or len(digest) != 64
                or any(c not in '0123456789abcdef' for c in digest)):
            raise ValueError('Invalid SHA256 for '+name)
        if not path.is_file() or sha(path) != digest:
            raise RuntimeError('Absent or modified: '+str(path))
        if require_size and 'bytes' not in record:
            raise ValueError('New delivery record requires byte size: '+name)
        if 'bytes' in record:
            size = record['bytes']
            if isinstance(size, bool) or not isinstance(size, int) or size < 0:
                raise ValueError('Invalid byte size for '+name)
            if path.stat().st_size != size:
                raise RuntimeError('Byte size differs: '+str(path))
    return len(records)


def predecessors():
    old = runpy.run_path(str(ROOT/'sandbox/stage3_spatial/coupled_20260923/verify_delivery.py'))
    old['predecessors']()
    previous = read(PREVIOUS/'delivery_manifest.json')
    overrides = {name: PREVIOUS/name for name in
                 ('docs/GEMINGA_COMMANDS.md', 'docs/implementation/stage3/CURRENT.md')}
    overrides[COUPLED_MANIFEST] = PREVIOUS/'delivery_manifest.json'
    count = verify(entries(previous), overrides=overrides)
    for key in ('references', 'retained_references'):
        count += verify(previous.get(key, []), overrides=overrides)
    return count


def require_flags(value, expected, label):
    for key, flag in expected.items():
        if value.get(key) is not flag:
            raise RuntimeError(f'{label}: {key} must be {flag}')


def scope():
    decision = read(DATA/'closure_decision.json')
    require_flags(decision, dict(development_closed=True,
        original_stage3_full_contract_complete=False, strict_numerical_admission=False,
        temporal_admission=False, full_D27_admission=False, kinetic_interface_admission=False,
        production_promotion=False, production_solver_changed=False), 'Development closure')
    if decision.get('status') != 'DEVELOPMENT_CLOSED_USER_AUTHORIZED_WITH_LIMITS':
        raise RuntimeError('Unexpected development closure status')
    if not decision.get('pending_before_stage4_5'):
        raise RuntimeError('Carried-over obligations are missing')
    review = read(DATA/decision['evidence'])
    if review.get('status') != 'VERIFIED_SIX_SNAPSHOT_CAMPAIGN_SCOPE_LIMITED' or review.get('issues') != []:
        raise RuntimeError('The full saved-data audit has not passed')
    if len(review.get('cases', [])) != 6:
        raise RuntimeError('Six reviewed spatial snapshots are required')
    require_flags(review['decision'], dict(development_closure_supported=True,
        stage3_complete_physical_dynamic_admission=False, temporal_admission=False,
        full_D27_admission=False, kinetic_interface_admission=False, production=False), 'Audit scope')
    raw = DATA/decision['completed_user_batch']['path']
    batch = read(raw/'summary.json')
    if batch.get('status') != 'RESERVOIR_INSTANTANEOUS_BATCH_PASS_SCOPE_LIMITED':
        raise RuntimeError('The user batch is incomplete')
    require_flags(batch, dict(stage3_closed=False, temporal_admission=False, production=False), 'Original raw batch')
    for folder in ('mixed_control', 'reservoir'):
        summary = read(raw/folder/'summary.json')
        require_flags(summary, dict(stage3_closed=False, temporal_admission=False, production=False), folder)
        if len(summary.get('cases', [])) != 6:
            raise RuntimeError('Incomplete case coverage: '+folder)
    research_path = ROOT/'docs/implementation/stage3_5/entry_contract.json'
    research = read(research_path)
    require_flags(research, dict(new_simulation_campaign_executed=False,
        final_parameter_ranges_adopted=False, new_physics_activated=False,
        production_promotion=False), 'Stage3.5 research scope')
    if not research.get('status', '').startswith('RESEARCH_OPEN'):
        raise RuntimeError('Stage3.5 must remain an open research stage')
    linked = (research_path.parent/research['stage3_development_closure']).resolve()
    if linked != (DATA/'closure_decision.json').resolve():
        raise RuntimeError('Stage3.5 references a different development closure')
    if not research.get('carryover_requirements_before_stage4_5'):
        raise RuntimeError('Stage3.5 dropped the pending implementation requirements')
    return decision, research


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    old_count = predecessors()
    value = read(MANIFEST)
    records = entries(value)
    current_count = verify(records, require_size=True)
    reference_count = 0
    for key in ('references', 'retained_references'):
        reference_count += verify(value.get(key, []))
    decision, research = scope()
    # Optional summary flags must not contradict the independently read contract.
    expected = dict(development_closed=True, original_stage3_full_contract_complete=False,
                    temporal_admission=False, full_D27_admission=False,
                    strict_numerical_admission=False, production_promotion=False)
    require_flags(value, {key: flag for key, flag in expected.items() if key in value}, 'Delivery metadata')
    external = Path('/home/jdiaz/GEMINGA_COMMANDS.md')
    if external.exists() and sha(external) != sha(ROOT/'docs/GEMINGA_COMMANDS.md'):
        raise RuntimeError('External Geminga command notebook differs from the repository')
    print(json.dumps(dict(verified=True, files=current_count, references=reference_count,
        coupled_predecessor_records=old_count, manifest_sha256=sha(MANIFEST),
        development_closed=decision['development_closed'], original_stage3_full_contract_complete=False,
        temporal_admission=False, full_D27_admission=False, production_promotion=False,
        next_stage=research['status'], physical_calculations_executed=0), indent=2))


if __name__ == '__main__':
    main()
