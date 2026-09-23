"""Verify stage-3 start delivery integrity; this does not run physical checks."""
from pathlib import Path
import argparse
import hashlib
import json

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT/'docs/implementation/stage3/iteration_20260923'
MANIFEST = DATA/'delivery_manifest.json'
REGISTRATION_SHA = '05fac6d65bbc7a3ab829d22ac3e7e0404cb670cc782bd059d3e5abdc50e09d35'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def verify_records(records):
    for record in records:
        path = (ROOT/record['path']).resolve()
        require(ROOT in path.parents and path.is_file(), 'Absent/outside root: '+record['path'])
        require(digest(path) == record['sha256'], 'Changed: '+record['path'])


def verify_evidence():
    require(digest(DATA/'registration.json') == REGISTRATION_SHA, 'Registration changed')
    retained = read(ROOT/'docs/implementation/stage2/delivery_manifest.json')
    verify_records(retained['files']+retained['retained_references'])
    summary = read(DATA/'pilot/summary.json')
    pilot = read(DATA/'pilot/weak_phase_thermal_n8.json')
    receipt = read(DATA/'pilot/weak_phase_thermal_n8.receipt.json')
    handoff = read(DATA/'handoff.json')
    require(summary['status'] == 'PILOT_PASS_FULL_CAMPAIGN_PENDING', 'Unexpected pilot scope')
    require(summary['stage3_closed'] is False and handoff['stage3_closed'] is False,
            'This delivery cannot close stage 3')
    require(handoff['long_run_launched_by_agent'] is False, 'Unexpected long run')
    sources = summary['source_hashes']
    for label, data in [('pilot', pilot), ('receipt', receipt), ('handoff', handoff)]:
        require(data['source_hashes'] == sources, 'Sources disagree: '+label)
    verify_records([dict(path=k, sha256=v) for k, v in sources.items()])
    for key, name in [('json_sha256', 'weak_phase_thermal_n8.json'),
                      ('arrays_sha256', 'weak_phase_thermal_n8.npz'),
                      ('initial_sha256', 'weak_phase_thermal_n8_initial.npz')]:
        require(digest(DATA/'pilot'/name) == receipt[key], 'Pilot receipt mismatch: '+name)
    tests = read(DATA/'checks/result.json')
    require(tests['status'] == 'PASS' and tests['exit_code'] == 0, 'Tests did not pass')
    verify_records([dict(path=k, sha256=v) for k, v in tests['tests'].items()])
    return sources


def payload():
    paths = []
    for folder in [DATA, ROOT/'sandbox/stage3_spatial']:
        paths.extend(p for p in folder.rglob('*') if p.is_file()
                     and '__pycache__' not in p.parts and p.suffix != '.pyc'
                     and p.name not in ('delivery_manifest.json', 'delivery_receipt.json'))
    paths.extend(ROOT/p for p in [
        'pysnspd/experimental/spatial_functional.py',
        'tests/test_experimental_spatial_functional.py',
        'tests/test_stage3_progress.py', 'tests/test_stage3_runner.py',
        'docs/implementation/stage3/CURRENT.md',
        'docs/GEMINGA_COMMANDS.md',
        'output/pdf/implementation/Informe_inicio_etapa_3_20260923.pdf',
    ])
    return sorted(set(paths))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true', help='Create manifest after evidence checks')
    args = parser.parse_args()
    sources = verify_evidence()
    if args.write:
        files = [dict(path=p.relative_to(ROOT).as_posix(), bytes=p.stat().st_size,
                      sha256=digest(p)) for p in payload()]
        references = dict(sources)
        for p in ['docs/implementation/stage2/delivery_manifest.json',
                  'docs/implementation/stage3/entry_contract.json']:
            references[p] = digest(ROOT/p)
        value = dict(schema='pysnspd.stage3.start.delivery.v1',
                     scope='Static eight-cell pilot only; full campaign pending user execution.',
                     stage3_closed=False, production_promotion=False,
                     files=files, retained_references=[dict(path=p, sha256=v)
                         for p, v in sorted(references.items())])
        MANIFEST.write_text(json.dumps(value, indent=2)+'\n', encoding='utf-8')
    manifest = read(MANIFEST)
    verify_records(manifest['files']+manifest['retained_references'])
    print(json.dumps(dict(verified=True, files=len(manifest['files']),
                         references=len(manifest['retained_references']),
                         manifest_sha256=digest(MANIFEST), stage3_closed=False), indent=2))


if __name__ == '__main__':
    main()
