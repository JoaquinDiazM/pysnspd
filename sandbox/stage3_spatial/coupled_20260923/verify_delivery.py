"""Read-only verification of the mixed-spatial delivery and archived predecessors."""
from pathlib import Path
import argparse
import hashlib
import json
import runpy

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT/'docs/implementation/stage3/coupled_20260923'
MANIFEST = DATA/'delivery_manifest.json'
read = lambda p: json.loads(p.read_text(encoding='utf-8'))
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()


def verify(records, overrides=None):
    for item in records:
        path = (overrides or {}).get(item['path'], ROOT/item['path'])
        if not path.is_file() or sha(path) != item['sha256']:
            raise RuntimeError('Absent or modified: '+str(path))


def predecessors():
    old = runpy.run_path(str(ROOT/'sandbox/stage3_spatial/ports_20260923/verify_delivery.py'))
    old['references']()
    before = DATA/'previous_delivery_exact'
    manifest = read(before/'delivery_manifest.json')
    overrides = {name: before/name for name in
                 ['docs/GEMINGA_COMMANDS.md', 'docs/implementation/stage3/CURRENT.md']}
    verify(manifest['files'], overrides)
    verify(manifest['references'])
    for name in ['raw/mixed_pilot', 'seeded_check', 'reservoir_postprocess_v2',
                 'reservoir_reproduction_geminga']:
        manifest = read(DATA/name/'manifest.json')
        sources = manifest.get('source_sha256', manifest.get('sources'))
        verify([dict(path=k, sha256=v) for k, v in sources.items()])
    if read(DATA/'raw/mixed_pilot/summary.json')['status'] != 'PILOT_INSTANTANEOUS_PASS_SCOPE_LIMITED':
        raise RuntimeError('Mixed pilot is not accepted')
    if read(DATA/'seeded_check/summary.json')['status'] != 'PASS_SEEDED_CAUSAL_CATALOG':
        raise RuntimeError('Corrected causal catalogue is not accepted')
    if read(DATA/'open_audit.json')['status'] != 'PASS_SAVED_OPEN_AUDIT':
        raise RuntimeError('Saved open cases have not passed review')
    checks = read(DATA/'checks.json')
    verify([dict(path=k, sha256=v) for k, v in checks['sources'].items()])
    for folder in ['reservoir_postprocess_v2', 'reservoir_reproduction_geminga']:
        summary = read(DATA/folder/'summary.json')
        if summary['status'] != 'PASS_SAVED_RESERVOIR_SNAPSHOTS_SCOPE_LIMITED':
            raise RuntimeError('Reservoir postprocess has not passed')
        for case in summary['cases']:
            for name, digest in case['artifacts_sha256'].items():
                if sha(DATA/folder/case['case_id']/name) != digest:
                    raise RuntimeError('Reservoir artifact changed: '+folder+'/'+name)
    failure = DATA/'reservoir_postprocess'
    for name, digest in read(failure/'manifest.json')['source_sha256'].items():
        if sha(failure/'sources'/name) != digest:
            raise RuntimeError('Historical preparation source changed: '+name)
    for case in read(DATA/'raw/mixed_pilot/summary.json')['cases']:
        folder = DATA/'raw/mixed_pilot'/(case['mesh']['id']+'_'+case['profile']['id'])
        for name, digest in case['artifacts_sha256'].items():
            if sha(folder/name) != digest:
                raise RuntimeError('Mixed artifact changed: '+str(folder/name))


def payload():
    files = []
    for folder in [DATA, ROOT/'sandbox/stage3_spatial/coupled_20260923']:
        files.extend(p for p in folder.rglob('*') if p.is_file() and '__pycache__' not in p.parts
                     and p.suffix != '.pyc' and p not in [MANIFEST, DATA/'delivery_receipt.json'])
    for module in ['local_energy_patch', 'monotone_energy_patch', 'seeded_count_catalog',
                   'mixed_spatial', 'spatial_dynamics', 'reservoir_spatial_dynamics']:
        files.append(ROOT/('pysnspd/experimental/'+module+'.py'))
        files.extend(ROOT.glob('tests/test_experimental_'+module+'*.py'))
    files.extend(ROOT/name for name in ['docs/GEMINGA_COMMANDS.md', 'docs/implementation/stage3/CURRENT.md',
        'output/pdf/implementation/Informe_empalme_y_acoplamiento_etapa_3_20260923.pdf'])
    return sorted(set(files))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    predecessors()
    if args.write:
        value = dict(schema='pysnspd.stage3.coupled.delivery.v1', stage3_closed=False,
                     files=[dict(path=p.relative_to(ROOT).as_posix(), sha256=sha(p),
                                 bytes=p.stat().st_size) for p in payload()])
        MANIFEST.write_text(json.dumps(value, indent=2)+'\n', encoding='utf-8', newline='\n')
    value = read(MANIFEST)
    verify(value['files'])
    external = Path('/home/jdiaz/GEMINGA_COMMANDS.md')
    if external.exists() and sha(external) != sha(ROOT/'docs/GEMINGA_COMMANDS.md'):
        raise RuntimeError('External command notebook differs')
    print(json.dumps(dict(verified=True, files=len(value['files']),
                         manifest_sha256=sha(MANIFEST), stage3_closed=False), indent=2))


if __name__ == '__main__':
    main()
