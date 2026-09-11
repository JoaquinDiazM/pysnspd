"""Write or verify the exact-byte stage-1 R2 delivery manifest (no calculation)."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / 'docs/implementation/stage1_r2/delivery_manifest.json'
RELEASE = '5ea0cd6a08d2cae414c82944f8230469ca5aa7d6'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def payload():
    files = []
    for folder in ('pysnspd/experimental', 'sandbox/stage1_catalog_r2',
                   'docs/implementation/stage1_r2'):
        files += [p for p in (ROOT/folder).rglob('*') if p.is_file()
                  and '__pycache__' not in p.parts and p != MANIFEST
                  and p.name != 'delivery_receipt.json' and p.suffix not in ('.pyc', '.log')
                  and (p.suffix != '.npz' or p.parent == ROOT/'docs/implementation/stage1_r2/catalogs')]
    files += list((ROOT/'tests').glob('test_experimental_*.py'))
    files.append(ROOT/'output/pdf/implementation/Informe_etapa_1_r2.pdf')
    return sorted(set(files))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true', help='create manifest after numerical and visual QA')
    args = parser.parse_args()
    if args.write:
        records = [{'path': p.relative_to(ROOT).as_posix(), 'bytes': p.stat().st_size,
                    'sha256': digest(p)} for p in payload()]
        MANIFEST.write_text(json.dumps({
            'schema': 'pysnspd.stage1_r2.delivery.v1', 'release_tag': 'v1.0.0',
            'release_commit_before_implementation': RELEASE,
            'previous_iteration_commit': 'd17d7c3fc3943699631126eaa636732bba56aabe',
            'scope': 'experimental data admission and uniform catalogue; no production activation',
            'files': records,
        }, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
    manifest = json.loads(MANIFEST.read_text(encoding='utf-8'))
    failures = []
    for item in manifest['files']:
        path = (ROOT/item['path']).resolve()
        if ROOT not in path.parents or not path.is_file():
            failures.append(item['path']+': missing or outside workspace')
        elif path.stat().st_size != item['bytes'] or digest(path) != item['sha256']:
            failures.append(item['path']+': content differs')
    result = {'verified': not failures, 'files': len(manifest['files']),
              'manifest_sha256': digest(MANIFEST), 'failures': failures}
    print(json.dumps(result, indent=2))
    if failures:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
