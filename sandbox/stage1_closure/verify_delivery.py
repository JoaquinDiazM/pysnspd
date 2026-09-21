"""Verify the exact-byte closure delivery and its retained R2 catalogue."""
from pathlib import Path
import argparse
import hashlib
import json

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT/'docs/implementation/stage1_closure'
MANIFEST = DATA/'delivery_manifest.json'
BASELINE = 'c0e07f9bc20ba2a75be96faecf6a5dfc4a58d815'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def payload():
    paths = []
    for folder in (DATA, ROOT/'sandbox/stage1_closure'):
        paths.extend(p for p in folder.rglob('*') if p.is_file()
                     and '__pycache__' not in p.parts
                     and p.name not in ('delivery_manifest.json', 'delivery_receipt.json')
                     and p.suffix not in ('.log', '.pyc'))
    paths.extend(ROOT/p for p in (
        'pysnspd/experimental/cell_validation.py',
        'pysnspd/experimental/material_preprocessing.py',
        'tests/test_experimental_cell_validation.py',
        'tests/test_experimental_material_preprocessing.py',
        'output/pdf/implementation/Informe_cierre_etapa_1.pdf',
    ) if (ROOT/p).is_file())
    return sorted(set(paths))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    if args.write:
        retained = json.loads((ROOT/'docs/implementation/stage1_r2/delivery_manifest.json').read_text())
        for record in retained['files']:
            if digest(ROOT/record['path']) != record['sha256']:
                raise RuntimeError('Historical R2 evidence changed: '+record['path'])
        files = [dict(path=p.relative_to(ROOT).as_posix(), bytes=p.stat().st_size,
                      sha256=digest(p)) for p in payload()]
        references = [dict(path=p, sha256=digest(ROOT/p)) for p in (
            'docs/implementation/stage1_r2/catalogs/occupation_catalog.npz',
            'docs/implementation/stage1_r2/catalog_admission.json',
            'docs/implementation/stage1_r2/delivery_manifest.json',
            'pysnspd/experimental/energy_catalog.py',
            'pysnspd/experimental/material_admission.py',
        )]
        MANIFEST.write_text(json.dumps(dict(
            schema='pysnspd.stage1_closure.delivery.v1', baseline_commit=BASELINE,
            immutable_release_commit='5ea0cd6a08d2cae414c82944f8230469ca5aa7d6',
            scope='Experimental catalogue use in cells and explicit material preprocessing; no production activation.',
            retained_references=references, files=files,
        ), indent=2)+'\n', encoding='utf-8')
    manifest = json.loads(MANIFEST.read_text())
    failures = []
    for record in manifest['files']+manifest['retained_references']:
        path = (ROOT/record['path']).resolve()
        if ROOT not in path.parents or not path.is_file():
            failures.append(record['path']+': absent/outside workspace')
        elif digest(path) != record['sha256']:
            failures.append(record['path']+': hash mismatch')
    result = dict(verified=not failures, files=len(manifest['files']),
                  references=len(manifest['retained_references']),
                  manifest_sha256=digest(MANIFEST), failures=failures)
    print(json.dumps(result, indent=2))
    if failures:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
