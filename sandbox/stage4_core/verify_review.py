"""Verify the result review and preserved predecessor bytes, without physics."""
from pathlib import Path
import hashlib
import json
import runpy
import sys

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT/'docs/implementation/stage4/review_20260923'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def verify(records, archive=None):
    for item in records:
        path = ROOT/item['path']
        if archive is not None and (archive/item['path']).is_file():
            path = archive/item['path']
        if sha(path) != item['sha256'] or path.stat().st_size != item['bytes']:
            raise RuntimeError('Changed delivery file: '+str(path))
    return len(records)


def main():
    old = runpy.run_path(str(ROOT/'sandbox/stage4_core/verify_delivery.py'))
    historical = old['historical_delivery']()
    predecessor = json.loads((ROOT/'docs/implementation/stage4/start_20260923/delivery_manifest.json').read_text())
    previous_count = verify(predecessor['files'], DATA/'previous_delivery_exact')
    manifest = json.loads((DATA/'delivery_manifest.json').read_text())
    assert manifest['stage4_complete'] is False and manifest['production_promotion'] is False
    count = verify(manifest['files'])
    for phase in ('local', 'spatial'):
        raw = DATA/'raw'/('stage4A_'+phase+'_20260923')
        identity = json.loads((raw/'identity.json').read_text())
        assert identity['plan_sha256'] == sha(ROOT/'docs/implementation/stage4/start_20260923/campaign_plan.json')
        for name, expected in identity['sources'].items():
            archived = DATA/'previous_delivery_exact'/name
            assert sha(archived if archived.is_file() else ROOT/name) == expected, name
        for case in json.loads((raw/'summary.json').read_text())['cases']:
            for record in case['files']:
                p = raw/record['path']
                assert sha(p) == record['sha256'] and p.stat().st_size == record['bytes']
    external = Path('/home/jdiaz/GEMINGA_COMMANDS.md')
    if sys.platform.startswith('linux') and external.exists():
        assert sha(external) == sha(ROOT/'docs/GEMINGA_COMMANDS.md')
    print(json.dumps(dict(verified=True, current_files=count, start_files=previous_count,
                         historical_delivery=historical, stage4_complete=False), indent=2))


if __name__ == '__main__':
    main()
