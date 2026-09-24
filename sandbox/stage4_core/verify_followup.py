"""Read-only stage4 follow-up check, including preserved historical deliveries."""
from contextlib import redirect_stdout
from functools import wraps
from io import StringIO
from pathlib import Path
from types import SimpleNamespace
import hashlib
import inspect
import json
import runpy
import sys

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'docs/implementation/stage4/followup_20260923'
ARCHIVE=DATA/'previous_delivery_exact'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding='utf8'))


def predecessor():
    archived={p.relative_to(ARCHIVE).as_posix():p for p in ARCHIVE.rglob('*') if p.is_file()}
    def load(path):
        ns=runpy.run_path(str(path)); gl=ns['main'].__globals__
        original=ns.get('verify')
        if original is not None:
            signature=inspect.signature(original)
            @wraps(original)
            def historical_verify(records,*args,**kwargs):
                bound=signature.bind(records,*args,**kwargs);bound.apply_defaults()
                if 'overrides' in signature.parameters:
                    lookup=Path(bound.arguments.get('root',ROOT)).resolve()
                    if lookup==ROOT:
                        bound.arguments['overrides']=archived | (bound.arguments.get('overrides') or {})
                    return original(*bound.args,**bound.kwargs)
                if 'archive' in signature.parameters:
                    previous=bound.arguments.get('archive')
                    for item in records:
                        name=item['path']
                        p=(Path(previous)/name if previous is not None and (Path(previous)/name).is_file()
                           else archived.get(name,ROOT/name))
                        if sha(p)!=item['sha256'] or p.stat().st_size!=item['bytes']:
                            raise RuntimeError('Changed historical delivery file: '+str(p))
                    return len(records)
                return original(records,*args,**kwargs)
            ns['verify']=gl['verify']=historical_verify
        ns['runpy']=gl['runpy']=SimpleNamespace(run_path=load)
        return ns
    path=ROOT/'sandbox/stage4_core/verify_review.py'
    previous=load(path);output=StringIO();arguments=sys.argv
    try:
        sys.argv=[str(path)]
        with redirect_stdout(output):previous['main']()
    finally:sys.argv=arguments
    return json.loads(output.getvalue())


def main():
    previous=predecessor()
    manifest=read(DATA/'delivery_manifest.json')
    assert manifest['stage4_complete'] is False and manifest['production_changed'] is False
    for item in manifest['files']:
        path=ROOT/item['path']
        assert sha(path)==item['sha256'] and path.stat().st_size==item['bytes'],str(path)
    raw=DATA/'raw/stage4A_core_followup_20260923';identity=read(raw/'identity.json')
    assert identity['plan_sha256']==sha(ROOT/'docs/implementation/stage4/review_20260923/next_campaign_plan.json')
    for name,expected in identity['sources'].items():
        path=ARCHIVE/name if (ARCHIVE/name).is_file() else ROOT/name
        assert sha(path)==expected,name
    for case in read(raw/'summary.json')['cases']:
        for item in case['files']:
            p=raw/item['path'];assert sha(p)==item['sha256'] and p.stat().st_size==item['bytes']
    reference=DATA/'raw/stage4_radial_reference_20260923'
    radial_identity=read(reference/'identity.json')
    assert radial_identity['plan_sha256']==sha(DATA/'radial_reference_plan.json')
    for name,expected in radial_identity['sources'].items(): assert sha(ROOT/name)==expected,name
    radial=read(reference/'summary.json')
    assert len(radial['task_records'])==518 and not radial['physical_core_admitted']
    for record in radial['task_records']:
        assert sha(reference/record['arrays'])==record['sha256'],record['id']
    for record in radial['cases'].values():
        assert sha(reference/record['fields_file'])==record['fields_sha256']
    pilot=DATA/'pilot/stage4_radial_pilot_20260923'
    assert read(pilot/'identity.json')['plan_sha256']==sha(DATA/'pilot/radial_pilot_plan.json')
    for record in read(pilot/'summary.json')['task_records']:
        assert sha(pilot/record['arrays'])==record['sha256']
    external=Path('/home/jdiaz/GEMINGA_COMMANDS.md')
    if sys.platform.startswith('linux') and external.exists():
        assert sha(external)==sha(ROOT/'docs/GEMINGA_COMMANDS.md')
    print(json.dumps(dict(verified=True,files=len(manifest['files']),previous=previous,
                         stage4_complete=False,production_changed=False),indent=2))


if __name__=='__main__':main()
