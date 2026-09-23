"""Verify stage 3 ports delivery and immutable predecessors, without physics."""
from pathlib import Path
import argparse
import hashlib
import json

ROOT=Path(__file__).resolve().parents[3]
DATA=ROOT/'docs/implementation/stage3/ports_20260923'
MANIFEST=DATA/'delivery_manifest.json'
read=lambda p:json.loads(p.read_text(encoding='utf-8'))
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()


def verify(records,root=ROOT,overrides=None):
    for record in records:
        name=record['path'];path=(root/name).resolve()
        if overrides and name in overrides:path=overrides[name]
        if not path.is_file() or sha(path)!=record['sha256']:
            raise RuntimeError('Absent or modified: '+str(path))


def references():
    stage2=read(ROOT/'docs/implementation/stage2/delivery_manifest.json')
    verify(stage2['files']+stage2['retained_references'])
    initial_snapshot=ROOT/'docs/implementation/stage3/review_20260923/previous_delivery_exact'
    initial=read(initial_snapshot/'delivery_manifest.json')
    verify(initial['files'],initial_snapshot);verify(initial['retained_references'])
    previous=read(DATA/'previous_delivery_exact/delivery_manifest.json')
    overrides={name:DATA/'previous_delivery_exact'/name for name in
               ['docs/GEMINGA_COMMANDS.md','docs/implementation/stage3/CURRENT.md']}
    verify(previous['files'],overrides=overrides);verify(previous['references'])
    pilot=read(DATA/'open_pilot/manifest.json');summary=read(DATA/'open_pilot/summary.json')
    if summary['status']!='PILOT_COMPLETE_SCOPE_LIMITED' or summary['stage3_closed']:
        raise RuntimeError('Unexpected open pilot scope')
    checks=read(DATA/'checks.json')
    if checks['status']!='PASS' or checks['passed_tests']!=168 or checks['passed_subtests']!=23:
        raise RuntimeError('New tests not passed')
    verify([dict(path=p,sha256=v) for p,v in checks['sources'].items()])
    values=dict(pilot['sources'])
    values.update(read(DATA/'electrical/electrical_diagnostics.json')['source_sha256'])
    # Runtime dependencies kept unchanged; include their exact bytes too.
    for name in ['pysnspd/experimental/cell_transport.py','pysnspd/experimental/cell_validation.py']:
        values[name]=sha(ROOT/name)
    records=[dict(path=p,sha256=v) for p,v in values.items()]
    verify(records)
    if read(DATA/'validation_review.json')['status']!='PASS_REVIEWED_SAVED_EVIDENCE':
        raise RuntimeError('Independent saved-data review not passed')
    return records


def payload():
    result=[]
    for folder in [DATA,ROOT/'sandbox/stage3_spatial/ports_20260923']:
        result.extend(p for p in folder.rglob('*') if p.is_file() and '__pycache__' not in p.parts
                      and p.suffix!='.pyc' and p not in [MANIFEST,DATA/'delivery_receipt.json'])
    result.extend(ROOT/('pysnspd/experimental/'+n+'.py') for n in
                  ['spatial_open','electrical_ports','superconducting_reservoir','reservoir_transport'])
    result.extend(ROOT/('tests/test_experimental_'+n+'.py') for n in
                  ['spatial_open','electrical_ports','superconducting_reservoir','reservoir_transport'])
    result.extend(ROOT/n for n in ['docs/GEMINGA_COMMANDS.md','docs/implementation/stage3/CURRENT.md',
        'output/pdf/implementation/Informe_avance_bordes_etapa_3_20260923.pdf'])
    return sorted(set(result))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write',action='store_true')
    args=parser.parse_args();refs=references()
    if args.write:
        value=dict(schema='pysnspd.stage3.ports.delivery.v1',stage3_closed=False,
            files=[dict(path=p.relative_to(ROOT).as_posix(),sha256=sha(p),bytes=p.stat().st_size) for p in payload()],
            references=refs)
        MANIFEST.write_text(json.dumps(value,indent=2)+'\n',encoding='utf-8',newline='\n')
    value=read(MANIFEST);verify(value['files']+value['references'])
    external=Path('/home/jdiaz/GEMINGA_COMMANDS.md')
    if external.exists() and sha(external)!=sha(ROOT/'docs/GEMINGA_COMMANDS.md'):
        raise RuntimeError('External Geminga command notebook differs')
    print(json.dumps(dict(verified=True,files=len(value['files']),references=len(value['references']),
                         manifest_sha256=sha(MANIFEST),stage3_closed=False),indent=2))


if __name__=='__main__':main()
