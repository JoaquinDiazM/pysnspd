"""Check exact files, historical evidence and nodal pilot; no physical solves."""
from pathlib import Path
import argparse
import hashlib
import json

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'docs/implementation/stage3/nodal_20260923'
REVIEW=ROOT/'docs/implementation/stage3/review_20260923'
MANIFEST=DATA/'delivery_manifest.json'
REGISTRATION_SHA='cedab95157aad073a3308aa3a9f322cfd887bbe1b45d9d594b503d05636ac0e8'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
read=lambda p:json.loads(p.read_text(encoding='utf-8'))


def verify(records,root=ROOT):
    for rec in records:
        p=(root/rec['path']).resolve()
        if root not in p.parents or not p.is_file() or sha(p)!=rec['sha256']:
            raise RuntimeError('Changed or absent: '+str(p))


def evidence():
    if sha(DATA/'registration.json')!=REGISTRATION_SHA:raise RuntimeError('Registration changed')
    older=read(ROOT/'docs/implementation/stage2/delivery_manifest.json')
    verify(older['files']+older['retained_references'])
    snapshot=REVIEW/'previous_delivery_exact'
    initial=read(snapshot/'delivery_manifest.json')
    verify(initial['files'],snapshot)
    verify(initial['retained_references'])
    summary=read(DATA/'pilot/summary.json');pilot=read(DATA/'pilot/weak_phase_thermal_n8.json')
    receipt=read(DATA/'pilot/weak_phase_thermal_n8.receipt.json')
    handoff=read(DATA/'handoff.json')
    if summary['status']!='PILOT_PASS_FULL_CAMPAIGN_PENDING':raise RuntimeError('Pilot not passed')
    sources=summary['source_hashes']
    for value in (pilot,receipt,handoff):
        if value['source_hashes']!=sources:raise RuntimeError('Nodal sources disagree')
    verify([dict(path=k,sha256=v) for k,v in sources.items()])
    for key,name in [('json_sha256','weak_phase_thermal_n8.json'),
                     ('arrays_sha256','weak_phase_thermal_n8.npz'),
                     ('initial_sha256','weak_phase_thermal_n8_initial.npz')]:
        if sha(DATA/'pilot'/name)!=receipt[key]:raise RuntimeError('Pilot file altered: '+name)
    tests=read(DATA/'checks/result.json')
    if tests['status']!='PASS' or tests['exit_code']!=0:raise RuntimeError('Focused checks failed')
    verify([dict(path=k,sha256=v) for k,v in tests['tests'].items()])
    if summary['stage3_closed'] or handoff['long_run_launched_by_agent']:
        raise RuntimeError('Unexpected scope or long run')
    return sources


def payload():
    files=[]
    for folder in [DATA,REVIEW,ROOT/'sandbox/stage3_spatial/review_20260923']:
        files.extend(p for p in folder.rglob('*') if p.is_file() and '__pycache__' not in p.parts
                     and p.suffix!='.pyc' and p not in [MANIFEST,DATA/'delivery_receipt.json'])
    files.extend(ROOT/p for p in [
        'pysnspd/experimental/spatial_nodal.py','tests/test_experimental_spatial_nodal.py',
        'tests/test_stage3_nodal_runner.py','sandbox/stage3_spatial/run_nodal_batch.py',
        'sandbox/stage3_spatial/run_nodal_checks.py','sandbox/stage3_spatial/verify_nodal_delivery.py',
        'sandbox/stage3_spatial/check_nodal_handoff.py',
        'docs/implementation/stage3/CURRENT.md','docs/GEMINGA_COMMANDS.md',
        'output/pdf/implementation/Informe_revision_espacial_etapa_3_20260923.pdf'])
    return sorted(set(files))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write',action='store_true')
    args=parser.parse_args();sources=evidence()
    if args.write:
        value=dict(schema='pysnspd.stage3.nodal.delivery.v1',stage3_closed=False,
                   files=[dict(path=p.relative_to(ROOT).as_posix(),sha256=sha(p),bytes=p.stat().st_size)
                          for p in payload()],
                   references=[dict(path=p,sha256=v) for p,v in sources.items()])
        MANIFEST.write_text(json.dumps(value,indent=2)+'\n',encoding='utf-8',newline='\n')
    value=read(MANIFEST);verify(value['files']+value['references'])
    print(json.dumps(dict(verified=True,files=len(value['files']),references=len(value['references']),
                         manifest_sha256=sha(MANIFEST),stage3_closed=False),indent=2))


if __name__=='__main__':main()
