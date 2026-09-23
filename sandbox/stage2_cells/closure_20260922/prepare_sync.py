"""Package only changed closure files; check all previously delivered task files.

The selected paths also form the explicit publication scope. Unrelated notebook
and meeting-document edits are not selected or changed.
"""
from pathlib import Path
import hashlib
import json
import subprocess
import zipfile

ROOT=Path(__file__).resolve().parents[3]
TARGET=ROOT/'tmp/stage2_final_20260922'
TARGET.mkdir(parents=True,exist_ok=True)
prior=json.loads((ROOT/'tmp/stage2_practical_review_20260922/sync_manifest.json').read_text())
previous={r['path']:r['sha256'] for r in prior['files']}
paths=set(previous)
for folder in ('docs/implementation/stage2/closure_20260922','sandbox/stage2_cells/closure_20260922'):
    paths.update(p.relative_to(ROOT).as_posix() for p in (ROOT/folder).rglob('*')
        if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc')
paths.update(('status.md','output/pdf/implementation/Informe_cierre_etapa_2.pdf'))
for name in paths-previous.keys():
    result=subprocess.run(['git','show','HEAD:'+name],cwd=ROOT,capture_output=True)
    previous[name]=hashlib.sha256(result.stdout).hexdigest() if result.returncode==0 else None
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
protected=json.loads((ROOT/'tmp/stage2_cells/preexisting_work.json').read_text())
assert not paths.intersection(protected),'Unrelated edits selected'
for name,digest in protected.items(): assert sha(ROOT/name)==digest,'Unrelated file changed: '+name
assert not subprocess.check_output(['git','diff','--cached','--name-only'],cwd=ROOT).strip(),'Index not empty'
records=[dict(path=p,sha256=sha(ROOT/p),previous_sha256=previous[p],
              bytes=(ROOT/p).stat().st_size,included=sha(ROOT/p)!=previous[p]) for p in sorted(paths)]
assert all(r['bytes']<100*1024*1024 for r in records),'Individual GitHub file exceeds limit'
manifest=dict(base_commit=prior['base_commit'],release_commit=prior['release_commit'],
    previous_notebook_sha256=previous['docs/GEMINGA_COMMANDS.md'],files=records)
(TARGET/'sync_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
(TARGET/'publication_paths.txt').write_text('\n'.join(sorted(paths))+'\n',encoding='utf-8')
with zipfile.ZipFile(TARGET/'working_delivery.zip','w',compression=zipfile.ZIP_DEFLATED,compresslevel=4) as bundle:
    bundle.write(TARGET/'sync_manifest.json','sync_manifest.json')
    for record in records:
        if record['included']: bundle.write(ROOT/record['path'],'payload/'+record['path'])
print(json.dumps(dict(selected_files=len(paths),changed_files=sum(r['included'] for r in records),
    preserved_unrelated_files=len(protected),zip_bytes=(TARGET/'working_delivery.zip').stat().st_size,
    maximum_file_bytes=max(r['bytes'] for r in records),zip_sha256=sha(TARGET/'working_delivery.zip'))))
