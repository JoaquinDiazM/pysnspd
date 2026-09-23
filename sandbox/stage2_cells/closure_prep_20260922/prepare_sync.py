"""Package this task's exact files; preserve the previous working delivery."""
from pathlib import Path
import hashlib,json,zipfile
ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).parent
TARGET=ROOT/'tmp/stage2_closure_prep_20260922';TARGET.mkdir(exist_ok=True)
prior=json.loads((ROOT/'tmp/stage2_time_pass_20260922/sync_manifest.json').read_text())
previous={r['path']:r['sha256'] for r in prior['files']}
paths=set(previous)
for folder in ('docs/implementation/stage2/closure_prep_20260922','sandbox/stage2_cells/closure_prep_20260922'):
    paths.update(p.relative_to(ROOT).as_posix() for p in (ROOT/folder).rglob('*')
        if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc')
paths.add('output/pdf/implementation/Informe_validacion_dos_celdas_etapa_2_20260922.pdf')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
manifest=dict(base_commit=prior['base_commit'],release_commit=prior['release_commit'],
    previous_notebook_sha256=previous['docs/GEMINGA_COMMANDS.md'],
    files=[dict(path=p,sha256=sha(ROOT/p),previous_sha256=previous.get(p)) for p in sorted(paths)])
(TARGET/'sync_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
with zipfile.ZipFile(TARGET/'working_delivery.zip','w',compression=zipfile.ZIP_DEFLATED,compresslevel=5) as bundle:
    bundle.write(TARGET/'sync_manifest.json','sync_manifest.json')
    for p in paths:bundle.write(ROOT/p,'payload/'+p)
print(json.dumps(dict(files=len(paths),zip_bytes=(TARGET/'working_delivery.zip').stat().st_size,
    zip_sha256=sha(TARGET/'working_delivery.zip'),scope='Working files only. No commit, push, long computation or production promotion.')))
