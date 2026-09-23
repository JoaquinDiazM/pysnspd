"""Package the reviewed working update with exact remote preconditions."""
from pathlib import Path
import hashlib,json,subprocess,zipfile
ROOT=Path(__file__).resolve().parents[3]
TARGET=ROOT/'tmp/stage2_practical_review_20260922';TARGET.mkdir(exist_ok=True)
prior=json.loads((ROOT/'tmp/stage2_closure_prep_20260922/sync_manifest.json').read_text())
previous={r['path']:r['sha256'] for r in prior['files']};paths=set(previous)
for folder in ('docs/implementation/stage2/practical_review_20260922','sandbox/stage2_cells/practical_review_20260922'):
    paths.update(p.relative_to(ROOT).as_posix() for p in (ROOT/folder).rglob('*')
        if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc')
for p in ('README.md','docs/implementation/stage2/README.md','docs/implementation/MODELO_VIGENTE.md',
    'docs/modelo_v0_4/actualizaciones/circuito_memoria_20260922.md',
    'output/pdf/implementation/Informe_revision_criterios_y_circuito_20260922.pdf'):
    paths.add(p)
    if p not in previous:
        git=subprocess.run(['git','show','HEAD:'+p],cwd=ROOT,capture_output=True)
        previous[p]=hashlib.sha256(git.stdout).hexdigest() if git.returncode==0 else None
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
manifest=dict(base_commit=prior['base_commit'],release_commit=prior['release_commit'],
    previous_notebook_sha256=previous['docs/GEMINGA_COMMANDS.md'],
    files=[dict(path=p,sha256=sha(ROOT/p),previous_sha256=previous.get(p)) for p in sorted(paths)])
(TARGET/'sync_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
with zipfile.ZipFile(TARGET/'working_delivery.zip','w',compression=zipfile.ZIP_DEFLATED,compresslevel=4) as bundle:
    bundle.write(TARGET/'sync_manifest.json','sync_manifest.json')
    for p in paths:bundle.write(ROOT/p,'payload/'+p)
print(json.dumps(dict(files=len(paths),zip_bytes=(TARGET/'working_delivery.zip').stat().st_size,
    zip_sha256=sha(TARGET/'working_delivery.zip'),scope='Working documentation, saved results and lightweight diagnostics. No production change, commit or push.')))
