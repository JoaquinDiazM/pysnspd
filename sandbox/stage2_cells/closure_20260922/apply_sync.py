"""Apply final closure update after verifying remote hashes, HEAD and release.

Only task files are updated. Every replaced file is backed up. No Git state is
changed, and no physical calculation is launched.
"""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess
import sys
import zipfile

ROOT=Path('/home/jdiaz/pysnspd').resolve()
NOTEBOOK=Path('/home/jdiaz/GEMINGA_COMMANDS.md')
sha=lambda b:hashlib.sha256(b).hexdigest()
with zipfile.ZipFile(Path(sys.argv[1])) as bundle:
    manifest=json.loads(bundle.read('sync_manifest.json'))
    assert subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()==manifest['base_commit']
    assert subprocess.check_output(['git','rev-parse','v1.0.0^{}'],cwd=ROOT,text=True).strip()==manifest['release_commit']
    assert not subprocess.check_output(['git','diff','--cached','--name-only'],cwd=ROOT).strip(),'Remote index not empty'
    assert sha(NOTEBOOK.read_bytes())==manifest['previous_notebook_sha256'],'External notebook changed'
    prepared=[]
    for record in manifest['files']:
        target=(ROOT/record['path']).resolve()
        assert target.is_relative_to(ROOT),'Path escapes repository'
        if target.exists():
            assert sha(target.read_bytes()) in (record['sha256'],record['previous_sha256']),record['path']
        else: assert record['previous_sha256'] is None,'Expected file missing: '+record['path']
        if record['included']:
            payload=bundle.read('payload/'+record['path'])
            assert sha(payload)==record['sha256'],record['path']
            prepared.append((target,payload))
    backup=ROOT/'tmp/stage2_cells/pre_final_closure_sync_20260922'
    assert not backup.exists(),'Backup exists; do not overwrite or repeat blindly'
    backup.mkdir(parents=True)
    for target,payload in prepared:
        if target.exists():
            dest=backup/target.relative_to(ROOT)
            dest.parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(target,dest)
        target.parent.mkdir(parents=True,exist_ok=True)
        target.write_bytes(payload)
    shutil.copy2(NOTEBOOK,backup/'external_GEMINGA_COMMANDS.md')
    NOTEBOOK.write_bytes((ROOT/'docs/GEMINGA_COMMANDS.md').read_bytes())
    for record in manifest['files']:
        assert sha((ROOT/record['path']).read_bytes())==record['sha256'],record['path']
    receipt=dict(status='FINAL_CLOSURE_WORKING_FILES_SYNCED',verified_files=len(manifest['files']),
        changed_files=len(prepared),backup=str(backup),base_commit=manifest['base_commit'],
        notebook_sha256=sha(NOTEBOOK.read_bytes()))
    (backup/'sync_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))
