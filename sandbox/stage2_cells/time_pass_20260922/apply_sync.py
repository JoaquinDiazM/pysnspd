"""Apply an exact-byte working delivery after checking all remote preconditions.

Does not commit, pull, reset or launch physics. Backs up every overwritten file.
The incoming zip is produced from the explicitly selected local task files.
"""
from pathlib import Path
import hashlib,json,shutil,subprocess,sys,zipfile
ROOT=Path('/home/jdiaz/pysnspd').resolve()
NOTEBOOK=Path('/home/jdiaz/GEMINGA_COMMANDS.md')
sha=lambda b:hashlib.sha256(b).hexdigest()
archive=Path(sys.argv[1])
with zipfile.ZipFile(archive) as bundle:
    manifest=json.loads(bundle.read('sync_manifest.json'))
    assert subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()==manifest['base_commit']
    assert subprocess.check_output(['git','rev-parse','v1.0.0^{}'],cwd=ROOT,text=True).strip()==manifest['release_commit']
    assert sha(NOTEBOOK.read_bytes())==manifest['previous_notebook_sha256'],'External notebook changed; preserve and review'
    prepared=[]
    for record in manifest['files']:
        destination=(ROOT/record['path']).resolve()
        assert destination.is_relative_to(ROOT),'Path outside repository'
        data=bundle.read('payload/'+record['path'])
        assert sha(data)==record['sha256'],'Corrupt incoming file'
        if destination.exists():
            assert sha(destination.read_bytes()) in (record['sha256'],record['previous_sha256']), 'Unexpected remote change: '+record['path']
        else:
            assert record['previous_sha256'] is None,'Expected tracked file missing: '+record['path']
        prepared.append((destination,data))
    backup=ROOT/'tmp/stage2_cells/pre_time_pass_sync_20260922'
    assert not backup.exists(),'Backup exists: do not overwrite or repeat blindly'
    backup.mkdir(parents=True)
    for destination,data in prepared:
        if destination.exists():
            copy=backup/destination.relative_to(ROOT);copy.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(destination,copy)
        destination.parent.mkdir(parents=True,exist_ok=True)
        destination.write_bytes(data)
    shutil.copy2(NOTEBOOK,backup/'external_GEMINGA_COMMANDS.md')
    NOTEBOOK.write_bytes((ROOT/'docs/GEMINGA_COMMANDS.md').read_bytes())
    for record in manifest['files']:assert sha((ROOT/record['path']).read_bytes())==record['sha256']
    receipt=dict(status='WORKING_FILES_SYNCED_NO_COMMIT_NO_PUSH',files=len(prepared),
        base_commit=manifest['base_commit'],backup=str(backup),
        notebook_sha256=sha(NOTEBOOK.read_bytes()),zip_sha256=sha(archive.read_bytes()))
    (backup/'sync_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))
