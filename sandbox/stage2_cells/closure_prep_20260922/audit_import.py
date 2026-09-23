"""Verify received bytes and preserve the unrelated local working files."""
from pathlib import Path
import hashlib,json,subprocess
ROOT=Path(__file__).resolve().parents[3]
DATA=ROOT/'docs/implementation/stage2/closure_prep_20260922'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
rows=[]
for line in (DATA/'remote_sha256.txt').read_text().splitlines():
    digest,remote=line.split(maxsplit=1)
    local=DATA/('raw_two/'+Path(remote).name if '/stage2_ssp_two_time_20260922/' in remote else Path(remote).name)
    rows.append(dict(remote=remote,path=local.relative_to(ROOT).as_posix(),sha256=digest,verified=sha(local)==digest,bytes=local.stat().st_size))
preserved=json.loads((ROOT/'tmp/stage2_cells/preexisting_work.json').read_text())
unrelated=[dict(path=p,verified=sha(ROOT/p)==s) for p,s in preserved.items()]
release=subprocess.check_output(['git','rev-parse','v1.0.0^{}'],cwd=ROOT,text=True).strip()
result=dict(status='PASS' if all(r['verified'] for r in rows+unrelated) and release=='5ea0cd6a08d2cae414c82944f8230469ca5aa7d6' else 'FAIL',
    imported_files=len(rows),imported_bytes=sum(r['bytes'] for r in rows),remote_imports=rows,
    unrelated_preexisting_files=unrelated,release_commit=release,physics_recomputed=False)
(DATA/'import_and_preservation_audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ('remote_imports','unrelated_preexisting_files')}))
if result['status']!='PASS':raise SystemExit(1)
