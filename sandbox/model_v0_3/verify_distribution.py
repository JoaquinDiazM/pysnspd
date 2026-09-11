"""Read-only integrity check for an unpacked delivery and its ZIP, on any host."""
from pathlib import Path
import json,hashlib,zipfile,platform
ROOT=Path(__file__).resolve().parents[2]
manifest=ROOT/'docs/modelo_v0_3/verificaciones/manifiesto_v0_3.json'
bundle=ROOT/'output/pdf/modelo_v0_3/pysnspd_documentos_v0_3.zip'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
m=json.loads(manifest.read_text(encoding='utf-8'))
failures=[]
for row in m['files']:
    p=ROOT/row['path']
    if not p.is_file() or sha(p)!=row['sha256']:failures.append(row['path'])
assert not failures,failures
with zipfile.ZipFile(bundle) as z:
    assert z.testzip() is None
    assert len(z.namelist())==len(m['files'])+1
    assert z.read(manifest.relative_to(ROOT).as_posix())==manifest.read_bytes()
    for row in m['files']:
        assert hashlib.sha256(z.read(row['path'])).hexdigest()==row['sha256'],row['path']
print(json.dumps({'status':'verified','host':platform.node(),'payload_files':len(m['files']),
 'manifest_sha256':sha(manifest),'zip_sha256':sha(bundle),'zip_bytes':bundle.stat().st_size}))
