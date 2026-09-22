"""Archive the user's rerun, then restore the exact historical static artifact."""
from pathlib import Path
import hashlib,json,subprocess
ROOT=Path('/home/jdiaz/pysnspd')
NAME='docs/implementation/stage2/resume_20260921/selected_continuum.json'
EXPECTED='4a30117d509ada1290021b38d64c10bbe0bd448f4aad19958e72e8d0c91db6d8'
sha=lambda data:hashlib.sha256(data).hexdigest()
def main():
    target=ROOT/NAME;actual=target.read_bytes()
    if sha(actual)!=EXPECTED:raise ValueError('Rerun artifact changed; do not overwrite')
    frozen=subprocess.check_output(['git','show','HEAD:'+NAME],cwd=ROOT)
    new,old=json.loads(actual),json.loads(frozen)
    changed=[k for k in new.keys()|old.keys() if new.get(k)!=old.get(k)]
    if changed!=['runtime_seconds']:raise ValueError('Unexpected changes beyond runtime')
    archive=ROOT/'docs/implementation/stage2/review_20260922/raw/diagnostics/selected_continuum.json'
    archive.parent.mkdir(parents=True,exist_ok=True)
    if archive.exists() and archive.read_bytes()!=actual:raise ValueError('Existing archive differs')
    archive.write_bytes(actual)
    target.write_bytes(frozen)
    print(json.dumps(dict(rerun_sha256=sha(actual),historical_sha256=sha(frozen),
        rerun_archived=True,only_changed_field='runtime_seconds')))
if __name__=='__main__':main()
