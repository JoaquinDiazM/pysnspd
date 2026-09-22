"""Apply the user-requested notebook cleanup after archiving its original bytes."""
from pathlib import Path
import hashlib,json
ROOT=Path('/home/jdiaz/pysnspd')
DATA=ROOT/'docs/implementation/stage2/review_20260922'
OLD='67b3556fe57507575568139357204d62d9bd4d74eddd68a26de1d74bcdbb03b9'
NEW='8647d223b12364fd8474c45dd54ff3d103ee7814a672367616ff041c68a01303'
sha=lambda data:hashlib.sha256(data).hexdigest()
def main():
    target=Path('/home/jdiaz/GEMINGA_COMMANDS.md')
    draft=ROOT/'tmp/stage2_cells/GEMINGA_COMMANDS_clean_draft.md'
    receipt=DATA/'command_notebook_cleanup.json'
    if receipt.exists():raise ValueError('Existing cleanup receipt is preserved')
    before=target.read_bytes();after=draft.read_bytes()
    if sha(before)!=OLD or sha(after)!=NEW:raise ValueError('Notebook/draft changed; review before applying')
    if b'screen' in after.lower():raise ValueError('Active notebook still contains session instructions')
    DATA.mkdir(parents=True,exist_ok=True)
    archive=DATA/'GEMINGA_COMMANDS_before_cleanup.md'
    if archive.exists() and archive.read_bytes()!=before:raise ValueError('Archive differs from current notebook')
    if not archive.exists():archive.write_bytes(before)
    temporary=target.with_name('GEMINGA_COMMANDS.review20260922.tmp')
    with temporary.open('xb') as stream:stream.write(after)
    temporary.replace(target)
    (ROOT/'docs/GEMINGA_COMMANDS.md').write_bytes(after)
    if target.read_bytes()!=after or archive.read_bytes()!=before:raise RuntimeError('Cleanup verification failed')
    record=dict(schema='pysnspd.command_notebook_cleanup.v1',status='APPLIED',
        requested_action='Remove session-management instructions and simplify the active command notebook.',
        old_sha256=OLD,new_sha256=NEW,original_archive=str(archive),original_bytes_preserved=True,
        before_lines=len(before.splitlines()),after_lines=len(after.splitlines()),
        active_screen_mentions=0,active_long_commands=1,long_command_executed=False,
        script_sha256=sha(Path(__file__).read_bytes()))
    receipt.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(record))
if __name__=='__main__':main()
