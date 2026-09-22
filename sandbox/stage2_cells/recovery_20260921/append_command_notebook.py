"""Preserve the existing Geminga notebook bytes and append the reviewed command."""
from pathlib import Path
import hashlib,json
ROOT=Path('/home/jdiaz/pysnspd')
OLD='2da748caccc7676eeb88557dcb5103cdf3b1c660914d174b01b47922c3ab2d7d'
sha=lambda data:hashlib.sha256(data).hexdigest()
def main():
    output=ROOT/'docs/implementation/stage2/recovery_20260921/command_notebook_update.json'
    if output.exists():raise ValueError('Existing receipt is not overwritten')
    target=Path('/home/jdiaz/GEMINGA_COMMANDS.md')
    old=target.read_bytes()
    if sha(old)!=OLD:raise ValueError('Notebook changed; review before appending')
    addition=(ROOT/'docs/implementation/stage2/recovery_20260921/command_addendum.md').read_bytes()
    new=old+b'\n'+addition
    with target.open('ab') as stream:stream.write(b'\n'+addition)
    if target.read_bytes()!=new:raise RuntimeError('Append verification failed')
    (ROOT/'docs/GEMINGA_COMMANDS.md').write_bytes(new)
    receipt=dict(previous_sha256=OLD,new_sha256=sha(new),previous_bytes=len(old),
        new_bytes=len(new),original_bytes_preserved=new[:len(old)]==old,
        authorization='Explicit user request to launch the long recovery in code_000 and detach.',
        command_addendum_sha256=sha(addition),script_sha256=sha(Path(__file__).read_bytes()))
    output.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt))
if __name__=='__main__':main()
