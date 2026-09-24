"""Hash Python source consistently across Windows/Linux newline conventions."""
import hashlib
from pathlib import Path

def source_sha(path):
    return hashlib.sha256(Path(path).read_bytes().replace(b'\r\n',b'\n')).hexdigest()
