"""Record a full, source-bound regression of the closure implementation.

Invoke under the external 240-second timeout documented for Geminga. This
runner never retries a timed-out command or starts a background calculation.
"""
from pathlib import Path
import hashlib
import json
import platform
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'docs/implementation/stage1_closure'


def hashes():
    paths = sorted((ROOT/'pysnspd/experimental').glob('*.py'))
    paths += sorted((ROOT/'tests').glob('test_experimental_*.py'))
    return {p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
            for p in paths}


def main():
    before = hashes()
    started = time.perf_counter()
    result = subprocess.run([sys.executable, '-m', 'pytest', '-q'], cwd=ROOT,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output = result.stdout.decode('utf-8', errors='replace')
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT/'pytest_final.log').write_text(output, encoding='utf-8')
    after = hashes()
    record = dict(command='python -m pytest -q', exit_code=result.returncode,
                  runtime_seconds=time.perf_counter()-started,
                  python=platform.python_version(), host=platform.node(),
                  platform=platform.platform(), hashes_before=before,
                  hashes_after=after, files_unchanged=before==after, output=output)
    (OUT/'regression_result.json').write_text(json.dumps(record, indent=2)+'\n', encoding='utf-8')
    print(output)
    if before != after:
        raise RuntimeError('Source or tests changed during the regression.')
    raise SystemExit(result.returncode)


if __name__ == '__main__':
    main()
