"""Run the full repository test suite and record source-bound evidence."""
from pathlib import Path
import hashlib
import json
import platform
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'docs/implementation/stage1_r2'
SOURCE = ROOT/'pysnspd/experimental/energy_catalog.py'


def main():
    digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    before = digest(SOURCE)
    started = time.perf_counter()
    result = subprocess.run([sys.executable, '-m', 'pytest', '-q'], cwd=ROOT,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    log = result.stdout.decode('utf-8', errors='replace')
    (OUT/'pytest_final.log').write_text(log, encoding='utf-8')
    record = dict(command='python -m pytest -q', exit_code=result.returncode,
                  runtime_seconds=time.perf_counter()-started, python=platform.python_version(),
                  platform=platform.platform(), source_sha256=before,
                  source_sha256_at_end=digest(SOURCE),
                  test_files_sha256={p.relative_to(ROOT).as_posix():digest(p)
                    for p in sorted((ROOT/'tests').glob('test_experimental_*.py'))},
                  output=log)
    (OUT/'regression_result.json').write_text(json.dumps(record, indent=2)+'\n', encoding='utf-8')
    print(log)
    if before != digest(SOURCE):
        raise RuntimeError('Source changed during verification.')
    raise SystemExit(result.returncode)


if __name__ == '__main__':
    main()
