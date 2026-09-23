"""Run the stage-3 nodal and historical tests with a bounded runtime and a saved receipt."""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[2]
TESTS=['tests/test_experimental_spatial_functional.py','tests/test_stage3_progress.py',
       'tests/test_stage3_runner.py', 'tests/test_experimental_spatial_nodal.py',
       'tests/test_stage3_nodal_runner.py']


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-root',required=True)
    args=parser.parse_args();out=Path(args.output_root)
    if out.exists():raise SystemExit('Choose a fresh output path.')
    out.mkdir(parents=True)
    command=[sys.executable,'-m','pytest','-q',*TESTS]
    started=time.perf_counter()
    print('Pruebas breves: funcional, progreso y reuso. Límite: 240 s.',flush=True)
    try:
        result=subprocess.run(command,cwd=ROOT,capture_output=True,text=True,timeout=240)
        status='PASS' if result.returncode==0 else 'FAIL'
        output=result.stdout+result.stderr;code=result.returncode
    except subprocess.TimeoutExpired as exc:
        status='INCOMPLETE_TIMEOUT';code=124
        output=''.join(v.decode(errors='replace') if isinstance(v,bytes) else v or '' for v in (exc.stdout,exc.stderr))
    (out/'pytest.log').write_text(output,encoding='utf-8')
    record=dict(status=status,exit_code=code,runtime_seconds=time.perf_counter()-started,
        command=command,tests={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in TESTS})
    (out/'result.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print(output,flush=True);print(json.dumps(record,indent=2),flush=True)
    raise SystemExit(code)


if __name__=='__main__':main()
