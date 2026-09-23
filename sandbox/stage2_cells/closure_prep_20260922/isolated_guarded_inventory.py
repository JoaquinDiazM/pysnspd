"""External environment and complete artifact inventory; no physical calls."""
from datetime import datetime, timezone
from pathlib import Path
import getpass
import hashlib
import json
import platform
import sys
import numpy as np
import scipy

ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'docs/implementation/stage2/closure_prep_20260922'
OUTPUT=BASE/'isolated_guarded_inventory.json'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    folders=['isolated_guarded','isolated_guarded_negative','isolated_guarded_negative_v2']
    rows=[]
    for folder in folders:
        if not (BASE/folder).is_dir():raise FileNotFoundError(folder)
        for path in sorted((BASE/folder).rglob('*')):
            if path.is_file():
                rows.append(dict(path=path.relative_to(ROOT).as_posix(),sha256=sha(path),bytes=path.stat().st_size))
    result=dict(schema='pysnspd.stage2.guarded-external-inventory.v1',
        recorded_at_utc=datetime.now(timezone.utc).isoformat(),
        purpose='Append environment provenance and every artifact hash without rewriting frozen QA results.',
        environment_observed_after_QA=dict(host=platform.node(),account=getpass.getuser(),
            python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
            executable=sys.executable,platform=platform.platform(),
            longdouble_bits=int(np.finfo(np.longdouble).bits),
            longdouble_mantissa_bits=int(np.finfo(np.longdouble).nmant),
            longdouble_min_exponent=int(np.finfo(np.longdouble).minexp)),
        roots=folders,artifact_count=len(rows),artifacts=rows,
        numerical_sources_sha256={path:sha(ROOT/path) for path in
            ['sandbox/stage2_cells/closure_prep_20260922/limited_ssp_guarded.py',
             'sandbox/stage2_cells/recovery_20260921/limited_ssp.py',
             'sandbox/stage2_cells/closure_prep_20260922/isolated_guarded_validation.py',
             'sandbox/stage2_cells/closure_prep_20260922/isolated_guarded_negative.py',
             'sandbox/stage2_cells/closure_prep_20260922/isolated_guarded_negative_v2.py']},
        inventory_script_sha256=sha(__file__),new_RHS_calls=0,new_integrations=0,
        limitation='This separate environment observation does not retroactively alter the original QA record. Original result and archive bytes remain unchanged.')
    with OUTPUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(dict(artifact_count=len(rows),environment=result['environment_observed_after_QA'],output=OUTPUT.as_posix())))


if __name__=='__main__':main()
