"""Read-only inventory of user-produced evidence; no trajectory is executed."""
from pathlib import Path
import hashlib,json
ROOT=Path('/home/jdiaz/pysnspd')
def main():
    mapping={'rk4':ROOT/'tmp/stage2_validation_20260922_015109',
             'ssp':ROOT/'tmp/stage2_recovery_20260921_ssp'}
    rows=[]
    for group,directory in mapping.items():
        for path in sorted(directory.iterdir()):
            if not path.is_file():raise ValueError('Unexpected nested directory')
            rows.append(dict(relative=group+'/'+path.name,source=str(path),bytes=path.stat().st_size,
                sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    for path in (ROOT/'docs/implementation/stage2/resume_20260921/selected_continuum.json',
                 ROOT/'tmp/rhs_smoothness_recheck.json',ROOT/'tmp/stage2_recovery_20260921_ssp.screen.log'):
        rows.append(dict(relative='diagnostics/'+path.name,source=str(path),bytes=path.stat().st_size,
            sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    print(json.dumps(dict(schema='pysnspd.stage2.user_results_inventory.v1',files=rows),indent=2))
if __name__=='__main__':main()
