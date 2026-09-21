"""Audit the untouched artifacts recovered from the user's GNU Screen session."""
from pathlib import Path
import hashlib,json,sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'sandbox/stage2_cells'))
from assess_time_refinement import load_run

DATA=ROOT/'docs/implementation/stage2/resume_20260921'
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    folder=DATA/'manual_reference'
    path=folder/'one_reference_probe.json'
    run=load_run(path);r=run['record']
    mismatch=[p for p,s in r['source_hashes'].items() if sha(ROOT/p)!=s]
    assert not mismatch,mismatch
    assert sha(ROOT/'docs/implementation/stage2/acceptance_criteria.json')==r['criteria_sha256']
    assert sha(ROOT/'docs/implementation/stage1_r2/catalogs/occupation_catalog.npz')==r['catalog_sha256']
    progress=[json.loads(line) for line in (folder/'one_reference_probe.progress.jsonl').read_text().splitlines()]
    assert progress[0]['event']=='START' and progress[-1]['event']=='STOP'
    artifacts=[dict(path=p.relative_to(ROOT).as_posix(),sha256=sha(p),bytes=p.stat().st_size)
               for p in sorted(folder.iterdir()) if p.is_file()]
    report=dict(schema='pysnspd.stage2.user_run_import.v1',status='VERIFIED_COMPLETED_REFERENCE',
        stage2_closed=False,source='User-run foreground command in detached screen3040582.code_000',
        remote_directory='/home/jdiaz/pysnspd/tmp/stage2_user_20260921_193252',
        retrieval='screen -S3040582.code_000 -Xhardcopy -h; existing files copied without changes',
        runtime_seconds=r['runtime_seconds'],rhs_calls=r['rhs_calls'],
        energy_ledger_scaled_max=r['energy_ledger_scaled_max'],
        instantaneous_residual_max=r['instantaneous_residual_max'],
        saved_times=run['times'].tolist(),source_mismatches=mismatch,
        amplitude_change=float(r['final']['amplitudes'][0]-r['initial']['amplitudes'][0]),
        eph_energy_fraction=abs(r['final']['electron_to_phonon'][0])/r['initial_excitation_energy'],
        scope='Local t=0.1 reference with only3 output points; not full activity, field/support or grid acceptance',
        artifacts=artifacts,script_sha256=sha(Path(__file__)))
    (DATA/'manual_reference_import.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='artifacts'},indent=2))

if __name__=='__main__':main()
