"""Check the user-executed fast trajectory against an existing normal reference."""
from pathlib import Path
import argparse,importlib.util,json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
def compare_preparation(fast,prior):
    if set(fast)!=set(prior):raise ValueError('Photon preparation arrays differ')
    phase=np.exp(1j*np.angle(np.vdot(fast['delta_bar'],prior['delta_bar'])))
    for k in fast:
        value=fast[k]*phase if k=='delta_bar' else fast[k]
        scale=max(float(np.max(np.abs(prior[k]))),1e-30)
        np.testing.assert_allclose(value,prior[k],rtol=1e-10,atol=1e-12*scale,err_msg='Different original physical preparation: '+k)

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--batch-catalog',type=Path,required=True);p.add_argument('--campaign-dir',type=Path,required=True)
    a=p.parse_args(argv);catalog=json.loads(a.batch_catalog.read_text());case=catalog['batches']['verify'][0]
    from pysnspd.analysis.fast_inventory import validate_run
    run=a.campaign_dir/'runs'/case['name'];result=validate_run(run,ROOT)
    if not result['estimated']:raise ValueError('Verification must exercise the early-HB estimate path')
    normal=Path(catalog['verification_normal_reference']);m=json.loads((normal/'manifest.json').read_text())
    observed=m['observation']['latency_s']*1e12;error=abs(result['tau_d_ps']-observed)
    if error>1.:raise ValueError(f'Existing-normal delay mismatch {error:.6g}ps exceeds1ps')
    with np.load(run/'source_initial_state.npz',allow_pickle=False) as fast,np.load(normal/'source_initial_state.npz',allow_pickle=False) as prior:
        compare_preparation({k:fast[k] for k in fast.files},{k:prior[k] for k in prior.files})
    s=importlib.util.spec_from_file_location('_v126_verified_driver',ROOT/'pipelines/56_v126_fast_batches.py');d=importlib.util.module_from_spec(s);s.loader.exec_module(d)
    receipt=dict(passed=True,execution_mode='fast',model_version='1.2.6',normal_reference=str(normal),fast_run=str(run),
                 normal_delay_ps=observed,estimated_delay_ps=result['tau_d_ps'],absolute_error_ps=error,implementation_sha256=d.implementation_pins(),
                 validation_scope='One user-executed same-fast checkpoint roundtrip and early-HB completion; not a full prospective scientific validation')
    target=Path(catalog['verification_receipt']);target.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2));return 0
if __name__=='__main__':raise SystemExit(main())
