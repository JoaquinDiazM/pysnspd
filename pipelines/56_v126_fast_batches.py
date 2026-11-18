"""Shared-budget v1.2.6 fast batches; production requires one user-run check."""
from pathlib import Path
import argparse,importlib.util,json,sys,hashlib
ROOT=Path(__file__).resolve().parents[1]
def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
S=load('_v126_fast_coordinator',ROOT/'pipelines/54_v126_fast_campaign.py')
D=load('_v126_fast_admission',ROOT/'pipelines/55_v126_shared_admission.py')
ENGINE=D.Admission(S,Path(__file__))
IMPLEMENTATION_FILES=('pipelines/53_v126_fast_photon.py','pipelines/54_v126_fast_campaign.py','pipelines/55_v126_shared_admission.py',
                      'pipelines/56_v126_fast_batches.py','pipelines/58_v126_verify_result.py','pysnspd/experimental/fast_latency.py','pysnspd/data/v126_fast_tail.json','pysnspd/analysis/fast_inventory.py')

def implementation_pins():return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in IMPLEMENTATION_FILES}

def main(argv=None):
    argv=list(sys.argv[1:] if argv is None else argv)
    early=argparse.ArgumentParser(add_help=False,allow_abbrev=False);early.add_argument('--batch-catalog',type=Path,default=S.CATALOG_PATH)
    selected,rest=early.parse_known_args(argv);S.set_catalog(selected.batch_catalog)
    if '--execute' in rest:
        batch=rest[rest.index('--batch')+1]
        if batch!='verify':
            receipt=Path(S.CATALOG['verification_receipt'])
            if not receipt.exists():raise RuntimeError('Run00_verify_v126.sh once before production; assistant must not launch it')
            proof=json.loads(receipt.read_text())
            if not proof.get('passed') or proof.get('implementation_sha256')!=implementation_pins():raise RuntimeError('User-run verification does not match this implementation')
    S.preview=ENGINE.preview;S.reservation=ENGINE.reservation
    return S.main(rest)

if __name__=='__main__':raise SystemExit(main())
