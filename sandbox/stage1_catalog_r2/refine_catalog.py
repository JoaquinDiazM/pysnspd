"""Insert explicitly diagnosed ratio nodes, retaining the existing catalogue."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
import numpy as np
from pysnspd.experimental.energy_catalog import OccupationEnergyCatalog,refine_occupation_catalog


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('seed',type=Path)
    parser.add_argument('--ratios',type=float,nargs='+',required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--report',type=Path,required=True)
    args=parser.parse_args()
    started=time.perf_counter()
    seed=OccupationEnergyCatalog.load(args.seed)
    seed_hash=hashlib.sha256(args.seed.read_bytes()).hexdigest()
    table=refine_occupation_catalog(seed,np.unique(np.r_[seed.gamma_ratio_axis,args.ratios]),reference_hash=seed_hash)
    table.save(args.output)
    report=dict(schema='pysnspd.stage1_r2.incremental_refinement.v1',runtime_seconds=time.perf_counter()-started,
                seed_sha256=seed_hash,catalog_sha256=hashlib.sha256(args.output.read_bytes()).hexdigest(),
                source_sha256=hashlib.sha256((ROOT/'pysnspd/experimental/energy_catalog.py').read_bytes()).hexdigest(),
                inserted_ratio_nodes=args.ratios,shape=list(table.excitation_energies.shape),
                admission_status='PENDING_INDEPENDENT_ASSESSMENT',production_connected=False)
    args.report.parent.mkdir(parents=True,exist_ok=True)
    args.report.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
