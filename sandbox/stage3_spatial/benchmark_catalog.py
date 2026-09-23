"""Bounded catalogue timing pilot; no transient or spatial admission claim."""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from pysnspd.experimental.energy_catalog import OccupationEnergyCatalog
from pysnspd.experimental.refined_cells import refined_count_catalog


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    args=parser.parse_args()
    output=Path(args.output)
    if output.exists(): raise SystemExit('Choose a fresh output; existing evidence is preserved.')
    started=time.perf_counter()
    source=ROOT/'docs/implementation/stage1_r2/catalogs/occupation_catalog.npz'
    base=OccupationEnergyCatalog.load(source)
    catalog=refined_count_catalog(base)
    rows=[]
    for index,(amplitude,gamma) in enumerate(((.9,0.),(.9,.005),(.93,.007)),1):
        before=time.perf_counter()
        print(f'[{index-1}/3] kernel delta={amplitude:g}, gamma={gamma:g}; elapsed={before-started:.2f}s',flush=True)
        energy,da,dg=catalog.energy_kernel(amplitude,gamma)
        rows.append(dict(amplitude=amplitude,gamma=gamma,seconds=time.perf_counter()-before,
                         finite=bool(np.isfinite([energy,da,dg]).all())))
        print(f'[{index}/3] complete; elapsed={time.perf_counter()-started:.2f}s',flush=True)
    data=dict(scope='Timing only; no spatial verdict',runtime_seconds=time.perf_counter()-started,
        electronic_nodes=len(catalog.count_nodes),host=platform.node(),python=platform.python_version(),
        numpy=np.__version__,catalog_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        amplitude_support=base.vacuum.delta_axis[[0,-1]].tolist(),
        gamma_support=base.vacuum.gamma_axis[[0,-1]].tolist(),
        delta0_J=base.vacuum.delta0_J,D_m2_s=base.vacuum.D_m2_s,
        N0_per_J_m3=base.vacuum.N0_per_J_m3,metadata=base.vacuum.metadata,rows=rows)
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:data[k] for k in ('scope','runtime_seconds','electronic_nodes','rows')}),flush=True)


if __name__=='__main__': main()
