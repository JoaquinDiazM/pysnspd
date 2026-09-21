"""Independent replay and closed-form moments of the explicit derived table."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time

import numpy as np

ROOT=Path(__file__).resolve().parents[4]
OUT=ROOT/'docs/implementation/stage1_closure'
sys.path.insert(0,str(ROOT))
from pysnspd.experimental.material_preprocessing import derive_phonon_shape


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-file',type=Path,default=ROOT/'tmp/stage1_closure/material_sources/nbn-a2f-ph.dat')
    args=parser.parse_args();started=time.perf_counter()
    measured=json.loads((OUT/'material_results.json').read_text())
    raw=np.loadtxt(args.source_file,comments='#')
    first={}
    for index,row in enumerate(raw):
        first.setdefault(float(row[0]),(index,row.copy()))
    expected=np.asarray([first[key][1] for key in sorted(first)])
    indices=np.array([first[key][0] for key in sorted(first)])
    bad=expected[:,2]<=0
    expected[bad,2]=0
    expected[bad|(expected[:,0]==0),1]=0
    saved=np.loadtxt(OUT/'material_nbn_shape_v1.csv',delimiter=',',skiprows=1)
    replay_error=float(np.max(abs(expected-saved[:,:3])))
    trace_equal=bool(np.array_equal(indices,saved[:,3]))
    x,a=expected[:,0],expected[:,1]
    dx=np.diff(x);lo=x[:-1];left=a[:-1];right=a[1:]
    slope=(right-left)/dx
    nonzero=lo>0
    terms=np.empty_like(dx)
    terms[~nonzero]=slope[~nonzero]*dx[~nonzero]
    logarithm=np.log1p(dx[nonzero]/lo[nonzero])
    terms[nonzero]=left[nonzero]*logarithm+slope[nonzero]*(dx[nonzero]-lo[nonzero]*logarithm)
    h=4.135667696 # only to compare the SAME explicitly conditional THz convention
    moments=dict(lambda_regular_origin=float(2*np.sum(terms)),
                 moment_E_power_0=float(h*np.sum(dx*(left+right)/2)),
                 moment_E_power_1=float(h*h*np.sum(dx*(lo*(left+right)/2+dx*(left+2*right)/6))),
                 moment_E_power_2=float(h**3*np.sum(dx*(lo*lo*(left+right)/2+lo*dx*(left+2*right)/3+dx*dx*(left+3*right)/12))))
    reference=measured['alpha_moments']['common_support_v1']
    moment_errors={key:abs(value-reference[key])/abs(reference[key]) for key,value in moments.items()}
    shape=derive_phonon_shape(args.source_file,source_url=measured['primary_sources']['sources']['nbn-a2f-ph.dat']['url'],
                             source_revision=measured['primary_sources']['revision'])
    q=(x[:-1]+x[1:])/2
    alpha,dos,ratio=shape.evaluate(q)
    identity=float(np.max(abs(alpha-dos*ratio)))
    finite=bool(np.all(np.isfinite(ratio)))
    no_units=all(shape.manifest[key] is None for key in ('certified_axis_unit','certified_DOS_unit',
                                                       'certified_normalization_basis','certified_density_m3'))
    passed=replay_error==0 and trace_equal and max(moment_errors.values())<1e-10 and identity<1e-12 and finite and no_units
    result=dict(schema='pysnspd.stage1_closure.independent_material_replay.v1',status='PASS' if passed else 'FAIL',
                interpretation='PASS verifies traceable numerical preprocessing only. It does not admit absolute units, normalization or global physical phonon kinetics.',
                source_sha256=sha(args.source_file),derived_csv_sha256=sha(OUT/'material_nbn_shape_v1.csv'),
                criteria_sha256=sha(OUT/'acceptance_criteria.json'),module_sha256=sha(ROOT/'pysnspd/experimental/material_preprocessing.py'),
                reviewer_source_sha256=sha(Path(__file__)),material_results_sha256=sha(OUT/'material_results.json'),
                summary=dict(replay_absolute_error=replay_error,source_row_trace_equal=trace_equal,
                             closed_form_moment_relative_errors=moment_errors,ratio_identity_abs_error=identity,
                             finite_derived_ratio=finite,absolute_metadata_remain_uncertified=no_units),
                conditional_closed_form_moments=moments,runtime_seconds=time.perf_counter()-started)
    (OUT/'review'/'material_replay.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({key:result[key] for key in ('status','summary','runtime_seconds')},indent=2))


if __name__=='__main__':
    main()
