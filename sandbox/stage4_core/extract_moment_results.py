"""Verify a completed moment-resolution run and retain compact review data.

No spectral/kinetic solves. All original maps are checked in place; only seven
integrated projection maps and metadata are copied into a fresh review folder.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(HERE))
from parallel_runtime import limit_thread_environment,linux_resources,resource_budget
limit_thread_environment()
import numpy as np


def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for piece in iter(lambda:stream.read(1024**2),b''):digest.update(piece)
    return digest.hexdigest()


def read(path):return json.loads(Path(path).read_text(encoding='utf8'))


def write(path,value):
    with Path(path).open('x',encoding='utf8') as stream:
        json.dump(value,stream,indent=2,allow_nan=False);stream.write('\n')


def check(job):
    phase,path,row=job
    actual=sha(path)
    if actual!=row['fields_sha256']:raise ValueError('Changed raw map: '+str(path))
    if read(path.with_suffix('.json'))!=row:raise ValueError('Map metadata disagrees with phase summary')
    with np.load(path) as a:
        if any(np.any(~np.isfinite(a[key])) for key in a.files):raise ValueError('Nonfinite raw map')
        extra={}
        if phase=='spectra':
            extra=dict(maximum_normalization_residual=float(np.max(abs(a['g']**2+a['f']*a['f_tilde']-1))),
                minimum_DOS=float(np.min(a['g'].real)))
            if extra['maximum_normalization_residual']>1e-7 or extra['minimum_DOS']< -1e-8:
                raise ValueError('Final frozen spectrum fails normalized causal chart')
    return dict(phase=phase,id=row['id'],path=str(path),sha256=actual,bytes=path.stat().st_size,
        arrays_finite=True,metadata_matches_summary=True,**extra)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();started=time.monotonic()
    if args.output.exists():raise FileExistsError('Fresh extraction output required')
    identity,summary,plan=(read(args.raw/name) for name in ('identity.json','summary.json','executed_plan.json'))
    if summary['status']!='MOMENT_RESOLUTION_CAMPAIGN_COMPLETE':raise ValueError('Campaign not complete')
    if sha(args.raw/'executed_plan.json')!=identity['plan_sha256']:raise ValueError('Plan changed')
    for name,expected in {**identity['sources'],**identity['inputs']}.items():
        if sha(ROOT/name)!=expected:raise ValueError('Frozen source/input changed: '+name)
    jobs=[];phases={}
    for phase,count in (('spectra',181),('kinetic',181),('projection',7)):
        result=read(args.raw/phase/'summary.json')
        if result['status']!='PHASE_COMPLETE' or len(result['records'])!=count:raise ValueError('Incomplete phase')
        if len({r['id'] for r in result['records']})!=count:raise ValueError('Duplicate record')
        phases[phase]=result
        jobs.extend((phase,args.raw/phase/r['fields_path'],r) for r in result['records'])
    if {r['id'] for r in phases['spectra']['records']}!={r['id'] for r in phases['kinetic']['records']}:
        raise ValueError('Spectral and kinetic queries differ')
    budget=resource_budget(linux_resources(),.9,.9,max_workers=4)
    original=os.sched_getaffinity(0)
    try:
        os.sched_setaffinity(0,set(budget['worker_affinity_cpus']+[budget['coordinator_cpu']]))
        with ThreadPoolExecutor(max_workers=budget['workers']) as pool:
            records=list(pool.map(check,jobs))
    finally:os.sched_setaffinity(0,original)
    args.output.mkdir(parents=True)
    for name in ('identity.json','summary.json','executed_plan.json','moment_comparisons.json','progress.jsonl'):
        shutil.copyfile(args.raw/name,args.output/name)
    copied={}
    for phase,result in phases.items():
        (args.output/phase).mkdir()
        for name in ('summary.json','resource_snapshot.json'):
            shutil.copyfile(args.raw/phase/name,args.output/phase/name)
        if phase=='projection':
            (args.output/phase/'fields').mkdir()
            for row in result['records']:
                path=Path(row['fields_path'])
                for suffix in ('.npz','.json'):
                    relative=Path(phase)/path.with_suffix(suffix)
                    shutil.copyfile(args.raw/relative,args.output/relative)
    for path in args.output.rglob('*'):
        if path.is_file():copied[path.relative_to(args.output).as_posix()]=sha(path)
    receipt=dict(status='ALL_RAW_MAPS_VERIFIED_NO_SOLVES',raw_directory=str(args.raw.resolve()),
        runtime_seconds=time.monotonic()-started,script_sha256=sha(__file__),budget=budget,
        sources_matching_identity=identity['sources'],inputs_matching_identity=identity['inputs'],
        maps_verified=records,copied_file_sha256=copied,physical_time_steps=0,new_spectral_solves=0,
        raw_spectral_and_kinetic_arrays='Remain on Geminga; all 362 field files checked by SHA-256, finite arrays and matching metadata',
        compact_projection_maps=7)
    write(args.output/'extraction_receipt.json',receipt)
    print(json.dumps(dict(status=receipt['status'],maps=len(records),raw_bytes=sum(r['bytes'] for r in records),
        runtime_seconds=receipt['runtime_seconds'])))


if __name__=='__main__':main()
