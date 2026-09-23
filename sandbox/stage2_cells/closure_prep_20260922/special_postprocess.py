"""Complete only postprocessing of saved special trajectories; never integrate.

Requires the existing reviewer dependency via PYTHONPATH=tmp/stage1_r2_review/deps.
The original failed import record and every trajectory remain unchanged.
"""
from pathlib import Path
import argparse
import json
import sys
import time
import special_cases as special


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True)
    parser.add_argument('--input-root',type=Path,required=True)
    parser.add_argument('--output-root',type=Path,required=True)
    args=parser.parse_args();plan=json.loads(args.plan.read_text())
    special.check_sources(plan)
    previous=json.loads((args.input_root/'special_assessment.json').read_text())
    if len(previous['cases'])!=3 or any(row['status']!='PASS' for row in previous['cases']):
        raise ValueError('All three stored trajectories must already have passed their own checks')
    if args.output_root.exists():raise ValueError('Preserve prior postprocessing; use a fresh output root')
    paths=[];hashes={}
    for task in plan['tasks']:
        path=args.input_root/(task['id']+'.json')
        data=json.loads(path.read_text())
        if data['status']!='COMPLETED_NOT_YET_ADJUDICATED':raise ValueError('Incomplete trajectory')
        if special.sha(path.with_suffix('.npz'))!=data['trajectory_sha256']:
            raise ValueError('Trajectory hash mismatch')
        for source,digest in data['source_hashes'].items():
            if special.sha(special.ROOT/source)!=digest:raise ValueError('Source changed:'+source)
        for p in (path,path.with_suffix('.npz'),path.with_name(path.stem+'_initial.npz')):
            hashes[str(p)]=special.sha(p)
        paths.append(path)
    args.output_root.mkdir(parents=True)
    started=time.perf_counter();result=dict(schema='pysnspd.stage2.special-postprocessing.v1',
        status='INCOMPLETE',new_trajectories=0,physical_sources_changed=False,
        stage2_admission=False,original_failed_record_sha256=special.sha(args.input_root/'special_assessment.json'),
        plan_sha256=special.sha(args.plan),source_sha256=special.sha(__file__),
        inputs=hashes,postprocessing={})
    output=args.output_root/'special_postprocess.json'
    special.save(output,result)
    try:
        result['postprocessing']['equilibrium_fields']=special.equilibrium_fields(
            paths[0],args.output_root/'special_equilibrium_fields.json')
        special.save(output,result)
        for name in ('fields','support'):
            target=args.output_root/f'special_boundary_{name}.json'
            special.run_command([sys.executable,str(special.ROOT/f'sandbox/stage2_cells/check_trajectory_{name}.py'),
                *map(str,paths[1:]),'--output',str(target)],target.with_suffix('.console.log'))
            item=json.loads(target.read_text())
            result['postprocessing']['boundary_'+name]={k:item[k] for k in
                ('status','maxima','worst_relative_upper_bound','runtime_seconds') if k in item}
            special.save(output,result)
            if item['status']!='PASS':raise RuntimeError('Failed boundary postprocessor:'+name)
        special.check_sources(plan)
        if any(special.sha(Path(p))!=digest for p,digest in hashes.items()):
            raise ValueError('Input changed during postprocessing')
        result['status']='PASS_POSTPROCESSING_FINITE_GAMMA_SUPPORT_SEPARATE'
        result['remaining']=['Equilibrium finite-Gamma support is evaluated by a separate root-owned diagnostic.',
            'No temporal-refinement or global stage2 admission is inferred.']
    except Exception as exc:
        result['status']='FAILED_OR_INCOMPLETE';result['reason']=repr(exc)
        raise
    finally:
        result['runtime_seconds']=time.perf_counter()-started;special.save(output,result)
        print(json.dumps({k:result[k] for k in ('status','runtime_seconds','new_trajectories')}),flush=True)


if __name__=='__main__':main()
