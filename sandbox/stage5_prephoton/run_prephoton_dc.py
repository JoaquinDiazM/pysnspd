"""Explicit user-run campaign: intrinsic bulk references and full-circuit dark holds."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from sandbox.stage4_core.parallel_runtime import limit_thread_environment,linux_resources
limit_thread_environment()
from sandbox.stage4_core.coupled_campaign import run_campaign,atomic_json,sha,allocate,validate_plan
from sandbox.stage4_core.dual_kwt_time import explicit_step_bound
from sandbox.stage5_prephoton.dc_hold import load_reference
from sandbox.stage5_prephoton.provenance import source_sha
import numpy as np

SCHEMA='pysnspd.stage5.prephoton_dc_campaign.v1'

def campaign(plan,cases):
    return dict(schema='pysnspd.stage4.coupled_campaign.v1',maximum_parallel_cases=plan['maximum_parallel_cases'],
        maximum_workers_per_case=plan['maximum_workers_per_case'],heartbeat_seconds=5.,cases=cases)

def reference_cases(plan):
    return [dict(id=x['id'],argv=['{python}','-u','sandbox/stage5_prephoton/prepare_dc_reference.py',
        '--plan',x['plan'],'--output','{output}','--workers','{workers}'],
        expected_outputs=['reference.npz','summary.json','manifest.json'],estimated_seconds=x['estimated_seconds'])
        for x in plan['references']]

def validate_workflow(plan):
    if plan.get('schema')!=SCHEMA:raise ValueError('Unsupported DC workflow')
    for path,digest in plan['source_sha256'].items():
        if source_sha(ROOT/path)!=digest:raise ValueError('Pinned source changed: '+path)
    for row in plan['references']:
        if sha(ROOT/row['plan'])!=row['sha256']:raise ValueError('Reference plan changed: '+row['id'])
    ids={x['id'] for x in plan['references']}
    if any(x['reference'] not in ids or not 0<x['step_multiplier']<=1 for x in plan['holds']):
        raise ValueError('Each hold needs a declared reference and valid time-step multiplier')
    if len({x['id'] for x in plan['holds']})!=len(plan['holds']):raise ValueError('Duplicate hold id')
    validate_plan(campaign(plan,reference_cases(plan)))

def hold_plan(plan,item,reference):
    """Only numerical step and omitted inductance are derived; physics is unchanged."""
    summary=json.loads((reference.parent/'summary.json').read_text())
    if not summary['admitted_for_dark_hold'] or sha(reference)!=summary['reference_sha256']:
        raise ValueError('Hold requires an admitted, unchanged DC reference')
    graph,gap,epsilon,_=load_reference(reference)
    result=json.loads(json.dumps(plan['hold_template']))
    # The old uniform Schur bound needs a real amplitude. Do not feed a q*x
    # complex boundary value into its amplitude-only formula.
    uniform=np.full(graph.n_nodes,summary['bulk']['gap_bulk_bar'],complex)
    bound=explicit_step_bound(graph,uniform,epsilon,dict(result,euler_safety=plan['euler_safety']))
    result['dt_ps']=item['step_multiplier']*min(plan['maximum_dt_ps'],bound['primary_step_ps'])
    result['bias_current_A']=summary['bulk']['current_target_A']
    result['circuit']['Lk_ext_H']=summary['fixed_external_inductance_H']
    result['reference_sha256']=sha(reference)
    result['source_sha256']=plan['source_sha256']
    result['source_hash_normalization']=plan['source_hash_normalization']
    result['step_estimate']=dict(bound,scope='Uniform symbol used only to choose a conservative initial step; biased-state acceptance requires observed stationarity and the explicit half-step comparison.')
    result['observation_cadence']=max(1,int(np.ceil(result['duration_ps']/result['dt_ps']/100)))
    return result

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True)
    parser.add_argument('--output-root',type=Path,required=True)
    parser.add_argument('--execute',action='store_true')
    args=parser.parse_args();plan=json.loads(args.plan.read_text());validate_workflow(plan)
    references=campaign(plan,reference_cases(plan));resources=linux_resources()
    if not args.execute:
        print(json.dumps(dict(status='DRY_RUN',references=references,holds=plan['holds'],
            allocation=allocate(references,resources),output_exists=args.output_root.exists(),
            note='No numerical calculation executed. Hold plans derive Lext and step from admitted references.'),indent=2))
        return 0
    output=args.output_root.resolve()
    if output.exists():raise FileExistsError('Use a fresh output root; no overwrite or implicit retry')
    output.mkdir(parents=True);atomic_json(output/'workflow_plan.json',plan)
    first=run_campaign(references,output/'references',resources=resources)
    admitted={x['id'] for x in first['cases'] if x['status']=='SUCCEEDED'}
    directory=output/'hold_plans';directory.mkdir()
    cases=[];blocked=[]
    for item in plan['holds']:
        if item['reference'] not in admitted:
            blocked.append(item['id']);continue
        ref=output/'references/cases'/item['reference']/'reference.npz'
        path=directory/(item['id']+'.json');atomic_json(path,hold_plan(plan,item,ref))
        cases.append(dict(id=item['id'],argv=['{python}','-u','sandbox/stage5_prephoton/dc_hold.py',
            '--plan',str(path),'--reference',str(ref),'--output','{output}','--workers','{workers}'],
            expected_outputs=['summary.json','history.json','fields.npz','manifest.json'],estimated_seconds=item['estimated_seconds']))
    atomic_json(output/'dependencies.json',dict(admitted_references=sorted(admitted),blocked_holds=blocked))
    second=run_campaign(campaign(plan,cases),output/'holds',resources=linux_resources()) if cases else None
    from sandbox.stage5_prephoton.analyze_prephoton_dc import analyze
    result=analyze(output)
    result['reference_execution_status']=first['status']
    result['hold_execution_status']=second['status'] if second else 'BLOCKED'
    atomic_json(output/'workflow_result.json',result);print(json.dumps(result),flush=True)
    return 0 if result['prephoton_dc_admitted'] else 2

if __name__=='__main__':raise SystemExit(main())
