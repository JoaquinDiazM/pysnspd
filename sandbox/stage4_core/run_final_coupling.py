"""Manual long-run handoff: stationary biased references, then coupled response.

Every stage is a foreground computation with durable logs, progress and ETA.
Independent cases share one <=90% resource budget. Failed references suppress
only their own dependent responses. Nothing is retried or overwritten.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from sandbox.stage4_core.coupled_campaign import run_campaign, atomic_json, sha, allocate, validate_plan
from sandbox.stage4_core.parallel_runtime import linux_resources, limit_thread_environment


def campaign(cases):
    return dict(schema='pysnspd.stage4.coupled_campaign.v1', maximum_parallel_cases=2,
        maximum_workers_per_case=13, heartbeat_seconds=5., cases=cases)


def reference_cases(plan):
    return [dict(id=item['id'], argv=['{python}', '-u', 'sandbox/stage4_core/biased_strip_reference.py',
        '--plan', item['plan'], '--output', '{output}', '--workers', '{workers}', '--execute'],
        expected_outputs=['reference.npz', 'manifest.json'], estimated_seconds=item['estimated_seconds'])
        for item in plan['references']]


def response_cases(plan, output, admitted):
    return [dict(id=item['id'], argv=['{python}', '-u', 'sandbox/stage4_core/biased_coupled_response.py',
        '--plan', item['plan'], '--reference', str(output/'references'/'cases'/item['reference']/'reference.npz'),
        '--output', '{output}', '--workers', '{workers}'],
        expected_outputs=['results.json', 'coupled_fields.npz', 'integrated_moments.npz', 'manifest.json'],
        estimated_seconds=item['estimated_seconds']) for item in plan['responses'] if item['reference'] in admitted]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--output-root', type=Path, required=True)
    parser.add_argument('--execute', action='store_true')
    args = parser.parse_args()
    limit_thread_environment()
    plan = json.loads(args.plan.read_text(encoding='utf8'))
    if plan.get('schema') != 'pysnspd.stage4.final_coupling.v1':
        raise ValueError('Unsupported workflow plan')
    for name, expected in plan.get('source_sha256', {}).items():
        if sha(ROOT/name) != expected:
            raise ValueError('Pinned source changed: '+name)
    for item in plan['references']+plan['responses']:
        if sha(ROOT/item['plan']) != item['sha256']:
            raise ValueError('Input plan changed: '+item['plan'])
    references = campaign(reference_cases(plan))
    validate_plan(references)
    reference_ids = {item['id'] for item in plan['references']}
    if any(item['reference'] not in reference_ids for item in plan['responses']):
        raise ValueError('Every response must name a declared reference')
    validate_plan(campaign(response_cases(plan,args.output_root.resolve(),reference_ids)))
    if not args.execute:
        print(json.dumps(dict(status='DRY_RUN', output_exists=args.output_root.exists(),
            references=references, responses=response_cases(plan, args.output_root, {x['id'] for x in plan['references']}),
            allocation=allocate(references, linux_resources()), note='No scientific calculation executed'), indent=2))
        return 0
    if args.output_root.exists():
        raise FileExistsError('A fresh output root is required; preserve previous output')
    args.output_root.mkdir(parents=True)
    atomic_json(args.output_root/'workflow_plan.json', plan)
    first = run_campaign(references, args.output_root/'references', resources=linux_resources())
    admitted = {row['id'] for row in first['cases'] if row['status'] == 'SUCCEEDED'}
    responses = response_cases(plan, args.output_root.resolve(), admitted)
    blocked = [row['id'] for row in plan['responses'] if row['reference'] not in admitted]
    atomic_json(args.output_root/'dependencies.json', dict(admitted_references=sorted(admitted), blocked_responses=blocked))
    if not responses:
        atomic_json(args.output_root/'workflow_result.json', dict(status='REFERENCE_STAGE_INCOMPLETE', blocked=blocked))
        return 2
    second = run_campaign(campaign(responses), args.output_root/'responses', resources=linux_resources())
    result = dict(status='COMPLETED' if second['status'] == 'SUCCEEDED' and not blocked else 'COMPLETED_WITH_INCOMPLETE_CASES',
        reference_status=first['status'], response_status=second['status'], blocked=blocked,
        meaning='Execution status only; compare physical observables and residuals before admitting the closure.')
    atomic_json(args.output_root/'workflow_result.json', result)
    print(json.dumps(result), flush=True)
    return 0 if result['status'] == 'COMPLETED' else 2


if __name__ == '__main__':
    raise SystemExit(main())
