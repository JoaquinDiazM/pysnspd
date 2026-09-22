"""Manual foreground diagnostic of one-cell SSP temporal convergence.

Default is a dry run. No calculation is launched during preparation. The user
executes the new frozen plan; a failure stops the sequence without any retry,
clipping, evidence overwrite, scheduler, screen session or background launch.
"""
from pathlib import Path
import argparse
from datetime import datetime, timezone
import importlib.metadata
import importlib.util
import json
import os
import platform
import re
import sys
import time

ROOT=Path(__file__).resolve().parents[3]
ORIGINAL='sandbox/stage2_cells/resume_20260921/run_acceptance_batch.py'
ADAPTER='sandbox/stage2_cells/recovery_20260921/run_limited_coupled.py'
INTEGRATOR='sandbox/stage2_cells/recovery_20260921/limited_ssp.py'
METHOD='ssprk3_common_flux_limited'
ASSESSORS={
    'limited_time':'sandbox/stage2_cells/recovery_20260921/assess_limited_time.py',
    'limited_pair':'sandbox/stage2_cells/recovery_20260921/assess_limited_pair.py',
    'grid':'sandbox/stage2_cells/assess_grid_refinement.py',
    'fields':'sandbox/stage2_cells/check_trajectory_fields.py',
    'support':'sandbox/stage2_cells/check_trajectory_support.py',
}
spec=importlib.util.spec_from_file_location('stage2_original_batch',ROOT/ORIGINAL)
old=importlib.util.module_from_spec(spec)
sys.modules[spec.name]=old
spec.loader.exec_module(old)
sha,write_new,workspace_path=old.sha,old.write_new,old.workspace_path
# Reuse the already tested foreground streaming/interrupt behavior unchanged.
stream_command=old.stream_command


def dependencies(task):
    return task.get('runs',[])+([task['reference']] if 'reference' in task else [])


def required_contract_paths(plan):
    paths={ORIGINAL,ADAPTER,INTEGRATOR,old.RUNNER,old.CRITERIA,old.CATALOG,
           'sandbox/stage2_cells/assess_time_refinement.py',
           Path(__file__).relative_to(ROOT).as_posix(),*ASSESSORS.values(),
           *[p.relative_to(ROOT).as_posix() for p in (ROOT/'pysnspd/experimental').glob('*.py')]}
    for task in plan['tasks']:
        if task['kind']=='existing_reference':
            path=workspace_path(task['path'])
            paths.update(p.relative_to(ROOT).as_posix() for p in (
                path,path.with_suffix('.npz'),path.with_name(path.stem+'_initial.npz')))
    return sorted(paths)


def validate_plan(plan,executing=False):
    if plan.get('schema')!='pysnspd.stage2.one-cell-time-plan.v1':
        raise ValueError('Unsupported recovery-plan schema')
    if executing and plan.get('status')!='READY_FOR_USER_FOREGROUND_RUN':
        raise ValueError('Only a reviewed and source-frozen plan may execute')
    seen=set()
    for task in plan.get('tasks',[]):
        identifier=task.get('id','');kind=task.get('kind')
        if not re.fullmatch(r'[a-z][a-z0-9_]*',identifier) or identifier in seen:
            raise ValueError('Unique lowercase task identifier required')
        if kind in ('trajectory','existing_reference'):
            parameters=task.get('parameters',{})
            if set(parameters)!=old.PARAMETERS:
                raise ValueError('Every frozen runner parameter must be explicit: '+identifier)
            if parameters['method']!=(METHOD if kind=='trajectory' else 'rk4'):
                raise ValueError('Unexpected integration-method contract: '+identifier)
            if parameters['case'] not in ('one','two') or parameters['scenario'] not in (
                    'driven','equilibrium','phonon_vacuum','sparse_electrons'):
                raise ValueError('Invalid physical case or scenario')
            for key in ('steps','phonon_nodes','face_order','reaction_order','reaction_outer_order','electron_refinement'):
                if type(parameters[key]) is not int or parameters[key]<1:
                    raise ValueError('Positive integer required: '+key)
            if any(isinstance(value,(dict,list,bool)) or value is None for value in parameters.values()):
                raise ValueError('Explicit scalar trajectory parameters required')
            if kind=='existing_reference':workspace_path(task['path'])
        elif kind in ASSESSORS:
            deps=dependencies(task)
            if not deps or any(name not in seen for name in deps):
                raise ValueError('Assessment dependencies must precede it: '+identifier)
            if kind in ('limited_time','grid') and 'reference' not in task:
                raise ValueError('An explicit reference is required')
            if kind=='limited_time' and len(task.get('runs',[]))<3:
                raise ValueError('Three temporal resolutions are required')
            if kind=='limited_pair' and (len(task.get('runs',[]))!=2 or 'reference' in task):
                raise ValueError('The supplemental diagnostic takes exactly two trajectories')
            if kind=='grid' and task.get('parameter') not in ('electron_refinement','phonon_nodes'):
                raise ValueError('Explicit grid parameter required')
            expected=['PASS_PAIR_DIAGNOSTIC'] if kind=='limited_pair' else ['PASS']
            if task.get('accepted_statuses')!=expected:
                raise ValueError('Assessment status must distinguish pair diagnostics from admission')
        else:raise ValueError('Unsupported task kind: '+str(kind))
        seen.add(identifier)
    if not seen:raise ValueError('Empty plan')
    contracts=plan.get('source_hashes',{})
    if executing or contracts:
        if not set(required_contract_paths(plan)).issubset(contracts):
            raise ValueError('Freeze every kernel, adapter, assessor and existing reference archive before execution')
        for name,digest in contracts.items():
            if sha(workspace_path(name))!=digest:raise ValueError('Frozen contract changed: '+name)
    return contracts


def command_for(task,output,completed,executable=sys.executable):
    if task['kind']=='existing_reference':return None
    if task['kind']=='trajectory':
        command=[executable,'-u',str(ROOT/ADAPTER)]
        for name,value in task['parameters'].items():
            command+=['--'+name.replace('_','-'),'rk4' if name=='method' else str(value)]
    else:
        command=[executable,'-u',str(ROOT/ASSESSORS[task['kind']]),
                 *[str(completed[name]) for name in task.get('runs',[])]]
        if 'reference' in task:command+=['--reference',str(completed[task['reference']])]
        if task['kind']=='grid':command+=['--parameter',task['parameter']]
    return command+['--output',str(output)]


def validate_trajectory(path,parameters,contracts):
    # Import only during execution; dry-run never imports the new integrator.
    sys.path.insert(0,str(ROOT/'sandbox/stage2_cells'))
    from assess_grid_refinement import read
    from run_coupled import setup
    import numpy as np
    loaded=read(path,ROOT);record=loaded['record']
    if {k:v for k,v in record['parameters'].items() if k!='output'}!=parameters:
        raise ValueError('Completed trajectory differs from the explicit parameter contract')
    sources={name:digest for name,digest in contracts.items()
             if name==old.RUNNER or name.startswith('pysnspd/experimental/')}
    if parameters['method']==METHOD:
        sources.update({name:contracts[name] for name in (ADAPTER,INTEGRATOR)})
        integration=record.get('time_integration',{})
        if integration.get('population_clipping') is not False or integration.get('posthoc_energy_repair') is not False:
            raise ValueError('Limited integration must explicitly exclude population clipping and energy repair')
    elif 'time_integration' in record:
        raise ValueError('Unrecognized reference integrator adapter')
    if record['source_hashes']!=sources:raise ValueError('Physical/integration source contract differs')
    keys=('case','phonon_nodes','infrared','face_order','escape','heating','scenario',
          'reaction_order','reaction_method','electron_refinement','reaction_layout',
          'reaction_outer_order','reaction_max_panel','reaction_max_energy_panel')
    setup_parameters={key:parameters[key] for key in keys}
    if setup_parameters['escape']=='inf':setup_parameters['escape']=float('inf')
    system,initial,_=setup(**setup_parameters)
    if not np.array_equal(loaded['arrays']['states'][0],initial):
        raise ValueError('Initial populations do not match the frozen physical setup')
    for name,expected in (('electron_count',system.catalog.count_nodes),
        ('electron_weights',system.catalog.count_weights),('phonon_energies',system.phonons.energies),
        ('phonon_capacities',system.phonons.capacities)):
        if not np.array_equal(loaded['arrays'][name],expected):raise ValueError('Grid mismatch: '+name)
    versions=dict(python=platform.python_version(),numpy=importlib.metadata.version('numpy'),
                  scipy=importlib.metadata.version('scipy'))
    if any(record.get(name)!=value for name,value in versions.items()):
        raise ValueError('Dependency versions differ from the reference run')
    limits=json.loads((ROOT/old.CRITERIA).read_text())['coupled_trajectories']
    if (record['energy_ledger_scaled_max']>limits['energy_ledger_scaled_max'] or
            record['instantaneous_residual_max']>limits['instantaneous_balance_scaled_max'] or
            record['minimum_electron']<0 or record['maximum_electron']>1 or record['minimum_phonon']<0):
        raise ValueError('Completed trajectory fails invariant/physical-population gates')
    return {str(p):sha(p) for p in (path,path.with_suffix('.npz'),path.with_name(path.stem+'_initial.npz'))}


def self_test():
    """Synthetic orchestration tests: no physical module or subprocess runs."""
    import contextlib
    import copy
    import io
    import tempfile
    from unittest.mock import patch
    started=time.perf_counter();checks={};module=sys.modules[__name__]
    plan=json.loads((ROOT/'docs/implementation/stage2/review_20260922/manual_time_plan.json').read_text())
    validate_plan(plan,False)
    task=next(item for item in plan['tasks'] if item['kind']=='trajectory')
    original_parameters=copy.deepcopy(task['parameters'])
    command=command_for(task,ROOT/'tmp/dummy_output.json',{})
    checks['entry_argument_rk4_but_registered_method_remains_explicit_SSP']=(
        command[command.index('--method')+1]=='rk4' and task['parameters']==original_parameters
        and original_parameters['method']==METHOD)
    draft=copy.deepcopy(plan);draft['status']='DRAFT_PENDING_SOURCE_FREEZE'
    try:validate_plan(draft,True)
    except ValueError:checks['unreviewed_unfrozen_plan_cannot_execute']=True
    else:checks['unreviewed_unfrozen_plan_cannot_execute']=False
    invalid=copy.deepcopy(plan)
    assessment=next(item for item in invalid['tasks'] if item['kind']=='limited_time')
    assessment['accepted_statuses']=['PASS_PAIR_DIAGNOSTIC']
    try:validate_plan(invalid,False)
    except ValueError:checks['three_level_assessment_cannot_accept_only_a_pair_diagnostic']=True
    else:checks['three_level_assessment_cannot_accept_only_a_pair_diagnostic']=False
    temporary_parent=ROOT/'tmp';temporary_parent.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='limited_batch_guards_',dir=temporary_parent) as folder:
        testroot=Path(folder)
        dry_output=testroot/'dry_output'
        args=['batch','--plan',str(ROOT/'docs/implementation/stage2/review_20260922/manual_time_plan.json'),
              '--output-root',str(dry_output),'--dry-run']
        before_integrator='limited_ssp' in sys.modules
        with patch.object(sys,'argv',args),patch.object(module,'stream_command') as child, \
                contextlib.redirect_stdout(io.StringIO()):
            main()
        checks['dry_run_has_no_child_output_directory_or_integrator_import']=(child.call_count==0
            and not dry_output.exists() and ('limited_ssp' in sys.modules)==before_integrator)
        fixture=testroot/'reference.json';write_new(fixture,{'status':'fixture'})
        first=dict(id='first',phase='fixture',kind='existing_reference',path=str(fixture),parameters={})
        fresh=dict(id='fresh',phase='fixture',kind='trajectory',parameters={})
        second=dict(id='second',phase='fixture',kind='trajectory',parameters={})
        def valid(path,parameters,contracts):return {str(path):sha(path)}
        def invoke(name,tasks,stream,validator=valid):
            plan_path=testroot/(name+'_plan.json');output=testroot/name
            if not plan_path.exists():write_new(plan_path,dict(tasks=tasks))
            argv=['batch','--plan',str(plan_path),'--output-root',str(output),'--execute']
            with patch.object(sys,'argv',argv),patch.object(module,'validate_plan',return_value={}), \
                 patch.object(module,'validate_trajectory',side_effect=validator), \
                 patch.object(module,'stream_command',side_effect=stream) as child, \
                 patch.object(importlib.metadata,'version',return_value='synthetic'), \
                 contextlib.redirect_stdout(io.StringIO()):
                try:main();failure=None
                except (ValueError,RuntimeError):failure=sys.exc_info()[1]
            return output,child.call_count,failure
        def never_launch(*args):raise AssertionError('No real child is permitted')
        output,calls,failure=invoke('reuse',[first],never_launch)
        receipt=(output/'first.receipt.json').read_bytes()
        _,calls_again,failure_again=invoke('reuse',[first],never_launch)
        checks['reference_reuse_verified_twice_without_reexecution_or_overwrite']=(
            failure is None and failure_again is None and calls==calls_again==0
            and (output/'first.receipt.json').read_bytes()==receipt)
        def invalid_reference(*args):raise ValueError('incomplete reference fixture')
        _,calls,failure=invoke('invalid_reference',[first,second],never_launch,invalid_reference)
        checks['invalid_reference_stops_before_any_child']=calls==0 and isinstance(failure,ValueError)
        def failed_child(*args):raise RuntimeError('synthetic child failure')
        output,calls,failure=invoke('child_fail',[fresh,second],failed_child)
        checks['failed_child_stops_without_retry_or_next_task']=(calls==1 and isinstance(failure,RuntimeError)
            and not (output/'second.json').exists() and 'BATCH_STOPPED' in (output/'batch_events.jsonl').read_text())
        assessment=dict(id='check',phase='fixture',kind='fields',runs=['first'],accepted_statuses=['PASS'])
        def failed_assessment(command,log):write_new(Path(command[-1]),dict(status='FAIL'))
        output,calls,failure=invoke('assessment_fail',[first,assessment,second],failed_assessment)
        checks['FAIL_even_with_zero_child_exit_stops_before_next_task']=(calls==1
            and isinstance(failure,ValueError) and not (output/'second.json').exists())
        checks['successful_offline_batch_does_not_issue_stage_admission']=(
            '"stage2_admission": false' in (testroot/'reuse'/'batch_events.jsonl').read_text())
    return dict(schema='pysnspd.stage2.manual-time-offline-selftest.v1',
        status='PASS' if all(checks.values()) else 'FAIL',checks=checks,
        scope='Synthetic guard fixtures only; no subprocess, physical RHS or time integration.',
        runner_sha256=sha(__file__),elapsed_seconds=time.perf_counter()-started)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path)
    parser.add_argument('--output-root',type=Path)
    parser.add_argument('--self-test',action='store_true')
    parser.add_argument('--self-test-output',type=Path)
    modes=parser.add_mutually_exclusive_group()
    modes.add_argument('--dry-run',action='store_true')
    modes.add_argument('--execute',action='store_true')
    args=parser.parse_args()
    if args.self_test:
        result=self_test()
        if args.self_test_output:write_new(workspace_path(args.self_test_output),result)
        print(json.dumps(result,indent=2))
        if result['status']!='PASS':raise SystemExit(1)
        return
    if args.plan is None or args.output_root is None:parser.error('--plan and --output-root are required')
    plan_path=workspace_path(args.plan);output_root=workspace_path(args.output_root)
    plan=json.loads(plan_path.read_text(encoding='utf-8'));contracts=validate_plan(plan,args.execute)
    completed={}
    if not args.execute:
        print('DRY RUN: no child process, RHS, new-integrator import, output directory or evidence modification.')
        print('Plan status:',plan.get('status'))
        for task in plan['tasks']:
            output=output_root/(task['id']+'.json')
            path=workspace_path(task['path']) if task['kind']=='existing_reference' else output
            print(json.dumps(dict(id=task['id'],phase=task['phase'],kind=task['kind'],
                existing_reference=str(path) if task['kind']=='existing_reference' else None,
                command=command_for(task,output,completed)),ensure_ascii=False))
            completed[task['id']]=path
        return
    for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[name]='1'
    signature=dict(schema='pysnspd.stage2.one-cell-time-manifest.v1',plan_sha256=sha(plan_path),
        source_hashes=contracts,python=platform.python_version(),
        numpy=importlib.metadata.version('numpy'),scipy=importlib.metadata.version('scipy'))
    manifest=output_root/'batch_manifest.json'
    if output_root.exists():
        if not manifest.exists() or json.loads(manifest.read_text())!=signature:
            raise ValueError('Output directory has another or missing contract')
    else:output_root.mkdir(parents=True);write_new(manifest,signature)
    started=time.perf_counter()
    with (output_root/'batch_events.jsonl').open('a',encoding='utf-8',buffering=1) as log:
        def emit(event,**data):
            row=dict(event=event,utc=datetime.now(timezone.utc).isoformat(),
                     elapsed_seconds=time.perf_counter()-started,**data)
            line=json.dumps(row,ensure_ascii=False);print(line,flush=True);log.write(line+'\n')
        emit('BATCH_START',tasks=len(plan['tasks']),plan=str(plan_path),contract_sha256=sha(manifest))
        try:
            for index,task in enumerate(plan['tasks']):
                validate_plan(plan,True)
                identifier=task['id'];output=output_root/(identifier+'.json')
                chosen=workspace_path(task['path']) if task['kind']=='existing_reference' else output
                receipt=output_root/(identifier+'.receipt.json')
                command=command_for(task,output,completed)
                contract=dict(task=task,command=command,batch_manifest_sha256=sha(manifest),
                    dependency_hashes={name:sha(completed[name]) for name in dependencies(task)})
                emit('TASK_START',index=index+1,total=len(plan['tasks']),task=identifier,
                     phase=task['phase'],output=str(chosen),command=command)
                saved=None
                if receipt.exists():
                    saved=json.loads(receipt.read_text())
                    if saved['contract']!=contract:raise ValueError('Existing receipt contract differs: '+identifier)
                    if any(sha(Path(name))!=digest for name,digest in saved['output_hashes'].items()):
                        raise ValueError('Previously completed output changed')
                if task['kind']=='existing_reference' or (task['kind']=='trajectory' and chosen.exists()):
                    hashes=validate_trajectory(chosen,task['parameters'],contracts)
                    emit('REUSED_COMPLETE_TRAJECTORY',task=identifier,hashes=hashes)
                elif saved:
                    if json.loads(output.read_text()).get('status') not in task['accepted_statuses']:
                        raise ValueError('Reused assessment is not an accepted result')
                    hashes=saved['output_hashes'];emit('REUSED_PASS_ASSESSMENT',task=identifier)
                else:
                    related=[output,output.with_suffix('.npz'),output.with_suffix('.integration.jsonl'),
                        output.with_name(output.stem+'_initial.npz'),output.with_suffix('.console.log')]
                    if any(path.exists() for path in related):
                        raise ValueError('Existing/unfinished evidence is preserved; no automatic retry')
                    stream_command(command,output.with_suffix('.console.log'))
                    if task['kind']=='trajectory':hashes=validate_trajectory(output,task['parameters'],contracts)
                    else:
                        result=json.loads(output.read_text())
                        if result.get('status') not in task['accepted_statuses']:
                            raise ValueError('Assessment did not pass: '+identifier)
                        hashes={str(output):sha(output)}
                if not receipt.exists():write_new(receipt,dict(contract=contract,output_hashes=hashes))
                completed[identifier]=chosen;emit('TASK_COMPLETE',task=identifier)
            emit('BATCH_COMPLETE',stage2_admission=False,
                 scope='The listed numerical checks only; final independent review remains required.')
        except BaseException as exc:
            emit('BATCH_STOPPED',exception=type(exc).__name__,reason=str(exc),
                 policy='No retry, fallback, clipping, overwrite or continuation after failure.')
            raise


if __name__=='__main__':main()
