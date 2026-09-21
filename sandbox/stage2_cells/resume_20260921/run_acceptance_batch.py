"""Foreground, user-run stage-2 batch; dry-run is the default.

This orchestrates the frozen runner. It never retries, clips populations,
overwrites evidence, launches background work, or changes the numerical model.
An existing complete run is reusable only after checking its full contract.
"""
from pathlib import Path
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import os
import platform
import re
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
RUNNER = 'sandbox/stage2_cells/run_coupled.py'
OBSERVER = 'sandbox/stage2_cells/run_reference_progress.py'
CRITERIA = 'docs/implementation/stage2/acceptance_criteria.json'
CATALOG = 'docs/implementation/stage1_r2/catalogs/occupation_catalog.npz'
PARAMETERS = {'case', 'steps', 'duration', 'phonon_nodes', 'infrared',
    'face_order', 'reaction_order', 'reaction_layout',
    'reaction_max_energy_panel', 'reaction_outer_order', 'reaction_max_panel',
    'electron_refinement', 'reaction_method', 'escape', 'heating', 'method',
    'scenario', 'rtol'}
ASSESSORS = {
    'time': 'sandbox/stage2_cells/assess_time_refinement.py',
    'grid': 'sandbox/stage2_cells/assess_grid_refinement.py',
    'fields': 'sandbox/stage2_cells/check_trajectory_fields.py',
    'support': 'sandbox/stage2_cells/check_trajectory_support.py',
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_new(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x', encoding='utf-8') as stream:
        stream.write(json.dumps(data, indent=2, allow_nan=False)+'\n')


def workspace_path(name):
    path = Path(name)
    path = (ROOT/path).resolve() if not path.is_absolute() else path.resolve()
    if not path.is_relative_to(ROOT.resolve()):
        raise ValueError('Path escapes the declared workspace: '+str(name))
    return path


def validate_plan(plan, executing=False):
    if plan.get('schema') != 'pysnspd.stage2.manual-batch-plan.v1':
        raise ValueError('Unsupported batch-plan schema')
    if executing and plan.get('status') != 'READY_FOR_MANUAL_RUN':
        raise ValueError('Draft plan cannot execute; unresolved numerical choices must be preregistered first.')
    contracts = plan.get('source_hashes', {})
    expected = {RUNNER, OBSERVER, CRITERIA, CATALOG, Path(__file__).relative_to(ROOT).as_posix(),
                *ASSESSORS.values(),
                *[p.relative_to(ROOT).as_posix() for p in (ROOT/'pysnspd/experimental').glob('*.py')]}
    if not expected.issubset(contracts):
        raise ValueError('The plan must freeze all kernels, runner, observer, assessors, criteria and catalogue.')
    for name, value in contracts.items():
        if sha(workspace_path(name)) != value:
            raise ValueError('Frozen source/input hash mismatch: '+name)
    seen = set()
    for task in plan.get('tasks', []):
        name = task.get('id', '')
        if not re.fullmatch(r'[a-z][a-z0-9_]*', name) or name in seen:
            raise ValueError('Task ids must be unique lowercase identifiers')
        kind = task.get('kind')
        if kind == 'trajectory':
            parameters = task.get('parameters', {})
            if set(parameters) != PARAMETERS:
                raise ValueError('Declare every runner parameter exactly once: '+name)
            if parameters['case'] not in ('one', 'two') or parameters['method'] not in ('rk4', 'dop853'):
                raise ValueError('Invalid case or method')
            if parameters['scenario'] not in ('driven', 'equilibrium', 'phonon_vacuum', 'sparse_electrons'):
                raise ValueError('Invalid scenario')
            for key in ('steps', 'phonon_nodes', 'face_order', 'reaction_order', 'reaction_outer_order'):
                if type(parameters[key]) is not int or parameters[key] < 1:
                    raise ValueError('Positive integer required for '+key)
            if type(parameters['electron_refinement']) is not int or parameters['electron_refinement'] < 0:
                raise ValueError('Nonnegative electron refinement required')
            for value in parameters.values():
                if isinstance(value, (dict, list, bool)) or value is None:
                    raise ValueError('Only scalar explicit trajectory parameters are accepted')
            if 'reuse_from' in task:
                workspace_path(task['reuse_from'])
        elif kind in ASSESSORS:
            dependencies = task.get('runs', []) + ([task['reference']] if 'reference' in task else [])
            if not dependencies or any(x not in seen for x in dependencies):
                raise ValueError('Assessment dependencies must precede the task: '+name)
            if kind in ('time', 'grid') and 'reference' not in task:
                raise ValueError('Time/grid assessment needs an explicit reference')
            if kind == 'grid' and task.get('parameter') not in ('electron_refinement', 'phonon_nodes'):
                raise ValueError('Declare the grid parameter')
        else:
            raise ValueError('Unsupported task kind: '+str(kind))
        seen.add(name)
    if not seen:
        raise ValueError('An empty plan cannot run')
    return contracts


def command_for(task, output, completed, executable=sys.executable):
    if task['kind'] == 'trajectory':
        command = [executable, '-u', str(ROOT/OBSERVER)]
        for name, value in task['parameters'].items():
            command += ['--'+name.replace('_', '-'), str(value)]
    else:
        command = [executable, '-u', str(ROOT/ASSESSORS[task['kind']]),
                   *[str(completed[name]) for name in task.get('runs', [])]]
        if 'reference' in task:
            command += ['--reference', str(completed[task['reference']])]
        if task['kind'] == 'grid':
            command += ['--parameter', task['parameter']]
    return command+['--output', str(output)]


def validate_trajectory(path, parameters, contracts):
    # The existing independent reader checks grids, all stored states, NPZ
    # hashes and their relation to snapshots. It performs no RHS evaluation.
    sys.path.insert(0, str(ROOT/'sandbox/stage2_cells'))
    from assess_grid_refinement import read
    from run_coupled import setup
    import numpy as np
    loaded = read(path, ROOT)
    record = loaded['record']
    actual = {key: value for key, value in record['parameters'].items() if key != 'output'}
    if actual != parameters:
        raise ValueError('Existing trajectory parameters differ: '+str(path))
    sources = {key: value for key, value in contracts.items()
               if key == RUNNER or key.startswith('pysnspd/experimental/')}
    if record.get('source_hashes') != sources:
        raise ValueError('Existing trajectory has a different source contract')
    setup_keys = ('case', 'phonon_nodes', 'infrared', 'face_order', 'escape', 'heating',
        'scenario', 'reaction_order', 'reaction_method', 'electron_refinement',
        'reaction_layout', 'reaction_outer_order', 'reaction_max_panel', 'reaction_max_energy_panel')
    setup_parameters = {key: parameters[key] for key in setup_keys}
    if setup_parameters['escape'] == 'inf':
        setup_parameters['escape'] = float('inf')
    system, expected_initial, _ = setup(**setup_parameters)
    if not np.array_equal(loaded['arrays']['states'][0], expected_initial):
        raise ValueError('Archived initial state differs from the frozen setup')
    for name, expected_grid in [('electron_count', system.catalog.count_nodes),
            ('electron_weights', system.catalog.count_weights),
            ('phonon_energies', system.phonons.energies),
            ('phonon_capacities', system.phonons.capacities)]:
        if not np.array_equal(loaded['arrays'][name], expected_grid):
            raise ValueError('Archived grid differs from the frozen setup: '+name)
    versions = {'python': platform.python_version(),
                'numpy': importlib.metadata.version('numpy'),
                'scipy': importlib.metadata.version('scipy')}
    if any(record.get(key) != value for key, value in versions.items()):
        raise ValueError('Existing trajectory was computed with different dependency versions')
    limits = json.loads((ROOT/CRITERIA).read_text())['coupled_trajectories']
    if (record['energy_ledger_scaled_max'] > limits['energy_ledger_scaled_max']
            or record['instantaneous_residual_max'] > limits['instantaneous_balance_scaled_max']
            or record['minimum_electron'] < 0 or record['maximum_electron'] > 1
            or record['minimum_phonon'] < 0):
        raise ValueError('Completed trajectory fails conservation/support gates; batch stopped.')
    return {str(p): sha(p) for p in [path, path.with_suffix('.npz'),
            path.with_name(path.stem+'_initial.npz')]}


def stream_command(command, log):
    environment = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    with log.open('x', encoding='utf-8', buffering=1) as stream:
        process = subprocess.Popen(command, cwd=ROOT, env=environment,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
        try:
            for line in process.stdout:
                print(line, end='', flush=True)
                stream.write(line)
            code = process.wait()
        except BaseException:
            process.terminate()
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
            raise
    if code:
        raise RuntimeError('Child command failed with exit code '+str(code)+'; no retry performed.')


def self_test():
    """Offline orchestration/guard checks; no RHS or child process is run."""
    import contextlib
    import io
    import tempfile
    from unittest.mock import patch
    started = time.perf_counter()
    checks = {}
    module = sys.modules[__name__]
    temporary_parent = (ROOT/'tmp').resolve()
    temporary_parent.mkdir(exist_ok=True)
    if not temporary_parent.is_relative_to(ROOT.resolve()):
        raise ValueError('Test directory must remain inside the workspace')
    with tempfile.TemporaryDirectory(prefix='batch_guard_', dir=temporary_parent) as folder:
        testroot = Path(folder)
        marker = testroot/'preserved.json'
        write_new(marker, {'original': True})
        original = marker.read_bytes()
        try:
            write_new(marker, {'replacement': True})
        except FileExistsError:
            checks['existing_evidence_never_overwritten'] = marker.read_bytes() == original

        # Exercise the real independent record reader at its early guards.
        for name, record in [
            ('incomplete', {'schema':'pysnspd.stage2.coupled_run.v1', 'status':'INCOMPLETE'}),
            ('source_mismatch', {'schema':'pysnspd.stage2.coupled_run.v1',
                'status':'COMPLETED_NOT_YET_ADJUDICATED', 'source_hashes':{RUNNER:'0'*64}})]:
            path = testroot/(name+'.json')
            write_new(path, record)
            try:
                validate_trajectory(path, {}, {})
            except ValueError as exc:
                checks[name+'_rejected_without_NPZ_or_RHS'] = (
                    'completed' in str(exc) if name == 'incomplete' else 'source hash mismatch' in str(exc))

        external = testroot/'existing_complete.json'
        write_new(external, {'status':'fixture-complete'})
        first = dict(id='first', kind='trajectory', parameters={}, reuse_from=str(external))
        def invoke(name, tasks, stream, validator):
            plan_path, output = testroot/(name+'_plan.json'), testroot/name
            if not plan_path.exists():
                write_new(plan_path, {'tasks':tasks})
            argv = ['batch', '--plan', str(plan_path), '--output-root', str(output), '--execute']
            with patch.object(sys, 'argv', argv), patch.object(module, 'validate_plan', return_value={}), \
                 patch.object(module, 'validate_trajectory', side_effect=validator), \
                 patch.object(module, 'stream_command', side_effect=stream) as child, \
                 patch.object(importlib.metadata, 'version', return_value='offline-fixture'), \
                 contextlib.redirect_stdout(io.StringIO()):
                try:
                    main()
                    failure = None
                except (ValueError, RuntimeError) as exc:
                    failure = exc
            return output, child.call_count, failure
        def valid_fixture(path, parameters, contracts):
            return {str(path):sha(path)}
        def never_launch(*args):
            raise AssertionError('No child process is permitted by this guard test')
        output, calls, failure = invoke('reuse', [first], never_launch, valid_fixture)
        receipt = (output/'first.receipt.json').read_bytes()
        output, calls_again, failure_again = invoke('reuse', [first], never_launch, valid_fixture)
        checks['complete_result_reused_twice_without_reexecution'] = (
            failure is None and failure_again is None and calls == calls_again == 0
            and (output/'first.receipt.json').read_bytes() == receipt)

        def invalid_fixture(path, parameters, contracts):
            raise ValueError('INCOMPLETE fixture')
        _, calls, failure = invoke('invalid_existing', [first], never_launch, invalid_fixture)
        checks['invalid_existing_result_stops_without_replacement_run'] = calls == 0 and isinstance(failure, ValueError)

        fresh = dict(id='first', kind='trajectory', parameters={})
        second = dict(id='second', kind='trajectory', parameters={})
        def child_failure(command, log):
            raise RuntimeError('Simulated child failure')
        output, calls, failure = invoke('failed_child', [fresh, second], child_failure, valid_fixture)
        checks['child_failure_stops_before_next_task_without_retry'] = (
            calls == 1 and isinstance(failure, RuntimeError)
            and not (output/'second.json').exists()
            and 'BATCH_STOPPED' in (output/'batch_events.jsonl').read_text())

        assessment = dict(id='assessment', kind='fields', runs=['first'])
        def failed_assessment(command, log):
            write_new(Path(command[-1]), {'status':'FAIL'})
        output, calls, failure = invoke('failed_assessment', [first, assessment, second], failed_assessment, valid_fixture)
        checks['FAIL_assessment_with_exit_zero_stops_before_next_task'] = (
            calls == 1 and isinstance(failure, ValueError) and not (output/'second.json').exists())
    return dict(schema='pysnspd.stage2.manual-batch-self-test.v1',
        status='PASS' if len(checks)==7 and all(checks.values()) else 'FAIL',
        scope='Offline synthetic guards only; no subprocess, time integration or RHS evaluation.',
        checks=checks, runner_sha256=sha(__file__), elapsed_seconds=time.perf_counter()-started)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path)
    parser.add_argument('--output-root', type=Path)
    parser.add_argument('--self-test', action='store_true')
    parser.add_argument('--self-test-output', type=Path)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--dry-run', action='store_true')
    group.add_argument('--execute', action='store_true', help='Manual foreground execution; can exceed five minutes.')
    args = parser.parse_args()
    if args.self_test:
        result = self_test()
        if args.self_test_output:
            write_new(workspace_path(args.self_test_output), result)
        print(json.dumps(result, indent=2))
        if result['status'] != 'PASS':
            raise SystemExit(1)
        return
    if args.plan is None or args.output_root is None:
        parser.error('--plan and --output-root are required outside --self-test')
    plan_path = workspace_path(args.plan)
    output_root = workspace_path(args.output_root)
    plan = json.loads(plan_path.read_text(encoding='utf-8'))
    contracts = validate_plan(plan, executing=args.execute)
    completed = {}
    if not args.execute:
        print('DRY RUN — no subprocess, RHS, output directory or evidence modification.', flush=True)
        print('Plan status:', plan.get('status'))
        for task in plan['tasks']:
            output = output_root/(task['id']+'.json')
            path = workspace_path(task['reuse_from']) if 'reuse_from' in task else output
            print(json.dumps(dict(id=task['id'], phase=task.get('phase'), kind=task['kind'],
                reuse_candidate=str(path) if 'reuse_from' in task else None,
                command=command_for(task, output, completed)), ensure_ascii=False))
            completed[task['id']] = path
        return
    for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        os.environ[name] = '1'
    signature = dict(schema='pysnspd.stage2.manual-batch-manifest.v1',
        plan_sha256=sha(plan_path), source_hashes=contracts,
        python=platform.python_version(), numpy=importlib.metadata.version('numpy'),
        scipy=importlib.metadata.version('scipy'))
    manifest = output_root/'batch_manifest.json'
    if output_root.exists():
        if not manifest.exists() or json.loads(manifest.read_text()) != signature:
            raise ValueError('Existing output directory has another or missing batch contract.')
    else:
        output_root.mkdir(parents=True)
        write_new(manifest, signature)
    events = output_root/'batch_events.jsonl'
    started = time.perf_counter()
    with events.open('a', encoding='utf-8', buffering=1) as log:
        def emit(event, **data):
            row = dict(event=event, utc=datetime.now(timezone.utc).isoformat(),
                       elapsed_seconds=time.perf_counter()-started, **data)
            line = json.dumps(row, ensure_ascii=False)
            print(line, flush=True); log.write(line+'\n')
        emit('BATCH_START', plan=str(plan_path), contract_sha256=sha(manifest), tasks=len(plan['tasks']))
        try:
            for index, task in enumerate(plan['tasks']):
                validate_plan(plan, executing=True)  # Abort if a frozen file changed during the batch.
                identifier = task['id']
                output = output_root/(identifier+'.json')
                receipt = output_root/(identifier+'.receipt.json')
                chosen = workspace_path(task['reuse_from']) if 'reuse_from' in task else output
                command = command_for(task, output, completed)
                dependencies = task.get('runs', []) + ([task['reference']] if 'reference' in task else [])
                dep_hashes = {name: sha(completed[name]) for name in dependencies}
                task_contract = dict(task=task, command=command, dependency_hashes=dep_hashes,
                                     batch_manifest_sha256=sha(manifest))
                emit('TASK_START', index=index+1, total=len(plan['tasks']), task=identifier,
                     phase=task.get('phase'), output=str(chosen))
                if receipt.exists():
                    saved = json.loads(receipt.read_text())
                    if saved['contract'] != task_contract:
                        raise ValueError('Existing task receipt contract differs: '+identifier)
                    for name, digest in saved['output_hashes'].items():
                        if sha(Path(name)) != digest:
                            raise ValueError('Previously completed output changed: '+name)
                if task['kind'] == 'trajectory' and chosen.exists():
                    hashes = validate_trajectory(chosen, task['parameters'], contracts)
                    emit('REUSED_COMPLETE_TRAJECTORY', task=identifier, hashes=hashes)
                elif receipt.exists():
                    if json.loads(output.read_text()).get('status') != 'PASS':
                        raise ValueError('A reused assessment is not PASS')
                    hashes = saved['output_hashes']
                    emit('REUSED_PASS_ASSESSMENT', task=identifier, hashes=hashes)
                else:
                    if 'reuse_from' in task:
                        raise ValueError('Explicit reuse candidate does not exist; no automatic replacement run.')
                    related = [output, output.with_suffix('.npz'), output.with_suffix('.progress.jsonl'),
                               output.with_name(output.stem+'_initial.npz'), output.with_suffix('.console.log')]
                    if any(path.exists() for path in related):
                        raise ValueError('Unfinished/existing task evidence is preserved; choose an amended plan and new output directory.')
                    stream_command(command, output.with_suffix('.console.log'))
                    if task['kind'] == 'trajectory':
                        hashes = validate_trajectory(output, task['parameters'], contracts)
                    else:
                        record = json.loads(output.read_text())
                        if record.get('status') != 'PASS':
                            raise ValueError('Assessment did not pass: '+identifier)
                        hashes = {str(output): sha(output)}
                    chosen = output
                if not receipt.exists():
                    write_new(receipt, dict(contract=task_contract, output_hashes=hashes))
                completed[identifier] = chosen
                emit('TASK_COMPLETE', task=identifier)
            emit('BATCH_COMPLETE', scope='Only the listed assessments; final stage admission still requires independent review.')
        except BaseException as exc:
            emit('BATCH_STOPPED', exception=type(exc).__name__, reason=str(exc),
                 policy='No retry, fallback, overwrite or continuation after a failed task.')
            raise


if __name__ == '__main__':
    main()
