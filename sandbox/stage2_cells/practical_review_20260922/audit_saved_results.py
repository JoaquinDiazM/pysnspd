"""Read-only provenance and descriptive sensitivity audit; never calls an RHS.

Writes new review artifacts only. Frozen raw evidence, sources and criteria are
not changed. Two-grid sensitivity and cell-relative descriptive norms are not
substitutes for the registered convergence/admission tests.
"""
from pathlib import Path
import csv
import hashlib
import json
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / 'sandbox/stage2_cells'),
               str(ROOT / 'sandbox/stage2_cells/closure_prep_20260922')]
import assess_grid_refinement as grid
import assess_time_refinement as temporal
import assess_guarded_time as guarded

OUT = ROOT / 'docs/implementation/stage2/practical_review_20260922'
RAW = OUT / 'raw'
PLAN = ROOT / 'docs/implementation/stage2/closure_prep_20260922/guarded_acceptance_plan.json'
FIELDS = ('amplitudes', 'excitation_energy', 'electron_energy', 'phonon_energy',
          'escape', 'input', 'electron_to_phonon', 'condensate_heat', 'transport')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def saved_path(remote):
    return RAW / Path(remote).name


def write_csv(path, rows):
    if path.exists():
        raise FileExistsError(path)
    with path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def metric(comparison, family, cell, t, actual, reference, population=False, capacity=1):
    a, b = np.asarray(actual), np.asarray(reference)
    absolute = float(np.sum(capacity*abs(a-b))) if population else float(np.max(abs(a-b)))
    scale = float(np.sum(capacity*abs(b))) if population else float(np.max(abs(b)))
    measured = temporal.measurement_error(absolute, scale)
    return dict(comparison=comparison, family=family, cell=cell, time=float(t),
                absolute=absolute, reference_norm=scale,
                relative_percent=None if measured['relative'] is None else 100*measured['relative'],
                original_norm_score=measured['assessment_error'],
                comparison_defined=measured['comparison_defined'],
                roundoff_limited=measured['roundoff_limited'])


def reduced(rows):
    groups = {}
    for row in rows:
        groups.setdefault((row['comparison'], row['family'], row['cell']), []).append(row)
    result = []
    for entries in groups.values():
        valid = [r for r in entries if r['original_norm_score'] is not None]
        worst = max(valid, key=lambda r: r['original_norm_score']) if valid else entries[0]
        result.append({**worst, 'undefined_checkpoints': len(entries)-len(valid),
                       'final_relative_percent': entries[-1]['relative_percent'],
                       'initial_relative_percent': entries[0]['relative_percent']})
    return result


def blocks(run):
    record, arrays = run['record'], run['arrays']
    cells = len(record['initial']['amplitudes'])
    ne, np_ = len(arrays['electron_count']), len(arrays['phonon_energies'])
    block = 1+ne+np_
    return arrays['states'][:, :cells*block].reshape(-1, cells, block), cells, ne


def compare_runs(actual, reference, name, mesh=False):
    rows = []
    a, ac, ne = blocks(actual)
    b, bc, nr = blocks(reference)
    assert ac == bc and ne == nr
    ta, tb = actual['arrays']['times'], reference['arrays']['times']
    names = grid.OBSERVABLES if mesh else FIELDS
    if mesh:
        ap = {k:v for k,v in actual['record']['parameters'].items() if k not in ('phonon_nodes','output')}
        bp = {k:v for k,v in reference['record']['parameters'].items() if k not in ('phonon_nodes','output')}
        assert ap == bp
        assert actual['record']['source_hashes'] == reference['record']['source_hashes']
        for key in ('catalog_sha256','criteria_sha256','material','python','numpy','scipy'):
            assert actual['record'][key] == reference['record'][key]
        for key in ('electron_count','electron_weights'):
            assert np.array_equal(actual['arrays'][key], reference['arrays'][key])
    else:
        guarded.compare(temporal.load_run(actual['path']), temporal.load_run(reference['path']))
    for i,t in enumerate(ta):
        j = int(np.argmin(abs(tb-t)))
        assert abs(tb[j]-t) <= 32*np.finfo(float).eps*max(1., abs(tb[-1]))
        sa, sb = actual['record']['snapshots'][i], reference['record']['snapshots'][j]
        for family in names:
            rows.append(metric(name, family, 'all', t, sa[family], sb[family]))
            if np.ndim(sa[family]):
                for c in range(ac):
                    rows.append(metric(name, family, str(c), t, sa[family][c], sb[family][c]))
        for family, low, high, capacity_key, factor in (
                ('electronic_population',1,1+ne,'electron_weights',4),
                ('phonon_population',1+ne,None,'phonon_capacities',1)):
            if mesh and family == 'phonon_population':
                x,y = actual['phonon_projection'][i],reference['phonon_projection'][j]
                family='phonon_65_bands'
                capacity=1
            else:
                capacity = factor*actual['arrays'][capacity_key]
                x,y = a[i,:,low:high],b[j,:,low:high]
            rows.append(metric(name, family, 'all', t, x, y, True, capacity))
            for c in range(ac):
                rows.append(metric(name, family, str(c), t, x[c], y[c], True, capacity))
        if mesh:
            for c in ['all', *map(str,range(ac))]:
                x,y=actual['electron_projection'][i],reference['electron_projection'][j]
                if c != 'all':x,y=x[int(c)],y[int(c)]
                rows.append(metric(name,'electron_129_bands',c,t,x,y,True))
    return rows


def phonon_difference(actual, reference, time_value):
    a,cells,ne=blocks(actual);b,_,_=blocks(reference)
    i=int(np.argmin(abs(actual['arrays']['times']-time_value)))
    j=int(np.argmin(abs(reference['arrays']['times']-time_value)))
    om=actual['arrays']['phonon_energies'];cap=actual['arrays']['phonon_capacities']
    assert np.array_equal(om,reference['arrays']['phonon_energies'])
    x,y=a[i,:,1+ne:],b[j,:,1+ne:]
    difference=cap*abs(x-y)
    mass=float(difference.sum())
    bands=[0,.005,.01,.025,.05,.1,.25,.5,1,2,4.0000000001]
    rows=[]
    for c in range(cells):
        for low,high in zip(bands[:-1],bands[1:]):
            keep=(om>=low)&(om<high)
            d=float(difference[c,keep].sum())
            rows.append(dict(time=float(actual['arrays']['times'][i]),cell=c,
                omega_low=low,omega_high=min(4.,high),mass_l1_difference=d,
                percent_of_global_mass_l1_difference=100*d/mass,
                signed_mass_difference=float(np.sum(cap[keep]*(x[c,keep]-y[c,keep]))),
                energy_weighted_l1_difference=float(np.sum(om[keep]*difference[c,keep])),
                reference_phonon_count=float(np.sum(cap[keep]*y[c,keep]))))
    summaries=[]
    for c in range(cells):
        arg=int(np.argmax(difference[c]));energy_ref=float(np.sum(om*cap*y[c]))
        summaries.append(dict(cell=c,mass_l1_difference=float(difference[c].sum()),
            percent_of_global_mass_l1_difference=100*float(difference[c].sum())/mass,
            reference_phonon_count=float(np.sum(cap*y[c])),
            relative_count_L1_percent=100*float(difference[c].sum())/float(np.sum(cap*y[c])),
            energy_weighted_L1_percent=100*float(np.sum(om*difference[c]))/energy_ref,
            signed_energy_difference=float(np.sum(om*cap*(x[c]-y[c]))),
            largest_mass_error_node=int(arg),largest_mass_error_omega=float(om[arg]),
            largest_mass_error_actual_occupation=float(x[c,arg]),
            largest_mass_error_reference_occupation=float(y[c,arg]),
            largest_node_share_of_global_mass_L1_percent=100*float(difference[c,arg])/mass))
    return rows,summaries


def main():
    started=time.perf_counter()
    outputs=['saved_results_audit.json','trajectory_audit.csv','comparison_metrics.csv',
             'comparison_checkpoints.csv','worst_phonon_error_bands.csv']
    if any((OUT/name).exists() for name in outputs):
        raise FileExistsError('Review evidence already exists; use a new review version')
    raw_hashes={p.name:sha(p) for p in sorted(RAW.iterdir()) if p.is_file()}
    manifest=read_json(RAW/'batch_manifest.json');plan=read_json(PLAN)
    assert sha(PLAN)==manifest['plan_sha256']
    assert plan['source_hashes']==manifest['source_hashes']
    for path,digest in manifest['source_hashes'].items():assert sha(ROOT/path)==digest,path
    tasks={t['id']:t for t in plan['tasks']}
    receipts=[]
    for path in sorted(RAW.glob('*.receipt.json')):
        r=read_json(path);contract=r['contract'];task=contract['task'];identifier=task['id']
        assert task==tasks[identifier],identifier
        assert contract['batch_manifest_sha256']==sha(RAW/'batch_manifest.json')
        for remote,digest in r['output_hashes'].items():
            assert sha(saved_path(remote))==digest,remote
        for dependency,digest in contract['dependency_hashes'].items():
            assert sha(RAW/(dependency+'.json'))==digest,dependency
        receipts.append(dict(task=identifier,kind=task['kind'],receipt_sha256=sha(path),
                             verified_output_count=len(r['output_hashes'])))
    events=[json.loads(line) for line in (RAW/'batch_events.jsonl').read_text().splitlines()]
    completed=[e['task'] for e in events if e['event']=='TASK_COMPLETE']
    assert sorted(completed)==sorted(r['task'] for r in receipts)
    criteria=read_json(ROOT/grid.CRITERIA)
    runs={};trajectory_rows=[]
    for path in sorted(RAW.glob('*.json')):
        record=read_json(path)
        if record.get('schema')!='pysnspd.stage2.coupled_run.v1':continue
        run=grid.read(path,ROOT)
        guarded.verify_contract({'record':record})
        runs[path.stem]=run
        state,cells,ne=blocks(run)
        p,n=state[:,:,1:1+ne],state[:,:,1+ne:]
        snaps=record['snapshots'];conserved=np.asarray([s['conserved'] for s in snaps])
        scale=max(1.,abs(snaps[0]['total']),float(np.sum(snaps[-1]['input'])))
        ledger=float(np.max(abs(conserved-conserved[0]))/scale)
        assert ledger==record['energy_ledger_scaled_max']
        assert float(p.min())==record['minimum_electron']
        assert float(p.max())==record['maximum_electron']
        assert float(n.min())==record['minimum_phonon']
        checks=temporal.invariant_checks(record,criteria)
        assert all(checks.values()),path.name
        for s in snaps:
            total=float(np.sum(s['electron_energy'])+np.sum(s['phonon_energy']))
            direct=total+float(np.sum(np.asarray(s['escape'])-np.asarray(s['input'])))
            assert total==s['total'] and abs(direct-s['conserved'])<1e-15
        stats=record['time_integration']['stats']
        params=record['parameters']
        trajectory_rows.append(dict(run=path.stem,cells=cells,scenario=params['scenario'],
            steps=params['steps'],phonon_nodes=params['phonon_nodes'],electron_nodes=ne,
            duration=params['duration'],runtime_seconds=record['runtime_seconds'],
            stored_checkpoints=len(snaps),ledger_scaled=ledger,
            instantaneous_residual_max_recorded=record['instantaneous_residual_max'],
            minimum_electron=float(p.min()),maximum_electron=float(p.max()),minimum_phonon=float(n.min()),
            invariants_pass=all(checks.values()),guard_extended_stages=stats['guarded_extended_steps'],
            guard_delegated_stages=stats['guarded_delegated_steps'],
            limited_event_evaluations=stats['limited_event_evaluations'],
            reaction_limited_evaluations=stats['reaction_limited_evaluations'],
            transport_limited_evaluations=stats.get('transport_limited_evaluations'),
            integrated_energy_weighted_flux_defect=stats['integrated_energy_weighted_flux_defect'],
            positive_values_rounded_to_zero=stats['guarded_positive_values_rounded_to_zero']))
    comparisons=[]
    for base in ('one','two'):
        for steps in (160,320,640):
            comparisons+=compare_runs(runs[f'{base}_guarded_{steps}'],runs[f'{base}_guarded_1280'],
                                      f'{base}_time_{steps}_vs_1280')
    for prefix in ('two_guarded','two_guarded_ph2049'):
        rows=compare_runs(runs[prefix+'_320'],runs[prefix+'_640'],prefix+'_pair_320_640')
        recorded=read_json(RAW/(prefix+'_pair_320_640.json'))
        assert recorded['reviewer_sha256']==sha(ROOT/'sandbox/stage2_cells/closure_prep_20260922/assess_guarded_pair.py')
        assert recorded['contract_sha256']==sha(guarded.__file__)
        for entry in recorded['runs']:assert sha(saved_path(entry['path']))==entry['sha256']
        formal=[r for r in rows if r['cell']=='all']
        assert max(r['original_norm_score'] for r in formal)==recorded['maximum_error']
        by_key={(r['time'],r['family']):r for r in formal}
        for point in recorded['checkpoints']:
            for family,error in point['errors'].items():
                r=by_key[point['time'],family]
                assert r['original_norm_score']==error['assessment_error']
                assert r['absolute']==error['absolute'] and r['reference_norm']==error['scale']
        comparisons+=rows
    comparisons+=compare_runs(runs['two_guarded_320'],runs['two_guarded_ph2049_320'],
                              'two_mesh_ph1025_vs_2049_steps320',mesh=True)
    summaries=reduced(comparisons)
    pair='two_guarded_ph2049_pair_320_640'
    worst=max((r for r in summaries if r['comparison']==pair and r['cell']=='all'),
              key=lambda r:r['original_norm_score'])
    bands,popdetail=phonon_difference(runs['two_guarded_ph2049_320'],
                                     runs['two_guarded_ph2049_640'],worst['time'])
    support={}
    for name in ('guarded_equilibrium_support','guarded_boundary_support',
                 'guarded_equilibrium_fields','guarded_boundary_fields','guarded_special_assessment'):
        r=read_json(RAW/(name+'.json'))
        support[name]={'status':r['status'],'sha256':sha(RAW/(name+'.json'))}
    assessment_statuses={}
    for name in ('one_guarded_time_assessment','two_guarded_time_assessment'):
        r=read_json(RAW/(name+'.json'))
        entries=[]
        for case in r['cases']:
            base='one' if name.startswith('one') else 'two'
            mine=[x for x in comparisons if x['comparison']==f'{base}_time_{case["steps"]}_vs_1280'
                  and x['cell']=='all' and x['time'] in r['shared_checkpoint_times']]
            assert max(x['original_norm_score'] for x in mine)==case['max_error']
            entries.append({'steps':case['steps'],'recorded_max_relative_percent':100*case['max_error']})
        assessment_statuses[name]={'status':r['status'],'cases':entries,
            'recorded_common_checkpoint_count':len(r['shared_checkpoint_times']),
            'note':'Reproduced at recorded common checkpoints; supplementary CSV also scans each coarse-run output.'}
    assert raw_hashes=={p.name:sha(p) for p in sorted(RAW.iterdir()) if p.is_file()}
    write_csv(OUT/'trajectory_audit.csv',trajectory_rows)
    write_csv(OUT/'comparison_metrics.csv',summaries)
    write_csv(OUT/'comparison_checkpoints.csv',comparisons)
    write_csv(OUT/'worst_phonon_error_bands.csv',bands)
    result=dict(schema='pysnspd.stage2.saved-results-audit.v1',
        status='AUDITED_SAVED_RESULTS_NOT_STAGE2_ADMISSION',
        scope='Read-only archive arithmetic; zero RHS evaluations, zero new trajectories, unchanged criteria and raw evidence.',
        script_sha256=sha(__file__),plan_sha256=sha(PLAN),raw_file_sha256=raw_hashes,
        manifest_sources_verified=len(manifest['source_hashes']),receipts=receipts,
        completed_task_count=len(completed),registered_task_count=len(tasks),
        stop_event=events[-1],trajectory_count=len(runs),
        instantaneous_balance_scope='Recorded initial/11-checkpoint RHS residuals inspected, not recomputed. Stored population support and snapshot/NPZ consistency checked at every saved time.',
        recorded_field_and_support_assessments=support,
        active_trajectory_field_and_upper_support_status='Not present in this stopped batch; population support is a separate checked property.',
        trajectory_checks=trajectory_rows,temporal_assessments=assessment_statuses,
        failed_pair_worst=worst,failed_pair_phonon_detail_at_worst_time=popdetail,
        formal_time_tolerance=criteria['coupled_trajectories']['finest_time_observable_relative_error_max'],
        supplemental_pair_budget=criteria['coupled_trajectories']['finest_time_observable_relative_error_max']/4,
        failed_pair_exceeds_supplemental_budget_factor=worst['original_norm_score']/(criteria['coupled_trajectories']['finest_time_observable_relative_error_max']/4),
        comparison_metrics=summaries,
        interpretation_limits=[
            'A two-resolution pair measures a difference, not an independently established error bound.',
            'Two phonon grids at320 steps measure combined mesh sensitivity and residual temporal error; no mesh convergence claim is made.',
            'Per-cell relative errors are descriptive; the registered test uses all-cell L1 population norms or all-cell maximum observable norms.',
            'The snapshot electron_energy includes a condensate/vacuum contribution; excitation_energy is reported separately.',
            'Percent values multiply dimensionless relative errors by100; ledger_scaled is its own frozen norm, not a signal error percentage.',
            'No measured device or circuit waveform appears in these synthetic one/two-cell trajectories.'
        ],runtime_seconds=time.perf_counter()-started)
    result['output_sha256']={name:sha(OUT/name) for name in outputs if name!='saved_results_audit.json'}
    (OUT/'saved_results_audit.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({k:result[k] for k in ('status','completed_task_count','trajectory_count','failed_pair_worst','runtime_seconds')},indent=2))


if __name__=='__main__':main()
