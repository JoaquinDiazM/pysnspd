"""Assess three true meshes at every common checkpoint, with file provenance.

Positive electron masses use 129 fixed count hats; phonons use 65 energy hats.
The first electron COUNT-coordinate moment is preserved. Physical electronic
energy is compared separately. No temporal interpolation or endpoint-only gate.
"""
from pathlib import Path
import argparse
import hashlib
import json
import tempfile
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
CRITERIA='docs/implementation/stage2/acceptance_criteria.json'
CATALOG='docs/implementation/stage1_r2/catalogs/occupation_catalog.npz'
OBSERVABLES=('amplitudes','electron_energy','excitation_energy','phonon_energy',
             'temperatures','forces','escape','input','electron_to_phonon',
             'condensate_heat','transport','total','conserved')
ROUND=100*np.finfo(float).eps


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def projected_mass(coordinates,populations,capacities,bands):
    x,p,c,b=(np.asarray(a,float) for a in (coordinates,populations,capacities,bands))
    if (x.ndim!=1 or c.shape!=x.shape or p.shape[-1:]!=x.shape or b.ndim!=1
            or np.any(~np.isfinite(x)) or np.any(~np.isfinite(p)) or np.any(~np.isfinite(c))
            or np.any(~np.isfinite(b)) or np.any(np.diff(x)<=0) or np.any(np.diff(b)<=0)
            or np.any(c<=0) or np.any(p<0) or x[0]<b[0] or x[-1]>b[-1]):
        raise ValueError('positive ordered finite grids and full projection support required')
    upper=np.searchsorted(b,x,side='left');lower=np.maximum(0,upper-1)
    widths=b[upper]-b[lower]
    fraction=np.divide(x-b[lower],widths,out=np.zeros_like(x),where=widths>0)
    rows=p.reshape(-1,len(x))
    projected=np.asarray([np.bincount(lower,weights=r*c*(1-fraction),minlength=len(b))
                          +np.bincount(upper,weights=r*c*fraction,minlength=len(b)) for r in rows])
    for expected,actual in ((np.sum(rows*c,axis=1),projected.sum(axis=1)),
                            (np.sum(rows*c*x,axis=1),projected@b)):
        if np.any(abs(actual-expected)>1e-12*np.maximum(1,abs(expected))):
            raise ValueError('population projection lost mass or first coordinate moment')
    return projected.reshape(*p.shape[:-1],len(b))


def source_file(root,name):
    path=(root/name).resolve()
    if not path.is_relative_to(root.resolve()):raise ValueError('source path escapes provenance root')
    return path


def read(path,source_root=ROOT):
    path=Path(path);record=json.loads(path.read_text(encoding='utf-8'))
    if record.get('schema')!='pysnspd.stage2.coupled_run.v1' or record.get('status')!='COMPLETED_NOT_YET_ADJUDICATED':
        raise ValueError('only completed coupled trajectories can be compared')
    sources=record.get('source_hashes')
    if not isinstance(sources,dict) or not sources:raise ValueError('source hashes are required')
    for name,expected in sources.items():
        if sha(source_file(source_root,name))!=expected:raise ValueError('source hash mismatch: '+name)
    for name,key in ((CRITERIA,'criteria_sha256'),(CATALOG,'catalog_sha256')):
        if sha(source_root/name)!=record.get(key):raise ValueError(key+' mismatch')
    npz_path=path.with_suffix('.npz');initial_path=path.with_name(path.stem+'_initial.npz')
    if sha(npz_path)!=record.get('trajectory_sha256'):raise ValueError('trajectory NPZ hash mismatch')
    if sha(initial_path)!=record.get('initial_conditions_sha256'):raise ValueError('initial NPZ hash mismatch')
    with np.load(npz_path,allow_pickle=False) as data:
        arrays={k:np.array(data[k],copy=True) for k in
                ('times','states','electron_count','electron_weights','phonon_energies','phonon_capacities')}
    if arrays['states'].ndim!=2 or arrays['states'].shape[0]<3:
        raise ValueError('trajectory requires initial, interior and final states')
    with np.load(initial_path,allow_pickle=False) as data:
        for key in ('electron_count','electron_weights','phonon_energies','phonon_capacities'):
            if not np.array_equal(arrays[key],data[key]):raise ValueError('initial/trajectory grid mismatch')
        if not np.array_equal(arrays['states'][0],data['state']):raise ValueError('initial population mismatch')
    t,y=arrays['times'],arrays['states'];params=record['parameters'];snapshots=record.get('snapshots',[])
    if (t.ndim!=1 or len(t)<3 or np.any(~np.isfinite(t)) or np.any(np.diff(t)<=0)
            or t[0]!=0 or t[-1]!=params['duration'] or len(t)!=params['steps']+1
            or not np.array_equal(t,np.asarray(record.get('times'))) or len(snapshots)!=len(t)
            or y.ndim!=2 or len(y)!=len(t) or np.any(~np.isfinite(y))):
        raise ValueError('times, steps, snapshots and trajectory disagree or lack interior checkpoints')
    if record['initial']!=snapshots[0] or record['final']!=snapshots[-1]:
        raise ValueError('initial/final summaries differ from checkpoint snapshots')
    x,w,om,cap=(arrays[k] for k in ('electron_count','electron_weights','phonon_energies','phonon_capacities'))
    for axis,weight in ((x,w),(om,cap)):
        if (axis.ndim!=1 or len(axis)<2 or weight.shape!=axis.shape or np.any(~np.isfinite(axis))
                or np.any(~np.isfinite(weight)) or np.any(axis<=0) or np.any(np.diff(axis)<=0) or np.any(weight<=0)):
            raise ValueError('invalid count or phonon quadrature')
    count=len(record['initial']['amplitudes']);block=1+len(x)+len(om)
    if count not in (1,2) or y.shape[1]!=count*block+4*count+1:
        raise ValueError('state layout does not match cell and grid sizes')
    if params['phonon_nodes']!=len(om) or params['case']!=('one' if count==1 else 'two'):
        raise ValueError('declared phonon/cell sizes differ from actual state')
    mesh=record.get('occupation_mesh')
    if mesh is not None and (mesh['electron_states']!=len(x) or mesh['refinement']!=params['electron_refinement']):
        raise ValueError('declared electronic mesh differs from actual state')
    blocks=y[:,:count*block].reshape(len(t),count,block)
    amplitudes=blocks[:,:,0];p=blocks[:,:,1:1+len(x)];n=blocks[:,:,1+len(x):]
    if np.any((p<0)|(p>1)) or np.any(n<0):raise ValueError('invalid stored Pauli/phonon populations')
    ledger=y[:,count*block:-1].reshape(len(t),count,4)
    direct={'amplitudes':amplitudes,'phonon_energy':np.sum(n*(cap*om),axis=-1),
            'escape':ledger[:,:,0],'input':ledger[:,:,1],'electron_to_phonon':ledger[:,:,2],
            'condensate_heat':ledger[:,:,3],'transport':y[:,-1]}
    for index,snapshot in enumerate(snapshots):
        for name in OBSERVABLES:
            if np.any(~np.isfinite(np.asarray(snapshot[name],float))):raise ValueError('nonfinite snapshot '+name)
        for name,actual in direct.items():
            expected=np.asarray(snapshot[name],float)
            if actual[index].shape!=expected.shape or np.any(abs(actual[index]-expected)>ROUND*np.maximum(1,abs(expected))):
                raise ValueError('snapshot/NPZ inconsistency: '+name)
    return {'record':record,'arrays':arrays,'path':path,'sha256':sha(path),
            'electron_projection':projected_mass(x,p,4*w,np.linspace(0,12,129)),
            'phonon_projection':projected_mass(om,n,cap,np.linspace(0,4,65))}


def error(value,expected,population=False):
    a,b=np.asarray(value,float),np.asarray(expected,float)
    if a.shape!=b.shape:raise ValueError('unmatched observable shapes')
    norm=(lambda x:float(np.sum(abs(x)))) if population else (lambda x:float(np.max(abs(x))))
    absolute,denominator=norm(a-b),norm(b);floor=ROUND*max(1.,denominator)
    if denominator<=floor:
        relative=None;gate=0. if absolute<=floor else None
        convention='degenerate reference: roundoff agreement only; no arbitrary unit denominator'
    else:
        relative=absolute/denominator;gate=0. if absolute<=floor else relative
        convention='projected mass L1/reference L1' if population else 'maximum difference/reference maximum magnitude'
    return {'absolute':absolute,'reference_norm':denominator,'relative':relative,
            'error_for_gate':gate,'roundoff_limited':absolute<=floor,'convention':convention}


def assess(paths,reference,parameter,source_root=ROOT):
    if parameter not in ('electron_refinement','phonon_nodes'):
        raise ValueError('only electron_refinement or phonon_nodes comparisons are supported')
    all_paths=[*map(Path,paths),Path(reference)]
    if len(paths)<2 or len({p.resolve() for p in all_paths})!=len(all_paths):
        raise ValueError('at least three distinct true mesh trajectories are required')
    runs=[read(p,source_root) for p in all_paths];ref=runs[-1];record=ref['record']
    baseline={k:v for k,v in record['parameters'].items() if k not in (parameter,'output')}
    values=[];sizes=[]
    target='electron_count' if parameter=='electron_refinement' else 'phonon_energies'
    held=('phonon_energies','phonon_capacities') if parameter=='electron_refinement' else ('electron_count','electron_weights')
    for run in runs:
        r=run['record'];params={k:v for k,v in r['parameters'].items() if k not in (parameter,'output')}
        if params!=baseline:raise ValueError('physical/time/method/tolerance parameters differ beyond '+parameter)
        if r['source_hashes']!=record['source_hashes']:raise ValueError('source sets/hashes differ')
        for key in ('criteria_sha256','catalog_sha256','material','python','numpy','scipy'):
            if r.get(key)!=record.get(key):raise ValueError('physical/provenance/dependency contract differs: '+key)
        if any(not np.array_equal(run['arrays'][key],ref['arrays'][key]) for key in held):
            raise ValueError('the unrefined population grid changed')
        values.append(r['parameters'][parameter]);sizes.append(len(run['arrays'][target]))
    if (any(not isinstance(v,(int,float)) or not np.isfinite(v) for v in values)
            or any(b<=a for a,b in zip(values[:-1],values[1:]))
            or any(b<=a for a,b in zip(sizes[:-1],sizes[1:]))):
        raise ValueError('refinement values AND actual grid sizes must increase coarse to reference')
    times=runs[0]['arrays']['times'];indices=[[] for _ in runs];common=[]
    tolerance=32*np.finfo(float).eps*max(1.,abs(times[-1]))
    for t in times:
        match=[]
        for run in runs:
            axis=run['arrays']['times'];i=int(np.argmin(abs(axis-t)))
            if abs(axis[i]-t)>tolerance:break
            match.append(i)
        if len(match)==len(runs):
            common.append(float(t))
            for dest,i in zip(indices,match):dest.append(i)
    if len(common)<3 or common[0]!=times[0] or common[-1]!=times[-1]:
        raise ValueError('common initial, final and interior checkpoints are required')
    criteria=json.loads((source_root/CRITERIA).read_text(encoding='utf-8'))
    limit=criteria['continuous_consistency']['relative_error_max']
    names=[*OBSERVABLES,'electron_129_bands','phonon_65_bands'];rows=[]
    for k,run in enumerate(runs[:-1]):
        checkpoints=[]
        for t,i,j in zip(common,indices[k],indices[-1]):
            errors={name:error(run['record']['snapshots'][i][name],record['snapshots'][j][name]) for name in OBSERVABLES}
            for name,key in (('electron_129_bands','electron_projection'),('phonon_65_bands','phonon_projection')):
                errors[name]=error(run[key][i],ref[key][j],True)
            checkpoints.append({'time':t,'errors':errors})
        metrics={}
        for name in names:
            entries=[r['errors'][name] for r in checkpoints];admissible=all(e['error_for_gate'] is not None for e in entries)
            maximum=max(e['error_for_gate'] for e in entries) if admissible else None
            metrics[name]={'initial':entries[0],'final':entries[-1],
                'maximum_error_for_gate':maximum,'maximum_absolute_error':max(e['absolute'] for e in entries),
                'all_roundoff_limited':all(e['roundoff_limited'] for e in entries),'passes':admissible and maximum<=limit}
        finite=all(m['maximum_error_for_gate'] is not None for m in metrics.values())
        rows.append({'path':run['path'].as_posix(),'sha256':run['sha256'],
            'trajectory_sha256':run['record']['trajectory_sha256'],'parameter_value':values[k],
            'actual_refined_grid_size':sizes[k],'checkpoints':checkpoints,'metrics':metrics,
            'errors':{name:m['maximum_error_for_gate'] for name,m in metrics.items()},
            'maximum_error':max(m['maximum_error_for_gate'] for m in metrics.values()) if finite else None,
            'precision_pass':all(m['passes'] for m in metrics.values())})
    convergence={}
    for name in names:
        series=[r['metrics'][name] for r in rows]
        convergence[name]=all(b['all_roundoff_limited'] or
            (a['maximum_error_for_gate'] is not None and b['maximum_error_for_gate'] is not None
             and b['maximum_error_for_gate']<a['maximum_error_for_gate']) for a,b in zip(series[:-1],series[1:]))
    convergent=all(convergence.values());passed=rows[-1]['precision_pass'] and convergent
    return {'schema':'pysnspd.stage2.grid-refinement-assessment.v2','status':'PASS' if passed else 'FAIL',
        'parameter':parameter,'rows':rows,'reference_path':ref['path'].as_posix(),
        'reference_sha256':ref['sha256'],'reference_trajectory_sha256':record['trajectory_sha256'],
        'reference_parameter_value':values[-1],'reference_actual_grid_size':sizes[-1],
        'common_checkpoint_times':common,'common_checkpoint_count':len(common),
        'population_basis':{'electrons':'129 fixed linear count hats on[0,12]; count and first COUNT-coordinate moment conserved, physical energy compared separately',
                            'phonons':'65 fixed linear energy hats on[0,4]; phonon number and energy conserved'},
        'relative_tolerance':limit,'criteria_sha256':record['criteria_sha256'],
        'maximum_finest_error':rows[-1]['maximum_error'],'convergent':convergent,
        'convergence_by_observable':convergence,'reviewer_sha256':sha(__file__),
        'provenance_verified':True,'same_physical_and_time_contract':True,
        'scope':'Numerical mesh convergence for the registered trajectory only; no long-transient or material admission.'}


def self_test():
    """Nonphysical fake files exercise validation, not an accepted trajectory."""
    started=time.perf_counter();checks={}
    with tempfile.TemporaryDirectory(prefix='grid_assessment_') as temp:
        root=Path(temp)
        for name,data in ((CRITERIA,(ROOT/CRITERIA).read_bytes()),(CATALOG,b'synthetic hash fixture'),('fixture.py',b'# synthetic source\n')):
            p=root/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(data)
        def fixture(level,epsilon,mode='electron'):
            count=4*level+1 if mode=='electron' else 5
            phonon_count=5 if mode=='electron' else 4*level+1
            x=np.linspace(1,1.001,count);w=np.full(count,.001/(count-1));w[[0,-1]]/=2
            om=np.linspace(1,1.001,phonon_count);cap=np.full(phonon_count,.001/(phonon_count-1));cap[[0,-1]]/=2
            times=np.linspace(0,1,5);block=1+count+phonon_count;states=np.zeros((5,block+5));snapshots=[]
            for i,t in enumerate(times):
                perturbation=epsilon*np.sin(np.pi*t)
                p=np.full(count,.2+.01*t+(perturbation if mode=='electron' else 0))
                n=np.full(phonon_count,.3+.02*t+(perturbation if mode=='phonon' else 0))
                states[i,:block]=np.r_[.7,p,n];ep=float(np.dot(4*w,p));ph=float(np.dot(cap*om,n))
                snapshots.append(dict(amplitudes=[.7],electron_energy=[ep],excitation_energy=[ep],phonon_energy=[ph],
                    temperatures=[.2],forces=[.1],escape=[0.],input=[0.],electron_to_phonon=[0.],condensate_heat=[0.],
                    transport=0.,total=ep+ph,conserved=ep+ph))
            path=root/f'{mode}_level_{level}.json';npz=path.with_suffix('.npz');initial=path.with_name(path.stem+'_initial.npz')
            grids=dict(electron_count=x,electron_weights=w,phonon_energies=om,phonon_capacities=cap)
            np.savez(npz,times=times,states=states,**grids);np.savez(initial,state=states[0],**grids)
            record=dict(schema='pysnspd.stage2.coupled_run.v1',status='COMPLETED_NOT_YET_ADJUDICATED',
                parameters=dict(case='one',steps=4,duration=1.,method='rk4',rtol=1e-9,
                    electron_refinement=level if mode=='electron' else 1,phonon_nodes=phonon_count,output=str(path)),
                source_hashes={'fixture.py':sha(root/'fixture.py')},criteria_sha256=sha(root/CRITERIA),catalog_sha256=sha(root/CATALOG),
                initial_conditions_sha256=sha(initial),trajectory_sha256=sha(npz),python='synthetic',numpy='synthetic',scipy='synthetic',
                material={'status':'NONPHYSICAL_ASSESSOR_TEST_FIXTURE'},times=times.tolist(),snapshots=snapshots,initial=snapshots[0],final=snapshots[-1])
            write(path,record);return path
        paths=[fixture(1,2e-5),fixture(2,5e-6),fixture(4,0.)]
        valid=assess(paths[:2],paths[-1],'electron_refinement',root)
        checks['three_true_meshes_and_all_checkpoints_pass']=valid['status']=='PASS' and valid['common_checkpoint_count']==5
        phonon_paths=[fixture(1,2e-5,'phonon'),fixture(2,5e-6,'phonon'),fixture(4,0.,'phonon')]
        checks['three_true_phonon_meshes_pass']=assess(phonon_paths[:2],phonon_paths[-1],'phonon_nodes',root)['status']=='PASS'
        checks['nonzero_error_with_degenerate_reference_is_not_false_PASS']=error([1e-6],[0.])['error_for_gate'] is None
        saved=paths[1].read_bytes();npz=paths[1].with_suffix('.npz');saved_npz=npz.read_bytes()
        initial=paths[1].with_name(paths[1].stem+'_initial.npz');saved_initial=initial.read_bytes()
        def rejects(change):
            change()
            try:assess(paths[:2],paths[-1],'electron_refinement',root)
            except (ValueError,FileNotFoundError):return True
            finally:paths[1].write_bytes(saved);npz.write_bytes(saved_npz);initial.write_bytes(saved_initial)
            return False
        def json_change(fn):
            r=json.loads(saved);fn(r);write(paths[1],r)
        checks['changed_source_hash_rejected']=rejects(lambda:json_change(lambda r:r['source_hashes'].update({'fixture.py':'0'*64})))
        checks['changed_NPZ_hash_rejected']=rejects(lambda:npz.write_bytes(saved_npz+b'tamper'))
        checks['different_time_method_rejected']=rejects(lambda:json_change(lambda r:r['parameters'].update(method='dop853')))
        checks['different_tolerance_rejected']=rejects(lambda:json_change(lambda r:r['parameters'].update(rtol=1e-4)))
        checks['different_physical_problem_rejected']=rejects(lambda:json_change(lambda r:r['material'].update(status='other input')))
        checks['nonincreasing_parameter_rejected']=rejects(lambda:json_change(lambda r:r['parameters'].update(electron_refinement=4)))
        def duplicate_actual_grid():
            coarse=json.loads(paths[0].read_text())
            npz.write_bytes(paths[0].with_suffix('.npz').read_bytes())
            initial.write_bytes(paths[0].with_name(paths[0].stem+'_initial.npz').read_bytes())
            r=json.loads(saved)
            for key in ('initial','final','snapshots'):r[key]=coarse[key]
            r['trajectory_sha256']=sha(npz);r['initial_conditions_sha256']=sha(initial);write(paths[1],r)
        checks['declared_refinement_without_real_grid_change_rejected']=rejects(duplicate_actual_grid)
        # Only an interior electron population changes; endpoints remain equal.
        # Refresh its hash to ensure the population comparison itself detects it.
        with np.load(npz) as data:arrays={k:data[k].copy() for k in data.files}
        arrays['states'][2,1:1+len(arrays['electron_count'])]+=.02
        np.savez(npz,**arrays);r=json.loads(saved);r['trajectory_sha256']=sha(npz);write(paths[1],r)
        bad=assess(paths[:2],paths[-1],'electron_refinement',root)
        pop=bad['rows'][1]['metrics']['electron_129_bands']
        checks['interior_population_error_detected_despite_equal_endpoints']=(bad['status']=='FAIL' and pop['initial']['error_for_gate']<1e-12
            and pop['final']['error_for_gate']<1e-12 and pop['maximum_error_for_gate']>1e-3)
    return {'schema':'pysnspd.stage2.grid-assessor-self-test.v1','status':'PASS' if all(checks.values()) else 'FAIL',
        'scope':'Nonphysical synthetic files test acceptance guards, not device physics.',
        'checks':checks,'runner_sha256':sha(__file__),'runtime_seconds':time.perf_counter()-started}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('runs',nargs='*',type=Path);parser.add_argument('--reference',type=Path)
    parser.add_argument('--parameter',choices=('electron_refinement','phonon_nodes'))
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--source-root',type=Path,default=ROOT)
    parser.add_argument('--self-test',action='store_true');args=parser.parse_args()
    try:
        if args.self_test:result=self_test()
        else:
            if args.reference is None or args.parameter is None:raise ValueError('reference and parameter are required')
            result=assess(args.runs,args.reference,args.parameter,args.source_root)
    except (ValueError,KeyError,OSError) as exc:
        result={'status':'INVALID_INPUT','reason':str(exc),'reviewer_sha256':sha(__file__)}
    write(args.output,result)
    print(json.dumps({k:result[k] for k in ('status','maximum_finest_error','convergent','reason','checks','runtime_seconds') if k in result},indent=2))
    if result['status']!='PASS':raise SystemExit(1)


if __name__=='__main__':main()
