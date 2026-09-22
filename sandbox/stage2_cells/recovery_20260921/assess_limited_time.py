"""Explicitly compare limited SSPRK3 with frozen, independently recorded runs.

Physical kernels, grids and initial conditions must be identical. Only the two
declared integration-source files may be additional to the frozen source table.
The unchanged stage-2 time criteria and all invariant guards are still applied.
"""
from pathlib import Path
import copy,hashlib,json,sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'sandbox/stage2_cells'))
import assess_time_refinement as original

EXTRA={'sandbox/stage2_cells/recovery_20260921/run_limited_coupled.py',
       'sandbox/stage2_cells/recovery_20260921/limited_ssp.py'}
CATALOG='docs/implementation/stage1_r2/catalogs/occupation_catalog.npz'
CRITERIA='docs/implementation/stage2/acceptance_criteria.json'
LIMITED_METHOD='ssprk3_common_flux_limited'
sha=lambda path:hashlib.sha256(Path(path).read_bytes()).hexdigest()
old_check=original.require_same_problem

def verify_contract(run):
    record=run['record'];sources=record['source_hashes']
    base={'sandbox/stage2_cells/run_coupled.py',
          *[p.relative_to(ROOT).as_posix() for p in (ROOT/'pysnspd/experimental').glob('*.py')]}
    method=record['parameters']['method']
    if method not in ('rk4','dop853',LIMITED_METHOD):
        raise ValueError('Unrecognized actual integration method')
    expected=base|EXTRA if method==LIMITED_METHOD else base
    if not isinstance(sources,dict) or set(sources)!=expected:
        raise ValueError('Only the complete frozen source table and two declared limited-integration sources are allowed')
    for p,value in sources.items():
        if sha(ROOT/p)!=value:raise ValueError('Recorded source differs from current source: '+p)
    for key,path in (('catalog_sha256',CATALOG),('criteria_sha256',CRITERIA)):
        if record.get(key)!=sha(ROOT/path):
            raise ValueError('Recorded '+key+' differs from the current frozen input')
    if method==LIMITED_METHOD:
        adapter=record.get('time_integration',{})
        stats=adapter.get('stats',{})
        if (adapter.get('population_clipping') is not False
                or adapter.get('posthoc_energy_repair') is not False
                or adapter.get('method')!='SSPRK3 with common donor-availability event-flux factors'
                or stats.get('method')!='experimental conservative-event-limited SSPRK3'
                or stats.get('status')!='COMPLETED_REQUIRES_CONVERGENCE_ASSESSMENT'
                or stats.get('completed_steps')!=record['parameters']['steps']):
            raise ValueError('Missing explicit limited-integration source contract')
    elif 'time_integration' in record:
        raise ValueError('Unrecognized integration adapter')
    result=copy.copy(run);result['record']=copy.deepcopy(record)
    result['record']['source_hashes']={p:v for p,v in sources.items() if p not in EXTRA}
    return result

def compare(actual,reference):
    old_check(verify_contract(actual),verify_contract(reference))
    for key in ('material','python','numpy','scipy'):
        if actual['record'].get(key)!=reference['record'].get(key):
            raise ValueError('Physical input/dependency contract differs: '+key)

def main():
    output=Path(sys.argv[sys.argv.index('--output')+1])
    if output.exists():raise SystemExit('Do not overwrite prior assessment.')
    original.require_same_problem=compare
    try:original.main()
    finally:
        original.require_same_problem=old_check
        if output.exists():
            record=json.loads(output.read_text())
            record['integration_contract_comparison']=dict(
                policy='Identical frozen physical kernels; explicit integration adapters verified by current source hashes.',
                allowed_extra_sources=sorted(EXTRA),assessor_sha256=sha(__file__),
                final_stage2_admission=False)
            output.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')

if __name__=='__main__':main()
