"""New guarded SSP temporal assessment with strict independent source checks.

Only the newly registered guarded integrator is admitted to this comparison.
Historical SSP/RK results are not substituted for a new guarded trajectory.
The numerical criteria and norm implementation remain those frozen in stage2.
"""
from pathlib import Path
from decimal import Decimal,InvalidOperation
import copy,hashlib,json,sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'sandbox/stage2_cells'))
import assess_time_refinement as original

EXTRA={'sandbox/stage2_cells/closure_prep_20260922/run_guarded_coupled.py',
       'sandbox/stage2_cells/closure_prep_20260922/limited_ssp_guarded.py',
       'sandbox/stage2_cells/recovery_20260921/limited_ssp.py'}
CATALOG='docs/implementation/stage1_r2/catalogs/occupation_catalog.npz'
CRITERIA='docs/implementation/stage2/acceptance_criteria.json'
GUARDED_METHOD='ssprk3_common_flux_guarded'
ADAPTER_DESCRIPTION='SSPRK3 with common conservative event-flux factors and subnormal arithmetic guard'
STATS_DESCRIPTION='experimental conservative-event-limited SSPRK3 with subnormal arithmetic guard'
sha=lambda path:hashlib.sha256(Path(path).read_bytes()).hexdigest()
old_check=original.require_same_problem


def verify_contract(run):
    record=run['record'];sources=record['source_hashes']
    base={'sandbox/stage2_cells/run_coupled.py',
          *[p.relative_to(ROOT).as_posix() for p in (ROOT/'pysnspd/experimental').glob('*.py')]}
    if record['parameters']['method']!=GUARDED_METHOD:
        raise ValueError('Only the newly registered guarded SSP integrator is accepted')
    if not isinstance(sources,dict) or set(sources)!=base|EXTRA:
        raise ValueError('Complete physical source table and exactly three guarded integration sources required')
    for path,digest in sources.items():
        if sha(ROOT/path)!=digest:raise ValueError('Recorded source changed: '+path)
    for key,path in (('catalog_sha256',CATALOG),('criteria_sha256',CRITERIA)):
        if record.get(key)!=sha(ROOT/path):raise ValueError('Recorded frozen input changed: '+key)
    adapter=record.get('time_integration',{});stats=adapter.get('stats',{})
    if (adapter.get('population_clipping') is not False
            or adapter.get('posthoc_energy_repair') is not False
            or adapter.get('method')!=ADAPTER_DESCRIPTION
            or stats.get('method')!=STATS_DESCRIPTION
            or stats.get('status')!='COMPLETED_REQUIRES_CONVERGENCE_ASSESSMENT'
            or stats.get('completed_steps')!=record['parameters']['steps']):
        raise ValueError('Missing or incomplete guarded integration metadata')
    counters=('guarded_delegated_steps','guarded_extended_steps','guarded_positive_values_rounded_to_zero')
    if any(type(stats.get(key)) is not int or stats[key]<0 for key in counters):
        raise ValueError('Explicit nonnegative guarded arithmetic counters required')
    if (stats[counters[0]]+stats[counters[1]]!=3*record['parameters']['steps']
            or stats.get('forward_euler_stages')!=3*record['parameters']['steps']):
        raise ValueError('Every SSP forward-Euler stage must have one recorded arithmetic branch')
    capability=stats.get('guarded_capability',{})
    if (stats.get('guarded_prototype') is not True or capability.get('extended') is not True
            or type(capability.get('longdouble_bits')) is not int
            or capability['longdouble_bits']<=64
            or capability.get('longdouble_mantissa_bits',0)<=52
            or capability.get('longdouble_min_exponent',0)>=-1022):
        raise ValueError('Recorded platform lacks the required extended arithmetic capability')
    for key in ('guarded_extent_defect_longdouble','guarded_energy_flux_defect_longdouble',
                'guarded_minimum_event_factor_longdouble'):
        try:
            value=Decimal(stats[key])
            valid=isinstance(stats[key],str) and value.is_finite() and value>=0
        except (InvalidOperation,KeyError,TypeError):valid=False
        if not valid or (key.endswith('factor_longdouble') and value>1):
            raise ValueError('Invalid extended arithmetic diagnostic: '+key)
    result=copy.copy(run);result['record']=copy.deepcopy(record)
    result['record']['source_hashes']={path:digest for path,digest in sources.items() if path not in EXTRA}
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
                policy='New guarded trajectories only; identical frozen physical kernels and exact declared integration sources verified by current hashes.',
                allowed_extra_sources=sorted(EXTRA),assessor_sha256=sha(__file__),
                final_stage2_admission=False)
            output.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':main()
