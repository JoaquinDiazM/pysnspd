"""Supplemental reverse-channel coverage of the frozen guarded prototype."""
from pathlib import Path
import argparse
import json
import time
import numpy as np
import isolated_checks as base
import limited_ssp_guarded as guarded

ROOT=base.ROOT
OUT=ROOT/'docs/implementation/stage2/closure_prep_20260922/isolated_guarded_negative'
PLAN=OUT/'isolated_registration.json'


def scalar(value):
    if isinstance(value,np.generic):return value.item()
    raise TypeError(type(value).__name__)


def write_new(path,value):
    encoded=json.dumps(value,indent=2,allow_nan=False,default=scalar)+'\n'
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x',encoding='utf-8') as f:f.write(encoded)


def sources():
    result=base.sources()
    for path in (Path(__file__),Path(guarded.__file__),ROOT/'sandbox/stage2_cells/recovery_20260921/test_limited_ssp.py'):
        result[path.relative_to(ROOT).as_posix()]=base.sha(path)
    return result


def registration():
    return dict(schema='pysnspd.stage2.guarded-negative-channel-registration.v1',
        status='PREREGISTERED_NOT_EXECUTED',source_hashes=sources(),
        reason='Independent review found that the prior n=0 toy did not exercise negative reaction rates or phonon consumption under the fallback.',
        scope='Four-state normal spectral oracle and six phonon nodes; same projected event law and guarded SSP map. Not selected-mesh or global time admission.',
        cases=dict(absorption=dict(p=[0.,.2,0.,0.],family='scattering',direction='backward'),
                   creation=dict(p=[0.,0.,0.,0.],family='recombination',direction='backward')),
        common=dict(amplitude=0.,gamma=0.,n=[1e-320,.05,.05,.05,.05,.05],
                    duration=.02,steps=20,rate_prefactor=50.,face=0.,
                    BGK=False,heating=False,escape=False,condensate_motion=False),
        gates=dict(first_step_guard_required=True,negative_rates_required=True,
                   strict_populations=True,energy_scaled_max=1e-7,count_scaled_max=1e-12,
                   minimum_saved_times=10,net_phonon_consumption=True,
                   clipping=False,energy_repair=False),
        execution=dict(wall_timeout_seconds=240,expected_seconds='under1',threads=1,max_attempts=1))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--register',action='store_true')
    args=parser.parse_args()
    if args.register:
        write_new(PLAN,registration());print('PREREGISTERED',PLAN);return
    if json.loads(PLAN.read_text(encoding='utf-8'))!=registration():raise ValueError('preregistration mismatch')
    if (OUT/'isolated_results.json').exists() or list(OUT.glob('*.npz')):raise FileExistsError('already attempted')
    guarded.require_extended()
    from test_limited_ssp import make_system
    started=time.perf_counter()
    output=dict(schema='pysnspd.stage2.guarded-negative-channel-results.v1',status='INCOMPLETE',
                registration_sha256=base.sha(PLAN),source_hashes=sources(),cases=[],
                scope=registration()['scope'])
    old=(base.OUT,base.PLAN,base.integrate,base.write_new)
    base.OUT,base.PLAN,base.integrate,base.write_new=OUT,PLAN,guarded.integrate,write_new
    try:
        template=make_system(cells=1,reaction_rate=50.,face=0.)
        cell=base.ElectronicCell(template.catalog,0.,0.)
        parent=template.reaction_events(cell)
        for name,case in registration()['cases'].items():
            mask=parent.recombination if name=='creation' else ~parent.recombination
            event=base.DirectionalEvents(parent,mask,'backward')
            system=base.fixed_system(template,[0.],[event])
            p=np.array(case['p']);n=np.array([1e-320,.05,.05,.05,.05,.05])
            initial=system.pack([0.],[p],[n])
            rates=event.rates(p,n)
            dp,dn=event.rhs_from_rates(rates)
            prechecks=dict(initial_guard_active=guarded.needs_extended(system,initial),
                           negative_rates=bool(np.any(rates<0)),no_positive_rates=bool(np.all(rates<=0)),
                           raw_phonon_consumption=bool(np.dot(system.phonons.capacities,dn)<0))
            if not all(prechecks.values()):raise ValueError('negative channel not exercised: '+name)
            record,states=base.save_case(name,system,initial,.02,20)
            final_p,final_n=system.unpack(states[-1])[1:]
            checks=dict(prechecks,
                fallback_executed=record['integrator_stats']['guarded_extended_steps']>0,
                net_phonon_consumption=record['final']['phonon_count']<record['initial']['phonon_count'],
                invariant_checks=all(record['gates'].values()),
                stage_accounting=(record['integrator_stats']['guarded_extended_steps']+
                                  record['integrator_stats']['guarded_delegated_steps']==60))
            output['cases'].append(dict(case=name,checks=checks,energy_error=record['energy_error'],
                count_error=record['count_error'],phonon_count_change=record['final']['phonon_count']-record['initial']['phonon_count'],
                electron_count_change=record['final']['electron_count']-record['initial']['electron_count'],
                guarded_steps=record['integrator_stats']['guarded_extended_steps'],
                first_phonon_initial=float(n[0]),first_phonon_final=float(final_n[0,0])))
            if not all(checks.values()):raise ValueError('reverse-channel check failed: '+name)
        output['status']='PASS_GUARDED_REVERSE_CHANNEL_MICROTESTS_ONLY'
    except Exception as exc:
        output['status']='FAIL_OR_INCOMPLETE';output['error']=repr(exc);raise
    finally:
        base.OUT,base.PLAN,base.integrate,base.write_new=old
        output['runtime_seconds']=time.perf_counter()-started
        write_new(OUT/'isolated_results.json',output)
        print(json.dumps(output,default=scalar),flush=True)


if __name__=='__main__':main()
