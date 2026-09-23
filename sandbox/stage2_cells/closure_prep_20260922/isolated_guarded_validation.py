"""One bounded, preregistered validation of the guarded SSP prototype."""
from pathlib import Path
import argparse
import json
import time
import numpy as np
import isolated_checks as base
import limited_ssp_guarded as guarded
import limited_ssp as original

ROOT=base.ROOT
OUT=ROOT/'docs/implementation/stage2/closure_prep_20260922/isolated_guarded'
PLAN=OUT/'isolated_registration.json'
PREVIOUS=ROOT/'docs/implementation/stage2/closure_prep_20260922/isolated_transport_diagnostic'


def scalar(value):
    if isinstance(value,np.generic):return value.item()
    raise TypeError(type(value).__name__)


def write_new(path,value):
    encoded=json.dumps(value,indent=2,allow_nan=False,default=scalar)+'\n'
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x',encoding='utf-8') as f:f.write(encoded)


def sources():
    result=base.sources()
    for path in [Path(__file__),Path(guarded.__file__),
                 ROOT/'sandbox/stage2_cells/recovery_20260921/test_limited_ssp.py']:
        result[path.relative_to(ROOT).as_posix()]=base.sha(path)
    return result


def registration():
    return dict(schema='pysnspd.stage2.guarded-prototype-registration.v1',status='PREREGISTERED_NOT_EXECUTED',
        source_hashes=sources(),capability_requirement='longdouble wider exponent AND precision than float64',
        guard_trigger='positive electron, hole or phonon capacity-weighted inventory <1024*float64.tiny',
        tasks=[dict(id='normal_bytewise',case='one selected630/1025 moving condensate',duration=.01,steps=5,
                    reference='frozen original SSP, same initial state; compare every macro-output bitwise'),
               dict(id='subnormal_budget_oracles',scope='longdouble scalar budgets, zero faces, and detection of phonon inventory'),
               dict(id='competing_channels_toy',scope='same physical event law, 4-state normal oracle with reaction+face; tiny populated tail'),
               dict(id='previous_rejected_stage',scope='guarded map only, identical saved FE input and h=.001'),
               dict(id='transport_fixed_fields',duration=.02,steps=20,amplitudes=[.5,1.],
                    initial='same registered exact p=0/1 band and phonon vacuum as failed original case'),
               dict(id='escape',duration=.2,tau=.7,steps=[10,20,40],
                    reference='exact exponential and bytewise comparison against every saved original escape state')],
        reference_hashes={p.relative_to(ROOT).as_posix():base.sha(p) for p in
            [PREVIOUS/'isolated_transport_rejected_state.npz',PREVIOUS/'isolated_transport_diagnosis.json',
             *[PREVIOUS/f'isolated_escape_{n}.npz' for n in (10,20,40)]]},
        gates=dict(strict_populations=True,energy_scaled_max=1e-7,count_scaled_max=1e-12,
                   smooth_error_reduction_min=6.,normal_branch_bytewise=True,clipping=False,energy_repair=False),
        execution=dict(threads=1,wall_timeout_seconds=240,expected_seconds='under60',max_attempts=1),
        scope='Prototype only. No global stage2 or production admission; new temporal and mesh validation remain required unless full-trajectory equivalence is measured.')


def run():
    if json.loads(PLAN.read_text(encoding='utf-8'))!=registration():raise ValueError('preregistration/source mismatch')
    if (OUT/'isolated_results.json').exists() or list(OUT.glob('*.npz')):raise FileExistsError('already attempted')
    guarded.require_extended()
    started=time.perf_counter()
    result=dict(schema='pysnspd.stage2.guarded-prototype-results.v1',status='INCOMPLETE',
                source_hashes=sources(),registration_sha256=base.sha(PLAN),capability=guarded.capability(),checks={},evidence={})
    checks=result['checks']
    old_out,old_plan,old_integrate,old_writer=base.OUT,base.PLAN,base.integrate,base.write_new
    base.OUT,base.PLAN,base.integrate,base.write_new=OUT,PLAN,guarded.integrate,write_new
    try:
        # Test both the actual moving-cell map and its integration history.
        old_system,initial,_=base.setup('one',heating=0.,escape=np.inf)
        new_system,new_initial,_=base.setup('one',heating=0.,escape=np.inf)
        if not np.array_equal(initial,new_initial):raise ValueError('normal oracle initial mismatch')
        t0=time.perf_counter()
        a,b,s0=original.integrate(old_system,initial,.01,5)
        c,d,s1=guarded.integrate(new_system,new_initial,.01,5)
        checks['normal_full_trajectory_bytewise']=(a.dtype==c.dtype and b.dtype==d.dtype and
                                                  a.shape==c.shape and b.shape==d.shape and
                                                  a.tobytes()==c.tobytes() and b.tobytes()==d.tobytes())
        checks['normal_branch_delegates_all_stages']=s1['guarded_delegated_steps']==15 and s1['guarded_extended_steps']==0
        result['evidence']['normal']=dict(runtime_seconds=time.perf_counter()-t0,stats=s1,
                                        max_state_difference=float(np.max(abs(b-d))))
        np.savez_compressed(OUT/'isolated_normal_equivalence.npz',times=a,original_states=b,guarded_states=d)

        # Underflow witnesses: computations here stay in longdouble until cast.
        LD=np.longdouble
        q=np.nextafter(0.,1.)
        available=LD(.04)*LD(70*q)
        demand=LD('6.178685682638e-311');h=LD(.001)
        theta=guarded._ratio(np.array([available],dtype=LD),np.array([demand],dtype=LD),h)[0]
        updated=LD(70*q)-h*theta*demand/LD(.04)
        checks['subnormal_inventory_not_rounded_before_budget']=available>0 and updated>=0 and updated<=LD(70*q)
        checks['zero_inventory_blocks_outflow']=guarded._ratio(np.array([0],dtype=LD),np.array([demand],dtype=LD),h)[0]==0
        result['evidence']['subnormal_scalar']=dict(initial=str(LD(70*q)),exact_inventory=str(available),
                                                  common_factor=str(theta),updated=str(updated),cast=float(updated))
        phonon_probe=new_initial.copy()
        phonon_probe[1+new_system.electron_size]=q
        checks['subnormal_phonon_inventory_selects_fallback']=guarded.needs_extended(new_system,phonon_probe)

        from test_limited_ssp import make_system,energy
        toy=make_system(reaction_rate=200.,face=200.)
        toy_state=toy.pack([0.,0.],[[.1,.4,.2,1e-320],[0.,0.,0.,0.]],np.zeros((2,6)))
        toy_stats={};obtained=guarded.step(toy,0.,toy_state,1.,toy_stats)
        _,tp,tn=toy.unpack(obtained)
        checks['competing_reaction_face_support']=np.all((tp>=0)&(tp<=1)) and np.all(tn>=0)
        checks['competing_reaction_face_energy']=abs(energy(toy,obtained)-energy(toy,toy_state))<=1e-12
        checks['competing_reaction_face_common_limiter_active']=toy_stats.get('guarded_extended_steps')==1 and toy_stats.get('reaction_limited_evaluations',0)>0 and toy_stats.get('transport_limited_evaluations',0)>0
        result['evidence']['competing_channels']=dict(stats=toy_stats,energy_difference=float(energy(toy,obtained)-energy(toy,toy_state)))

        template,_,_=base.setup('one',heating=0.,escape=np.inf)
        parent=template.reaction_events(base.ElectronicCell(template.catalog,.5,0.))
        empty=base.DirectionalEvents(parent,np.zeros(len(parent.omega),bool),'disabled')
        del parent
        fixed=base.fixed_system(template,[.5,1.],[empty,empty],transport=True)
        with np.load(PREVIOUS/'isolated_transport_rejected_state.npz',allow_pickle=False) as saved:
            failed_input=np.array(saved['stage_input'],copy=True)
            transport_initial=np.array(saved['initial_state'],copy=True)
            old_rejected=np.array(saved['rejected_state'],copy=True)
        probe_stats={};probe=guarded.step(fixed,.017,failed_input,.001,probe_stats)
        _,p,n=fixed.unpack(probe)
        checks['previous_rejected_stage_physical']=np.all((p>=0)&(p<=1)) and np.all(n>=0)
        checks['previous_rejected_stage_uses_extended']=probe_stats.get('guarded_extended_steps')==1
        result['evidence']['previous_rejected_stage']=dict(stats=probe_stats,
            witness_old_value=float(old_rejected[1+297]),witness_guarded_value=float(probe[1+297]),
            energy_change=float(base.measures(fixed,probe)['total']-base.measures(fixed,failed_input)['total']))
        np.savez_compressed(OUT/'isolated_stage_correction.npz',stage_input=failed_input,original_rejected=old_rejected,guarded_state=probe)

        record,states=base.save_case('transport',fixed,transport_initial,.02,20)
        checks['transport_complete_registered_case']=record['status']=='PASS'
        checks['transport_extended_branch_exercised']=record['integrator_stats']['guarded_extended_steps']>0
        result['evidence']['transport']=dict(energy_error=record['energy_error'],count_error=record['count_error'],
                                           stats=record['integrator_stats'],runtime_seconds=record['runtime_seconds'])

        escape=base.fixed_system(template,[.5],[empty],escape=.7)
        cell=base.ElectronicCell(template.catalog,.5,0.)
        p=cell.fermi_dirac(.2)
        n0=escape.bath_phonons+.02*np.exp(-escape.phonons.energies)
        expected=escape.bath_phonons+(n0-escape.bath_phonons)*np.exp(-.2/.7)
        errors=[];bytewise=[]
        for steps in (10,20,40):
            initial=escape.pack([.5],[p],[n0])
            record,states=base.save_case(f'escape_{steps}',escape,initial,.2,steps)
            with np.load(PREVIOUS/f'isolated_escape_{steps}.npz',allow_pickle=False) as old:
                bytewise.append(bool(states.dtype==old['states'].dtype and states.shape==old['states'].shape
                                     and states.tobytes()==old['states'].tobytes()))
            errors.append(float(np.dot(escape.phonons.capacities,abs(escape.unpack(states[-1])[2][0]-expected))))
        reductions=[errors[i]/errors[i+1] for i in (0,1)]
        checks['escape_all_stored_states_bytewise']=all(bytewise)
        checks['escape_third_order']=min(reductions)>=6.
        result['evidence']['escape']=dict(errors_capacity_L1=errors,reductions=reductions,bytewise_by_steps=bytewise)
        if not all(checks.values()):raise ValueError('prototype guard check failed')
        result['status']='PASS_GUARDED_PROTOTYPE_ONLY_NOT_GLOBALLY_ADMITTED'
    except Exception as exc:
        result['status']='FAIL_OR_INCOMPLETE';result['error']=repr(exc)
        raise
    finally:
        base.OUT,base.PLAN,base.integrate,base.write_new=old_out,old_plan,old_integrate,old_writer
        result['runtime_seconds']=time.perf_counter()-started
        write_new(OUT/'isolated_results.json',result)
        print(json.dumps(dict(status=result['status'],runtime_seconds=result['runtime_seconds'],checks=checks),default=scalar),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--register',action='store_true')
    args=parser.parse_args()
    if args.register:
        write_new(PLAN,registration());print('PREREGISTERED',PLAN)
    else:run()


if __name__=='__main__':main()
