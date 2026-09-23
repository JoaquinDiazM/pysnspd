"""Observe a rejected transport stage and separately test smooth SSP escape.

The observer does not change inputs, outputs, floating point operations, fluxes
or exceptions of the frozen map. A reproduced rejection is never an admission.
"""
from pathlib import Path
import argparse
from decimal import Decimal, localcontext
import json
import time
import numpy as np
import isolated_checks as base
import limited_ssp as ssp
from pysnspd.experimental.cell_transport import NativeTransportEvents

OUT=base.OUT/'isolated_transport_diagnostic'
PLAN=OUT/'isolated_registration.json'


def sha(p):
    return base.sha(p)


def scalar(value):
    if isinstance(value,np.generic):
        return value.item()
    raise TypeError(type(value).__name__)


def write_new(path,value):
    encoded=json.dumps(value,indent=2,allow_nan=False,default=scalar)+'\n'
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x',encoding='utf-8') as f:
        f.write(encoded)


def finite_or_label(value):
    value=float(value)
    return value if np.isfinite(value) else str(value)


def sources():
    result=base.sources()
    result[Path(__file__).relative_to(base.ROOT).as_posix()]=sha(__file__)
    return result


def registration():
    return dict(schema='pysnspd.stage2.isolated-observational-diagnosis.v1',
        status='PREREGISTERED_NOT_EXECUTED',source_hashes=sources(),
        failed_result_sha256=sha(base.ROOT/'docs/implementation/stage2/closure_prep_20260922/isolated_recovery/isolated_results.json'),
        transport=dict(amplitudes=[.5,1.],gammas=[0.,0.],electron_states=630,phonon_nodes=1025,
                       duration=.02,steps=20,filled_count_band=[1.1,1.5],
                       populations='left band exactly1, left exterior0; right electrons0; all phonons0',
                       BGK=False,condensate_motion=False,heating=False,reactions=False,escape=False,
                       diffusion_over_length_squared=.2,face_order=2),
        purpose='Capture the first rejected state and donor inventories; no bypass, repair or acceptance.',
        escape=dict(duration=.2,tau=.7,steps=[10,20,40],reference='exact affine exponential',
                    expected_order=3,minimum_error_reduction=6.),
        execution=dict(wall_timeout_seconds=240,expected_seconds='under30',threads=1,max_attempts=1),
        reuses='No approved reaction trajectory is repeated. Only the known failing transport is observed, as explicitly authorized.',
        no_map_or_kernel_changes=True)


def observe_transport(template,empty):
    started=time.perf_counter()
    system=base.fixed_system(template,[.5,1.],[empty,empty],transport=True)
    band=((template.catalog.count_nodes>1.1)&(template.catalog.count_nodes<1.5)).astype(float)
    initial=system.pack([.5,1.],[band,np.zeros_like(band)],np.zeros((2,1025)))
    capture={}
    calls=0
    actual_step=ssp.step
    actual_unpack=system.unpack

    def unpack_observer(state):
        try:
            return actual_unpack(state)
        except ValueError:
            if 'rejected_state' not in capture:
                capture['rejected_state']=np.array(state,copy=True)
            raise

    def step_observer(sys,t,z,h,stats=None):
        nonlocal calls
        calls+=1
        capture['stage_input']=np.array(z,copy=True)
        capture['stage_time']=float(t)
        capture['step']=float(h)
        capture['FE_call']=calls
        capture['macrostep']=(calls-1)//3+1
        capture['stage_in_macrostep']=(calls-1)%3+1
        return actual_step(sys,t,z,h,stats)

    system.unpack=unpack_observer
    ssp.step=step_observer
    stats={}
    try:
        ssp.integrate(system,initial,.02,20,stats=stats)
        status='KNOWN_REJECTION_NOT_REPRODUCED_NO_ADMISSION'
    except ValueError as exc:
        if 'rejected_state' not in capture:
            raise
        status='REPRODUCED_SUPPORT_REJECTION_NOT_ADMITTED'
        capture['exception']=repr(exc)
    finally:
        ssp.step=actual_step
        system.unpack=actual_unpack
    result=dict(schema='pysnspd.stage2.isolated-transport-diagnosis.v1',status=status,
                registration_sha256=sha(PLAN),source_hashes=sources(),stats_before_rejection=stats)
    if 'rejected_state' in capture:
        before,after=capture.pop('stage_input'),capture.pop('rejected_state')
        a,p,n=system.unpack(before)
        cells=[base.ElectronicCell(system.catalog,a[i],0.) for i in range(2)]
        face=NativeTransportEvents(*cells,.2,2)
        rates=face.rates(p)
        capacity=4*system.catalog.count_weights
        h=capture['step']
        rows=[]; diagnostics={}; particle=[];holes=[];indices=[];coefficients=[]
        for side in (0,1):
            ix=[face.indices[side,:,k] for k in (0,1)]
            coefficient=[(2*side-1)*face.barycentric[side,:,k] for k in (0,1)]
            demand_p,demand_h=ssp._demands(ix,coefficient,rates,system.electron_size)
            with np.errstate(over='ignore',divide='ignore',invalid='ignore',under='ignore'):
                rp=ssp._ratio(capacity*p[side],demand_p,h)
                rh=ssp._ratio(capacity*(1-p[side]),demand_h,h)
            indices.append(ix);coefficients.append(coefficient);particle.append(rp);holes.append(rh)
            diagnostics[f'demand_particle_{side}']=demand_p
            diagnostics[f'demand_hole_{side}']=demand_h
            diagnostics[f'available_particle_{side}']=capacity*p[side]
            diagnostics[f'available_hole_{side}']=capacity*(1-p[side])
        theta=np.ones_like(rates)
        for side in (0,1):
            theta=ssp._limit(theta,indices[side],coefficients[side],rates,particle[side],holes[side])
        extent=theta*rates
        for side in (0,1):
            start=side*system.block_size+1
            obtained=after[start:start+system.electron_size]
            dp=h*ssp._scatter(indices[side],coefficients[side],extent,system.electron_size)/capacity
            for index in np.flatnonzero((obtained<0)|(obtained>1)):
                index=int(index)
                incident=(face.indices[side,:,0]==index)|(face.indices[side,:,1]==index)
                avail=diagnostics[f'available_particle_{side}'][index]
                demand=diagnostics[f'demand_particle_{side}'][index]
                with np.errstate(all='ignore'):
                    raw=np.float64(avail)/np.float64(h*demand)
                with localcontext() as context:
                    context.prec=1200
                    exact_add=Decimal.from_float(float(p[side,index]))+Decimal.from_float(float(dp[index]))
                    exact_inventory=Decimal.from_float(float(capacity[index]))*Decimal.from_float(float(p[side,index]))
                rows.append(dict(cell=side,electron_index=index,count_coordinate=float(system.catalog.count_nodes[index]),
                    p_before=float(p[side,index]),p_rejected=float(obtained[index]),capacity=float(capacity[index]),
                    particle_inventory=float(avail),particle_inventory_exact_decimal=str(exact_inventory),
                    h_times_particle_demand=float(h*demand),particle_demand=float(demand),
                    raw_ratio=finite_or_label(raw),particle_ratio=float(particle[side][index]),
                    hole_inventory=float(diagnostics[f'available_hole_{side}'][index]),
                    hole_demand=float(diagnostics[f'demand_hole_{side}'][index]),hole_ratio=float(holes[side][index]),
                    min_incident_factor=float(theta[incident].min()) if np.any(incident) else None,
                    increment=float(dp[index]),exact_add_of_float_inputs=str(exact_add),
                    reproduced_float_sum=float(p[side,index]+dp[index]),
                    below_smallest_normal=bool(abs(obtained[index])<np.finfo(float).tiny)))
        ph_bad=[]
        for side in (0,1):
            start=side*system.block_size+1+system.electron_size
            ph=after[start:start+system.phonon_size]
            ph_bad.extend(dict(cell=side,index=int(k),value=float(ph[k])) for k in np.flatnonzero(ph<0))
        target=OUT/'isolated_transport_rejected_state.npz'
        if target.exists():raise FileExistsError(target)
        np.savez_compressed(target,initial_state=initial,stage_input=before,rejected_state=after,
                            electron_count=system.catalog.count_nodes,electron_weights=system.catalog.count_weights,
                            face_rates=rates,face_theta=theta,face_indices=face.indices,
                            face_barycentric=face.barycentric,**diagnostics)
        result.update(capture)
        result.update(rejected_electrons=rows,rejected_phonons=ph_bad,
                      state_sha256=sha(target),observation_only=True,population_clipping=False)
    result['runtime_seconds']=time.perf_counter()-started
    write_new(OUT/'isolated_transport_diagnosis.json',result)
    print(json.dumps({k:result[k] for k in ('status','runtime_seconds')},default=scalar),flush=True)
    return result


def escape_order(template,empty):
    # This independent operator was never reached by either previous batch.
    original_out,original_plan,original_writer=base.OUT,base.PLAN,base.write_new
    base.OUT,base.PLAN,base.write_new=OUT,PLAN,write_new
    started=time.perf_counter()
    try:
        system=base.fixed_system(template,[.5],[empty],escape=.7)
        cell=base.ElectronicCell(system.catalog,.5,0.)
        p=cell.fermi_dirac(.2)
        n0=system.bath_phonons+.02*np.exp(-system.phonons.energies)
        reference=system.bath_phonons+(n0-system.bath_phonons)*np.exp(-.2/.7)
        errors=[]
        for steps in (10,20,40):
            initial=system.pack([.5],[p],[n0])
            record,states=base.save_case(f'escape_{steps}',system,initial,.2,steps)
            errors.append(float(np.dot(system.phonons.capacities,abs(system.unpack(states[-1])[2][0]-reference))))
        ratios=[errors[i]/errors[i+1] for i in (0,1)]
        result=dict(schema='pysnspd.stage2.ssp-escape-order.v1',status='PASS' if min(ratios)>=6. else 'FAIL',
                    errors_capacity_L1=errors,reductions=ratios,minimum_reduction=6.,
                    registration_sha256=sha(PLAN),source_hashes=sources(),runtime_seconds=time.perf_counter()-started,
                    scope='Smooth isolated escape only; no transport or complete stage2 admission.')
        write_new(OUT/'isolated_escape_order.json',result)
        print(json.dumps(result),flush=True)
        if result['status']!='PASS':raise ValueError('SSP escape order failed')
    finally:
        base.OUT,base.PLAN,base.write_new=original_out,original_plan,original_writer


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--register',action='store_true')
    args=parser.parse_args()
    if args.register:
        write_new(PLAN,registration()); print('PREREGISTERED',PLAN); return
    if json.loads(PLAN.read_text(encoding='utf-8'))!=registration():raise ValueError('source/preregistration mismatch')
    if list(OUT.glob('*.npz')) or (OUT/'isolated_transport_diagnosis.json').exists():raise FileExistsError('diagnosis already attempted')
    template,_,_=base.setup('one',heating=0.,escape=np.inf)
    parent=template.reaction_events(base.ElectronicCell(template.catalog,.5,0.))
    empty=base.DirectionalEvents(parent,np.zeros(len(parent.omega),bool),'disabled')
    del parent
    observe_transport(template,empty)
    escape_order(template,empty)


if __name__=='__main__':
    main()
