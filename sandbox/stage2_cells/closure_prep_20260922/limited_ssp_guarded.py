"""Experimental extended-inventory fallback; no global numerical admission.

Ordinary steps delegate to the frozen SSP map without changing its arithmetic.
Only a positive particle/hole inventory near float64 underflow selects the
extended branch. Its COMMON event extents are accumulated before rounding the
population result to float64. There is no clipping or energy repair.
"""
from __future__ import annotations
from pathlib import Path
import sys
from types import FunctionType
import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'recovery_20260921'))
import limited_ssp as original
from pysnspd.experimental.cell_validation import ElectronicCell
from pysnspd.experimental.cell_transport import NativeTransportEvents

LD=np.longdouble
TRIGGER=LD(1024)*LD(np.finfo(np.float64).tiny)
EXTENDED=(np.finfo(LD).tiny<np.finfo(np.float64).tiny
          and np.finfo(LD).eps<np.finfo(np.float64).eps)


def require_extended():
    if not EXTENDED:
        raise RuntimeError('This Geminga prototype requires longdouble with genuinely wider exponent and mantissa than float64.')


def capability():
    return dict(extended=bool(EXTENDED),longdouble_bits=int(np.finfo(LD).bits),
                longdouble_mantissa_bits=int(np.finfo(LD).nmant),
                longdouble_min_exponent=int(np.finfo(LD).minexp),
                trigger_inventory=str(TRIGGER))


def needs_extended(system,state):
    _,p,n=system.unpack(state)
    ec=np.asarray(4*system.catalog.count_weights,dtype=LD)
    pc=np.asarray(system.phonons.capacities,dtype=LD)
    for values,capacity in ((p,ec),(1-p,ec),(n,pc)):
        inventory=np.asarray(values,dtype=LD)*capacity
        if np.any((inventory>0)&(inventory<TRIGGER)):
            return True
    return False


def _demands(indices,coefficients,rates,size):
    particle=np.zeros(size,dtype=LD);hole=np.zeros(size,dtype=LD)
    for ix,coefficient in zip(indices,coefficients):
        signed=np.asarray(coefficient,dtype=LD)*rates
        np.add.at(particle,ix,np.where(signed<0,-signed,LD(0)))
        np.add.at(hole,ix,np.where(signed>0,signed,LD(0)))
    return particle,hole


def _scatter(indices,coefficients,rates,size):
    result=np.zeros(size,dtype=LD)
    for ix,coefficient in zip(indices,coefficients):
        np.add.at(result,ix,np.asarray(coefficient,dtype=LD)*rates)
    return result


def _ratio(available,demand,h):
    result=np.ones_like(available,dtype=LD)
    active=demand>0
    result[active]=np.minimum(LD(1),LD(original.SAFETY)*available[active]/(h*demand[active]))
    if np.any(~np.isfinite(result)) or np.any(result<0):
        raise FloatingPointError('Invalid extended availability; no repair applied.')
    return result


def _limit(theta,indices,coefficients,rates,particle,hole=None):
    for ix,coefficient in zip(indices,coefficients):
        signed=np.asarray(coefficient,dtype=LD)*rates
        theta=np.minimum(theta,np.where(signed<0,particle[ix],LD(1)))
        if hole is not None:
            theta=np.minimum(theta,np.where(signed>0,hole[ix],LD(1)))
    return theta


def _electron_slots(event):
    if hasattr(event,'electron_lower_i'):
        ix=[event.electron_lower_i,event.electron_upper_i,event.electron_lower_j,event.electron_upper_j]
        sign=np.where(event.recombination,LD(-1),LD(1))
        coefficient=[sign*np.asarray(event.electron_beta_lower_i,dtype=LD),
                     sign*np.asarray(event.electron_beta_upper_i,dtype=LD),
                     -np.asarray(event.electron_beta_lower_j,dtype=LD),
                     -np.asarray(event.electron_beta_upper_j,dtype=LD)]
    else:
        ix=[event.index_i,event.index_j]
        coefficient=[np.where(event.recombination,LD(-1),LD(1)),-np.ones_like(event.omega,dtype=LD)]
    coefficient=[np.array(v,dtype=LD,copy=True) for v in coefficient]
    for k in range(1,len(ix)):
        for j in range(k):
            same=ix[k]==ix[j]
            coefficient[j]+=np.where(same,coefficient[k],LD(0))
            coefficient[k][same]=0
    return ix,coefficient


def _long_sum(stats,key,value):
    stats[key]=str(LD(stats.get(key,'0'))+LD(value))


def _flux(stats,rates,theta,energy,h,kind):
    # Preserve ordinary diagnostic fields, and retain underflow-scale information
    # as decimal strings rather than forcing it into an unrepresentable float64.
    original._record_flux(stats,rates,theta,energy,h,kind)
    weight=LD(stats.get('_stage_weight',1.))
    defect=abs(rates)*(1-theta)
    _long_sum(stats,'guarded_extent_defect_longdouble',weight*h*np.sum(defect,dtype=LD))
    _long_sum(stats,'guarded_energy_flux_defect_longdouble',weight*h*np.dot(np.asarray(energy,dtype=LD),defect))
    if len(theta):
        old=LD(stats.get('guarded_minimum_event_factor_longdouble','1'))
        stats['guarded_minimum_event_factor_longdouble']=str(min(old,theta.min()))


def _extended_step(system,t,z,h,stats):
    h=LD(h)
    amplitudes,p,n=system.unpack(z)
    original._add(stats,'forward_euler_stages',1)
    system.rhs_calls+=1
    cells=[ElectronicCell(system.catalog,amplitudes[i],system.gammas[i]) for i in range(system.cell_count)]
    capacity=np.asarray([4*c.weights for c in cells],dtype=LD)
    phcap=np.asarray(system.phonons.capacities,dtype=LD)
    base=np.asarray(z,dtype=LD).copy()
    base_ledger=base[system.population_size:-1].reshape(system.cell_count,4)
    pbase=np.empty_like(p,dtype=LD);nbase=np.empty_like(n,dtype=LD)
    weight=LD(stats.get('_stage_weight',1.))
    for i,cell in enumerate(cells):
        temperature=cell.equivalent_temperature(p[i]);force=cell.moments(p[i])[1]
        motion=system.mobility.amplitude_response(amplitudes[i],force,temperature)
        heating=cell.heating(p[i],system.external_powers[i]+motion.heat,system.bath_temperature)
        fraction=h/LD(system.tau_kin);escape_fraction=h/LD(system.tau_escape)
        if fraction>1 or escape_fraction>1:
            raise ValueError('regular baseline exceeds its convex interval')
        pbase[i]=(1-fraction)*np.asarray(p[i],dtype=LD)+fraction*np.asarray(cell.fermi_dirac(temperature),dtype=LD)+h*np.asarray(heating,dtype=LD)
        nbase[i]=(1-escape_fraction)*np.asarray(n[i],dtype=LD)+escape_fraction*np.asarray(system.bath_phonons,dtype=LD)
        escaping=(np.asarray(n[i],dtype=LD)-np.asarray(system.bath_phonons,dtype=LD))/LD(system.tau_escape)
        start=i*system.block_size
        base[start]=LD(amplitudes[i])+h*LD(motion.velocity)
        base[start+1:start+1+system.electron_size]=pbase[i]
        base[start+1+system.electron_size:start+system.block_size]=nbase[i]
        power=np.dot(phcap*np.asarray(system.phonons.energies,dtype=LD),escaping)
        base_ledger[i]+=h*np.asarray([power,system.external_powers[i],0.,motion.heat],dtype=LD)
        original._add(stats,'integrated_external_input',float(weight*h*LD(system.external_powers[i])))
        original._add(stats,'integrated_condensate_heat',float(weight*h*LD(motion.heat)))
    if (np.any(~np.isfinite(base)) or np.any(pbase<0) or np.any(pbase>1) or np.any(nbase<0)
            or np.any(base[0:system.population_size:system.block_size]<0)):
        raise ValueError('extended regular baseline outside physical domain')
    face=None;face_rates=None
    face_particle=np.zeros_like(p,dtype=LD);face_hole=np.zeros_like(p,dtype=LD)
    if system.cell_count==2:
        face=NativeTransportEvents(*cells,system.diffusion_over_length_squared,system.face_order)
        face_rates=np.asarray(face.rates(p),dtype=LD)
        for side in (0,1):
            ix=[face.indices[side,:,k] for k in (0,1)]
            coefficient=[LD(2*side-1)*np.asarray(face.barycentric[side,:,k],dtype=LD) for k in (0,1)]
            face_particle[side],face_hole[side]=_demands(ix,coefficient,face_rates,system.electron_size)
    ratios_p=[];ratios_h=[]
    result=base.copy();ledger=result[system.population_size:-1].reshape(system.cell_count,4)
    for i,cell in enumerate(cells):
        event=system.reaction_events(cell)
        rates=np.asarray(event.rates(p[i],n[i]),dtype=LD)
        ix,coefficient=_electron_slots(event)
        demand_p,demand_h=_demands(ix,coefficient,rates,system.electron_size)
        pix=[event.phonon_lower,event.phonon_upper]
        pc=[np.asarray(event.beta_lower,dtype=LD),np.asarray(event.beta_upper,dtype=LD)]
        demand_n,_=_demands(pix,pc,rates,system.phonon_size)
        rp=_ratio(capacity[i]*pbase[i],demand_p+face_particle[i],h)
        rh=_ratio(capacity[i]*(1-pbase[i]),demand_h+face_hole[i],h)
        rn=_ratio(phcap*nbase[i],demand_n,h)
        ratios_p.append(rp);ratios_h.append(rh)
        theta=_limit(np.ones_like(rates,dtype=LD),ix,coefficient,rates,rp,rh)
        theta=_limit(theta,pix,pc,rates,rn)
        actual=theta*rates
        start=i*system.block_size
        result[start+1:start+1+system.electron_size]+=h*_scatter(ix,coefficient,actual,system.electron_size)/capacity[i]
        result[start+1+system.electron_size:start+system.block_size]+=h*_scatter(pix,pc,actual,system.phonon_size)/phcap
        ledger[i,2]+=h*np.dot(np.asarray(event.omega,dtype=LD),actual)
        _flux(stats,rates,theta,event.omega,h,'reaction')
        del event,rates,ix,coefficient,pix,pc,theta,actual
    if face is not None:
        theta=np.ones_like(face_rates,dtype=LD)
        for side in (0,1):
            ix=[face.indices[side,:,k] for k in (0,1)]
            coefficient=[LD(2*side-1)*np.asarray(face.barycentric[side,:,k],dtype=LD) for k in (0,1)]
            theta=_limit(theta,ix,coefficient,face_rates,ratios_p[side],ratios_h[side])
        actual=theta*face_rates
        for side in (0,1):
            ix=[face.indices[side,:,k] for k in (0,1)]
            coefficient=[LD(2*side-1)*np.asarray(face.barycentric[side,:,k],dtype=LD) for k in (0,1)]
            start=side*system.block_size+1
            result[start:start+system.electron_size]+=h*_scatter(ix,coefficient,actual,system.electron_size)/capacity[side]
        result[-1]+=h*np.dot(np.asarray(face.energies,dtype=LD),actual)
        _flux(stats,face_rates,theta,face.energies,h,'transport')
    # Validate extended populations before rounding; rounding is the ordinary
    # dtype conversion, never a sign correction or projection onto support.
    blocks=result[:system.population_size].reshape(system.cell_count,system.block_size)
    ep=blocks[:,1:1+system.electron_size];pn=blocks[:,1+system.electron_size:]
    if np.any(~np.isfinite(result)) or np.any(ep<0) or np.any(ep>1) or np.any(pn<0):
        raise FloatingPointError('extended event update left support; no repair applied')
    obtained=np.asarray(result,dtype=np.float64)
    system.unpack(obtained)
    original._add(stats,'guarded_extended_steps',1)
    original._add(stats,'guarded_positive_values_rounded_to_zero',int(np.count_nonzero((result>0)&(obtained==0))))
    return obtained


def step(system,t,z,h,stats=None):
    require_extended()
    stats={} if stats is None else stats
    if not np.isfinite(h) or h<=0:
        raise ValueError('positive finite step required')
    if not needs_extended(system,z):
        original._add(stats,'guarded_delegated_steps',1)
        return original.step(system,t,z,h,stats)
    return _extended_step(system,t,z,h,stats)


# Same bytecode/defaults as the original integrator, with a PRIVATE global
# dictionary. The frozen module and every other caller remain unchanged.
_globals=dict(original.integrate.__globals__)
_globals['step']=step
_integrate=FunctionType(original.integrate.__code__,_globals,
                        name='integrate_with_guarded_step',argdefs=original.integrate.__defaults__)
_integrate.__kwdefaults__=dict(original.integrate.__kwdefaults__)


def integrate(system,initial,duration,steps,*,callback=None,stats=None):
    require_extended()
    times,states,record=_integrate(system,initial,duration,steps,callback=callback,stats=stats)
    record['method']='experimental conservative-event-limited SSPRK3 with subnormal arithmetic guard'
    for name in ('guarded_delegated_steps','guarded_extended_steps','guarded_positive_values_rounded_to_zero'):
        record.setdefault(name,0)
    for name in ('guarded_extent_defect_longdouble','guarded_energy_flux_defect_longdouble'):
        record.setdefault(name,'0')
    record.setdefault('guarded_minimum_event_factor_longdouble','1')
    record['guarded_prototype']=True
    record['guarded_capability']=capability()
    record['guarded_scope']='Subnormal-inventory correction prototype; requires renewed time/mesh validation unless trajectory equivalence is measured.'
    return times,states,record
