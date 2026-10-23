"""Historical K20/Allmaras fit inputs in the current v1.2.3 native equations.

Isolated execution adapter. Does not edit imported files or change current runs.
Fixed D replaces the depairing fit. Gamma60 independently rescales the phonon
DOS, and an exactly integrated 40nm hotbelt replaces Gaussian deposition.
This is neither the historical 1D solver nor a reproduction of its IRF fit.
"""
from __future__ import annotations
import argparse
import copy
import importlib.util
import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'

REFERENCE = dict(reference_id='K20_Supplementary_Table1_A20_pp95_98',
    source_url='https://eprints.lancs.ac.uk/id/eprint/140252/3/Binder1.pdf#page=33',
    thesis_url='https://thesis.caltech.edu/13748/08/Allmaras_Thesis_Final.pdf#page=107',
    sheet_resistance_ohm=608., diffusion_m2_s=5e-5, gamma=60., retained_fraction=.667,
    tau_ee_Tc_ps=6.,tau_ep_Tc_ps=24.7,tau_escape_ps=20.,bath_K=2.,
    hotbelt_length_nm=40.,width_nm=80.,film_thickness_nm=7.,Tc_K=8.65,
    reported_Fano_sigma_1550_meV=92.,reported_nonuniformity_sigma_1550_meV=40.,
    fluctuations_sampled=False,
    scope='Deterministic mean-energy inputs; current 2D BCS/KWT/native kinetics and thesis circuit retained; not an exact historical solver or fitted-IRF reproduction')


def clip_polygon(poly, normal, offset):
    import numpy as np
    if len(poly)<3:
        return np.empty((0,2))
    result=[]
    for start,end in zip(poly,np.roll(poly,-1,axis=0)):
        a=float(start@normal-offset);b=float(end@normal-offset)
        ina=a<=0.;inb=b<=0.
        if ina:result.append(start)
        if ina!=inb:result.append(start+(end-start)*(a/(a-b)))
    return np.asarray(result).reshape(-1,2)


def hotbelt_weights(geometry, volumes, thickness):
    """Exact intersection of clipped Voronoi cells with the 40nm rectangle."""
    import numpy as np
    xy=geometry.coordinates_m*1e9
    center=.5*(xy[:,0].min()+xy[:,0].max());W=geometry.width_m*1e9
    adjacency=[[] for _ in xy]
    for i,j in geometry.graph.edges:adjacency[i].append(j);adjacency[j].append(i)
    areas=[]
    for i,site in enumerate(xy):
        poly=np.array([[center-20,-W/2],[center+20,-W/2],[center+20,W/2],[center-20,W/2]])
        for j in adjacency[i]:
            poly=clip_polygon(poly,xy[j]-site,.5*float((xy[j]+site)@(xy[j]-site)))
            if len(poly)<3:break
        if len(poly)<3:areas.append(0.);continue
        local=poly-site
        areas.append(.5*abs(local[:,0]@np.roll(local[:,1],-1)-local[:,1]@np.roll(local[:,0],-1))*1e-18)
    active_volume=np.asarray(areas)*thickness
    exact=40e-9*geometry.width_m*thickness
    if not np.isclose(active_volume.sum(),exact,rtol=2e-10,atol=0):
        raise ArithmeticError('Hotbelt cell intersections do not cover the exact 40nm volume')
    if np.any(active_volume>np.asarray(volumes)*(1+2e-9)):
        raise ArithmeticError('Hotbelt intersection exceeds a native dual-cell volume')
    cuts=np.r_[geometry.left_nodes,geometry.right_nodes]
    if np.any(active_volume[cuts]>0):raise ArithmeticError('Hotbelt touches a longitudinal cut')
    return active_volume/exact


def configure(source_runner, threads):
    os.environ['NUMBA_NUM_THREADS']=str(threads)
    path=Path(source_runner).resolve()
    if path != (ROOT/'pipelines/13_dimensional_photon_run.py').resolve():
        raise ValueError('The current v1.2.3 producer is required')
    spec=importlib.util.spec_from_file_location('_allmaras_current_native',path)
    R=importlib.util.module_from_spec(spec);spec.loader.exec_module(R)
    import numpy as np
    from scipy.constants import Boltzmann as kB, Planck, speed_of_light
    from pysnspd.experimental.effective_material import EffectivePhononMaterial
    original_calibration=R.OLD.calibration_result
    original_material=R.B.allmaras_acoustic_material
    original_build=R.OLD.build_system
    original_pins=R.source_pins
    original_write=R.write_json

    def fixed_calibration(**kwargs):
        if kwargs['sheet_resistance_ohm']!=608.:raise ValueError('Historical control requires Rs608')
        c=original_calibration(**kwargs)
        ref=c['reference'];D=ref['diffusion_m2_s'];dep=ref['depairing_current_A']
        if D!=5e-5:raise ArithmeticError('Historical diffusion reference changed')
        c.update(diffusion_m2_s=D,N0_per_J_m3=ref['N0_per_J_m3'],ell0_m=ref['ell0_m'],
            target_depairing_current_A=dep,achieved_depairing_current_A=dep,
            relative_target_residual=0.,requested_current_fraction=kwargs['requested_current_A']/dep,
            requested_current_below_depairing=kwargs['requested_current_A']<dep,
            historical_parameter_control=True,measured_depairing_fit_enforced=False)
        c['calibration'].update(method='Historical fixed-D input; target field stores predicted bulk for native consistency checks only',
            diffusion_is_fitted_effective_parameter=False,
            sheet_resistance_status='Historical K20 simulation input, not an independent wire measurement')
        c['assumptions'].update(target='No measured-depairing fit; report the current closure prediction',
            scope=REFERENCE['scope'])
        return c

    def gamma_material(**kwargs):
        material=original_material(**kwargs)
        info=copy.deepcopy(material.manifest)
        old=float(info['gamma_A20']);scale=old/60.
        for key in ('g_ph_coefficient_per_J3_m3','acoustic_B_J_per_m3_K4',
                    'target_capacity_at_reference_J_per_m3_K','heat_capacity_at_reference_J_per_m3_K',
                    'total_modes_per_m3'):
            if key in info:info[key]*=scale
        info.update(material_id='K20_historical_gamma60_control',status='HISTORICAL_PARAMETER_CONTROL',
            gamma_A20=60.,gamma_independent=60.,gamma_before_override=old,
            historical_gamma_DOS_multiplier=scale,
            acoustic_speed_alone_determines_DOS=False,
            limitation='Gamma60 overrides the acoustic DOS amplitude independently of nominal sound speed; not an acoustic-consistent material inference',
            historical_reference=REFERENCE)
        return EffectivePhononMaterial(material.energy_J,material.alpha2F,
            material.g_ph_per_J_m3*scale,material.N0_per_J_m3,info)

    def prepare(functional,kinetics,gap,p0,geometry,*,wavelength_nm,sigma_nm,retained_fraction,bath_temperature_K):
        count=functional.graph.n_nodes;bath=np.full(count,bath_temperature_K)
        n0=R.B.bose_matrix(kinetics.phonons.energies,functional.Eu_J,bath)
        baseline=functional.evaluate(gap,p0,temperature_K=bath).node_energy_J
        baseline=baseline+np.sum(kinetics.phonon_energy_weights_J*n0,axis=1)
        weights=hotbelt_weights(geometry,functional.volumes_m3,functional.thickness_m)
        photon=Planck*speed_of_light/(wavelength_nm*1e-9);retained=photon*retained_fraction
        target=retained*weights;active=target>0;minimum_ratio=float('inf')
        def values(T):
            nonlocal minimum_ratio
            p=functional.fermi_dirac(gap,T);n=R.B.bose_matrix(kinetics.phonons.energies,functional.Eu_J,T)
            thermal=functional.thermal_data(gap,p,T)
            ce=thermal['capacity_effective_J_K'];cfd=thermal['capacity_FD_J_K']
            cph=np.sum(kinetics.phonon_energy_weights_J*n*(n+1)*
                kinetics.phonons.energies[None,:]*functional.Eu_J/(kB*T[:,None]**2),axis=1)
            if np.any(~np.isfinite(ce+cph)) or np.any(ce<=0) or np.any(cfd<=0) or np.any(ce+cph<=0):
                raise ArithmeticError('Historical preparation left the positive-capacity domain')
            minimum_ratio=min(minimum_ratio,float(np.min(ce/cfd)))
            energy=functional.evaluate(gap,p,temperature_K=T).node_energy_J
            return energy+np.sum(kinetics.phonon_energy_weights_J*n,axis=1)-baseline,p,n
        values(bath);lo=bath.copy();hi=np.where(active,np.maximum(2*bath,2.),bath)
        for _ in range(32):
            added,_,_=values(hi);need=(added<target)&active
            if not np.any(need):break
            hi[need]*=2
        else:raise ArithmeticError('Historical hotbelt energy bracket failed')
        for _ in range(60):
            mid=(lo+hi)/2;added,_,_=values(mid);below=added<target
            lo=np.where(below&active,mid,lo);hi=np.where((~below)&active,mid,hi)
        T=np.where(active,(lo+hi)/2,bath);added,p,n=values(T);defect=added-target
        if max(abs(float(defect.sum())),float(np.max(abs(defect))))>1e-10*retained:
            raise ArithmeticError('Historical hotbelt energy normalization failed')
        center=[float(.5*(geometry.coordinates_m[:,0].min()+geometry.coordinates_m[:,0].max())),0.]
        receipt=dict(wavelength_nm=wavelength_nm,photon_energy_J=photon,retained_fraction=retained_fraction,
            retained_energy_J=retained,profile='uniform_hotbelt',hotbelt_length_nm=40.,
            sigma_nm=40/np.sqrt(12),sigma_interpretation='Longitudinal RMS of ideal40nm top-hat; not Gaussian',
            center_m=center,initialization='Historical common FD/BE 40nm transverse hotbelt; current native DC gap retained',
            time_origin='Effective postcascade transfer; no extra fitted clock shift',
            energy_functional='Current U_BCS+G-sum(Te*dG/dTe)+phonons at fixed initial DC gap',
            energy_deposited_J=float(added.sum()),energy_defect_J=float(defect.sum()),
            relative_total_energy_defect=float(abs(defect.sum())/retained),
            maximum_initial_temperature_K=float(T.max()),minimum_preparation_effective_capacity_ratio=minimum_ratio,
            positive_capacity_checked_throughout_preparation=True,capacity_floor_used=False,
            source_must_not_be_applied_again=True,cut_deposition_J=0.,historical_reference=REFERENCE)
        return p,n,T,target,receipt

    def build(args,output,*,resume=False):
        result=original_build(args,output,resume=resume)
        result[-1].update(historical_parameter_control=True,historical_reference=REFERENCE,
            diffusion_is_fitted=False,gamma_override=60.,
            calibration_scope='Historical fixed D; target field is current-model predicted depairing, not a fit to measured depairing')
        return result

    def pins():
        values=original_pins();values[str(Path(__file__).relative_to(ROOT))]=R.B.sha256(__file__)
        return values
    def write(path,value):
        if Path(path).name=='manifest.json':
            value.update(historical_parameter_control=True,historical_reference=REFERENCE)
        return original_write(path,value)
    R.OLD.calibration_result=fixed_calibration
    R.B.allmaras_acoustic_material=gamma_material
    R.OLD.prepare_photon=prepare;R.OLD.build_system=build
    R.source_pins=pins;R.write_json=write
    return R


def main(argv=None):
    early=argparse.ArgumentParser(add_help=False)
    early.add_argument('--source-runner',type=Path,default=ROOT/'pipelines/13_dimensional_photon_run.py')
    early.add_argument('--threads',type=int,default=int(os.environ.get('NUMBA_NUM_THREADS','3')))
    options,remaining=early.parse_known_args(argv)
    R=configure(options.source_runner,options.threads)
    p=R.parser();p.description=__doc__
    p.set_defaults(sheet_resistance_ohm=608.,retained_fraction=.667,tau_ee_Tc_ps=6.,
        tau_ep_Tc_ps=24.7,tau_escape_ps=20.,bath_K=2.,sigma_nm=40/(12**.5),latency_only=True)
    a=p.parse_args(remaining+['--threads',str(options.threads)])
    for key in ('sheet_resistance_ohm','retained_fraction','tau_ee_Tc_ps','tau_ep_Tc_ps','tau_escape_ps','bath_K'):
        if getattr(a,key)!=REFERENCE[key]:raise ValueError('Historical fixed input changed: '+key)
    if not a.latency_only or a.recover_from_run is not None:raise ValueError('Historical control requires fresh/own-resumed latency-only runs')
    c=R.OLD.calibration_result(bulk_factory=R.ResummedThermalFunctional,sheet_resistance_ohm=608.,
        target_depairing_current_A=38.8*80/120*1e-6,width_m=80e-9,requested_current_A=a.current_uA*1e-6,
        bath_temperature_K=2.,critical_temperature_K=8.65,count_cutoff_bar=a.xi_cutoff,count_order=a.xi_order)
    a.target_depairing_uA=c['target_depairing_current_A']*1e6
    a.historical_parameter_control=True;a.historical_gamma=60.;a.initial_profile='uniform_hotbelt';a.hotbelt_length_nm=40.
    try:return R.run(a)
    except Exception as exc:
        path=Path(a.output_dir)/'manifest.json'
        if path.exists():
            m=json.loads(path.read_text())
            if m.get('status')=='PREPARATION':
                m.update(status='PREPARATION_FAILURE',error=type(exc).__name__+': '+str(exc));R.write_json(path,m)
        raise


if __name__=='__main__':raise SystemExit(main())
