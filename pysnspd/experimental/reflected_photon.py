"""One retained photon with a reflecting transverse postcascade preparation.

The Neumann image kernel is an explicit preparation hypothesis, not a second
photon, a dynamical source, or proof of microscopic sidewall reflectivity.
"""
import math
import numpy as np
from scipy.constants import Boltzmann as kB, Planck, speed_of_light, elementary_charge


def reflecting_shape(xy,width_m,sigma_m,y0_m):
    xy=np.asarray(xy,float)
    if (not np.all(np.isfinite(xy)) or not np.isfinite(width_m+sigma_m+y0_m)
        or min(width_m,sigma_m)<=0 or abs(y0_m)>width_m/2):
        raise ValueError('Finite positive width/sigma and impact inside the transverse strip required')
    u=xy[:,1]+width_m/2;u0=y0_m+width_m/2
    image_radius=max(2,int(math.ceil(9*sigma_m/(2*width_m)))+1)
    transverse=np.zeros(len(xy))
    for k in range(-image_radius,image_radius+1):
        transverse+=np.exp(-.5*((u-u0+2*k*width_m)/sigma_m)**2)
        transverse+=np.exp(-.5*((u+u0+2*k*width_m)/sigma_m)**2)
    x0=.5*(xy[:,0].min()+xy[:,0].max())
    shape=np.exp(-.5*((xy[:,0]-x0)/sigma_m)**2)*transverse
    return shape,dict(profile='reflecting_gaussian_images',center_m=[float(x0),float(y0_m)],
        impact_y_nm=float(y0_m*1e9),image_radius=image_radius,
        transverse_kernel='sum_k[g(u-u0+2kw)+g(u+u0+2kw)], u=y+w/2',
        boundary_hypothesis='Perfectly reflecting transverse postcascade preparation (Neumann); substrate escape unchanged',
        microscopic_sidewall_reflectivity_measured=False,virtual_images_are_additional_photons=False,
        continuum_edge_limit='Coincident primary and image: doubled interior half-Gaussian at abs(y0)=w/2',
        source_reference='https://doi.org/10.1007/s00020-023-02749-0')


def prepare_reflected(functional,kinetics,gap,p0,geometry,*,wavelength_nm,
        sigma_nm,retained_fraction,bath_temperature_K,bose_matrix,impact_y_nm=0.,
        deposited_energy_eV=None):
    bath=np.full(functional.graph.n_nodes,bath_temperature_K)
    n0=bose_matrix(kinetics.phonons.energies,functional.Eu_J,bath)
    baseline=functional.evaluate(gap,p0,temperature_K=bath).node_energy_J
    baseline=baseline+np.sum(kinetics.phonon_energy_weights_J*n0,axis=1)
    shape,receipt=reflecting_shape(geometry.coordinates_m,geometry.width_m,sigma_nm*1e-9,impact_y_nm*1e-9)
    cuts=np.r_[geometry.left_nodes,geometry.right_nodes];shape[cuts]=0.
    weighted=shape*functional.volumes_m3
    if not np.isfinite(weighted.sum()) or weighted.sum()<=0:raise ValueError('No finite interior photon support')
    photon=Planck*speed_of_light/(wavelength_nm*1e-9)
    retained=photon*retained_fraction if deposited_energy_eV is None else deposited_energy_eV*elementary_charge
    if not np.isfinite(retained) or not 0<retained<=photon:raise ValueError('Require 0 < retained photon energy <= hc/lambda')
    target=retained*weighted/weighted.sum();active=target>0;minimum_ratio=math.inf
    def values(T):
        nonlocal minimum_ratio
        p=functional.fermi_dirac(gap,T);n=bose_matrix(kinetics.phonons.energies,functional.Eu_J,T)
        thermal=functional.thermal_data(gap,p,T)
        ce=thermal['capacity_effective_J_K'];cfd=thermal['capacity_FD_J_K']
        cph=np.sum(kinetics.phonon_energy_weights_J*n*(n+1)*
            (kinetics.phonons.energies[None,:]*functional.Eu_J)/(kB*T[:,None]**2),axis=1)
        if np.any(~np.isfinite(ce+cph)) or np.any(ce<=0) or np.any(cfd<=0) or np.any(ce+cph<=0):
            raise ArithmeticError('Reflected preparation left positive-capacity domain; no floor')
        minimum_ratio=min(minimum_ratio,float(np.min(ce/cfd)))
        energy=functional.evaluate(gap,p,temperature_K=T).node_energy_J
        return energy+np.sum(kinetics.phonon_energy_weights_J*n,axis=1)-baseline,p,n
    values(bath);lo=bath.copy();hi=np.where(active,np.maximum(2*bath,2.),bath)
    for _ in range(32):
        added,_,_=values(hi);need=(added<target)&active
        if not np.any(need):break
        hi[need]*=2
    else:raise ArithmeticError('No positive-capacity reflected source bracket')
    for _ in range(60):
        mid=(lo+hi)/2;added,_,_=values(mid);below=added<target
        lo=np.where(below&active,mid,lo);hi=np.where((~below)&active,mid,hi)
    T=np.where(active,(lo+hi)/2,bath);added,p,n=values(T);defect=added-target
    if max(abs(float(defect.sum())),float(np.max(abs(defect))))>1e-10*retained:
        raise ArithmeticError('Reflected source violated retained native internal-energy balance')
    receipt.update(wavelength_nm=wavelength_nm,photon_energy_J=photon,retained_fraction=retained/photon,
        retained_energy_J=retained,deposited_energy_eV=retained/elementary_charge,sigma_nm=sigma_nm,
        initialization='Common FD/BE once at fixed original DC gap; virtual image shapes summed BEFORE the energy inversion',
        time_origin='Effective postcascade transfer; no extra cascade delay',
        energy_functional='U_BCS+G-sum(Te*dG/dTe)+phonons at fixed DC gap',
        energy_deposited_J=float(added.sum()),energy_defect_J=float(defect.sum()),
        maximum_nodal_energy_defect_J=float(np.max(abs(defect))),
        relative_total_energy_defect=float(abs(defect.sum())/retained),
        maximum_initial_temperature_K=float(T.max()),minimum_preparation_effective_capacity_ratio=minimum_ratio,
        positive_capacity_checked_throughout_preparation=True,capacity_floor_used=False,
        cut_deposition_J=float(target[cuts].sum()),source_must_not_be_applied_again=True,
        old_prepared_photon_state_reused=False,candidate_physical_admission=False,
        deposition_quadrature='Existing native nodal dual-volume quadrature; exact discrete total-energy normalization')
    return p,n,T,target,receipt
