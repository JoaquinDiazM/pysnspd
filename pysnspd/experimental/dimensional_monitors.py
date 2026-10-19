"""Accepted-state guards for finite domain, electrical scope and long recovery.

Bounds are declared engineering diagnostics, not measured device ratings.
A bipolar readout pulse is allowed; only global branch reversal is excluded.
Latching is a possible physical outcome, never mislabeled a solver defect.
"""
from dataclasses import dataclass
import numpy as np
from .dimensional_transition import population_distance


@dataclass(frozen=True)
class MonitorConfig:
    collar_nm: float = 25.
    collar_temperature_rise_K: float = .5
    collar_gap_fraction: float = .9
    collar_confirmation_ps: float = 2.
    current_limit_uA: float = 28.453333333333333
    device_voltage_limit_mV: float = 302.31666666666666
    readout_voltage_limit_mV: float = 1.4226666666666667
    capacitor_voltage_limit_mV: float = 430.
    latching_earliest_ns: float = 5.
    latching_window_ns: float = 1.
    latching_plateau_fraction: float = .01
    latching_normal_length_nm: float = 20.
    recovery_earliest_ns: float = 2.
    recovery_confirmation_ns: float = 1.
    recovery_temperature_rise_K: float = .02
    recovery_current_fraction: float = .01
    recovery_gap_fraction: float = .01
    recovery_native_energy_fraction: float = .01
    recovery_voltage_uV: float = 1.


def longest_hot_length_nm(mask,volumes,width,thickness):
    length=volumes/(width*thickness)*1e9
    largest=current=0.
    for hot,dx in zip(mask,length):
        current=current+dx if hot else 0.;largest=max(largest,current)
    return largest


def metrics(system,geometry,state,rhs,dc,config,reference_J,extended=False):
    f,k=system.functional,system.kinetics
    Ib,Is,vc=state.circuit_SI;IDC=system.circuit_scales[1]
    vout=system.circuit.R_load_ohm*(Ib-Is)
    out=dict(Ib_uA=float(Ib*1e6),Is_uA=float(Is*1e6),vc_mV=float(vc*1e3),
        Vout_mV=float(vout*1e3),Vdev_mV=float(rhs.electrical.Vdev_V*1e3),
        current_error_fraction=float(max(abs(Ib-IDC),abs(Is-IDC))/IDC),
        max_temperature_K=float(np.max(rhs.temperature_K)),minimum_temperature_K=float(np.min(rhs.temperature_K)),
        gap_amplitude_error=float(np.max(abs(abs(state.delta_bar)/abs(dc.delta_bar)-1))),
        hot_length_nm=longest_hot_length_nm(abs(state.delta_bar)<.2*abs(dc.delta_bar),
            f.volumes_m3,geometry.width_m,f.thickness_m) if extended else 0.,
        film_power_residual_W=float(rhs.film_power_residual_W),
        boundary_alarm=False,collar_max_Te_K=float(k.bath_temperature_K),
        collar_min_gap_ratio=1.,collar_phonon_hot=False,recovered=False)
    if extended:
        x=geometry.coordinates_m[:,0];dist=np.minimum(x-x.min(),x.max()-x)
        collar=(dist>1e-16)&(dist<=config.collar_nm*1e-9+1e-16)
        if not np.any(collar):raise ValueError("Boundary collar contains no free nodes")
        out["collar_max_Te_K"]=float(np.max(rhs.temperature_K[collar]))
        out["collar_min_gap_ratio"]=float(np.min(abs(state.delta_bar[collar])/abs(dc.delta_bar[collar])))
        actual=np.sum(k.phonon_energy_weights_J[collar]*state.n[collar],axis=1)
        limit=np.sum(k.phonon_energy_weights_J[collar]*k.bose(
            k.bath_temperature_K+config.collar_temperature_rise_K),axis=1)
        out["collar_phonon_hot"]=bool(np.any(actual>limit))
        out["boundary_alarm"]=bool(out["collar_max_Te_K"]>k.bath_temperature_K+config.collar_temperature_rise_K or
            out["collar_min_gap_ratio"]<config.collar_gap_fraction or out["collar_phonon_hot"])
    small_voltage=max(abs(out["Vout_mV"]),abs(out["Vdev_mV"]),abs(out["vc_mV"]))<=config.recovery_voltage_uV*1e-3
    if (small_voltage and out["current_error_fraction"]<=config.recovery_current_fraction and
        out["max_temperature_K"]<=k.bath_temperature_K+config.recovery_temperature_rise_K and
        out["minimum_temperature_K"]>=max(0.,k.bath_temperature_K-config.recovery_temperature_rise_K) and
        out["gap_amplitude_error"]<=config.recovery_gap_fraction):
        ph=np.sum(k.phonon_energy_weights_J*state.n,axis=1)
        low=np.sum(k.phonon_energy_weights_J*k.bose(max(0.,k.bath_temperature_K-config.recovery_temperature_rise_K)),axis=1)
        high=np.sum(k.phonon_energy_weights_J*k.bose(k.bath_temperature_K+config.recovery_temperature_rise_K),axis=1)
        out["phonon_recovery_moment_ok"]=bool(np.all((ph>=low)&(ph<=high)))
        out["native_energy_fraction"]=population_distance(system,state,dc,reference_J)
        out["recovered"]=(out["phonon_recovery_moment_ok"] and
                          out["native_energy_fraction"]<=config.recovery_native_energy_fraction)
    return out


def observe(t_ps,phase,values,config,tracker,bias_uA,sheet_resistance_ohm,width_nm):
    """Update JSON-checkpointable monitors only on accepted official states."""
    if any(not np.isfinite(values[k]) for k in ("Ib_uA","Is_uA","vc_mV","Vout_mV","Vdev_mV")):
        return "CIRCUIT_NONFINITE_GUARD","Nonfinite accepted electrical observable"
    if min(values["Ib_uA"],values["Is_uA"])<=bias_uA*1e-6:
        return "CIRCUIT_POLARITY_GUARD","Global branch current reached its positive polarity margin"
    if max(values["Ib_uA"],values["Is_uA"])>config.current_limit_uA:
        return "CIRCUIT_OVERDRIVE_GUARD","Global current exceeds the declared uA ceiling"
    for field,limit in (("Vdev_mV",config.device_voltage_limit_mV),
                        ("Vout_mV",config.readout_voltage_limit_mV),
                        ("vc_mV",config.capacitor_voltage_limit_mV)):
        if abs(values[field])>limit:
            return "CIRCUIT_OVERVOLTAGE_GUARD",field+" exceeds "+str(limit)+" mV"
    if phase!="1d":return None,None
    if values["boundary_alarm"]:
        if tracker.get("boundary_since_ps") is None:tracker["boundary_since_ps"]=t_ps
        if t_ps-tracker["boundary_since_ps"]>=config.collar_confirmation_ps:
            return "POSSIBLE_DOMAIN_ERROR","Thermal/suppression front persisted in the free outer collar"
    else:tracker["boundary_since_ps"]=None
    if t_ps>=config.recovery_earliest_ns*1e3 and values["recovered"]:
        if tracker.get("recovery_since_ps") is None:tracker["recovery_since_ps"]=t_ps
        if t_ps-tracker["recovery_since_ps"]>=config.recovery_confirmation_ns*1e3:
            return "RECOVERY_CONFIRMED","Circuit, gap, temperatures and native populations returned within declared margins"
    else:tracker["recovery_since_ps"]=None
    if t_ps<tracker.get("next_plateau_ps",0.):return None,None
    tracker["next_plateau_ps"]=t_ps+10.
    sample=[t_ps]+[values[k] for k in ("Ib_uA","Is_uA","Vdev_mV","max_temperature_K","hot_length_nm")]
    rows=tracker.setdefault("plateau_samples",[]);rows.append(sample)
    window=config.latching_window_ns*1e3
    while len(rows)>2 and rows[1][0]<t_ps-window:rows.pop(0)
    if t_ps<config.latching_earliest_ns*1e3 or t_ps-rows[0][0]<window:return None,None
    a=np.asarray(rows,float)
    scales=np.array([bias_uA,bias_uA,max(1.,abs(values["Vdev_mV"])),
                     max(8.65,values["max_temperature_K"])])
    plateau=bool(np.all(np.ptp(a[:,1:5],axis=0)<=config.latching_plateau_fraction*scales))
    minimum_resistance=sheet_resistance_ohm*config.latching_normal_length_nm/width_nm
    # V_mV/I_uA * 1000 gives ohm. AC readout can have Vout=0 in a latched DC state.
    resistance=values["Vdev_mV"]/values["Is_uA"]*1e3
    if (plateau and values["hot_length_nm"]>=config.latching_normal_length_nm and
        values["Is_uA"]>.05*bias_uA and resistance>=minimum_resistance):
        return "POSSIBLE_LATCHING","Persistent resistive normal segment and electrical/thermal plateau; physical diagnosis pending"
    return None,None
