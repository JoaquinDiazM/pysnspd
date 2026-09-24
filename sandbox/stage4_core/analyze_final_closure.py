"""Reconstruct physical observables from the saved final stage-4 campaign.

Only postprocessing is performed. No simulation, fitted correction, tolerance
change or automatic admission is hidden in this analysis. Complex responses
use peak phasors exp(-i omega t), normalized by v_dev=e V_dev/(kB Tc).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.constants import k, e
from scipy.linalg import eigvalsh

ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT/'docs/implementation/stage4/closure_20260924'


def complex_value(value):
    return np.asarray(value['real'])+1j*np.asarray(value['imag'])


def encoded(value):
    return {'real':float(np.real(value)), 'imag':float(np.imag(value))}


def load_arrays(path):
    with np.load(path) as data:
        return {name:data[name].copy() for name in data.files}


def analyze(data):
    out = dict(schema='pysnspd.stage4.final_physics_review.v1',
        scope='Postprocessing of the weak harmonic campaign; no nonlinear energy or photon acceptance',
        definitions=dict(
            response='Peak complex phasors, exp(-i omega t), per v_dev=e V_dev/(kB Tc)',
            field_rms='Area-weighted RMS over non-contact nodes',
            field_comparison='Area-weighted L2 over all nodes; ordinary Euclidean L2 for edge-current coefficients; denominator is refined field norm',
            scalar_comparison='Absolute complex or real change divided by the magnitude of the refined observable',
            crosssection='Cuts halfway between every consecutive distinct nodal x coordinate; sum oriented edge currents crossing each cut',
            crosssection_max='Maximum complex departure from the summed left contact current divided by its magnitude',
            crosssection_rms='Interval-width-weighted RMS complex departure from the interval-width-weighted mean over all cuts, divided by magnitude of that mean',
            radial_heat='Omega^2 sum_i m_i gamma_i |radial_i|^2/(16 R_sheet), W per V_dev_peak^2; 1/2 from peak-to-cycle-average',
            terminal_power='(Re Y_left + Re Y_right)/4 for antisymmetric terminal voltages +/- V_dev/2'),
        references={}, responses={}, comparison={})
    for path in sorted((data/'references/cases').iterdir()):
        metrics = json.loads((path/'summary.json').read_text(encoding='utf8'))['metrics']
        fields = load_arrays(path/'reference.npz')
        gap, xy, mass = abs(fields['delta_bar']), fields['coordinates_bar'], fields['area_weights']
        center = (xy[:,0].min()+xy[:,0].max())/2
        out['references'][path.name] = {key:metrics[key] for key in (
            'gap_mass_rms_relative','minimum_gap_relative','phase_advance_rad',
            'phase_branch_integer','maximum_phase_jump_on_path_rad','phase_path_minimum_gap_relative',
            'maximum_free_current_divergence','maximum_absolute_link_current_bar','maximum_noether_residual','sweep')}
        out['references'][path.name].update(gap_bar_min=float(min(gap)),gap_bar_max=float(max(gap)),
            gap_mass_mean=float(np.sum(mass*gap)/mass.sum()),
            gap_central_mean=float(gap[np.abs(xy[:,0]-center)<1].mean()),
            central_region_definition='Arithmetic node mean in |x-x_center|/ell0 < 1',
            node_count=len(gap),coordinate_range_bar=np.ptp(xy,axis=0).tolist())
    for path in sorted((data/'responses/cases').iterdir()):
        result = json.loads((path/'results.json').read_text(encoding='utf8'))
        plan = json.loads((path/'manifest.json').read_text(encoding='utf8'))['plan']
        tc, rsheet = plan['Tc_K'],plan['sheet_resistance_ohm']
        fields = load_arrays(path/'coupled_fields.npz')
        reference_id = 'phase2' if path.name.startswith('phase2') else 'phase8'
        reference = load_arrays(data/'references/cases'/reference_id/'reference.npz')
        mass, gap = fields['area_weights'],fields['gap_reference']
        xy, edges = fields['coordinates_bar'],fields['edges']
        free = np.ones(len(gap),bool);free[reference['boundary_nodes']] = False

        def rms(value):
            return float(np.sqrt(np.sum(mass[free]*abs(value[free])**2)/mass[free].sum()))

        current = fields['current_response']
        tail,head = edges.T
        distinct_x = np.unique(xy[:,0])
        interval_width = np.diff(distinct_x)
        cuts_x = distinct_x[:-1]+interval_width/2
        cuts = np.array([np.dot((xy[tail,0]<=x).astype(float)-(xy[head,0]<=x).astype(float),current) for x in cuts_x])
        cut_mean = np.sum(interval_width*cuts)/interval_width.sum()
        current_unit = k*tc/(2*e*rsheet)
        left_boundary = np.zeros(len(gap),bool);right_boundary = np.zeros(len(gap),bool)
        boundary = reference['boundary_nodes'];midpoint = .5*(xy[:,0].min()+xy[:,0].max())
        left_boundary[boundary] = xy[boundary,0]<midpoint
        right_boundary[boundary] = xy[boundary,0]>midpoint
        current_left = np.dot(left_boundary[tail].astype(float)-left_boundary[head].astype(float),current)
        current_right = -np.dot(right_boundary[tail].astype(float)-right_boundary[head].astype(float),current)
        yright,yleft = current_right/(2*rsheet),current_left/(2*rsheet)
        impedance = complex_value(result['impedance_ohm'])
        readout = complex_value(result['readout_peak_V'])
        circuit_state = complex_value(result['circuit_state_peak'])
        device_voltage = impedance*circuit_state[1]
        gamma_resolved = fields['gamma_resolved']
        gamma_kwt = gamma_resolved+fields['gamma_residual']
        relative_mobility = eigvalsh(gamma_resolved,gamma_kwt)
        phase_response = fields['angular_response']/abs(gap)
        angular_frequency = result['angular_frequency_per_ps']
        actual_factor = device_voltage*e/(k*tc)
        row = {key:result[key] for key in (
            'omega','mode_count','reference_dc_current_A','resolved_inductance_H','external_inductance_H',
            'admittance_S','impedance_ohm','readout_peak_V','full_field_residuals',
            'resolved_mobility_norm_over_kwt','minimum_residual_mobility_eigenvalue',
            'film_port_power_per_voltage_squared_W','candidate_radial_heat_per_voltage_squared_W',
            'minimum_static_projected_stiffness')}
        row.update(right_terminal_admittance_S=encoded(yright),
            terminal_current_mismatch_relative=float(abs(current_left-current_right)/abs(current_left)),
            power_both_terminals_per_V2=float((yleft.real+yright.real)/4),
            frequency_GHz=angular_frequency*1e3/(2*np.pi),period_ps=2*np.pi/angular_frequency,
            reactive_inductance_nH=-impedance.imag/(angular_frequency*1e12)*1e9,
            device_voltage_peak_V=encoded(device_voltage),device_voltage_magnitude_V=abs(device_voltage),
            readout_magnitude_V=abs(readout),radial_rms_per_vbar=rms(fields['radial_response']),
            radial_max_per_vbar=float(np.max(abs(fields['radial_response']))),
            phase_rms_per_vbar=rms(phase_response),potential_rms_per_vbar=rms(fields['potential_response']),
            relative_radial_response_rms_per_vbar=rms(fields['radial_response']/abs(gap)),
            source_drive_peak_V=plan['bias_voltage_peak_V'],
            radial_rms_at_declared_source_drive=rms(fields['radial_response']*actual_factor),
            crosssection_current_max_deviation_from_left_relative=float(np.max(abs(cuts-current_left))/abs(current_left)),
            crosssection_cut_count=len(cuts_x),
            crosssection_current_rms_deviation_from_mean_relative=float(np.sqrt(np.sum(interval_width*abs(cuts-cut_mean)**2)/interval_width.sum())/abs(cut_mean)),
            current_crosssection_min_abs_uA_per_vbar=float(abs(cuts).min()*current_unit*1e6),
            current_crosssection_max_abs_uA_per_vbar=float(abs(cuts).max()*current_unit*1e6),
            resolved_relative_mobility_eigenvalue_min=float(relative_mobility.min()),
            resolved_relative_mobility_eigenvalue_max=float(relative_mobility.max()),
            radial_heat_over_port_power=result['candidate_radial_heat_per_voltage_squared_W']/result['film_port_power_per_voltage_squared_W'])
        out['responses'][path.name] = row

    base = data/'responses/cases/phase8_w020';refined = data/'responses/cases/phase8_refined'
    base_result = json.loads((base/'results.json').read_text(encoding='utf8'))
    refined_result = json.loads((refined/'results.json').read_text(encoding='utf8'))
    base_fields = load_arrays(base/'coupled_fields.npz');refined_fields = load_arrays(refined/'coupled_fields.npz')
    mass = refined_fields['area_weights']
    for name in ('admittance_S','impedance_ohm','readout_peak_V'):
        out['comparison'][name+'_complex_relative_change'] = abs(complex_value(base_result[name])-complex_value(refined_result[name]))/abs(complex_value(refined_result[name]))
    for name in ('resolved_inductance_H','candidate_radial_heat_per_voltage_squared_W','film_port_power_per_voltage_squared_W'):
        out['comparison'][name+'_relative_change'] = abs(base_result[name]-refined_result[name])/abs(refined_result[name])
    for name in ('radial_response','angular_response','potential_response','current_response'):
        weights = mass if base_fields[name].shape==mass.shape else np.ones(len(base_fields[name]))
        out['comparison'][name+'_weighted_L2_relative_change'] = float(np.sqrt(np.sum(weights*abs(base_fields[name]-refined_fields[name])**2)/np.sum(weights*abs(refined_fields[name])**2)))
        out['comparison'][name+'_max_error_over_refined_max'] = float(np.max(abs(base_fields[name]-refined_fields[name]))/np.max(abs(refined_fields[name])))
    delta_real_y = abs(complex_value(base_result['admittance_S']).real-complex_value(refined_result['admittance_S']).real)
    out['comparison']['ReY_relative_change'] = delta_real_y/abs(complex_value(refined_result['admittance_S']).real)
    out['comparison']['ReY_change_over_total_refined_Y'] = delta_real_y/abs(complex_value(refined_result['admittance_S']))
    out['comparison']['gamma_resolved_over_kwt_max_all_cases'] = max(row['resolved_mobility_norm_over_kwt'] for row in out['responses'].values())
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root',type=Path,default=DEFAULT/'data')
    parser.add_argument('--output',type=Path,default=DEFAULT/'physics_review.json')
    args = parser.parse_args()
    result = analyze(args.data_root)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf8')
    print(json.dumps({'output':str(args.output),'references':len(result['references']),
        'responses':len(result['responses']),'comparison':result['comparison']},indent=2))


if __name__ == '__main__':
    main()
