"""Assess saved DC preparation and preservation data; never run a field solver.

Missing/failed cases remain visible and cannot pass the campaign. Voltage
contrasts use R_load*abs(I_DC), not the relative error of a tiny dark signal.
Generated analysis files may be refreshed; the input data are never changed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np


SCHEMA = 'pysnspd.stage5.prephoton_dc_analysis.v1'
HBAR = 1.054571817e-34
KB = 1.380649e-23


def _read(path):
    def reject_constant(value):
        raise ValueError('Nonfinite JSON constant: '+value)
    with Path(path).open(encoding='utf8') as stream:
        return json.load(stream, parse_constant=reject_constant)


def _write(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)+'\n', encoding='utf8')
    temporary.replace(path)


def _number(value, name, *, positive=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(name+' must be a finite number')
    if positive and value <= 0:
        raise ValueError(name+' must be positive')
    return float(value)


def _arrays(path, names):
    with np.load(path, allow_pickle=False) as data:
        result = {name: np.asarray(data[name]).copy() for name in names}
    if any(not np.all(np.isfinite(value)) for value in result.values()):
        raise ValueError('Nonfinite saved field in '+str(path))
    return result


def _manifest(directory):
    path = directory/'manifest.json'
    if path.exists():
        status = _read(path).get('status')
        if status in ('FAILED', 'RUNNING', 'BLOCKED'):
            raise ValueError('Manifest is '+status)


def _reference(directory):
    summary = _read(directory/'summary.json')
    _manifest(directory)
    arrays = _arrays(directory/'reference.npz', ('delta_bar', 'coordinates_bar', 'area_weights', 'boundary_nodes'))
    n = arrays['delta_bar'].size
    if arrays['delta_bar'].shape != (n,) or arrays['coordinates_bar'].shape != (n, 2):
        raise ValueError('Invalid reference field or coordinate shape')
    mass = arrays['area_weights']
    if mass.shape != (n,) or np.any(mass <= 0):
        raise ValueError('Reference needs positive nodal areas')
    boundary = arrays['boundary_nodes']
    if boundary.dtype.kind not in 'iu' or np.any(boundary < 0) or np.any(boundary >= n):
        raise ValueError('Invalid boundary indices')
    free = np.ones(n, dtype=bool)
    free[boundary] = False
    if not np.any(free):
        raise ValueError('No interior nodes in reference')
    if 'reference_sha256' in summary:
        actual = hashlib.sha256((directory/'reference.npz').read_bytes()).hexdigest()
        if actual != summary['reference_sha256']:
            raise ValueError('Reference hash differs from its saved summary')
    bulk = summary['bulk']
    gap_bulk = _number(bulk['gap_bulk_bar'], 'gap_bulk_bar', positive=True)
    target = _number(bulk['current_target_A'], 'current_target_A')
    if target == 0:
        raise ValueError('This campaign requires nonzero current')
    current = _number(bulk['current_reference_A'], 'current_reference_A')
    normalized = float(np.average(abs(arrays['delta_bar'][free]), weights=mass[free])/gap_bulk)
    parameters = summary['plan']
    ell = math.sqrt(HBAR*_number(parameters['diffusion_m2_s'], 'D', positive=True)/
                    (2*KB*_number(parameters['Tc_K'], 'Tc', positive=True)))
    inductance = {key: _number(summary[key], key, positive=True) for key in (
        'fixed_external_inductance_H', 'resolved_equilibrium_inductance_H', 'total_equilibrium_inductance_H')}
    if not np.isclose(inductance['fixed_external_inductance_H']+inductance['resolved_equilibrium_inductance_H'],
                      inductance['total_equilibrium_inductance_H'], rtol=1e-9, atol=0):
        raise ValueError('Saved exterior/resolved inductances do not add to total')
    added = _number(parameters['added_inductance_H'], 'added_inductance_H', positive=True)
    if inductance['fixed_external_inductance_H'] < added*(1-1e-10):
        raise ValueError('Fixed exterior omits part of the declared added inductor')
    admitted = (summary.get('stationary') is True and summary.get('admitted_for_dark_hold') is True
                and summary.get('status') == 'DC_REFERENCE_ADMITTED')
    record = dict(available=True, admitted=admitted, status=summary.get('status'),
        current_target_A=target, current_reference_A=current,
        mean_gap_over_bulk=normalized, gap_bulk_bar=gap_bulk,
        length_nm=float(np.ptp(arrays['coordinates_bar'][:, 0])*ell*1e9),
        width_nm=float(np.ptp(arrays['coordinates_bar'][:, 1])*ell*1e9),
        node_count=n, bulk=bulk, **inductance)
    return record, dict(summary=summary, arrays=arrays, ell0_m=ell, free=free)


def _hold(directory, reference):
    summary = _read(directory/'summary.json')
    _manifest(directory)
    history = _read(directory/'history.json')
    if not isinstance(history, list) or len(history) < 2:
        raise ValueError('Need initial and final accepted-time observations')
    names = ('time_ps', 'gap_drift_relative', 'current_cut_mean_A', 'V_out_V', 'terminal_voltage_V')
    series = {key: np.asarray([_number(row[key], key) for row in history]) for key in names}
    times = series['time_ps']
    duration = _number(summary['duration_ps'], 'duration_ps', positive=True)
    if times[0] != 0 or np.any(np.diff(times) <= 0) or not np.isclose(times[-1], duration, rtol=1e-10, atol=1e-12):
        raise ValueError('History must cover the whole declared horizon with increasing times')
    arrays = _arrays(directory/'fields.npz', ('gap', 'initial_gap', 'time_ps', 'area_weights', 'coordinates_bar'))
    mass = arrays['area_weights']
    n = mass.size
    if mass.shape != (n,) or np.any(mass <= 0) or arrays['coordinates_bar'].shape != (n, 2):
        raise ValueError('Invalid hold spatial weights or coordinates')
    if arrays['initial_gap'].shape != (n,) or arrays['gap'].shape != (len(arrays['time_ps']), n):
        raise ValueError('Invalid hold snapshot shapes')
    snapshots = arrays['time_ps']
    if len(snapshots) < 2 or snapshots[0] != 0 or np.any(np.diff(snapshots) <= 0) or not np.isclose(snapshots[-1], duration):
        raise ValueError('Snapshots must cover the complete horizon')
    source = reference['arrays']
    if source['coordinates_bar'].shape != arrays['coordinates_bar'].shape or not np.allclose(
            source['coordinates_bar'], arrays['coordinates_bar'], rtol=0, atol=1e-12):
        raise ValueError('Hold geometry differs from declared reference')
    if not np.allclose(source['delta_bar'], arrays['initial_gap'], rtol=1e-11, atol=1e-13):
        raise ValueError('Hold initial gap differs from declared reference')
    current = _number(summary['initial_discrete_current_A'], 'initial_discrete_current_A')
    target = _number(summary['requested_current_A'], 'requested_current_A')
    circuit = summary['circuit']
    rload = _number(circuit['R_load_ohm'], 'R_load_ohm', positive=True)
    if current == 0 or target == 0:
        raise ValueError('Nonzero DC current is needed for dark-voltage normalization')
    scale = rload*abs(current)
    xy = arrays['coordinates_bar']
    width_bar = float(np.ptp(xy[:, 1]))
    if width_bar <= 0:
        raise ValueError('Need finite strip width for the normal-state voltage scale')
    normal_resistance = _number(reference['summary']['plan']['sheet_resistance_ohm'], 'sheet resistance', positive=True)*float(np.ptp(xy[:, 0]))/width_bar
    device_scale = normal_resistance*abs(current)
    gap_bulk = reference['summary']['bulk']['gap_bulk_bar']
    free = reference['free']
    gap_mean = np.average(abs(arrays['gap'][:, free]), weights=mass[free], axis=1)/gap_bulk
    checks = summary.get('checks', summary.get('gates', {}))
    if not isinstance(checks, dict) or not checks or any(type(v) is not bool for v in checks.values()):
        raise ValueError('Hold summary needs explicit boolean gates')
    admitted = (summary.get('status') == 'COMPLETED' and summary.get('dc_hold_passed') is True
                and all(checks.values()))
    maxima = summary.get('maxima', summary.get('maximum', {}))
    if not isinstance(maxima, dict):
        raise ValueError('Saved maxima must be a mapping')
    maximum_sources = {}

    def maximum(key, sampled_values):
        sampled = float(np.max(sampled_values))
        if key not in maxima:
            maximum_sources[key] = 'SAMPLED_HISTORY_ONLY_FALLBACK'
            return sampled
        value = _number(maxima[key], 'maxima.'+key)
        if value < 0 or value < sampled*(1-1e-10):
            raise ValueError('Saved every-step maximum is smaller than a sampled value: '+key)
        maximum_sources[key] = 'SUMMARY_MAXIMUM_OVER_EVERY_ACCEPTED_STEP'
        return value

    gap_maximum = maximum('gap_drift_relative', series['gap_drift_relative'])
    readout_maximum = maximum('dark_readout_absolute_V', abs(series['V_out_V']))
    device_maximum = maximum('device_voltage_absolute_V', abs(series['terminal_voltage_V']))
    record = dict(available=True, admitted=admitted, status=summary.get('status'), gates=checks,
        initial_discrete_current_A=current, current_target_A=target, duration_ps=duration,
        dt_ps=_number(summary['dt_ps'], 'dt_ps', positive=True), dark_scale_V=scale,
        maximum_gap_drift_relative=gap_maximum,
        maximum_absolute_dark_readout_V=readout_maximum,
        maximum_dark_readout_over_Rload_I=readout_maximum/scale,
        maximum_absolute_device_voltage_V=device_maximum,
        device_voltage_scale_V=device_scale,
        maximum_device_voltage_over_Rnormal_I=device_maximum/device_scale,
        device_voltage_scale_definition='R_sheet*(resolved length/width)*abs(initial discrete current); normal-state strip voltage scale, not an experimental trigger',
        initial_mean_gap_over_bulk=float(gap_mean[0]), final_mean_gap_over_bulk=float(gap_mean[-1]),
        final_current_cut_mean_A=float(series['current_cut_mean_A'][-1]),
        energy_balance_relative_to_condensation_scale=summary.get('energy_balance_relative_to_condensation_scale'),
        energy_normalization_definition='U0*sum(area_weights*abs(initial_delta_bar)^2), not F_normal-F_superconducting; legacy summary key calls it condensation_scale_J',
        maximum_sources=maximum_sources, maxima=maxima,
        history_sampling=summary.get('history_sampling', 'Saved observation times only; full-step sampling not asserted'),
        saved_history_observations=len(history), accepted_steps=summary.get('accepted_steps'))
    return record, dict(summary=summary, series=series, arrays=arrays,
                        mean_gap_over_bulk=gap_mean, reference=reference)


def _same_physics(a, b):
    for key in ('T_K', 'Tc_K', 'diffusion_m2_s', 'sheet_resistance_ohm', 'width_m', 'current_A'):
        if not np.isclose(a['summary']['plan'][key], b['summary']['plan'][key], rtol=1e-10, atol=0):
            raise ValueError('Comparison changes physical parameter '+key)


def _contrast(series_a, series_b, field, scale=1.):
    grid = np.union1d(series_a['time_ps'], series_b['time_ps'])
    overlap = grid[(grid >= max(series_a['time_ps'][0], series_b['time_ps'][0])) &
                   (grid <= min(series_a['time_ps'][-1], series_b['time_ps'][-1]))]
    if len(overlap) < 2:
        raise ValueError('No accepted-time overlap for comparison')
    delta = (np.interp(overlap, series_a['time_ps'], series_a[field])-
             np.interp(overlap, series_b['time_ps'], series_b[field]))
    return float(np.max(abs(delta))/scale)


def _compare(kind, identifiers, refs, holds, ref_data, hold_data, hold_sources, limits):
    if kind not in ('mesh', 'length', 'time'):
        raise ValueError('Unknown comparison kind: '+kind)
    if not isinstance(identifiers, list) or len(identifiers) != 2:
        raise ValueError('A comparison must declare two case ids')
    a, b = identifiers
    ra, rb = hold_sources.get(a, a), hold_sources.get(b, b)
    if ra not in ref_data or rb not in ref_data or a not in hold_data or b not in hold_data:
        raise ValueError('Comparison lacks complete reference/hold data')
    _same_physics(ref_data[ra], ref_data[rb])
    ha, hb = holds[a], holds[b]
    if not np.isclose(ha['duration_ps'], hb['duration_ps'], rtol=1e-10, atol=1e-12):
        raise ValueError('Comparison horizons differ; overlap is not full-horizon admission')
    if kind == 'time' and ra != rb:
        raise ValueError('Time-step comparison must use the same spatial reference')
    if kind == 'time' and np.isclose(ha['dt_ps'], hb['dt_ps'], rtol=1e-10, atol=0):
        raise ValueError('Time-step comparison needs distinct time steps')
    same_length = np.isclose(refs[ra]['length_nm'], refs[rb]['length_nm'], rtol=1e-8, atol=0)
    if kind == 'mesh' and not same_length:
        raise ValueError('Mesh comparison changes the physical domain length')
    if kind == 'length' and same_length:
        raise ValueError('Length comparison needs distinct physical lengths')
    ca, cb = hold_data[a]['summary']['circuit'], hold_data[b]['summary']['circuit']
    for key in ('R_load_ohm', 'R_bias_ohm', 'L_bias_H', 'C_couple_F'):
        if not np.isclose(ca[key], cb[key], rtol=1e-10, atol=0):
            raise ValueError('Comparison changes circuit component '+key)
    if kind != 'length' and not np.isclose(ca['Lk_ext_H'], cb['Lk_ext_H'], rtol=1e-8, atol=0):
        raise ValueError('Mesh/time comparison changes fixed exterior inductance')
    current_scale = abs(refs[ra]['current_target_A'])
    voltage_scale = ca['R_load_ohm']*current_scale
    sa, sb = hold_data[a]['series'], hold_data[b]['series']
    ga = dict(time_ps=hold_data[a]['arrays']['time_ps'], gap=hold_data[a]['mean_gap_over_bulk'])
    gb = dict(time_ps=hold_data[b]['arrays']['time_ps'], gap=hold_data[b]['mean_gap_over_bulk'])
    measures = dict(
        reference_mean_gap_over_bulk_difference=abs(refs[ra]['mean_gap_over_bulk']-refs[rb]['mean_gap_over_bulk']),
        reference_current_difference_over_target=abs(refs[ra]['current_reference_A']-refs[rb]['current_reference_A'])/current_scale,
        trajectory_mean_gap_over_bulk_difference=_contrast(ga, gb, 'gap'),
        trajectory_gap_drift_difference=_contrast(sa, sb, 'gap_drift_relative'),
        trajectory_current_difference_over_target=_contrast(sa, sb, 'current_cut_mean_A', current_scale),
        trajectory_dark_voltage_difference_over_Rload_target=_contrast(sa, sb, 'V_out_V', voltage_scale))
    gates = {
        key: value <= (limits['gap_comparison_relative'] if 'gap' in key else
            limits['normalized_dark_readout_difference'] if 'voltage' in key else
            limits['current_comparison_relative']) for key, value in measures.items()}
    cases_pass = all((refs[ra]['admitted'], refs[rb]['admitted'], ha['admitted'], hb['admitted']))
    return dict(kind=kind, cases=identifiers, available=True, admitted=cases_pass and all(gates.values()),
        gates=gates, measures=measures, current_normalization_A=current_scale,
        voltage_normalization_V=voltage_scale, limits=limits,
        interpolation='Linear interpolation on union of saved accepted-time observations within full common horizon; snapshot times for area-weighted gap. Saved observations can be sparse (about 101 times); this is not an every-step trajectory contrast. Individual case gates retain the runner every-accepted-step maxima.',
        reference_gap_definition='Difference of area-weighted interior |Delta|/Delta_bulk, same physical DC current.',
        diagnostics_not_gated=dict(trajectory_absolute_device_voltage_difference_V=_contrast(sa, sb, 'terminal_voltage_V'),
            meaning='Vdev contrast retained explicitly; no threshold invented beyond the declared campaign acceptance.'),
        voltage_definition='Absolute difference of signed V_out divided by R_load*abs(target current); no relative error of near-zero signal.')


def _label(name, record):
    return f"{record['current_target_A']*1e6:g} µA; L={record['length_nm']:g} nm; N={record['node_count']}"


def _plots(output, result, ref_data, hold_data, plan):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False,
                         'savefig.dpi': 180, 'figure.facecolor': 'white'})
    directory = output/'figures'
    directory.mkdir(exist_ok=True)
    generated = []
    if ref_data:
        names = list(ref_data)
        delta0 = float(plan.get('hold_template', {}).get('delta0_over_kBTc', 1.764))
        fig, axes = plt.subplots(len(names), 1, figsize=(8.4, 2.4*len(names)), squeeze=False)
        for ax, name in zip(axes[:, 0], names):
            data = ref_data[name]; arrays = data['arrays']; record = result['references'][name]
            x = arrays['coordinates_bar'][:, 0]*data['ell0_m']*1e9
            x = x-.5*(x.max()+x.min())
            values = abs(arrays['delta_bar'])/delta0
            ax.scatter(x, values, s=3, alpha=.2, color='#2573ad', label='Todos los nodos; sin resta de referencia')
            bins = np.linspace(x.min(), x.max()+1e-10, 25)
            labels = np.digitize(x, bins)-1
            bx, by = [], []
            for index in range(len(bins)-1):
                mask = labels == index
                if np.any(mask):
                    bx.append(np.average(x[mask], weights=arrays['area_weights'][mask]))
                    by.append(np.average(values[mask], weights=arrays['area_weights'][mask]))
            ax.plot(bx, by, color='#ca6b28', label='Promedio por bin x, ponderado por área')
            ax.axhline(record['gap_bulk_bar']/delta0, color='black', ls='--', lw=1, label='Bulk con la misma corriente')
            ax.set(xlabel='x desde el centro [nm]', ylabel=r'$|\Delta|/\Delta_0$', title=_label(name, record))
            ax.ticklabel_format(axis='y', useOffset=False)
            ax.legend(fontsize=7, loc='best')
        fig.suptitle(f'Preparación DC sin fotón. Δ₀ = {delta0:g} kB Tc; fases qx, no contactos metálicos', fontsize=11)
        fig.tight_layout(rect=(0, 0, 1, .98))
        path = directory/'01_gap_profiles.png'; fig.savefig(path); plt.close(fig)
        generated.append(dict(path=path.relative_to(output).as_posix(), quantity='Initial |Delta|/Delta0 at all nodes and area-weighted x bins, versus current-carrying bulk. Delta0 is the declared energy scale.'))

        selected = [name for name in plan.get('comparisons', {}).get('mesh', []) if name in ref_data]
        if not selected:
            selected = names[:2]
        fig, axes = plt.subplots(len(selected), 1, figsize=(8.4, 2.8*len(selected)), squeeze=False)
        all_values = np.concatenate([abs(ref_data[name]['arrays']['delta_bar'])/result['references'][name]['gap_bulk_bar'] for name in selected])
        lo, hi = float(all_values.min()), float(all_values.max())
        if hi-lo < 1e-10:
            lo -= 1e-6; hi += 1e-6
        for ax, name in zip(axes[:, 0], selected):
            data = ref_data[name]; arrays = data['arrays']; record = result['references'][name]
            xy = arrays['coordinates_bar']*data['ell0_m']*1e9
            xy = xy-.5*(xy.max(axis=0)+xy.min(axis=0))
            values = abs(arrays['delta_bar'])/record['gap_bulk_bar']
            image = ax.scatter(xy[:, 0], xy[:, 1], c=values, s=12, vmin=lo, vmax=hi, cmap='viridis', rasterized=True)
            ax.set(aspect='equal', xlabel='x desde el centro [nm]', ylabel='y desde el centro [nm]', title=_label(name, record))
            fig.colorbar(image, ax=ax, label=r'$|\Delta|/\Delta_{\rm bulk}(I_{DC})$', format='%.6f')
        fig.suptitle('Mapas iniciales DC: valores nodales; misma escala de color; referencia no restada', fontsize=11)
        fig.tight_layout(rect=(0, 0, 1, .97))
        path = directory/'02_gap_maps.png'; fig.savefig(path); plt.close(fig)
        generated.append(dict(path=path.relative_to(output).as_posix(), quantity='Initial nodal gap magnitude divided by each current-carrying bulk reference; common color limits, no spatial interpolation.'))

        fig, axes = plt.subplots(1, 2, figsize=(10, 5.6))
        positions = np.arange(len(names))
        added = np.asarray([float(ref_data[name]['summary']['plan']['added_inductance_H'])*1e9 for name in names])
        exterior = np.asarray([result['references'][name]['fixed_external_inductance_H']*1e9 for name in names])
        resolved = np.asarray([result['references'][name]['resolved_equilibrium_inductance_H']*1e9 for name in names])
        omitted = exterior-added
        for ax, include_added in zip(axes, (True, False)):
            base = added if include_added else np.zeros(len(names))
            if include_added:
                ax.bar(positions, added, color='#4b7392', label='Inductor exterior añadido: estimación publicada')
            ax.bar(positions, omitted, bottom=base, color='#e0a456', label='Tramo activo no resuelto: referencia fija')
            ax.bar(positions, resolved, bottom=base+omitted, color='#9bc9ad', label='Tramo resuelto: contribución de referencia')
            ax.set_xticks(positions, [f"{result['references'][name]['current_target_A']*1e6:g} µA\n{result['references'][name]['length_nm']:g} nm\nN={result['references'][name]['node_count']}" for name in names], fontsize=7)
            ax.set(ylabel='Inductancia de referencia [nH]', title='Rama completa contabilizada' if include_added else 'Sólo el tramo activo de 5 µm')
            ax.set_ylim(bottom=0)
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, fontsize=7, loc='upper center', bbox_to_anchor=(.5, .93), ncol=2)
        fig.suptitle('Partición inicial constante: exterior = añadido + activo no resuelto; no Lk(t)', y=.995, fontsize=11)
        fig.tight_layout(rect=(0, 0, 1, .81))
        path = directory/'04_fixed_inductance.png'; fig.savefig(path, bbox_inches='tight'); plt.close(fig)
        generated.append(dict(path=path.relative_to(output).as_posix(), quantity='Fixed initial inductance accounting in nH: added 96 nH estimate, omitted active bulk and resolved equilibrium contribution. Not a dynamic inductance measurement; tapers/connections unquantified.'))
    if hold_data:
        fig, axes = plt.subplots(3, 1, figsize=(9, 10), sharex=True)
        for index, (name, data) in enumerate(hold_data.items()):
            series = data['series']
            reference_name = result['holds'][name]['reference']
            label = _label(reference_name, result['references'][reference_name])+f"; dt={result['holds'][name]['dt_ps']:.3g} ps"
            spacing = max(1, len(series['time_ps'])//9)
            style = dict(label=label, linestyle=('-', '--', ':', '-.')[index % 4],
                         marker=('o', 's', '^', 'x', 'd')[index % 5], markersize=3,
                         markevery=(index % spacing, spacing), linewidth=1.3)
            axes[0].plot(series['time_ps'], 100*series['gap_drift_relative'], **style)
            axes[1].plot(series['time_ps'], 1e6*series['terminal_voltage_V'], **style)
            axes[2].plot(series['time_ps'], 1e9*series['V_out_V'], **style)
        axes[0].set(ylabel='Máx. cambio de |Δ| / media inicial [%]', title='Norma máxima nodal de la diferencia respecto a t=0')
        axes[1].set(ylabel='Vdev = φL − φR [µV]', title='Caída firmada del segmento resuelto; sin filtrado de la lectura')
        axes[2].set(xlabel='Tiempo físico [ps]', ylabel='Vout [nV]', title='Voltaje firmado de carga; línea basal DC ideal = 0 V')
        axes[1].axhline(0, color='black', lw=.8, ls=':')
        axes[2].axhline(0, color='black', lw=.8, ls=':')
        axes[0].legend(fontsize=7, ncol=1)
        fig.suptitle('Conservación DC sin fotón: observaciones guardadas; sin umbral experimental', fontsize=11)
        fig.tight_layout(rect=(0, 0, 1, .97))
        path = directory/'03_dark_hold_history.png'; fig.savefig(path); plt.close(fig)
        generated.append(dict(path=path.relative_to(output).as_posix(), quantity='Saved accepted-time observations, possibly sparse; lines connect stored observations, not every integration step. At each saved time: max nodal gap-magnitude drift from t=0 normalized by initial area mean in percent; signed device voltage phiL-phiR in microvolts; signed load voltage Vout in nanovolts. Ideal DC baseline is zero. Full-step maxima are reported separately from summary.json. No photon or experimental trigger.'))
    return generated


def _markdown(result):
    lines = ['# Preparación DC anterior al fotón', '',
        '**Estado: '+result['status']+'.**', '', result['scope'], '',
        '## Referencias con corriente', '',
        '| Caso | Datos | Admitido | I medida [µA] | Media interior / bulk | L exterior fija [nH] |',
        '|---|---|---|---:|---:|---:|']
    for name, record in result['references'].items():
        if record['available']:
            lines.append(f"| {name} | disponibles | {record['admitted']} | {record['current_reference_A']*1e6:.7g} | {record['mean_gap_over_bulk']:.8g} | {record['fixed_external_inductance_H']*1e9:.7g} |")
        else:
            lines.append(f"| {name} | incompletos: {record['reason']} | False | — | — | — |")
    lines += ['', '## Conservación temporal', '',
        '| Caso | Admitido | Horizonte [ps] | Deriva máxima gap [%] | Máx. abs(Vout) [nV] | Vout / (RL abs(Iref)) |',
        '|---|---|---:|---:|---:|---:|']
    for name, record in result['holds'].items():
        if record['available']:
            lines.append(f"| {name} | {record['admitted']} | {record['duration_ps']:.7g} | {record['maximum_gap_drift_relative']*100:.6g} | {record['maximum_absolute_dark_readout_V']*1e9:.6g} | {record['maximum_dark_readout_over_Rload_I']:.6g} |")
        else:
            lines.append(f"| {name} | False: {record['reason']} | — | — | — | — |")
    lines += ['', 'Los máximos de estas tablas se leen de summary.json, que acumula todos los pasos aceptados. '
        'Cuando esa entrada falta, se usa únicamente el máximo de las observaciones guardadas y se marca '
        'como tal a continuación; no se atribuye cobertura temporal completa a ese reemplazo.', '']
    for name, record in result['holds'].items():
        if record['available']:
            sampled = [key for key, value in record['maximum_sources'].items() if value == 'SAMPLED_HISTORY_ONLY_FALLBACK']
            if sampled:
                lines.append('- '+name+': máximos sólo muestreados para '+', '.join('`'+key+'`' for key in sampled)+'.')
    lines += ['', 'La lectura filtrada puede ser pequeña aunque el segmento mantenga una caída numérica. '
        'Por eso Vdev se informa por separado, sin ocultarlo en Vout.', '',
        '| Caso | Máx. abs(Vdev) [µV] | Escala Rnormal abs(Iref) [mV] | Cociente adimensional |',
        '|---|---:|---:|---:|']
    for name, record in result['holds'].items():
        if record['available']:
            lines.append(f"| {name} | {record['maximum_absolute_device_voltage_V']*1e6:.7g} | {record['device_voltage_scale_V']*1e3:.7g} | {record['maximum_device_voltage_over_Rnormal_I']:.7g} |")
    lines += ['', 'Rnormal = R□ L/W para el segmento resuelto. Este cociente indica la escala de la '
        'caída respecto al estado normal; no es un umbral experimental ni un criterio de pase nuevo. '
        'El balance energético térmico usa U0 Σ área |Δ(0)|² como normalización numérica '
        '(las áreas y Δ son adimensionales); esa escala no es Fnormal − Fsuperconductor.', '']
    lines += ['', '## Comparaciones declaradas', '',
        'Se comparan igual corriente, material y horizonte. Las diferencias de gap se expresan respecto al bulk; '
        'el voltaje usa la escala RL |I|, sin dividir por una señal oscura casi nula. '
        'Las curvas se interpolan sobre los instantes guardados (pueden ser unas 101 observaciones), '
        'no sobre cada paso del integrador; los pases individuales usan además sus máximos de todos los pasos. '
        'Los límites son los registrados en workflow_plan.json.', '']
    for name, record in result['comparisons'].items():
        lines.append('- **'+name+'**: '+('admitida' if record['admitted'] else 'no admitida')+'.')
        if record['available']:
            for key, value in record['measures'].items():
                lines.append(f'  - `{key}` = {value:.7g}; pase = {record["gates"][key]}.')
        else:
            lines.append('  - '+record['reason'])
    lines += ['', '## Figuras físicas', '']
    for item in result['figures']:
        lines += ['!['+item['quantity']+']('+item['path']+')', '', item['quantity'], '']
    if result['errors']:
        lines += ['## Incidencias', '']+['- '+item for item in result['errors']]+['']
    lines += ['## Alcance de la decisión', '',
        'Este resultado evalúa una preparación DC interior y su conservación térmica. '
        'Los bordes continúan el superconductor con corriente, la fuente es DC y la inductancia exterior es fija. '
        'No hay fotón, fuente AC ni distribución de posiciones longitudinales. '
        'El potencial usa el cierre óhmico heredado y las poblaciones permanecen térmicas. '
        'El pase no demuestra el balance de calor no lineal ni la cinética fotónica completa. '
        'Los datos ausentes o fallidos no se reemplazan ni se consideran pases.', '']
    return '\n'.join(lines)


def analyze(output: Path) -> dict:
    """Generate reproducible derived reports; incomplete input never admits DC."""
    output = Path(output)
    if not output.is_dir():
        raise FileNotFoundError('Existing campaign output directory required')
    result = dict(schema=SCHEMA, status='DC_PREPARATION_NOT_ADMITTED', prephoton_dc_admitted=False,
        references={}, holds={}, comparisons={}, figures=[], errors=[],
        scope='Photon-free current-carrying interior strip; fixed exterior inductance and full CM DC circuit. Thermal preservation only, not nonlinear photon heat admission.',
        photon=False, ac_source=False, nonlinear_heat_closed=False, production_changed=False)
    ref_data, hold_data = {}, {}
    try:
        plan = _read(output/'workflow_plan.json')
        for group in ('references', 'holds'):
            names = [row['id'] for row in plan[group]]
            if not names or len(set(names)) != len(names):
                raise ValueError('Need unique nonempty '+group+' ids')
            if any(Path(name).name != name or name in ('.', '..') for name in names):
                raise ValueError('Case ids must be simple directory names')
        reference_ids = {row['id'] for row in plan['references']}
        if any(row['reference'] not in reference_ids for row in plan['holds']):
            raise ValueError('Every hold must name a declared reference')
        limits = {key: _number(plan['acceptance'][key], key, positive=True) for key in (
            'gap_comparison_relative', 'current_comparison_relative', 'normalized_dark_readout_difference')}
        result['acceptance'] = limits
    except Exception as exc:
        result['errors'].append('Invalid workflow plan: '+str(exc))
        _write(output/'analysis.json', result)
        (output/'analysis.md').write_text(_markdown(result), encoding='utf8')
        return result
    for item in plan['references']:
        name = item['id']
        try:
            record, data = _reference(output/'references/cases'/name)
            result['references'][name] = record; ref_data[name] = data
        except Exception as exc:
            result['references'][name] = dict(available=False, admitted=False, reason=str(exc))
    hold_sources = {item['id']: item['reference'] for item in plan['holds']}
    for item in plan['holds']:
        name = item['id']
        try:
            if item['reference'] not in ref_data:
                raise ValueError('Declared reference data unavailable: '+item['reference'])
            record, data = _hold(output/'holds/cases'/name, ref_data[item['reference']])
            record['admitted'] = record['admitted'] and result['references'][item['reference']]['admitted']
            record['reference'] = item['reference']
            result['holds'][name] = record; hold_data[name] = data
        except Exception as exc:
            result['holds'][name] = dict(available=False, admitted=False, reason=str(exc))
    comparisons = plan.get('comparisons', {})
    missing_comparisons = {'mesh', 'length', 'time'}-set(comparisons)
    if missing_comparisons:
        result['errors'].append('Required comparisons not declared: '+', '.join(sorted(missing_comparisons)))
    for kind, identifiers in comparisons.items():
        try:
            result['comparisons'][kind] = _compare(kind, identifiers, result['references'], result['holds'],
                ref_data, hold_data, hold_sources, limits)
        except Exception as exc:
            result['comparisons'][kind] = dict(available=False, admitted=False, cases=identifiers, reason=str(exc))
    admitted = (not result['errors'] and all(record['admitted'] for group in
        ('references', 'holds', 'comparisons') for record in result[group].values()))
    result['prephoton_dc_admitted'] = bool(admitted)
    result['status'] = 'DC_PREPARATION_ADMITTED' if admitted else 'DC_PREPARATION_NOT_ADMITTED'
    try:
        result['figures'] = _plots(output, result, ref_data, hold_data, plan)
    except Exception as exc:
        result['errors'].append('Figure generation failed; numerical verdict retained: '+str(exc))
    _write(output/'analysis.json', result)
    (output/'analysis.md').write_text(_markdown(result), encoding='utf8')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-root', type=Path, required=True)
    args = parser.parse_args()
    result = analyze(args.output_root)
    print(json.dumps(result, ensure_ascii=False, allow_nan=False), flush=True)
    return 0 if result['prephoton_dc_admitted'] else 2


if __name__ == '__main__':
    raise SystemExit(main())
