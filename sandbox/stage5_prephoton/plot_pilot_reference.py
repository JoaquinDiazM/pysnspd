"""Plot the static intrinsic-DC pilot without duplicating its spectral archive.

The first panel separates current depairing from the zero-current equilibrium
gap. The second resolves the small discretization residual around the biased
bulk. All nodes are plotted; longitudinal means use original dual-cell areas.
This artifact contains no temporal or photon admission.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import Boltzmann, elementary_charge

from pysnspd.experimental.bulk_current_reference import solve_bulk_reference, physical_scales


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def area_binned_profile(x, values, weights, bins=20):
    limits = np.linspace(x.min(), x.max(), bins+1)
    ids = np.minimum(np.searchsorted(limits, x, side='right')-1, bins-1)
    means, centers = [], []
    for index in range(bins):
        active = ids == index
        if np.any(active):
            means.append(float(np.average(values[active], weights=weights[active])))
            centers.append(float(np.average(x[active], weights=weights[active])))
    return np.asarray(centers), np.asarray(means)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--summary', type=Path, required=True)
    parser.add_argument('--output', type=Path,
        default=ROOT/'docs/implementation/stage5/prephoton_dc_20260924/pilot')
    parser.add_argument('--remote-reference',
        default='/home/jdiaz/scratch/prephoton_dc_pilot_20260924/reference.npz')
    parser.add_argument('--runtime-seconds', type=float, default=107.965)
    parser.add_argument('--workers', type=int, default=12)
    args = parser.parse_args()
    summary = json.loads(args.summary.read_text(encoding='utf8'))
    digest = sha(args.reference)
    if digest != summary['reference_sha256']:
        raise ValueError('Reference data hash does not match its summary')
    if summary.get('physical_time_steps') != 0:
        raise ValueError('This static-pilot renderer expects no physical time steps')
    plan = summary['plan']
    with np.load(args.reference, allow_pickle=False) as data:
        gap = abs(data['delta_bar'])
        mass = data['area_weights'].copy()
        coordinates = data['coordinates_bar'].copy()
        boundary = data['boundary_nodes'].copy()
        frequencies = data['epsilon_bar'].copy()
    count = len(frequencies)
    zero = solve_bulk_reference(plan['T_K']/plan['Tc_K'], count, 0.)
    scales = physical_scales(**{key: plan[key] for key in (
        'Tc_K', 'diffusion_m2_s', 'sheet_resistance_ohm')})
    bulk = summary['bulk']['gap_bulk_bar']
    x = (coordinates[:, 0]-.5*(coordinates[:, 0].min()+coordinates[:, 0].max()))*scales['ell0_m']*1e9
    normalized = gap/zero.gap_bar
    residual_percent = 100*(gap/bulk-1)
    centers, mean_gap = area_binned_profile(x, normalized, mass)
    _, mean_residual = area_binned_profile(x, residual_percent, mass)
    bulk_ratio = bulk/zero.gap_bar
    nominal_uA = plan['current_A']*1e6
    measured_uA = summary['bulk']['current_reference_A']*1e6
    plt.rcParams.update({'font.size':10.5, 'axes.titlesize':12, 'axes.labelsize':11,
        'axes.spines.top':False, 'axes.spines.right':False, 'figure.facecolor':'white',
        'savefig.facecolor':'white', 'font.family':'DejaVu Sans'})
    fig, axes = plt.subplots(2, 1, figsize=(10.4, 8.0), sharex=True,
        gridspec_kw={'height_ratios':[1, 1.1]}, constrained_layout=False)
    fig.subplots_adjust(left=.105, right=.98, bottom=.25, top=.81, hspace=.47)
    fig.suptitle('Una ventana interior: depresión por corriente y piso espacial uniforme',
        x=.105, ha='left', fontsize=15, fontweight='bold', y=.975)
    fig.text(.105, .915,
        f'Equilibrio estático sin fotón · T = {plan["T_K"]:g} K · ancho = 80 nm · longitud = {np.ptp(x):.0f} nm\n'
        f'I nominal = {nominal_uA:.2f} µA; corriente discreta = {measured_uA:.4f} µA · {len(gap)} nodos · {count} Matsubara',
        fontsize=10.5, va='center')
    ax = axes[0]
    ax.axhline(1., color='#727985', lw=1.4, ls=':', label='Equilibrio sin corriente, mismo T y corte')
    ax.axhline(bulk_ratio, color='#2a735a', lw=1.6, ls='--', label='Bulk uniforme con I nominal')
    ax.scatter(x, normalized, s=8, color='#3477ab', alpha=.30, linewidths=0,
        label='Todos los nodos de la malla')
    ax.plot(centers, mean_gap, color='#c15b17', lw=1.9, label='Media de 20 franjas, ponderada por área')
    ax.scatter(x[boundary], normalized[boundary], marker='s', s=16, facecolors='none',
        edgecolors='#175b49', linewidths=.9, label='Bordes: mismo bulk con corriente')
    ax.set_ylabel(r'$|\Delta|/\Delta_0(T)$')
    ax.set_title(f'A. La corriente reduce el piso un {100*(1-bulk_ratio):.2f} % respecto de I = 0', loc='left', pad=10)
    ax.set_ylim(min(.952, normalized.min()-.002), 1.003)
    ax.grid(axis='y', alpha=.16)
    ax.legend(loc='upper right', fontsize=8.6, frameon=True, framealpha=.96,
        bbox_to_anchor=(.985, .91), ncol=2)
    ax = axes[1]
    ax.axhline(0, color='#2a735a', lw=1.5, ls='--')
    ax.scatter(x, residual_percent, s=9, color='#3477ab', alpha=.35, linewidths=0)
    ax.plot(centers, mean_residual, color='#c15b17', lw=1.9)
    ax.scatter(x[boundary], residual_percent[boundary], marker='s', s=18,
        facecolors='none', edgecolors='#175b49', linewidths=.9)
    ax.set_title('B. Ampliación: dispersión de la malla alrededor del bulk con corriente', loc='left', pad=10)
    ax.set_ylabel(r'$100\,(|\Delta|/\Delta_{\rm bulk}-1)$ [%]')
    ax.set_xlabel('Posición longitudinal x respecto del centro de la ventana [nm]')
    ax.grid(alpha=.16)
    ax.set_xlim(x.min()-2, x.max()+2)
    extent = max(abs(residual_percent.min()), abs(residual_percent.max()))
    ax.set_ylim(-1.2*extent, 1.2*extent)
    bulk_record = summary['bulk']
    fig.text(.105, .102,
        f'Desviación RMS respecto del bulk = {100*bulk_record["gap_relative_rms_vs_bulk"]:.4f} %; '
        f'rango pico a pico = {100*bulk_record["gap_peak_to_peak_relative"]:.4f} % (nodos interiores).\n'
        f'Variación máxima entre cortes de corriente = {100*bulk_record["crosscut_max_relative_error"]:.4f} % de I nominal.\n'
        r'$\Delta_0(T)$: solución homogénea sin corriente a 0,9 K, suma finita de 256 Matsubara; no es $\Delta(T=0)$.',
        fontsize=9.6, va='center', linespacing=1.55)
    fig.text(.105, .028,
        'Sólo referencia estacionaria. No certifica aún conservación temporal, respuesta al fotón ni calor no térmico.',
        fontsize=9.2, color='#60646d')
    args.output.mkdir(parents=True, exist_ok=True)
    figure_path = args.output/'perfil_dc_intrinseco.png'
    fig.savefig(figure_path, dpi=180)
    plt.close(fig)
    receipt = dict(schema='pysnspd.stage5.prephoton_static_pilot_receipt.v1',
        status=summary['status'], admitted_for_dark_hold=summary['admitted_for_dark_hold'],
        scope='Static current-carrying intrinsic-bulk reference only; physical_time_steps=0; long DC hold remains independent evidence',
        physical_time_steps=0, photon=False, ac_source=False,
        runtime_seconds=args.runtime_seconds, numerical_workers=args.workers,
        runtime_provenance='Observed Geminga pilot wall time and admitted worker allocation, supplied with the static pilot receipt',
        nodes=len(gap), positive_matsubara_frequencies=count,
        optimization_sweeps=len(summary['history']),
        source=dict(remote_reference=args.remote_reference, reference_sha256=digest,
            raw_archive_copied_into_repository=False, summary_sha256=sha(args.summary),
            mesh_path=plan['mesh'], mesh_sha256=plan['mesh_sha256']),
        physical_case={key:plan[key] for key in ('T_K','Tc_K','diffusion_m2_s',
            'sheet_resistance_ohm','width_m','current_A','active_length_m','added_inductance_H')},
        resolved_length_m=summary['resolved_length_m'], q_bar=summary['q_bar'],
        bulk=bulk_record, stationary_gap_mass_rms_relative=summary['metrics']['gap_mass_rms_relative'],
        maximum_spectral_residual=summary['metrics']['maximum_spectral_residual'],
        gap_zero_current_same_temperature_kBTc=zero.gap_bar,
        gap_zero_current_same_temperature_meV=zero.gap_bar*Boltzmann*plan['Tc_K']/elementary_charge*1e3,
        bulk_gap_over_zero_current_gap=bulk_ratio,
        fixed_external_inductance_H=summary['fixed_external_inductance_H'],
        resolved_equilibrium_inductance_H=summary['resolved_equilibrium_inductance_H'],
        total_equilibrium_inductance_H=summary['total_equilibrium_inductance_H'],
        figure=dict(path=figure_path.name, sha256=sha(figure_path),
            upper='All nodal gap magnitudes divided by zero-current equilibrium gap at the same T and finite Matsubara count; no subtraction',
            lower='100*(nodal gap magnitude/current-carrying continuum bulk gap-1), in percent',
            mean='20 non-overlapping equal-width x bins; within each bin the mean x and gap use original dual-cell area weights, including end cells',
            coordinates='x measured from the center of the resolved window, in nm',
            boundary='End-node gap magnitudes equal the same current-carrying continuum bulk'),
        renderer=dict(path='sandbox/stage5_prephoton/plot_pilot_reference.py',sha256=sha(__file__)),
        reproduction_command='python sandbox/stage5_prephoton/plot_pilot_reference.py '
            '--reference '+args.remote_reference+' --summary '+str(Path(args.remote_reference).parent/'summary.json')+
            ' --output docs/implementation/stage5/prephoton_dc_20260924/pilot')
    (args.output/'receipt.json').write_text(json.dumps(receipt, indent=2, ensure_ascii=False, allow_nan=False)+'\n', encoding='utf8')
    print(json.dumps(dict(figure=str(figure_path), receipt=str(args.output/'receipt.json'),
        zero_current_gap=zero.gap_bar, bulk_gap_ratio=bulk_ratio), indent=2))


if __name__ == '__main__':
    main()
