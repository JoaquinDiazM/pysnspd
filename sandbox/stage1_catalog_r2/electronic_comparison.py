"""Compare the legacy and experimental electronic spectra without a PRE run.

The Simon input is a PHONON DOS and deliberately never appears on these axes.
Identical physical fields isolate solver differences from regulator changes.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pysnspd.usadel.solver import usadel_dos
from pysnspd.experimental.energy_catalog import retarded_spectrum


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--eta', type=float, default=1e-8)
    args = parser.parse_args()
    started = time.perf_counter()
    out = ROOT/'docs/implementation/stage1_r2'
    figures = out/'figures'
    figures.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.size': 11, 'axes.spines.top': False,
                         'axes.spines.right': False, 'savefig.dpi': 180})
    blue, orange, teal = '#2463a6', '#d05f17', '#008b86'
    delta, old_eta = .72, .001
    energy = np.unique(np.r_[np.linspace(0, 1.8, 1001),
                             np.linspace(.68, .76, 801)])
    fig, axes = plt.subplots(2, 2, figsize=(10.2, 5.8),
                             gridspec_kw={'height_ratios': [2.4, 1]},
                             constrained_layout=True)
    rows = []
    for col, gamma in enumerate([0., .2]):
        old = usadel_dos(energy, delta_J=delta, gamma_J=gamma, eta_J=old_eta)
        new = retarded_spectrum(energy, delta=delta, gamma=gamma, eta=old_eta)[0].real
        diff = new-old
        axes[0, col].plot(energy, old, color=blue, lw=2.5, label='Legado de pySNSPD')
        axes[0, col].plot(energy, new, color=orange, lw=1.4, ls='--',
                          marker='o', markevery=90, ms=3.8, mfc='white',
                          label='Experimental R2')
        axes[0, col].set_title(rf'$|\Delta|/\Delta_0={delta:g}$; $\Gamma/\Delta_0={gamma:g}$')
        axes[0, col].set_ylabel('DOS electrónica / DOS normal')
        axes[0, col].legend(fontsize=11)
        axes[1, col].plot(energy, diff, color=teal, lw=1.1)
        axes[1, col].axhline(0, color='#718096', lw=.7)
        axes[1, col].set_ylabel('Nueva - legado')
        axes[1, col].set_xlabel(r'Energía electrónica $E/\Delta_0$')
        axes[1, col].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
        if np.max(abs(diff)) == 0:
            axes[1, col].set_yticks([0])
            axes[1, col].text(.5, .70, 'Diferencia cero en todas las muestras',
                               ha='center', transform=axes[1, col].transAxes, fontsize=10)
        rows.append({'delta': delta, 'gamma': gamma, 'eta': old_eta,
                     'samples': len(energy), 'max_abs_DOS_difference': float(np.max(abs(diff))),
                     'max_scaled_DOS_difference': float(np.max(abs(diff)/np.maximum(1, abs(new))))})
    fig.suptitle(r'Control: mismo campo y mismo regulador $\eta/\Delta_0=10^{-3}$', fontsize=14)
    for ext in ('png', 'pdf'):
        fig.savefig(figures/f'electronic_same_equation.{ext}')
    plt.close(fig)

    # Gamma=0 has an analytic retarded branch. A logarithmic horizontal axis
    # resolves both regulators without discarding the sharp peak at E=Delta.
    offsets = np.geomspace(1e-10, .2, 1600)
    eta_exponent = int(np.floor(np.log10(args.eta)))
    eta_mantissa = args.eta / 10.**eta_exponent
    eta_label = (rf'10^{{{eta_exponent}}}' if np.isclose(eta_mantissa, 1.)
                 else rf'{eta_mantissa:g}\times 10^{{{eta_exponent}}}')
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.15), constrained_layout=True)
    for sign, ax, name in [(1, axes[0], 'Por encima del borde'),
                           (-1, axes[1], 'Por debajo del borde')]:
        for eta, color, label, style in [
                (old_eta, blue, r'Anterior: $\eta/\Delta_0=10^{-3}$', '-'),
                (args.eta, orange, rf'R2: $\eta/\Delta_0={eta_label}$', '--')]:
            rho = retarded_spectrum(delta+sign*offsets, delta=delta,
                                     gamma=0., eta=eta)[0].real
            ax.loglog(offsets, rho, color=color, lw=2, ls=style, label=label)
        ax.set_xlabel(r'Distancia al borde $|E-|\Delta||/\Delta_0$')
        ax.set_ylabel('DOS electrónica / DOS normal')
        ax.set_title(name)
        ax.grid(which='major', alpha=.18)
        ax.legend(fontsize=11)
    fig.suptitle('Efecto de reducir el regulador numérico $\\eta$\n'
                 r'$|\Delta|/\Delta_0=0.72$; $\Gamma=0$; misma ecuación de Usadel', fontsize=13)
    for ext in ('png', 'pdf'):
        fig.savefig(figures/f'electronic_regulator_edge.{ext}')
    plt.close(fig)
    record = {
        'schema': 'pysnspd.stage1_r2.electronic_comparison.v1',
        'quantity': 'Electronic DOS, Re(c), normalized by the normal electronic DOS',
        'not_the_Simon_input': 'Simon nbn-a2f-ph.dat contains a phonon DOS and alpha2F, not these electronic spectra.',
        'same_equation_controls': rows,
        'regulator_comparison': {'delta': delta, 'gamma': 0, 'legacy_eta': old_eta,
                                  'r2_eta': args.eta,
                                  'interpretation': 'Numerical causal regulator only; no fitted physical lifetime or new material DOS is claimed.'},
        'source_sha256': {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                          for name in ['pysnspd/usadel/solver.py', 'pysnspd/experimental/energy_catalog.py',
                                       'sandbox/stage1_catalog_r2/electronic_comparison.py']},
        'elapsed_seconds': time.perf_counter()-started,
    }
    (out/'electronic_comparison.json').write_text(json.dumps(record, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
