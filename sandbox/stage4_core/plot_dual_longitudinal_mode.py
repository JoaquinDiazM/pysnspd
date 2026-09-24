"""Plot the saved invariant-mode control; no spectrum or trajectory is solved."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-root', type=Path, required=True)
    parser.add_argument('--output-root', type=Path, required=True)
    args = parser.parse_args()
    receipt = json.loads((args.input_root/'receipt.json').read_text())
    with np.load(args.input_root/'evolution.npz') as saved:
        data = {name:saved[name].copy() for name in saved.files}
    args.output_root.mkdir(parents=True, exist_ok=True)
    t, a = data['times_ps'], data['amplitudes']
    energy = data['energy_over_kBTc']
    temperature = receipt['T_K']/receipt['Tc_K']
    q = np.exp(-energy/temperature)
    chi = 2*q/(temperature*(1+q)**2)
    plt.rcParams.update({'font.size':10, 'axes.titlesize':11, 'axes.labelsize':10,
        'legend.fontsize':8, 'savefig.facecolor':'white'})
    fig, ax = plt.subplots(2, 2, figsize=(12, 8), layout='constrained')
    coords = data['coordinates_bar']
    points = ax[0, 0].scatter(coords[:, 0], coords[:, 1], c=data['mode'],
        s=5, cmap='viridis', vmin=0, vmax=1, rasterized=True)
    fig.colorbar(points, ax=ax[0, 0], label=r'Modo espacial $v(X)$, máximo = 1')
    ax[0, 0].set(xlabel=r'$x/\ell_0$', ylabel=r'$y/\ell_0$',
        title='A. Preparación suave sobre todos los nodos duales', aspect='equal')
    scale = 100/receipt['gap_over_kBTc']
    ax[0, 1].plot(t, scale*a[:, 0], label='Euler: paso principal', color='#2461A4')
    ax[0, 1].plot(data['refined_times_ps'], scale*data['refined_amplitudes'][:, 0],
        '--', label='Euler: mitad de paso', color='#E3A323')
    ax[0, 1].plot(t[-1], scale*data['exact_final'][0], 'o', ms=5,
        color='#B1293D', label='Exponencial matricial: valor final')
    ax[0, 1].set(xlabel='Tiempo físico (ps)', ylabel=r'$100\,\delta|\Delta|/|\Delta_{eq}|$ (%)',
        title='B. Cambio del gap en el máximo del modo')
    ax[0, 1].legend()
    difference = a[:, 0]-data['refined_amplitudes'][::2, 0]
    inset = ax[0, 1].inset_axes([.58, .12, .38, .26])
    inset.plot(t, 1e6*difference/receipt['gap_over_kBTc'], color='#6D3A7B')
    inset.set_title('Diferencia entre pasos', fontsize=8)
    inset.set_ylabel('ppm del gap eq.', fontsize=7)
    inset.tick_params(labelsize=7)
    inset.grid(alpha=.15)
    ax[1, 0].plot(energy, -.5*chi*a[0, 1:], 'o-', label='Inicio', color='#2461A4')
    ax[1, 0].plot(energy, -.5*chi*a[-1, 1:], 's--', label=f'Final: {t[-1]:g} ps', color='#B1293D')
    ax[1, 0].axhline(0, color='.5', lw=.6)
    ax[1, 0].set(xlabel=r'Energía $E/(k_BT_c)$', ylabel=r'$\delta f(E)= -\chi(E)y(E)/2$',
        title='C. Cambio de ocupación en el máximo del modo')
    ax[1, 0].ticklabel_format(axis='y', style='sci', scilimits=(-3, 3))
    ax[1, 0].legend()
    available = data['availability']/data['availability'][0]
    loss = data['integrated_losses']/data['availability'][0]
    ax[1, 1].plot(t, available, label=r'Disponibilidad $A(t)/A(0)$', color='#2461A4')
    ax[1, 1].plot(t, loss[:, 0], label='Pérdida KWT acumulada / A(0)', color='#E3A323')
    ax[1, 1].plot(t, loss[:, 1], label='Pérdida por transporte / A(0)', color='#589157')
    ax[1, 1].plot(t, available+loss.sum(axis=1), '--', label='Suma (balance esperado: 1)', color='#B1293D')
    ax[1, 1].set(xlabel='Tiempo físico (ps)', ylabel='Fracción de la disponibilidad inicial',
        title='D. Energía libre de la perturbación y sus pérdidas')
    ax[1, 1].legend()
    ax[1, 1].text(.03, .06, f'Transporte final: {loss[-1, 1]:.2e} A(0)',
        transform=ax[1, 1].transAxes, fontsize=8, color='#386A38')
    for axis in ax.ravel():
        axis.grid(alpha=.15)
    fig.suptitle('Intercambio longitudinal débil, sin fotón ni corriente aplicada\n'
        r'$\delta\Delta/(k_BT_c)=x(t)v(X)$; $\delta h_L(E,X)=\chi(E)y(E,t)v(X)$', fontsize=13)
    fig.savefig(args.output_root/'longitudinal_exchange.png', dpi=180)
    fig.savefig(args.output_root/'longitudinal_exchange.pdf')
    plt.close(fig)
    text = ('La figura usa el modo propio más suave del Laplaciano de la malla dual de '
        f"{receipt['nodes']} nodos, con contactos de extremo fijos y bordes laterales naturales. "
        'A muestra su valor normalizado, no una distribución fotónica. B muestra la '
        'variación real del gap en el punto donde el modo alcanza uno, respecto del '
        'gap uniforme de equilibrio. C representa incrementos de la ocupación '
        'electrónica respecto de Fermi–Dirac a 0,9 K; las líneas sólo unen las siete '
        'energías de la cuadratura de control. D usa disponibilidad cuadrática de la '
        'perturbación y pérdidas integradas, todas divididas por A(0). No representa '
        'energía interna total ni calor entregado a una población de fonones. La '
        'referencia matricial de B sólo contrasta el Euler del sector lineal; no '
        'es otro solver del detector. No se ajustó una temperatura a la población.\n')
    (args.output_root/'longitudinal_exchange_caption.md').write_text(text, encoding='utf-8')
    print(json.dumps({'plots':'longitudinal_exchange.png/.pdf', 'new_physics_solves':0}))


if __name__ == '__main__':
    main()
