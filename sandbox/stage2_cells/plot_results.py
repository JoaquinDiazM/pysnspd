"""Results plots with explicit scales, separate residuals and readable legends."""
from pathlib import Path
import argparse
import json
import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.special import expit
from pysnspd.experimental.energy_catalog import energy_at_count_batch

BLUE, ORANGE, TEAL, GRAY = '#2166ac', '#d66024', '#00897b', '#5d6470'


def load(path):
    return json.loads(path.read_text(encoding='utf-8'))


def polish(ax):
    ax.spines[['top','right']].set_visible(False)
    ax.grid(alpha=.17)


def save(fig, path):
    fig.savefig(path.with_suffix('.png'), dpi=190, facecolor='white')
    fig.savefig(path.with_suffix('.pdf'), facecolor='white')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--one', type=Path, required=True)
    parser.add_argument('--two', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=ROOT/'docs/implementation/stage2/figures')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.size':10, 'axes.titlesize':11, 'legend.fontsize':8.5})
    one, two = load(args.one), load(args.two)
    fig, axes = plt.subplots(2,2,figsize=(10,6.4),layout='constrained')
    for record, label, colors, style in ((one,'Aislada',[GRAY],'--'),
                                        (two,'Acopladas',[BLUE,ORANGE],'-')):
        t = np.asarray(record['times'])
        states = record['snapshots']
        for cell, color in enumerate(colors):
            suffix = '' if label == 'Aislada' else f' {cell+1}'
            axes[0,0].plot(t,[v['amplitudes'][cell] for v in states],style,color=color,label=label+suffix)
            axes[0,1].plot(t,[v['temperatures'][cell] for v in states],style,color=color)
        baseline = states[0]['conserved']
        error = np.asarray([abs(v['conserved']-baseline) for v in states])
        axes[1,1].plot(t,error,color=colors[0],linestyle=style,label=label)
    t=np.asarray(two['times']); states=two['snapshots']
    for name, label, color in (('condensate_heat','Condensado → electrones',ORANGE),
                               ('electron_to_phonon','Electrones → fonones',BLUE),
                               ('escape','Fonones → baño',TEAL)):
        axes[1,0].plot(t,[sum(v[name]) for v in states],label=label,color=color)
    axes[1,0].plot(t,[v['transport'] for v in states],label='Transporte: celda 1 → 2',color=GRAY,linestyle='--')
    axes[0,0].set(ylabel=r'$|\Delta|/\Delta_0$',title='Amplitud: recuperación con cinética activa')
    axes[0,1].set(ylabel=r'$k_B T_E/\Delta_0$',title='Temperatura equivalente de la energía electrónica')
    axes[1,0].set(ylabel=r'Energía / $(N_0\Delta_0^2)$',title='Transferencias acumuladas en las dos celdas')
    axes[1,1].set(ylabel=r'Error absoluto / $(N_0\Delta_0^2)$',title='Balance con fuente y escape contabilizados')
    axes[1,1].set_yscale('symlog',linthresh=1e-14)
    axes[1,1].axhline(1e-7,color='#b2182b',linestyle=':',label='Límite de aceptación')
    for ax in axes.ravel():
        polish(ax); ax.set_xlabel(r'$t/t_{ref}$  (escala de ensayo: 1 ps)')
    for index in ((0,0),(1,0),(1,1)): axes[index].legend(loc='best')
    save(fig,args.output/'coupled_dynamics')

    with np.load(args.one.with_suffix('.npz')) as data:
        x=data['electron_count']; om=data['phonon_energies']; y=data['states']
    ne=len(x); p0=y[0,1:1+ne]; pf=y[-1,1:1+ne]
    n0=y[0,1+ne:1+ne+len(om)]; nf=y[-1,1+ne:1+ne+len(om)]
    a0=one['initial']['amplitudes'][0]; af=one['final']['amplitudes'][0]
    e0=energy_at_count_batch(x,delta=a0,gamma=0.,eta=1e-8)
    ef=energy_at_count_batch(x,delta=af,gamma=0.,eta=1e-8)
    fd=expit(-ef/one['final']['temperatures'][0])
    fig,axes=plt.subplots(2,2,figsize=(10,6.2),layout='constrained')
    axes[0,0].plot(e0,p0,color=BLUE,label='Distribución inicial')
    axes[0,0].plot(ef,pf,color=ORANGE,label='Distribución final p')
    axes[0,0].plot(ef,fd,color=GRAY,linestyle='--',label='Fermi-Dirac con la misma energía')
    axes[0,0].set(xlim=(0,3),ylabel='Ocupación electrónica p',title='Una celda aislada: distribución electrónica')
    axes[1,0].plot(ef,pf-fd,color=TEAL)
    axes[1,0].axhline(0,color=GRAY,linewidth=.8)
    axes[1,0].set(xlim=(0,3),ylabel=r'$p-f_{FD}(T_E)$',title='La temperatura equivalente no sustituye p')
    axes[0,1].plot(om,n0,color=BLUE,label='Distribución inicial')
    axes[0,1].plot(om,nf,color=ORANGE,label='Distribución final n')
    axes[0,1].set(yscale='log',ylabel='Ocupación fonónica n',title='Reacciones redistribuyen los fonones')
    axes[1,1].plot(om,nf-n0,color=TEAL)
    axes[1,1].axhline(0,color=GRAY,linewidth=.8)
    axes[1,1].set(ylabel=r'$n_{final}-n_{inicial}$',title='Cambio de ocupación fonónica')
    for row in range(2):
        axes[row,0].set_xlabel(r'$E/\Delta_0$ (cada curva en su espectro)')
        axes[row,1].set_xlabel(r'$\Omega/\Delta_0$')
    for ax in axes.ravel(): polish(ax)
    axes[0,0].legend(); axes[0,1].legend()
    save(fig,args.output/'resolved_populations')


if __name__ == '__main__':
    main()
