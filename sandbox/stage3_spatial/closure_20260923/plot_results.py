"""Plot saved stage-3 measurements; never evaluate the physical model."""
from pathlib import Path
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT/'docs/implementation/stage3/closure_20260923'
RAW = DATA/'raw/stage3_reservoir_full_20260923'
OUT = DATA/'figures'


def main():
    OUT.mkdir(exist_ok=True)
    path = RAW/'reservoir/summary.json'
    result = json.loads(path.read_text(encoding='utf-8'))
    assert result['status'] == 'PASS_SAVED_RESERVOIR_SNAPSHOTS_SCOPE_LIMITED'
    cases = {c['case_id']: c for c in result['cases']}
    h = [cases[f'm{i}_helix'] for i in range(3)]
    p = [cases[f'm{i}_smooth_perturbation'] for i in range(3)]
    x = np.arange(3)
    plt.rcParams.update({'font.size': 11, 'axes.spines.top': False,
                         'axes.spines.right': False, 'axes.titlepad': 8})
    blue, orange = '#246786', '#b95828'
    fig, axes = plt.subplots(1, 2, figsize=(7, 3.1))
    heat = np.array([c['reservoir']['condensate_rate_bar'] for c in h])
    normal = 100*np.array([max(abs(np.array(c['terminal']['adjacent_internal_normal_current_over_Iref']))) for c in h])
    for ax, values, title, label in zip(axes, [heat, normal],
        ['Calor de la hélice', 'Corriente normal adyacente'],
        [r'$Q_\Delta$ integrado (normalizado)', r'$|I_n|/I_{ref}$ (%)']):
        ax.semilogy(x, values, 'o-', color=blue)
        ax.set_xticks(x, ['180', '90', '45'])
        ax.set_xlabel('Elemento longitudinal (nm)', fontsize=10)
        ax.set_title(title, fontsize=11)
        ax.set_ylabel(label)
        ax.grid(axis='y', alpha=.25)
    fig.subplots_adjust(left=.12, right=.98, top=.86, bottom=.24, wspace=.55)
    fig.text(.5, .015, 'Tres mallas; el flujo normal externo completo de D.27 sigue pendiente.', ha='center', fontsize=10)
    fig.savefig(OUT/'01_boundary_mesh.png', dpi=300)
    plt.close(fig)
    q = np.array([c['reservoir']['condensate_rate_bar'] for c in p])
    v = np.array([c['reservoir']['maximum_material_speed_bar'] for c in p])
    excess = np.array([b['potential']['Vdev_V']-a['potential']['Vdev_V'] for a,b in zip(h,p)])
    observables = {'QDelta':q, 'maximum_material_speed':v, 'excess_Vdev':excess}
    relative = 100*np.array([abs((a[-2]-a[-1])/a[-1]) for a in observables.values()])
    fig, axes = plt.subplots(1, 2, figsize=(7, 3.1))
    axes[0].plot(x, q/q[-1], 'o-', label=r'$Q_\Delta$', color=blue)
    axes[0].plot(x, v/v[-1], 's-', label=r'$v_{material,max}$', color=orange)
    axes[0].set_xticks(x, ['180', '90', '45'])
    axes[0].set_xlabel('Tamaño de elemento (nm)')
    axes[0].set_ylabel('Valor / valor en malla fina')
    axes[0].set_title('Respuesta a la perturbación', fontsize=11)
    axes[0].legend(fontsize=9)
    axes[0].grid(axis='y', alpha=.25)
    labels = [r'$Q_\Delta$', r'$v_{material,max}$', r'$V_{dev,p}-V_{dev,h}$']
    axes[1].barh(labels, relative, color=[blue, orange, '#458b7c'])
    axes[1].set_xlim(0, max(relative)*1.35)
    for i, value in enumerate(relative):
        axes[1].text(value+.02, i, f'{value:.3f} %', va='center', fontsize=10)
    axes[1].set_xlabel('Diferencia media-fina (%)')
    axes[1].set_title('Cambio media-fina', fontsize=11)
    axes[1].grid(axis='x', alpha=.25)
    fig.subplots_adjust(left=.11, right=.98, top=.84, bottom=.26, wspace=.65)
    fig.text(.5, .012, 'Instantáneas: las diferencias no constituyen cotas de error del transiente.', ha='center', fontsize=10)
    fig.savefig(OUT/'02_perturbation_mesh.png', dpi=300)
    plt.close(fig)
    metrics = dict(source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        reservoir_helix_heat=heat.tolist(), adjacent_normal_relative_percent=normal.tolist(),
        perturbation_heat=q.tolist(), maximum_material_speed=v.tolist(),
        excess_Vdev_V=excess.tolist(), last_mesh_difference_percent=dict(zip(observables,relative)),
        interpretation='Saved instantaneous sensitivity only; no time or universal spatial error certificate.')
    (OUT/'figure_data.json').write_text(json.dumps(metrics, indent=2)+'\n', encoding='utf-8')


if __name__ == '__main__':
    main()
