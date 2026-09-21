"""Render convergence and a fixed-field profile; no dynamics or metric edits."""
from pathlib import Path
import hashlib
import json
import sys
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'docs/implementation/stage2'
sys.path.insert(0,str(ROOT))
from transport_diagnostics import population,FIELDS
from pysnspd.experimental.energy_catalog import OccupationEnergyCatalog,retarded_spectrum
from pysnspd.experimental.cell_validation import ElectronicCell
from pysnspd.experimental.refined_cells import refined_count_catalog


def profile():
    """Evaluate the existing geometric rate density on common interior energies."""
    catalogue=ROOT/'docs/implementation/stage1_r2/catalogs/occupation_catalog.npz'
    base=OccupationEnergyCatalog.load(catalogue)
    meshes=[base,refined_count_catalog(base,1,bulk_spacing=.02)]
    fields=FIELDS['gapped']; query=np.linspace(1.4,3.8,801)
    diffusion=[]
    for amplitude,gamma in fields:
        c,s=retarded_spectrum(query,delta=amplitude,gamma=gamma,eta=base.eta)
        diffusion.append(c.real*c.real-s.imag*s.imag)
    left,right=diffusion
    conductance_density=4*.03*2*left*right/(left+right)
    reference=conductance_density*(population(query,'nonthermal',0)-population(query,'nonthermal',1))
    curves=[]
    for mesh in meshes:
        products=[]
        for side,field in enumerate(fields):
            cell=ElectronicCell(mesh,*field)
            p=population(cell.energies,'nonthermal',side)
            upper=np.searchsorted(cell.energies,query); lower=upper-1
            if np.any(lower<0) or np.any(upper>=len(cell.energies)):
                raise ValueError('plot would extrapolate')
            fraction=(query-cell.energies[lower])/(cell.energies[upper]-cell.energies[lower])
            la=(1-fraction)*np.log(p[lower])+fraction*np.log(p[upper])
            lb=(1-fraction)*np.log1p(-p[lower])+fraction*np.log1p(-p[upper])
            products.append((np.exp(la),np.exp(lb)))
        (al,bl),(ar,br)=products
        curves.append(conductance_density*(al*br-ar*bl))
    return query,reference,curves


def main():
    started=time.perf_counter()
    paths=[OUT/'pilot_initial/transport_original_results.json',OUT/'transport_results.json']
    old,new=[json.loads(p.read_text(encoding='utf-8')) for p in paths]
    plt.rcParams.update({'font.size':11,'axes.labelsize':11,'axes.titlesize':12,'legend.fontsize':10.5})
    fig=plt.figure(figsize=(10,6.8),layout='constrained')
    grid=fig.add_gridspec(2,2,width_ratios=[1,1.2])
    convergence=fig.add_subplot(grid[:,0]); curve=fig.add_subplot(grid[0,1]); residual=fig.add_subplot(grid[1,1],sharex=curve)
    summaries=new['mesh_summary']
    counts=[r['electron_nodes'][0] for r in summaries]
    for key,label,marker in [('maximum_RHS_relative_error','Respuesta de población','o'),('maximum_power_relative_error','Potencia transferida','s')]:
        convergence.loglog(counts,[100*r[key] for r in summaries],marker+'-',label=label)
    convergence.plot(180,100*old['native_max_relative_RHS_error'],'x',ms=10,mew=2,color='#a64e36',label='R2 histórico: respuesta')
    convergence.set(xlabel='Estados electrónicos por celda',ylabel='Error máximo de los seis casos [%]',title='a) Convergencia del transporte',ylim=(8e-5,1.4),xlim=(145,2500))
    convergence.set_xticks([180,*counts],[str(n) for n in [180,*counts]],rotation=25); convergence.minorticks_off()
    convergence.legend(loc='upper right')
    convergence.axhline(.1,color='#555555',ls='--',lw=1.2)
    convergence.text(195,.072,'Límite: 0,1 %',color='#555555')
    query,reference,curves=profile(); scale=float(np.max(reference))
    curve.semilogy(query,reference/scale,color='black',lw=2,label='Referencia continua')
    for values,color,style,label in zip(curves,['#a64e36','#226e96'],['--',':'],['R2: 180 estados','Complementario: 526']):
        curve.semilogy(query,values/scale,color=color,ls=style,lw=1.8,label=label)
        residual.plot(query,100*(values-reference)/reference,color=color,ls=style,lw=1.4,label=label)
    curve.set(ylabel='Flujo J(E) / máximo de referencia',title='b) Perfil no térmico; brechas distintas')
    curve.legend(loc='lower left'); curve.tick_params(labelbottom=False)
    residual.axhline(0,color='black',lw=1)
    local_max=100*float(np.max(np.abs((curves[1]-reference)/reference)))
    residual.text(.98,.06,f'526: máximo local {local_max:.4f} %',transform=residual.transAxes,
                  ha='right',color='#226e96',bbox={'facecolor':'white','edgecolor':'none','alpha':.9})
    residual.set(xlabel='Energía E/Δ₀',ylabel='Diferencia local respecto a Jref [%]',title='c) Residuos de las curvas superpuestas')
    for ax in [convergence,curve,residual]:
        ax.grid(axis='y',alpha=.2)
    figure=OUT/'figures/transport_convergence'
    fig.savefig(figure.with_suffix('.png'),dpi=180)
    fig.savefig(figure.with_suffix('.pdf'))
    plt.close(fig)
    sources=[ROOT/'pysnspd/experimental/energy_catalog.py',ROOT/'pysnspd/experimental/refined_cells.py',ROOT/'sandbox/stage2_cells/transport_diagnostics.py']
    curves_path=OUT/'transport_profile_curves.json'
    curves_path.write_text(json.dumps({'schema':'pysnspd.stage2.transport-profile.v1',
        'field':FIELDS['gapped'],'profile':'registered nonthermal','energy_over_Delta0':query.tolist(),
        'reference_flux_density':reference.tolist(),'native180_flux_density':curves[0].tolist(),
        'complementary526_flux_density':curves[1].tolist(),
        'complementary526_maximum_local_relative_error_percent':local_max,
        'scope':'Illustration on common interior energies; admission metrics remain in transport_results.json.'},indent=2)+'\n',encoding='utf-8')
    manifest={'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'inputs':{p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'sources':{p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        'profile_curves_sha256':hashlib.sha256(curves_path.read_bytes()).hexdigest(),
        'outputs':{figure.with_suffix(s).relative_to(ROOT).as_posix():hashlib.sha256(figure.with_suffix(s).read_bytes()).hexdigest() for s in ['.png','.pdf']},
        'runtime_seconds':time.perf_counter()-started}
    (OUT/'transport_figure_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':main()
