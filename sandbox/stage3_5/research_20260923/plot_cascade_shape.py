"""Plot only the extracted historical A20 percentiles, not a new cascade."""
from pathlib import Path
import json,math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[3]
DATA=ROOT/'docs/implementation/stage3_5/research_20260923'

def main():
    v=json.loads((DATA/'cascade_digitization.json').read_text(encoding='utf8'))
    fig,ax=plt.subplots(figsize=(7.3,3.8))
    for i,row in enumerate(v['rows']):
        vals=[];err=[]
        for key,q in [('R50',.5),('R90',.9)]:
            factor=math.sqrt(-2*math.log(1-q))
            val=row[key]['radius_nm']/factor
            lo,hi=row[key]['extraction_envelope_nm']
            vals.append(val);err.append([val-lo/factor,hi/factor-val])
        ax.errorbar(np.array([0,1])+(i-.5)*.06,vals,yerr=np.array(err).T,
            fmt='o-',capsize=4,color=['#207c91','#b75b26'][i],
            label=f"A20: {row['time_ps']:.4f} ps")
    ax.set_xticks([0,1],['Usando el radio que encierra 50 %','Usando el radio que encierra 90 %'])
    ax.set_ylabel('Ancho gaussiano equivalente s (nm)')
    ax.set_title('Si el perfil fuera gaussiano, ambos anchos coincidirían')
    ax.set_xlim(-.25,1.25);ax.set_ylim(1.25,1.86);ax.grid(axis='y',alpha=.2)
    ax.legend(loc='upper left',fontsize=9)
    ax.spines[['top','right']].set_visible(False)
    fig.text(.5,.015,'Energía total de A20, figura 2.7(a). Barras: lectura de imagen, no incertidumbre física.',ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.05,1,1])
    target=DATA/'figures/04_cascade_shape.png';fig.savefig(target,dpi=220);plt.close(fig)
    print(target)

if __name__=='__main__':main()
