"""Lightweight prescribed-resistance comparison, not a detector transient.

Fixed before evaluating: Ib0=20uA, Rload=50ohm, Rbias=10kohm,
Lbias=1uH, external series L=10nH, C=100pF, Vbias=Rbias*Ib0.
Rdev=0/150/0ohm during0-1/1-4/4-16ns. Constant-section affine
matrix exponentials serve as the reference. No temporal tolerance is fitted.
"""
from pathlib import Path
import hashlib,json,time
import numpy as np
from scipy.linalg import expm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[3]
DATA=ROOT/'docs/implementation/stage2/practical_review_20260922'
FIG=DATA/'figures';FIG.mkdir(parents=True,exist_ok=True)
I0=20e-6;RL=50.;Rb=1e4;Lb=1e-6;L=10e-9;C=100e-12;Vb=I0*Rb
ts=L/RL;vs=I0*RL

def generator(resistance):
    # Dimensionless state [Ib/I0, Is/I0, vc/(I0*RL), 1].
    matrix=np.zeros((4,4))
    matrix[0]=[ -ts*(Rb+RL)/Lb,ts*RL/Lb,-ts*RL/Lb,ts*Rb/Lb ]
    matrix[1]=[ ts*RL/L,-ts*(RL+resistance)/L,ts*RL/L,0 ]
    matrix[2]=[ts/(C*RL),-ts/(C*RL),0,0]
    return matrix

started=time.perf_counter();states=[np.array([1.,1.,0.,1.])];times=[0.];resistances=[0.];simple=[1.]
for left,right,resistance in [(0.,1e-9,0.),(1e-9,4e-9,150.),(4e-9,16e-9,0.)]:
    number=int(round((right-left)/5e-12));h=(right-left)/number
    propagator=expm(generator(resistance)*h/ts)
    equilibrium=RL/(RL+resistance);decay=np.exp(-(RL+resistance)*h/L)
    for j in range(number):
        states.append(propagator@states[-1]);times.append(left+(j+1)*h)
        resistances.append(resistance);simple.append(equilibrium+(simple[-1]-equilibrium)*decay)
x=np.array(states);t=np.array(times);rd=np.array(resistances);is_old=I0*np.array(simple)
ib=I0*x[:,0];is_=I0*x[:,1];vc=vs*x[:,2];out=RL*(ib-is_);vnode=vc+out;vdev=rd*is_
energy=.5*Lb*ib**2+.5*L*is_**2+.5*C*vc**2
power_residual=[]
for z,r in zip(x,rd):
    rhs=generator(r)@z/ts
    dib=I0*rhs[0];dis=I0*rhs[1];dvc=vs*rhs[2]
    b,s,c=I0*z[0],I0*z[1],vs*z[2]
    udot=Lb*b*dib+L*s*dis+C*c*dvc
    terms=[udot,r*s*s,Rb*b*b,RL*(b-s)**2,Vb*b]
    power_residual.append(abs(sum(terms[:-1])-terms[-1])/max(sum(abs(q) for q in terms),1e-30))
result=dict(status='MEASURED_PRESCRIBED_CIRCUIT_IDENTITY_NOT_ADMISSION',runtime_seconds=time.perf_counter()-started,
    parameters=dict(I0_A=I0,Rload_ohm=RL,Rbias_ohm=Rb,Lbias_H=Lb,Lseries_external_H=L,C_F=C,Vbias_V=Vb),
    peak_full_V=float(out.max()),peak_simple_V=float((RL*(I0-is_old)).max()),
    minimum_full_V=float(out.min()),final_full_V=float(out[-1]),
    max_instantaneous_power_residual_scaled=float(max(power_residual)),
    source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    method='Exact affine matrix-exponential propagation on each constant prescribed-resistance interval.',
    scope='Topology diagnostic with identical reference current and series inductance. Synthetic voltage-bias value; Lseries10nH is an external test element, not an identified partition of the thesis total. No film, photon, detection latency, readout fit or production change.')
(DATA/'circuit_diagnostic.json').write_text(json.dumps(result,indent=2)+'\n')
np.savetxt(DATA/'circuit_diagnostic.csv',np.column_stack([t,ib,is_,vc,out,vnode,vdev,energy,RL*(I0-is_old)]),delimiter=',',header='t_s,Ib_A,Is_A,vc_V,Vload_V,Vnode_V,Vdev_V,Ucircuit_J,Vload_simple_V',comments='')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
fig,axes=plt.subplots(1,2,figsize=(10,3.4),layout='constrained')
axes[0].plot(t*1e9,out*1e3,color='#007e80',lw=2,label='Topología de la memoria (prueba)')
axes[0].plot(t*1e9,RL*(I0-is_old)*1e3,color='#cf791b',ls='--',lw=1.8,label='Fuente ideal + carga')
axes[0].axhline(0,color='#666666',lw=.7);axes[0].set(xlabel='Tiempo (ns)',ylabel='Voltaje en carga (mV)',title='Red completa frente a fuente ideal')
axes[0].legend(fontsize=8)
axes[1].plot(t*1e9,vc*1e3,label='Voltaje del capacitor',color='#6c539e')
axes[1].plot(t*1e9,vnode*1e3,label='Voltaje del nodo',color='#007e80')
axes[1].plot(t*1e9,out*1e3,label='Voltaje en la carga',color='#3274b2',ls=':')
axes[1].set(xlabel='Tiempo (ns)',ylabel='Voltaje (mV)',title='Tres voltajes con funciones distintas');axes[1].legend(fontsize=8)
for ax in axes:ax.grid(alpha=.2)
fig.savefig(FIG/'circuit_comparison.png');plt.close(fig)
print(json.dumps(result))
