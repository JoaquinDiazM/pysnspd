"""Cancellation-aware differences of the same alpha=0 thermal graph action.

These identities compare nearby nonlinear stationary spectral states. They
do not approximate the gap dependence or replace the underlying functional.
"""
import numpy as np


def spectral_difference(graph,d0,u0,d,u,epsilon):
    n=graph.n_nodes
    d0,u0,d,u=(np.asarray(value,complex) for value in (d0,u0,d,u))
    if any(value.shape!=(n,) or np.any(~np.isfinite(value)) for value in (d0,u0,d,u)) or not np.isfinite(epsilon) or epsilon<=0:
        raise ValueError('Finite full-node states and positive Matsubara energy required')
    gap_change=d-d0;du=u-u0;s0=np.hypot(1.,abs(u0));s=np.hypot(1.,abs(u))
    g0=1/s0;g=1/s
    dg=-g*g0*(2*np.real(np.conj(u0)*du)+abs(du)**2)/(s+s0)
    f0=g0*u0;df=g*du+u0*dg
    mass=graph.area_weights;tail,head=graph.edges.T;c=graph.conductance
    spectral_node=np.sum(mass*(-2*epsilon*dg-2*np.real(np.conj(d0)*df+np.conj(gap_change)*(f0+df))))
    edgef=f0[tail]-f0[head];deltaf=df[tail]-df[head]
    edgeg=g0[tail]-g0[head];deltag=dg[tail]-dg[head]
    spectral_edge=np.sum(c*(2*np.real(np.conj(edgef)*deltaf)+abs(deltaf)**2+2*edgeg*deltag+deltag**2))
    amplitude_change=2*np.real(np.conj(d0)*gap_change)+abs(gap_change)**2
    current_change=2*c*np.imag(np.conj(df[tail])*f0[head]+np.conj(f0[tail])*df[head]+np.conj(df[tail])*df[head])
    return dict(f_difference=df,g_difference=dg,
        renormalized_energy_difference=float(spectral_node+spectral_edge+np.dot(mass,amplitude_change)/epsilon),
        gap_gradient_difference=2*mass*(gap_change/epsilon-df),current_difference=current_change)
