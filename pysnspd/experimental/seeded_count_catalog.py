"""Accelerate the causal inverse using a predictor, retaining the exact equations.

A local polynomial supplies ONLY initial trial energies. Every returned energy
is corrected by bracketed Newton against the original causal state count, with
the same tolerance and bracket as energy_at_count_batch. Forces use implicit
derivatives of that count. An inaccurate/nonmonotone predictor is never used as
the productive energy; state ordering is checked after solving the count.
"""
from dataclasses import dataclass
import numpy as np
from numpy.polynomial.chebyshev import chebval2d
from .energy_catalog import energy_at_count,retarded_spectrum_batch,retarded_spectrum


def invert_count_seeded(count,delta,gamma,eta,predictor):
    x=np.asarray(count,float);seed=np.asarray(predictor,float)
    if (x.ndim!=1 or seed.shape!=x.shape or np.any(~np.isfinite(x)) or np.any(x<=0)
            or not np.isfinite(delta+gamma+eta) or delta<=0 or gamma<=0 or eta<=0):
        raise ValueError('positive finite weak-field scales and count coordinates required')
    lo=x.copy();hi=energy_at_count(x,delta=delta,gamma=0,eta=eta)
    use=np.isfinite(seed)&(seed>lo)&(seed<hi)
    energy=np.where(use,seed,hi)
    active=np.ones(len(x),bool);history=[]
    for iteration in range(80):
        if not np.any(active):
            if np.any(np.diff(energy)<=0) or np.any(~np.isfinite(energy)):
                raise FloatingPointError('corrected causal energies lost finite state ordering')
            return energy,dict(iterations=iteration,spectrum_points=sum(history),active_counts=history,
                               seed_bracket_replacements=int(np.count_nonzero(~use)))
        indices=np.flatnonzero(active);history.append(len(indices));e=energy[indices]
        c,s=retarded_spectrum_batch(e,delta=delta,gamma=gamma,eta=eta)
        w=1j*delta/s
        actual=np.real(w-1j*gamma*delta*delta/(2*w*w))
        residual=actual-x[indices]
        hi[indices]=np.where(residual>=0,e,hi[indices])
        lo[indices]=np.where(residual<0,e,lo[indices])
        step=residual/c.real;tolerance=4e-14*np.maximum(1,abs(e))
        done=(abs(step)<=tolerance)|((hi[indices]-lo[indices])<=tolerance)
        active[indices[done]]=False;candidate=e-step
        safe=(candidate>lo[indices])&(candidate<hi[indices])&np.isfinite(candidate)
        finish=done&safe;energy[indices[finish]]=candidate[finish]
        candidate=np.where(safe,candidate,(lo[indices]+hi[indices])/2)
        energy[indices[~done]]=candidate[~done]
    raise FloatingPointError('seeded causal inverse did not converge within80iterations')


@dataclass(frozen=True)
class SeededCountCatalog:
    source: object
    predictor: object

    def __post_init__(self):
        for key in ['count_nodes','count_weights']:
            if not np.array_equal(getattr(self.source,key),getattr(self.predictor,key)):
                raise ValueError('predictor and source count quadratures differ')
        if self.source.eta!=self.predictor.eta:raise ValueError('causal regulators differ')

    @property
    def vacuum(self):return self.source.vacuum
    @property
    def count_nodes(self):return self.source.count_nodes
    @property
    def count_weights(self):return self.source.count_weights
    @property
    def eta(self):return self.source.eta

    def energy_kernel_with_diagnostics(self,amplitude,gamma):
        a,g=float(amplitude),float(gamma);patch=self.predictor
        if not (np.isfinite(a) and np.isfinite(g) and
                patch.amplitude_bounds[0]<=a<=patch.amplitude_bounds[1] and
                patch.gamma_bounds[0]<=g<=patch.gamma_bounds[1]):
            raise ValueError('field outside the registered predictor box')
        self.vacuum.evaluate(a,g)
        x=2*(a-patch.amplitude_bounds[0])/np.diff(patch.amplitude_bounds)[0]-1
        y=2*(g-patch.gamma_bounds[0])/np.diff(patch.gamma_bounds)[0]-1
        # Deliberately bypass predictor.energy_kernel: no predictor force,
        # positivity or state-order claim is needed for an internal trial.
        seed=chebval2d(x,y,patch.coefficients)
        energy,diagnostics=invert_count_seeded(self.count_nodes,a,g,self.eta,seed)
        c,s=retarded_spectrum(energy,delta=a,gamma=g,eta=self.eta)
        da=s.imag/c.real;dg=-s.real*da
        values=(energy,da,dg)
        if np.any(~np.isfinite(values)):raise FloatingPointError('nonfinite corrected kernels')
        return values,diagnostics

    def energy_kernel(self,amplitude,gamma):
        return self.energy_kernel_with_diagnostics(amplitude,gamma)[0]

    def evaluate(self,amplitude,gamma,occupation):
        raw=np.asarray(occupation)
        if np.iscomplexobj(raw):raise ValueError('occupation must be real')
        p=np.asarray(raw,float)
        if p.shape!=self.count_nodes.shape or np.any(~np.isfinite(p)) or np.any((p<0)|(p>1)):
            raise ValueError('physical occupation required on unchanged native count grid')
        return tuple(float(u+4*np.dot(self.count_weights*k,p)) for u,k in
                     zip(self.vacuum.evaluate(amplitude,gamma),self.energy_kernel(amplitude,gamma)))
