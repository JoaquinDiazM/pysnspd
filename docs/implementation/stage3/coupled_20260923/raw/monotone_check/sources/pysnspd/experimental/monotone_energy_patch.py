"""Positive-increment energy patch; forces derive from its reconstructed E.

Interpolate logarithms of E[0] and adjacent positive energy increments. Summing
their exponentials preserves state ordering without sorting, clipping, adding
states or remapping populations. Accuracy is a separate measured requirement.
"""
from pathlib import Path
import numpy as np
from numpy.polynomial.chebyshev import chebval2d,chebvander2d
from .local_energy_patch import LocalEnergyPatch


class MonotoneEnergyPatch(LocalEnergyPatch):
    representation='log_positive_energy_increments_v1'

    @classmethod
    def build(cls,source,amplitude_bounds,gamma_bounds,degree=6,*,on_sample=None):
        if type(degree) is not int or degree<2:raise ValueError('integer degree >=2 required')
        raw=np.asarray([amplitude_bounds,gamma_bounds])
        if np.iscomplexobj(raw):raise ValueError('field bounds must be real')
        bounds=np.asarray(raw,float)
        if bounds.shape!=(2,2) or np.any(~np.isfinite(bounds)) or np.any(np.diff(bounds,axis=1)<=0):
            raise ValueError('two increasing field intervals required')
        nodes=np.cos(np.pi*np.arange(degree+1)/degree)
        xx,yy=np.meshgrid(nodes,nodes,indexing='ij')
        aa=bounds[0].mean()+np.diff(bounds[0])[0]*xx.ravel()/2
        gg=bounds[1].mean()+np.diff(bounds[1])[0]*yy.ravel()/2
        aa=np.where(xx.ravel()==1,bounds[0,1],np.where(xx.ravel()==-1,bounds[0,0],aa))
        gg=np.where(yy.ravel()==1,bounds[1,1],np.where(yy.ravel()==-1,bounds[1,0],gg))
        values=[]
        for i,(a,g) in enumerate(zip(aa,gg)):
            energy=source.energy_kernel(float(a),float(g))[0]
            increments=np.diff(np.r_[0.,energy])
            if np.any(~np.isfinite(increments)) or np.any(increments<=0):
                raise ValueError('source must have positive strictly ordered energies')
            values.append(np.log(increments))
            if on_sample is not None:on_sample(i+1,len(aa))
        vander=chebvander2d(xx.ravel(),yy.ravel(),[degree,degree])
        coeff=np.linalg.solve(vander,np.asarray(values)).reshape(degree+1,degree+1,-1)
        return cls(source,tuple(bounds[0]),tuple(bounds[1]),coeff)

    def energy_kernel(self,amplitude,gamma):
        a,g=float(amplitude),float(gamma)
        if not (np.isfinite(a) and np.isfinite(g)
                and self.amplitude_bounds[0]<=a<=self.amplitude_bounds[1]
                and self.gamma_bounds[0]<=g<=self.gamma_bounds[1]):
            raise ValueError('field outside the registered local energy patch; no extrapolation')
        self.vacuum.evaluate(a,g)
        x=2*(a-self.amplitude_bounds[0])/np.diff(self.amplitude_bounds)[0]-1
        y=2*(g-self.gamma_bounds[0])/np.diff(self.gamma_bounds)[0]-1
        polynomials=[chebval2d(x,y,c) for c in (self.coefficients,self._da,self._dg)]
        with np.errstate(over='raise',invalid='raise',under='ignore'):
            increments=np.exp(polynomials[0])
            values=(np.cumsum(increments),np.cumsum(increments*polynomials[1]),
                    np.cumsum(increments*polynomials[2]))
        if (np.any(~np.isfinite(values)) or np.any(increments<=0)
                or np.any(np.diff(values[0])<=0)):
            raise FloatingPointError('positive increments lost resolved finite ordering; no repair')
        return values

    def save(self,path):
        with Path(path).open('xb') as stream:
            np.savez_compressed(stream,coefficients=self.coefficients,count_nodes=self.count_nodes,
                count_weights=self.count_weights,amplitude_bounds=self.amplitude_bounds,
                gamma_bounds=self.gamma_bounds,eta=self.eta,representation=self.representation,
                scales=[self.vacuum.delta0_J,self.vacuum.N0_per_J_m3,self.vacuum.D_m2_s])

    @classmethod
    def load(cls,path,source):
        """Caller verifies source/artifact hashes; this adds representation guards."""
        with np.load(path,allow_pickle=False) as data:
            if 'representation' not in data or data['representation'].item()!=cls.representation:
                raise ValueError('wrong energy representation')
        return super().load(path,source)
