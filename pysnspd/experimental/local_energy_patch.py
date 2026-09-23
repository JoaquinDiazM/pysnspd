"""Conservative local interpolation of E(x,a,Gamma) for weak spatial trials.

Both field derivatives come from the SAME tensor Chebyshev polynomial as the
energy. The native state-count quadrature and exact vacuum are retained. This
is a bounded numerical approximation, not a new material or a global catalogue.
Queries outside its declared field box fail; no extrapolation or clipping.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.polynomial.chebyshev import chebder, chebval2d, chebvander2d


@dataclass(frozen=True)
class LocalEnergyPatch:
    source: object
    amplitude_bounds: tuple[float, float]
    gamma_bounds: tuple[float, float]
    coefficients: np.ndarray

    def __post_init__(self):
        for name in ('amplitude_bounds','gamma_bounds'):
            raw=np.asarray(getattr(self,name))
            if np.iscomplexobj(raw):raise ValueError('field bounds must be real')
            bounds=np.asarray(raw,float)
            if bounds.shape!=(2,) or np.any(~np.isfinite(bounds)) or bounds[0]>=bounds[1]:
                raise ValueError('two finite increasing bounds required')
            object.__setattr__(self,name,tuple(map(float,bounds)))
        if np.iscomplexobj(self.coefficients):raise ValueError('coefficients must be real')
        coeff=np.array(self.coefficients,float,copy=True)
        if (coeff.ndim!=3 or min(coeff.shape[:2])<3 or coeff.shape[2]!=len(self.source.count_nodes)
                or np.any(~np.isfinite(coeff))):
            raise ValueError('finite tensor coefficients and matching state count required')
        coeff.setflags(write=False);object.__setattr__(self,'coefficients',coeff)
        da=chebder(coeff,axis=0)*(2/np.diff(self.amplitude_bounds)[0])
        dg=chebder(coeff,axis=1)*(2/np.diff(self.gamma_bounds)[0])
        da.setflags(write=False);dg.setflags(write=False)
        object.__setattr__(self,'_da',da);object.__setattr__(self,'_dg',dg)
        for a in self.amplitude_bounds:
            for g in self.gamma_bounds:self.vacuum.evaluate(a,g)

    @property
    def count_nodes(self):return self.source.count_nodes

    @property
    def count_weights(self):return self.source.count_weights

    @property
    def vacuum(self):return self.source.vacuum

    @property
    def eta(self):return self.source.eta

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
        # Assign exact endpoints to avoid a construction-only roundoff overshoot.
        aa=np.where(xx.ravel()==1,bounds[0,1],np.where(xx.ravel()==-1,bounds[0,0],aa))
        gg=np.where(yy.ravel()==1,bounds[1,1],np.where(yy.ravel()==-1,bounds[1,0],gg))
        values=[]
        for i,(a,g) in enumerate(zip(aa,gg)):
            values.append(source.energy_kernel(float(a),float(g))[0])
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
        values=tuple(chebval2d(x,y,c) for c in (self.coefficients,self._da,self._dg))
        if (np.any(~np.isfinite(values)) or np.any(values[0]<=0)
                or np.any(np.diff(values[0])<=0)):
            raise FloatingPointError('local interpolation lost finite monotone positive energies')
        return values

    def evaluate(self,amplitude,gamma,occupation):
        raw=np.asarray(occupation)
        if np.iscomplexobj(raw):raise ValueError('occupation must be real')
        p=np.asarray(raw,float)
        if p.shape!=self.count_nodes.shape or np.any(~np.isfinite(p)) or np.any((p<0)|(p>1)):
            raise ValueError('occupation must match count nodes and Pauli bounds')
        return tuple(float(u+4*np.dot(self.count_weights*k,p)) for u,k in
                     zip(self.vacuum.evaluate(amplitude,gamma),self.energy_kernel(amplitude,gamma)))

    def save(self,path):
        path=Path(path)
        with path.open('xb') as stream:
            np.savez_compressed(stream,coefficients=self.coefficients,count_nodes=self.count_nodes,
                count_weights=self.count_weights,amplitude_bounds=self.amplitude_bounds,
                gamma_bounds=self.gamma_bounds,eta=self.eta,
                scales=[self.vacuum.delta0_J,self.vacuum.N0_per_J_m3,self.vacuum.D_m2_s])

    @classmethod
    def load(cls,path,source):
        """Restore numerical data; caller must verify source/artifact SHA manifest.

        Array, regulator and scale compatibility checks below do not identify
        the source physics implementation. The registered runner verifies its
        frozen source hashes and the patch-file hash before calling this method.
        """
        with np.load(path,allow_pickle=False) as data:
            if (not np.array_equal(data['count_nodes'],source.count_nodes)
                or not np.array_equal(data['count_weights'],source.count_weights)
                or data['eta'].item()!=source.eta
                or not np.array_equal(data['scales'],[source.vacuum.delta0_J,source.vacuum.N0_per_J_m3,source.vacuum.D_m2_s])):
                raise ValueError('patch source quadrature, regulator or physical scales mismatch')
            return cls(source,tuple(data['amplitude_bounds']),tuple(data['gamma_bounds']),data['coefficients'])
