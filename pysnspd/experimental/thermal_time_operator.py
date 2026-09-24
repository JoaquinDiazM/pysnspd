"""Reuse saved thermal spectral Jacobians at alpha=0 for full-node directions.

The admitted thermal core has no Peierls link field. This helper does not infer
or silently discard a nonzero vector potential from a different experiment.
"""
import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import splu


class SavedSpectralMode:
    def __init__(self,graph,arrays,epsilon,temperature_ratio):
        self.graph=graph;self.epsilon=float(epsilon);self.weight=2*np.pi*temperature_ratio
        self.u=np.asarray(arrays['u'],complex).copy();self.f=np.asarray(arrays['f'],complex).copy();self.g=np.asarray(arrays['g'],float).copy()
        self.components=np.asarray(arrays['free_real_components'],int).copy()
        matrix=csc_matrix((arrays['jacobian_data'],arrays['jacobian_indices'],arrays['jacobian_indptr']),shape=tuple(arrays['jacobian_shape']))
        self.matrix=matrix;self.factor=splu(matrix)

    def apply(self,direction):
        v=np.asarray(direction,complex)
        if v.shape!=(self.graph.n_nodes,) or np.any(~np.isfinite(v)):raise ValueError('Finite full-node direction required')
        allowed=np.zeros(2*len(v),bool);allowed[self.components]=True
        cart=np.column_stack((v.real,v.imag)).ravel()
        if np.any(cart[~allowed]!=0):raise ValueError('Saved fixed spectral contacts cannot move')
        rhs=np.repeat(self.graph.area_weights,2)*cart
        ducart=np.zeros(2*len(v));ducart[self.components]=self.factor.solve(rhs[self.components])
        du=ducart[::2]+1j*ducart[1::2]
        df=self.g*du-self.g**3*self.u*np.real(np.conj(self.u)*du)
        tail,head=self.graph.edges.T
        current=self.weight*2*self.graph.conductance*np.imag(np.conj(df[tail])*self.f[head]+np.conj(self.f[tail])*df[head])
        hessian=2*self.weight*self.graph.area_weights*(v/self.epsilon-df)
        residual=float(np.max(abs(self.matrix@ducart[self.components]-rhs[self.components])))
        return hessian,current,residual


class FullNodeCoordinates:
    """Mass-weighted Cartesian coordinates; every unconstrained node is retained."""
    def __init__(self,graph,gap_scale):
        self.graph=graph;self.free=np.ones(graph.n_nodes,bool);self.free[graph.boundary_nodes]=False
        self.scale=np.sqrt(graph.area_weights[self.free])/gap_scale
        self.size=2*np.count_nonzero(self.free)

    def encode(self,field):
        z=np.asarray(field,complex)
        if z.shape!=(self.graph.n_nodes,) or np.any(z[~self.free]!=0):raise ValueError('Full-node field with zero fixed contacts required')
        value=z[self.free]*self.scale
        return np.column_stack((value.real,value.imag)).ravel()

    def decode(self,value):
        raw=np.asarray(value,float)
        if raw.shape!=(self.size,):raise ValueError('Incorrect full-node coordinate dimension')
        z=np.zeros(self.graph.n_nodes,complex);z[self.free]=(raw[::2]+1j*raw[1::2])/self.scale
        return z
