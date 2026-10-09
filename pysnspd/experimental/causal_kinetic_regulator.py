"""Causal R/A/K completion of the SAME numerical spectral eta.

This regulator is algebraic and follows H, rather than a physical reservoir.
Its Keldysh entry is SigmaR star H - H star SigmaA. It restores consistency
between spectra solved with B+eta*tau3 and the kinetic equation. In particular
its two physical trace projections vanish for diagonal H and constant eta;
its full anomalous matrix entries are never labeled heating.
The physical B and collision source must remain explicitly independent.
"""
from __future__ import annotations
from dataclasses import replace
import warnings
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import MatrixRankWarning,spsolve
from . import adiabatic_matrix_usadel as ad
from . import matrix_kinetic_rhs as kinetic

CONTRACT='Numerical causal completion SigmaR=eta*tau3, SigmaA=-eta*tau3, SigmaK=SigmaR star H-H star SigmaA; follows the same H; no physical collision bath or heating'


def _eta(eta):
    raw=np.asarray(eta)
    if raw.ndim or np.iscomplexobj(raw) or raw.dtype.kind=='b' or not np.isfinite(raw) or raw<0:
        raise ValueError('Numerical eta must be a finite nonnegative real scalar')
    return float(raw)


def self_energy_triad(H,eta,*,kappa_s):
    eta=_eta(eta);kappa=kinetic._positive(kappa_s,'kappa_s')
    SR=ad.Jet.of(np.broadcast_to(eta*ad.TAU3,H.value.shape));SA=SR.scale(-1)
    SK=ad.star(SR,H,kappa)-ad.star(H,SA,kappa)
    return SR,SA,SK


def regulator_term(R,H,eta,*,kappa_s):
    """Full matrix local term, without node mass: not a physical source.

    The exact first-gradient associative reduction is
    [SigmaR,R]_star star H - H star [SigmaA,A]_star.
    No full residual is projected away; only the known causal eta is completed.
    """
    SR,SA,_=self_energy_triad(H,eta,kappa_s=kappa_s)
    return (ad.star(ad.commutator(SR,R,kappa_s),H,kappa_s)
            -ad.star(H,ad.commutator(SA,ad.advanced(R),kappa_s),kappa_s))


def kinetic_residual(graph,R,H,B_without_eta,*,numerical_eta,kappa_s,alpha=None,alpha_time_per_s=None):
    raw,flux,K=ad.kinetic_residual(graph,R,H,B_without_eta,alpha=alpha,
        alpha_time=alpha_time_per_s,kappa=kappa_s)
    correction=regulator_term(R,H,numerical_eta,kappa_s=kappa_s).scale(graph.area_weights/2)
    return raw-correction,flux,K,correction


def assemble_storage(graph,R,B_without_eta,*,numerical_eta,kappa_s,alpha=None):
    M=kinetic.assemble_storage(graph,R,B_without_eta,alpha=alpha,kappa_s=kappa_s)
    R=kinetic._stationary_jet(R);n=graph.n_nodes
    blocks=np.empty((n,2,2))
    for component,basis in enumerate((ad.I2,ad.TAU3)):
        H=ad.Jet.of(np.zeros_like(R.value),time=np.broadcast_to(basis,R.value.shape))
        correction=regulator_term(R,H,numerical_eta,kappa_s=kappa_s).scale(-graph.area_weights/2)
        blocks[:,:,component]=ad.projections(correction.first)
    nodes=np.arange(n);rows=[];cols=[];values=[]
    for out in range(2):
        for component in range(2):
            rows.append(2*nodes+out);cols.append(2*nodes+component);values.append(blocks[:,out,component])
    extra=coo_matrix((np.concatenate(values),(np.concatenate(rows),np.concatenate(cols))),shape=M.shape).tocsr()
    answer=M+extra;answer.eliminate_zeros();return answer


def solve_distribution_rhs(graph,R,h,h_energy,B_without_eta,*,numerical_eta,kappa_s,
        projected_source,source_contract,alpha=None,alpha_time_per_s=None,
        fixed_nodes=(),fixed_rhs_per_s=None,require_pauli_support=True):
    """Physical source is explicit; causal eta completion is retained, not heat.

    R/RE/Rt/R1 MUST come from spectral B_without_eta+numerical_eta*tau3.
    This consistency is a caller contract; the full Nambu residual is retained
    so spectral root/tangent errors are visible. Eta-backed subgap storage is
    not by itself physical population support: eta convergence remains required.
    """
    eta=_eta(numerical_eta);kappa=kinetic._positive(kappa_s,'kappa_s');n=graph.n_nodes
    h=kinetic._real_fields(h,n,'h');hE=kinetic._real_fields(h_energy,n,'h_energy')
    source=kinetic._real_fields(projected_source,n,'projected_source')
    if not isinstance(source_contract,str) or not source_contract.strip():raise ValueError('Explicit physical source/self-energy contract required')
    if np.any(abs((B_without_eta.value[:,0,0]-B_without_eta.value[:,1,1]).real)>1e-12):
        raise ValueError('Numerical eta is forbidden in physical kinetic B')
    margin=1-np.max(abs(np.column_stack((h[:,0]+h[:,1],h[:,0]-h[:,1]))))
    if require_pauli_support and margin<-32*np.finfo(float).eps:raise ValueError('hL +/- hT leave Pauli support; no clipping applied')
    M=assemble_storage(graph,R,B_without_eta,numerical_eta=eta,kappa_s=kappa,alpha=alpha)
    H0=ad.distribution(h[:,0],h[:,1],energy_L=hE[:,0],energy_T=hE[:,1])
    base,_,_,_=kinetic_residual(graph,R,H0,B_without_eta,numerical_eta=eta,kappa_s=kappa,alpha=alpha,alpha_time_per_s=alpha_time_per_s)
    frozen=ad.projections(base.value);temporal=ad.projections(base.first);rhs=(source-frozen-temporal).ravel()
    fixed=np.asarray(fixed_nodes)
    if fixed.size==0:fixed=np.empty(0,int)
    if (fixed.dtype.kind not in 'iu' or fixed.ndim!=1 or np.any(fixed<0) or np.any(fixed>=n) or len(np.unique(fixed))!=len(fixed)):
        raise ValueError('Unique integer prescribed nodes required')
    value=np.zeros((n,2))
    if len(fixed):
        if fixed_rhs_per_s is None:raise ValueError('Explicit time law of prescribed distributions required')
        value[fixed]=kinetic._real_fields(fixed_rhs_per_s,len(fixed),'fixed_rhs_per_s')
    elif fixed_rhs_per_s is not None and np.asarray(fixed_rhs_per_s).size:raise ValueError('Boundary RHS without prescribed nodes')
    free=np.ones(n,bool);free[fixed]=False;components=np.flatnonzero(np.repeat(free,2));vector=value.ravel()
    if len(components):
        with warnings.catch_warnings():
            warnings.simplefilter('error',MatrixRankWarning)
            try:vector[components]=spsolve(M[components][:,components],(rhs-M@vector)[components])
            except MatrixRankWarning as error:raise FloatingPointError('Singular temporal storage; no DOS floor or regularization') from error
    if np.any(~np.isfinite(vector)):raise FloatingPointError('Nonfinite kinetic RHS; no repair applied')
    H=ad.distribution(h[:,0],h[:,1],energy_L=hE[:,0],energy_T=hE[:,1],time_L=value[:,0],time_T=value[:,1])
    residual,flux,_,_=kinetic_residual(graph,R,H,B_without_eta,numerical_eta=eta,kappa_s=kappa,alpha=alpha,alpha_time_per_s=alpha_time_per_s)
    projected=ad.projections(residual.value+residual.first)-source;error=M@vector-rhs
    return kinetic.MatrixKineticRHS(value.copy(),M,frozen,temporal,source.copy(),projected,projected[fixed].copy(),
        ad.projections(flux.value),ad.projections(flux.first),residual.value+residual.first,
        float(np.max(abs(projected[free]),initial=0)),float(np.max(abs(error[components]),initial=0)),
        float(margin),kappa,source_contract+'; '+CONTRACT,
        'Causally completed numerical eta on supplied first-gradient spectrum; no detector admission')
