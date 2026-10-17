"""Generic local graph gradient G_i(rho_i,L_i,T_i), SI and analytic geometry.

Variational graph geometry shared by the thermal-gradient closures.
Wirtinger force: delta G = 2 Re<F,delta d> + G_alpha delta alpha + G_T delta T.
L_i = .5 sum_{edges incident i} c_e |d_t-exp(-i alpha_e)d_h|^2.
"""
from dataclasses import dataclass
import numpy as np
from scipy.constants import elementary_charge as electron_charge,hbar
from scipy.sparse import coo_matrix,diags
from scipy.special import expit
from pysnspd.experimental.effective_bcs_functional import EffectiveBCSFunctional,EffectiveBCSDirection
from pysnspd.experimental.thermal_stiffness import ThermalStiffnessFunctional,ThermalStiffnessEvaluation

@dataclass(frozen=True)
class GradientPartials:
    G_J:np.ndarray
    G_rho_J:np.ndarray
    G_L:np.ndarray
    G_rho_rho_J:np.ndarray
    G_rho_L:np.ndarray
    G_L_L_per_J:np.ndarray
    G_T_J_K:np.ndarray
    G_rho_T_J_K:np.ndarray
    G_L_T_per_K:np.ndarray
    G_T_T_J_K2:np.ndarray

    def matching(self,n):
        values={}
        for key,value in self.__dict__.items():
            a=np.broadcast_to(np.asarray(value,float),(n,))
            if np.any(~np.isfinite(a)):raise ValueError("Finite matching gradient partials required: "+key)
            values[key]=a
        return GradientPartials(**values)

@dataclass(frozen=True)
class GradientGeometry:
    rho:np.ndarray
    L_J:np.ndarray
    phase:np.ndarray
    difference:np.ndarray
    product:np.ndarray
    edges:np.ndarray
    c_J:np.ndarray

def gradient_geometry(d,edges,c_J,alpha=None):
    d=np.asarray(d,complex);edges=np.asarray(edges,int);c=np.asarray(c_J,float)
    alpha=np.zeros(len(edges)) if alpha is None else np.asarray(alpha,float)
    if d.ndim!=1 or edges.ndim!=2 or edges.shape[1]!=2 or c.shape!=(len(edges),) or alpha.shape!=(len(edges),):
        raise ValueError("Matching graph fields required")
    if np.any(~np.isfinite(d)) or np.any(~np.isfinite(c)) or np.any(c<0) or np.any(~np.isfinite(alpha)):
        raise ValueError("Finite graph fields and nonnegative edge coefficients required")
    tail,head=edges.T;phase=np.exp(-1j*alpha)
    difference=d[tail]-phase*d[head];product=np.conj(d[tail])*phase*d[head]
    edge=.5*c*abs(difference)**2
    L=np.bincount(tail,edge,minlength=len(d))+np.bincount(head,edge,minlength=len(d))
    return GradientGeometry(abs(d)**2,L,phase,difference,product,edges,c)

@dataclass(frozen=True)
class GraphGradientEvaluation:
    geometry:GradientGeometry
    partials:GradientPartials
    free_energy_J:float
    internal_energy_J:float
    entropy_J_K:np.ndarray
    node_internal_energy_J:np.ndarray
    force_J:np.ndarray
    edge_current_A:np.ndarray
    link_derivative_alpha_J:np.ndarray

def gradient_evaluate(d,T,geometry,partials):
    d=np.asarray(d,complex);T=np.asarray(T,float);n=len(d)
    if T.shape!=(n,) or np.any(T<=0) or np.any(~np.isfinite(T)):raise ValueError("Finite positive nodal T required")
    g=partials.matching(n);tail,head=geometry.edges.T;c=geometry.c_J
    ke=(g.G_L[tail]+g.G_L[head])/2
    force=g.G_rho_J*d
    np.add.at(force,tail,c*ke*geometry.difference)
    np.add.at(force,head,-c*ke*np.conj(geometry.phase)*geometry.difference)
    link=-2*c*ke*np.imag(geometry.product)
    current=-2*electron_charge/hbar*link
    internal=g.G_J-T*g.G_T_J_K
    return GraphGradientEvaluation(geometry,g,float(np.sum(g.G_J)),float(np.sum(internal)),
        -g.G_T_J_K,internal,force,current,link)

def gradient_geometry_jacobians(d,geometry):
    """JR=d rho/d cart, JL=d L/d cart, JS=d Im(conj dt phase dh)/d cart."""
    d=np.asarray(d,complex);n=len(d);m=len(geometry.edges)
    tail,head=geometry.edges.T;phase=geometry.phase;c=geometry.c_J
    cols=np.column_stack((2*tail,2*tail+1,2*head,2*head+1))
    w=geometry.difference;headw=-np.conj(phase)*w
    lvalues=c[:,None]*np.column_stack((w.real,w.imag,headw.real,headw.imag))
    rows=np.vstack((np.repeat(tail[:,None],4,axis=1),np.repeat(head[:,None],4,axis=1)))
    JL=coo_matrix((np.tile(lvalues.ravel(),2),(rows.ravel(),np.tile(cols.ravel(),2))),shape=(n,2*n)).tocsr()
    JR=coo_matrix((np.column_stack((2*d.real,2*d.imag)).ravel(),
        (np.repeat(np.arange(n),2),np.arange(2*n))),shape=(n,2*n)).tocsr()
    rhead=phase*d[head];rtail=np.conj(phase)*d[tail]
    # Head gradient is R(alpha)^T[-Im(dt),Re(dt)].
    svalues=np.column_stack((rhead.imag,-rhead.real,-rtail.imag,rtail.real))
    JS=coo_matrix((svalues.ravel(),(np.repeat(np.arange(m),4),cols.ravel())),shape=(m,2*n)).tocsr()
    return JR,JL,JS

def weighted_phase_laplacian(geometry,weights,n):
    tail,head=geometry.edges.T;phase=geometry.phase;w=np.asarray(weights,float)
    rot=np.stack((np.stack((phase.real,-phase.imag),axis=1),
                  np.stack((phase.imag,phase.real),axis=1)),axis=1)
    blocks=(np.broadcast_to(np.eye(2),(len(w),2,2))*w[:,None,None],
            -w[:,None,None]*rot,-w[:,None,None]*np.swapaxes(rot,1,2),
            np.broadcast_to(np.eye(2),(len(w),2,2))*w[:,None,None])
    rownodes=(tail,tail,head,head);colnodes=(tail,head,tail,head)
    rows=np.concatenate([np.repeat(2*rr,4)+np.tile([0,0,1,1],len(w)) for rr in rownodes])
    cols=np.concatenate([np.repeat(2*cc,4)+np.tile([0,1,0,1],len(w)) for cc in colnodes])
    vals=np.concatenate([block.ravel() for block in blocks])
    return coo_matrix((vals,(rows,cols)),shape=(2*n,2*n)).tocsr()

def gradient_hessian_current_jacobian(d,evaluation):
    """Fixed-T analytic half real Hessian and current Jacobian."""
    geom=evaluation.geometry;g=evaluation.partials;n=len(d);tail,head=geom.edges.T
    JR,JL,JS=gradient_geometry_jacobians(d,geom)
    H=.5*(JR.T@diags(g.G_rho_rho_J)@JR+
       JR.T@diags(g.G_rho_L)@JL+JL.T@diags(g.G_rho_L)@JR+
       JL.T@diags(g.G_L_L_per_J)@JL)
    ke=(g.G_L[tail]+g.G_L[head])/2
    H=H+diags(np.repeat(g.G_rho_J,2))+weighted_phase_laplacian(geom,geom.c_J*ke,n)
    average=coo_matrix((np.full(2*len(tail),.5),
       (np.repeat(np.arange(len(tail)),2),np.column_stack((tail,head)).ravel())),
       shape=(len(tail),n)).tocsr()
    dGL=diags(g.G_rho_L)@JR+diags(g.G_L_L_per_J)@JL
    JI=(4*electron_charge/hbar*diags(geom.c_J)@
        (diags(np.imag(geom.product))@average@dGL+diags(ke)@JS))
    return H.tocsr(),JI.tocsr()

def gradient_alpha_jacobians(d,evaluation):
    """Fixed-T d(force_cart)/dalpha and d(current)/dalpha, analytically.

    The mixed force Jacobian is half the real G_{cart,alpha} Hessian.
    Current-alpha Jacobian equals -(2e/hbar)*G_{alpha,alpha}.
    """
    geom=evaluation.geometry;g=evaluation.partials;n=len(d);tail,head=geom.edges.T
    JR,JL,JS=gradient_geometry_jacobians(d,geom);m=len(tail)
    JLa=coo_matrix((-np.repeat(geom.c_J*np.imag(geom.product),2),
        (np.column_stack((tail,head)).ravel(),np.repeat(np.arange(m),2))),shape=(n,m)).tocsr()
    ke=(g.G_L[tail]+g.G_L[head])/2
    force_alpha=.5*(JR.T@diags(g.G_rho_L)@JLa+
        JL.T@diags(g.G_L_L_per_J)@JLa)-JS.T@diags(geom.c_J*ke)
    current_alpha=(-2*electron_charge/hbar)*(JLa.T@diags(g.G_L_L_per_J)@JLa+
        diags(2*geom.c_J*ke*np.real(geom.product)))
    return force_alpha.tocsr(),current_alpha.tocsr()

@dataclass(frozen=True)
class GraphGradientDirection:
    free_energy_direction_J:float
    internal_energy_direction_J:float
    G_T_direction_J_K:np.ndarray
    rho_direction:np.ndarray
    L_direction_J:np.ndarray
    force_direction_J:np.ndarray
    current_direction_A:np.ndarray

def gradient_direction(d,T,dd,dT,evaluation,alpha_direction=None):
    d=np.asarray(d,complex);T=np.asarray(T,float);dd=np.asarray(dd,complex);dT=np.asarray(dT,float)
    geom=evaluation.geometry;g=evaluation.partials;tail,head=geom.edges.T;c=geom.c_J
    da=np.zeros(len(c)) if alpha_direction is None else np.asarray(alpha_direction,float)
    if dd.shape!=d.shape or dT.shape!=T.shape or da.shape!=(len(c),) or np.any(~np.isfinite(dd)) or np.any(~np.isfinite(dT)) or np.any(~np.isfinite(da)):
        raise ValueError("Matching finite gradient direction required")
    drho=2*np.real(np.conj(d)*dd)
    dw=dd[tail]-geom.phase*dd[head]+1j*da*geom.phase*d[head]
    edge_Ldot=c*np.real(np.conj(geom.difference)*dw)
    dL=np.bincount(tail,edge_Ldot,minlength=len(d))+np.bincount(head,edge_Ldot,minlength=len(d))
    dGr=g.G_rho_rho_J*drho+g.G_rho_L*dL+g.G_rho_T_J_K*dT
    dGL=g.G_rho_L*drho+g.G_L_L_per_J*dL+g.G_L_T_per_K*dT
    ke=(g.G_L[tail]+g.G_L[head])/2;dke=(dGL[tail]+dGL[head])/2
    df=dGr*d+g.G_rho_J*dd
    common=c*(dke*geom.difference+ke*dw)
    np.add.at(df,tail,common)
    np.add.at(df,head,-np.conj(geom.phase)*(common+1j*da*c*ke*geom.difference))
    dproduct=(np.conj(dd[tail])*geom.phase*d[head]+np.conj(d[tail])*geom.phase*dd[head]-1j*da*geom.product)
    di=4*electron_charge/hbar*c*(dke*np.imag(geom.product)+ke*np.imag(dproduct))
    dGT=g.G_rho_T_J_K*drho+g.G_L_T_per_K*dL+g.G_T_T_J_K2*dT
    fixed=float(np.dot(g.G_rho_J,drho)+np.dot(g.G_L,dL))
    dG=fixed+float(np.dot(g.G_T_J_K,dT))
    dU=fixed-float(np.dot(T,dGT))
    return GraphGradientDirection(dG,dU,dGT,drho,dL,df,di)

class GeneralThermalGradientFunctional(ThermalStiffnessFunctional):
    """Native local BCS plus arbitrary differentiable G(rho,L,T).

    gradient_partials(rho,L_J,T_K) must return GradientPartials.
    Caller supplies the scientific constitutive law; this class only chains it.
    """
    def __init__(self,*args,gradient_partials,closure_id="general_local_gradient_candidate",**kwargs):
        super().__init__(*args,**kwargs)
        self.gradient_partials=gradient_partials
        self.closure_id=closure_id

    def _local_bcs(self):
        obj=object.__new__(EffectiveBCSFunctional);obj.__dict__.update(self.__dict__)
        obj.edge_energy_scale_J=np.zeros_like(self.edge_energy_scale_J)
        return obj

    def gradient_evaluation(self,d,T,alpha=None):
        geometry=gradient_geometry(d,self.graph.edges,self.edge_energy_scale_J,alpha)
        data=self.gradient_partials(geometry.rho,geometry.L_J,T)
        return gradient_evaluate(d,T,geometry,data)

    def thermal_data(self,d,p,temp=None,alpha=None):
        d,p,alpha=self._state(d,p,alpha);T=self._temperature(d,p,temp)
        # Native differential moments are independent of the gradient law.
        en=self.energies_bar(d);theta=T/self.Tc_K;u=self.node_energy_scale_J;w=self.weights_bar
        fd=expit(-en/theta[:,None])
        C=4*u*np.sum(w*en**2*fd*(1-fd),axis=1)/(self.Tc_K*theta**2)
        B=2*u*np.sum(w*((p-fd)/en+fd*(1-fd)/theta[:,None]),axis=1)
        td=dict(capacity_FD_J_K=C,B_J=B,epsilon_J=4*u[:,None]*w*en,
            temperature_K=T,fd=fd)
        ge=self.gradient_evaluation(d,T,alpha);g=ge.partials
        td.update(L_J=ge.geometry.L_J,gradient_partials=g,
            G_T_J_per_K=g.G_T_J_K,G_TT_J_per_K2=g.G_T_T_J_K2,
            G_Trho_J_per_K=g.G_rho_T_J_K,G_TL_per_K=g.G_L_T_per_K,
            capacity_effective_J_K=td['capacity_FD_J_K']-T*g.G_T_T_J_K2)
        return td

    def evaluate(self,gapbar,p,alpha=None,*,temperature_K=None):
        d,p,alpha=self._state(gapbar,p,alpha);T=self._temperature(d,p,temperature_K)
        local=self._local_bcs().evaluate(d,p,alpha)
        ge=self.gradient_evaluation(d,T,alpha);g=ge.partials
        force=local.gap_force_J_per_gapbar+ge.force_J
        node=local.node_energy_J+ge.node_internal_energy_J
        div=np.zeros(len(d));tail,head=self.graph.edges.T
        np.add.at(div,tail,ge.edge_current_A);np.add.at(div,head,-ge.edge_current_A)
        return ThermalStiffnessEvaluation(float(node.sum()),local.local_energy_J,ge.internal_energy_J,
            node,local.electron_energies_J,force,2*force,ge.edge_current_A,div,
            ge.link_derivative_alpha_J,local.local_coefficient_bar,g.G_L,
            ge.free_energy_J,ge.entropy_J_K,g.G_T_J_K,g.G_T_T_J_K2,ge.geometry.L_J,T,
            "G(rho,L,Te) free gradient; Ugrad=G-Te*G_T; F=d(U_BCS+G)/dd* at fixed Te")

    def hessian_and_current_jacobian(self,gapbar,p,alpha=None,*,temperature_K=None):
        d,p,alpha=self._state(gapbar,p,alpha);T=self._temperature(d,p,temperature_K)
        local,_=self._local_bcs().hessian_and_current_jacobian(d,p,alpha)
        H,JI=gradient_hessian_current_jacobian(d,self.gradient_evaluation(d,T,alpha))
        return local+H,JI

    def direction(self,gapbar,p,gap_direction,p_direction=None,alpha=None,alpha_direction=None,
                  *,temperature_K=None,temperature_direction_K=None):
        d,p,alpha=self._state(gapbar,p,alpha);T=self._temperature(d,p,temperature_K)
        dd=np.asarray(gap_direction,complex);dp=np.zeros_like(p) if p_direction is None else np.asarray(p_direction,float)
        da=np.zeros_like(alpha) if alpha_direction is None else np.asarray(alpha_direction,float)
        td=self.thermal_data(d,p,T,alpha);drho=2*np.real(np.conj(d)*dd)
        dT=((td['B_J']*drho+np.sum(td['epsilon_J']*dp,axis=1))/td['capacity_FD_J_K']
             if temperature_direction_K is None else np.broadcast_to(np.asarray(temperature_direction_K,float),T.shape))
        local=self._local_bcs().direction(d,p,dd,dp,alpha,da)
        ge=self.gradient_evaluation(d,T,alpha);dg=gradient_direction(d,T,dd,dT,ge,da)
        return EffectiveBCSDirection(local.energy_direction_J+dg.internal_energy_direction_J,
            local.gap_force_direction_J_per_gapbar+dg.force_direction_J,
            dg.current_direction_A,local.electron_energy_direction_J)
