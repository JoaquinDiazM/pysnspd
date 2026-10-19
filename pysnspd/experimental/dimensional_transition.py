"""Conservative, guarded 2D-to-1D restriction of the neutral BCS/KWT model.

No hotbelt constitutive law is added. The continuum restriction is dy=0,
dV=W*df*dx, face area=W*df; the native energy grid and circuit are unchanged.
The finite-volume remap and its small energy reconciliation are explicit
approximation operations, accepted only under recorded tolerances.
"""
from dataclasses import dataclass
import numpy as np
from scipy.optimize import brentq
from scipy.special import expit
from .thermal_spatial_usadel import ThermalGraph
from .effective_bcs_preparation import EffectiveBCSGeometry
from .effective_coupled import EffectiveState
from .thermal_coupled import ThermalCoupledSystem


@dataclass(frozen=True)
class TransitionTolerances:
    gap_relative: float = .03
    populations_energy_relative: float = .05
    current_pattern_relative: float = .03
    projection_energy_relative: float = .005
    reconciliation_energy_relative: float = .01
    shadow_voltage_V: float = 2e-6
    shadow_current_relative: float = .005
    shadow_energy_relative: float = .01
    cut_energy_relative: float = .02


def chain_geometry(length_m, width_m, spacing_m, ell0_m):
    count = max(3, int(round(length_m/spacing_m))+1)
    x = np.linspace(0., length_m, count)
    dx = np.diff(x)
    dual = np.r_[dx[0]/2, (dx[:-1]+dx[1:])/2, dx[-1]/2]
    edges = np.column_stack((np.arange(count-1), np.arange(1,count)))
    xy = np.column_stack((x,np.zeros(count)))
    fixed = np.array([0,count-1])
    graph = ThermalGraph(width_m*dual/ell0_m**2, edges, width_m/dx,
                         xy/ell0_m, fixed)
    return EffectiveBCSGeometry(graph, xy, np.empty((0,3),int), dx,
        np.full(count-1,width_m), fixed[:1], fixed[1:], ell0_m,
        length_m,width_m,"generated:1D-finite-volume",
        dict(dimension=1,length_m=length_m,width_m=width_m,
             nodes=count,spacing_m=float(dx[0]),volume_rule="W*df*dual_dx"))


def volume_restriction(source_x, source_volume, target_volume):
    """Positive monotone FV remap, exact constant/volume marginals.

    Cumulative-volume bins are ordered by longitudinal coordinate. This is a
    declared remap, not an exact transverse integral on an irregular mesh.
    """
    v = np.asarray(source_volume,float); w = np.asarray(target_volume,float)
    if np.any(v<=0) or np.any(w<=0) or not np.isclose(v.sum(),w.sum(),rtol=2e-10,atol=0):
        raise ValueError("Incompatible positive physical volumes")
    order = np.argsort(source_x,kind="stable")
    a = np.r_[0.,np.cumsum(v[order])]
    b = np.r_[0.,np.cumsum(w)]
    # Remove accumulated roundoff at the identical physical right cut.
    a[-1]=b[-1]
    overlap=np.maximum(0.,np.minimum(b[1:,None],a[None,1:])-
                       np.maximum(b[:-1,None],a[None,:-1]))
    result=np.zeros_like(overlap);result[:,order]=overlap/w[:,None]
    if not np.allclose(result.sum(axis=1),1.,rtol=1e-10,atol=1e-12):
        raise ArithmeticError("Restriction does not reproduce constants")
    if not np.allclose(w@result,v,rtol=1e-9,atol=1e-40):
        raise ArithmeticError("Restriction is not volume conservative")
    return result


def film_energy(system,state):
    f,k=system.functional,system.kinetics
    temp=k.equivalent_temperature_K(state.delta_bar,state.p)
    return float(f.evaluate(state.delta_bar,state.p,temperature_K=temp).energy_J+
                 np.sum(k.phonon_energy_weights_J*state.n))


def population_distance(system,a,b,scale):
    f,k=system.functional,system.kinetics
    eps=4*f.node_energy_scale_J[:,None]*f.weights_bar*f.energies_bar(a.delta_bar)
    return float((np.sum(eps*abs(a.p-b.p))+
                  np.sum(k.phonon_energy_weights_J*abs(a.n-b.n)))/scale)


def resolved_winding(delta,triangles,dc_gap,floor=1e-3):
    if len(triangles)==0:return 0
    z=delta[triangles]
    phase=np.angle(np.roll(z,-1,axis=1)*np.conj(z)).sum(axis=1)
    good=np.min(abs(z),axis=1)>floor*float(np.mean(abs(dc_gap)))
    return int(np.count_nonzero(good & (abs(phase)>np.pi)))


class DimensionalProjection:
    def __init__(self,system2,geometry2,dc2,system1,geometry1,dc1):
        self.s2,self.g2,self.dc2=system2,geometry2,dc2
        self.s1,self.g1,self.dc1=system1,geometry1,dc1
        self.R=volume_restriction(geometry2.coordinates_m[:,0],
            system2.functional.volumes_m3,system1.functional.volumes_m3)
        self.Edc2=film_energy(system2,dc2);self.Edc1=film_energy(system1,dc1)

    def lift(self,state):
        x=self.g2.coordinates_m[:,0];xx=self.g1.coordinates_m[:,0]
        ratio=state.delta_bar/self.dc1.delta_bar
        d=self.dc2.delta_bar*(np.interp(x,xx,ratio.real)+1j*np.interp(x,xx,ratio.imag))
        p=np.column_stack([np.interp(x,xx,state.p[:,i]) for i in range(state.p.shape[1])])
        n=np.column_stack([np.interp(x,xx,state.n[:,i]) for i in range(state.n.shape[1])])
        return EffectiveState(d,p,n,state.circuit_SI.copy())

    def current_pattern_error(self,state):
        s,g=self.s2,self.g2;f=s.functional
        T=s.kinetics.equivalent_temperature_K(state.delta_bar,state.p)
        ev=f.evaluate(state.delta_bar,state.p,temperature_K=T)
        electrical=s.ports.solve(ev.edge_current_A,state.circuit_SI[1])
        i,j=g.graph.edges.T
        expected=state.circuit_SI[1]/g.width_m*g.dual_face_lengths_m*(
            g.coordinates_m[j,0]-g.coordinates_m[i,0])/g.edge_lengths_m
        weights=np.zeros(len(i));active=g.graph.conductance>0
        weights[active]=s.ports.sheet_resistance_ohm/g.graph.conductance[active]
        residual=electrical.I_total_A-expected
        normal_resistance=s.ports.sheet_resistance_ohm*g.length_m/g.width_m
        return float(np.sqrt(np.dot(weights,residual*residual)/
            (s.circuit_scales[1]**2*normal_resistance)))

    def errors(self,original,reduced,scale):
        lifted=self.lift(reduced);v=self.s2.functional.volumes_m3
        dg=abs(lifted.delta_bar-original.delta_bar)
        return dict(gap_relative=float(np.sqrt(np.dot(v,dg*dg)/
            np.dot(v,abs(self.dc2.delta_bar)**2))),
            populations_energy_relative=population_distance(self.s2,original,lifted,scale),
            film_energy_relative=abs((film_energy(self.s2,original)-self.Edc2)-
                (film_energy(self.s1,reduced)-self.Edc1))/scale,
            current_pattern_relative=self.current_pattern_error(original),
            winding_count=resolved_winding(original.delta_bar,self.g2.triangles,self.dc2.delta_bar))

    def project(self,state,scale,tolerances):
        # Exact DC maps to exact DC: protect the matched discrete references.
        if (np.array_equal(state.delta_bar,self.dc2.delta_bar) and
            np.array_equal(state.p,self.dc2.p) and np.array_equal(state.n,self.dc2.n)):
            result=self.dc1.copy();result.circuit_SI=state.circuit_SI.copy()
            return result,dict(raw_energy_relative=0.,reconciliation_energy_relative=0.,
                               activity_shift=0.,**self.errors(state,result,scale))
        f,k=self.s1.functional,self.s1.kinetics
        d=self.dc1.delta_bar*(self.R@(state.delta_bar/self.dc2.delta_bar))
        p=self.R@state.p;n=self.R@state.n
        # Same rigid superconducting cuts; preserve their phases, never reset CM.
        for old,new in zip(self.s2.cuts,self.s1.cuts):
            relative_phase=np.angle(np.mean(state.delta_bar[old]/self.dc2.delta_bar[old]))
            d[new]=self.dc1.delta_bar[new]*np.exp(1j*relative_phase)
            p[new]=self.dc1.p[new];n[new]=self.dc1.n[new]
        raw=EffectiveState(d,p,n,state.circuit_SI.copy())
        target=film_energy(self.s2,state)-self.Edc2+self.Edc1
        before=film_energy(self.s1,raw)
        raw_error=abs(before-target)/scale
        if raw_error>tolerances.projection_energy_relative:
            raise ArithmeticError("Projection energy tolerance exceeded: "+str(raw_error))
        free=self.s1.free;en=f.energies_bar(d)
        ep=np.zeros_like(p);mask=(p>0)&(p<1)
        ep[mask]=np.log(p[mask])-np.log1p(-p[mask])
        bn=np.full_like(n,np.inf);occupied=n>0
        bn[occupied]=np.log1p(n[occupied])-np.log(n[occupied])
        omega=k.phonons.energies
        def shifted(lam):
            pp=p.copy();nn=n.copy()
            vals=ep[free]+lam*en[free]
            pp[free]=np.where(mask[free],expit(vals),p[free])
            act=bn[free]-lam*omega
            if np.any(act<=0):raise ArithmeticError("Bose activity reached its pole")
            nn[free]=np.exp(-act)/(-np.expm1(-act))
            return EffectiveState(d,pp,nn,state.circuit_SI.copy())
        def residual(lam):return (film_energy(self.s1,shifted(lam))-target)/scale
        shift=0.
        if raw_error>1e-12:
            sign=1. if before<target else -1.
            limit=np.min(bn[free]/omega)
            hi=sign*1e-5;lo=0.;r0=residual(0.)
            for _ in range(35):
                if hi>0 and np.isfinite(limit):hi=min(hi,.8*limit)
                try:rr=residual(hi)
                except (ArithmeticError,ValueError):
                    hi*=.5;rr=residual(hi)
                    if rr*r0>0:raise ArithmeticError("No supported energy reconciliation bracket")
                if rr*r0<=0:break
                hi*=2
            else:raise ArithmeticError("No energy reconciliation bracket")
            shift=brentq(residual,min(lo,hi),max(lo,hi),xtol=1e-13,rtol=1e-12)
        result=shifted(shift)
        correction=population_distance(self.s1,raw,result,scale)
        if correction>tolerances.reconciliation_energy_relative:
            raise ArithmeticError("Energy reconciliation activity too large: "+str(correction))
        if abs(film_energy(self.s1,result)-target)/scale>2e-9:
            raise ArithmeticError("Projection energy was not conserved")
        return result,dict(raw_energy_relative=raw_error,
            reconciliation_energy_relative=correction,activity_shift=float(shift),
            **self.errors(state,result,scale))


class ObservedThermalSystem(ThermalCoupledSystem):
    """Expose already-evaluated accepted midpoint powers without another RHS."""
    def evaluate(self,*args,**kwargs):
        rhs=super().evaluate(*args,**kwargs);self.last_observed_rhs=rhs
        self.last_observed_delta=args[0].delta_bar
        return rhs

    def midpoint_step(self,state,h_ps,**kwargs):
        result,receipt=super().midpoint_step(state,h_ps,**kwargs)
        rhs=self.last_observed_rhs
        f=self.functional
        eps=4*f.node_energy_scale_J[:,None]*f.weights_bar*f.energies_bar(self.last_observed_delta)
        rawpower=np.sum(eps*rhs.kinetics.dp_dt_s,axis=1)+rhs.kinetics.phonon_power_W
        powers=[float(np.sum(rawpower[nodes])) for nodes in self.cuts]
        receipt["absolute_cut_energy_increment_J"]=h_ps*1e-12*sum(abs(x) for x in powers)
        receipt["cut_left_power_W"],receipt["cut_right_power_W"]=powers
        midcm=(state.circuit_SI+result.circuit_SI)/2
        receipt["absolute_device_work_increment_J"]=h_ps*1e-12*abs(midcm[1]*rhs.electrical.Vdev_V)
        return result,receipt


def observed(system):
    # Same constructor/data/physics, only receipt instrumentation differs.
    system.__class__=ObservedThermalSystem
    return system


def combine_receipts(a,b):
    out=dict(b)
    for key in ("external_energy_increment_J","limited_event_energy_J","totalabs_event_energy_J",
                "absolute_cut_energy_increment_J","absolute_device_work_increment_J","absolute_external_energy_increment_J"):
        out[key]=a.get(key,0.)+b.get(key,0.)
    out["accepted_local_energy_fraction"]=max(a["accepted_local_energy_fraction"],b["accepted_local_energy_fraction"])
    out["iterations"]=a.get("iterations",0)+b.get("iterations",0)
    return out


def checked_step(system,state,rhs,h,*,tolerance,max_iterations,energy_tolerance,
                 reference,absolute_work):
    new,receipt=system.midpoint_step(state,h,tolerance=tolerance,max_iterations=max_iterations)
    rr=system.evaluate(new)
    dw=receipt["external_energy_increment_J"]
    defect=abs(rr.energy_J-rhs.energy_J-dw)
    receipt["accepted_local_energy_fraction"]=defect/(reference+absolute_work+abs(dw))
    receipt["absolute_external_energy_increment_J"]=abs(dw)
    if receipt["accepted_local_energy_fraction"]>energy_tolerance:
        raise ArithmeticError("Independent local energy acceptance failed")
    return new,rr,receipt


def doubled_step(system,state,rhs,h,*,time_rtol,time_atol,**kwargs):
    full,_,_=checked_step(system,state,rhs,h,**kwargs)
    half,rhalf,a=checked_step(system,state,rhs,h/2,**kwargs)
    second,rsecond,b=checked_step(system,half,rhalf,h/2,**kwargs)
    f=system.functional;reference=kwargs["reference"]
    gap=float(np.max(abs(second.delta_bar-full.delta_bar))/
              (time_atol+time_rtol*np.max(abs(state.delta_bar))))
    cm=float(np.max(abs(second.circuit_SI-full.circuit_SI)/
        (time_atol*system.circuit_scales+time_rtol*np.maximum(
            abs(second.circuit_SI),system.circuit_scales))))
    native=population_distance(system,second,full,reference)/(time_atol+time_rtol)
    error=max(gap,cm,native)/3.
    if not np.isfinite(error) or error>1:
        raise ArithmeticError("1D step-doubling temporal error: "+str(error))
    receipt=combine_receipts(a,b);receipt["temporal_error"]=error
    return second,rsecond,receipt
