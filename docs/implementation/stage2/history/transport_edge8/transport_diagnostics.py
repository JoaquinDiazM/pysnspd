"""Independent fixed-energy transport consistency and native-moment diagnostics.

No time trajectory or production module is edited. A denser reference cell is
used only to diagnose the interpolation limit, never to label the native R2
180-node response as passing when its own error fails the contract.
"""
from __future__ import annotations
import argparse
from dataclasses import dataclass
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import sys
import time

import numpy as np
from scipy.special import expit

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from pysnspd.experimental.cell_validation import ElectronicCell,occupation_array
from pysnspd.experimental.energy_catalog import OccupationEnergyCatalog,retarded_spectrum
import pysnspd.experimental.cell_transport as current_module

FIELDS={"normal":((0.,0.),(0.,0.)),
        "gapped":((.55,0.),(.95,0.)),
        "gapless":((.35,.65),(.6,.8))}


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")


@dataclass
class DenseReferenceCell:
    """Interpolation study over the same energy hull, not a new R2 catalogue."""
    base: ElectronicCell
    factor: int

    def __post_init__(self):
        native=self.base.energies
        fractions=np.arange(self.factor)/self.factor
        self.energies=np.r_[(native[:-1,None]+np.diff(native)[:,None]*fractions).ravel(),native[-1]]
        # Count control volumes are diagnostic only. In weighted-L1 comparison
        # these positive capacities cancel; they are not a dynamics replacement.
        count=np.interp(self.energies,native,self.base.catalog.count_nodes)
        bounds=np.r_[count[0],(count[:-1]+count[1:])/2,count[-1]]
        self.weights=np.diff(bounds)
        self.catalog=self.base.catalog; self.amplitude=self.base.amplitude; self.gamma=self.base.gamma

    def validate(self,p): return occupation_array(p,self.energies.shape)


def population(energy,label,side):
    energy=np.asarray(energy)
    if label=="thermal":return expit(-energy/(.35 if side==0 else .12))
    if label=="nonthermal":
        if side==0:
            logodds=-energy/.3+1.5*np.exp(-.5*((energy-2.5)/.4)**2)
        else:
            logodds=-energy/.18-.8*np.exp(-.5*((energy-3.)/.55)**2)
        return expit(logodds)
    raise ValueError(label)


def reference_geometry(cells,order):
    """Geometry/spectrum shared by profiles; independent of the event class."""
    lo=max(c.energies[0] for c in cells); hi=min(c.energies[-1] for c in cells)
    edges=np.concatenate([np.array([lo,hi])]+[c.energies for c in cells])
    # The source operator is not reused: add known spectral edges explicitly.
    for cell in cells:
        if cell.amplitude>cell.gamma:
            edge=(cell.amplitude**(2/3)-cell.gamma**(2/3))**1.5
            edges=np.r_[edges,edge]
    edges=np.unique(edges[(edges>=lo)&(edges<=hi)])
    z,w=np.polynomial.legendre.leggauss(order)
    q=(edges[:-1,None]+np.diff(edges)[:,None]*(z+1)/2).ravel()
    weights=(np.diff(edges)[:,None]*w/2).ravel()
    diffusion=[]
    for cell in cells:
        c,s=retarded_spectrum(q,delta=cell.amplitude,gamma=cell.gamma,eta=cell.catalog.eta)
        dl=np.real(c)**2-np.imag(s)**2
        if np.any(dl<0):raise ValueError("reference diffusion coefficient is negative")
        diffusion.append(dl)
    dl,dr=diffusion
    shared=np.divide(2*dl*dr,dl+dr,out=np.zeros_like(q),where=(dl+dr)>0)
    return q,weights,shared


def reference(cells,label,coefficient,order,geometry=None):
    """Direct DL(f_L-f_R) integral and independent weak nodal deposition."""
    q,weights,shared=reference_geometry(cells,order) if geometry is None else geometry
    flux=4*coefficient*weights*shared*(population(q,label,0)-population(q,label,1))
    rhs=[]
    for side,cell in enumerate(cells):
        upper=np.searchsorted(cell.energies,q,side="left"); lower=upper-1
        if np.any(lower<0) or np.any(upper>=len(cell.energies)):raise ValueError("reference would extrapolate")
        fraction=(q-cell.energies[lower])/(cell.energies[upper]-cell.energies[lower])
        count_rate=np.zeros_like(cell.energies)
        np.add.at(count_rate,lower,(1-fraction)*flux)
        np.add.at(count_rate,upper,fraction*flux)
        rhs.append((2*side-1)*count_rate/(4*cell.weights))
    return rhs,float(np.dot(q,flux)),float(np.sum(np.abs(flux)))


def population_norm(cells,left,right):
    numerator=sum(float(np.dot(4*c.weights,np.abs(a-b))) for c,a,b in zip(cells,left,right))
    denominator=sum(float(np.dot(4*c.weights,np.abs(b))) for c,b in zip(cells,right))
    return {"absolute_weighted_L1":numerator,"reference_weighted_L1":denominator,
            "relative_weighted_L1":numerator/denominator if denominator>0 else None}


def algebra(operator,cells):
    p=[population(c.energies,"nonthermal",i) for i,c in enumerate(cells)]
    dp=operator.rhs(p)
    counts=[float(np.dot(4*c.weights,r)) for c,r in zip(cells,dp)]
    energies=[float(np.dot(4*c.weights*c.energies,r)) for c,r in zip(cells,dp)]
    entropy=sum(float(np.dot(4*c.weights*np.log((1-f)/f),r)) for c,f,r in zip(cells,p,dp))
    fd=[]
    for theta in [.08,.2,.5]:
        peq=[expit(-c.energies/theta) for c in cells]
        rhs=operator.rhs(peq)
        fd.append({"temperature_bar":theta,"weighted_RHS_L1":sum(float(np.dot(4*c.weights,np.abs(r))) for c,r in zip(cells,rhs)),
                   "absolute_power":abs(operator.transferred_power(peq))})
    bounds=[]
    for kind in ["empty_left_full_right","alternating_faces"]:
        if kind=="empty_left_full_right":pop=[np.zeros_like(cells[0].energies),np.ones_like(cells[1].energies)]
        else:
            pop=[np.full_like(c.energies,.25) for c in cells]
            for state in pop:state[::5]=0.;state[1::7]=1.
        rhs=operator.rhs(pop)
        inward=[r[f==0] for f,r in zip(pop,rhs)]+[-r[f==1] for f,r in zip(pop,rhs)]
        minimum=min(float(np.min(v)) for v in inward if len(v))
        bounds.append({"case":kind,"minimum_inward_population_derivative":minimum})
    return {"event_energy_residual_absolute":float(np.max(np.abs(operator.event_energy_residual()))),
        "QP_count_balance_scaled":abs(sum(counts))/max(1.,sum(abs(v) for v in counts)),
        "native_energy_balance_scaled":abs(sum(energies))/max(1.,sum(abs(v) for v in energies)),
        "entropy_production":entropy,"common_FD":fd,"Pauli_boundaries":bounds}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog",type=Path,default=ROOT/"docs/implementation/stage1_r2/catalogs/occupation_catalog.npz")
    parser.add_argument("--output",type=Path,default=ROOT/"docs/implementation/stage2/transport_results.json")
    parser.add_argument("--operator-file",type=Path)
    parser.add_argument("--fields",nargs="+",choices=tuple(FIELDS),default=list(FIELDS))
    parser.add_argument("--refinement-factors",nargs="+",type=int,default=[1,2,4])
    parser.add_argument("--face-orders",nargs="+",type=int,default=[1,2,4])
    parser.add_argument("--complementary-refinements",nargs="+",type=int,
                        help="Use actual refined causal count catalogues instead of diagnostic energy subdivisions")
    parser.add_argument("--bulk-spacing",type=float,default=.02)
    args=parser.parse_args(); started=time.perf_counter()
    source=Path(current_module.__file__) if args.operator_file is None else args.operator_file
    if args.operator_file is None:
        Operator=current_module.NativeTransportEvents
    else:
        name="pysnspd.experimental._transport_review_snapshot"
        spec=importlib.util.spec_from_file_location(name,source); saved=importlib.util.module_from_spec(spec)
        sys.modules[name]=saved; spec.loader.exec_module(saved);Operator=saved.NativeTransportEvents
    criteria=ROOT/"docs/implementation/stage2/acceptance_criteria.json"
    if args.complementary_refinements is not None:
        from pysnspd.experimental import refined_cells
        if any(level<1 for level in args.complementary_refinements):
            raise ValueError("complementary refinements must be positive")
    registration={"schema":"pysnspd.stage2.transport-scenarios.v2","operator_sha256":sha(source),
        "catalog_sha256":sha(args.catalog),"criteria_sha256":sha(criteria),
        "fields":{name:FIELDS[name] for name in args.fields},
        "refinement_factors":args.refinement_factors if args.complementary_refinements is None else None,
        "complementary_refinements":args.complementary_refinements,
        "complementary_source_sha256":None if args.complementary_refinements is None else sha(refined_cells.__file__),
        "complementary_bulk_spacing":args.bulk_spacing if args.complementary_refinements else None,
        "face_orders":args.face_orders,"diffusion_over_length_squared":.03,
        "profiles":{"thermal":"FD left Tbar=.35, right=.12",
            "nonthermal":"logodds left=-E/.3+1.5 exp[-.5((E-2.5)/.4)^2]; right=-E/.18-.8 exp[-.5((E-3)/.55)^2]"},
        "reference_orders":[24,48],"native_population_gate_relative":1e-3,
        "reference_refinement_budget_relative":2e-4,
        "dense_reference_scope":"Linear subdivisions of each native energy interval are diagnostic only. With --complementary-refinements, each level instead has its own positive count Gauss rule, exact causal energies and consistent energy derivatives. This does not change the historical native180 failure.",
        "time_dynamics":False}
    write(args.output.with_name(args.output.stem+"_scenarios.json"),registration)
    catalog=OccupationEnergyCatalog.load(args.catalog)
    rows=[]; exact=[]
    for name in args.fields:
        base=[ElectronicCell(catalog,*field) for field in FIELDS[name]]
        levels=args.refinement_factors if args.complementary_refinements is None else args.complementary_refinements
        for factor in levels:
            if args.complementary_refinements is None:
                cells=base if factor==1 else [DenseReferenceCell(c,factor) for c in base]
                state_mesh="native_R2" if factor==1 else "diagnostic_energy_subdivision"
            else:
                complementary=refined_cells.refined_count_catalog(catalog,factor,bulk_spacing=args.bulk_spacing)
                cells=[ElectronicCell(complementary,*field) for field in FIELDS[name]]
                state_mesh="complementary_causal_count"
            exact.append({"field":name,"state_mesh":state_mesh,"electron_factor":factor,
                          **algebra(Operator(*cells,gauss_order=4),cells)})
            geometries={order:reference_geometry(cells,order) for order in [24,48]}
            for profile in ["thermal","nonthermal"]:
                low=reference(cells,profile,.03,24,geometries[24]); high=reference(cells,profile,.03,48,geometries[48])
                refcheck=population_norm(cells,low[0],high[0])
                for order in args.face_orders:
                    operator=Operator(*cells,gauss_order=order)
                    p=[population(c.energies,profile,i) for i,c in enumerate(cells)]
                    rhs=operator.rhs(p); power=operator.transferred_power(p)
                    comparison=population_norm(cells,rhs,high[0])
                    rows.append({"field":name,"profile":profile,"electron_factor":factor,"state_mesh":state_mesh,
                        "electron_nodes":[len(c.energies) for c in cells],"face_order":order,"event_count":len(operator.energies),
                        "shared_support":list(operator.shared_support),"RHS_comparison":comparison,
                        "power":power,"reference_power":high[1],"power_relative_error":abs(power-high[1])/abs(high[1]),
                        "reference_RHS_refinement":refcheck,"reference_power_refinement_relative":abs(low[1]-high[1])/abs(high[1]),
                        "passes_population_0_1_percent":comparison["relative_weighted_L1"]<=1e-3,
                        "passes_power_0_1_percent":abs(power-high[1])/abs(high[1])<=1e-3})
    native=[r for r in rows if r["state_mesh"]=="native_R2" and r["face_order"]==max(args.face_orders)]
    selected=[r for r in rows if r["electron_factor"]==max(levels) and r["face_order"]==max(args.face_orders)]
    summaries=[]
    for level in levels:
        group=[r for r in rows if r["electron_factor"]==level and r["face_order"]==max(args.face_orders)]
        summaries.append({"electron_factor":level,"electron_nodes":group[0]["electron_nodes"],
            "state_mesh":group[0]["state_mesh"],
            "maximum_RHS_relative_error":max(r["RHS_comparison"]["relative_weighted_L1"] for r in group),
            "maximum_power_relative_error":max(r["power_relative_error"] for r in group),
            "precision_status":"PASS" if all(r["passes_population_0_1_percent"] and r["passes_power_0_1_percent"] for r in group) else "FAIL"})
    result={"schema":"pysnspd.stage2.transport-review.v2","source_path":str(source),
        "source_sha256":sha(source),"runner_sha256":sha(__file__),"catalog_sha256":sha(args.catalog),
        "criteria_sha256":sha(criteria),"registration_sha256":sha(args.output.with_name(args.output.stem+"_scenarios.json")),
        "complementary_source_sha256":registration["complementary_source_sha256"],
        "algebra":exact,"continuous_comparisons":rows,
        "native_status":("PASS" if all(r["passes_population_0_1_percent"] and r["passes_power_0_1_percent"] for r in native) else "FAIL") if native else "NOT_RETESTED_HISTORICAL_FAILURE_PRESERVED",
        "native_max_relative_RHS_error":max((r["RHS_comparison"]["relative_weighted_L1"] for r in native),default=None),
        "native_max_relative_power_error":max((r["power_relative_error"] for r in native),default=None),
        "selected_precision_status":"PASS" if all(r["passes_population_0_1_percent"] and r["passes_power_0_1_percent"] for r in selected) else "FAIL",
        "mesh_summary":summaries,
        "reference_max_relative_RHS_refinement":max(r["reference_RHS_refinement"]["relative_weighted_L1"] for r in rows),
        "runtime_seconds":time.perf_counter()-started,"platform":platform.platform(),"python":platform.python_version()}
    write(args.output,result)
    print(json.dumps({key:result[key] for key in ["native_status","selected_precision_status","mesh_summary","reference_max_relative_RHS_refinement","runtime_seconds"]},indent=2))


if __name__=="__main__":main()
