"""Append-only phase-aware HDF5 plus atomic accepted-state checkpoints."""
import json
from pathlib import Path
import h5py
import numpy as np

TRACE=("time_ps","Vout_V","Vdev_V","Ib_A","Is_A","vc_V","energy_J",
       "energy_defect_J","external_J","absolute_external_J","cut_absolute_J",
       "device_absolute_J","domain_fraction","escape_power_W","cuts_power_W",
       "step_ps","iterations","temporal_error")


class PhaseArchive:
    def __init__(self,path,systems,geometries,dc_states,*,resume=False,counts=None):
        self.path=Path(path)
        self.file=h5py.File(self.path,"r+" if resume else "w",libver="latest")
        if not resume:
            self.file.attrs["schema"]="pysnspd.dimensional_run.v1"
            self.file.attrs["model_version"]="1.2.1"
            self.file.attrs["voltage_scope"]="Intrinsic Rload*(Ib-Is); no amplifier/filter"
            self.file.attrs["phase_policy"]="2d official before handoff; shadow1d independent; 1d official after handoff"
            for name,system in systems.items():
                g=geometries[name];dc=dc_states[name];f=system.functional;k=system.kinetics
                q=self.file.create_group(name)
                for key,val in dict(coordinates_m=g.coordinates_m,triangles=g.triangles,
                    edges=g.graph.edges,volumes_m3=f.volumes_m3,dc_delta_bar=dc.delta_bar,
                    dc_p=dc.p,dc_n=dc.n,xi_bar=f.xi_bar,xi_weights=f.weights_bar,
                    phonon_energies_bar=k.phonons.energies).items():q.create_dataset(key,data=val)
                q.attrs["energy_scale_J"]=f.Eu_J
                for field in ("Lk_ext_H","R_bias_ohm","L_bias_H","R_load_ohm","C_couple_F","V_bias_V"):
                    q.attrs[field]=getattr(system.circuit,field)
                for field in ("core_begin_m","core_end_m"):
                    if field in g.receipt:q.attrs[field]=g.receipt[field]
                q.create_dataset("trace",shape=(0,),maxshape=(None,),
                    dtype=[(key,"f8") for key in TRACE],chunks=(256,))
                n=g.graph.n_nodes
                for key,shape,dtype in (
                    ("snapshot_time_ps",(),float),("delta_bar",(n,),complex),
                    ("p",(n,k.electron_size),float),("n",(n,k.phonon_size),float),
                    ("Te_K",(n,),float),("u_excess_J_m3",(n,),float)):
                    q.create_dataset(key,shape=(0,)+shape,maxshape=(None,)+shape,
                        chunks=(1,)+shape,dtype=dtype,compression="lzf")
        elif counts is not None:
            for name,sizes in counts.items():
                q=self.file[name]
                for key,size in sizes.items():
                    ds=q[key]
                    if len(ds)<size:raise ValueError("Archive is shorter than committed checkpoint")
                    ds.resize(size,axis=0)
        self.file.flush();self.file.swmr_mode=True
        self.reference={}
        for name,system in systems.items():
            f,k=system.functional,system.kinetics;dc=dc_states[name]
            T=k.equivalent_temperature_K(dc.delta_bar,dc.p)
            self.reference[name]=f.evaluate(dc.delta_bar,dc.p,temperature_K=T).node_energy_J+np.sum(k.phonon_energy_weights_J*dc.n,axis=1)

    def append(self,name,row):
        q=self.file[name]["trace"];i=len(q);q.resize(i+1,axis=0)
        q[i]=tuple(float(row.get(key,0.)) for key in TRACE)

    def snapshot(self,name,system,state,rhs,t):
        f,k=system.functional,system.kinetics
        values=dict(snapshot_time_ps=t,delta_bar=state.delta_bar,p=state.p,n=state.n,
            Te_K=rhs.temperature_K,u_excess_J_m3=(rhs.evaluation.node_energy_J+
                np.sum(k.phonon_energy_weights_J*state.n,axis=1)-self.reference[name])/f.volumes_m3)
        q=self.file[name];i=len(q["snapshot_time_ps"])
        if i and abs(float(q["snapshot_time_ps"][-1])-t)<1e-10:return
        for key,val in values.items():
            ds=q[key];ds.resize(i+1,axis=0);ds[i]=val

    def counts(self):
        return {name:{key:len(ds) for key,ds in q.items() if key in
            ("trace","snapshot_time_ps","delta_bar","p","n","Te_K","u_excess_J_m3")}
            for name,q in self.file.items()}

    def close(self):self.file.close()
