"""v1.2.4 spatial photons or authenticated full-1D continuation of v1.2.3 bases.

The native evolution sources remain unchanged. This entry adds reflecting
postcascade preparation and a fresh-directory latency-to-recovery handoff.
"""
from __future__ import annotations
import argparse,importlib.util,json,os,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'

def configure(source_runner,threads):
    os.environ['NUMBA_NUM_THREADS']=str(threads)
    if Path(source_runner).resolve()!=(ROOT/'pipelines/13_dimensional_photon_run.py').resolve():raise ValueError('Current native runner required')
    spec=importlib.util.spec_from_file_location('_spatial_v124_native',source_runner)
    R=importlib.util.module_from_spec(spec);spec.loader.exec_module(R)
    from pysnspd.experimental.reflected_photon import prepare_reflected
    from pysnspd.experimental import latency_continuation as recovery
    from pysnspd.experimental import dimensional_recovery as previous_recovery
    original_pins=R.source_pins;original_write=R.write_json;archive_class=R.PhaseArchive
    original_build=R.OLD.build_system
    def build(args,output,*,resume=False):
        # The legacy DC receipt includes scheduling fields and every parser
        # attribute. Keep its physical comparison exact while allowing the
        # explicitly requested3-thread scheduling and new source-only inputs.
        clean=argparse.Namespace(**{k:v for k,v in vars(args).items() if k not in
            ('impact_y_nm','deposited_energy_eV','source_boundary')})
        identity_function=R.B.parameter_identity
        if resume:
            saved=json.loads((Path(output)/'candidate_physics_identity.json').read_text())['parameter_identity']
            current=identity_function(clean)
            for key in set(saved)|set(current):
                if key not in ('threads','physical_cores') and saved.get(key)!=current.get(key):
                    raise ValueError('Saved setup physical input differs: '+key)
            R.B.parameter_identity=lambda _:saved
        try:return original_build(clean,output,resume=resume)
        finally:R.B.parameter_identity=identity_function
    def pins():
        result=original_pins()
        for rel in sorted(recovery.NEW_SOURCES):result[rel]=R.B.sha256(ROOT/rel)
        return result
    def write(path,value):
        if Path(path).name=='manifest.json':
            value['model_version']='1.2.4'
            value['release_scope']='Native equations unchanged; reflecting source or authenticated saved-latency continuation'
        return original_write(path,value)
    class Archive(archive_class):
        def __init__(self,*args,**kwargs):
            super().__init__(*args,**kwargs);self.file.attrs.modify('model_version','1.2.4')
    R.source_pins=pins;R.write_json=write;R.PhaseArchive=Archive;R.OLD.build_system=build
    previous_recovery.prepare_continuation=recovery.prepare_latency_continuation
    return R,prepare_reflected

def main(argv=None):
    early=argparse.ArgumentParser(add_help=False)
    early.add_argument('--source-runner',type=Path,default=ROOT/'pipelines/13_dimensional_photon_run.py')
    early.add_argument('--threads',type=int,default=3)
    opt,rest=early.parse_known_args(argv);R,prepare=configure(opt.source_runner,opt.threads)
    p=R.parser();p.description=__doc__
    p.set_defaults(sheet_resistance_ohm=608.,retained_fraction=.9,tau_ee_Tc_ps=5.,tau_escape_ps=20.)
    p.add_argument('--impact-y-nm',type=float,default=0.)
    p.add_argument('--deposited-energy-eV',type=float)
    p.add_argument('--source-boundary',choices=('reflecting','preserved'),default='reflecting')
    a=p.parse_args(rest+['--threads',str(opt.threads)])
    if a.recover_from_run is not None or a.resume:
        origin=a.recover_from_run if a.recover_from_run is not None else a.output_dir
        m=json.loads((Path(origin)/'manifest.json').read_text())
        # Adopt all physical/numerical base inputs; CLI supplies only scheduling,
        # session wall budget and output. Recovery authentication rejects changes.
        preserve={k:getattr(a,k) for k in ('output_dir','recover_from_run','resume','threads','cpu_list','physical_cores','max_1d_wall_hours','duration_ps')}
        for key,value in m['parameter_identity'].items():setattr(a,key,value)
        for key,value in preserve.items():setattr(a,key,value)
        if a.recover_from_run is not None:
            a.latency_only=False;a.impact_y_nm=0.;a.deposited_energy_eV=None;a.source_boundary='preserved'
    if a.source_boundary=='preserved' and not (a.recover_from_run is not None or a.resume):p.error('Preserved preparation only applies to saved-state continuation')
    if a.deposited_energy_eV is not None:
        from scipy.constants import elementary_charge,Planck,speed_of_light
        photon=Planck*speed_of_light/(a.wavelength_nm*1e-9)
        a.retained_fraction=a.deposited_energy_eV*elementary_charge/photon
        if not 0<a.retained_fraction<=1:p.error('Retained energy must lie in (0,hc/lambda]')
    if not -40<=a.impact_y_nm<=40:p.error('Impact y must be within the80nm strip')
    if a.source_boundary=='reflecting':
        original=R.OLD.prepare_photon
        def reflected(f,k,d,p0,g,**kwargs):
            return prepare(f,k,d,p0,g,**kwargs,bose_matrix=R.B.bose_matrix,
                impact_y_nm=a.impact_y_nm,deposited_energy_eV=a.deposited_energy_eV)
        R.OLD.prepare_photon=reflected
    try:return R.run(a)
    except Exception as exc:
        path=Path(a.output_dir)/'manifest.json'
        if path.exists():
            m=json.loads(path.read_text());m.update(run_complete=False,error=type(exc).__name__+': '+str(exc))
            if a.recover_from_run is not None or m.get('continuation'):
                # Preserve committed2D and any saved diagnostic shadow on replay
                # failures. Never relabel a failed overlap as admitted1D.
                import h5py,numpy as np
                if 'archive_counts' not in m:
                    previous=Path(a.output_dir)/'pre_update_manifest.json'
                    if previous.exists():m={**json.loads(previous.read_text()),**m}
                    checkpoint=Path(a.output_dir)/'checkpoint_latest.npz'
                    if checkpoint.exists():
                        with np.load(checkpoint,allow_pickle=False) as z:meta=json.loads(str(z['metadata_json']))
                        for key in ('archive_counts','continuation','time_ps','phase'):m[key]=meta[key]
                if (Path(a.output_dir)/'data.h5').exists() and 'archive_counts' in m:
                    with h5py.File(Path(a.output_dir)/'data.h5','r',swmr=True) as f:
                        m['archive_counts']={name:{key:len(f[name][key]) for key in m['archive_counts'][name]} for name in m['archive_counts']}
                m.update(status='POSSIBLE_DIMENSIONAL_ERROR',shadow_error=m['error'])
            else:m['status']='PREPARATION_FAILURE'
            R.write_json(path,m)
        raise

if __name__=='__main__':raise SystemExit(main())
