"""Benchmark and independently check a bounded spectral energy interpolant."""
from pathlib import Path
import argparse,hashlib,json,sys,time
import numpy as np

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from progress import Progress
from pysnspd.experimental.energy_catalog import OccupationEnergyCatalog
from pysnspd.experimental.refined_cells import refined_count_catalog
from pysnspd.experimental.monotone_energy_patch import MonotoneEnergyPatch
from pysnspd.experimental.spatial_functional import PeriodicSpatialFunctional


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def save(path,value):Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def main():
    p=argparse.ArgumentParser();p.add_argument('--output-root',required=True)
    p.add_argument('--registration',default='docs/implementation/stage3/coupled_20260923/monotone_patch_registration.json')
    args=p.parse_args();registration=ROOT/args.registration;reg=json.loads(registration.read_text(encoding='utf-8'))
    output=ROOT/args.output_root
    if output.exists():raise FileExistsError('output exists; no overwrite')
    output.mkdir(parents=True);start=time.perf_counter()
    files=['pysnspd/experimental/'+n+'.py' for n in ['local_energy_patch','monotone_energy_patch','refined_cells','energy_catalog','spatial_functional']]
    files+=['sandbox/stage3_spatial/coupled_20260923/check_monotone_patch.py',args.registration,reg['catalogue']]
    save(output/'manifest.json',dict(registration=reg,sources={f:sha(ROOT/f) for f in files}))
    progress=Progress(len(reg['degrees'])+len(reg['fractional_probes']),'CatÃ¡logo local')
    try:
        source=refined_count_catalog(OccupationEnergyCatalog.load(ROOT/reg['catalogue']),refinement=reg['refinement'])
        models={};patches={};build_seconds={}
        for degree in reg['degrees']:
            progress.start_task('construcciÃ³n grado '+str(degree));begin=time.perf_counter()
            patch=MonotoneEnergyPatch.build(source,reg['amplitude_bounds'],reg['gamma_bounds'],degree,
                on_sample=lambda i,n:progress.checkpoint(f'espectro {i}/{n}'))
            patch.save(output/f'patch_degree{degree}.npz');patches[degree]=patch
            build_seconds[degree]=time.perf_counter()-begin
            models[degree]=PeriodicSpatialFunctional(patch,360e-9,840e-18,5);progress.advance()
        direct_model=PeriodicSpatialFunctional(source,360e-9,840e-18,5)
        records=[];direct_times=[];patch_times=[]
        for index,(fa,fg) in enumerate(reg['fractional_probes']):
            progress.start_task(f'contraste causal independiente {index+1}')
            a=reg['amplitude_bounds'][0]+fa*np.diff(reg['amplitude_bounds'])[0]
            g=reg['gamma_bounds'][0]+fg*np.diff(reg['gamma_bounds'])[0]
            begin=time.perf_counter();exact=np.asarray(source.energy_kernel(a,g));direct_times.append(time.perf_counter()-begin)
            occupied=1/(np.exp(np.minimum(exact[0]/reg['theta'],700))+1)+.001*np.exp(-((source.count_nodes-.035)/.01)**2)
            d=1j*a*np.sqrt(g*direct_model.gap_ratio)*(a*a+.01)/(a*a)
            reference=direct_model.principal_symbol(a,d,occupied)
            point=dict(amplitude=float(a),gamma=float(g),degrees={})
            for degree,patch in patches.items():
                begin=time.perf_counter()
                for _ in range(reg['benchmark_repetitions']):approx=np.asarray(patch.energy_kernel(a,g))
                per_call=(time.perf_counter()-begin)/reg['benchmark_repetitions']
                if degree==reg['selected_degree']:patch_times.append(per_call)
                delta=approx-exact;w=source.count_weights
                kernel_error=np.sum(abs(delta)*w,axis=1)/np.sum(abs(exact)*w,axis=1)
                moment_error=abs(delta@(w*occupied))/np.sum(abs(exact)*w*occupied,axis=1)
                symbol=models[degree].principal_symbol(a,d,occupied)
                relative_matrix=float(np.linalg.norm(symbol.matrix-reference.matrix,2)/np.linalg.norm(reference.matrix,2))
                sign_margin=float(min(symbol.eigenvalues[0],reference.eigenvalues[0])-symbol.uncertainty-reference.uncertainty-
                                  np.linalg.norm(symbol.matrix-reference.matrix,2))
                metrics=dict(pointwise_relative_energy=float(np.max(abs(delta[0]/exact[0]))),
                    weighted_relative_kernel_L1=float(np.max(kernel_error)),
                    moment_error_over_absolute_weighted_reference=float(np.max(moment_error)),
                    relative_principal_matrix=relative_matrix)
                passed=all(metrics[k]<=v for k,v in reg['limits'].items()) and sign_margin>0
                point['degrees'][str(degree)]=dict(metrics=metrics,kernel_errors=kernel_error.tolist(),
                    moment_errors=moment_error.tolist(),sign_margin=sign_margin,passed=bool(passed),
                    average_kernel_seconds=per_call)
            records.append(point);save(output/'probes.json',records);progress.advance()
        candidate=str(reg['selected_degree']);passed=all(x['degrees'][candidate]['passed'] for x in records)
        result=dict(status='PASS_LOCAL_PATCH_SAMPLED' if passed else 'LOCAL_PATCH_NOT_ACCEPTED',
            runtime_seconds=time.perf_counter()-start,build_seconds=build_seconds,
            selected_degree=reg['selected_degree'],probes=records,
            max_candidate_errors={k:max(x['degrees'][candidate]['metrics'][k] for x in records) for k in reg['limits']},
            direct_median_kernel_seconds=float(np.median(direct_times)),
            patch_median_kernel_seconds=float(np.median(patch_times)),
            kernel_speedup=float(np.median(direct_times)/np.median(patch_times)),
            speedup_scope='Kernel lookup only, not whole RHS or trajectory speedup.',
            stage3_closed=False,production=False,patch_sha256=sha(output/f'patch_degree{candidate}.npz'))
        save(output/'summary.json',result);progress.finish(success=passed)
        print(json.dumps({k:v for k,v in result.items() if k not in ['probes','build_seconds']},indent=2))
        if not passed:raise SystemExit(2)
    except Exception as error:
        save(output/'failure.json',dict(exception=type(error).__name__,reason=str(error),runtime_seconds=time.perf_counter()-start))
        raise


if __name__=='__main__':main()
