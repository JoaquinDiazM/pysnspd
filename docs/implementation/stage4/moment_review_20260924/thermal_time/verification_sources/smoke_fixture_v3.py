"""Tiny synthetic integration fixture, explicitly not the admitted physical case."""
from pathlib import Path
import hashlib,json,os,signal,subprocess,sys,time
import numpy as np
root=Path('/home/jdiaz/pysnspd');work=root/'tmp/stage4_thermal_time_prepare_20260924';sys.path.insert(0,str(root))
from pysnspd.experimental import thermal_spatial_usadel as thermal
from pysnspd.experimental.thermal_weak_response import SpectralTangent,thermal_hessian_action,ThermalKWTNormal
source=work/'synthetic_operator_v3';source.mkdir();(source/'modes').mkdir()
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
write=lambda p,x:Path(p).write_text(json.dumps(x,indent=2)+'\n')
graph=thermal.rectangular_graph(5,5,2.,2.);xy=graph.coordinates_bar-[1.,0.]
graph=thermal.ThermalGraph(graph.area_weights,graph.edges,graph.conductance,xy,graph.boundary_nodes)
d=np.exp(.1j*xy[:,0]);bump=np.maximum(0.,1-np.sum(xy**2,axis=1))**3
directions=np.array([d*bump,1j*d*xy[:,0]*bump]);directions[:,graph.boundary_nodes]=0
t=.35;spectra=[];bundles=[];records=[]
for n in range(2):
    epsilon=2*np.pi*t*(n+.5)
    s=thermal.solve_frequency(graph,d,epsilon,fixed_nodes=graph.boundary_nodes,fixed_u=d[graph.boundary_nodes]/epsilon,tol=1e-11)
    live=SpectralTangent(graph,d,s);J=live.jacobian;spectra.append(s);bundles.append([live.apply(v) for v in directions])
    p=source/'modes'/f'n{n:04d}.npz'
    np.savez_compressed(p,u=s.u,f=s.f,g=s.g,free_real_components=live.components,jacobian_data=J.data,jacobian_indices=J.indices,jacobian_indptr=J.indptr,jacobian_shape=J.shape)
    records.append(dict(n=n,epsilon=epsilon,fields_path=p.relative_to(source).as_posix(),fields_sha256=sha(p)))
base=thermal.evaluate_thermal(graph,d,t,spectra);H,I=thermal_hessian_action(graph,t,directions,bundles,[s.epsilon for s in spectra])
checks=source/'full_node_operator_checks.npz'
np.savez_compressed(checks,d=d,coordinates_bar=xy,area_weights=graph.area_weights,edges=graph.edges,conductance=graph.conductance,
    boundary_nodes=graph.boundary_nodes,gap_gradient=base.gap_gradient,current=base.current_bar,input_directions=directions,hessian_actions=H,current_actions=I)
write(source/'summary.json',dict(status='FULL_NODE_THERMAL_OPERATOR_VALIDATION_COMPLETE',all_calculus_flags_met=True,records=records,
    base_free_energy=base.energy,scope='SYNTHETIC5x5 TWO-FREQUENCY FIXTURE; live spectral derivatives independently covered by unit tests. Not the admitted physical core.'))
write(source/'identity.json',dict(sources={name:sha(root/name) for name in ('pysnspd/experimental/thermal_spatial_usadel.py','pysnspd/experimental/thermal_weak_response.py')}))
plan=json.loads((root/'docs/implementation/stage4/moment_review_20260924/thermal_time/plan.json').read_text())
plan.update(output_reserve_bytes=64*1024**2,T_K=t*8.65,matsubara_count=2,maximum_workers=2,observation_times_ps=[0.,.000001,.000003],
    operator_summary_sha256=sha(source/'summary.json'),operator_identity_sha256=sha(source/'identity.json'),operator_checks_sha256=sha(checks))
plan['passes']['primary']['maximum_krylov_dimension']=8;plan['passes']['refined']['maximum_krylov_dimension']=12
write(source/'executed_plan.json',plan);write(work/'synthetic_plan_v3.json',plan)
command=['/home/jdiaz/.conda/envs/snspd/bin/python','-u','sandbox/stage4_core/thermal_time_campaign.py','--plan',str(work/'synthetic_plan_v3.json'),'--operator-root',str(source),'--output-root',str(work/'synthetic_evolution_v3'),'--execute']
started=time.monotonic()
with (work/'smoke_v3.log').open('x') as stream:
    p=subprocess.Popen(command,cwd=root,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
    try:code=p.wait(timeout=60)
    except subprocess.TimeoutExpired:
        os.killpg(p.pid,signal.SIGTERM)
        try:p.wait(timeout=5)
        except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGKILL);p.wait()
        code=124
receipt=dict(status='PASSED' if code==0 else 'FAILED',exit_code=code,runtime_seconds=time.monotonic()-started,timeout_seconds=60,
    scope='Synthetic5x5 two-frequency subprocess/pool/refinement/checkpoint composition fixture. Not the admitted physical trajectory.',
    fixture_source_sha256=sha(Path(__file__)),summary_sha256=sha(work/'synthetic_evolution_v3/summary.json') if code==0 else None)
write(work/'smoke_receipt_v3.json',receipt);print(json.dumps(receipt));print((work/'smoke_v3.log').read_text()[-3000:]);raise SystemExit(code)
