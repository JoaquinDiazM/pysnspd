"""Reproduce static stage4A audit without re-running the user batches.

Checks every recorded result hash, compares analytic prescribed Gamma, separates
outer-boundary algebraic heat, and performs three direct center-symbol queries.
This changes no original output or physical solver and advances no trajectory.
Typical center-query runtime is under 30 seconds. Output must be new.
"""
from pathlib import Path
import argparse, hashlib, json, collections, sys, time
import numpy as np

R = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--raw-root', type=Path, default=R/'docs/implementation/stage4/review_20260923/raw',
                    help='Parent of stage4A_local_20260923 and stage4A_spatial_20260923')
parser.add_argument('--output-dir', type=Path, required=True)
args = parser.parse_args()
B, O = args.raw_root, args.output_dir
O.mkdir(parents=True, exist_ok=True)
if (O/'independent_audit.json').exists():
    raise FileExistsError('Choose a new output directory; previous audit is preserved')
plan = R/'docs/implementation/stage4/start_20260923/campaign_plan.json'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
archive_root=R/'docs/implementation/stage4/review_20260923/previous_delivery_exact'
def source_lookup(path, digest):
    current=R/path
    if current.exists() and sha(current)==digest:
        return 'current'
    archived=archive_root/path
    return 'archived' if archived.exists() and sha(archived)==digest else 'unverified'

audit={'schema':'pysnspd.stage4A.independent_review.v1','date':'2026-09-23','meaning':'Audit of completed prescribed static fields; no new physical trajectory, no production admission','stage4_complete':False,'photon':False,'identity':{},'local':{},'spatial':{},'new_light_diagnostic':{},'recommendations':[]}
results={}
for phase in ['local','spatial']:
 d=B/f'stage4A_{phase}_20260923'; s=json.loads((d/'summary.json').read_text());i=json.loads((d/'identity.json').read_text());rs={p['id']:json.loads((d/p['result']).read_text()) for p in s['cases']};results[phase]=rs
 missing=[];bad=[];nfile=0
 for case in s['cases']:
  for f in case['files']:
   p=d/f['path'];nfile+=1
   if not p.exists():missing.append(f['path'])
   elif sha(p)!=f['sha256'] or p.stat().st_size!=f['bytes']:bad.append(f['path'])
 audit['identity'][phase]={'cases_expected':len(i['task_ids']),'cases_completed':len(rs),'task_ids_match':list(rs)==i['task_ids'],'plan_hash_matches':sha(plan)==i['plan_sha256'],'source_hashes_match_current':all(sha(R/p)==h for p,h in i['sources'].items()),'source_hashes_current_or_archived':{p:source_lookup(p,h) for p,h in i['sources'].items()},'all_source_hashes_verified':all(source_lookup(p,h)!='unverified' for p,h in i['sources'].items()),'source_files':len(i['sources']),'output_files_checked':nfile,'missing':missing,'hash_or_size_mismatch':bad,'batch_runtime_seconds':s['runtime_seconds'],'summary_sha256':sha(d/'summary.json'),'identity_sha256':sha(d/'identity.json')}
loc=results['local'];audit['local']['verdict_counts']=dict(collections.Counter(x['verdict'] for x in loc.values()));audit['local']['negative_controls']=[{'id':k,'lambda_min':x['symbol_eigenvalues'][0],'uncertainty':x['symbol_uncertainty'],'population':x['population'],'delta_reg':x['delta_reg']} for k,x in loc.items() if not x['symbol_stable']]
audit['local']['fd_max_absolute_error_final_step']=max(x['fd_identity']['absolute_errors'][-1] for x in loc.values());audit['local']['kwt_max_absolute_power_identity_residual']=max(abs(m['power_identity_residual']) for x in loc.values() for m in x['mobility'].values())
anchors=[]
for control in ['current','small']:
 a=loc['local_'+control+'_synthetic_fixed_count_d0.1']
 for resolution in ['count_refined','eta_halved']:
  b=loc['anchor_'+control+'_'+resolution]; fa=np.array(a['force_cartesian']);fb=np.array(b['force_cartesian'])
  anchors.append({'control':control,'resolution':resolution,'energy_absolute_change':abs(b['energy_density_bar']-a['energy_density_bar']),'force_relative_change':float(np.linalg.norm(fb-fa)/np.linalg.norm(fa)),'lambda_min_absolute_change':abs(b['symbol_eigenvalues'][0]-a['symbol_eigenvalues'][0]),'equivalent_temperature_K_absolute_change':abs(b['equivalent_temperature_K']-a['equivalent_temperature_K'])})
audit['local']['spectral_anchors']=anchors
material=json.loads((B/'stage4A_spatial_20260923/identity.json').read_text())['material'];volunit=material['width_m']*material['thickness_m']*material['ell0_m'];gap=2.1065301756880164e-22/(1.380649e-23*material['Tc_K'])
for key,x in results['spatial'].items():
 p=B/'stage4A_spatial_20260923'/key/'fields.npz';f=np.load(p);m=f['quadrature_volume_m3']/volunit;shape=f['shape'];ix=np.indices(shape);bnd=((ix[0]==0)|(ix[1]==0)|(ix[0]==shape[0]-1)|(ix[1]==shape[1]-1)).ravel();gd=f['cartesian_force']/m[:,None];xy=f['coordinates_m']/material['ell0_m'];xx=xy[:,0]-160e-9/material['ell0_m']/2;yy=xy[:,1];amp=abs(f['delta']);px=.07-.0005*xx*np.sin(yy/5)*np.exp(-(xx/10)**2);py=.005*np.cos(yy/5)*np.exp(-(xx/10)**2);exact=(amp*amp/(amp*amp+x['delta_reg']**2))**2*(px*px+py*py)/gap;diff=f['gamma']-exact;center=int(np.argmin(xx*xx+yy*yy));negative=f['symbol_eigenvalues'][:,0]<-f['symbol_uncertainty']
 q={'nodes':x['nodes'],'runtime_seconds':x['runtime_seconds'],'energy_bar':x['energy_bar'],'gamma_relative_mass_weighted_L2_error_vs_analytic_profile':float(np.sqrt(np.dot(m,diff*diff)/np.dot(m,exact*exact))),'gamma_max_absolute_error_vs_analytic_profile':float(max(abs(diff))),'center_gamma_discrete':float(f['gamma'][center]),'center_gamma_analytic':float(exact[center]),'negative_symbols':int(sum(negative)),'negative_coordinates_nm':(f['coordinates_m'][negative]*1e9).tolist(),'negative_lambda_values':f['symbol_eigenvalues'][negative,0].tolist(),'negative_local_uncertainties':f['symbol_uncertainty'][negative].tolist(),'regions':{}}
 for mask,label in [(np.ones(len(m),bool),'all'),(bnd,'boundary'),(~bnd,'interior')]:
  q['regions'][label]={'nodes':int(sum(mask)),'mass_bar':float(m[mask].sum()),'force_density_L2':float(np.sqrt(np.sum(m[mask]*np.sum(gd[mask]**2,axis=1)))),'force_density_RMS':float(np.sqrt(np.sum(m[mask]*np.sum(gd[mask]**2,axis=1))/m[mask].sum())),'instantaneous_unconstrained_heat':{k.removeprefix('heat_'):float(m[mask]@f[k][mask]) for k in f.files if k.startswith('heat_')}}
 audit['spatial'][key]=q
comparisons={}
for profile in ['smooth','suppressed']:
 a=audit['spatial'][f'2d_{profile}_d0.1_4x2'];b=audit['spatial'][f'2d_{profile}_d0.1_8x4'];ar=a['regions'];br=b['regions'];comparisons[profile]={'energy_relative_change':(b['energy_bar']-a['energy_bar'])/abs(a['energy_bar']),'interior_force_L2_relative_change':br['interior']['force_density_L2']/ar['interior']['force_density_L2']-1,'Korzh_interior_heat_relative_change':br['interior']['instantaneous_unconstrained_heat']['Korzh_reference']/ar['interior']['instantaneous_unconstrained_heat']['Korzh_reference']-1,'Korzh_boundary_heat_ratio':br['boundary']['instantaneous_unconstrained_heat']['Korzh_reference']/ar['boundary']['instantaneous_unconstrained_heat']['Korzh_reference'],'Korzh_boundary_heat_fraction_fine':br['boundary']['instantaneous_unconstrained_heat']['Korzh_reference']/br['all']['instantaneous_unconstrained_heat']['Korzh_reference']}
audit['mesh_comparisons']=comparisons
audit['norm_convention']={'force_density':'F_i=G_i/m_i, where G_i is the integrated Cartesian nodal gradient and m_i is quadrature mass','L2':'sqrt(sum_i m_i |F_i|^2)','RMS':'L2 / sqrt(sum_i m_i)','boundary_mask':'All four outer GLL node rows; internal element boundaries remain interior','caution':'Post-processing only, not an implemented boundary condition or projected transient; heat is stored algebraic response before constraints','edge_current_caution':'maximum_current_A is an integrated graph-edge current, not a mesh-invariant current density or whole-section current'}
audit['new_light_diagnostic']={'purpose':'Isolate underresolved spatial derivative from constitutive principal symbol at prescribed suppressed center','analytic_field':'z=(0.95-0.91 exp[-(x/2.5)^2-(y/2.5)^2]) exp[i(0.07x+0.025 sin(y/5) exp[-(x/10)^2])], centered dimensionless x,y','analytic_center':{'z':[.04,0.],'derivatives':[[0.,.0028],[0.,.0002]]},'population':'synthetic_fixed_count','spectral_resolution':'reference','method':'Existing controls.catalogue, population, RectangularSpatialFunctional.principal_symbol with analytic center derivatives, no time advance','results':[], 'kernel_runtime_seconds':[],'conclusion':'The coarse delta=.05 negative center is not evidence for a negative continuum symbol of that prescribed field; its gamma is about 109.78 times the analytic value. Do not increase delta to mask spatial underresolution.'}
# Direct constitutive reference with the exact derivatives of the input field.
sys.path.insert(0, str(R))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from controls import catalogue, population
from pysnspd.experimental.rectangular_spatial import RectangularSpatialFunctional
p = json.loads(plan.read_text(encoding='utf8'))
cat = catalogue(p, 'reference')
occupation = population(cat, 'synthetic_fixed_count')
z = .04+0j
derivatives = np.array([.0028j, .0002j])
center_results, elapsed = [], []
for delta in (.05, .1, .2):
    started = time.monotonic()
    model = RectangularSpatialFunctional(cat, 160e-9, 80e-9, 7e-9,
        elements_x=1, elements_y=1, degree=2, delta_regularizer_bar=delta)
    q = np.imag(np.conj(z)*derivatives)/(abs(z)**2+delta**2)
    symbol = model.principal_symbol(z, derivatives, occupation, spatial_dimensions=2)
    center_results.append(dict(delta_reg=delta, gamma=float(q@q/model.gap_ratio),
        lambda_min=float(symbol.eigenvalues[0]), eigenvalues=symbol.eigenvalues.tolist(),
        uncertainty=float(symbol.uncertainty), positive=bool(symbol.stable)))
    elapsed.append(time.monotonic()-started)
audit['new_light_diagnostic']['results'] = center_results
audit['new_light_diagnostic']['kernel_runtime_seconds'] = elapsed
audit['new_light_diagnostic']['script_sha256'] = sha(Path(__file__))
audit['spatial_identity_residuals'] = {
    'max_noether_absolute':max(x['noether_max'] for x in results['spatial'].values()),
    'max_kwt_power_absolute':max(abs(m['identity_residual_bar']) for x in results['spatial'].values() for m in x['mobility'].values()),
    'field_fd_max_absolute':max(abs(c['field_derivative']-c['expected_field']) for x in results['spatial'].values() for c in x['identity_checks']),
    'gauge_fd_max_absolute':max(abs(c['gauge_derivative']-c['expected_gauge']) for x in results['spatial'].values() for c in x['identity_checks']),
    'field_fd_max_relative':max(abs(c['field_derivative']-c['expected_field'])/abs(c['expected_field']) for x in results['spatial'].values() for c in x['identity_checks']),
    'gauge_fd_max_relative':max(abs(c['gauge_derivative']-c['expected_gauge'])/abs(c['expected_gauge']) for x in results['spatial'].values() for c in x['identity_checks']),
}
audit['recommendations']=['Accept identity/source reproduction and constitutive directional derivative/KWT algebra diagnostics; do not claim full stage4 completion.','Resolve prescribed narrow-core gradients with analytic derivative comparison and spatial refinement, including delta=.05 and .1; four independent local high-gradient negative controls remain valid exclusions.','Declare actual boundary constraints before interpreting total KWT heating or evolving. Current unconstrained boundary reactions dominate heat; removing boundary samples in plots alone does not implement those constraints.','Use volume-weighted interior force and interior dissipation plus the analytic-profile Gamma error to separate interior convergence from boundary reactions; do not use Euclidean integrated force or max edge current as mesh-invariant observables.','After boundary implementation and resolved spatial field, use a bounded smooth weak non-photon transient/reference for time/energy accounting before a core transient. Physical photon or jitter remains excluded.']
audit['cost_estimate']={'measured_no_FD_561_nodes_seconds':[517.7688461309299,517.966562661808],'no_FD_2145_nodes_16x8_degree4_rough_seconds':float(517.87*2145/561),'estimate_qualification':'Linear direct-spectrum query estimate only; actual cost can change with cache and spectral support. User-run command, not launched by audit.'}
(O/'independent_audit.json').write_text(json.dumps(audit,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf8')
print(json.dumps({'identity':audit['identity'],'mesh_comparisons':comparisons,'spatial_identity_residuals':audit['spatial_identity_residuals'],'center_results':center_results},indent=2))
