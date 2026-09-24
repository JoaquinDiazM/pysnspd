"""Compact evidence and static-overlap check for the frozen charge resolvent."""
from pathlib import Path
import argparse,hashlib,json
import numpy as np


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text())


def analyze(source,moments,output):
    if output.exists():raise FileExistsError('Fresh compact evidence directory required')
    summary=read(source/'summary.json');records=summary['records'];prior=read(moments/'projection/summary.json')['records']
    comparisons=read(source/'comparisons.json');checks=[];phase_rows=[];certificate=[];compact={}
    wardmax=dynmax=solvermax=0.
    for row in records:
        path=source/row['fields_path']
        if sha(path)!=row['fields_sha256']:raise ValueError('Changed response array')
        certificate.append(dict(path=row['fields_path'],sha256=row['fields_sha256'],bytes=path.stat().st_size))
        for probe in row['probes'].values():
            for name in ('full','projected'):
                wardmax=max(wardmax,probe['norms'][name]['ward_max'])
                dynmax=max(dynmax,probe['norms'][name]['integrated_dynamic_residual_free_max'])
            solvermax=max(solvermax,max(p['full_equation_free_max'] for p in probe['per_energy']))
        static=next(r for r in records if r['case_id']==row['case_id'] and r['eta_relative']==row['eta_relative'] and r['nu']==0)
        with np.load(path) as a,np.load(source/static['fields_path']) as b:
            xy,m,c,edges=a['coordinates_bar'],a['area_weights'],a['conductance'],a['edges']
            core=np.linalg.norm(xy,axis=1)<=4;ecore=np.linalg.norm((xy[edges[:,0]]+xy[edges[:,1]])/2,axis=1)<=4
            for probe in ('radial','angular'):
                for key in ('current','energy_flux','phase_torque_density'):
                    old=b[probe+'_full_'+key];new=a[probe+'_full_'+key]
                    mask,weight=(ecore,1/c) if key!='phase_torque_density' else (core,m)
                    denom=np.sum(weight[mask]*abs(old[mask])**2)
                    gain=np.sum(weight[mask]*old[mask].conj()*new[mask])/denom if denom>1e-24 else None
                    phase_rows.append(dict(case_id=row['case_id'],eta_relative=row['eta_relative'],nu=row['nu'],probe=probe,moment=key,
                        projected_gain_real=None if gain is None else float(gain.real),
                        projected_gain_imag=None if gain is None else float(gain.imag),
                        phase_degrees=None if gain is None else float(np.angle(gain,deg=True))))
            if row['nu']==0:
                grid='fine50' if row['energy_nodes']==50 else 'base31'
                oldrow=next(r for r in prior if r['case_id']==row['case_id'] and r['eta_relative']==row['eta_relative'] and r['grid']==grid)
                with np.load(moments/'projection'/oldrow['fields_path']) as old:
                    differences=[]
                    for probe in ('radial','angular'):
                        for name in ('zero_hT','full','projected'):
                            prefix=probe+'_'+name+'_'
                            for key in ('current','energy_flux','phase_torque_density','amplitude_force_density'):
                                differences.append(float(np.max(abs(a[prefix+key]-old[prefix+key]))))
                            force=a[prefix+'force_cartesian'];previous=old[prefix+'force']
                            differences.extend([float(np.max(abs(force[:,0]-previous.real))),float(np.max(abs(force[:,1]-previous.imag)))])
                    checks.append(dict(case_id=row['case_id'],eta_relative=row['eta_relative'],
                        baseline_fields=oldrow['fields_path'],baseline_sha256=oldrow['fields_sha256'],maximum_absolute_difference=max(differences)))
            if row['case_id']=='radial_65_N256' and row['eta_relative']==.01 and row['nu'] in (0.,.001,.1):
                prefix='nu'+str(row['nu']).replace('.','p')+'_'
                for name in ('full','projected'):
                    for key in ('force_cartesian','current','energy_flux','phase_torque_density','amplitude_force_density'):
                        compact[prefix+name+'_'+key]=a['angular_'+name+'_'+key]
                for key in ('coordinates_bar','d','area_weights','edges','conductance','boundary_nodes'):compact[key]=a[key]
    result=dict(status='CHARGE_FREQUENCY_EVIDENCE_COMPLETE',runtime_seconds=summary['runtime_seconds'],jobs=len(records),
        maximum_ward_residual=wardmax,maximum_integrated_dynamic_charge_residual=dynmax,
        maximum_per_energy_charge_solve_residual=solvermax,static_overlap=checks,phasor_alignment=phase_rows,
        full_array_certificate=certificate,full_array_bytes=sum(r['bytes'] for r in certificate),
        slow_anchor_ranges=[],raw_remote_directory=str(source),
        interpretation='Fixed-spectrum leading-adiabatic charge resolvent. Slow anchors test static charge elimination; high frequencies are exploratory, not a detector bandwidth certificate.')
    for nu in (.0001,.001,.01,.1,1.):
        for moment in ('current','phase_torque_density','energy_flux'):
            vals=[r['relative_dynamic_change'] for r in comparisons if r['nu']==nu and r['probe']=='angular' and r['moment']==moment]
            result['slow_anchor_ranges'].append(dict(nu=nu,moment=moment,minimum=min(vals),maximum=max(vals)))
    output.mkdir(parents=True)
    (output/'analysis.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    np.savez_compressed(output/'representative_angular_phasors.npz',**compact)
    (output/'compact_identity.json').write_text(json.dumps(dict(analysis_source_sha256=sha(Path(__file__)),
        summary_sha256=sha(source/'summary.json'),comparisons_sha256=sha(source/'comparisons.json'),
        compact_sha256=sha(output/'representative_angular_phasors.npz'),
        compact_scope='Angular radial65 eta.01 at nu=0,.001,.1; all24 full maps remain certified on Geminga',
        static_prior_projection_summary_sha256=sha(moments/'projection/summary.json')),indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('full_array_certificate','phasor_alignment')},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--moments',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();analyze(a.source,a.moments,a.output)
