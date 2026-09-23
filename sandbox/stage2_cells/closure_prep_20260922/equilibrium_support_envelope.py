"""Conservative analytic envelope: no new quadrature or dynamics.

For the causal retarded solution, c=u/sqrt(u^2-a^2),
u=E+i eta+i Gamma c. At E>a, Re(u)>=E and
|c|^2<=|u|^2/(|u|^2-a^2)<=E^2/(E^2-a^2).
The scattering factor NN'-RR'<=NN' is therefore bounded by the
BCS bound at the lowest tail energy, even at finite eta and Gamma.
The remaining Debye polynomial is integrated analytically. Maxima of the
actual piecewise-linear populations bound the geometric interpolation.
"""
from pathlib import Path
import argparse,hashlib,json,sys,time
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT),str(ROOT/'sandbox/stage2_cells')]
import numpy as np
from run_coupled import setup
from pysnspd.experimental.cell_validation import ElectronicCell
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('trajectory',type=Path);parser.add_argument('--quadrature-record',type=Path,required=True)
    parser.add_argument('--registration',type=Path,required=True);parser.add_argument('--output',type=Path)
    parser.add_argument('--prepare',action='store_true');args=parser.parse_args()
    record=json.loads(args.trajectory.read_text());quad=json.loads(args.quadrature_record.read_text())
    assert quad['trajectory_sha256']==sha(args.trajectory)
    assert record['trajectory_sha256']==sha(args.trajectory.with_suffix('.npz'))
    registration=dict(schema='pysnspd.stage2.equilibrium-envelope-registration.v1',
        source_sha256=sha(__file__),trajectory_sha256=sha(args.trajectory),
        archive_sha256=record['trajectory_sha256'],quadrature_record_sha256=sha(args.quadrature_record),
        physical_sources=record['source_hashes'],relative_omission_limit=1e-3,
        method='Analytic causal spectral envelope, maxima of actual interpolation populations, exact Debye polynomial integral.',
        reason='The prior16/32quadrature bound was tiny but did not converge relative to itself; it is retained as FAIL. This broader exact envelope has no quadrature error and does not change tolerances.',
        reference_note='Uses the same gross in-support activity measured in the preserved quadrature record only as the normalization denominator.',
        no_new_RHS=True,no_new_trajectories=True)
    if args.prepare:
        with args.registration.open('x',encoding='utf-8') as f:json.dump(registration,f,indent=2)
        print('ANALYTIC_ENVELOPE_REGISTERED');return
    assert json.loads(args.registration.read_text())==registration
    for name,digest in record['source_hashes'].items():assert sha(ROOT/name)==digest
    assert args.output is not None and not args.output.exists()
    started=time.perf_counter();par=record['parameters']
    system,_,debye=setup(par['case'],par['phonon_nodes'],par['infrared'],par['face_order'],float(par['escape']),par['heating'],par['scenario'],par['reaction_order'],par['reaction_method'],par['electron_refinement'],par['reaction_layout'],par['reaction_outer_order'],par['reaction_max_panel'],par['reaction_max_energy_panel'])
    with np.load(args.trajectory.with_suffix('.npz'),allow_pickle=False) as file:times,states=file['times'],file['states']
    rows=[];cut=debye.cutoff_energy_bar;lo=debye.infrared_cutoff_bar
    for time_value,state in zip(times,states):
        amplitudes,p,n=system.unpack(state)
        for i,a in enumerate(amplitudes):
            cell=ElectronicCell(system.catalog,a,system.gammas[i]);emax=float(cell.energies[-1]);emin=emax-cut
            if not emin>a:raise ValueError('The envelope requires all tail energies strictly above the gap scale')
            pmax=float(max(np.interp(emin,cell.energies,p[i]),np.max(p[i][cell.energies>=emin])))
            nmax=float(np.max(n[i]));coherence=emin*emin/(emin*emin-a*a)
            factor=system.rate_prefactor*coherence*pmax*nmax*debye.lambda_eph/cut**2
            bounds=dict(number=factor*(cut**4-lo**4)/4,power=factor*(cut**5-lo**5)/5)
            saved=next(r for r in quad['rows'] if r['time']==float(time_value) and r['cell']==i)
            denom=saved['gross_activity']
            relative={k:bounds[k]/denom[k] for k in bounds}
            if any(denom[k]<=0 for k in bounds):raise ValueError('Undefined gross normalization')
            if any(bounds[k]<system.rate_prefactor*saved['bounds'][-1][k] for k in bounds):
                raise ValueError('Envelope does not dominate the sampled quadrature')
            rows.append(dict(time=float(time_value),cell=i,gamma=system.gammas[i],amplitude=float(a),
                minimum_tail_energy=emin,maximum_linear_electron_population=pmax,
                maximum_linear_phonon_population=nmax,coherence_upper_bound=coherence,
                absolute_rate_upper_bounds=bounds,relative_upper_bounds=relative))
    worst=max(value for row in rows for value in row['relative_upper_bounds'].values())
    result=dict(status='PASS' if worst<=1e-3 else 'FAIL',worst_relative_upper_bound=worst,
        relative_limit=1e-3,samples=len(rows),distinct_times_per_cell=len(np.unique(times)),
        proof='docs/implementation/stage2/closure_prep_20260922/EQUILIBRIUM_SUPPORT_ENVELOPE.md',
        quadrature_refinement_required=False,reason='All integrations in the envelope are exact polynomial primitives; population maxima at native knots/endpoints are exact for the declared linear interpolants.',
        finite_eta_included_in_causal_bound=True,rows=rows,registration_sha256=sha(args.registration),
        source_sha256=sha(__file__),runtime_seconds=time.perf_counter()-started,new_trajectories=0,
        scope='Absorption to unresolved high-energy electronic states; no arbitrary incoming external population is assumed or bounded.',stage2_admission=False)
    with args.output.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2,allow_nan=False)
    print(json.dumps({k:result[k] for k in ('status','worst_relative_upper_bound','samples','runtime_seconds')}))
    if result['status']!='PASS':raise SystemExit(1)

if __name__=='__main__':main()
