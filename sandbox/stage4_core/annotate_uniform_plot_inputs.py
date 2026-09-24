"""Add derived plotting metadata without changing archived numerical records."""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.constants import hbar, k as KB


def decode(value):
    return np.asarray(value['real'])+1j*np.asarray(value['imag'])


def encode(value):
    value=np.asarray(value)
    return dict(real=value.real.tolist(),imag=value.imag.tolist())


def annotate(case, output):
    case,output=Path(case),Path(output)
    source=case/'results.json'; manifest_path=case/'manifest.json'
    raw=json.loads(source.read_text(encoding='utf8'))
    manifest=json.loads(manifest_path.read_text(encoding='utf8'))
    derived=copy.deepcopy(raw)
    E0=KB*manifest['plan']['Tc_K']; d=manifest['d']; t=manifest['plan']['T_K']/manifest['plan']['Tc_K']
    gx,gw=leggauss(24)
    errors=[]
    for row in derived['records']:
        if row['kind']!='port':continue
        # CM.4 independently reconstructs the missing recorded source amplitude
        # from the SAVED exp(-iwt) three-state response. No numerical model rerun.
        ib,detector,vc=decode(row['circuit_state_peak'])
        omega=row['omega']*E0/hbar
        drive=(1e4+50.-1j*omega*1e-6)*ib-50.*detector+vc
        error=abs(drive-1e-6)
        if error>1e-17:
            raise ValueError('Saved circuit state does not recover the documented 1 microvolt source')
        errors.append(error)
        row['bias_voltage_peak_V']=1e-6
        row['time_convention']='exp(-iwt)'
        row['derived_bias_reconstruction_peak_V']=encode(drive)
        center=row['energy_cutoff']; Omega=row['omega']; eta=row['eta']
        energy=center+Omega*gx/2
        z=eta-1j*energy
        root=np.sqrt(z*z+d*d);root=np.where(root.real<0,-root,root)
        rho=(z/root).real
        endpoint=1-.5*np.dot(gw,np.tanh(energy/(2*t))*rho)
        row['neutral_finite_cutoff_target']=encode(complex(endpoint))
    derived['derived_metadata']=dict(
        raw_results_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        raw_manifest_sha256=hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
        executed_source_sha256=manifest['sources']['sandbox/stage4_core/coupled_response.py'],
        annotation_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        source_amplitude_basis='Known hardcoded 1e-6 V call, independently recovered by CM.4 from every stored circuit state',
        maximum_reconstructed_drive_error_V=max(errors),
        extra_inputs=dict(Rb_ohm=1e4,Lb_H=1e-6,RL_ohm=50.,provenance='CM thesis defaults used by the executed circuit call'),
        gauge_tail='Analytic uniform BCS endpoint integral at the recorded eta, frequency and cutoff; no fitted normalization',
        original_numeric_fields_changed=False)
    output.mkdir(parents=True,exist_ok=True)
    target=output/'results_for_plot.json'
    target.write_text(json.dumps(derived,indent=2,allow_nan=False)+'\n',encoding='utf8')
    (output/'derived_metadata.json').write_text(json.dumps(derived['derived_metadata'],indent=2)+'\n',encoding='utf8')
    return target


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    print(annotate(args.case,args.output))
