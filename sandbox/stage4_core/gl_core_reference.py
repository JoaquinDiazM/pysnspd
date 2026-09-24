"""Small static thermal-core correspondence check against the GL tanh solution.

Independent of experimental spatial assembly and retarded spectral catalogue.
No device transient, vortex solver, stochastic barrier or physical run is used.
"""
from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq
from scipy.special import zeta


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT/'docs/implementation/stage4/followup_20260923/gl_core_reference.json'
KAPPA = np.pi/4


def tail(power, count, t):
    return float(zeta(power, count+.5)/(2*np.pi*t)**power)


def gap_bracket(d, t, count):
    values = np.asarray(d, dtype=float)
    eps = 2*np.pi*t*(np.arange(count)+.5)
    a = values[..., None]**2
    energy = np.sqrt(eps**2+a)
    summed = np.sum(a/(eps*energy*(eps+energy)), axis=-1)
    aa = values**2
    summed += aa/2*tail(3,count,t)-3*aa**2/8*tail(5,count,t)+5*aa**3/16*tail(7,count,t)
    return np.log(t)+2*np.pi*t*summed


def bcs_potential(d, t, count):
    """V(d)=ln(t)d^2+2*pi*t sum[d^2/eps-2*(E-eps)], stably evaluated."""
    values = np.asarray(d, dtype=float)
    eps = 2*np.pi*t*(np.arange(count)+.5)
    a = values[..., None]**2
    energy = np.sqrt(eps**2+a)
    summed = np.sum(a*a/(eps*(energy+eps)**2), axis=-1)
    aa = values**2
    summed += aa**2/4*tail(3,count,t)-aa**3/8*tail(5,count,t)+5*aa**4/64*tail(7,count,t)
    return np.log(t)*values**2+2*np.pi*t*summed


def evaluate(t, *, count, order, include_profiles):
    alpha = -math.log(t)
    beta = float(7*zeta(3, 1)/(16*np.pi**2*t*t))
    d_gl = math.sqrt(alpha/(2*beta))
    d_bcs = float(brentq(lambda d: float(gap_bracket(d,t,count)), .001, 2.))
    xi_gl = math.sqrt(KAPPA/alpha)
    length = math.sqrt(2)*xi_gl
    condensation_gl = alpha**2/(4*beta)
    e_gl_exact = 8*math.sqrt(2)/3*condensation_gl*xi_gl
    points, weights = leggauss(order)
    y = 8*points
    weights = weights*8*length
    tanh = np.tanh(y)
    sech2 = 1-tanh*tanh
    def profile(gap):
        return gap*tanh, gap/length*sech2, -2*gap/length**2*sech2*tanh
    gl, gl_prime, gl_second = profile(d_gl)
    bcs, bcs_prime, bcs_second = profile(d_bcs)
    # Analytic positive GL excess prevents cancellation near the endpoints.
    e_gl_quad = float(np.dot(weights, KAPPA*gl_prime**2+beta*(gl**2-d_gl**2)**2))
    bcs_excess = bcs_potential(bcs,t,count)-float(bcs_potential(d_bcs,t,count))
    e_bcs = float(np.dot(weights, KAPPA*bcs_prime**2+bcs_excess))
    gl_force = -2*alpha*gl+4*beta*gl**3-2*KAPPA*gl_second
    bcs_force = 2*bcs*gap_bracket(bcs,t,count)-2*KAPPA*bcs_second
    force_scale = 2*alpha*d_bcs
    row = dict(T_over_Tc=t, matsubara_count=count, spatial_quadrature_order=order,
        dimensionless_half_domain_in_profile_widths=8., kappa=KAPPA,
        alpha=alpha, quartic_beta=beta, gap_GL_kBTc=d_gl, gap_BCS_kBTc=d_bcs,
        gap_BCS_over_GL_minus_one=d_bcs/d_gl-1,
        xi_GL_over_ell0=xi_gl, profile_width_over_ell0=length,
        condensation_energy_GL_bar=condensation_gl,
        core_energy_GL_exact_bar=e_gl_exact, core_energy_GL_quadrature_bar=e_gl_quad,
        core_energy_BCS_prescribed_bar=e_bcs,
        prescribed_energy_BCS_over_GL_minus_one=e_bcs/e_gl_exact-1,
        GL_energy_quadrature_relative_defect=abs(e_gl_quad/e_gl_exact-1),
        force_scale=force_scale,
        GL_maximum_force_on_quadrature=float(np.max(abs(gl_force))),
        BCS_residual_force_L2_scaled=float(np.sqrt(np.dot(weights,bcs_force*bcs_force)/(2*8*length))/force_scale),
        BCS_gap_equation_residual=float(gap_bracket(d_bcs,t,count)))
    if include_profiles:
        yy = np.linspace(-4,4,161)
        th = np.tanh(yy); s2 = 1-th*th
        zz = d_bcs*th; zpp = -2*d_bcs/length**2*s2*th
        f = 2*zz*gap_bracket(zz,t,count)-2*KAPPA*zpp
        row['profiles'] = dict(x_over_profile_width=yy.tolist(),
            GL_amplitude_kBTc=(d_gl*th).tolist(), BCS_prescribed_amplitude_kBTc=zz.tolist(),
            BCS_force_scaled=(f/force_scale).tolist(), GL_force_scaled=np.zeros_like(yy).tolist())
        row['BCS_maximum_sampled_residual_force_scaled'] = float(np.max(abs(f))/force_scale)
    return row


def calculate():
    started = time.monotonic()
    temperatures = (.95,.98,.99)
    coarse = [evaluate(t,count=256,order=160,include_profiles=False) for t in temperatures]
    fine = [evaluate(t,count=512,order=320,include_profiles=True) for t in temperatures]
    return dict(schema='pysnspd.stage4.gl_core_reference.v1',
        status='COMPUTED_STATIC_THERMAL_LIMIT_NOT_KORZH_CORE_ADMISSION',
        references=[
            'https://doi.org/10.1103/PhysRev.164.498',
            'https://arxiv.org/html/2506.18130v1',
            'docs/modelo_v0_4/C_condensado_corriente_y_senal_v0_4.md C.9 and C.29-C.31'],
        physical_time_steps=0, production_changed=False, Korzh_reference_temperature_changed=False,
        independent_of_experimental_spatial_assembly_and_retarded_catalogue=True,
        ensemble='Prescribed temperature, equilibrium Matsubara free energy, zero current; not frozen counts',
        coordinates='x/ell0, ell0=sqrt(hbar*D/(2*kB*Tc)); d=Delta/(kB*Tc) is signed real',
        energy_unit='N0*(kB*Tc)^2*A*ell0 for a longitudinal profile uniform across area A',
        candidate='Integral[kappa*(d_prime)^2 + V_BCS(d)-V_BCS(d_BCS)] with kappa=pi/4',
        gl_reference='Quartic thermal expansion with alpha=-ln(t), beta=7*zeta(3)/(16*pi^2*t^2); same kappa',
        profile_convention='d(x)=d_bulk*tanh[x/(sqrt(2)*xi_GL)]; each theory uses its own equilibrium bulk amplitude',
        candidate_profile_is_stationary_solution=False,
        zero_current_GL_profile_is_exact_stationary_solution=True,
        finite_domain='Integrate x/profile_width in [-8,8], retain whole-line analytic GL reference',
        tail='Gap inverse frequency powers3,5,7; potential powers3,5,7; analytic Hurwitz-zeta correction',
        checks=dict(coarse_count=256,fine_count=512,coarse_quadrature=160,fine_quadrature=320,
            maximum_BCS_energy_relative_difference=max(abs(a['core_energy_BCS_prescribed_bar']/b['core_energy_BCS_prescribed_bar']-1) for a,b in zip(coarse,fine)),
            maximum_GL_exact_energy_relative_defect=max(r['GL_energy_quadrature_relative_defect'] for r in fine),
            interpretation='One combined cutoff/quadrature repeat; numerical stability, not a physical uncertainty estimate'),
        results=fine, coarse_results=coarse,
        restrictions=[
            'Temperatures are separate mathematical correspondence controls, not new device operating points',
            'No finite-current or low-temperature core/boundary/vortex admission',
            'Real field implies q_delta=0; this does not validate the phase regularizer',
            'Candidate profile is prescribed, not a computed BCS saddle or barrier',
            'GL leading behavior is checked, not arbitrary-temperature microscopic gradient physics',
            'No admission of1D reduction for the80nm strip or prediction of dark-count rate',
            'No KWT calibration or electron-phonon heat partition validation'],
        runtime_seconds=time.monotonic()-started)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=OUTPUT)
    args=parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Preserve previous evidence; choose a new output path')
    result=calculate()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf8')
    print(json.dumps(dict(output=str(args.output),runtime_seconds=result['runtime_seconds'],
        checks=result['checks'],results=[{k:v for k,v in r.items() if k!='profiles'} for r in result['results']]),indent=2))


if __name__=='__main__':
    main()
