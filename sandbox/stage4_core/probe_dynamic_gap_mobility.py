"""Uniform moving-gap coefficient with kinetic response, not a damping fit.

This bounded analytic quadrature compares holding h=tanh(E/2T) artificially
fixed with the collisionless adiabatic response determined by the longitudinal
kinetic equation. Eta is only a regulator. No physical relaxation is inferred.
"""
from __future__ import annotations
import os
for _key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_key, '1')
import argparse
import json
from pathlib import Path
import time
import numpy as np
from numpy.polynomial.legendre import leggauss


def evaluate(gap, temperature, eta, order=64):
    """Return coefficient G1/(m*kappa*d_dot), positive-energy convention."""
    edges = sorted(set([0., temperature, 3*temperature, 8*temperature,
        max(0., gap-10*eta), max(0., gap-eta), gap, gap+eta,
        gap+10*eta, max(8*gap, gap+20*eta)]))
    nodes, weights = leggauss(order)
    energy = np.concatenate([(left+right)/2+(right-left)*nodes/2 for left, right in zip(edges[:-1], edges[1:])])
    weight = np.concatenate([(right-left)*weights/2 for left, right in zip(edges[:-1], edges[1:])])
    z = eta-1j*energy
    root = np.sqrt(z*z+gap*gap)
    g, f = z/root, gap/root
    f_gap = z*z/(root**3)
    f_energy = 1j*gap*z/(root**3)
    h_energy = .5/temperature/np.cosh(energy/(2*temperature))**2
    # From rho*hLdot+Im(f)*d_dot*hL_E=0, with d_dot=1.
    h_time = -f.imag/g.real*h_energy
    held_integrand = 2*f_gap.real*h_energy
    moving_integrand = 2*(f_gap.real*h_energy-f_energy.real*h_time)
    temporal_residual = g.real*h_time+f.imag*h_energy
    return {'eta': float(eta), 'eta_over_gap': float(eta/gap), 'order': int(order),
        'energy_nodes': int(len(energy)),
        'held_distribution_coefficient': float(weight@held_integrand),
        'kinetically_evolved_coefficient': float(weight@moving_integrand),
        'absolute_integral_evolved': float(weight@abs(moving_integrand)),
        'maximum_kinetic_residual': float(np.max(abs(temporal_residual))),
        'minimum_DOS': float(np.min(g.real))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--gap', type=float, default=1.76392592757718)
    parser.add_argument('--temperature-ratio', type=float, default=.9/8.65)
    args = parser.parse_args()
    start = time.monotonic()
    records = [evaluate(args.gap, args.temperature_ratio, fraction*args.gap, order)
               for fraction in (.1, .05, .02, .01, .005, .002) for order in (32, 64)]
    result = {'schema': 'stage4_uniform_moving_gap_coefficient_v1',
        'scope': 'Uniform real-gap collisionless adiabatic response; no circuit, no spatial bias, no physical damping calibration',
        'normalization': 'G_radial^(1)/(node_area*kappa*d_gap_dt)',
        'gap': args.gap, 'temperature_over_Tc': args.temperature_ratio,
        'records': records, 'runtime_seconds': time.monotonic()-start,
        'conclusion': 'The held-h coefficient is not an intrinsic mobility. Use the response satisfying the kinetic equation; regulator dependence cannot define a physical relaxation time.'}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
