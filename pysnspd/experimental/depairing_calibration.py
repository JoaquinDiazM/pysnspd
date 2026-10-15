"""Provisional diffusion calibration against a stated homogeneous depairing target.

This utility never changes production defaults or an existing geometry. Its tiny
ThermalGraph is only a constructor placeholder for continuum ``bulk`` queries;
it is not a two-node device calculation. A future functional must be supplied
explicitly and re-evaluated before claiming calibration of that closure.
"""
from __future__ import annotations

import math
from typing import Callable

import numpy as np
from scipy.constants import Boltzmann as kB, elementary_charge as e, hbar
from scipy.optimize import brentq

from .effective_bcs_functional import EffectiveBCSFunctional, bcs_count_quadrature
from .thermal_spatial_usadel import ThermalGraph

K20_SOURCE = "https://eprints.lancs.ac.uk/id/eprint/140252/3/Binder1.pdf"
K20_DOI = "https://doi.org/10.1038/s41566-020-0589-x"
K20_REFERENCE_WIDTH_M = 120e-9
K20_REFERENCE_DEPAIRING_CURRENT_A = 38.8e-6


def _real(value, name, *, zero=False):
    if isinstance(value, (bool, np.bool_)) or np.iscomplexobj(value):
        raise ValueError(f"{name} must be a real finite number")
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be a real finite number") from exc
    if not math.isfinite(result) or (result < 0 if zero else result <= 0):
        raise ValueError(f"{name} must be {'nonnegative' if zero else 'positive'} and finite")
    return result


def calibration_result(*, sheet_resistance_ohm=425., width_m=80e-9,
                       target_depairing_current_A=None, requested_current_A=None,
                       thickness_m=7e-9, critical_temperature_K=8.65,
                       bath_temperature_K=.9, reference_diffusion_m2_s=5e-5,
                       count_cutoff_bar=32., count_order=6,
                       bulk_factory: Callable | None = None):
    """Return JSON-ready calibrated SI parameters and their scientific scope.

    The default target scales K20's 38.8 uA at 120 nm linearly with width;
    that extrapolation is an explicit inference, not an 80 nm measurement.
    At fixed sheet resistance the native closure has I_dep proportional to
    1/(Rs*sqrt(D)), N0=1/(2*e**2*D*Rs*d), ell0=sqrt(hbar*D/(2*kB*Tc)).
    The inferred D is checked by rebuilding and evaluating the functional.

    ``bulk_factory`` accepts the same keyword arguments as
    ``EffectiveBCSFunctional`` (including ``graph``) and returns an object with
    ``bulk(temperature_K=..., current_A=0., width_m=...)``. It is called afresh
    for each D, with updated N0 and ell0. If its continuum scaling differs,
    a bracketed root recalibrates D against that actual callback. Only bulk DC
    is checked; this does not validate geometry, dynamics or measured latency.
    """
    rs = _real(sheet_resistance_ohm, "sheet_resistance_ohm")
    width = _real(width_m, "width_m")
    thickness = _real(thickness_m, "thickness_m")
    tc = _real(critical_temperature_K, "critical_temperature_K")
    tb = _real(bath_temperature_K, "bath_temperature_K")
    dref = _real(reference_diffusion_m2_s, "reference_diffusion_m2_s")
    if tb >= tc:
        raise ValueError("bath_temperature_K must be below critical_temperature_K")
    inferred_target = target_depairing_current_A is None
    target = (K20_REFERENCE_DEPAIRING_CURRENT_A*width/K20_REFERENCE_WIDTH_M
              if inferred_target else _real(target_depairing_current_A,
                                             "target_depairing_current_A"))
    requested = (None if requested_current_A is None else
                 _real(requested_current_A, "requested_current_A", zero=True))
    if requested is not None and (requested >= target or math.isclose(requested, target, rel_tol=8*np.finfo(float).eps, abs_tol=0.)):
        raise ValueError("requested_current_A must be strictly below the depairing target")
    xi, weights = bcs_count_quadrature(count_cutoff_bar, count_order)
    factory = EffectiveBCSFunctional if bulk_factory is None else bulk_factory
    if not callable(factory):
        raise ValueError("bulk_factory must be callable")
    samples = {}

    def sample(diffusion):
        diffusion = _real(diffusion, "trial diffusion_m2_s")
        if diffusion not in samples:
            n0 = 1/(2*e**2*diffusion*rs*thickness)
            ell = math.sqrt(hbar*diffusion/(2*kB*tc))
            graph = ThermalGraph(np.array([.5, .5]), np.array([[0, 1]]),
                                 np.array([1.]), np.array([[0., 0.], [1., 0.]]))
            model = factory(graph=graph, energy_scale_J=kB*tc,
                            N0_per_J_m3=n0, ell0_m=ell, thickness_m=thickness,
                            count_nodes=xi.copy(), count_weights=weights.copy(),
                            gradient_law="bcs_radial", bath_temperature_K=tb,
                            diffusion_m2_s=diffusion)
            result = model.bulk(temperature_K=tb, current_A=0., width_m=width)
            current = _real(result.depairing_current_A, "bulk depairing_current_A")
            samples[diffusion] = (current, n0, ell)
        return samples[diffusion]

    reference, n0ref, ellref = sample(dref)
    guess = _real(dref*(reference/target)**2, "scaling-estimated diffusion_m2_s")
    achieved, n0, ell = sample(guess)
    calibrated = guess
    method = "analytic continuum scaling, independently re-evaluated bulk"
    if not math.isclose(achieved, target, rel_tol=2e-8, abs_tol=0.):
        if bulk_factory is None:
            raise ArithmeticError("native bulk violated the stated D scaling")
        center = math.log(guess)
        for power in range(1, 6):
            radius = power*math.log(4.)
            lo, hi = center-radius, center+radius
            if (sample(math.exp(lo))[0]-target)*(sample(math.exp(hi))[0]-target) <= 0:
                break
        else:
            raise ValueError("custom bulk_factory depairing target could not be bracketed")
        root = brentq(lambda value: math.log(sample(math.exp(value))[0]/target),
                      lo, hi, xtol=2e-11, rtol=2e-11, maxiter=32)
        calibrated = math.exp(root)
        achieved, n0, ell = sample(calibrated)
        method = "bracketed log-D root of supplied continuum bulk_factory"
    residual = (achieved-target)/target
    if not math.isfinite(residual) or abs(residual) > 2e-8:
        raise ArithmeticError("calibrated bulk does not reproduce depairing target")
    if requested is not None and requested >= achieved:
        raise ValueError("requested_current_A must be strictly below achieved depairing current")
    rs_status = ("conditional inference 340 ohm/square divided by RRR 0.8"
                 if rs == 425. else "K20 historical simulation parameter"
                 if rs == 608. else "explicit supplied sheet-resistance hypothesis")
    return dict(
        schema="pysnspd.depairing_calibration.v1",
        diffusion_m2_s=calibrated, N0_per_J_m3=n0, ell0_m=ell,
        sheet_resistance_ohm=rs, width_m=width, thickness_m=thickness,
        critical_temperature_K=tc, bath_temperature_K=tb,
        target_depairing_current_A=target, achieved_depairing_current_A=achieved,
        relative_target_residual=residual, requested_current_A=requested,
        requested_current_fraction=None if requested is None else requested/achieved,
        requested_current_below_depairing=None if requested is None else True,
        reference=dict(diffusion_m2_s=dref, depairing_current_A=reference,
                       N0_per_J_m3=n0ref, ell0_m=ellref),
        quadrature=dict(cutoff_bar=float(count_cutoff_bar), order=count_order,
                        electronic_states=len(xi), finite_native_BCS=True),
        calibration=dict(method=method, bulk_evaluations=len(samples),
                         factory=getattr(factory, "__module__", type(factory).__module__)
                         +"."+getattr(factory, "__qualname__", type(factory).__qualname__),
                         continuum=True, finite_phase_step_m=None,
                         target_is_width_scaled_inference=inferred_target,
                         diffusion_is_fitted_effective_parameter=True,
                         sheet_resistance_status=rs_status,
                         physical_admission=False,
                         status="PROVISIONAL_CONTINUUM_BULK_CALIBRATION",
                         requires_verification_with_new_closure=bulk_factory is None,
                         closure_supplied_explicitly=bulk_factory is not None),
        geometry=dict(bulk_constructor_nodes=2, device_mesh_validated=False,
                      rebuild_physical_geometry_with_calibrated_ell0=True,
                      reuse_old_normalized_mesh_or_ell0=False),
        assumptions=dict(
            scaling="Idep proportional to width/(Rs*sqrt(D)) at fixed Tc,Tb and finite quadrature",
            N0_convention="single-spin normal DOS inferred with sigma=2*e^2*N0*D=1/(Rs*d)",
            target="38.8 uA at 120 nm, linearly scaled to requested width" if inferred_target
                   else "explicit caller-supplied target; K20 width inference not imposed",
            future_closure="A thermal stiffness sharing the same bath-DC functional may retain this calibration; verify actual new closure",
            sigma_nm="independent preparation input; never recomputed by this utility",
            scope="Homogeneous continuum bulk only; no finite mesh, photon or latency admission"),
        sources=dict(
            K20=dict(url=K20_SOURCE, doi=K20_DOI, page_convention="one-based PDF pages",
                     depairing_locator="PDF13, Methods: Depairing current estimation; 120 nm, 38.8 uA",
                     sheet_resistance_locator="PDF9-10, Nanofabrication: 340 ohm/square at room temperature and RRR0.8; 425 is conditional inference",
                     simulation_locator="PDF33, Supplementary Table1: Rs608 ohm/square and D0.5 cm^2/s"),
            closure="pysnspd.experimental.effective_bcs_functional: finite bcs_count_quadrature and continuum bulk()"))