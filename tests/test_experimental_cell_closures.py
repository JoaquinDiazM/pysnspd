"""Independent D.11--13 and Debye identities, not detector calibration."""
from dataclasses import replace
import hashlib
import json
import numpy as np
import pytest
from scipy.integrate import quad

from pysnspd.experimental.cell_closures import (
    CellScales, KWTMobility, DebyePhonons, ConditionalPhononShape,
    PhononQuadrature, bose_occupation,
)

HBAR=1.054571817e-34
KB=1.380649e-23
SCALES=CellScales(2.106e-22,1.1e47,8.65,1.83,1.0)


def reference_SI(z,force,theta,scales):
    """Solve the dimensional Cartesian equation directly, then convert once."""
    delta=scales.delta0_J*np.asarray(z)
    tmax=max(theta*scales.delta0_J/KB,scales.Tb_K)
    tau=1e-12/((tmax/scales.Tc_K)/.50+(tmax/scales.Tc_K)**3/2.47)
    tau0=np.pi*HBAR/(8*KB*scales.Tc_K)
    A0=scales.N0_per_J_m3*np.sqrt((1+tmax/scales.Tc_K)/2)
    R=np.sqrt(1+4*np.dot(delta,delta)*tau*tau/HBAR**2)
    temporal=A0*tau0/R*(np.eye(2)+4*tau*tau/HBAR**2*np.outer(delta,delta))
    dimensional_force=scales.N0_per_J_m3*scales.delta0_J*np.asarray(force)/2
    velocity=np.linalg.solve(temporal,-dimensional_force)
    rho_dot=2*np.dot(delta,velocity)
    heat=2*A0*tau0/R*(np.dot(velocity,velocity)+tau*tau/HBAR**2*rho_dot*rho_dot)
    return velocity*scales.t_ref_s/scales.delta0_J, heat*scales.t_ref_s/scales.energy_density_scale_J_m3


@pytest.mark.parametrize("amplitude",[0.,.2,.8,1.5])
@pytest.mark.parametrize("temperature",[0.,.12,.6,1.2])
def test_cartesian_equation_and_heat_match_independent_SI(amplitude,temperature):
    z=np.array([.6,.8])*amplitude; force=np.array([.7,-.45])
    response=KWTMobility(SCALES).tensor_response(z,force,temperature)
    velocity,heat=reference_SI(z,force,temperature,SCALES)
    np.testing.assert_allclose(response.velocity,velocity,rtol=1e-11,atol=1e-13)
    assert response.heat==pytest.approx(heat,rel=1e-11,abs=1e-13)
    assert response.heat==pytest.approx(-np.dot(force,response.velocity),rel=1e-11)
    assert response.heat>=0 and np.min(np.linalg.eigvalsh(response.mobility))>0


def test_radial_factor_two_tangential_eigenvalue_and_zero_continuity():
    closure=KWTMobility(SCALES)
    for amplitude in [0.,1e-10,.35,1.3]:
        radial=closure.amplitude_response(amplitude,.8,.3)
        vector=closure.tensor_response([amplitude,0],[.8,0],.3)
        np.testing.assert_allclose(vector.velocity,[radial.velocity,0],rtol=1e-13,atol=0)
        assert radial.heat==pytest.approx(vector.heat,rel=1e-13)
        # Independently derive the radial/tangent ratio from the temporal rank-one matrix.
        coeff=closure.coefficients(amplitude,.3)
        expected=1+4*(SCALES.delta0_J*amplitude*coeff["taupsi_ps"]*1e-12/HBAR)**2
        assert vector.mobility[1,1]/vector.mobility[0,0]==pytest.approx(expected,rel=1e-13)
    zero=closure.tensor_response([0,0],[0,0],0)
    np.testing.assert_array_equal(zero.velocity,[0,0]); assert zero.heat==0
    assert zero.Tmob_K==SCALES.Tb_K and zero.R==1


def test_rotation_covariance_without_phase_division():
    phi=.723; rotation=np.array([[np.cos(phi),-np.sin(phi)],[np.sin(phi),np.cos(phi)]])
    model=KWTMobility(SCALES); z=np.array([.4,.8]); force=np.array([-.7,.2])
    base=model.tensor_response(z,force,.15); turned=model.tensor_response(rotation@z,rotation@force,.15)
    np.testing.assert_allclose(turned.velocity,rotation@base.velocity,rtol=1e-12,atol=1e-13)
    np.testing.assert_allclose(turned.mobility,rotation@base.mobility@rotation.T,rtol=1e-12,atol=1e-13)
    assert turned.heat==pytest.approx(base.heat,rel=1e-12)


def test_reference_time_conversion_and_temperature_floor():
    first=KWTMobility(SCALES).amplitude_response(.8,1.2,.1)
    second=KWTMobility(replace(SCALES,t_ref_ps=3)).amplitude_response(.8,1.2,.1)
    assert second.velocity/3==pytest.approx(first.velocity,rel=1e-13)
    assert second.heat/3==pytest.approx(first.heat,rel=1e-13)
    assert second.taupsi_ps==first.taupsi_ps
    assert SCALES.temperature_to_K(SCALES.temperature_bar(2.5))==pytest.approx(2.5)
    with pytest.raises(TypeError):
        KWTMobility(SCALES).amplitude_response(.8,1,.1,tau_kin=1.)


@pytest.mark.parametrize("name",["delta0_J","N0_per_J_m3","Tc_K","Tb_K","t_ref_ps"])
@pytest.mark.parametrize("value",[0.,-1.,np.nan,np.inf])
def test_invalid_scales_are_rejected(name,value):
    with pytest.raises(ValueError): replace(SCALES,**{name:value})


@pytest.mark.parametrize("amplitude,force,temperature",[(-1,1,.1),(1,np.nan,.1),(1,1,-1),(1,1,np.inf)])
def test_invalid_KWT_input(amplitude,force,temperature):
    with pytest.raises(ValueError): KWTMobility(SCALES).amplitude_response(amplitude,force,temperature)


@pytest.mark.parametrize("low",[.02,.01,.005])
def test_Debye_Lobatto_modes_lambda_and_energy_moment_in_both_unit_systems(low):
    nat=10*SCALES.N0_per_J_m3*SCALES.delta0_J
    model=DebyePhonons(SCALES,4.,nat,.1,"unit-test synthetic Debye",low)
    grid=model.quadrature(5); ratio=low/4
    assert grid.energies[0]==low and grid.energies[-1]==4
    assert np.all(grid.capacities>0)
    np.testing.assert_allclose(grid.capacities,grid.weights*grid.dos,rtol=0,atol=0)
    assert np.sum(grid.capacities)==pytest.approx(3*nat/(SCALES.N0_per_J_m3*SCALES.delta0_J)*(1-ratio**3),rel=1e-11)
    mode_SI=np.sum(grid.weights*model.dos_SI(grid.energies*SCALES.delta0_J))*SCALES.delta0_J
    assert mode_SI==pytest.approx(3*nat*(1-ratio**3),rel=1e-11)
    assert 2*np.dot(grid.weights,model.alpha2F(grid.energies)/grid.energies)==pytest.approx(.1*(1-ratio**2),rel=1e-11)
    moment=9*nat/(SCALES.N0_per_J_m3*SCALES.delta0_J)*4/4*(1-ratio**4)
    assert np.dot(grid.capacities,grid.energies)==pytest.approx(moment,rel=1e-11)
    assert model.metadata()["infrared_omitted_fractions"]["lambda"]==ratio**2
    np.testing.assert_array_equal(model.alpha2F([0,low/2,4.1]),[0,0,0])


@pytest.mark.parametrize("temperature",[.03,.12,.6,2.,10.])
def test_Debye_thermal_energy_capacity_against_adaptive_quadrature(temperature):
    model=DebyePhonons(SCALES,4,10*SCALES.N0_per_J_m3*SCALES.delta0_J,.1,"thermal synthetic",.005)
    grid=model.quadrature(129)
    coefficient=90/4**3
    def u_integrand(e): return coefficient*e**3*np.exp(-e/temperature)/(-np.expm1(-e/temperature))
    def c_integrand(e): return coefficient*e**4/temperature**2*np.exp(-e/temperature)/(-np.expm1(-e/temperature))**2
    u=quad(u_integrand,.005,4,epsabs=1e-13,epsrel=1e-12)[0]
    c=quad(c_integrand,.005,4,epsabs=1e-13,epsrel=1e-12)[0]
    assert grid.thermal_energy(temperature)==pytest.approx(u,rel=1e-8)
    assert grid.thermal_capacity(temperature)==pytest.approx(c,rel=1e-8)


def test_Bose_zero_and_classical_limits_and_explicit_no_zero_grid():
    model=DebyePhonons(SCALES,4,10*SCALES.N0_per_J_m3*SCALES.delta0_J,.1,"limits synthetic",.01)
    grid=model.quadrature(9)
    assert grid.thermal_energy(0)==grid.thermal_capacity(0)==0
    t=1e10
    assert grid.thermal_energy(t)/t==pytest.approx(np.sum(grid.capacities),rel=1e-8)
    assert grid.thermal_capacity(t)==pytest.approx(np.sum(grid.capacities),rel=1e-8)
    with pytest.raises(ValueError): bose_occupation([0,1],.1)
    with pytest.raises(ValueError): replace(model,infrared_cutoff_bar=0).quadrature(9)


@pytest.mark.parametrize("updates",[{"lambda_eph":0},{"atom_density_m3":-1},{"synthetic_label":""},{"infrared_cutoff_bar":4}])
def test_invalid_Debye_declarations(updates):
    model=DebyePhonons(SCALES,4,1e26,.1,"synthetic",.01)
    with pytest.raises(ValueError): replace(model,**updates)


def make_conditional_fixture(tmp_path):
    data=np.array([[0,0,0,0],[1,.1,2,1],[2,.2,3,2],[3,0,0,3]])
    path=tmp_path/"derived.csv"; np.savetxt(path,data,delimiter=",",header="axis,alpha,dos,source_row",comments="")
    manifest={"schema":"pysnspd.derived-phonon-shape.v1","derived_csv_sha256":hashlib.sha256(path.read_bytes()).hexdigest(),"source_sha256":"f"*64,"source_url":"https://example.org/source","operations":[{"name":"synthetic_fixture"}]}
    mp=tmp_path/"manifest.json"; mp.write_text(json.dumps(manifest))
    return path,mp


def test_conditional_mapping_is_explicit_hash_bound_and_not_SI_admission(tmp_path):
    path,manifest=make_conditional_fixture(tmp_path)
    shape=ConditionalPhononShape(path,manifest,source_axis_to_energy_bar=2,dos_ordinate_to_dos_bar=5,synthetic_label="explicit shape comparison",infrared_cutoff_bar=.1)
    np.testing.assert_allclose(shape.alpha2F([0,2,4,6,7]),[0,.1,.2,0,0])
    np.testing.assert_allclose(shape.dos_bar([0,2,4,6,7]),[0,10,15,0,0])
    assert not shape.metadata()["SI_admitted"]
    assert shape.metadata()["dos_ordinate_to_dos_bar"]==5
    path.write_bytes(path.read_bytes()+b"\n")
    with pytest.raises(ValueError,match="hash"):
        ConditionalPhononShape(path,manifest,source_axis_to_energy_bar=2,dos_ordinate_to_dos_bar=5,synthetic_label="same")


def test_invalid_phonon_populations_are_never_clipped():
    grid=PhononQuadrature(np.array([.1,1]),np.ones(2),np.ones(2),{})
    for p in [[-1e-16,0],[0,np.inf],[0]]:
        with pytest.raises(ValueError): grid.energy(p)


def test_finite_quadrature_factors_cannot_create_infinite_capacities():
    with pytest.raises(ValueError,match="capacities"):
        PhononQuadrature(np.array([.1,1]),np.full(2,1e308),np.full(2,1e308),{})


@pytest.mark.parametrize("axis_factor,dos_factor",[(1e308,1.),(1.,1e308)])
def test_conditional_shape_mapping_rejects_derived_overflow(tmp_path,axis_factor,dos_factor):
    path,manifest=make_conditional_fixture(tmp_path)
    with pytest.raises(ValueError,match="nonfinite"):
        ConditionalPhononShape(path,manifest,source_axis_to_energy_bar=axis_factor,
            dos_ordinate_to_dos_bar=dos_factor,synthetic_label="invalid numerical mapping")
