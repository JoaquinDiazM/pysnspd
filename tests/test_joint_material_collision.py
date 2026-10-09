"""Physical/sign tests of joint material Born building blocks."""
import numpy as np
import pytest
from scipy.constants import hbar, Boltzmann
from pysnspd.experimental import adiabatic_matrix_usadel as ad
from pysnspd.experimental.joint_material_collision import (
    InelasticCalibrators, ElectronPhononBornQuadrature,collision_source,
    weak_born_first_order, require_general_phonon_feedback,
    require_general_bgk_self_energy,GeneralPhononFeedbackRequired)

EU=Boltzmann*8.65


def normal(energy,temperature=4.):
    energy=np.asarray(energy)
    R=ad.Jet.of(np.broadcast_to(ad.TAU3,energy.shape+(2,2)).copy())
    H=ad.distribution(np.tanh(energy/(2*Boltzmann*temperature)),np.zeros(energy.shape))
    return R,H


def build(energies,temperature=4.):
    omega=EU*np.array([.11,.31,.77])
    weights=EU*np.array([.10,.21,.34])
    op=ElectronPhononBornQuadrature(omega,weights,np.array([.001,.004,.011]),
        energy_unit_J=EU,shift_contract="Exact independent E +/- Omega evaluations")
    target=np.asarray(energies)
    rp,hp=normal(target[None,:]+omega[:,None],temperature)
    rm,hm=normal(target[None,:]-omega[:,None],temperature)
    n=1/np.expm1(omega/(Boltzmann*temperature))
    return op,rp,hp,rm,hm,n


def test_allmaras_nominal_free_reference_scales():
    nominal=InelasticCalibrators()
    assert nominal.tau_ee_Tc_s==5e-12
    assert nominal.tau_ep_Tc_s==24.7e-12
    assert nominal.tau0_s==pytest.approx(1.86840e-9,rel=2e-5)
    assert nominal.tau_sc_s(8.65,8.65)==pytest.approx(1/(1/5e-12+1/24.7e-12))
    varied=InelasticCalibrators(7e-12,30e-12)
    assert varied.tau_sc_s(17.3,8.65)==pytest.approx(1/(2/7e-12+8/30e-12))


def test_normal_equilibrium_detailed_balance_and_fdt():
    energy=EU*np.array([-.8,-.1,.2,1.2])
    op,rp,hp,rm,hm,n=build(energy)
    sigma=op.evaluate(rp,hp,rm,hm,n)
    r,H=normal(energy)
    result=collision_source(r,H,sigma,np.ones(len(energy)))
    assert np.max(abs(result.projected_source))<64*np.finfo(float).eps*np.max(abs(sigma.retarded_J.value))/EU
    h=H.value[:,0,0]
    assert np.allclose(sigma.keldysh_J.value,
        (sigma.retarded_J.value-sigma.advanced_J.value)*h[:,None,None],
        atol=2e-39,rtol=2e-13)
    assert max(sigma.symmetry_errors_J().values())<1e-38
    assert np.all(sigma.retarded_J.value[:,0,0].imag<0)
    assert np.all(sigma.retarded_generator.value[:,0,0].real>0)


def test_normal_nonequilibrium_golden_rule_sign_and_prefactor():
    energy=EU*np.array([-.6,0.,.4,1.1])
    op,rp,hp,rm,hm,n=build(energy)
    rng=np.random.default_rng(41)
    plus=rng.uniform(-.9,.9,(3,4));minus=rng.uniform(-.9,.9,(3,4))
    h=rng.uniform(-.9,.9,4)
    hp,hm=ad.distribution(plus,0.),ad.distribution(minus,0.)
    r,_=normal(energy);H=ad.distribution(h,0.)
    sigma=op.evaluate(rp,hp,rm,hm,n)
    area=np.array([1.,2.,3.,4.])
    result=collision_source(r,H,sigma,area)
    bracket=(2*n[:,None]+1)*(plus+minus-2*h)+h*(plus-minus)
    expected=np.pi/hbar*np.sum((op.quadrature_weights_J*op.alpha2F)[:,None]*bracket,axis=0)
    assert np.allclose(ad.projections(result.collision_rate_per_s.value)[:,0],
                       expected,atol=2e-5,rtol=1e-14)
    assert np.max(abs(ad.projections(result.collision_rate_per_s.value)[:,1]))==0
    assert np.allclose(result.projected_source[:,0],-hbar/(2*EU)*area*expected,rtol=1e-14)


def test_charge_imbalance_mode_not_discarded():
    energy=EU*np.array([-.6,0.,.4,1.1])
    op,rp,hp,rm,hm,n=build(energy)
    hp=ad.distribution(np.full((3,4),.15),np.full((3,4),.21))
    hm=ad.distribution(np.full((3,4),-.13),np.full((3,4),-.09))
    sigma=op.evaluate(rp,hp,rm,hm,n)
    r,_=normal(energy)
    result=collision_source(r,ad.distribution(np.full(4,.31),np.full(4,-.17)),sigma,np.ones(4))
    assert np.max(abs(result.projected_source[:,1]))>1e-6
    assert max(sigma.symmetry_errors_J().values())<1e-38


def test_superconducting_anomalous_entries_are_retained():
    energy=EU*np.array([.3,1.2])
    op,rp,hp,rm,hm,n=build(energy)
    angle=np.array([[.31+.1j,.27+.2j]]*3)
    R=np.zeros((3,2,2,2),complex)
    R[:,:,0,0]=np.cosh(angle);R[:,:,1,1]=-np.cosh(angle)
    R[:,:,0,1]=np.sinh(angle);R[:,:,1,0]=-np.sinh(angle)
    rp,rm=ad.Jet.of(R),ad.Jet.of(R*.93)
    hp=ad.distribution(np.full((3,2),.2),.1)
    hm=ad.distribution(np.full((3,2),-.3),-.1)
    sigma=op.evaluate(rp,hp,rm,hm,n)
    assert np.max(abs(sigma.retarded_J.value[:,0,1]))>0
    assert max(sigma.symmetry_errors_J().values())<2e-38


def test_phonon_time_derivative_and_first_are_retained():
    op,rp,hp,rm,hm,n=build(EU*np.array([.2,.8]))
    ndot=np.array([.4,-.3,.2])*1e10
    base=op.evaluate(rp,hp,rm,hm,n,phonon_time_per_s=ndot)
    dt=1e-16
    plus=op.evaluate(rp,hp,rm,hm,n+dt*ndot)
    minus=op.evaluate(rp,hp,rm,hm,n-dt*ndot)
    assert np.allclose(base.retarded_J.time,
        (plus.retarded_J.value-minus.retarded_J.value)/(2*dt),rtol=2e-9,atol=1e-21)
    shifted=op.evaluate(ad.Jet.of(rp.value,first=rp.value*.03),hp,rm,hm,n)
    assert np.max(abs(shifted.retarded_J.first))>0
    assert max(shifted.symmetry_errors_J().values())<1e-38


def test_global_gauge_covariance():
    op,rp,hp,rm,hm,n=build(EU*np.array([.2,.8]))
    rng=np.random.default_rng(44)
    Rp=ad.Jet.of(rng.normal(size=(3,2,2,2))+1j*rng.normal(size=(3,2,2,2)))
    Rm=ad.Jet.of(rng.normal(size=(3,2,2,2))+1j*rng.normal(size=(3,2,2,2)))
    phase=ad.links(np.full((3,2),.71))
    original=op.evaluate(Rp,hp,Rm,hm,n)
    moved=op.evaluate(ad.transport(phase,Rp),hp,ad.transport(phase,Rm),hm,n)
    p=phase.value[0]
    for name in ("retarded_J","advanced_J","keldysh_J"):
        assert np.allclose(getattr(moved,name).value,
            p@getattr(original,name).value@np.swapaxes(p.conj(),-1,-2),
            atol=2e-38,rtol=2e-14)


def test_weak_born_counts_collision_only_at_first_order():
    energy=EU*np.array([.2,.8])
    op,rp,hp,rm,hm,n=build(energy)
    sigma=op.evaluate(rp,hp,rm,hm,n)
    R,H=normal(energy)
    result=weak_born_first_order(R,H,sigma,np.ones(2))
    assert np.max(abs(result.spectral_rhs_first))==0  # diagonal normal commute
    assert np.max(abs(result.projected_kinetic_rhs_first))<64*np.finfo(float).eps*np.max(abs(sigma.retarded_J.value))/EU
    assert result.max_self_energy_over_Eu>0
    # Existing fast spectral time derivatives cannot leak into weak source.
    moved=ad.Jet.of(R.value,time=R.value*1e12,first=R.value*.1)
    assert np.array_equal(weak_born_first_order(moved,H,sigma,np.ones(2)).kinetic_rhs_first,
                          result.kinetic_rhs_first)


def test_general_feedback_and_bgk_are_explicitly_guarded():
    op,rp,hp,rm,hm,n=build(EU*np.array([.2,.8]))
    sigma=op.evaluate(rp,hp,rm,hm,n)
    with pytest.raises(GeneralPhononFeedbackRequired,match="phonon bubble"):
        require_general_phonon_feedback(sigma)
    with pytest.raises(GeneralPhononFeedbackRequired,match="BGK"):
        require_general_bgk_self_energy(charge_imbalance=[0.,.1],moving_spectrum=False)
    with pytest.raises(GeneralPhononFeedbackRequired,match="BGK"):
        require_general_bgk_self_energy(charge_imbalance=0.,moving_spectrum=True)


def test_invalid_support_not_repaired():
    op,rp,hp,rm,hm,n=build(EU*np.array([.2,.8]))
    with pytest.raises(ValueError,match="nonnegative"):op.evaluate(rp,hp,rm,hm,-n)
    with pytest.raises(ValueError,match="kappa_s"):op.evaluate(rp,hp,rm,hm,n,kappa_s=1.)


def test_phonon_bubble_normal_matches_explicit_emission_absorption():
    from pysnspd.experimental.joint_material_collision import phonon_polarization_rhs
    omega=EU*np.array([.2,.5])
    lower=EU*np.array([-.8,-.1,.3,.9])
    electron_weight=EU*np.array([.11,.23,.19,.17])
    Rh,Hh=normal(lower[None,:]+omega[:,None],7.)
    Rl,Hl=normal(np.broadcast_to(lower,(2,4)),7.)
    n=1/np.expm1(omega/(Boltzmann*2.))
    alpha=np.array([.001,.006]);dos=np.array([1e65,3e65]);N0=2e47
    result=phonon_polarization_rhs(Rh,Hh,Rl,Hl,
        electron_weights_J=electron_weight,phonon_energy_J=omega,
        phonon_weights_J=EU*np.array([.1,.3]),alpha2F=alpha,
        g_ph_per_J_m3=dos,N0_per_J_m3=N0,phonon_occupation=n,energy_unit_J=EU)
    fh=(1-Hh.value[:,:,0,0].real)/2
    fl=(1-Hl.value[:,:,0,0].real)/2
    rates=(n[:,None]+1)*fh*(1-fl)-n[:,None]*fl*(1-fh)
    expected=4*np.pi*N0*alpha/(hbar*dos)*np.sum(electron_weight*rates,axis=1)
    assert np.allclose(result.dn_dt_per_s,expected,rtol=3e-14)
    assert np.all(result.dn_dt_per_s>0)


def test_phonon_bubble_thermal_detailed_balance_with_anomalous_fields():
    from pysnspd.experimental.joint_material_collision import phonon_polarization_rhs
    omega=EU*np.array([.2,.5]);lower=EU*np.array([-.8,-.1,.3,.9])
    Rh,Hh=normal(lower[None,:]+omega[:,None],4.)
    Rl,Hl=normal(np.broadcast_to(lower,(2,4)),4.)
    # Any anomalous spectral fields with scalar thermal H share FDT;
    # the pointwise thermal cancellation does not require normal DOS.
    x=np.array([[.17+.03j,.27+.1j,.31+.2j,.24+.04j]]*2)
    def matrix(z):
        out=np.zeros(z.shape+(2,2),complex)
        out[:,:,0,0]=np.cos(z);out[:,:,1,1]=-np.cos(z)
        out[:,:,0,1]=np.sin(z);out[:,:,1,0]=np.sin(z)
        return ad.Jet.of(out)
    Rh,Rl=matrix(x),matrix(x*.83)
    n=1/np.expm1(omega/(Boltzmann*4.))
    result=phonon_polarization_rhs(Rh,Hh,Rl,Hl,
        electron_weights_J=np.full(4,EU*.2),phonon_energy_J=omega,
        phonon_weights_J=EU*np.array([.1,.3]),alpha2F=np.array([.001,.006]),
        g_ph_per_J_m3=np.array([1e65,3e65]),N0_per_J_m3=2e47,
        phonon_occupation=n,energy_unit_J=EU)
    assert np.max(abs(result.bubble_integral_J))/EU<1e-14


def test_phonon_bubble_energy_ward_for_independent_full_nambu_pair():
    from pysnspd.experimental.joint_material_collision import phonon_polarization_rhs
    # One paired event proves the algebra before any E-cut convergence claim.
    omega=EU*.42;wp=EU*.13;we=EU*.17;alpha=.004;N0=2e47;dos=1e65
    low=EU*.31;high=low+omega
    rng=np.random.default_rng(111)
    Rlo=ad.Jet.of(rng.normal(size=(1,2,2))+1j*rng.normal(size=(1,2,2)))
    Rhi=ad.Jet.of(rng.normal(size=(1,2,2))+1j*rng.normal(size=(1,2,2)))
    Hlo,Hhi=ad.distribution(np.array([.2]),np.array([.11])),ad.distribution(np.array([.4]),np.array([-.13]))
    zero=ad.Jet.of(np.zeros((1,1,2,2)))
    wrap=lambda j:ad.Jet.of(j.value[None,...])
    op=ElectronPhononBornQuadrature([omega],[wp],[alpha],energy_unit_J=EU,
        shift_contract="One exact high-low pair for algebraic Ward identity")
    slo=op.evaluate(wrap(Rhi),wrap(Hhi),zero,zero,[.3])
    shi=op.evaluate(zero,zero,wrap(Rlo),wrap(Hlo),[.3])
    Ilo=collision_source(Rlo,Hlo,slo,np.ones(1)).collision_rate_per_s.value
    Ihi=collision_source(Rhi,Hhi,shi,np.ones(1)).collision_rate_per_s.value
    electron_power=-N0*we/4*(low*np.trace(Ilo[0]).real+high*np.trace(Ihi[0]).real)
    bubble=phonon_polarization_rhs(wrap(Rhi),wrap(Hhi),wrap(Rlo),wrap(Hlo),
        electron_weights_J=[we],phonon_energy_J=[omega],phonon_weights_J=[wp],
        alpha2F=[alpha],g_ph_per_J_m3=[dos],N0_per_J_m3=N0,
        phonon_occupation=[.3],energy_unit_J=EU)
    phonon_power=float(bubble.phonon_power_density_W_per_m3)
    assert abs(electron_power+phonon_power)<2e-14*max(abs(electron_power),abs(phonon_power))


def test_superconducting_particle_hole_pair_charge_ward():
    """Charge Ward on normalized local spectra, including hT and PHS partners."""
    omega=EU*.42;wp=EU*.13
    low=np.array([.31,-.73]);high=low+.42
    def spectrum(E):
        z=.04-1j*np.asarray(E);den=np.sqrt(z*z+.4**2)
        value=(z/den)[...,None,None]*ad.TAU3+(.4/den)[...,None,None]*np.array([[0,1],[1,0]])
        return ad.Jet.of(value)
    distribution=lambda E:ad.distribution(np.tanh(E/.7),.05*np.exp(-E**2/.49))
    Rl,Rh=spectrum(low),spectrum(high)
    Hl,Hh=distribution(low),distribution(high)
    wrap=lambda j:ad.Jet.of(j.value[None,...])
    zero=ad.Jet.of(np.zeros((1,2,2,2)))
    op=ElectronPhononBornQuadrature([omega],[wp],[.004],energy_unit_J=EU,
        shift_contract="Exact paired positive/negative energies, identical local gap")
    sl=op.evaluate(wrap(Rh),wrap(Hh),zero,zero,[.3])
    sh=op.evaluate(zero,zero,wrap(Rl),wrap(Hl),[.3])
    Il=collision_source(Rl,Hl,sl,np.ones(2)).collision_rate_per_s.value
    Ih=collision_source(Rh,Hh,sh,np.ones(2)).collision_rate_per_s.value
    charge=lambda I:np.trace(ad.TAU3@I,axis1=-2,axis2=-1).real
    scale=max(np.max(abs(charge(Il))),np.max(abs(charge(Ih))))
    assert abs(np.sum(charge(Il)+charge(Ih)))<2e-12*scale


def test_superconducting_pair_coherence_matches_frozen_material_events():
    """The same R/A/K vertex yields S:NN-R2R2 and R:NN+R2R2."""
    def spectrum(energy):
        z=.02-1j*energy;root=np.sqrt(z*z+.4**2)
        return (z*ad.TAU3+.4*np.array([[0,1],[1,0]]))/root
    low,high=spectrum(.7),spectrum(1.3)
    negative=spectrum(-.7)
    advanced=lambda R:-ad.TAU3@R.conj().T@ad.TAU3
    Qlo,Qhi,Qnegative=(R-advanced(R) for R in (low,high,negative))
    Nlo,Nhi=low[0,0].real,high[0,0].real
    Flo,Fhi=low[0,1].imag,high[0,1].imag
    assert np.trace(Qhi@Qlo).real/8==pytest.approx(Nhi*Nlo-Fhi*Flo,rel=2e-14)
    assert np.trace(Qhi@Qnegative).real/8==pytest.approx(Nhi*Nlo+Fhi*Flo,rel=2e-14)


def material_snapshot_inputs():
    E=EU*np.array([-1.2,-.3,.3,1.2])
    we=EU*np.array([.21,.13,.13,.21])
    omega=EU*np.array([.2,.7])
    wp=EU*np.array([.15,.31])
    nodes=2
    # Node-dependent electronic states but phonons independent of electron E.
    thermal=lambda values:normal(values,4.)
    R,H=thermal(np.broadcast_to(E[:,None],(4,nodes)))
    plus,minus=E[None,:]+omega[:,None],E[None,:]-omega[:,None]
    Rp,Hp=thermal(np.broadcast_to(plus[:,:,None],(2,4,nodes)))
    Rm,Hm=thermal(np.broadcast_to(minus[:,:,None],(2,4,nodes)))
    n=1/np.expm1(omega/(Boltzmann*4.))
    return dict(electron_energy_J=E,electron_weights_J=we,
        phonon_energy_J=omega,phonon_weights_J=wp,alpha2F=[.001,.005],
        g_ph_per_J_m3=[1e65,3e65],N0_per_J_m3=2e47,energy_unit_J=EU,
        retarded=R,distribution=H,retarded_plus=Rp,distribution_plus=Hp,
        retarded_minus=Rm,distribution_minus=Hm,
        shifted_energy_plus_J=plus,shifted_energy_minus_J=minus,
        phonon_occupation=n,area_weights=[1.,1.7],
        shift_contract="Exact analytic normal queries on signed symmetric E grid")


def test_one_material_provider_returns_both_paired_sectors():
    from pysnspd.experimental.joint_material_collision import evaluate_born_material_snapshot
    values=material_snapshot_inputs()
    result=evaluate_born_material_snapshot(**values)
    assert result.retarded_generator_first.shape==(4,2,2,2)
    assert result.full_nambu_source.shape==(4,2,2,2)
    assert result.projected_source.shape==(4,2,2)
    assert result.phonon_rate_per_s.shape==(2,2)
    assert result.shift_error_max_J==0
    assert result.general_material_provider_complete
    assert result.self_energy.phonon_feedback_complete
    require_general_phonon_feedback(result.self_energy)
    assert not result.energy_admission and not result.temporal_admission
    scale=np.max(abs(result.retarded_generator_first))
    assert np.max(abs(result.projected_source))<64*np.finfo(float).eps*scale*1.7
    assert np.max(abs(result.phonon_rate_per_s))<1e-20


def test_material_provider_retains_independent_finite_cut_balance_defect():
    from pysnspd.experimental.joint_material_collision import (
        evaluate_born_material_snapshot,phonon_polarization_rhs)
    values=material_snapshot_inputs()
    values["phonon_occupation"]=np.array([[.7,.9],[.4,.6]])
    result=evaluate_born_material_snapshot(**values)
    # Coarse finite E quadrature is deliberately NOT invariant under +/-Omega.
    # A real defect must be returned, rather than heat adjusted to close it.
    assert np.max(abs(result.collision_energy_defect_W_per_m3))>1e3
    assert np.array_equal(result.collision_energy_defect_W_per_m3,
        result.electron_collision_power_density_W_per_m3+
        result.phonon_collision_power_density_W_per_m3)
    assert not result.energy_admission


@pytest.mark.parametrize("change",["wrong_shift","positive_only","unpaired_weights","pauli"])
def test_material_provider_rejects_unpaired_or_unphysical_queries(change):
    from pysnspd.experimental.joint_material_collision import evaluate_born_material_snapshot
    values=material_snapshot_inputs()
    if change=="wrong_shift":
        values["shifted_energy_plus_J"]=values["shifted_energy_plus_J"].copy()
        values["shifted_energy_plus_J"][0,0]+=EU*.01
        message="E \\+/- Omega"
    elif change=="positive_only":
        values["electron_energy_J"]=abs(values["electron_energy_J"])
        message="ordered signed"
    elif change=="unpaired_weights":
        values["electron_weights_J"]=values["electron_weights_J"].copy()
        values["electron_weights_J"][0]*=1.1
        message="symmetric E and weights"
    else:
        values["distribution"]=ad.distribution(np.full((4,2),1.02),np.zeros((4,2)))
        message="Pauli support"
    with pytest.raises(ValueError,match=message):
        evaluate_born_material_snapshot(**values)


def test_complete_material_provider_retains_superconducting_charge_and_anomalous_sources():
    from pysnspd.experimental.joint_material_collision import evaluate_born_material_snapshot
    values=material_snapshot_inputs()
    phase=np.array([.21,.64])
    gap=.4*np.exp(1j*phase)
    def field(energies):
        energies=np.asarray(energies)
        x=energies[...,None]/EU
        z=.03-1j*x
        root=np.sqrt(z*z+abs(gap)**2)
        R=np.zeros(x.shape[:-1]+(2,2,2),complex)
        R[...,0,0]=z/root;R[...,1,1]=-z/root
        R[...,0,1]=gap/root;R[...,1,0]=gap.conj()/root
        H=ad.distribution(np.broadcast_to(np.tanh(x/.8),R.shape[:-2]),
            np.broadcast_to(.03*np.exp(-x*x/.4),R.shape[:-2]))
        return ad.Jet.of(R),H
    R,H=field(values["electron_energy_J"])
    Rp,Hp=field(values["shifted_energy_plus_J"])
    Rm,Hm=field(values["shifted_energy_minus_J"])
    values.update(retarded=R,distribution=H,retarded_plus=Rp,distribution_plus=Hp,
        retarded_minus=Rm,distribution_minus=Hm,
        phonon_occupation=np.array([[.4,.6],[.2,.3]]))
    result=evaluate_born_material_snapshot(**values)
    assert np.max(abs(result.retarded_generator_first[...,0,1]))>1e-6
    assert np.max(abs(result.full_nambu_source[...,0,1]))>1e-6
    assert np.max(abs(result.projected_source[...,1]))>1e-6
    assert np.max(abs(result.phonon_rate_per_s))>0
    assert max(result.self_energy.symmetry_errors_J().values())<1e-37
