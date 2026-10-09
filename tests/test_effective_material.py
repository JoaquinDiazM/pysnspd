import unittest
import tempfile
from pathlib import Path
import numpy as np
from scipy.constants import Boltzmann as kB,hbar,elementary_charge as e
from pysnspd.experimental.effective_material import allmaras_acoustic_material,acoustic_anchored_material,EffectivePhononMaterial
from pysnspd.experimental.kinetic_events import ElectronPhononEvents

class EffectiveMaterialTest(unittest.TestCase):
    def setUp(self):self.m=allmaras_acoustic_material(nodes=513)
    def test_acoustic_heat_and_derivative_in_SI(self):
        B=self.m.manifest['acoustic_B_J_per_m3_K4']
        for T in (.9,4.325,8.65,17.3):
            self.assertAlmostEqual(self.m.heat_capacity_J_per_m3_K(T)/(B*T**3),1.,delta=2e-10)
            h=T*1e-4
            fd=(self.m.energy_density_J_per_m3(T+h)-self.m.energy_density_J_per_m3(T-h))/(2*h)
            self.assertAlmostEqual(fd/self.m.heat_capacity_J_per_m3_K(T),1.,delta=2e-8)
        self.assertEqual(self.m.energy_density_J_per_m3(0),0)
    def test_Einstein_single_spin_and_A20_gamma(self):
        self.assertEqual(self.m.N0_per_J_m3,1/(2*e**2*5e-5*608*7e-9))
        self.assertAlmostEqual(self.m.manifest['gamma_A20']/self.m.manifest['gamma_independent'],1.,delta=1e-14)
    def test_normal_collision_rate_seconds_matches_tau0_without_fit(self):
        Eu=kB*8.65;omega=np.array([.07,.5,1.,5.,30.]);A,_,_=self.m.evaluate(omega*Eu)
        expected=2*omega**2/self.m.manifest['tau0_s']
        np.testing.assert_allclose(2*np.pi/hbar*A*Eu,expected,rtol=2e-15)
        delta=1.764*Eu
        for tref in (1e-14,1e-12,1e-9):
            native=self.m.event_rate_prefactor(delta,tref)/4*A*Eu/delta/tref
            np.testing.assert_allclose(native,expected,rtol=2e-15)
    def events(self,tref):
        delta=1.764*kB*8.65
        pp=np.linspace(.1,8.,80);weights=np.full_like(pp,.1)*delta
        grid=self.m.phonon_grid(pp*delta,weights,delta)
        ep=np.array([.8,1.3,2.,3.5]);ew=np.array([.1,.2,.15,.25])
        return ElectronPhononEvents(ep,ew,np.array([.1,.2,.3,.1]),grid,
            rate_prefactor=self.m.event_rate_prefactor(delta,tref),label='test supplied spectrum; not detector admission'),delta
    def test_shared_native_events_detailed_balance_and_energy(self):
        events,delta=self.events(1e-12);temp=kB*7./delta
        p=1/(1+np.exp(events.electron_energies/temp));n=1/np.expm1(events.phonons.energies/temp)
        rates=events.rates(p,n)
        self.assertLess(np.max(abs(rates)),1e-16)
        p=p+np.array([.04,.03,.02,.01]);dp,dn=events.rhs(p,n)
        electronic=float(np.dot(events.electron_capacities*events.electron_energies,dp))
        phononic=float(np.dot(events.phonons.capacities*events.phonons.energies,dn))
        self.assertGreater(abs(electronic),1e-7)
        self.assertLess(abs(electronic+phononic),2e-15*abs(electronic))
    def test_native_rates_convert_to_identical_seconds(self):
        a,_=self.events(1e-12);b,_=self.events(3e-12)
        p=np.array([.2,.17,.11,.07]);n=np.full(len(a.phonons.energies),.03)
        for x,y in zip(a.rhs(p,n),b.rhs(p,n)):np.testing.assert_allclose(x/1e-12,y/3e-12,rtol=4e-15)
    def test_arbitrary_shape_ordinate_cancels_from_acoustic_anchor(self):
        energy=np.linspace(0,kB*100,301);shape=(energy/(kB*10))**2
        def build(f):return acoustic_anchored_material(energy,shape,f*shape,reference_temperature_K=8.65,
            sound_speed_m_per_s=4912,N0_per_J_m3=self.m.N0_per_J_m3,shape_contract='arbitrary unit shape',alpha2F_contract='explicit dimensionless A')
        a,b=build(1.),build(17.)
        np.testing.assert_allclose(a.g_ph_per_J_m3,b.g_ph_per_J_m3,rtol=4e-15)
        self.assertFalse(a.manifest['mode_count_normalized'])
    def test_normal_phonon_scattering_and_unordered_pair_factors_from_primary(self):
        # Delta=Eu so all variables have the source's kB Tc normalization.
        # Exact event brackets, unit energy weights and normal coherence.
        Eu=kB*8.65;tref=1e-12;pp=np.array([1.,2.,3.,4.])
        phonons=self.m.phonon_grid(pp*Eu,np.ones(4)*Eu,Eu)
        events=ElectronPhononEvents(np.array([1.,2.]),np.ones(2),np.zeros(2),phonons,
            rate_prefactor=self.m.event_rate_prefactor(Eu,tref),label='normal source factor test')
        gamma=self.m.manifest['gamma_A20'];tau=self.m.manifest['tau0_s']
        # One scattering event: source Eq4 has 2 gamma/tau0.
        scattering=np.flatnonzero(~events.recombination)[0]
        omega=events.omega[scattering];q=np.flatnonzero(pp==omega)[0]
        actual=events.coefficients[scattering]/phonons.capacities[q]/tref
        self.assertAlmostEqual(actual/(2*gamma/tau),1.,delta=1e-14)
        # Off-diagonal unordered pair accounts for two entries in the source's
        # ordered pair integral; the diagonal has one entry and half event rate.
        for event in np.flatnonzero(events.recombination):
            q=np.flatnonzero(pp==events.omega[event])[0]
            multiplicity=1 if events.index_i[event]==events.index_j[event] else 2
            actual=events.coefficients[event]/phonons.capacities[q]/tref
            self.assertAlmostEqual(actual/(multiplicity*gamma/tau),1.,delta=1e-14)

    def test_roundtrip_and_support_no_silent_extrapolation(self):
        with tempfile.TemporaryDirectory() as td:
            self.m.save(Path(td)/'m.npz',Path(td)/'m.json');loaded=EffectivePhononMaterial.load(Path(td)/'m.npz')
            np.testing.assert_allclose(loaded.evaluate(kB*np.array([1.,7.]))[0],self.m.evaluate(kB*np.array([1.,7.]))[0],rtol=0,atol=0)
        with self.assertRaises(ValueError):self.m.evaluate(kB*665)
        with self.assertRaises(ValueError):allmaras_acoustic_material(tau0_s=0)

if __name__=='__main__':unittest.main()
