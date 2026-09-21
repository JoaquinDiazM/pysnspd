"""Homogeneous one/two-cell kinetic tests of D.11--D.20.

The two superconducting amplitudes are local variables; the connection carries
only quasiparticle diffusion. This stage deliberately does not discretize the
spatial condensate functional, phase, electrostatic potential, or a circuit.
All energies and forces use the SAME native catalogue quadrature.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .cell_validation import ElectronicCell
from .cell_transport import NativeTransportEvents
from .kinetic_events import ElectronPhononEvents, ProjectedElectronPhononEvents, PhononGrid


@dataclass
class CoupledCellSystem:
    catalog: object
    phonons: PhononGrid
    mobility: object
    gammas: tuple[float, ...] = (0.0,)
    tau_kin: float = .7
    tau_escape: float = 15.0
    diffusion_over_length_squared: float = .2
    face_order: int = 2
    external_powers: tuple[float, ...] = (.01,)
    rate_prefactor: float = 1.0
    reaction_quadrature_order: int = 2
    reaction_method: str = 'projected'
    reaction_layout: str = 'resolved_panels'
    reaction_outer_order: int = 2
    reaction_max_panel: float = .5
    reaction_max_energy_panel: float = .125

    def __post_init__(self):
        if len(self.gammas) not in (1, 2) or len(self.external_powers) != len(self.gammas):
            raise ValueError("one or two cells and one external power per cell required")
        if any(not np.isfinite(v) or v < 0 for v in self.external_powers):
            raise ValueError("external heating powers must be finite and nonnegative")
        for name in ("tau_kin", "tau_escape"):
            value = getattr(self, name)
            if value <= 0 or np.isnan(value):
                raise ValueError(name+" must be positive (infinity disables the term)")
        self.electron_size = len(self.catalog.count_nodes)
        self.phonon_size = len(self.phonons.energies)
        self.cell_count = len(self.gammas)
        self.block_size = 1+self.electron_size+self.phonon_size
        self.population_size = self.cell_count*self.block_size
        # Four independently integrated diagnostic ledgers per cell:
        # escape, external input, net e->ph transfer, condensate dissipation.
        # Last entry is the net left->right electronic transport.
        self.state_size = self.population_size+4*self.cell_count+1
        self.bath_temperature = self.mobility.scales.bath_temperature_bar
        self.bath_phonons = 1/np.expm1(self.phonons.energies/self.bath_temperature)
        self.rhs_calls = 0

    def reaction_events(self, cell):
        if self.reaction_method == 'projected':
            return ProjectedElectronPhononEvents.from_cell(
                cell, self.phonons, rate_prefactor=self.rate_prefactor,
                quadrature_order=self.reaction_quadrature_order,
                integration_layout=self.reaction_layout,
                outer_order=self.reaction_outer_order,
                max_count_panel=self.reaction_max_panel,
                max_energy_panel=self.reaction_max_energy_panel)
        if self.reaction_method == 'native_pairs':
            return ElectronPhononEvents.from_cell(cell, self.phonons,
                                                  rate_prefactor=self.rate_prefactor)
        raise ValueError('unknown reaction quadrature method')

    def pack(self, amplitudes, electrons, phonons):
        if len(amplitudes) != self.cell_count:
            raise ValueError("incorrect cell count")
        state = np.zeros(self.state_size)
        for i in range(self.cell_count):
            start = i*self.block_size
            state[start] = amplitudes[i]
            state[start+1:start+1+self.electron_size] = electrons[i]
            state[start+1+self.electron_size:start+self.block_size] = phonons[i]
        self.unpack(state)
        return state

    def unpack(self, state):
        state = np.asarray(state, float)
        if state.shape != (self.state_size,) or np.any(~np.isfinite(state)):
            raise ValueError("one finite complete cell state required")
        amplitudes, electrons, phonons = [], [], []
        for i in range(self.cell_count):
            start = i*self.block_size
            amplitudes.append(state[start])
            p = state[start+1:start+1+self.electron_size]
            n = state[start+1+self.electron_size:start+self.block_size]
            if np.any(p < 0) or np.any(p > 1) or np.any(n < 0):
                raise ValueError("time stage outside Pauli/phonon support; no clipping applied")
            electrons.append(p); phonons.append(n)
        return np.asarray(amplitudes), np.asarray(electrons), np.asarray(phonons)

    def ledgers(self, state):
        values = np.asarray(state)[self.population_size:]
        return values[:-1].reshape(self.cell_count, 4), float(values[-1])

    def snapshot(self, state):
        amplitudes, p, n = self.unpack(state)
        electronic, phononic, temperatures, forces, excitations = [], [], [], [], []
        for i in range(self.cell_count):
            cell = ElectronicCell(self.catalog, amplitudes[i], self.gammas[i])
            energy, force, _ = cell.moments(p[i])
            electronic.append(energy); forces.append(force)
            excitations.append(cell.excitation_energy(p[i]))
            phononic.append(np.dot(self.phonons.capacities*self.phonons.energies, n[i]))
            temperatures.append(cell.equivalent_temperature(p[i]))
        ledger, transport = self.ledgers(state)
        return dict(amplitudes=amplitudes, electron_energy=np.asarray(electronic),
                    excitation_energy=np.asarray(excitations),
                    phonon_energy=np.asarray(phononic), temperatures=np.asarray(temperatures),
                    forces=np.asarray(forces), escape=ledger[:, 0], input=ledger[:, 1],
                    electron_to_phonon=ledger[:, 2], condensate_heat=ledger[:, 3],
                    transport=transport, total=float(np.sum(electronic)+np.sum(phononic)),
                    conserved=float(np.sum(electronic)+np.sum(phononic)+np.sum(ledger[:, 0]-ledger[:, 1])),
                    minimum_electron=float(p.min()), maximum_electron=float(p.max()),
                    minimum_phonon=float(n.min()))

    def rhs(self, time, state):
        self.rhs_calls += 1
        amplitudes, p, n = self.unpack(state)
        result = np.zeros_like(state)
        ledger = result[self.population_size:-1].reshape(self.cell_count, 4)
        cells, derivatives = [], []
        for i in range(self.cell_count):
            cell = ElectronicCell(self.catalog, amplitudes[i], self.gammas[i])
            cells.append(cell)
            event = self.reaction_events(cell)
            dp, dn = event.rhs(p[i], n[i])
            eph_power = float(np.dot(self.phonons.capacities*self.phonons.energies, dn))
            temperature = cell.equivalent_temperature(p[i])
            force = cell.moments(p[i])[1]
            motion = self.mobility.amplitude_response(amplitudes[i], force, temperature)
            power = self.external_powers[i]+motion.heat
            dp += cell.heating(p[i], power, self.bath_temperature)
            if np.isfinite(self.tau_kin):
                dp += (cell.fermi_dirac(temperature)-p[i])/self.tau_kin
            escaping = (n[i]-self.bath_phonons)/self.tau_escape
            dn -= escaping
            start = i*self.block_size
            result[start] = motion.velocity
            result[start+1+self.electron_size:start+self.block_size] = dn
            derivatives.append(dp)
            ledger[i] = (np.dot(self.phonons.capacities*self.phonons.energies, escaping),
                         self.external_powers[i], eph_power, motion.heat)
        if self.cell_count == 2:
            face = NativeTransportEvents(*cells, self.diffusion_over_length_squared,
                                         self.face_order)
            exchange = face.rhs(p)
            for i in range(2):
                derivatives[i] += exchange[i]
            result[-1] = face.transferred_power(p)
        for i, dp in enumerate(derivatives):
            start = i*self.block_size
            result[start+1:start+1+self.electron_size] = dp
        if np.any(~np.isfinite(result)):
            raise FloatingPointError("nonfinite coupled derivative")
        return result

    def instantaneous_energy_residual(self, state):
        """Differentiate the same native energy, independently of saved ledgers."""
        amplitudes, p, _ = self.unpack(state)
        derivative = self.rhs(0.0, state)
        value = 0.0
        for i in range(self.cell_count):
            start = i*self.block_size
            cell = ElectronicCell(self.catalog, amplitudes[i], self.gammas[i])
            dp = derivative[start+1:start+1+self.electron_size]
            dn = derivative[start+1+self.electron_size:start+self.block_size]
            value += cell.moments(p[i])[1]*derivative[start]
            value += 4*np.dot(cell.weights*cell.energies, dp)
            value += np.dot(self.phonons.capacities*self.phonons.energies, dn)
        ledgers = derivative[self.population_size:-1].reshape(self.cell_count, 4)
        return float(value+np.sum(ledgers[:, 0]-ledgers[:, 1]))
