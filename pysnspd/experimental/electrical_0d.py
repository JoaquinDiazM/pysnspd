"""v1.2.5 electrical-only phase with the resolved DC inductance restored.

This is an explicit reduction, not thermal equilibration. Native distributions
are archived, never reset or represented by fictitious cold spatial fields.
"""
from dataclasses import asdict, dataclass
import numpy as np
from scipy.linalg import expm, eig


@dataclass(frozen=True)
class AdmissionConfig:
    earliest_ns: float = 12.
    overlap_ns: float = 1.
    voltage_error_uV: float = 1.
    current_error_fraction: float = .005
    gap_error_fraction: float = .03
    normal_current_fraction: float = .001
    maximum_snapshot_gap_ps: float = 10.

    def __post_init__(self):
        if any(not np.isfinite(v) or v <= 0 for v in asdict(self).values()):
            raise ValueError('All electrical admission bounds must be positive and finite')


class InductiveCircuit:
    """Exact linear three-state circuit; polynomial lift integrates work separately.

    State deviations preserve precision near DC. The additional 9 products
    delta_i*delta_j and one work state integrate source-minus-Joule work without
    subtracting the DC source and bias-resistor powers numerically.
    """
    def __init__(self, circuit, film_inductance_H):
        self.c = dict(circuit)
        required = ('Lk_ext_H', 'R_bias_ohm', 'L_bias_H', 'R_load_ohm',
                    'C_couple_F', 'V_bias_V')
        if set(self.c) != set(required):
            raise ValueError('Exact six-parameter thesis circuit required')
        if any(not np.isfinite(self.c[k]) or self.c[k] <= 0 for k in required):
            raise ValueError('Positive finite circuit parameters required')
        if not np.isfinite(film_inductance_H) or film_inductance_H <= 0:
            raise ValueError('Positive resolved-film DC inductance required')
        self.Lfilm = float(film_inductance_H)
        c = self.c
        self.Ltotal = c['Lk_ext_H'] + self.Lfilm
        Rb, RL, Lb, C = (c[k] for k in
                        ('R_bias_ohm', 'R_load_ohm', 'L_bias_H', 'C_couple_F'))
        I0 = c['V_bias_V']/Rb
        self.equilibrium = np.array([I0, I0, 0.])
        self.A = np.array([[-(Rb+RL)/Lb, RL/Lb, -1/Lb],
                           [RL/self.Ltotal, -RL/self.Ltotal, 1/self.Ltotal],
                           [1/C, -1/C, 0.]])*1e-12  # time is ps
        self.H = np.array([Lb, self.Ltotal, C])
        self.K = np.zeros((13, 13))
        self.K[:3, :3] = self.A
        self.K[3:12, 3:12] = np.kron(self.A, np.eye(3)) + np.kron(np.eye(3), self.A)
        self.K[12, 0] = -c['V_bias_V']*1e-12
        Q = np.zeros((3, 3)); Q[0, 0] = -(Rb+RL); Q[1, 1] = -RL
        Q[0, 1] = Q[1, 0] = RL
        self.K[12, 3:12] = Q.ravel()*1e-12
        self._propagators = {}

    def initialize(self, state):
        y = np.asarray(state, float)
        if y.shape != (3,) or not np.isfinite(y).all():
            raise ValueError('Finite (Ib, Is, vc) state required')
        d = y-self.equilibrium
        return np.r_[d, np.outer(d, d).ravel(), 0.]

    def advance(self, lifted, h_ps):
        if not np.isfinite(h_ps) or h_ps <= 0:
            raise ValueError('Positive finite time step required')
        if h_ps not in self._propagators:
            if len(self._propagators) > 8:
                self._propagators.clear()
            self._propagators[h_ps] = expm(self.K*h_ps)
        z = self._propagators[h_ps] @ lifted
        if not np.isfinite(z).all():
            raise ArithmeticError('Nonfinite 0D state')
        return z

    def state(self, z):
        return self.equilibrium+z[:3]

    def energy_deviation(self, z):
        return float(np.dot(self.H*self.equilibrium, z[:3]) +
                     .5*np.dot(self.H, np.diag(z[3:12].reshape(3, 3))))

    def outputs(self, z):
        y = self.state(z)
        return dict(Ib_A=float(y[0]), Is_A=float(y[1]), vc_V=float(y[2]),
                    Vout_V=float(self.c['R_load_ohm']*(z[0]-z[1])),
                    Vdev_V=float(self.Lfilm*(self.A@z[:3])[1]*1e12))

    def electrical_distance(self, z):
        o = self.outputs(z); I0 = self.equilibrium[0]
        return float(max(abs(z[0])/(.001*I0), abs(z[1])/(.001*I0),
                         abs(z[2])/1e-6, abs(o['Vout_V'])/1e-6))


def qualify(circuit, trace, snapshot_times_ps, gap, dc_gap, manifest,
            config=AdmissionConfig()):
    """Require a persistent SC interval and independently evolved circuit overlap.

    A passing electrical overlap deliberately does not require bath temperature.
    All physical failure/energy/domain stops remain disqualifying.
    """
    t = np.asarray(trace['time_ps']); end = float(t[-1])
    reasons = []
    if end < config.earliest_ns*1000:
        reasons.append('Earlier than the declared electrical-reduction scope')
    start = end-config.overlap_ns*1000
    if len(t) < 3 or not np.isfinite(trace.view('f8')).all() or np.any(np.diff(t) <= 0):
        raise ValueError('Finite monotone accepted native trace required')
    if t[0] > start:
        reasons.append('Insufficient native overlap history')
    last = manifest.get('latest_monitor') or {}
    if last.get('boundary_alarm', True) or last.get('hot_length_nm', 1.) != 0:
        reasons.append('Native domain/hot-segment monitor is not clear')
    ident = manifest['parameter_identity']
    if manifest['domain_fraction'] > ident['guard_cut_energy_relative']:
        reasons.append('Native cut-energy guard failed')
    if manifest['max_global_energy_fraction'] > ident['global_energy_tolerance']:
        reasons.append('Native global energy guard failed')
    st = np.asarray(snapshot_times_ps)
    g = np.abs(np.asarray(gap)); dc = np.abs(np.asarray(dc_gap))
    if st.ndim != 1 or len(st) < 3 or g.shape != (len(st), len(dc)):
        raise ValueError('Consistent native snapshots required')
    if not np.isfinite(g).all() or not np.isfinite(dc).all() or np.any(dc <= 0):
        raise ValueError('Finite positive native gap reference required')
    if (st[0] > start or abs(st[-1]-end) > 1e-7 or
            np.any(np.diff(st) <= 0) or
            np.max(np.diff(st)) > config.maximum_snapshot_gap_ps):
        reasons.append('Native gap snapshots do not cover the complete overlap')
    scale = .5*(g[:, 0]/dc[0] + g[:, -1]/dc[-1])
    if np.any(scale <= 0):
        raise ValueError('Cold-current boundary amplitude vanished')
    gap_error = float(np.max(abs(g/(dc[None, :]*scale[:, None])-1)))
    if gap_error > config.gap_error_fraction:
        reasons.append('Full-domain gap differs from the current-matched SC branch')
    keep = t >= start
    times = np.r_[start, t[keep]]
    times = np.unique(times)
    y0 = np.array([np.interp(start, t, trace[k]) for k in ('Ib_A', 'Is_A', 'vc_V')])
    eigenvalues, vectors = eig(circuit.A)
    coefficients = np.linalg.solve(vectors, y0-circuit.equilibrium)
    deviations = (vectors @ (coefficients[:, None]*np.exp(
        eigenvalues[:, None]*(times-start)[None, :]))).T.real
    errors = dict(voltage_error_V=0., capacitor_error_V=0., current_relative=0.,
                  film_voltage_residual_V=0., normal_current_fraction=0.)
    for now, deviation in zip(times, deviations):
        actual = np.array([np.interp(now, t, trace[k]) for k in ('Ib_A', 'Is_A', 'vc_V')])
        vd = float(np.interp(now, t, trace['Vdev_V']))
        vout = circuit.c['R_load_ohm']*(actual[0]-actual[1])
        didt = (actual[2]+vout-vd)/circuit.c['Lk_ext_H']
        predicted = circuit.equilibrium+deviation
        o = dict(Vout_V=circuit.c['R_load_ohm']*(deviation[0]-deviation[1]), vc_V=predicted[2])
        errors['voltage_error_V'] = max(errors['voltage_error_V'], abs(o['Vout_V']-vout))
        errors['capacitor_error_V'] = max(errors['capacitor_error_V'], abs(o['vc_V']-actual[2]))
        errors['current_relative'] = max(errors['current_relative'],
            float(np.max(abs(predicted[:2]-actual[:2]))/circuit.equilibrium[0]))
        errors['film_voltage_residual_V'] = max(errors['film_voltage_residual_V'],
                                               abs(vd-circuit.Lfilm*didt))
        Rnormal = ident['sheet_resistance_ohm']*ident['length_1d_nm']/80.
        if min(actual[:2]) <= 0:
            reasons.append('Native global branch polarity failed'); break
        errors['normal_current_fraction'] = max(errors['normal_current_fraction'], abs(vd)/(Rnormal*actual[1]))
    for key in ('voltage_error_V', 'capacitor_error_V', 'film_voltage_residual_V'):
        if errors[key] > config.voltage_error_uV*1e-6:
            reasons.append(key+' exceeds the declared electrical margin')
    if errors['current_relative'] > config.current_error_fraction:
        reasons.append('Independent circuit current overlap failed')
    if errors['normal_current_fraction'] > config.normal_current_fraction:
        reasons.append('Normal-current scale separation failed')
    return dict(admitted=not reasons, reasons=reasons, config=asdict(config),
                overlap_start_ps=start, overlap_end_ps=end, checks=len(times),
                maximum_gap_error_fraction=gap_error, errors=errors,
                thermal_admission=False, thermal_audit_pending=True,
                last_native_temperature_K=last.get('max_temperature_K'),
                interpretation='Electrical surrogate only; no thermal reset or retrigger certification')
