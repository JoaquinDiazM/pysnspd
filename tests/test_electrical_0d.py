"""Independent circuit/energy checks and authenticated phase lifecycle."""
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import h5py
import numpy as np
import pytest
from scipy.integrate import solve_ivp
from pysnspd.experimental.electrical_0d import AdmissionConfig, InductiveCircuit, qualify
from pysnspd.experimental import electrical_continuation as E

C = dict(Lk_ext_H=101.71486210109221e-9, R_bias_ohm=1e4,
         L_bias_H=1e-6, R_load_ohm=50., C_couple_F=1e-10, V_bias_V=.155)
LFILM = 1.428715525273054e-9


def test_circuit_matches_independent_ode_and_preserves_inductive_voltage():
    model = InductiveCircuit(C, LFILM)
    y = np.array([15.504e-6, 17e-6, 35e-6]); z = model.initialize(y)
    first = model.outputs(z)
    assert first['Vdev_V'] != 0
    def rhs(t, state):
        Ib, Is, vc, work = state
        dIb = (.155-1e4*Ib-vc-50*(Ib-Is))/1e-6
        dIs = (vc+50*(Ib-Is))/(C['Lk_ext_H']+LFILM)
        dvc = (Ib-Is)/1e-10
        d = Ib-15.5e-6
        power = -.155*d-1e4*d*d-50*(Ib-Is)**2
        return np.array([dIb, dIs, dvc, power])*1e-12
    sol = solve_ivp(rhs, (0., 2000.), np.r_[y, 0.], method='DOP853',
                    rtol=2e-11, atol=[1e-18, 1e-18, 1e-18, 1e-30], max_step=5.)
    final = model.advance(z, 2000.)
    np.testing.assert_allclose(model.state(final), sol.y[:3, -1], rtol=2e-10, atol=1e-17)
    assert abs(final[12]-sol.y[3, -1]) < 1e-26
    assert abs(model.energy_deviation(final)-model.energy_deviation(z)-final[12]) < 1e-28


def test_dc_equilibrium_and_precision_at_100ns():
    model = InductiveCircuit(C, LFILM)
    z = model.initialize(model.equilibrium)
    assert np.array_equal(model.advance(z, 88000.), z)
    assert model.outputs(z)['Vdev_V'] == 0
    z = model.initialize([15.504e-6, 17e-6, 35e-6])
    final = model.advance(z, 88000.)
    assert model.electrical_distance(final) < 1e-7
    assert np.any(final[:3] != 0)  # no forced reset at finite time


@pytest.mark.parametrize('bad', [0., -1., float('nan')])
def test_bad_inductance_and_admission_margins_rejected(bad):
    with pytest.raises(ValueError):
        InductiveCircuit(C, bad)
    with pytest.raises(ValueError):
        AdmissionConfig(overlap_ns=bad)


@pytest.fixture
def native(tmp_path, monkeypatch):
    monkeypatch.setattr(E, 'SCRATCH', tmp_path)
    source = tmp_path/'native/runs/case'; source.mkdir(parents=True)
    model = InductiveCircuit(C, LFILM)
    times = np.arange(12000., 13202., 2.)
    z0 = model.initialize([15.504e-6, 17e-6, 35e-6])
    dtype = [(k, 'f8') for k in E.FIELDS+('regime_id',)]
    tr = np.zeros(len(times), dtype=dtype)
    for j, t in enumerate(times):
        z = z0 if j == 0 else model.advance(z0, t-times[0])
        o = model.outputs(z)
        for key, val in o.items():
            tr[key][j] = val
        tr['time_ps'][j] = t
        tr['energy_J'][j] = 1e-16 + model.energy_deviation(z)-model.energy_deviation(z0)
        tr['external_J'][j] = z[12]
    gap = np.ones((len(times), 3), complex); p = np.zeros((3, 2)); n = np.zeros((3, 2))
    pins = {k: E.sha(E.ROOT/k) for k in E.NATIVE_SOURCES | E.V124_SOURCES}
    ident = dict(sheet_resistance_ohm=608., length_1d_nm=1000.,
        guard_cut_energy_relative=.02, global_energy_tolerance=.01, local_energy_tolerance=2.5e-6)
    counts = {'1d': dict(trace=len(times), snapshot_time_ps=len(times),
                        delta_bar=len(times), p=len(times), n=len(times))}
    cm = np.array([tr[k][-1] for k in ('Ib_A', 'Is_A', 'vc_V')])
    meta = dict(phase='1d', time_ps=times[-1], source_pins_sha256=pins,
                parameter_identity=ident, archive_counts=counts, terminal_diagnostic=False,
                initial_energy_J=1e-16, external_J=tr['external_J'][-1], absolute_external_J=1e-17)
    E.atomic_npz(source/'checkpoint_latest.npz', delta_bar=gap[-1], p=p, n=n,
                 circuit_SI=cm, metadata_json=np.array(json.dumps(meta)))
    E.atomic_npz(source/'dc1d.npz', delta_bar=gap[-1])
    with h5py.File(source/'data.h5', 'w', libver='latest') as f:
        q = f.create_group('1d'); q.create_dataset('trace', data=tr)
        q.create_dataset('snapshot_time_ps', data=times); q.create_dataset('delta_bar', data=gap)
        q.create_dataset('dc_delta_bar', data=gap[-1]); q.create_dataset('p', data=np.broadcast_to(p, (len(times), 3, 2)))
        q.create_dataset('n', data=np.broadcast_to(n, (len(times), 3, 2)))
    m = dict(meta, model_version='1.2.4', physical_kernel_version='1.1.3', status='WALL_BUDGET_CHECKPOINTED',
        error=None, run_complete=False, cold_recovery={'active': False},
        checkpoint_sha256=E.sha(source/'checkpoint_latest.npz'), dc1d_sha256=E.sha(source/'dc1d.npz'),
        latest_monitor=dict(boundary_alarm=False, hot_length_nm=0., max_temperature_K=2.4),
        inductance_partition=dict(added_96nH_preserved=True, partition_frozen_from_DC=True,
            bulk_differential_H_per_m=LFILM/1e-6, length_1d_m=1e-6),
        circuit_by_phase={'1d': C}, domain_fraction=.001, max_global_energy_fraction=.001,
        energy_defect_J=0., source={'retained_energy_J': 1e-19}, observation={'crossing_time_s': 33e-12},
        monitor_config=dict(current_limit_uA=28.45, device_voltage_limit_mV=300.,
            readout_voltage_limit_mV=1.4, capacitor_voltage_limit_mV=310.))
    E.write_json(source/'manifest.json', m)
    return source, tmp_path/'out'


def test_phase_preserves_native_states_ledger_archive_and_resume(native):
    source, out = native
    before = {p.name: p.read_bytes() for p in source.iterdir()}
    plan = E.prepare(source, out, final_time_ns=14.)
    assert plan['admission']['admitted']
    assert plan['admission']['thermal_admission'] is False
    pause = E.execute(out, stop_requested=lambda: True)
    assert pause['status'] == 'INTERRUPTED_CHECKPOINTED'
    done = E.execute(out, resume=True)
    assert done['status'] == 'ELECTRICAL_HORIZON_COMPLETE' and done['phase'] == '0d'
    assert done['max_phase_energy_fraction'] < 1e-10
    with h5py.File(out/'data.h5', 'r') as f:
        assert isinstance(f.get('1d', getlink=True), h5py.ExternalLink)
        tr = f['0d/trace'][:]
        assert tr['time_ps'][0] == plan['initial_time_ps'] and tr['time_ps'][-1] == 14000.
        assert np.all(np.diff(tr['time_ps']) > 0) and np.all(tr['thermal_evolved'] == 0)
        assert tr['energy_J'][0] == plan['native_end']['energy_J']
        assert tr['external_J'][0] == plan['native_external_J']
        assert np.any(tr['Vdev_V'] != 0)
    completed = {p.name: p.read_bytes() for p in out.iterdir() if p.is_file() and p.name != '.execution.lock'}
    E.execute(out, resume=True)
    assert completed == {p.name: p.read_bytes() for p in out.iterdir() if p.is_file() and p.name != '.execution.lock'}
    assert before == {p.name: p.read_bytes() for p in source.iterdir()}


@pytest.mark.parametrize('damage', ['hash', 'version', 'status', 'cold', 'terminal', 'archive'])
def test_incompatible_or_damaged_native_checkpoint_refused(native, damage):
    source, out = native; p = source/'manifest.json'; m = json.loads(p.read_text())
    if damage == 'hash':
        m['checkpoint_sha256'] = 'bad'
    if damage == 'version':
        m['model_version'] = '1.2.2'
    if damage == 'status':
        m['status'] = 'POSSIBLE_DOMAIN_ERROR'
    if damage == 'cold':
        m['cold_recovery']['active'] = True
    if damage == 'terminal':
        with np.load(source/'checkpoint_latest.npz') as z:
            arrays = {k: z[k].copy() for k in z.files}
        meta = json.loads(str(arrays['metadata_json'].item())); meta['terminal_diagnostic'] = True
        arrays['metadata_json'] = np.array(json.dumps(meta)); E.atomic_npz(source/'checkpoint_latest.npz', **arrays)
        m['checkpoint_sha256'] = E.sha(source/'checkpoint_latest.npz')
    if damage == 'archive':
        with h5py.File(source/'data.h5', 'r+') as f:
            f['1d/delta_bar'][-1] = .5
    E.write_json(p, m)
    with pytest.raises(ValueError):
        E.prepare(source, out, final_time_ns=14.)


def test_gap_or_electrical_overlap_failures_do_not_clear_physics(native):
    source, out = native
    with h5py.File(source/'data.h5', 'r+') as f:
        d = f['1d/delta_bar'][:]; d[10:-1, 1] = .5; f['1d/delta_bar'][:] = d
    p = E.prepare(source, out, final_time_ns=14.)
    assert not p['admission']['admitted']
    with pytest.raises(ValueError, match='not admitted'):
        E.execute(out)


def test_exact_v123_catalog_is_supported(native):
    source, out = native; p = source/'manifest.json'; m = json.loads(p.read_text())
    m['model_version'] = '1.2.3'
    m['source_pins_sha256'] = {k: v for k, v in m['source_pins_sha256'].items() if k in E.NATIVE_SOURCES}
    with np.load(source/'checkpoint_latest.npz') as z:
        arrays = {k: z[k].copy() for k in z.files}
    meta = json.loads(str(arrays['metadata_json'].item())); meta['source_pins_sha256'] = m['source_pins_sha256']
    arrays['metadata_json'] = np.array(json.dumps(meta)); E.atomic_npz(source/'checkpoint_latest.npz', **arrays)
    m['checkpoint_sha256'] = E.sha(source/'checkpoint_latest.npz'); E.write_json(p, m)
    assert E.prepare(source, out, final_time_ns=14.)['admission']['admitted']


def test_forged_prepared_initial_state_is_refused(native):
    source, out = native
    plan = E.prepare(source, out, final_time_ns=14.)
    plan['initial_circuit_SI'][1] *= .9
    E.write_json(out/'phase_plan.json', plan)
    with pytest.raises(ValueError, match='state/ledger'):
        E.execute(out)


def test_single_cli_requires_an_admitted_execution_cpu():
    import importlib.util
    path = E.ROOT/'pipelines/25_v125_electrical_phase.py'
    spec = importlib.util.spec_from_file_location('_v125_cli_guard', path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    with pytest.raises(SystemExit) as caught:
        module.main(['--from-run', '/unused', '--output-dir', '/unused2', '--execute'])
    assert caught.value.code == 2
