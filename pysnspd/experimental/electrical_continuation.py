"""Authenticated, append-only v1.2.3/v1.2.4 1D -> v1.2.5 0D handoff.

Native files are read only. The new archive links their committed history and
stores 0D in its own group. Thermal arrays remain in the native checkpoint.
"""
import contextlib
from dataclasses import asdict
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import time

import h5py
import numpy as np
from .electrical_0d import AdmissionConfig, InductiveCircuit, qualify

ROOT = Path(__file__).resolve().parents[2]
SCRATCH = Path('/home/jdiaz/scratch')
NATIVE_SOURCES = frozenset('''pysnspd/experimental/effective_bcs_functional.py
pysnspd/experimental/effective_bcs_preparation.py
pysnspd/experimental/effective_bcs_kinetics.py
pysnspd/experimental/effective_bcs_ports.py
pysnspd/experimental/effective_coupled.py
pysnspd/experimental/effective_material.py
pysnspd/experimental/electrical_ports.py
pysnspd/experimental/thermal_spatial_usadel.py
pysnspd/experimental/kinetic_events.py
pysnspd/experimental/trigger_window.py
pysnspd/experimental/effective_bcs_kinetics_numba.py
pipelines/05_effective_photon_run.py
sandbox/stage4_core/parallel_runtime.py
pipelines/09_thermal_stiffness_photon_run.py
pysnspd/experimental/depairing_calibration.py
pysnspd/experimental/thermal_stiffness.py
pysnspd/experimental/thermal_coupled.py
pysnspd/experimental/general_thermal_gradient.py
pysnspd/experimental/finite_depairing_sums.py
pysnspd/experimental/resummed_thermal.py
pipelines/13_dimensional_photon_run.py
pysnspd/experimental/dimensional_transition.py
pysnspd/experimental/dimensional_archive.py
pysnspd/experimental/dimensional_geometry.py
pysnspd/experimental/dimensional_monitors.py
pysnspd/experimental/dimensional_recovery.py
pysnspd/experimental/dimensional_1d.py
pysnspd/experimental/dimensional_cold.py'''.split())
V124_SOURCES = frozenset(('pipelines/18_spatial_photon_run.py',
    'pysnspd/experimental/latency_continuation.py', 'pysnspd/experimental/reflected_photon.py'))
PHASE_SOURCES = ('pysnspd/experimental/electrical_0d.py',
                 'pysnspd/experimental/electrical_continuation.py',
                 'pipelines/25_v125_electrical_phase.py')
FIELDS = ('time_ps', 'Ib_A', 'Is_A', 'vc_V', 'Vout_V', 'Vdev_V',
          'delta_Ib_A', 'delta_Is_A', 'delta_vc_V', 'energy_J', 'external_J',
          'energy_defect_J', 'phase_work_J', 'phase_energy_defect_J',
          'electrical_distance', 'thermal_evolved')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(2**20), b''):
            h.update(block)
    return h.hexdigest()


def write_json(path, value):
    p = Path(path); tmp = p.with_name(p.name+'.tmp')
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')
    os.replace(tmp, p)


def atomic_npz(path, **arrays):
    path = Path(path); tmp = path.with_name(path.name+'.tmp')
    with tmp.open('wb') as f:
        np.savez_compressed(f, **arrays)
    os.replace(tmp, path)


def scratch(path):
    p = Path(path).expanduser().resolve()
    if not p.is_relative_to(SCRATCH.resolve()):
        raise ValueError('All phase products must remain in scratch')
    return p


def phase_pins():
    return {k: sha(ROOT/k) for k in PHASE_SOURCES}


def refuse_live(source):
    for p in Path('/proc').iterdir():
        if not p.name.isdigit() or int(p.name) == os.getpid():
            continue
        try:
            words = (p/'cmdline').read_bytes().decode(errors='replace').split('\0')
            if '--output-dir' in words and Path(words[words.index('--output-dir')+1]).resolve() == source:
                raise ValueError('Native source has an active writer; preserved')
        except OSError:
            pass


@contextlib.contextmanager
def source_lock(source):
    # Multiple 0D readers may coexist; every existing campaign writer takes EX.
    path = source.parent.parent/'.execution.lock'
    with path.open('a+') as f:
        try:
            fcntl.flock(f, fcntl.LOCK_SH | fcntl.LOCK_NB)
        except BlockingIOError:
            raise ValueError('Native campaign is active; source preserved')
        refuse_live(source)
        yield


def authenticate(source, config=AdmissionConfig()):
    source = scratch(source); raw = (source/'manifest.json').read_bytes(); m = json.loads(raw)
    expected = NATIVE_SOURCES | (V124_SOURCES if m.get('model_version') == '1.2.4' else set())
    if m.get('model_version') not in ('1.2.3', '1.2.4') or set(m.get('source_pins_sha256', {})) != expected:
        raise ValueError('Only exact v1.2.3/v1.2.4 native producer catalogs are compatible')
    if any(sha(ROOT/k) != v for k, v in m['source_pins_sha256'].items()):
        raise ValueError('Native producer source changed')
    if (m.get('phase') != '1d' or m.get('error') is not None or
        m.get('status') not in ('WALL_BUDGET_CHECKPOINTED', 'INTERRUPTED_CHECKPOINTED',
                               'DEVELOPMENT_DURATION_REACHED', 'FINAL_TIME_REACHED', 'RECOVERY_CONFIRMED') or
        (m.get('cold_recovery') or {}).get('active')):
        raise ValueError('Intact paused/completed native1D required; failures/cold reconstructions preserved')
    if sha(source/'checkpoint_latest.npz') != m['checkpoint_sha256'] or sha(source/'dc1d.npz') != m['dc1d_sha256']:
        raise ValueError('Native checkpoint/DC checksum mismatch')
    with np.load(source/'checkpoint_latest.npz', allow_pickle=False) as z:
        arrays = {k: z[k].copy() for k in z.files}
    meta = json.loads(str(arrays['metadata_json'].item()))
    for key in ('phase', 'time_ps', 'source_pins_sha256', 'parameter_identity', 'archive_counts'):
        if meta[key] != m[key]:
            raise ValueError('Native manifest/checkpoint mismatch: '+key)
    if meta.get('terminal_diagnostic'):
        raise ValueError('Native terminal diagnostic cannot be promoted')
    with h5py.File(source/'data.h5', 'r', swmr=True) as f:
        for phase, sizes in m['archive_counts'].items():
            for key, n in sizes.items():
                if len(f[phase][key]) < n:
                    raise ValueError('Archive shorter than committed counts')
        n = m['archive_counts']['1d']['snapshot_time_ps']; q = f['1d']
        st = q['snapshot_time_ps'][:n]
        for key in ('delta_bar', 'p', 'n'):
            if not np.array_equal(q[key][n-1], arrays[key]):
                raise ValueError('Final archive/native checkpoint differs: '+key)
        if abs(st[-1]-meta['time_ps']) > 1e-8:
            raise ValueError('Native final saved snapshot is missing')
        tr = q['trace'][:m['archive_counts']['1d']['trace']]
        if (abs(tr['time_ps'][-1]-meta['time_ps']) > 1e-8 or
            not np.array_equal(np.array([tr[k][-1] for k in ('Ib_A', 'Is_A', 'vc_V')]), arrays['circuit_SI'])):
            raise ValueError('Final accepted circuit/native checkpoint mismatch')
        if np.any(tr['regime_id'] != 0):
            raise ValueError('Only full-native 1D evidence may admit the electrical phase')
        first = max(0, int(np.searchsorted(st, meta['time_ps']-1000*config.overlap_ns, side='right'))-1)
        gap = q['delta_bar'][first:n]; dc = q['dc_delta_bar'][:]
        partition = m['inductance_partition']
        if not partition.get('added_96nH_preserved') or not partition.get('partition_frozen_from_DC'):
            raise ValueError('Authenticated frozen DC inductance partition required')
        Lfilm = partition['bulk_differential_H_per_m']*partition['length_1d_m']
        circuit = InductiveCircuit(m['circuit_by_phase']['1d'], Lfilm)
        admission = qualify(circuit, tr, st[first:n], gap, dc, m, config)
        signature = dict(manifest_sha256=hashlib.sha256(raw).hexdigest(),
            checkpoint_sha256=m['checkpoint_sha256'], dc1d_sha256=m['dc1d_sha256'],
            native_trace_sha256=hashlib.sha256(tr.tobytes()).hexdigest(),
            archive_counts=m['archive_counts'], native_archive_size=(source/'data.h5').stat().st_size)
        native_end = {k: float(tr[k][-1]) for k in tr.dtype.names}
    if (source/'manifest.json').read_bytes() != raw:
        raise ValueError('Native source changed during authentication')
    return m, meta, arrays, circuit, admission, signature, native_end


def prepare(source, destination, final_time_ns=100., sample_ps=10., config=AdmissionConfig()):
    source, destination = scratch(source), scratch(destination)
    if destination == source or destination.is_relative_to(source):
        raise ValueError('Separate 0D phase directory required')
    if not np.isfinite(final_time_ns) or not np.isfinite(sample_ps) or not 0 < sample_ps <= 50:
        raise ValueError('Finite horizon and positive sampling interval required')
    m, meta, arrays, circuit, admission, signature, end = authenticate(source, config)
    if final_time_ns*1000 <= meta['time_ps']:
        raise ValueError('0D horizon must be later than the native checkpoint')
    plan = dict(schema='pysnspd.electrical_phase_plan.v1', model_version='1.2.5',
        native_source_run=str(source), native_identity=m['parameter_identity'],
        native_signature=signature, native_source_pins_sha256=m['source_pins_sha256'],
        source_pins_sha256=phase_pins(), final_time_ns=final_time_ns, sample_ps=sample_ps,
        admission_config=asdict(config), admission=admission, native_end=end,
        initial_time_ps=meta['time_ps'], initial_circuit_SI=arrays['circuit_SI'].tolist(),
        circuit=m['circuit_by_phase']['1d'], film_inductance_H=circuit.Lfilm,
        initial_energy_J=meta['initial_energy_J'], native_external_J=meta['external_J'],
        native_absolute_external_J=meta['absolute_external_J'],
        native_energy_defect_J=m['energy_defect_J'],
        native_retained_energy_J=m['source']['retained_energy_J'],
        native_global_energy_tolerance=m['parameter_identity']['global_energy_tolerance'],
        native_local_energy_tolerance=m['parameter_identity']['local_energy_tolerance'],
        frozen_thermal_state=dict(checkpoint_sha256=m['checkpoint_sha256'],
            time_ps=meta['time_ps'], latest_monitor=m['latest_monitor'],
            thermal_audit_pending=True, thermal_reset_certified=False))
    destination.mkdir(parents=True, exist_ok=True)
    p = destination/'phase_plan.json'
    if p.exists() and json.loads(p.read_text()) != plan:
        raise ValueError('Existing 0D plan differs; preserve it and choose a fresh output')
    write_json(p, plan)
    return plan


def execute(destination, *, resume=False, stop_requested=lambda: False):
    """Execute only an admitted prepared plan. Checkpoints commit flushed prefixes."""
    out = scratch(destination); plan = json.loads((out/'phase_plan.json').read_text())
    source = Path(plan['native_source_run'])
    if phase_pins() != plan['source_pins_sha256']:
        raise ValueError('0D evolution sources changed since preparation')
    if not plan['admission']['admitted']:
        raise ValueError('Electrical 0D handoff not admitted: '+str(plan['admission']['reasons']))
    with source_lock(source), (out/'.execution.lock').open('a+') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise ValueError('This 0D output already has a writer')
        # Re-authenticate immediately before any write to the new phase archive.
        verified = authenticate(source, AdmissionConfig(**plan['admission_config']))
        if verified[5] != plan['native_signature']:
            raise ValueError('Native source changed after preparation')
        inherited = dict(initial_time_ps=verified[1]['time_ps'],
            initial_circuit_SI=verified[2]['circuit_SI'].tolist(), native_end=verified[6],
            initial_energy_J=verified[1]['initial_energy_J'],
            native_external_J=verified[1]['external_J'],
            native_absolute_external_J=verified[1]['absolute_external_J'],
            native_energy_defect_J=verified[0]['energy_defect_J'],
            native_retained_energy_J=verified[0]['source']['retained_energy_J'],
            native_identity=verified[0]['parameter_identity'])
        if any(plan.get(k) != v for k, v in inherited.items()):
            raise ValueError('Prepared native state/ledger was changed')
        circuit = verified[3]; z0 = circuit.initialize(plan['initial_circuit_SI'])
        initial_deviation = circuit.energy_deviation(z0)
        checkpoint = out/'checkpoint_latest.npz'; new = not checkpoint.exists()
        if not new:
            if not resume:
                raise ValueError('Existing 0D checkpoint: use --resume')
            m = json.loads((out/'manifest.json').read_text())
            if sha(checkpoint) != m['checkpoint_sha256']:
                raise ValueError('0D checkpoint checksum mismatch')
            if m['status'] == 'ELECTRICAL_HORIZON_COMPLETE':
                with h5py.File(out/'data.h5', 'r', swmr=True) as archive:
                    count = m['archive_counts']['0d']['trace']
                    if len(archive['0d/trace']) != count or archive['0d/trace'][-1]['time_ps'] != m['time_ps']:
                        raise ValueError('Completed 0D archive integrity mismatch')
                return m
            if m['status'] not in ('RUNNING', 'INTERRUPTED_CHECKPOINTED'):
                raise ValueError('Diagnostic 0D stop preserved; no blind resume')
            with np.load(checkpoint, allow_pickle=False) as f:
                z = f['lifted_state'].copy(); meta = json.loads(str(f['metadata_json'].item()))
            if meta['plan_sha256'] != sha(out/'phase_plan.json') or meta['time_ps'] != m['time_ps']:
                raise ValueError('0D checkpoint plan/time mismatch')
            t = meta['time_ps']; count = meta['trace_count']
            if count != m['archive_counts']['0d']['trace']:
                raise ValueError('0D checkpoint/manifest trace counts differ')
        else:
            if (out/'data.h5').exists():
                raise ValueError('Uncommitted output archive preserved')
            shutil.copy2(source/'manifest.json', out/'native_manifest.json')
            shutil.copy2(source/'checkpoint_latest.npz', out/'native_checkpoint.npz')
            shutil.copy2(source/'dc1d.npz', out/'dc1d.npz')
            z = z0; t = plan['initial_time_ps']; count = 0
            meta = dict(plan_sha256=sha(out/'phase_plan.json'), phase='0d', time_ps=t,
                        trace_count=0, max_phase_energy_fraction=0.)
            m = dict(schema='pysnspd.electrical_phase.v1', model_version='1.2.5',
                physical_kernel_version=verified[0]['physical_kernel_version'], phase='0d',
                run_complete=False, status='RUNNING', source_pins_sha256=plan['source_pins_sha256'],
                parameter_identity=plan['native_identity'], native_source_pins_sha256=plan['native_source_pins_sha256'],
                native_reference=plan['native_signature'], native_source_run=str(source),
                admission=plan['admission'], frozen_thermal_state=plan['frozen_thermal_state'],
                initial_time_ps=t, target_time_ps=plan['final_time_ns']*1000,
                circuit_by_phase={'0d': dict(circuit.c, resolved_film_inductance_H=circuit.Lfilm,
                                           total_detector_inductance_H=circuit.Ltotal)},
                source=verified[0]['source'], observation=verified[0]['observation'],
                energy_scope='Frozen unresolved energy offset plus lumped electrical energy; native thermal ledger ends at handoff',
                thermal_evolved=False, thermal_reset_certified=False,
                provenance=dict(native_phase_unmodified=True, photon_redeposited=False,
                    native_arrays_preserved=True, native_history_linked_read_only=True))
        f = h5py.File(out/'data.h5', 'w' if new else 'r+', libver='latest')
        try:
            if new:
                f.attrs.update(schema='pysnspd.electrical_phase.v1', model_version='1.2.5',
                    thermal_policy='No native field reconstruction or thermal extrapolation after 1D handoff')
                for phase in plan['native_signature']['archive_counts']:
                    f[phase] = h5py.ExternalLink(str(source/'data.h5'), '/'+phase)
                q = f.create_group('0d')
                q.attrs['start_time_ps'] = t
                ds = q.create_dataset('trace', shape=(0,), maxshape=(None,),
                    dtype=[(k, 'f8') for k in FIELDS], chunks=(256,))
            else:
                ds = f['0d/trace']
                if len(ds) < count:
                    raise ValueError('0D archive shorter than committed checkpoint')
                ds.resize(count, axis=0)  # discard only this phase's uncommitted tail
                if count and abs(ds[-1]['time_ps']-t) > 1e-8:
                    raise ValueError('0D archive/checkpoint timestamp mismatch')
            f.flush(); f.swmr_mode = True

            def row():
                nonlocal count
                o = circuit.outputs(z); de = circuit.energy_deviation(z)-initial_deviation
                defect = de-z[12]
                scale = plan['native_retained_energy_J'] + plan['native_absolute_external_J'] + abs(z[12])
                fraction = abs(defect)/scale
                meta['max_phase_energy_fraction'] = max(meta['max_phase_energy_fraction'], fraction)
                if fraction > plan['native_local_energy_tolerance']:
                    raise ArithmeticError('Independent lumped-phase energy ledger failed')
                if abs(plan['native_energy_defect_J']+defect)/scale > plan['native_global_energy_tolerance']:
                    raise ArithmeticError('Inherited plus reduced global energy guard failed')
                if min(o['Ib_A'], o['Is_A']) <= 0:
                    raise ArithmeticError('0D global circuit polarity guard failed')
                limits = verified[0]['monitor_config']
                if max(o['Ib_A'], o['Is_A'])*1e6 > limits['current_limit_uA']:
                    raise ArithmeticError('0D inherited current guard failed')
                for key, limit in [('Vdev_V','device_voltage_limit_mV'),
                                   ('Vout_V','readout_voltage_limit_mV'),
                                   ('vc_V','capacitor_voltage_limit_mV')]:
                    if abs(o[key])*1000 > limits[limit]:
                        raise ArithmeticError('0D inherited voltage guard failed: '+key)
                values = dict(o, time_ps=t, delta_Ib_A=z[0], delta_Is_A=z[1], delta_vc_V=z[2],
                    energy_J=plan['native_end']['energy_J']+de,
                    external_J=plan['native_external_J']+z[12],
                    energy_defect_J=plan['native_energy_defect_J']+defect,
                    phase_work_J=z[12], phase_energy_defect_J=defect,
                    electrical_distance=circuit.electrical_distance(z), thermal_evolved=0.)
                ds.resize(count+1, axis=0); ds[count] = tuple(values[k] for k in FIELDS); count += 1
                m['latest_electrical'] = {k: float(v) for k, v in values.items()}

            def save(status):
                f.flush(); meta.update(time_ps=t, trace_count=count)
                atomic_npz(checkpoint, lifted_state=z, circuit_SI=circuit.state(z),
                           circuit_delta_SI=z[:3], metadata_json=np.array(json.dumps(meta)))
                m.update(status=status, run_complete=status == 'ELECTRICAL_HORIZON_COMPLETE',
                    time_ps=t, checkpoint_sha256=sha(checkpoint), error=None,
                    archive_counts=dict(plan['native_signature']['archive_counts'], **{'0d': {'trace': count}}),
                    max_phase_energy_fraction=meta['max_phase_energy_fraction'])
                write_json(out/'manifest.json', m)

            if new:
                row(); save('RUNNING')
            target = plan['final_time_ns']*1000; next_commit = t+1000
            try:
                while t < target-1e-8:
                    if stop_requested():
                        save('INTERRUPTED_CHECKPOINTED'); return m
                    h = min(plan['sample_ps'], target-t)
                    z = circuit.advance(z, h); t += h; row()
                    if t >= next_commit:
                        save('RUNNING'); next_commit = t+1000
                save('ELECTRICAL_HORIZON_COMPLETE')
            except Exception as exc:
                # The last committed checkpoint stays resumable; diagnostics are explicit.
                m.update(status='ELECTRICAL_DIAGNOSTIC_STOP', run_complete=False,
                         error=type(exc).__name__+': '+str(exc))
                write_json(out/'manifest.json', m)
                raise
        finally:
            f.close()
        if sha(source/'checkpoint_latest.npz') != plan['native_signature']['checkpoint_sha256']:
            raise ValueError('Native source changed during electrical continuation')
        return m
