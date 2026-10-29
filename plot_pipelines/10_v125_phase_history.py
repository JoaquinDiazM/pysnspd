"""Inspect all committed 2D/1D/0D phases; never reconstruct native 0D fields."""
import argparse
import json
import os
from pathlib import Path
import sys
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pysnspd.experimental.electrical_continuation import sha, scratch


def read_history(run):
    run = scratch(run); m = json.loads((run/'manifest.json').read_text())
    if m.get('model_version') != '1.2.5' or m.get('phase') != '0d':
        raise ValueError('Use the native report entry for historical 2D/1D-only cases')
    plan = json.loads((run/'phase_plan.json').read_text())
    source = Path(plan['native_source_run'])
    for file, key in [('manifest.json', 'manifest_sha256'), ('checkpoint_latest.npz', 'checkpoint_sha256'),
                      ('dc1d.npz', 'dc1d_sha256')]:
        if sha(source/file) != plan['native_signature'][key]:
            raise ValueError('Linked native source changed; committed history must be reauthenticated')
    if sha(run/'checkpoint_latest.npz') != m['checkpoint_sha256']:
        raise ValueError('0D checkpoint changed')
    trace = {}; temperatures = []
    with h5py.File(run/'data.h5', 'r', swmr=True) as f:
        for phase, sizes in m['archive_counts'].items():
            trace[phase] = f[phase]['trace'][:sizes['trace']]
            if phase in ('2d', '1d'):
                n = sizes['snapshot_time_ps']
                temperatures.append((f[phase]['snapshot_time_ps'][:n], f[phase]['Te_K'][:n].max(axis=1)))
    return m, trace, temperatures


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-dir', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(argv); output = scratch(a.output); output.parent.mkdir(parents=True, exist_ok=True)
    m, traces, temps = read_history(a.run_dir)
    fig, ax = plt.subplots(3, 2, figsize=(11.7, 8.3), layout='constrained')
    for phase, color in [('2d', '#174e75'), ('1d', '#398e76'), ('0d', '#df672b')]:
        if phase not in traces:
            continue
        tr = traces[phase]; t = tr['time_ps']/1000
        for panel, key, scale, label in [(ax[0, 0], 'Vout_V', 1e3, 'Intrinsic Vout (mV)'),
                (ax[0, 1], 'Is_A', 1e6, 'Detector current Is (uA)'),
                (ax[1, 0], 'vc_V', 1e3, 'Capacitor voltage (mV)'),
                (ax[1, 1], 'Vdev_V', 1e6, 'Film terminal voltage (uV)'),
                (ax[2, 1], 'energy_defect_J', 1., 'Inherited + phase ledger defect (J)')]:
            panel.plot(t, tr[key]*scale, color=color, label=phase.upper())
            panel.set_ylabel(label)
    for t, T in temps:
        ax[2, 0].plot(t/1000, T, color='#174e75')
    ax[2, 0].set_ylabel('Saved native maximum effective Te (K)')
    ax[2, 0].set_title('No thermal evolution after the 1D endpoint')
    ax[1, 1].set_yscale('symlog', linthresh=1.)
    start = m['initial_time_ps']/1000
    for panel in ax.ravel():
        panel.axvline(start, ls=':', color='gray'); panel.grid(alpha=.2)
        panel.set_xlabel('Global time since 2D start (ns)')
    ax[0, 0].legend(); ax[2, 1].legend()
    fig.suptitle('pysnspd 1.2.5: native history + admitted electrical-only 0D phase\n'
                 'Frozen DC film inductance retained; thermal audit pending')
    fig.savefig(output); plt.close(fig)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
