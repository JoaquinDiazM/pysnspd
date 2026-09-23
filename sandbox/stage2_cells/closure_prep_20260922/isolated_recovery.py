"""Recover a serializer failure without repeating its completed trajectory.

The original adapter, SSP implementation, partial JSON and NPZ stay immutable.
Only NumPy scalar serialization, provenance and reuse/output handling differ.
"""
from pathlib import Path
import argparse
import json
import time
import numpy as np
import isolated_checks as legacy

ORIGINAL_OUT = legacy.OUT
OUT = ORIGINAL_OUT/'isolated_recovery'
PLAN = OUT/'isolated_registration.json'
original_sources = legacy.sources
original_registration = legacy.registration


def scalar(value):
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError('only NumPy scalars receive additional serialization: '+type(value).__name__)


def write_new(path, value):
    # Finish validation before opening a file, so a serialization failure does
    # not leave a misleading truncated JSON result. No numerical value changes.
    encoded=json.dumps(value, indent=2, allow_nan=False, default=scalar)+'\n'
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x', encoding='utf-8') as f:
        f.write(encoded)


def sources():
    values=original_sources()
    values[Path(__file__).relative_to(legacy.ROOT).as_posix()]=legacy.sha(__file__)
    return values


def registration():
    result=original_registration()
    result['schema']='pysnspd.stage2.isolated-serializer-recovery.v1'
    result['isolated']['cases']=['absorption','recombination','creation','transport']
    result['recovery']=dict(reason='Original run stopped after emission integration on np.bool_ JSON serialization.',
        original_artifacts={p.relative_to(legacy.ROOT).as_posix():legacy.sha(p) for p in
            (ORIGINAL_OUT/'isolated_registration.json', ORIGINAL_OUT/'isolated_emission.npz',
             ORIGINAL_OUT/'isolated_emission.json', ORIGINAL_OUT/'isolated_results.json',
             ORIGINAL_OUT/'isolated_run.log')},
        emission_policy='Recompute state invariants from saved NPZ only; do not integrate or invent lost limiter statistics.',
        adapter_changes=['JSON default converts np.generic scalars with .item().',
                         'Serialize completely before exclusive file creation.',
                         'Separate output directory and preregistration.',
                         'Reuse emission and run only cases never attempted.'],
        numerical_adapter_and_SSP_unchanged=True,
        no_prior_timeout_or_physical_failure=True)
    return result


legacy.sources=sources
legacy.registration=registration
legacy.write_new=write_new
legacy.OUT=OUT
legacy.PLAN=PLAN


def recover_emission():
    started=time.perf_counter()
    previous=json.loads((ORIGINAL_OUT/'isolated_registration.json').read_text(encoding='utf-8'))
    for name, value in previous['source_hashes'].items():
        if legacy.sha(legacy.ROOT/name)!=value:
            raise ValueError('emission source changed: '+name)
    source=ORIGINAL_OUT/'isolated_emission.npz'
    with np.load(source,allow_pickle=False) as saved:
        arrays={name:np.array(saved[name],copy=True) for name in saved.files}
    system, _, _=legacy.setup('one', heating=0., escape=np.inf)
    for name, expected in (('electron_count',system.catalog.count_nodes),
                           ('electron_weights',system.catalog.count_weights),
                           ('phonon_energies',system.phonons.energies),
                           ('phonon_capacities',system.phonons.capacities)):
        if not np.array_equal(arrays[name],expected):
            raise ValueError('saved emission grid mismatch')
    times, states=arrays['times'], arrays['states']
    if not np.array_equal(times,np.linspace(0.,.02,21)):
        raise ValueError('saved emission times differ from registration')
    band=((system.catalog.count_nodes>1.1)&(system.catalog.count_nodes<1.5)).astype(float)
    if not np.array_equal(states[0],system.pack([.5],[band],[np.zeros(1025)])):
        raise ValueError('saved emission initial state differs from registration')
    rows=[legacy.measures(system,state) for state in states]
    a0,p0,n0=system.unpack(states[0])
    scale=max(1.,abs(rows[0]['total']))
    energy=max(abs(r['total']-rows[0]['total'])/scale for r in rows)
    count=max(abs(r['electron_count']-rows[0]['electron_count']) for r in rows)/max(1.,rows[0]['electron_count'])
    eph=max(abs(r['phonon_energy']-rows[0]['phonon_energy']-r['eph'])/scale for r in rows)
    no_other_channels=all(r['external']==r['heat']==r['escape']==r['transport']==0. for r in rows)
    populations=all(np.all((system.unpack(s)[1]>=0)&(system.unpack(s)[1]<=1))
                    and np.all(system.unpack(s)[2]>=0) for s in states)
    fields=all(np.array_equal(system.unpack(s)[0],a0) for s in states)
    delta_ph=rows[-1]['phonon_energy']-rows[0]['phonon_energy']
    gates=dict(energy=energy<=1e-7,count=count<=1e-12,shared_energy_ledger=eph<=1e-12,
               physical_populations=populations,fields_fixed=fields,
               expected_emission=delta_ph>100*np.finfo(float).eps,
               disabled_channels=no_other_channels,sample_count=len(times)>=10)
    value=dict(schema='pysnspd.stage2.isolated-emission-recovered.v1',
        status='PASS_SAVED_STATE_INVARIANTS_ONLY' if all(gates.values()) else 'FAIL',
        gates=gates,energy_error=energy,count_error=count,eph_ledger_error=eph,
        initial=rows[0],final=rows[-1],samples=rows,times=times.tolist(),
        source_trajectory=source.relative_to(legacy.ROOT).as_posix(),trajectory_sha256=legacy.sha(source),
        original_registration_sha256=legacy.sha(ORIGINAL_OUT/'isolated_registration.json'),
        registration_sha256=legacy.sha(PLAN),source_hashes=previous['source_hashes'],
        reintegrated=False,new_RHS_evaluations=0,
        lost_statistics=['integrator limiter counters','instantaneous residual and event-geometry numerical values'],
        limitation='Original partial JSON is preserved but is not valid evidence of the unrecorded statistics. Recovered verdict checks all saved states only.',
        runtime_seconds=time.perf_counter()-started)
    write_new(OUT/'isolated_emission_recovered.json',value)
    print(json.dumps(dict(case='emission_reused',status=value['status'],energy_error=energy,count_error=count)),flush=True)
    if not all(gates.values()):
        raise ValueError('saved emission invariant failed')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--register',action='store_true')
    args=parser.parse_args()
    if args.register:
        write_new(PLAN,registration())
        print('PREREGISTERED',PLAN)
        return
    if json.loads(PLAN.read_text(encoding='utf-8'))!=registration():
        raise ValueError('recovery preregistration/source mismatch')
    if list(OUT.glob('isolated_*.npz')) or (OUT/'isolated_results.json').exists() or (OUT/'isolated_emission_recovered.json').exists():
        raise FileExistsError('recovery outputs already exist; no overwrite or reexecution')
    recover_emission()
    legacy.run()


if __name__=='__main__':
    main()
