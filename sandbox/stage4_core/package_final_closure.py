"""Extract the compact stage-4 review set; never rerun a physical calculation."""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DEST = ROOT / 'docs/implementation/stage4/closure_20260924/data'
REMOTE = '/home/jdiaz/scratch/stage4_final_coupling_20260924'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    args = parser.parse_args()
    source = args.source.resolve()
    if source == DEST.resolve():
        raise ValueError('Raw source must differ from compact destination')
    records = []
    for path in sorted(source.rglob('*')):
        if not path.is_file():
            continue
        relative = path.relative_to(source)
        record = dict(path=relative.as_posix(), source_sha256=digest(path),
                      source_bytes=path.stat().st_size,
                      raw_path=REMOTE + '/' + relative.as_posix())
        target = DEST / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        if path.name == 'reference.npz':
            keys = ['delta_bar', 'coordinates_bar', 'edges', 'conductance',
                    'area_weights', 'boundary_nodes', 'phase_path_nodes', 'current_bar']
            with np.load(path) as data:
                np.savez_compressed(target, **{key: data[key] for key in keys})
            record.update(storage='review_subset_not_solver_restart', retained_arrays=keys)
        elif path.name == 'summary.json':
            data = json.loads(path.read_text())
            metrics = {key: value for key, value in data['metrics'].items() if key != 'modes'}
            selected = {key: data[key] for key in
                        ('status', 'criterion_met', 'completed_sweeps', 'gap_reference_bar', 'elapsed_seconds')
                        if key in data}
            selected['metrics'] = metrics
            # No duplicated per-frequency or iteration logs in the review packet.
            for key in list(selected):
                if isinstance(selected[key], list) and len(selected[key]) > 20:
                    selected.pop(key)
            target.write_text(json.dumps(selected, indent=2) + '\n', encoding='utf8')
            record['storage'] = 'compact_summary'
        elif path.name == 'integrated_moments.npz':
            record['storage'] = 'raw_only_on_geminga'
        else:
            shutil.copy2(path, target)
            record['storage'] = 'exact_copy'
        if record['storage'] != 'raw_only_on_geminga':
            record.update(review_sha256=digest(target), review_bytes=target.stat().st_size)
        records.append(record)
    manifest = dict(schema='pysnspd.stage4.compact_evidence.v1',
                    source_commit='c80c0f8612d8fd5a2b390938ca5ed4b9d58d3caa',
                    raw_root=REMOTE, physics_solves_performed=0,
                    note='Compact NPZ references support plotting, not solver restart. Full stationary Matsubara fields, integrated response moments and original outputs remain in scratch on Geminga; individual real-energy perturbation fields were not saved.',
                    files=records)
    DEST.mkdir(parents=True, exist_ok=True)
    (DEST.parent / 'data_manifest.json').write_text(json.dumps(manifest, indent=2)+'\n', encoding='utf8')
    print(json.dumps(dict(files=len(records), review_bytes=sum(r.get('review_bytes', 0) for r in records))))


if __name__ == '__main__':
    main()
