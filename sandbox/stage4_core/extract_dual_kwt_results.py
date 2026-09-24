"""Copy compact completed dual KWT results; never evolve a physical state."""
from pathlib import Path
import argparse,hashlib,json,shutil


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    summary=json.loads((a.source/'summary.json').read_text())
    if summary.get('completed_horizon_ps')!=1.0:
        raise ValueError('A completed 1 ps campaign is required')
    if a.output.exists():raise FileExistsError('Preserve existing extraction')
    names=['identity.json','executed_plan.json','summary.json','step_history.json','refinement.json','progress.jsonl']
    names+=sorted(q.name for q in a.source.glob('observation_*') if q.suffix in ('.json','.npz'))
    a.output.mkdir(parents=True)
    for name in names:shutil.copyfile(a.source/name,a.output/name)
    observations=[]
    for q in sorted(a.output.glob('observation_*.json')):
        row=json.loads(q.read_text());assert sha(a.output/row['fields_path'])==row['fields_sha256']
        observations.append(row['time_ps'])
    receipt=dict(source=str(a.source.resolve()),status='COMPLETED_RESULTS_EXTRACTED',
        physical_solves=0,observations_ps=observations,
        checkpoints_retained_remotely=len(list(a.source.glob('checkpoint_*.npz'))),
        source_hashes={name:sha(a.source/name) for name in names},
        extractor_sha256=sha(Path(__file__)))
    (a.output/'extraction_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(dict(status=receipt['status'],files=len(names),physical_solves=0)))


if __name__=='__main__':main()
