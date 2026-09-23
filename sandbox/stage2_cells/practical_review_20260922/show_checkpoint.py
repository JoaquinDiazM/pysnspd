"""Show the received checkpoint, optionally checking preserved bytes; no RHS."""
from pathlib import Path
import argparse,hashlib,json
ROOT=Path(__file__).resolve().parents[3];DATA=ROOT/'docs/implementation/stage2/practical_review_20260922'
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--verify',action='store_true');args=parser.parse_args()
data=json.loads((DATA/'saved_results_audit.json').read_text())
if args.verify:
    for name,digest in data['raw_file_sha256'].items():
        assert hashlib.sha256((DATA/'raw'/name).read_bytes()).hexdigest()==digest,name
    for name,digest in data['output_sha256'].items():
        assert hashlib.sha256((DATA/name).read_bytes()).hexdigest()==digest,name
pair=json.loads((DATA/'raw/two_guarded_ph2049_pair_320_640.json').read_text())
print(json.dumps(dict(completed_tasks=data['completed_task_count'],trajectories=data['trajectory_count'],
    all_trajectory_invariants_pass=all(r['invariants_pass'] for r in data['trajectory_checks']),
    pair_status=pair['status'],pair_percent=100*pair['maximum_error'],auxiliary_budget_percent=100*pair['tolerance'],
    numerical_stage2_admission=False,static3a='Separate preparation allowed; not implemented',
    raw_hashes_checked=args.verify,new_long_run_required_now=False),indent=2))
