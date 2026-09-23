"""Print the final development-stage closure; optionally verify saved evidence.

No numerical kernels are imported and no trajectories are executed.
"""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT/'docs/implementation/stage2/closure_20260922'
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--verify', action='store_true')
args = parser.parse_args()
decision = json.loads((DATA/'closure_decision.json').read_text(encoding='utf-8'))
admission = json.loads((ROOT/'docs/implementation/stage2/stage2_admission.json').read_text(encoding='utf-8'))
for key in ('development_stage_closed', 'accepted_for_stage3_development',
            'numerical_admission', 'production_promotion'):
    assert decision[key] == admission[key], 'Conflicting active decisions: '+key
if args.verify:
    subprocess.run([sys.executable, str(ROOT/'sandbox/stage2_cells/practical_review_20260922/show_checkpoint.py'),
                    '--verify'], check=True, stdout=subprocess.DEVNULL)
    subprocess.run([sys.executable, str(ROOT/'sandbox/stage2_cells/verify_delivery.py')],
                    check=True, stdout=subprocess.DEVNULL)
    for record in decision['evidence']:
        assert hashlib.sha256((ROOT/record['path']).read_bytes()).hexdigest()==record['sha256'],record['path']
print(json.dumps(dict(status=decision['status'],
    development_stage_closed=decision['development_stage_closed'],
    accepted_for_stage3_development=decision['accepted_for_stage3_development'],
    strict_numerical_certificate=decision['strict_certificate']['status'],
    production_promotion=decision['production_promotion'],
    trajectories=decision['evidence_summary']['completed_trajectories'],
    completed_tasks=decision['evidence_summary']['completed_tasks'],
    auxiliary_pair_status=decision['evidence_summary']['auxiliary_failure']['status'],
    hashes_verified=args.verify, new_rhs_evaluations=0,
    final_report='output/pdf/implementation/Informe_cierre_etapa_2.pdf',
    next_stage='docs/implementation/stage3/entry_contract.json',
    new_long_run_required_now=False), indent=2))
