"""Independent finite-Gamma field checks on the new guarded equilibrium run."""
from pathlib import Path
import argparse,json,sys
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'tmp/stage1_r2_review/deps'),str(Path(__file__).parent)]
from assess_guarded_time import verify_contract,sha
from special_cases import equilibrium_fields


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('runs',nargs=1,type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();path=args.runs[0]
    if args.output.exists():raise ValueError('Preserve prior field assessment')
    record=json.loads(path.read_text());verify_contract({'record':record})
    if record['parameters']['scenario']!='equilibrium':raise ValueError('Equilibrium trajectory required')
    result=equilibrium_fields(path,args.output)
    print(json.dumps(result),flush=True)
    if result['status']!='PASS':raise SystemExit(1)


if __name__=='__main__':main()
