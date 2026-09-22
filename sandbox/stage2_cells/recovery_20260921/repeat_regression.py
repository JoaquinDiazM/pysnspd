"""Fresh regression record for the recovery; preserve previous test evidence."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'sandbox/stage2_cells'))
import run_regression
if __name__=='__main__':
    run_regression.OUT=ROOT/'docs/implementation/stage2/recovery_20260921'
    if (run_regression.OUT/'regression_result.json').exists():raise ValueError('Preserve prior regression evidence')
    run_regression.main()
