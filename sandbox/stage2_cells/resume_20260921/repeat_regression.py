"""Fresh repository regression; preserve the previous stage-2 test record."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'sandbox/stage2_cells'))
import run_regression
if __name__=='__main__':
    run_regression.OUT=ROOT/'docs/implementation/stage2/resume_20260921'
    run_regression.main()
