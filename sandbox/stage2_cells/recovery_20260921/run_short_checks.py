"""Three small new SSP trajectories; the caller caps this whole script at 240 s."""
from pathlib import Path
import subprocess,sys
ROOT=Path(__file__).resolve().parents[3]
DATA=ROOT/'docs/implementation/stage2/recovery_20260921'
def main():
    paths=[]
    for steps in (10,20,40):
        output=DATA/f'one_ssp_short_{steps}.json'
        log=output.with_suffix('.console.log')
        if log.exists():raise ValueError('Do not overwrite short-check evidence')
        with log.open('x') as stream:
            subprocess.run([sys.executable,'-u',str(Path(__file__).with_name('run_limited_coupled.py')),
                '--case','one','--steps',str(steps),'--duration','0.1','--escape','inf',
                '--heating','0.0','--method','rk4','--output',str(output)],
                cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT,check=True)
        paths.append(str(output))
    subprocess.run([sys.executable,str(Path(__file__).with_name('assess_limited_time.py')),
        *paths,'--reference',str(ROOT/'docs/implementation/stage2/resume_20260921/manual_reference/one_reference_probe.json'),
        '--output',str(DATA/'limited_short_time_assessment.json')],cwd=ROOT,check=True)
if __name__=='__main__':main()
