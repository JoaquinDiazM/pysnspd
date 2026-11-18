"""v1.2.6 complete fast-mode photon run; frozen native physics, sparse outputs."""
from pathlib import Path
import importlib.util,json,sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
spec=importlib.util.spec_from_file_location('_v126_frozen_energy',ROOT/'pipelines/35_spatial_energy_entry.py')
BASE=importlib.util.module_from_spec(spec);spec.loader.exec_module(BASE)
MODEL=ROOT/'pysnspd/data/v126_fast_tail.json'

def configure(source_runner,threads):
    from pysnspd.experimental.fast_latency import TailModel,install
    R,prepare=BASE.configure(source_runner,threads)
    install(R,TailModel(MODEL),'pipelines/53_v126_fast_photon.py')
    return R,prepare

def dispatch(argv,validate=False):
    # Reject foreign checkpoints BEFORE entering18's legacy exception writer.
    # Otherwise even a deliberately refused resume would annotate its origin.
    if '--resume' in argv or '--recover-from-run' in argv:
        import argparse
        early=argparse.ArgumentParser(add_help=False);early.add_argument('--output-dir',type=Path);early.add_argument('--recover-from-run',type=Path)
        a,_=early.parse_known_args(list(argv))
        if a.recover_from_run is not None:raise ValueError('Normal-to-fast continuation is forbidden; no recover-from-run in fast mode')
        if a.output_dir is not None and (a.output_dir/'manifest.json').exists():
            from pysnspd.experimental.fast_latency import PROTOCOL
            saved=json.loads((a.output_dir/'manifest.json').read_text())
            if saved.get('execution_mode')!='fast' or saved.get('fast_protocol')!=PROTOCOL:
                raise ValueError('Normal-to-fast continuation is forbidden')
    previous=BASE.A.configure
    def selected(source_runner,threads):
        R,prepare=configure(source_runner,threads)
        if validate:
            def validate_only(a):
                R.OLD.validate_args(a);R.fast_validate(a);return dict(vars(a))
            R.run=validate_only
        return R,prepare
    BASE.A.configure=selected
    try:return BASE.A.main(list(argv))
    finally:BASE.A.configure=previous

def validate_arguments(argv):return dispatch(argv,validate=True)

def main(argv=None):
    argv=list(sys.argv[1:] if argv is None else argv)
    if '--validate-only' in argv:
        argv.remove('--validate-only');a=validate_arguments(argv)
        print(json.dumps(dict(status='FAST_ARGUMENTS_VALIDATED',execution_mode=a['execution_mode'],current_uA=a['current_uA'],energy_eV=a['deposited_energy_eV'],y_nm=a['impact_y_nm'])))
        return 0
    return dispatch(argv)

if __name__=='__main__':raise SystemExit(main())
