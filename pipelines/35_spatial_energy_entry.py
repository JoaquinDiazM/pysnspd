"""Positive finite wavelength CLI extension over frozen native1.2.4 evolution.

The original13/18 entrypoints and native equations remain byte-identical.
New trajectories authenticate this adapter as an additional source pin.
"""
from __future__ import annotations
import argparse,importlib.util,json,math,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
ENTRY='pipelines/35_spatial_energy_entry.py'
spec=importlib.util.spec_from_file_location('_energy_frozen_spatial',ROOT/'pipelines/18_spatial_photon_run.py')
A=importlib.util.module_from_spec(spec);spec.loader.exec_module(A)
_native_configure=A.configure

def positive_finite(value):
    value=float(value)
    if not math.isfinite(value) or value<=0:raise argparse.ArgumentTypeError('Wavelength must be positive and finite')
    return value

def configure(source_runner,threads):
    R,prepare=_native_configure(source_runner,threads)
    old_parser=R.parser;old_pins=R.source_pins;old_validate=R.OLD.validate_args
    def parser():
        p=old_parser()
        action=p._option_string_actions['--wavelength-nm']
        action.choices=None;action.type=positive_finite
        return p
    def pins():
        return dict(old_pins(),**{ENTRY:R.B.sha256(ROOT/ENTRY)})
    def validate(args):
        positive_finite(args.wavelength_nm)
        # In the frozen validator wavelength is used only by its two-value enum.
        # Validate every other field on a copy with a legacy enum representative;
        # the real700nm value is never changed in the simulation or its identity.
        legacy=argparse.Namespace(**vars(args));legacy.wavelength_nm=775.
        old_validate(legacy)
        return args
    R.parser=parser;R.source_pins=pins;R.OLD.validate_args=validate
    return R,prepare

def _dispatch(argv,validate=False):
    previous=A.configure
    def selected(source_runner,threads):
        R,prepare=configure(source_runner,threads)
        if validate:
            def validate_run(args):
                R.OLD.validate_args(args)
                return dict(vars(args))
            R.run=validate_run
        return R,prepare
    A.configure=selected
    try:return A.main(argv)
    finally:A.configure=previous

def validate_arguments(argv):
    """Use the real argument/energy checks; never construct a physical system."""
    return _dispatch(list(argv),validate=True)

def main(argv=None):
    argv=list(sys.argv[1:] if argv is None else argv)
    if '--validate-only' in argv:
        argv.remove('--validate-only');a=validate_arguments(argv)
        print(json.dumps(dict(status='ARGUMENTS_VALIDATED',current_uA=a['current_uA'],
            wavelength_nm=a['wavelength_nm'],deposited_energy_eV=a['deposited_energy_eV'],
            retained_fraction=a['retained_fraction']),sort_keys=True))
        return 0
    return _dispatch(argv)

if __name__=='__main__':raise SystemExit(main())
