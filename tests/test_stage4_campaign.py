"""Fast checks of the manual campaign's no-duplicate-work and retention paths."""
from pathlib import Path
import importlib.util
import json
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('stage4_controls_test',ROOT/'sandbox/stage4_core/controls.py')
module=importlib.util.module_from_spec(spec)
import sys
sys.modules[spec.name]=module
spec.loader.exec_module(module)


class TestCampaignControls(unittest.TestCase):
    def test_nonanchor_skips_repeated_fd_but_keeps_full2D_diagnostics(self):
        plan=json.loads((ROOT/'docs/implementation/stage4/start_20260923/campaign_plan.json').read_text())
        case=dict(phase='spatial',profile='suppressed',population='vacuum',delta_reg=.2,
                  elements_x=1,elements_y=1,degree=2,finite_differences=False)
        cat=module.catalogue(plan,'reference')
        with tempfile.TemporaryDirectory() as folder:
            result=module.spatial_case(plan,case,cat,lambda *args:None,Path(folder))
            self.assertEqual(result['identity_checks'],[])
            self.assertTrue((Path(folder)/'fields.npz').is_file())
        self.assertEqual(result['nodes'],9)
        self.assertEqual(result['positive_symbols']+result['negative_symbols']+result['unresolved_symbols'],9)
        self.assertEqual(len(result['mobility']),3)
        self.assertEqual(cat.calls,0)  # Empty populations need no costly spectra.

if __name__=='__main__':unittest.main()
