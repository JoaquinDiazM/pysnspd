import numpy as np

from pysnspd.gtdgl.pytdgl_like.currents import _rms


def test_rms_helper_zero_and_nonzero():
    assert _rms(np.zeros(4)) == 0.0
    assert np.isclose(_rms(np.array([3.0, 4.0])), np.sqrt(12.5))
