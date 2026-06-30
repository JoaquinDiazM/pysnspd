from __future__ import annotations

import numpy as np

from pysnspd.gtdgl.pytdgl_like.usadel_current import _interp_nd_axiswise


def _slow_reference(table, delta_axis, q_axis, delta, q):
    out = []
    for d, qq in zip(delta, q):
        # Bilinear interpolation with q oddness for a positive q axis.
        qsign = np.sign(qq)
        qa = abs(float(qq))
        d = float(np.clip(d, delta_axis[0], delta_axis[-1]))
        qa = float(np.clip(qa, q_axis[0], q_axis[-1]))
        id1 = int(np.searchsorted(delta_axis, d, side="right"))
        iq1 = int(np.searchsorted(q_axis, qa, side="right"))
        id1 = min(max(id1, 1), len(delta_axis) - 1)
        iq1 = min(max(iq1, 1), len(q_axis) - 1)
        id0 = id1 - 1
        iq0 = iq1 - 1
        ad = (d - delta_axis[id0]) / (delta_axis[id1] - delta_axis[id0])
        aq = (qa - q_axis[iq0]) / (q_axis[iq1] - q_axis[iq0])
        val = (
            (1 - ad) * (1 - aq) * table[id0, iq0]
            + ad * (1 - aq) * table[id1, iq0]
            + (1 - ad) * aq * table[id0, iq1]
            + ad * aq * table[id1, iq1]
        )
        out.append(qsign * val)
    return np.asarray(out)


def test_vectorized_delta_q_interpolation_matches_reference_and_is_odd():
    delta_axis = np.array([1.0, 2.0, 3.0])
    q_axis = np.array([0.0, 10.0, 20.0, 30.0])
    table = np.array(
        [
            [0.0, 1.0, 2.0, 3.0],
            [0.0, 2.0, 4.0, 6.0],
            [0.0, 3.0, 6.0, 9.0],
        ],
        dtype=float,
    )
    delta = np.array([1.2, 1.7, 2.5, 2.9, 1.2])
    q = np.array([5.0, -12.0, 22.0, -29.0, 0.0])

    got = _interp_nd_axiswise(
        table,
        [("delta", delta_axis, delta), ("q", q_axis, q)],
    )
    expected = _slow_reference(table, delta_axis, q_axis, delta, q)
    assert np.allclose(got, expected)
    assert got[1] < 0.0
    assert got[3] < 0.0
    assert got[-1] == 0.0


def test_vectorized_interpolation_supports_one_point_delta_axis():
    delta_axis = np.array([2.0])
    q_axis = np.array([0.0, 10.0, 20.0])
    table = np.array([[0.0, 4.0, 8.0]])
    delta = np.array([2.0, 2.0, 2.0])
    q = np.array([-5.0, 0.0, 15.0])
    got = _interp_nd_axiswise(
        table,
        [("delta", delta_axis, delta), ("q", q_axis, q)],
    )
    assert np.allclose(got, [-2.0, 0.0, 6.0])
