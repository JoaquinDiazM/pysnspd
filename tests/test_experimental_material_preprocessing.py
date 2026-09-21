import hashlib
import numpy as np
import pytest
from pysnspd.experimental.material_preprocessing import derive_phonon_shape


def write(tmp_path, rows):
    path = tmp_path / "source.dat"
    np.savetxt(path, rows)
    return path


def derive(path, **kwargs):
    return derive_phonon_shape(path, source_url="https://example.org/primary.dat",
                              source_revision="fixed-test-revision", **kwargs)


def test_derivation_keeps_stable_first_duplicates_and_exact_source(tmp_path):
    path = write(tmp_path, [[2, 7, -2], [1, 3, 4], [0, 1e-40, 1e-20], [1, 9, 8], [3, 2, 0]])
    raw = path.read_bytes()
    result = derive(path)
    np.testing.assert_array_equal(result.axis, [0, 1, 2, 3])
    np.testing.assert_array_equal(result.alpha2F, [0, 3, 0, 0])
    np.testing.assert_array_equal(result.phdos, [1e-20, 4, 0, 0])
    np.testing.assert_array_equal(result.source_rows, [2, 1, 0, 4])
    assert result.manifest["source_sha256"] == hashlib.sha256(raw).hexdigest()
    assert path.read_bytes() == raw
    assert result.manifest["operations"][1]["dropped_source_rows_zero_based"] == [3]
    assert result.manifest["renormalization_factor"] == 1
    assert result.manifest["certified_DOS_unit"] is None
    assert not result.axis.flags.writeable


def test_common_support_has_bounded_one_sided_ratio_and_zero_outside(tmp_path):
    result = derive(write(tmp_path, [[0, 0, 0], [1, 6, 2], [2, 9, -1], [3, 4, -2]]))
    q = np.r_[np.linspace(.01, 1.99, 20), 2, 2.5, 3]
    alpha, dos, ratio = result.evaluate(q)
    np.testing.assert_allclose(ratio[:-3], 3, rtol=2e-14)
    np.testing.assert_array_equal(ratio[-3:], [0, 0, 0])
    np.testing.assert_allclose(alpha, dos * ratio)


def test_floor_is_explicit_and_does_not_change_large_values(tmp_path):
    path = write(tmp_path, [[0, 1e-8, 1e-12], [1, 5, 2], [2, 6, 1e-9]])
    minimal = derive(path)
    floored = derive(path, dos_floor=1e-10)
    np.testing.assert_array_equal(minimal.phdos, [1e-12, 2, 1e-9])
    np.testing.assert_array_equal(floored.phdos, [0, 2, 1e-9])
    assert floored.manifest["operations"][2]["floor_original_ordinate_units"] == 1e-10


@pytest.mark.parametrize("rows", [
    [[0, 1, 1]], [[0, 1, 1], [0, 2, 2]], [[-1, 1, 1], [0, 1, 1]],
    [[0, -1, 1], [1, 1, 1]], [[0, 1, 1], [1, np.nan, 1]], [[0, 1], [1, 2]],
])
def test_invalid_source_not_silently_repaired(tmp_path, rows):
    with pytest.raises(ValueError):
        derive(write(tmp_path, rows))


@pytest.mark.parametrize("floor", [-1, np.nan, np.inf])
def test_invalid_floor(tmp_path, floor):
    with pytest.raises(ValueError):
        derive(write(tmp_path, [[0, 0, 0], [1, 1, 1]]), dos_floor=floor)


@pytest.mark.parametrize("query", [-1, 2, np.inf, np.nan])
def test_extrapolation_and_nonfinite_queries_fail(tmp_path, query):
    result = derive(write(tmp_path, [[0, 0, 0], [1, 1, 1]]))
    with pytest.raises(ValueError):
        result.evaluate(query)


def test_finite_extreme_ordinates_cannot_silently_overflow_ratio(tmp_path):
    result = derive(write(tmp_path, [[0, 0, 0], [1, 1e300, 1e-300]]))
    with pytest.raises(ValueError, match="not representable"):
        result.evaluate(1.)
