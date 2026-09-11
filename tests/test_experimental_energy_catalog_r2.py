"""Small R2 regressions for the finite-eta causal batch APIs.

These checks target previously rejected edge roots and the narrower eta=1e-8
regime. They do not replace catalog field-grid, quadrature or transient gates.
"""
import json
import hashlib
from dataclasses import replace
from unittest.mock import patch
import numpy as np
import pytest

from pysnspd.experimental.energy_catalog import (
    OccupationEnergyCatalog, build_occupation_catalog, build_vacuum_catalog,
    energy_at_count_batch, retarded_spectrum_batch, refine_occupation_catalog,
)


ETA = 1e-8


def _assert_causal_solution(energy, c, s, delta, gamma):
    """Check the unsquared physical equation, not the polynomial used to solve it."""
    z = np.asarray(energy) + 1j * ETA
    lhs, rhs = delta * c, (gamma * c - 1j * z) * s
    scale = np.maximum(1, np.maximum(abs(lhs), abs(rhs)))
    assert np.all(np.isfinite(c)) and np.all(np.isfinite(s))
    assert np.all(c.real >= 0)
    assert np.max(abs(lhs-rhs) / scale) < 1e-8
    assert np.max(abs(c*c+s*s-1) / np.maximum(1, abs(c)**2+abs(s)**2)) < 1e-10
    if delta > 0:
        # The retarded auxiliary energy must stay in the upper half plane.
        assert np.all((1j * delta * c / s).imag > 0)


@pytest.mark.parametrize("delta,gamma,energy", [
    # Independently recorded scalar rejection, spectral_crosscheck.json.
    (1.2, 1e-6, 1.1998406167419657),
    # Recorded batch rejection in the compensated-quadrature pilot.
    (1.5, 8.25829e-8, 1.4999690988106482),
    # Nearly repeated quartic roots at a much smaller depairing parameter.
    (1.2, 1e-12, 1.2),
])
def test_recorded_edge_cases_are_causal_batch_solutions(delta, gamma, energy):
    energies = energy + ETA * np.array([-2., -1., 0., 1., 2.])
    c, s = retarded_spectrum_batch(energies, delta=delta, gamma=gamma, eta=ETA)
    _assert_causal_solution(energies, c, s, delta, gamma)
    # Scalar and array dispatch of the same batch API must select one branch.
    cs, ss = retarded_spectrum_batch(energy, delta=delta, gamma=gamma, eta=ETA)
    assert cs.shape == ss.shape == ()
    np.testing.assert_allclose([cs, ss], [c[2], s[2]], rtol=1e-10, atol=1e-12)


@pytest.mark.parametrize("gamma", [1.2*(1-1e-8), 1.2, 1.2*(1+1e-8)])
def test_gap_closing_boundary_keeps_causal_normalized_solution(gamma):
    energy = np.r_[0, np.geomspace(1e-12, 4., 27)]
    c, s = retarded_spectrum_batch(energy, delta=1.2, gamma=gamma, eta=ETA)
    _assert_causal_solution(energy, c, s, 1.2, gamma)
    # At E=0, particle-hole symmetry makes both functions real.
    assert abs(c[0].imag) < 1e-10
    assert abs(s[0].imag) < 1e-10


def test_bcs_edge_and_inverse_analytic_limits_at_small_eta():
    delta = 1.2
    energy = np.array([[0., delta-ETA, delta], [delta+ETA, 2., 4.]])
    z = energy + 1j*ETA
    root = np.sqrt(z*z-delta*delta)
    c, s = retarded_spectrum_batch(energy, delta=delta, gamma=0, eta=ETA)
    np.testing.assert_allclose(c, z/root, rtol=2e-15, atol=0)
    np.testing.assert_allclose(s, 1j*delta/root, rtol=2e-15, atol=0)
    assert c.shape == s.shape == energy.shape
    _assert_causal_solution(energy, c, s, delta, 0)
    count = np.array([[0., ETA/10, ETA], [10*ETA, .2, 3.]])
    reference = count*np.sqrt(1+delta**2/(count**2+ETA**2))
    inverse = energy_at_count_batch(count, delta=delta, gamma=0, eta=ETA)
    np.testing.assert_allclose(inverse, reference, rtol=2e-15, atol=0)
    assert inverse.shape == count.shape


@pytest.mark.parametrize("gamma", [0., 1e-12, 1e-6, 1.2])
def test_normal_batch_limit_is_exact_for_any_depairing(gamma):
    count = np.array([[0., 1e-11, .02], [.3, 1., 12.]])
    energy = energy_at_count_batch(count, delta=0, gamma=gamma, eta=ETA)
    c, s = retarded_spectrum_batch(energy, delta=0, gamma=gamma, eta=ETA)
    np.testing.assert_array_equal(energy, count)
    np.testing.assert_array_equal(c, np.ones_like(count, dtype=complex))
    np.testing.assert_array_equal(s, np.zeros_like(count, dtype=complex))


@pytest.mark.parametrize("delta,gamma", [
    (1.2, 1e-12), (1.2, 1e-6), (1.2, 1.2*(1-1e-8)),
    (1.2, 1.2), (1.2, 1.2*(1+1e-8)), (.08, 1.2),
])
def test_inverse_preserves_order_support_and_count_at_small_eta(delta, gamma):
    count = np.r_[0., np.geomspace(1e-10, 12., 29)]
    energy = energy_at_count_batch(count, delta=delta, gamma=gamma, eta=ETA)
    bcs_upper = count*np.sqrt(1+delta**2/(count**2+ETA**2))
    assert energy[0] == 0
    assert np.all(np.diff(energy) > 0)
    assert np.all(energy >= count)
    assert np.all(energy <= bcs_upper)
    c, s = retarded_spectrum_batch(energy, delta=delta, gamma=gamma, eta=ETA)
    _assert_causal_solution(energy, c, s, delta, gamma)
    # Exact state-count primitive, independently recomputed from returned c,s.
    w = 1j*delta/s
    reconstructed = np.real(w-1j*gamma*delta*delta/(2*w*w))
    assert np.max(abs(reconstructed-count)/np.maximum(1, count)) < 1e-8
    # The API accepts arbitrary input order and must not sort user arrays.
    permutation = np.random.default_rng(401).permutation(len(count))
    shuffled = energy_at_count_batch(count[permutation], delta=delta, gamma=gamma, eta=ETA)
    np.testing.assert_allclose(shuffled, energy[permutation], rtol=1e-11, atol=1e-13)


@pytest.mark.parametrize("api", [retarded_spectrum_batch, energy_at_count_batch])
@pytest.mark.parametrize("values,parameters", [
    ([-1., 0.], {}), ([0., np.nan], {}), ([0., np.inf], {}),
    ([0., 1.], {"delta": -1.}), ([0., 1.], {"delta": np.nan}),
    ([0., 1.], {"gamma": -1.}), ([0., 1.], {"gamma": np.inf}),
    ([0., 1.], {"eta": 0.}), ([0., 1.], {"eta": -1e-8}),
    ([0., 1.], {"eta": np.nan}),
])
def test_invalid_batch_inputs_raise_instead_of_repair(api, values, parameters):
    kwargs = dict(delta=1.2, gamma=1e-6, eta=ETA)
    kwargs.update(parameters)
    with pytest.raises(ValueError):
        api(values, **kwargs)


def _small_vacuum():
    return build_vacuum_catalog(np.array([.5, .65, .8, .95]),
                               np.array([0., 2e-14, 4e-14, 8e-14]),
                               Tc_K=8.65, N0_per_J_m3=1e47, D_m2_s=1.581e-4,
                               analytic=True)


@pytest.mark.parametrize("direct_energy", [False, True])
def test_compensated_transformed_builder_roundtrip_and_nodal_gamma_slope(tmp_path, direct_energy):
    scale = 1e-12
    table = build_occupation_catalog(
        _small_vacuum(), count_nodes=np.array([1e-9, 1e-7, 1e-5, .03, .2, 1., 3.]),
        count_weights=np.full(7, .1), eta=ETA, inverse_method="batch",
        compensated_gamma=True, gamma_coordinate_scale=scale,
        gamma_offset_threshold=1e-10, gamma_offset_order=8, direct_energy=direct_energy,
    )
    path = table.save(tmp_path/"compensated_transformed.npz")
    with np.load(path, allow_pickle=False) as archive:
        metadata = json.loads(str(archive["metadata_json"]))
        version = 'v4' if direct_energy else 'v3'
        assert metadata["schema"] == f"pysnspd.experimental.occupation_catalog.{version}"
        assert metadata["gamma_coordinate_scale"] == scale
        offset_name = 'gamma_energy_offsets' if direct_energy else 'gamma_log_offsets'
        assert offset_name in archive.files
    restored = OccupationEnergyCatalog.load(path)
    assert restored.gamma_coordinate_scale == scale
    restored_offsets = getattr(restored, offset_name)
    np.testing.assert_array_equal(restored_offsets, getattr(table, offset_name))
    np.testing.assert_array_equal(restored_offsets[:, 0, :], 0.)
    assert not restored_offsets.flags.writeable
    if direct_energy:
        np.testing.assert_array_equal(restored.delta_energy_derivatives, table.delta_energy_derivatives)
        assert not restored.delta_energy_derivatives.flags.writeable
    assert not restored.excitation_energies.flags.writeable
    population = .2*np.exp(-((table.count_nodes-.1)/.3)**2)
    np.testing.assert_allclose(restored.evaluate(.72, 3e-14, population),
                               table.evaluate(.72, 3e-14, population), rtol=1e-13, atol=1e-13)
    # A transformed coordinate must retain the derivative with respect to
    # physical Gamma at a field node, rather than accidentally return d/dy.
    delta, gamma = float(restored.vacuum.delta_axis[1]), float(restored.vacuum.gamma_axis[2])
    c, s = retarded_spectrum_batch(restored.excitation_energies[1, 2],
                                   delta=delta, gamma=gamma, eta=ETA)
    np.testing.assert_allclose(restored.energy_kernel(delta, gamma)[2],
                               -s.real*s.imag/c.real, rtol=2e-10, atol=1e-10)


def test_compensated_transform_retains_sub_ulp_change_and_chain_rule():
    # Deliberately synthetic interpolation fixture, NOT a material spectrum:
    # log increment = log(base_i) + .2*delta + beta_i*log1p(Gamma/scale).
    # Its Gamma changes are below/near double-precision ULP, but its derivative
    # is finite. This exact function supplies an independent chain-rule oracle.
    vacuum = _small_vacuum()
    scale, beta = 1e-12, np.array([-1e-15, -2e-15, -3e-15])
    base = np.array([.1, .3, .6])
    y = np.log1p(vacuum.gamma_axis/scale)
    offsets = np.broadcast_to(y[None, :, None]*beta, (4, 4, 3)).copy()
    increments = base*np.exp(.2*vacuum.delta_axis[:, None, None]+offsets)
    energies = np.cumsum(increments, axis=-1)
    derivatives = np.cumsum(increments*beta/(vacuum.gamma_axis[None, :, None]+scale), axis=-1)
    table = OccupationEnergyCatalog(vacuum, np.array([.1, .4, 1.]), np.ones(3),
                                    energies, ETA, derivatives, offsets, scale)
    delta, gamma = .72, 3e-14
    expected_inc = base*np.exp(.2*delta+beta*np.log1p(gamma/scale))
    expected_energy = np.cumsum(expected_inc)
    expected_gamma_derivative = np.cumsum(expected_inc*beta/(gamma+scale))
    actual = table.energy_kernel(delta, gamma)
    np.testing.assert_allclose(actual[0], expected_energy, rtol=3e-15, atol=0)
    np.testing.assert_allclose(actual[1], .2*expected_energy, rtol=3e-13, atol=0)
    np.testing.assert_allclose(actual[2], expected_gamma_derivative, rtol=3e-13, atol=0)


def _ratio_analytic_energy(count, coefficient, delta, gamma, logarithmic):
    """Synthetic E and derivatives written directly in physical (delta,Gamma).

    The added energy is a polynomial in r=Gamma/delta, or in log1p(r/.2),
    anchored to zero at Gamma=0. It is an interpolation oracle, not a DOS model.
    """
    factor = np.sqrt(1+delta*delta/(count*count+ETA*ETA))
    bcs = count*factor
    bcs_delta = count*delta/((count*count+ETA*ETA)*factor)
    if logarithmic:
        coordinate = np.log1p(gamma/(.2*delta))
        correction = coefficient*((delta*delta+delta)*coordinate+2*coordinate**2)
        gamma_force = coefficient*(delta*delta+delta+4*coordinate)/(.2*delta+gamma)
        delta_force = coefficient*((2*delta+1)*coordinate
            -(delta*delta+delta+4*coordinate)*gamma/(delta*(.2*delta+gamma)))
    else:
        correction = coefficient*(delta*gamma+gamma+2*gamma*gamma/(delta*delta))
        gamma_force = coefficient*(delta+1+4*gamma/(delta*delta))
        delta_force = coefficient*(gamma-4*gamma*gamma/(delta**3))
    return bcs+correction, bcs_delta+delta_force, gamma_force, correction


def _synthetic_ratio_catalog(logarithmic=False):
    vacuum = build_vacuum_catalog(np.array([.4, .65, .9, 1.15]),
                                 np.array([0., .4, .8, 1.2]),
                                 Tc_K=8.65, N0_per_J_m3=1e47, D_m2_s=1.581e-4,
                                 analytic=True)
    count = np.array([.2, .7, 1.5, 3.])
    coefficients = np.array([.01, .02, .03, .04])
    ratio = np.array([0., .15, .8, 1.7, 3.2])
    data = [_ratio_analytic_energy(count[None, :], coefficients[None, :], delta,
                                  delta*ratio[:, None], logarithmic)
            for delta in vacuum.delta_axis]
    energies, derivatives_d, derivatives_g, corrections = [
        np.array([row[index] for row in data]) for index in range(4)]
    table = OccupationEnergyCatalog(
        vacuum=vacuum, count_nodes=count, count_weights=np.array([.1, .2, .3, .4]),
        excitation_energies=energies, eta=ETA, gamma_energy_derivatives=derivatives_g,
        gamma_coordinate_scale=.2 if logarithmic else None,
        delta_energy_derivatives=derivatives_d, gamma_energy_offsets=corrections,
        gamma_ratio_axis=ratio,
    )
    return table, coefficients


@pytest.mark.parametrize("logarithmic", [False, True])
def test_ratio_representation_restores_physical_forces_from_analytic_energy(logarithmic):
    table, coefficients = _synthetic_ratio_catalog(logarithmic)
    population = np.array([.15, .04, .2, .07])
    for delta, gamma in ((.53, .09), (.77, .73), (1.07, 1.13)):
        energy, force_d, force_g, _ = _ratio_analytic_energy(
            table.count_nodes, coefficients, delta, gamma, logarithmic)
        np.testing.assert_allclose(table.energy_kernel(delta, gamma),
                                   [energy, force_d, force_g], rtol=2e-12, atol=2e-13)
        expected = np.array(table.vacuum.evaluate(delta, gamma)) + 4*np.array([
            np.dot(kernel*population, table.count_weights)
            for kernel in (energy, force_d, force_g)])
        np.testing.assert_allclose(table.evaluate(delta, gamma, population),
                                   expected, rtol=2e-12, atol=2e-13)


def test_ratio_v5_roundtrip_keeps_internal_axis_separate_from_physical_support(tmp_path):
    table, _ = _synthetic_ratio_catalog(logarithmic=True)
    path = table.save(tmp_path/'ratio_v5.npz')
    with np.load(path, allow_pickle=False) as archive:
        metadata = json.loads(str(archive['metadata_json']))
        assert metadata['schema'] == 'pysnspd.experimental.occupation_catalog.v5'
        assert metadata['array_units']['gamma_ratio_axis'] == 'dimensionless Gamma/absDelta'
        assert len(archive['gamma_ratio_axis']) == 5
        assert len(archive['gamma_axis']) == 4
        assert archive['gamma_axis'][-1] == 1.2
    loaded = OccupationEnergyCatalog.load(path)
    np.testing.assert_array_equal(loaded.gamma_ratio_axis, table.gamma_ratio_axis)
    np.testing.assert_array_equal(loaded.vacuum.gamma_axis, table.vacuum.gamma_axis)
    assert not loaded.gamma_ratio_axis.flags.writeable
    assert loaded.gamma_coordinate_scale == .2
    np.testing.assert_allclose(loaded.energy_kernel(.73, .27),
                               table.energy_kernel(.73, .27), rtol=1e-14, atol=1e-14)


def test_ratio_internal_coverage_does_not_allow_physical_extrapolation():
    table, _ = _synthetic_ratio_catalog()
    # The internal ratio grid represents points beyond physical Gamma=1.2.
    assert table.gamma_ratio_axis[-1]*table.vacuum.delta_axis[-1] > 1.2
    for delta, gamma in ((.4, 1.2), (1.15, 1.2), (.4, 0.)):
        assert np.all(np.isfinite(table.energy_kernel(delta, gamma)))
    for delta, gamma in ((.8, np.nextafter(1.2, np.inf)), (.8, 1.3),
                         (.8, -.01), (.39, .1), (1.16, .1), (0., 1.3)):
        with pytest.raises(ValueError):
            table.energy_kernel(delta, gamma)
    normal = table.energy_kernel(0., 1.2)
    np.testing.assert_array_equal(normal[0], table.count_nodes)
    np.testing.assert_array_equal(normal[1:], np.zeros((2, len(table.count_nodes))))
    with pytest.raises(ValueError, match='full physical Gamma support'):
        replace(table, gamma_ratio_axis=table.gamma_ratio_axis*.5)


@pytest.fixture(scope='module')
def small_physical_ratio_seed():
    vacuum = build_vacuum_catalog(np.array([.4, .65, .9, 1.15]),
                                 np.array([0., .4, .8, 1.2]),
                                 Tc_K=8.65, N0_per_J_m3=1e47, D_m2_s=1.581e-4,
                                 analytic=True)
    return build_occupation_catalog(
        vacuum, count_nodes=np.array([.05, .2, .8, 2.]), count_weights=np.full(4, .25),
        eta=ETA, inverse_method='batch', compensated_gamma=True, direct_energy=True,
        gamma_coordinate_scale=.02, gamma_ratio_axis=np.array([0., .5, 1.2, 3.2]),
    )


def test_incremental_refinement_preserves_old_data_and_solves_only_insertions(
        small_physical_ratio_seed, tmp_path):
    seed = small_physical_ratio_seed
    seed_path = seed.save(tmp_path/'seed.npz')
    seed_hash = hashlib.sha256(seed_path.read_bytes()).hexdigest()
    inserted = np.array([.3, .9])
    ratio = np.sort(np.r_[seed.gamma_ratio_axis, inserted])
    with patch('pysnspd.experimental.energy_catalog.energy_at_count_batch',
               wraps=energy_at_count_batch) as inverse_spy:
        refined = refine_occupation_catalog(seed, ratio, reference_hash=seed_hash)
    # Insertions are above the tiny-Gamma integration threshold, so exactly one
    # inverse is needed for each new (amplitude,ratio) pair, and none for old nodes.
    assert inverse_spy.call_count == len(seed.vacuum.delta_axis)*len(inserted)
    actual_fields = []
    for call in inverse_spy.call_args_list:
        np.testing.assert_array_equal(call.args[0], seed.count_nodes)
        assert call.kwargs['eta'] == seed.eta
        actual_fields.append((call.kwargs['delta'], call.kwargs['gamma']))
    expected_fields = [(delta, delta*r) for delta in seed.vacuum.delta_axis for r in inserted]
    np.testing.assert_array_equal(actual_fields, expected_fields)
    old_positions = np.searchsorted(ratio, seed.gamma_ratio_axis)
    for name in ('excitation_energies', 'delta_energy_derivatives',
                 'gamma_energy_derivatives', 'gamma_energy_offsets'):
        np.testing.assert_array_equal(getattr(refined, name)[:, old_positions, :],
                                      getattr(seed, name))
    np.testing.assert_array_equal(refined.count_nodes, seed.count_nodes)
    np.testing.assert_array_equal(refined.count_weights, seed.count_weights)
    np.testing.assert_array_equal(refined.vacuum.gamma_axis, seed.vacuum.gamma_axis)
    assert refined.vacuum.metadata['refinement_seed_sha256'] == seed_hash
    np.testing.assert_array_equal(refined.vacuum.metadata['refinement_new_ratio_nodes'], inserted)
    loaded = OccupationEnergyCatalog.load(refined.save(tmp_path/'refined.npz'))
    assert loaded.vacuum.metadata['refinement_seed_sha256'] == seed_hash
    assert loaded.vacuum.metadata['refinement_policy'] == refined.vacuum.metadata['refinement_policy']
    np.testing.assert_array_equal(loaded.gamma_ratio_axis, refined.gamma_ratio_axis)
    np.testing.assert_array_equal(loaded.gamma_energy_offsets, refined.gamma_energy_offsets)


def test_incremental_refinement_rejects_removed_nodes_or_changed_limits_before_solving(
        small_physical_ratio_seed):
    seed = small_physical_ratio_seed
    old = seed.gamma_ratio_axis
    invalid_axes = (
        np.sort(np.r_[old[[0, 2, 3]], .7]),  # remove the interior r=.5 node
        np.r_[old, old[-1]+.1],              # extend the upper limit
        np.r_[old[:-1], old[-1]-.1],         # contract the upper limit
        np.r_[.01, old[1:]],                 # move the Gamma=0 anchor
    )
    with patch('pysnspd.experimental.energy_catalog.energy_at_count_batch',
               side_effect=AssertionError('an invalid refinement solved a spectrum')):
        for axis in invalid_axes:
            with pytest.raises(ValueError, match='retain every old node and the same support'):
                refine_occupation_catalog(seed, axis)
