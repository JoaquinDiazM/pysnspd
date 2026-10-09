"""Pinned FD/BE snapshot inputs for the single stage5.1 Schur system.

Saved numeric/source/reference identities remain exact. Historical producer
code identities are preserved as provenance, not forged into current pins:
updated spectral code must revalidate every unchanged R0 against its current
physical residual, normalization and causal DOS before reuse. No root is
recomputed, photon deposited again, or physical time advanced by this loader.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.constants import Boltzmann as kB, elementary_charge as e, hbar
from pysnspd.experimental import adiabatic_matrix_usadel as ad
from pysnspd.experimental import retarded_spatial_usadel as ret
from pysnspd.experimental.bulk_current_reference import solve_bulk_reference, bulk_current_tangent
from pysnspd.experimental.effective_material import EffectivePhononMaterial
from pysnspd.experimental.electrical_ports import ThesisCircuitParameters
from pysnspd.experimental.joint_evolution import BulkCutLaw, JointSnapshotSystem, JointSpectralSlice
from pysnspd.experimental.joint_nambu_rows import JointNambuRows
from pysnspd.experimental.nambu_moments import NambuMoments
from sandbox.stage4_core.parallel_runtime import linux_resources
from sandbox.stage5_rhs.run_matrix_rhs_audit import load_inputs, energy_rule

ROOT = Path(__file__).resolve().parents[2]
SOURCE_CONTRACT = ('Actual prepared initial common local FD/BE state only: local Markov e-ph '
    'detailed balance zero at transfer; zero collisionless electron source, no eta bath, '
    'no KWT heat, no repeated photon. General conserving Born R/A/K and phonon evolution '
    'remain part of stage5.1 continuation before a general time evolution.')
REFERENCE_CONTRACT = ('Same instantaneous DC gap, phase, links and local electrochemical/Wilson energy coordinate; frozen thermal '
    'K0 is subtracted on the identical real-energy quadrature; independent finite F256 '
    'thermal force and Matsubara current preserve the original admitted DC reference.')


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def _verify_pin(path, expected, ledger):
    path = Path(path)
    actual = sha(path)
    if actual != expected:
        raise ValueError('Pinned numeric/receipt file changed: '+str(path))
    ledger[str(path)] = actual
    return path


def _producer_provenance(mapping, origin):
    result = []
    for filename, saved in mapping.items():
        path = ROOT/filename
        actual = sha(path) if path.is_file() else None
        result.append(dict(origin=origin, path=filename, saved_sha256=saved,
            current_sha256=actual, unchanged=actual == saved,
            policy='provenance only; unchanged numerical R0 independently revalidated'))
    return result


def _validate_cached_R0(graph, gap, alpha, R0, energy_J, eta_J, Eu,
                        *, chart_tolerance=2e-8, matrix_tolerance=1e-7,
                        normalization_tolerance=1e-10, causal_tolerance=1e-8):
    """Revalidate a saved root with the currently loaded physical equations.

    Chart residual divides by m_i*max(1,|d_i|); matrix residual divides by m_i.
    The two tolerances have different normalizations and are reported. These
    are unchanged-root reuse criteria, never a new detector/energy admission.
    """
    R0 = np.asarray(R0, complex)
    if R0.shape != (graph.n_nodes, 2, 2) or np.any(~np.isfinite(R0)):
        raise ValueError('Finite full R0 per current graph node required')
    denominator = 1+R0[:, 0, 0]
    if np.any(abs(denominator) < 1e-14):
        raise ValueError('Saved R0 is outside the admitted stereographic chart')
    a, b = R0[:, 0, 1]/denominator, R0[:, 1, 0]/denominator
    z = float(eta_J)/Eu-1j*float(energy_J)/Eu
    observed = ret.evaluate(graph, gap, z, a, b, alpha)
    f, ft, g = ret.fields(a, b)
    reconstructed = np.empty_like(R0)
    reconstructed[:, 0, 0], reconstructed[:, 1, 1] = g, -g
    reconstructed[:, 0, 1], reconstructed[:, 1, 0] = f, ft
    if not np.allclose(reconstructed, R0, rtol=1e-12, atol=1e-12):
        raise ValueError('Saved full R0 does not match its stereographic representation')
    free = np.ones(graph.n_nodes, bool)
    free[graph.boundary_nodes] = False
    chart = float(np.max(abs(observed.residual[free]/
        (graph.area_weights[free]*np.maximum(1., abs(gap[free])))[:, None]), initial=0.))
    B = ad.generator(gap, energy_J/Eu, eta=eta_J/Eu)
    residual, normalization = ad.spectral_residual(graph, ad.Jet.of(R0), B, alpha=alpha)
    matrix = float(np.max(abs(residual.value[free]/graph.area_weights[free, None, None]), initial=0.))
    normal = float(np.max(abs(normalization.value), initial=0.))
    minimum_DOS = float(np.min(R0[:, 0, 0].real))
    if chart > chart_tolerance or matrix > matrix_tolerance or normal > normalization_tolerance:
        raise ValueError('Unchanged cached R0 fails the current spectral equation/normalization')
    if minimum_DOS < -causal_tolerance:
        raise ValueError('Unchanged cached R0 has negative causal DOS; no clipping')
    return dict(chart_scaled_max=chart, matrix_per_mass_max=matrix,
        normalization_max=normal, minimum_DOS=minimum_DOS,
        chart_tolerance=chart_tolerance, matrix_tolerance=matrix_tolerance,
        normalization_tolerance=normalization_tolerance, causal_tolerance=causal_tolerance)


def load_joint_initial(initial_path):
    """Return all pinned source arrays plus original DC inputs and R0 fields."""
    path = Path(initial_path)
    receipt = json.loads((path/'receipt.json').read_text())
    if (receipt.get('schema') != 'pysnspd.stage5.thermalized_photon_initialization.v1'
            or not receipt.get('source_energy_already_in_populations')
            or not receipt.get('source_must_not_be_applied_again')
            or not receipt.get('tail_screen_passed')):
        raise ValueError('An admitted-tail, single-absorption FD/BE initializer is required')
    ledger = {str(path/'receipt.json'): sha(path/'receipt.json')}
    _verify_pin(path/'initial_state.npz', receipt['fields_sha256'], ledger)
    with np.load(path/'initial_state.npz', allow_pickle=False) as archive:
        data = {key: archive[key].copy() for key in archive.files}
    if any(np.any(~np.isfinite(value)) for value in data.values()):
        raise ValueError('Nonfinite prepared source fields')
    if (not bool(data['source_energy_already_in_populations']) or np.any(data['h_T'] != 0)
            or np.any(data['phi_neutral_V'] != 0)):
        raise ValueError('This factory requires the actual neutral common-temperature FD/BE transfer')
    cache = Path(receipt['spectral_cache'])
    _verify_pin(cache/'identity.json', receipt['cached_identity_sha256'], ledger)
    identity = json.loads((cache/'identity.json').read_text())
    plan = identity['plan']
    reference = Path(plan['reference'])
    _verify_pin(reference, receipt['reference_sha256'], ledger)
    if identity['reference_sha256'] != receipt['reference_sha256']:
        raise ValueError('Initial and cached DC reference identities differ')
    observer = Path(plan['thermal_observer'])
    # These metadata files determine admitted electrical/reference inputs.
    # Record their current identity as well as the historical numeric pins:
    # any later change must invalidate a stage5.1 checkpoint resume.
    for metadata in (reference.parent/'summary.json',observer/'manifest.json'):
        ledger[str(metadata)]=sha(metadata)
    for filename, key in [('summary.json', 'observer_summary_sha256'),
            ('fields.npz', 'observer_fields_sha256'),
            ('spectral_receipts/mode_00000.npz', 'observer_lowest_mode_sha256')]:
        _verify_pin(observer/filename, identity[key], ledger)
    w = load_inputs(plan)
    if w['observed']['count'] != 256 or w['summary']['plan']['matsubara_count'] != 256:
        raise ValueError('The current photon campaign requires its original finite F256 reference')
    if not np.allclose(data['delta_DC_J'], w['gap']*w['Eu'], rtol=1e-12, atol=0):
        raise ValueError('Prepared photon and pinned spectral/DC gaps differ')
    source = Path(receipt['source_path'])
    _verify_pin(source/'receipt.json', receipt['source_sha256'], ledger)
    source_receipt = json.loads((source/'receipt.json').read_text())
    _verify_pin(source/'source_and_dc_fields.npz', source_receipt['fields_sha256'], ledger)
    if (source_receipt['reference_sha256'] != receipt['reference_sha256']
            or source_receipt['material_sha256'] != receipt['material_sha256']):
        raise ValueError('Prepared source material/reference mismatch')
    material_path = _verify_pin(receipt['material_path'], receipt['material_sha256'], ledger)
    material_object = EffectivePhononMaterial.load(material_path)
    material = material_object.manifest
    for key, expected in [('Tc_K', w['summary']['plan']['Tc_K']),
            ('diffusivity_m2_s', w['summary']['plan']['diffusion_m2_s']),
            ('sheet_resistance_ohm', w['summary']['plan']['sheet_resistance_ohm'])]:
        if not np.isclose(material[key], expected, rtol=1e-12):
            raise ValueError('Effective material differs from the admitted DC device')
    if not np.allclose(data['volume_m3'], w['area']*material['film_thickness_m'], rtol=1e-12, atol=0):
        raise ValueError('Source and DC finite-volume ownership differ')
    energy, weights = energy_rule(plan, w['Eu'])
    energy *= w['Eu']
    producer = _producer_provenance(identity['sources'], 'base spectral cache')
    extension = receipt.get('static_spectral_extension')
    if extension is not None:
        ext = Path(extension['path'])
        _verify_pin(ext/'identity.json', extension['identity_sha256'], ledger)
        _verify_pin(ext/'summary.json', extension['summary_sha256'], ledger)
        extended = json.loads((ext/'identity.json').read_text())
        completed = json.loads((ext/'summary.json').read_text())
        if (completed['status'] != 'STATIC_CAUSAL_R0_EXTENSION_COMPLETED'
                or extended['base_identity_sha256'] != receipt['cached_identity_sha256']
                or extended['base_reference_sha256'] != receipt['reference_sha256']
                or extended['base_observer_fields_sha256'] != identity['observer_fields_sha256']):
            raise ValueError('Static UV extension is not the unchanged original DC cache')
        extE, extW = energy_rule(extended['plan'], w['Eu'])
        extE *= w['Eu']
        if (not np.array_equal(extE, extended['added_energy_J'])
                or not np.array_equal(extW, extended['added_weights_J'])
                or extE[0] <= energy[-1]
                or extended['plan']['energy_panel_edges_J'][0] != plan['energy_panel_edges_J'][-1]):
            raise ValueError('Static extension changed or overlaps the original quadrature')
        energy, weights = np.r_[energy, extE], np.r_[weights, extW]
        producer.extend(_producer_provenance(extended['sources'], 'static UV extension'))
    if (not np.array_equal(data['electron_energy_J'], energy)
            or not np.array_equal(data['electron_weights_J'], weights)):
        raise ValueError('Prepared physical-energy quadrature differs from its pinned construction')
    expected_shape = (len(energy), w['graph'].n_nodes)
    if data['h_L'].shape != expected_shape or data['h_T'].shape != expected_shape:
        raise ValueError('Complete hL/hT fields per energy/node required')
    T = data['common_temperature_K']
    if T.shape != (w['graph'].n_nodes,) or np.any(T <= 0):
        raise ValueError('Positive prepared node temperatures are required for initial FD derivatives')
    if not np.allclose(data['h_L'], np.tanh(energy[:, None]/(2*kB*T)), rtol=1e-13, atol=1e-15):
        raise ValueError('Initial fields do not implement the declared common-temperature FD preparation')
    if np.any(abs(data['h_L']+data['h_T']) > 1+32*np.finfo(float).eps) or np.any(abs(data['h_L']-data['h_T']) > 1+32*np.finfo(float).eps):
        raise ValueError('Initial Pauli support failed; no clipping')
    entries = list(receipt['cached_root_sha256'].items())
    if len(entries) != len(energy):
        raise ValueError('A complete energy-to-root identity map is required')
    roots, validation, paths = [], [], []
    for index, (filename, digest) in enumerate(entries):
        root_path = _verify_pin(filename, digest, ledger)
        with np.load(root_path, allow_pickle=False) as archive:
            R0 = archive['R0'].copy()
            if 'energy_J' in archive and float(archive['energy_J']) != float(energy[index]):
                raise ValueError('Saved static-root energy mapping changed')
            if 'eta_J' in archive and float(archive['eta_J']) != receipt['eta_J']:
                raise ValueError('Saved static-root eta changed')
        if not np.array_equal(R0[:, 0, 0].real, data['rho_energy_node'][index]):
            raise ValueError('Initial source-to-root DOS mapping changed')
        metrics = _validate_cached_R0(w['graph'], w['gap'], w['alpha'], R0,
            energy[index], receipt['eta_J'], w['Eu'])
        validation.append(dict(energy_index=index, energy_J=float(energy[index]), **metrics))
        roots.append(R0); paths.append(str(root_path))
    with np.load(reference, allow_pickle=False) as archive:
        current_DC = archive['current_bar']*w['Eu']/(2*e*w['summary']['plan']['sheet_resistance_ohm'])
    with np.load(observer/'fields.npz', allow_pickle=False) as archive:
        thermal_current = archive['current_A'].copy()
    current_roundoff = 64*np.finfo(float).eps*max(float(np.max(abs(thermal_current))), float(np.max(abs(current_DC))))
    if not np.allclose(thermal_current, current_DC, rtol=1e-12, atol=current_roundoff):
        raise ValueError('Independent thermal Matsubara current differs from the original DC reference')
    producer.extend(_producer_provenance(receipt.get('source_code_sha256', {}), 'FD/BE initializer'))
    producer.extend(_producer_provenance(source_receipt.get('source_sha256', {}), 'Gaussian source'))
    provenance = dict(pinned_files_sha256=ledger, producer_history=producer,
        changed_producer_count=sum(not item['unchanged'] for item in producer),
        current_R0_validation=validation, spectral_roots_repeated=0,
        source_reapplied=False, physical_time_steps=0,
        producer_changes_are_not_numeric_identity_changes=True,
        reference_contract=REFERENCE_CONTRACT, source_contract=SOURCE_CONTRACT,
        thermal_force_modes=256, full_energy_count=len(energy),
        current_Matsubara_DC_max_difference_A=float(np.max(abs(thermal_current-current_DC))))
    data.update(receipt=receipt, identity=identity, inputs=w, roots_R0=np.asarray(roots),
        root_paths=paths, material=material, effective_material=material_object,
        thermal_edge_current_A=thermal_current, provenance=provenance)
    return data


def _resident_bytes():
    for line in Path('/proc/self/status').read_text().splitlines():
        if line.startswith('VmRSS:'):
            return int(line.split()[1])*1024
    return 0


def build_initial_system(initial_path, energy_indices=None, resource_memory_limit_bytes=None):
    """Construct the full cached 112-energy system; subsets are tests only.

    Cache construction is serial. The caller controls CPU affinity and the
    outer process budget. Available-memory admission is refreshed here; a
    caller limit can reduce but never raise the 90% ceiling. No affine matrix
    or long solve is launched by this function.
    """
    resources = linux_resources()
    ceiling = int(.9*resources['available_memory_bytes'])
    if resource_memory_limit_bytes is not None:
        if isinstance(resource_memory_limit_bytes, bool) or not isinstance(resource_memory_limit_bytes, (int, np.integer)) or resource_memory_limit_bytes <= 0:
            raise ValueError('Positive integer memory ceiling required')
        ceiling = min(ceiling, int(resource_memory_limit_bytes))
    data = load_joint_initial(initial_path)
    w = data['inputs']; graph = w['graph']; Eu = w['Eu']; p = w['summary']['plan']
    total = len(data['electron_energy_J'])
    if energy_indices is None:
        selected = np.arange(total)
    else:
        selected = np.asarray(energy_indices)
        if (selected.ndim != 1 or selected.dtype.kind not in 'iu' or not len(selected)
                or np.any(selected < 0) or np.any(selected >= total)
                or np.any(np.diff(selected) <= 0)):
            raise ValueError('Test subset must contain ordered unique valid energy indices')
    q = w['q']
    bulk_reference = solve_bulk_reference(p['T_K']/p['Tc_K'], p['matsubara_count'], q)
    if not np.allclose(abs(w['gap'][graph.boundary_nodes]), bulk_reference.gap_bar, rtol=1e-10, atol=1e-12):
        raise ValueError('Saved cuts do not match the original homogeneous current-carrying bulk')
    tangent = bulk_current_tangent(bulk_reference, width_m=p['width_m'], Tc_K=p['Tc_K'],
        diffusion_m2_s=p['diffusion_m2_s'], sheet_resistance_ohm=p['sheet_resistance_ohm'])
    law = BulkCutLaw(bulk_reference.gap_bar, q, tangent.gap_bar_derivative_per_A,
        tangent.q_bar_derivative_per_A,
        'Analytic tangent of the SAME original finite-N256 bulk DC branch; no new inductance or fitted rate')
    x = graph.coordinates_bar[:, 0]
    left = np.flatnonzero(np.isclose(x, x.min(), rtol=0, atol=1e-10))
    right = np.flatnonzero(np.isclose(x, x.max(), rtol=0, atol=1e-10))
    if not np.array_equal(np.sort(np.r_[left, right]), np.sort(graph.boundary_nodes)):
        raise ValueError('Only the two complete longitudinal superconducting cuts may be prescribed')
    obs = NambuMoments(data['electron_energy_J'][selected], data['electron_weights_J'][selected],
        data['volume_m3'], data['material']['N0_per_J_m3'], Eu, Eu/(2*e*p['sheet_resistance_ohm']))
    rows = JointNambuRows(obs, graph.edges, left_nodes=left, right_nodes=right)
    cut_current = w['summary']['metrics']['contact_reactions']
    IDC = .5*obs.I0*(cut_current['left']['outgoing_link_current_bar']
        -cut_current['right']['outgoing_link_current_bar'])
    circuit = ThesisCircuitParameters(Lk_ext_H=w['summary']['fixed_external_inductance_H'],
        V_bias_V=1e4*IDC)
    slices, Kth, Fth, cache_bytes = [], [], [], 0
    last_reserve = 64*1024**2
    for index in selected:
        if _resident_bytes()+last_reserve > ceiling:
            raise MemoryError('Current available/caller memory budget cannot admit another spectral cache')
        h = np.column_stack((data['h_L'][index], data['h_T'][index]))
        hE = np.column_stack(((1-h[:, 0]**2)*Eu/(2*kB*data['common_temperature_K']),
            np.zeros(graph.n_nodes)))
        R0 = data['roots_R0'][index]
        sl = JointSpectralSlice(graph, w['gap'], w['alpha'], R0,
            data['electron_energy_J'][index], data['receipt']['eta_J'], Eu,
            h, hE, graph.boundary_nodes, law, projected_source=np.zeros_like(h),
            source_contract=SOURCE_CONTRACT)
        slices.append(sl)
        owned = sl.spectral.estimated_storage_bytes+sl.storage.data.nbytes+sl.storage.indices.nbytes+sl.storage.indptr.nbytes
        owned += sl.tangent.jacobian.data.nbytes+sl.tangent.jacobian.indices.nbytes+sl.tangent.jacobian.indptr.nbytes
        for factor in (sl.storage_factor, sl.tangent.factor):
            if factor is not None:
                for matrix in (factor.L, factor.U):
                    owned += matrix.data.nbytes+matrix.indices.nbytes+matrix.indptr.nbytes
        cache_bytes += int(owned)
        last_reserve = max(64*1024**2, 2*int(owned))
        if _resident_bytes() > ceiling:
            raise MemoryError('Spectral cache exceeds the admitted memory ceiling')
        h0 = np.tanh(data['electron_energy_J'][index]/(2*kB*p['T_K']))
        H0 = ad.distribution(np.full(graph.n_nodes, h0), np.zeros(graph.n_nodes))
        R = ad.Jet.of(R0)
        _, flux, K = ad.kinetic_residual(graph, R, H0,
            ad.generator(w['gap'], data['electron_energy_J'][index]/Eu), alpha=w['alpha'])
        Kth.append(K.value); Fth.append(flux.value)
    system = JointSnapshotSystem(graph, w['gap'], w['alpha'], slices, rows,
        K_thermal_zero=np.asarray(Kth), thermal_gradient_J_per_gapbar=w['Fgradient'],
        thermal_edge_current_A=data['thermal_edge_current_A'], thermal_edge_flux_zero=np.asarray(Fth),
        circuit=circuit, circuit_state=np.array([IDC, IDC, 0.]), bulk_law=law,
        reference_contract=REFERENCE_CONTRACT)
    system.initial_data = data
    system.provenance = dict(data['provenance'], selected_energy_indices=selected.tolist(),
        full_energy_quadrature=bool(np.array_equal(selected, np.arange(total))),
        subset_test_only=not np.array_equal(selected, np.arange(total)),
        full_snapshot_admitted=False, energy_admitted=False, temporal_admitted=False,
        cached_operator_array_bytes=cache_bytes, process_resident_bytes=_resident_bytes(),
        memory_limit_bytes=ceiling, original_fixed_external_inductance_H=circuit.Lk_ext_H,
        circuit_initial_state_SI=system.y.tolist(), bulk_continuation_q_bar=q,
        bulk_cut_law=law.__dict__, source_reapplied=False, physical_time_steps=0)
    return system
