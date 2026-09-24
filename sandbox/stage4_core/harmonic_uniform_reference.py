"""Exact uniform references and Laplacian modes for weak harmonic diagnostics.

No kinetic closure, deposition, collision rate or photon preparation is chosen.
The spatial basis remains a basis on the 2D graph; truncating it is an explicit
linear-response approximation, never evidence that the detector is 1D.
Finite-frequency branches below use exp(-i omega t), Omega=hbar*omega/(kB Tc).
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
import sys
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import eigsh, spsolve

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pysnspd.experimental.energy_catalog import K_B_J_K, HBAR_J_S, E_CHARGE_C
from pysnspd.experimental.thermal_spatial_usadel import _fixed_indices
from pysnspd.experimental.thesis_circuit_harmonic import partition_inductance

TAU3 = np.diag([1., -1.]).astype(complex)


@dataclass(frozen=True)
class UniformReference:
    temperature_ratio: float
    gap_bar: float
    matsubara: np.ndarray
    gap_equation_residual: float
    phase_stiffness_bar: float


def uniform_reference(temperature_ratio=.9/8.65, matsubara_count=256):
    """Stationary gap of the SAME finite renormalized Matsubara action."""
    t = float(temperature_ratio)
    if not np.isfinite(t) or not 0 < t < 1 or type(matsubara_count) is not int or matsubara_count < 1:
        raise ValueError('0<T/Tc<1 and a positive integer frequency count required')
    epsilon = 2*np.pi*t*(np.arange(matsubara_count)+.5)
    def equation(d):
        return np.log(t)+2*np.pi*t*np.sum(1/epsilon-1/np.hypot(epsilon, d))
    if equation(8.) <= 0:
        raise ValueError('Finite sum has no admitted gap below 8 kB Tc; increase count')
    gap = brentq(equation, 0., 8., xtol=1e-14)
    stiffness = 4*np.pi*t*np.sum(gap*gap/(epsilon*epsilon+gap*gap))
    return UniformReference(t, float(gap), epsilon, float(equation(gap)), float(stiffness))


def thermal_hessian_symbols(eigenvalues, reference):
    """H_radial/phase v = M v times the returned symbols.

    Exact uniform Schur symbols of the finite thermal action. They do not
    define finite-frequency KWT mobility or a kinetic collision time.
    """
    lam = np.asarray(eigenvalues, float)
    if np.any(~np.isfinite(lam)) or np.any(lam < 0):
        raise ValueError('Nonnegative finite graph eigenvalues required')
    epsilon = reference.matsubara
    g = epsilon/np.hypot(epsilon, reference.gap_bar)
    denominator = epsilon+lam[..., None]*g
    radial = 2*np.log(reference.temperature_ratio)+4*np.pi*reference.temperature_ratio*np.sum(
        1/epsilon-g**3/denominator, axis=-1)
    # Use the stationarity identity to avoid subtracting two large terms.
    phase = 4*np.pi*reference.temperature_ratio*np.sum(
        g*g*lam[..., None]/(epsilon*denominator), axis=-1)
    return radial, phase


@dataclass(frozen=True)
class UniformNambuBranch:
    energy_bar: float
    eta_bar: float
    kind: str
    gap_bar: complex
    root: complex
    matrix: np.ndarray
    generator: np.ndarray


def uniform_branch(energy_bar, gap_bar, eta_bar, kind='R'):
    """Causal R or corresponding A, with B=b R/A and positive-real b.

    eta is a numerical contour regulator here, not an admitted collision bath.
    The advanced branch is constructed separately; a response at the same
    harmonic frequency is NOT obtained by taking a dagger of retarded deltaR.
    """
    energy, eta, d = float(energy_bar), float(eta_bar), complex(gap_bar)
    if not all(np.isfinite(x) for x in (energy, eta, d)) or eta <= 0 or kind not in ('R', 'A'):
        raise ValueError('Finite energy/gap, positive eta and branch R or A required')
    z = eta-1j*energy
    root = np.sqrt(z*z+abs(d)**2)
    if root.real < 0:
        root = -root
    D = np.array([[0., d], [np.conj(d), 0.]], complex)
    if kind == 'R':
        generator = D+z*TAU3
    else:
        generator = D-np.conj(z)*TAU3
        root = np.conj(root)
    matrix = generator/root
    return UniformNambuBranch(energy, eta, kind, d, complex(root), matrix, generator)


def shifted_branches(energy_bar, Omega_bar, gap_bar, eta_bar, *, kinds=('R', 'R')):
    """Branches at E +/- Omega/2; Omega=hbar*omega/E0, exp(-i omega t)."""
    omega = float(Omega_bar)
    if not np.isfinite(omega) or len(kinds) != 2:
        raise ValueError('Finite dimensionless angular frequency and two branches required')
    return (uniform_branch(energy_bar+omega/2, gap_bar, eta_bar, kinds[0]),
            uniform_branch(energy_bar-omega/2, gap_bar, eta_bar, kinds[1]))


@dataclass(frozen=True)
class ModalSpectralResponse:
    matrix: np.ndarray
    denominator: complex
    normalization_residual_max: float
    equation_residual_max: float


def modal_response(plus, minus, eigenvalue, delta_generator):
    """Closed form on L v=lambda Mv, retaining all four Nambu entries.

    X = (deltaB - Rplus deltaB Rminus)/(bplus+bminus+2 lambda).
    It solves both mixed normalization and the complete projected matrix
    equation, including retarded/advanced mixed pairs as algebraic references.
    It is not by itself a Keldysh occupation or self-consistency closure.
    """
    lam = float(eigenvalue)
    dB = np.asarray(delta_generator, complex)
    if not np.isfinite(lam) or lam < 0 or dB.shape != (2, 2) or np.any(~np.isfinite(dB)):
        raise ValueError('Nonnegative eigenvalue and finite 2x2 forcing required')
    if plus.gap_bar != minus.gap_bar or plus.eta_bar != minus.eta_bar:
        raise ValueError('Both energies must use the same uniform gap and regulator')
    Rp, Rm, Bp, Bm = plus.matrix, minus.matrix, plus.generator, minus.generator
    denominator = plus.root+minus.root+2*lam
    if abs(denominator) == 0:
        raise ValueError('Exact modal singularity; no artificial denominator floor')
    X = (dB-Rp@dB@Rm)/denominator
    norm = Rp@X+X@Rm
    residual = -.5*(lam*(Rp@X-X@Rm)+Bp@X-X@Bm+dB@Rm-Rp@dB)
    return ModalSpectralResponse(X, complex(denominator), float(np.max(abs(norm))),
                                 float(np.max(abs(residual))))


def graph_laplacian(graph):
    tail, head = graph.edges.T
    c = graph.conductance
    return coo_matrix((np.r_[c, c, -c, -c],
        (np.r_[tail, head, tail, head], np.r_[tail, head, head, tail])),
        shape=(graph.n_nodes, graph.n_nodes)).tocsr()


def end_contacts(graph, axis=0):
    """Identify only minimum/maximum coordinate contacts; sides stay natural."""
    if graph.coordinates_bar is None or axis not in (0, 1):
        raise ValueError('Planar coordinates and axis 0 or 1 required')
    position = graph.coordinates_bar[:, axis]
    lo, hi = float(position.min()), float(position.max())
    if hi <= lo:
        raise ValueError('Distinct end planes required')
    atol = 64*np.finfo(float).eps*max(1., abs(lo), abs(hi), hi-lo)
    left = np.flatnonzero(abs(position-lo) <= atol)
    right = np.flatnonzero(abs(position-hi) <= atol)
    if not len(left) or not len(right) or np.intersect1d(left, right).size:
        raise ValueError('Distinct nonempty terminal node sets required')
    return left, right


@dataclass(frozen=True)
class GraphModes:
    eigenvalues: np.ndarray
    vectors: np.ndarray
    fixed_nodes: np.ndarray
    free_nodes: np.ndarray
    relative_residuals: np.ndarray
    mass_orthogonality_error: float


def dirichlet_modes(graph, count, *, fixed_nodes=None, tolerance=1e-10):
    """M-orthonormal lowest modes of the full graph, zero on fixed contacts.

    Defaults to two end planes even if graph.boundary_nodes also lists sides.
    This distinction is material for a rectangular graph built by the generic
    thermal helper, whose default boundary list contains the whole perimeter.
    """
    if fixed_nodes is None:
        fixed_nodes = np.unique(np.r_[end_contacts(graph)])
    fixed = _fixed_indices(fixed_nodes, graph.n_nodes)
    if not len(fixed):
        raise ValueError('Fixed contacts are required; no zero-mode gauge is inferred')
    free_mask = np.ones(graph.n_nodes, bool); free_mask[fixed] = False
    free = np.flatnonzero(free_mask)
    if type(count) is not int or not 1 <= count <= len(free) or not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError('Positive mode count within free-node dimension and finite tolerance required')
    lap = graph_laplacian(graph)[free][:, free]
    mass = graph.area_weights[free]
    if len(free) <= 64 or count >= len(free)-1:
        inverse_sqrt = 1/np.sqrt(mass)
        dense = lap.toarray()*inverse_sqrt[:, None]*inverse_sqrt[None, :]
        values, vectors = eigh(dense, subset_by_index=(0, count-1))
        vectors = vectors*inverse_sqrt[:, None]
    else:
        values, vectors = eigsh(lap, k=count, M=diags(mass), sigma=0., which='LM',
                                v0=np.linspace(1., 2., len(free)), tol=tolerance)
        order = np.argsort(values); values, vectors = values[order], vectors[:, order]
    if np.any(values <= 0):
        raise ValueError('Dirichlet graph lacks positive modes; check connectivity and contacts')
    full = np.zeros((graph.n_nodes, count))
    for j in range(count):
        vector = vectors[:, j]
        vector /= np.sqrt(np.dot(mass, vector*vector))
        if vector[np.argmax(abs(vector))] < 0:
            vector *= -1
        full[free, j] = vector
    residual = lap@full[free]-mass[:, None]*full[free]*values
    scales = np.linalg.norm(mass[:, None]*full[free]*values, axis=0)
    relative = np.linalg.norm(residual, axis=0)/scales
    orthogonality = full.T@(graph.area_weights[:, None]*full)-np.eye(count)
    return GraphModes(values, full, fixed, free, relative, float(np.max(abs(orthogonality))))


def harmonic_lift(graph, fixed_nodes, fixed_values):
    """Unique scalar graph-harmonic extension of real or complex contacts."""
    fixed = _fixed_indices(fixed_nodes, graph.n_nodes)
    values = np.asarray(fixed_values)
    if not len(fixed) or values.shape != (len(fixed),) or np.any(~np.isfinite(values)):
        raise ValueError('One finite value at each explicit fixed contact required')
    free = np.setdiff1d(np.arange(graph.n_nodes), fixed)
    lap = graph_laplacian(graph)
    lift = np.zeros(graph.n_nodes, dtype=np.result_type(values, float))
    lift[fixed] = values
    if len(free):
        lift[free] = spsolve(lap[free][:, free], -(lap[free][:, fixed]@lift[fixed]))
    if np.any(~np.isfinite(lift)):
        raise ValueError('Contact lift is singular; graph connectivity must be checked')
    return lift


def project_field(graph, modes, field):
    """Mass projection coefficients and omitted norm, after contact lift removal.

    field is a full-node scalar or complex perturbation vanishing at the fixed
    contacts. The error is the mass L2 norm; a large omitted component remains
    visible and is not silently removed from the input field.
    """
    values = np.asarray(field)
    if values.shape != (graph.n_nodes,) or np.any(~np.isfinite(values)):
        raise ValueError('One finite scalar/complex field value per graph node required')
    if np.max(abs(values[modes.fixed_nodes]), initial=0.) > 1e-13*max(1., np.max(abs(values), initial=0.)):
        raise ValueError('Subtract an explicit harmonic contact lift before modal projection')
    coefficient = modes.vectors.T@(graph.area_weights*values)
    reconstruction = modes.vectors@coefficient
    omitted = np.sqrt(np.dot(graph.area_weights, abs(values-reconstruction)**2))
    total = np.sqrt(np.dot(graph.area_weights, abs(values)**2))
    return dict(coefficients=coefficient, reconstructed=reconstruction,
                omitted_mass_norm=float(omitted), original_mass_norm=float(total))


def korzh_scales():
    """Selected K20 FITTED reference, not a new measurement or confidence range."""
    Tc, Tb, diffusion, sheet, thickness, width = 8.65, .9, 5e-5, 608., 7e-9, 80e-9
    sigma = 1/(sheet*thickness)
    n0 = sigma/(2*E_CHARGE_C**2*diffusion)
    energy = K_B_J_K*Tc
    ell = np.sqrt(HBAR_J_S*diffusion/(2*energy))
    time = HBAR_J_S/(2*energy)
    U0 = n0*energy**2*thickness*ell**2
    I0 = 2*E_CHARGE_C/HBAR_J_S*U0
    return dict(Tc_K=Tc, Tb_K=Tb, D_m2_s=diffusion, R_sheet_ohm=sheet,
        thickness_m=thickness, width_m=width, conductivity_S_m=sigma,
        N0_per_J_m3=n0, E0_J=energy, ell0_m=float(ell), tD_s=time,
        graph_energy_unit_J=float(U0), graph_current_unit_A=float(I0),
        potential_v_equals_ephi_over_E0_unit_V=energy/E_CHARGE_C,
        inherited_v_equals_2ephi_over_E0_unit_V=energy/(2*E_CHARGE_C),
        finite_frequency_Omega_definition='hbar*omega/E0 = 2*tD*omega',
        current_definition='I_tail_to_head=(2e/hbar)*U0*(-partial(Fbar)/partial(alpha))',
        N0_convention='single spin; sigma=2*e^2*N0*D',
        scope='K20 fitted reference, no material collision rate selected')


def uniform_inductance_reference(graph, reference, *, total_reference_H=10e-9):
    """Differential zero-current inductance from the uniform thermal graph.

    Both gap and spectral reservoir phases follow the stationary branch. The
    scalar harmonic phase profile has left=0, right=1 and natural side edges.
    Finite-sum stiffness is retained; no London value is silently substituted.
    """
    left, right = end_contacts(graph)
    fixed = np.unique(np.r_[left, right])
    free = np.setdiff1d(np.arange(graph.n_nodes), fixed)
    lap = graph_laplacian(graph)
    phase = np.zeros(graph.n_nodes); phase[right] = 1.
    if len(free):
        phase[free] = spsolve(lap[free][:, free], -(lap[free][:, fixed]@phase[fixed]))
    tail, head = graph.edges.T
    geometric_flux = graph.conductance*(phase[head]-phase[tail])
    divergence = np.zeros(graph.n_nodes)
    np.add.at(divergence, tail, geometric_flux)
    np.add.at(divergence, head, -geometric_flux)
    geometric_conductance = float(np.sum(divergence[left]))
    if geometric_conductance <= 0 or not np.all(np.isfinite(phase)):
        raise ValueError('Positive finite two-terminal graph conductance required')
    scales = korzh_scales()
    current_per_phase = scales['graph_current_unit_A']*reference.phase_stiffness_bar*geometric_conductance
    partition = partition_inductance(total_reference_H, 1/current_per_phase)
    return dict(**asdict(partition), phase_profile=phase,
        geometric_conductance=geometric_conductance,
        current_per_phase_A_per_rad=float(current_per_phase),
        phase_stiffness_bar=reference.phase_stiffness_bar,
        free_current_residual_bar=float(np.max(abs(divergence[free]), initial=0.)),
        left_nodes=left, right_nodes=right,
        scope='Zero-current uniform stationary differential branch, same port planes')
