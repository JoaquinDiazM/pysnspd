"""Independent arithmetic and operator checks for learning fragments E05-E07."""
from math import expm1, isclose, pi
import numpy as np


def main():
    h = 6.62607015e-34
    kb = 1.380649e-23
    classical_energy = 0.5 * (2*pi)**2 * (1e-3)**2
    ratio = classical_energy / h
    thermal_ratio = h*5e9 / (kb*0.020)
    occupancy = 1/expm1(thermal_ratio)
    assert isclose(classical_energy, 1.97e-5, rel_tol=0.003)
    assert isclose(ratio, 2.98e28, rel_tol=0.001)
    assert isclose(occupancy, 6.16e-6, rel_tol=0.003)

    stiffness = np.array([[2., -1.], [-1., 2.]])
    eigenvalues, _ = np.linalg.eigh(stiffness)
    omega = np.sqrt(eigenvalues)
    counts, edges = np.histogram(omega, bins=[0, 1.5, 2])
    heights = counts / np.diff(edges)
    assert np.allclose(eigenvalues, [1, 3])
    assert np.allclose(heights, [2/3, 2])
    assert np.dot(heights, np.diff(edges)) == 2
    # The frequency-to-energy Jacobian preserves both bin counts.
    hbar = h/(2*pi)
    energy_edges = edges*hbar
    energy_heights = heights/hbar
    assert np.allclose(energy_heights*np.diff(energy_edges), counts)

    m, spring, x = 0.2, 8.0, 0.030
    force = -spring*x
    # The finite difference of V independently recovers the analytic force.
    dx = 1e-6
    V = lambda x: 0.5*spring*x*x
    force_fd = -(V(x+dx)-V(x-dx))/(2*dx)
    assert isclose(force_fd, force, rel_tol=1e-10)
    assert isclose(force/m, -1.2)

    c = np.array([[0., 1.], [0., 0.]])
    number = c.T @ c
    assert np.allclose(number, np.diag([0, 1]))
    assert np.allclose(c @ c.T + number, np.eye(2))
    # Jordan-Wigner operators verify the full pair Hamiltonian, including signs.
    identity = np.eye(2)
    c1 = np.kron(c, identity)
    c2 = np.kron(np.diag([1, -1]), c)
    assert np.allclose(c1@c2 + c2@c1, 0)
    xi, delta = 2., 1.+0.3j
    bdg = np.array([[xi, delta], [delta.conjugate(), -xi]])
    psi = [c1, c2.T]
    expanded = sum(psi[i].conj().T @ psi[j] * bdg[i, j]
                   for i in range(2) for j in range(2))
    expected = xi*(c1.T@c1+c2.T@c2-np.eye(4))
    expected = expected.astype(complex)
    expected += delta*c1.T@c2.T + delta.conjugate()*c2@c1
    assert np.allclose(expanded, expected)
    assert np.allclose(expected, expected.conj().T)
    assert np.allclose(np.linalg.eigvalsh([[2, 1], [1, -2]]),
                       [-np.sqrt(5), np.sqrt(5)])
    print(f"E05: E_classical={classical_energy:.9g} J; E/hnu={ratio:.9g}")
    print(f"E05: hnu/kT={thermal_ratio:.9g}; nbar={occupancy:.9g}")
    print(f"E05: frequencies={omega.tolist()}; histogram heights={heights.tolist()}")
    print(f"E06: force={force:.4g} N; acceleration={force/m:.4g} m/s^2")
    print("E06/E07: single-mode algebra, two-mode anticommutation, pair signs, Hermiticity and BdG eigenvalues passed")


if __name__ == '__main__':
    main()
