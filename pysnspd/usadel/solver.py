"""Uniform dirty-limit Usadel solver placeholders."""


def solve_usadel_matsubara(delta, q, temperature, material_params):
    """Solve the Matsubara Usadel equation for a given ``delta`` and ``q``."""
    return 0


def solve_usadel_dos(delta, q, energy_grid, material_params):
    """Compute the quasiparticle density of states ``rho(E; |Delta|, q)``."""
    return 0


def compute_supercurrent_from_usadel(delta, q, temperature, material_params):
    """Compute the uniform Usadel supercurrent for the local current-carrying state."""
    return 0
