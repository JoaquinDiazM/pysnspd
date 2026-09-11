"""Independent, lightweight v0.4 kinetic diagnostics. No production solver imports.

Run from the repository root. The DAT file is audited as supplied, never repaired.
The collision trajectory uses an explicitly synthetic normal spectrum, not NbN rates.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from numpy.polynomial.legendre import leggauss
from scipy.integrate import cumulative_trapezoid, quad, solve_ivp
from scipy.optimize import brentq
from scipy.special import expit


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "docs/modelo_v0_4"
FIG = OUT / "figuras"
VER = OUT / "verificaciones"
for directory in (FIG, VER):
    directory.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                     "axes.spines.right": False, "savefig.dpi": 190,
                     "figure.constrained_layout.use": True})


def savefig(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"{name}.{ext}")
    plt.close(fig)


def grid(n=256, xmax=30.0):
    z, w = leggauss(n)
    return (z + 1) * xmax / 2, w * xmax / 2


def temperature(E, w, p):
    energy = np.dot(w * E, p)
    if energy <= 0:
        return 0.0
    upper = np.dot(w * E, np.full_like(p, .5))
    if energy >= upper:
        raise ValueError("No positive-temperature FD inverse on the finite grid.")
    return brentq(lambda T: np.dot(w * E, expit(-E / T)) - energy,
                  1e-8, 1e8, xtol=2e-14)


def material_audit(path):
    a = np.loadtxt(path)
    nu, alpha, phdos = a.T
    pos = np.maximum(phdos, 0)
    neg = np.minimum(phdos, 0)
    h_mev_thz = 4.135667696
    signed = np.trapezoid(phdos, nu)
    positive = np.trapezoid(pos, nu)
    negative = np.trapezoid(neg, nu)
    lambda_ep = 2 * np.trapezoid(alpha[1:] / nu[1:], nu[1:])
    cumulative = cumulative_trapezoid(phdos, nu, initial=0)
    fig, ax = plt.subplots(1, 2, figsize=(10.0, 3.8))
    ax[0].plot(nu, phdos, label="Columna original (encabezado: estados/THz)")
    ax[0].fill_between(nu, phdos, 0, where=phdos < 0, color="#d36b56", alpha=.6,
                       label="Peso negativo: se registra, no se corrige")
    ax[0].set(xlabel=r"Frecuencia $\nu$ [THz]", ylabel="DOS tal como está tabulada")
    ax[0].legend(fontsize=8, loc="upper left")
    ax[1].plot(nu, cumulative, label="Integral con unidad declarada")
    ax[1].plot(nu, cumulative * h_mev_thz, "--", label="Hipótesis estados/meV (sin adoptar)")
    ax[1].axhline(3, color="black", ls=":", label="3 modos por átomo, si DOS completa")
    ax[1].set(xlabel=r"Frecuencia máxima [THz]", ylabel="Número integrado de modos")
    ax[1].legend(fontsize=8, loc="lower right")
    fig.suptitle("Auditoría de entrada NbN: la cabecera no determina una DOS admisible")
    savefig(fig, "B_06_auditoria_dos")
    np.savetxt(VER / "B_04_material_audit.csv", np.c_[nu, alpha, phdos, cumulative],
               delimiter=",", header="nu_THz,alpha2F,phdos_as_supplied,cumulative_modes", comments="")
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "rows": len(a), "nu_THz": [float(nu[0]), float(nu[-1])],
            "declared_header": path.read_text().splitlines()[0].strip(),
            "signed_integral_modes": float(signed), "positive_integral_modes": float(positive),
            "negative_integral_modes": float(negative),
            "negative_to_positive_fraction": float(-negative / positive),
            "if_per_meV_signed_integral": float(signed * h_mev_thz),
            "if_per_meV_positive_integral": float(positive * h_mev_thz),
            "lambda_ep": float(lambda_ep), "minimum_phdos": float(phdos.min()),
            "admission": "REJECTED for absolute phonon energy/rates until units, mode count and negative tail are resolved",
            "repairs_applied": []}


def packet_and_bgk():
    # New quadrature and an independent adaptive integral, not the v0.3 Gaussians code.
    # Same occupations as v0.3: u_qp/(4*N0*Delta0**2)=.025, hence u_qp=.1 here.
    energy_target = .1
    centers = [1.2, 6.0]
    widths = [.045, .2]
    records, snapshots = [], []
    for n in (256, 512, 1024):
        x, w = grid(n)
        E = np.sqrt(x*x + 1)
        ps, force = [], []
        for center, width in zip(centers, widths):
            shape = np.exp(-.5*((E-center)/width)**2)
            p = shape * energy_target / (4*np.dot(w*E, shape))
            ps.append(p)
            force.append(4*np.dot(w/E, p))
        T = temperature(E, w, ps[0])
        pf = expit(-E/T)
        force_eq = 4*np.dot(w/E, pf)
        ts = np.linspace(0, 6, 161)
        ft = []
        energy_errors, bounds = [], []
        for p in ps:
            pt = pf[:, None] + (p-pf)[:, None]*np.exp(-ts)[None, :]
            ft.append(4*(w/E)@pt)
            energy_errors.append(float(np.max(abs(4*(w*E)@pt-energy_target))))
            bounds.append([float(pt.min()), float(pt.max())])
        records.append({"nodes": n, "forces": force, "ratio": force[0]/force[1],
                        "T_E": T, "equilibrium_force": force_eq,
                        "max_energy_error": max(energy_errors), "bounds": bounds})
        snapshots = (E, ps, pf, ts, ft, force_eq)
    independent = []
    for center, width in zip(centers, widths):
        def shape(x):
            return np.exp(-.5*((np.sqrt(1+x*x)-center)/width)**2)
        energy_i = quad(lambda x: 4*np.sqrt(1+x*x)*shape(x), 0, 30,
                        points=[np.sqrt(center*center-1)], epsabs=1e-12, epsrel=1e-12)[0]
        force_i = quad(lambda x: 4/np.sqrt(1+x*x)*shape(x), 0, 30,
                       points=[np.sqrt(center*center-1)], epsabs=1e-12, epsrel=1e-12)[0]
        independent.append(energy_target*force_i/energy_i)
    E, ps, pf, ts, ft, feq = snapshots
    fig, ax = plt.subplots(1, 2, figsize=(10.0, 3.8))
    for p, center in zip(ps, centers):
        ax[0].plot(E, p, label=f"Paquete E*={center:g}; misma energía")
    ax[0].plot(E, pf, "k--", label="FD de igual energía")
    ax[0].set(xlim=(1, 7), xlabel=r"$E/\Delta_0$", ylabel="Ocupación")
    ax[0].legend(fontsize=8)
    for tau in (.2, 1., 5.):
        ax[1].plot(ts, np.exp(-ts/tau), label=rf"$\tau_{{\rm kin}}/t_{{\rm ref}}={tau:g}$")
    ax[1].set(xlabel=r"$t/t_{\rm ref}$", ylabel="Diferencia de fuerzas / diferencia inicial")
    ax[1].legend(fontsize=8)
    fig.suptitle("Se retiene f(E): igual energía conserva diferencias de fuerza durante la relajación")
    savefig(fig, "B_07_espectro_y_relajacion")
    np.savetxt(VER / "B_04_bgk_trajectory.csv", np.c_[ts, ft[0], ft[1]],
               delimiter=",", header="t_over_tau,force_low_packet,force_high_packet", comments="")
    return {"quadratures": records, "independent_quad_forces": independent,
            "normalization": "N0=Delta0=1; total u_qp=.1, i.e. u_qp/(4*N0*Delta0^2)=.025 as v0.3",
            "tau_kin_over_t_ref_sweep": [.2, 1., 5.],
            "normalized_force_difference_at_t_ref": [float(np.exp(-1/tau)) for tau in (.2, 1., 5.)],
            "last_vs_independent_max_relative": float(np.max(abs(np.array(records[-1]["forces"])/independent-1))),
            "meaning": "Conservative model relaxation in units of tau_kin; no material time or latency prediction"}


def heat_direction(E, w, p, bath=.08, power=1.0):
    T = max(temperature(E, w, p), bath)
    # Logarithmic rescaling prevents underflow at a large gap/bath ratio.
    log_weight = np.log(E)-np.logaddexp(0, E/T)
    seed = np.exp(log_weight-log_weight.max())
    direction = seed*(1-p)
    return power*direction/(4*np.dot(w*E, direction)), T


def heating_checks():
    cases, refinements = {}, []
    for n in (96, 192, 384):
        x, w = grid(n, 24.)
        E = np.sqrt(1+x*x)
        p0 = np.zeros_like(x)
        pthermal = expit(-E/.3)
        ppacket = .08*np.exp(-.5*((E-4)/.35)**2)
        current = {}
        vacuum_H, _ = heat_direction(E, w, p0)
        for name, p in [("vacuum", p0), ("thermal", pthermal), ("nonthermal", ppacket),
                        ("near_vacuum_FD_seed", 1e-12*pthermal),
                        ("near_vacuum_packet_seed", 1e-12*ppacket)]:
            H, Tstar = heat_direction(E, w, p)
            tangent = E*p*(1-p)
            old_denom = 4*np.dot(w*E, tangent)
            old = tangent/old_denom if old_denom > 0 else None
            current[name] = {"power_error": abs(float(4*np.dot(w*E, H))-1),
                             "minimum_source": float(H.min()), "maximum_source": float(H.max()),
                             "T_star": Tstar, "old_defined": old is not None,
                             "distance_old_l1": None if old is None else float(np.dot(w, abs(H-old))),
                             "new_distance_to_vacuum_l1": float(np.dot(w, abs(H-vacuum_H))),
                             "new_force_per_power": float(4*np.dot(w/E, H)),
                             "old_force_per_power": None if old is None else float(4*np.dot(w/E, old))}
        times = np.linspace(0, 2, 81)
        sol = solve_ivp(lambda t, p: heat_direction(E, w, p, power=.05)[0],
                        (0, 2), p0, method="DOP853", t_eval=times,
                        rtol=2e-9, atol=2e-12)
        if not sol.success:
            raise RuntimeError(sol.message)
        energy = 4*(w*E)@sol.y
        force = 4*(w/E)@sol.y
        refinements.append({"nodes": n, "max_energy_residual": float(max(abs(energy-.05*times))),
                            "minimum_p": float(sol.y.min()), "maximum_p": float(sol.y.max()),
                            "final_force": float(force[-1]), "nfev": sol.nfev})
        cases = current
    fig, ax = plt.subplots(1, 2, figsize=(10.0, 3.8))
    for name, p, label in [("vacuum", p0, "Vacío: fuente 0.4 finita"),
                            ("thermal", pthermal, "FD: coincide con tangente térmica"),
                            ("nonthermal", ppacket, "No térmico: cierre declarado")]:
        H, _ = heat_direction(E, w, p)
        ax[0].plot(x, 4*E*H, label=label)
    ax[0].set(xlim=(0, 5.5), xlabel=r"Coordenada de estados $x/\Delta_0$", ylabel="Peso energético por unidad de x")
    ax[0].legend(fontsize=8)
    ax[1].plot(times, energy, label="Energía integrada desde p=0")
    ax[1].plot(times, .05*times, "k--", label="Energía suministrada P t")
    ax[1].set(xlabel="Tiempo [u. a.]", ylabel="Energía QP [u. a.]")
    ax[1].legend(fontsize=8)
    fig.suptitle("Calor 0.4: energía exacta y fuente definida; la forma espectral sigue siendo un cierre")
    savefig(fig, "B_08_calor_espectral")
    np.savetxt(VER / "B_04_heat_trajectory.csv", np.c_[times, energy, force],
               delimiter=",", header="time,qp_energy,qp_force", comments="")
    # At exactly p=1 the heat source vanishes. Inward directions at p=0 and p=1
    # establish continuous Pauli invariance wherever the global denominator > 0.
    return {"cases_last_grid": cases, "refinement": refinements,
            "final_force_relative_change": abs(refinements[-1]["final_force"]/refinements[-2]["final_force"]-1),
            "old_v03_vacuum_denominator": 0.,
            "entropy_claim": "No universal entropy positivity claimed for inverted occupations",
            "physical_status": "Phenomenological deposition of positive dissipation, not microscopic Joule collision integral"}


def collision_model(n, escape=True):
    # Synthetic normal DOS, midpoint electrons; exact equal-energy reaction grid.
    h = 6/n
    E = (np.arange(n)+.5)*h
    om = np.arange(1, 2*n+1)*h
    we = np.full(n, 4*h)  # N0=1, rho=1
    gp = 9*10*om**2/12**3  # Ni=10, normalized Debye DOS over [0,12]
    wp = gp*h
    si, sj = np.triu_indices(n, k=1)
    sk = sj-si-1
    ri, rj = np.triu_indices(n, k=0)
    rk = ri+rj
    alpha = .03*(om/12)**2
    ks = 8*np.pi*alpha[sk]*h*h
    kr = 8*np.pi*alpha[rk]*h*h*np.where(ri == rj, .5, 1.)
    bath = .3
    peq = expit(-E/bath)
    nb = 1/np.expm1(om/bath)
    cells = np.tile(np.r_[peq, nb], (2, 1))
    shape = np.exp(-.5*((om-6)/1.3)**2)
    cells[0, n:] += .5*shape/np.dot(wp*om, shape)
    initial = np.r_[cells.ravel(), 0.]
    nc = n+2*n

    def rhs(t, y):
        c = y[:-1].reshape(2, nc)
        dc = np.zeros_like(c)
        escaped_power = 0.
        for k in range(2):
            p, phonon = c[k, :n], c[k, n:]
            bs = p[sj]*(1-p[si])*(phonon[sk]+1)-p[si]*(1-p[sj])*phonon[sk]
            br = p[ri]*p[rj]*(phonon[rk]+1)-(1-p[ri])*(1-p[rj])*phonon[rk]
            rs, rr = ks*bs, kr*br
            counts = np.bincount(si, rs, n)-np.bincount(sj, rs, n)
            counts -= np.bincount(ri, rr, n)+np.bincount(rj, rr, n)
            phcounts = np.bincount(sk, rs, 2*n)+np.bincount(rk, rr, 2*n)
            dc[k, :n] = counts/we
            dc[k, n:] = phcounts/wp
            if escape:
                loss = (phonon-nb)/8
                dc[k, n:] -= loss
                escaped_power += np.dot(wp*om, loss)
        # Same face count flux, opposite signs; B13 normal coefficient D=.03.
        face_counts = 4*.03*h*(c[0, :n]-c[1, :n])
        dc[0, :n] -= face_counts/we
        dc[1, :n] += face_counts/we
        return np.r_[dc.ravel(), escaped_power]

    times = np.linspace(0, 20, 161)
    sol = solve_ivp(rhs, (0, 20), initial, t_eval=times, method="DOP853", rtol=2e-9, atol=1e-12)
    if not sol.success:
        raise RuntimeError(sol.message)
    c = sol.y[:-1].reshape(2, nc, -1)
    ue = np.einsum("i,kit->kt", we*E, c[:, :n])
    up = np.einsum("i,kit->kt", wp*om, c[:, n:])
    total = ue.sum(axis=0)+up.sum(axis=0)+sol.y[-1]
    energy0 = float(total[0])
    # Equilibrium RHS uses exact Bose/Fermi identities and a uniform spatial state.
    y_eq = np.r_[np.tile(np.r_[peq, nb], 2), 0.]
    equilibrium = float(np.max(abs(rhs(0, y_eq))))
    result = {"electron_bins": n, "phonon_bins": 2*n, "energy_step": h,
              "max_relative_energy_residual": float(max(abs(total-energy0))/energy0),
              "minimum_f": float(c[:, :n].min()), "maximum_f": float(c[:, :n].max()),
              "minimum_n": float(c[:, n:].min()), "equilibrium_rhs_max": equilibrium,
              "final_qp_excess": float(ue[:, -1].sum()-ue[:, 0].sum()),
              "final_escape": float(sol.y[-1, -1]), "nfev": sol.nfev,
              "escape": escape}
    return result, (times, ue, up, sol.y[-1], total)


def trajectory_checks():
    rs, series = [], None
    for n in (12, 24, 48):
        r, series = collision_model(n)
        rs.append(r)
    closed, _ = collision_model(24, escape=False)
    t, ue, up, esc, total = series
    fig, ax = plt.subplots(1, 2, figsize=(10., 3.8))
    ax[0].plot(t, ue.sum(axis=0)-ue[:, 0].sum(), label="Exceso electrónico")
    ax[0].plot(t, up.sum(axis=0)-up[:, 0].sum()+.5, label="Exceso fonónico")
    ax[0].plot(t, esc, label="Escape acumulado")
    ax[0].axhline(.5, color="black", ls="--", label="Energía depositada inicial")
    ax[0].set(xlabel="Tiempo [u. a.]", ylabel="Energía [u. a.]")
    ax[0].legend(fontsize=8)
    ax[1].plot(t, ue[0]-ue[0, 0], label="Celda de la burbuja")
    ax[1].plot(t, ue[1]-ue[1, 0], label="Celda vecina")
    ax[1].set(xlabel="Tiempo [u. a.]", ylabel="Exceso electrónico por celda")
    ax[1].legend(fontsize=8)
    fig.suptitle("Prueba 0.4: reacciones, difusión y escape en dos celdas; espectro normal sintético")
    savefig(fig, "B_09_trayectoria_conservativa")
    np.savetxt(VER / "B_04_collision_trajectory.csv", np.c_[t, ue.T, up.T, esc, total],
               delimiter=",", header="time,ue_0,ue_1,uph_0,uph_1,escape,total", comments="")
    return {"refinement": rs, "closed_test": closed,
            "last_relative_qp_change": abs(rs[-1]["final_qp_excess"]/rs[-2]["final_qp_excess"]-1),
            "last_relative_escape_change": abs(rs[-1]["final_escape"]/rs[-2]["final_escape"]-1),
            "model": "normal rho=1, N0=1, hbar=1, Ni=10, alpha2F=.03*(Omega/12)^2, Debye cutoff=12, D=.03, tau_escape=8",
            "scope": "Finite conservative two-cell reaction network; not a complete SNSPD transient or convergence proof for moving-gap Usadel transport"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--material", type=Path, required=True)
    args = parser.parse_args()
    report = {"version": "0.4", "host": platform.node(), "python": platform.python_version(),
              "numpy": np.__version__, "scipy": scipy.__version__, "matplotlib": matplotlib.__version__,
              "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "material_audit": material_audit(args.material),
              "spectral_vs_energy_and_bgk": packet_and_bgk(),
              "heat": heating_checks(), "collision_transport_escape": trajectory_checks()}
    # Admission of the material is a separate expected failure, not hidden by
    # passing structural checks on an unrelated synthetic spectrum.
    checks = {
        "independent_packet_quadrature": report["spectral_vs_energy_and_bgk"]["last_vs_independent_max_relative"] < 1e-9,
        "heating_refinement": report["heat"]["final_force_relative_change"] < 1e-6,
        "heating_energy_and_bounds": all(r["max_energy_residual"] < 1e-9 and r["minimum_p"] >= -1e-10 and r["maximum_p"] <= 1+1e-10 for r in report["heat"]["refinement"]),
        "collision_energy_and_bounds": all(r["max_relative_energy_residual"] < 1e-9 and r["minimum_f"] >= -1e-10 and r["maximum_f"] <= 1+1e-10 and r["minimum_n"] >= -1e-10 for r in report["collision_transport_escape"]["refinement"]),
        "collision_refinement": report["collision_transport_escape"]["last_relative_qp_change"] < 1e-3 and report["collision_transport_escape"]["last_relative_escape_change"] < 1e-3,
    }
    report["structural_acceptance"] = checks
    (VER / "B_verificaciones_v04.json").write_text(json.dumps(report, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
    if not all(checks.values()):
        raise AssertionError(checks)
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
