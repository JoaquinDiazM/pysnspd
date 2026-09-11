"""Independent uniform Usadel audit for A v0.4; no production imports.

Solve the positive renormalized Matsubara energy with bracketed Newton,
instead of the angle iteration and analytic tail used in v0.3. Raw sums
are extrapolated in cutoff independently. All quantities below use k_B Tc
as energy unit and Q=q sqrt(hbar D/(2 k_B Tc)).
"""
from pathlib import Path
from functools import lru_cache
import json
import platform
import time
import numpy as np
from scipy.optimize import brentq, minimize_scalar
import scipy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'docs/modelo_v0_4'
FIG, DATA = OUT / 'figuras', OUT / 'verificaciones'
for directory in (FIG, DATA):
    directory.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({'font.size': 11, 'axes.spines.top': False,
                     'axes.spines.right': False, 'savefig.dpi': 210,
                     'font.family': 'DejaVu Sans'})
COLORS = ['#146a8a', '#bf6334', '#59894b', '#80639c']


@lru_cache(maxsize=100)
def matsubara(T, n):
    return np.pi*T*(2*np.arange(n)+1)


def spectral_matsubara(T, delta, gamma, n):
    eps = matsubara(T, n)
    lo, hi = eps.copy(), eps+gamma
    u = hi.copy()
    for iterations in range(60):
        r = np.hypot(u, delta)
        c = u/r
        f = u-eps-gamma*c
        lo = np.where(f < 0, u, lo)
        hi = np.where(f >= 0, u, hi)
        df = 1-gamma*delta**2/r**3
        candidate = u-f/df
        candidate = np.where((candidate < lo) | (candidate > hi),
                             (lo+hi)/2, candidate)
        step = np.max(np.abs(candidate-u)/np.maximum(1, np.abs(u)))
        u = candidate
        if step < 3e-15:
            break
    r = np.hypot(u, delta)
    s, c = delta/r, u/r
    residual = np.max(np.abs(delta*c-(eps+gamma*c)*s))
    if residual > 2e-11*max(1, delta):
        raise AssertionError(('spectral branch', T, delta, gamma, residual))
    return eps, s, c, u, r, float(residual)


def state(T, delta, Q, n=2048):
    gamma = Q**2
    eps, s, c, u, r, residual = spectral_matsubara(T, delta, gamma, n)
    # Stable identities: r-eps=delta^2/(r+u)+gamma*c, avoiding large subtraction.
    r_minus_eps = delta**2/(r+u)+gamma*c
    terms = s*s*(r_minus_eps**2/eps+eps*s*s/(1+c)**2+gamma)
    f = delta**2*np.log(T)+2*np.pi*T*np.sum(terms)
    force = 2*(delta*np.log(T)+2*np.pi*T*np.sum(s*r_minus_eps/eps))
    conjugate_Q = 4*np.pi*T*Q*np.sum(s*s)
    return np.array([f, force, conjugate_Q]), residual


def derivative(fun, x, h):
    return (fun(x-2*h)-8*fun(x-h)+8*fun(x+h)-fun(x+2*h))/(12*h)


def extrapolate(ns, values):
    # Midpoint Matsubara sums admit an inverse-cutoff asymptotic expansion.
    # The intercept is compared with a disjointly shifted cutoff window.
    ns, values = np.asarray(ns, float), np.asarray(values, float)
    return np.polynomial.polynomial.polyfit(ns[0]/ns, values, len(ns)-1)[0]


def converged(T, delta, Q):
    ns = np.array([256, 512, 1024, 2048, 4096])
    vals = np.array([state(T, delta, Q, int(n))[0] for n in ns])
    return extrapolate(ns, vals)


def eq_gap(T, Q):
    if converged(T, 1e-6, Q)[1] > 0:
        return 0.
    return brentq(lambda d: converged(T, d, Q)[1], 1e-6, 2.4,
                  xtol=3e-12)


def retarded(E, delta, gamma, eta):
    z = E+1j*eta
    if gamma == 0:
        u = z
    elif E > 5*delta:
        u = z+1j*gamma
        for _ in range(30):
            w = np.sqrt(u*u-delta*delta)
            if w.imag < 0:
                w = -w
            step = (u-z-1j*gamma*u/w)/(1+1j*gamma*delta**2/w**3)
            u -= step
            if abs(step) < 2e-14*max(1, abs(u)):
                break
    else:
        # Squaring produces spurious roots; screen them using the original equation.
        roots = np.roots([1, -2*z, z*z-delta*delta+gamma*gamma,
                         2*z*delta*delta, -z*z*delta*delta])
        candidates = []
        for candidate in roots:
            w = np.sqrt(candidate*candidate-delta*delta)
            if w.imag < 0:
                w = -w
            if candidate.imag > 0 and candidate.real > -1e-10:
                score = abs(candidate-z-1j*gamma*candidate/w)
                candidates.append((score, candidate))
        if not candidates:
            raise AssertionError(('no retarded root', E, gamma, eta))
        u = min(candidates, key=lambda row: row[0])[1]
    w = np.sqrt(u*u-delta*delta)
    if w.imag < 0:
        w = -w
    c, s = u/w, 1j*delta/w
    residual = abs(delta*c-(gamma*c-1j*z)*s)/max(1, abs(delta*c), abs((gamma*c-1j*z)*s))
    return c, s, residual


def save(fig, name):
    for extension in ('png', 'pdf'):
        fig.savefig(FIG/f'{name}.{extension}', bbox_inches='tight')
    plt.close(fig)


def main():
    started = time.perf_counter()
    points = [(T, d, Q) for T in (.05, .15, .4, .8, .98)
              for d in (.25, .8, 1.5, 2.) for Q in (0., .15, .45, .75)]
    steps = np.array([1e-2, 3e-3, 1e-3, 3e-4, 1e-4])
    derivatives = []
    spectral_max = 0.
    for T, delta, Q in points:
        exact, residual = state(T, delta, Q)
        spectral_max = max(spectral_max, residual)
        for h in steps:
            force = derivative(lambda d: state(T, d, Q)[0][0], delta, h)
            current = derivative(lambda q: state(T, delta, q)[0][0], Q, h)
            mixed_a = derivative(lambda q: state(T, delta, q)[0][1], Q, h)
            mixed_b = derivative(lambda d: state(T, d, Q)[0][2], delta, h)
            derivatives.append([T, delta, Q, h,
                                abs(force-exact[1])/max(1, abs(exact[1])),
                                abs(current-exact[2])/max(1, abs(exact[2])),
                                abs(mixed_a-mixed_b)/max(1, abs(mixed_a), abs(mixed_b))])
    derivatives = np.array(derivatives)
    np.savetxt(DATA/'A_derivadas_malla_v04.csv', derivatives, delimiter=',',
               header='T_over_Tc,delta_kBTc,Q,h,force_scaled_error,current_scaled_error,mixed_scaled_error', comments='')
    max_by_h = np.array([np.max(derivatives[derivatives[:, 3] == h, 4:], axis=0) for h in steps])
    chosen = max_by_h[2]
    assert np.max(chosen) < 2e-7, ('derivative errors at h=0.001', chosen)

    ns = np.array([128, 256, 512, 1024, 2048, 4096, 8192])
    convergence, raw_errors = [], []
    for T, delta, Q in points:
        values = np.array([state(T, delta, Q, int(n))[0] for n in ns])
        reference = extrapolate(ns[-5:], values[-5:])
        earlier = extrapolate(ns[-6:-1], values[-6:-1])
        diff = np.abs(reference-earlier)/np.maximum(1, np.abs(reference))
        convergence.append([T, delta, Q, *diff])
        raw_errors.append(np.abs(values-reference)/np.maximum(1, np.abs(reference)))
    convergence, raw_errors = np.array(convergence), np.array(raw_errors)
    assert np.max(convergence[:, 3:]) < 2e-8, ('cutoff-window dependence', np.max(convergence[:, 3:]))
    np.savetxt(DATA/'A_convergencia_malla_v04.csv', convergence, delimiter=',',
               header='T_over_Tc,delta_kBTc,Q,F_extrapolation_change,X_extrapolation_change,Pi_extrapolation_change', comments='')
    np.savetxt(DATA/'A_convergencia_cortes_v04.csv', np.column_stack([ns, np.max(raw_errors, axis=0)]), delimiter=',',
               header='n_matsubara,max_F_scaled_error,max_X_scaled_error,max_Pi_scaled_error', comments='')
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.0), layout='constrained')
    for i, label in enumerate(('Fuerza de amplitud', 'Corriente conjugada', 'Derivadas mixtas')):
        axes[0].loglog(steps, np.maximum(max_by_h[:, i], 1e-16), 'o-', color=COLORS[i], label=label)
    axes[0].set(xlabel='Paso de diferencia finita (adimensional)', ylabel='Error máximo escalado',
                title='80 estados: no basta un solo paso')
    axes[0].legend(fontsize=8)
    for i, label in enumerate(('Energía', 'Fuerza', 'Corriente conjugada')):
        axes[1].loglog(ns, np.maximum(np.max(raw_errors, axis=0)[:, i], 1e-16), 'o-', color=COLORS[i], label=label)
    axes[1].set(xlabel='Número de frecuencias de Matsubara', ylabel='Error de la suma sin cola',
                title='El corte crudo exige convergencia')
    axes[1].legend(fontsize=8)
    save(fig, 'A_04_auditoria_independiente')

    spectral_rows, edge_rows = [], []
    energies = np.unique(np.r_[np.linspace(0, 3, 601), 1.])
    etas = [1e-2, 1e-3, 1e-4, 1e-5]
    ratios = [0., .1, .4, .8, 1.1]
    max_retarded, min_dos, max_normalization = 0., 1., 0.
    for ratio in ratios:
        for eta in etas:
            for E in energies:
                c, s, res = retarded(E, 1., ratio, eta)
                max_retarded = max(max_retarded, res)
                min_dos = min(min_dos, c.real)
                max_normalization = max(max_normalization, abs(c*c+s*s-1))
                spectral_rows.append([ratio, eta, E, c.real, s.real, s.imag, res])
        if 0 < ratio < 1:
            result = minimize_scalar(lambda v: -v*(1-ratio/np.sqrt(1-v*v)),
                                     bounds=(0, np.sqrt(1-ratio**2)), method='bounded',
                                     options={'xatol': 1e-13})
            predicted = (1-ratio**(2/3))**1.5
            leaks = [retarded(.5*predicted, 1., ratio, eta)[0].real for eta in etas]
            edge_rows.append({'Gamma_over_delta': ratio, 'edge_from_stationary_equation': float(-result.fun),
                              'edge_formula': predicted, 'edge_absolute_error': abs(-result.fun-predicted),
                              'subgap_dos_at_half_edge': leaks})
    spectral_rows = np.array(spectral_rows)
    assert max_retarded < 1e-8 and min_dos > -1e-10 and max_normalization < 1e-8
    assert max(e['edge_absolute_error'] for e in edge_rows) < 1e-10
    assert all(e['subgap_dos_at_half_edge'][-1] < .002*e['subgap_dos_at_half_edge'][0] for e in edge_rows)
    np.savetxt(DATA/'A_espectros_causales_v04.csv', spectral_rows, delimiter=',',
               header='Gamma_over_delta,eta_over_delta,E_over_delta,rho,N2,R2,equation_scaled_residual', comments='')
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.0), layout='constrained')
    for i, ratio in enumerate((0., .1, .4, .8)):
        use = spectral_rows[(spectral_rows[:, 0] == ratio) & (spectral_rows[:, 1] == 1e-4)]
        axes[0].plot(use[:, 2], use[:, 3], color=COLORS[i], label=f'Γ / |Δ| = {ratio:g}')
    axes[0].set(xlabel=r'$E/|\Delta|$', ylabel='DOS normalizada', ylim=(0, 4), xlim=(0, 2.3),
                title='Misma amplitud, diferentes espectros')
    axes[0].legend(fontsize=8)
    for i, row in enumerate(edge_rows):
        axes[1].loglog(etas, row['subgap_dos_at_half_edge'], 'o-', color=COLORS[i],
                       label=f'Γ / |Δ| = {row["Gamma_over_delta"]:g}')
    axes[1].set(xlabel=r'Regulador $\eta/|\Delta|$', ylabel=r'DOS en $E=E_g/2$',
                title='El regulador genera una cola')
    axes[1].legend(fontsize=8)
    save(fig, 'A_05_espectro_y_regulador')

    # Recompute the two earlier illustrations with this independent solver.
    T = .9/8.65
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), layout='constrained')
    ds = np.linspace(.15, 2.15, 150)
    for i, Q in enumerate((0., .55)):
        axes[0].plot(ds, [converged(T, d, Q)[0] for d in ds], color=COLORS[i], label=f'Q = {Q:.2f}')
        deq = eq_gap(T, Q)
        axes[0].plot(deq, converged(T, deq, Q)[0], 'o', color=COLORS[i])
    axes[0].set(xlabel=r'$|\Delta|/(k_BT_c)$', ylabel=r'$\delta f/[N_0(k_BT_c)^2]$', title='Pendiente en amplitud: fuerza')
    axes[0].legend(fontsize=8)
    qs = np.linspace(0, .65, 120)
    f0 = converged(T, 1.4, 0)[0]
    axes[1].plot(qs, [converged(T, 1.4, q)[0]-f0 for q in qs], color=COLORS[0])
    p = converged(T, 1.4, .43)
    xx = np.linspace(.28, .57, 20)
    axes[1].plot(xx, p[0]-f0+p[2]*(xx-.43), '--', color=COLORS[1], label='Tangente: corriente conjugada')
    axes[1].plot(.43, p[0]-f0, 'o', color=COLORS[1])
    axes[1].set(xlabel='$Q$ (amplitud fija)', ylabel='Incremento de energía normalizado', title='Pendiente en flujo: corriente')
    axes[1].legend(fontsize=8)
    save(fig, 'A_01_energia_comun')
    ratios_grid = np.linspace(0, 1.2, 401)
    gap_grid = np.maximum(0, 1-ratios_grid**(2/3))**1.5
    fig, ax = plt.subplots(figsize=(7.8, 3.6), layout='constrained')
    ax.plot(ratios_grid, gap_grid, color=COLORS[0], lw=2.5, label=r'$E_g/|\Delta|$')
    ax.axhline(1, color=COLORS[1], ls='--', label=r'$|\Delta|$ como referencia')
    ax.fill_between(ratios_grid, gap_grid, 1, color=COLORS[0], alpha=.1)
    ax.set(xlabel=r'$\Gamma/|\Delta|$', ylabel='Energía normalizada', title='Amplitud y borde espectral son distintos', ylim=(-.03, 1.12))
    ax.legend(fontsize=9)
    save(fig, 'A_03_gap_espectral')
    report = {'document':'A', 'revision':'0.4', 'host':platform.node(),
              'python':platform.python_version(), 'numpy':np.__version__, 'scipy':scipy.__version__,
              'scope':'Catálogo uniforme estático; sin transitorio, ajuste material ni imports de producción o de v0.3.',
              'independent_method':'Newton acotado para energía Matsubara renormalizada; sumas crudas y extrapolación en corte; raíces cuárticas filtradas con ecuación retardada original.',
              'point_count':len(points), 'T_over_Tc':[.05,.15,.4,.8,.98], 'delta_kBTc':[.25,.8,1.5,2.], 'Q':[0.,.15,.45,.75],
              'derivative_steps':steps.tolist(), 'derivative_error_scaling':'abs(numerical-reference)/max(1,abs(reference)); mixed uses max of both derivatives',
              'max_errors_by_step':max_by_h.tolist(), 'chosen_step':1e-3,
              'max_errors_at_chosen_step':dict(zip(['amplitude','current','mixed'],chosen.tolist())),
              'cutoffs':ns.tolist(), 'extrapolation_windows':[ns[-6:-1].tolist(),ns[-5:].tolist()],
              'max_extrapolation_window_change':np.max(convergence[:,3:],axis=0).tolist(),
              'max_matsubara_absolute_residual':spectral_max,
              'spectral_state_count':len(spectral_rows), 'max_retarded_scaled_residual':max_retarded,
              'max_retarded_normalization_residual':max_normalization, 'min_dos':min_dos,
              'spectral_edge_checks':edge_rows, 'runtime_seconds':time.perf_counter()-started}
    (DATA/'A_auditoria_independiente_v04.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2,ensure_ascii=False))


if __name__ == '__main__':
    main()
