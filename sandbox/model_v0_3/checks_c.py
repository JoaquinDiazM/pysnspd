"""Reproducible, lightweight illustrative checks for document C revision 0.3.

No production solver, material tables, or detector transient is used.
Requires numpy, scipy and matplotlib. Outputs under docs/modelo_v0_3.
"""
from pathlib import Path
import csv
import json
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'docs/modelo_v0_3'
FIG = OUT / 'figuras'
DATA = OUT / 'verificaciones'
FIG.mkdir(parents=True, exist_ok=True)
DATA.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 10,
    'axes.labelsize': 10, 'axes.titlesize': 11,
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.titleweight': 'bold', 'legend.frameon': False,
    'grid.alpha': .20, 'savefig.dpi': 200,
})
BLUE, RED, GOLD, GREEN = '#235789', '#c4503c', '#b27a19', '#278473'

def save_csv(name, columns, rows):
    with (DATA / name).open('w', encoding='utf-8', newline='') as f:
        w = csv.writer(f)
        w.writerow(columns)
        w.writerows(rows)

def save_plot(fig, name):
    fig.savefig(FIG / name, facecolor='white', bbox_inches='tight')
    fig.savefig((FIG / name).with_suffix('.pdf'), facecolor='white', bbox_inches='tight')
    plt.close(fig)

def quadrature(order):
    nodes, weights = leggauss(order)
    return (nodes+1)*20, weights*20

X, W = quadrature(512)
X_FINE, W_FINE = quadrature(1024)

def bcs_u(temp, gap, x=X, w=W):
    en = np.hypot(x, gap)
    exp_minus = np.exp(-en / temp)
    return float(np.dot(w, en * exp_minus/(1+exp_minus)))

def bcs_cap(temp, gap):
    en = np.hypot(X, gap)
    em = np.exp(-en/temp)
    return float(np.dot(W, en**2/temp**2*em/(1+em)**2))

# Figure C.1: A conditional coordinate map, complementary to B's u(T) plots.
gaps = np.linspace(0, 1.2, 145)
levels = [.02, .06, .15, .30]
curves = []
fig, ax = plt.subplots(figsize=(8.7, 4.7), layout='constrained')
for lev, color in zip(levels, [BLUE, RED, GOLD, GREEN]):
    temps = np.array([brentq(lambda tt: bcs_u(tt, gg)-lev, .01, 2., xtol=1e-13) for gg in gaps])
    curves.append(temps)
    ax.plot(gaps, temps, color=color, lw=2.7 if lev == .06 else 1.7,
            label=rf'$u_{{\rm qp}}/(4N_0\Delta_{{\rm ref}}^2)={lev:.2f}$')
point_gaps = [.2, .9]
point_temps = [brentq(lambda tt: bcs_u(tt, gg)-.06, .01, 2.) for gg in point_gaps]
for gg, tt in zip(point_gaps, point_temps):
    ax.scatter(gg, tt, s=58, color=RED, edgecolors='white', zorder=4)
ax.annotate(f'Misma energía\nTemperaturas: {point_temps[0]:.3f} y {point_temps[1]:.3f}',
            xy=(.9, point_temps[1]), xytext=(.47, .12),
            arrowprops={'arrowstyle': '->', 'color': RED}, color=RED, fontsize=10)
ax.set(xlabel=r'Amplitud instantánea $|\Delta|/\Delta_{\rm ref}$',
       ylabel=r'Temperatura equivalente $k_BT_E/\Delta_{\rm ref}$',
       title='La inversa energía–temperatura depende del espectro', xlim=(0, 1.2), ylim=(.07, .87))
ax.grid(True)
ax.legend(loc='upper left', ncol=2, fontsize=9)
save_plot(fig, 'C_01_coordenadas_energia.png')
save_csv('C_coordenadas_energia.csv', ['gap_Dref']+[f'T_for_u_{v:g}' for v in levels],
         zip(gaps, *curves))
quad_err = max(abs(bcs_u(tt, gg)-bcs_u(tt, gg, X_FINE, W_FINE)) for gg in [0., .2, .9, 1.2] for tt in [.1, .3, .7])
cap_err = max(abs((bcs_u(tt+1e-5, gg)-bcs_u(tt-1e-5, gg))/(2e-5)-bcs_cap(tt, gg))
              for gg in [0., .2, .9] for tt in [.15, .3, .7])
assert quad_err < 1e-9, quad_err
assert cap_err < 1e-8, cap_err

def rk4(rhs, y0, end, step):
    n = int(round(end/step))
    ts = np.linspace(0., end, n+1)
    ys = np.empty((n+1, len(y0)))
    ys[0] = y0
    for k in range(n):
        t, y = ts[k], ys[k]
        a = rhs(t, y)
        b = rhs(t+step/2, y+step*a/2)
        c = rhs(t+step/2, y+step*b/2)
        d = rhs(t+step, y+step*c)
        ys[k+1] = y+step*(a+2*b+2*c+d)/6
    return ts, ys

# Figure C.2: A thermodynamic toy with a known energy/free-energy relation.
def order_energy(d): return (1-d*d)**2/4
def coeff(d): return 1-.4*d*d
def energy(d, temp): return order_energy(d)+coeff(d)*temp*temp/2
def force(d, temp): return d*(d*d-1+.4*temp*temp)

d0, temp0, step = .35, .2, .004
u0 = energy(d0, temp0)
def temp_from_energy(d):
    rad = 2*(u0-order_energy(d))/coeff(d)
    if rad <= 0: raise ValueError('Energy outside admissible positive-temperature range')
    return np.sqrt(rad)

def energy_rhs(t, y):
    d = y[0]
    return np.array([-force(d, temp_from_energy(d))])

def thermal_rhs(t, y):
    d, temp = y
    ff = force(d, temp)
    dd = -ff
    dt = (ff*ff+.8*d*temp*temp*dd)/(coeff(d)*temp)
    return np.array([dd, dt])

t, de = rk4(energy_rhs, [d0], 12., step)
_, thermal = rk4(thermal_rhs, [d0, temp0], 12., step)
_, fine = rk4(thermal_rhs, [d0, temp0], 12., step/2)
_, naive = rk4(lambda tt, yy: np.array([-force(yy[0], temp0)]), [d0], 12., step)
d = de[:, 0]
temp = np.array([temp_from_energy(v) for v in d])
uo = order_energy(d)
uq = coeff(d)*temp*temp/2
thermal_energy = energy(thermal[:, 0], thermal[:, 1])
un = energy(naive[:, 0], temp0)
entropy = coeff(d)*temp
entres = []
for dd, tt in thermal:
    derivatives = thermal_rhs(0., [dd, tt])
    entres.append(tt*(-.8*dd*tt*derivatives[0]+coeff(dd)*derivatives[1])-force(dd, tt)**2)
cell_check = {
    'model': 'Illustrative polynomial free energy C.20b, dimensionless',
    'initial': {'d': d0, 'temperature': temp0, 'energy': float(u0)},
    'final_at_t12': {'d': float(d[-1]), 'temperature': float(temp[-1]), 'entropy': float(entropy[-1])},
    'max_energy_drift_thermal_coordinates': float(np.max(abs(thermal_energy-u0))),
    'max_thermal_vs_energy_coordinates': float(max(np.max(abs(thermal[:, 0]-d)), np.max(abs(thermal[:, 1]-temp)))),
    'max_change_halving_step': float(np.max(abs(thermal-fine[::2]))),
    'max_entropy_identity_residual': float(np.max(np.abs(entres))),
    'initial_entropy': float(entropy[0]),
    'min_entropy_increment': float(np.min(np.diff(entropy))),
    'unbalanced_final_energy_loss_fraction': float((u0-un[-1])/u0),
}
assert cell_check['max_energy_drift_thermal_coordinates'] < 1e-8, cell_check
assert cell_check['max_thermal_vs_energy_coordinates'] < 1e-7, cell_check
assert cell_check['max_change_halving_step'] < 1e-7, cell_check
assert cell_check['min_entropy_increment'] > -1e-12, cell_check
fig, ax = plt.subplots(1, 3, figsize=(11.0, 3.7), layout='constrained')
ax[0].plot(t, d, color=BLUE, label='Amplitud')
ax[0].plot(t, temp, color=RED, label='Temperatura')
ax[0].set(title='El condensado se recupera', ylabel='Valor adimensional')
ax[0].legend(fontsize=9)
ax[1].plot(t, uo, color=BLUE, label='Orden')
ax[1].plot(t, uq, color=RED, label='Excitaciones')
ax[1].plot(t, uo+uq, color=GREEN, lw=2, label='Total')
ax[1].set(title='La energía se transfiere', ylabel='Energía adimensional')
ax[1].legend(fontsize=9)
ax[2].plot(t, thermal_energy/u0, color=GREEN, lw=2, label='Balance completo')
ax[2].plot(t, un/u0, color=GOLD, ls='--', label='Temperatura fija')
ax[2].set(title='Prueba de conservación', ylabel=r'$u(t)/u(0)$', ylim=(0, 1.09))
ax[2].legend(fontsize=9)
for aa in ax: aa.set_xlabel('Tiempo ilustrativo'); aa.grid(True)
save_plot(fig, 'C_02_celda_balance.png')
save_csv('C_celda_balance.csv', ['t', 'd_energy', 'T_energy', 'u_order', 'u_exc', 'd_thermal', 'T_thermal', 'u_thermal', 'entropy', 'u_naive'],
         zip(t, d, temp, uo, uq, thermal[:, 0], thermal[:, 1], thermal_energy, entropy, un))

# Figure C.3: Prescribed resistance, not a computed detector event.
def resistance(tt): return 4*(1-np.exp(-tt/.08))*np.exp(-tt/.5)
step_c = .002
circuits, circuit_data = [], []
fig, axes = plt.subplot_mosaic([['r', 'r'], ['i', 'v']], figsize=(9.0, 6.3),
                               height_ratios=[.65, 1], layout='constrained')
for lam, color in zip([.25, 1., 4.], [BLUE, RED, GREEN]):
    def rhs(tt, yy): return (1-(1+resistance(tt))*yy)/lam
    tt, yy = rk4(rhs, [1.], 6., step_c)
    tf, yf = rk4(rhs, [1.], 6., step_c/2)
    i, r = yy[:, 0], resistance(tt)
    vp, vpatch = 1-i, r*i
    # Independent quadrature checks the integrated circuit power identity.
    power_net = vp-r*i*i-vp*vp
    energy_change = lam*(i[-1]**2-i[0]**2)/2
    power_integral = np.trapezoid(power_net, tt)
    peak = int(np.argmax(vp))
    patch_peak = int(np.argmax(vpatch))
    item = {
        'lambda': lam, 'peak_port_time': float(tt[peak]), 'peak_port_voltage': float(vp[peak]),
        'peak_patch_time': float(tt[patch_peak]), 'peak_patch_voltage': float(vpatch[patch_peak]),
        'max_current_change_halving_step': float(np.max(abs(i-yf[::2, 0]))),
        'integrated_energy_balance_absolute_residual': float(abs(power_integral-energy_change)),
    }
    assert item['max_current_change_halving_step'] < 1e-7, item
    assert item['integrated_energy_balance_absolute_residual'] < 2e-5, item
    circuits.append(item)
    axes['i'].plot(tt, i, color=color, label=rf'$\lambda={lam:g}$')
    axes['v'].plot(tt, vp, color=color, label=rf'Puerto, $\lambda={lam:g}$')
    axes['v'].scatter(tt[peak], vp[peak], color=color, s=22)
    if lam == 1.:
        axes['v'].plot(tt, vpatch, color=RED, ls='--', lw=1.3, label=r'Región, $\lambda=1$')
    circuit_data.extend(zip(np.full_like(tt, lam), tt, r, i, vpatch, vp))
axes['r'].plot(tt, resistance(tt), color=GOLD, lw=2)
axes['r'].set(title='Entrada prescrita común a los tres circuitos', ylabel=r'$R_{\rm patch}/R_L$')
axes['i'].set(title='La corriente se desvía a la carga', ylabel=r'$I/I_b$')
axes['v'].set(title='Región y puerto responden distinto', ylabel=r'$V/(I_bR_L)$')
for aa in axes.values():
    aa.set_xlim(0, 4)
    aa.set_xlabel(r'Tiempo $\zeta=tR_L/L_{\rm ref}$')
    aa.grid(True)
axes['i'].legend(fontsize=9)
axes['v'].legend(fontsize=8.5)
save_plot(fig, 'C_03_circuito_prescrito.png')
save_csv('C_circuito_prescrito.csv', ['lambda', 'time_zeta', 'r_patch', 'current', 'voltage_patch', 'voltage_port'], circuit_data)

summary = {
    'revision': '0.3', 'date': '2026-09-08',
    'scope': 'BCS coordinate map; illustrative cell; prescribed circuit. No production detector transient.',
    'bcs': {'quadrature_order': 512, 'quadrature_upper_x_Dref': 40,
            'max_absolute_change_order_1024': quad_err, 'max_heat_capacity_derivative_error': cap_err,
            'equal_energy': .06, 'equal_energy_gaps': point_gaps, 'equal_energy_temperatures': point_temps},
    'cell': cell_check,
    'circuit': circuits,
    'diffusion_C24': {'D_m2_s': 1.581e-4, 'epsilon': .01, 'times_ps': [5, 20, 50],
                      'lengths_nm': [float(np.sqrt(4*1.581e-4*v*1e-12*np.log(100))*1e9) for v in [5, 20, 50]]},
    'gradient_C28': {'K_ratio': 1/1.764, 'length_ratio': float(np.sqrt(1/1.764))},
}
(DATA / 'C_checks.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
print(json.dumps(summary, indent=2, ensure_ascii=False))
