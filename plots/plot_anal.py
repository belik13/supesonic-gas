# plot_anal.py — карта существования, зависимости углов, сравнение моделей газа.
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':11,'font.family':'DejaVu Sans','figure.dpi':125})

RES = '../results/'

# 1. карта областей существования (заливка по граничным кривым)
b = np.genfromtxt(RES + 'exist_bounds.csv', delimiter=',', names=True)
M = b['M']; mach = b['dth_mach']; vac = b['dth_vac']
fig, ax = plt.subplots(figsize=(8.6, 5.2)); ytop = 100; ybot = -175
ax.fill_between(M, mach, ytop, color='#c9c9ef', label='нет регулярного решения (режим Маха)')
ax.fill_between(M, 0,    mach, color='#fb8b6c', label='скачок – скачок')
ax.fill_between(M, vac,  0,    color='#9ecae1', label='разрежение – разрежение')
ax.fill_between(M, ybot, vac,  color='#dddddd', label='вакуум (решение отсутствует)')
ax.plot(M, mach, 'k-', lw=2.0); ax.plot(M, vac, 'k--', lw=2.0); ax.axhline(0, color='gray', lw=0.6)
ax.set_xlim(M.min(), M.max()); ax.set_ylim(ybot, ytop)
ax.set_xlabel('число Маха набегающих потоков M'); ax.set_ylabel('Δθ = θ_L − θ_R, град')
ax.set_title('Области существования и типы конфигураций распада разрыва')
ax.text(5.5,  72, 'нет регулярного решения (Маха)', ha='center', fontsize=9)
ax.text(5.5,  28, 'скачок – скачок', ha='center', fontsize=9)
ax.text(5.5, -28, 'разрежение – разрежение', ha='center', fontsize=9)
ax.text(5.5,-120, 'вакуум (нет решения)', ha='center', fontsize=9)
plt.tight_layout(); plt.savefig('07 Карта областей существования.png'); plt.close()

# 2. углы при симметричном схождении vs alpha
d = np.genfromtxt(RES + 'ang_sym.csv', delimiter=',', names=True); m = d['sigma_c'] > 0
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
ax[0].plot(d['alpha'][m], d['pstar_c'][m], 'b-',  label='γ = const')
ax[0].plot(d['alpha'][m], d['pstar_v'][m], 'r--', label='c_p(T)')
ax[0].set_xlabel('полуугол схождения α, град'); ax[0].set_ylabel('p*/p')
ax[0].grid(alpha=.3); ax[0].legend(); ax[0].set_title('Давление в зоне взаимодействия')
ax[1].plot(d['alpha'][m], d['sigma_c'][m], 'b-',  label='γ = const')
ax[1].plot(d['alpha'][m], d['sigma_v'][m], 'r--', label='c_p(T)')
ax[1].set_xlabel('полуугол схождения α, град'); ax[1].set_ylabel('σ, град')
ax[1].grid(alpha=.3); ax[1].legend(); ax[1].set_title('Угол наклона косого скачка')
plt.tight_layout(); plt.savefig('08 Зависимость от угла схождения.png'); plt.close()

# 3. контактный разрыв и скачок vs перепад давления
d = np.genfromtxt(RES + 'ang_pratio.csv', delimiter=',', names=True)
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
ax[0].plot(d['ratio'], d['thetaStar'], 'm-'); ax[0].axhline(0, color='gray', lw=.7)
ax[0].set_xlabel('p_L / p_R'); ax[0].set_ylabel('θ* (наклон контактного разрыва), град')
ax[0].grid(alpha=.3); ax[0].set_title('Наклон контактного разрыва')
ax[1].plot(d['ratio'], d['sigma_shock'], 'b-')
ax[1].set_xlabel('p_L / p_R'); ax[1].set_ylabel('σ, град')
ax[1].grid(alpha=.3); ax[1].set_title('Угол скачка (со стороны меньшего давления)')
plt.tight_layout(); plt.savefig('09 Зависимость от перепада давления.png'); plt.close()

# 4. vs M
d = np.genfromtxt(RES + 'ang_vsM.csv', delimiter=',', names=True)
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
ax[0].plot(d['M'], d['pstar_c'], 'b-',  label='γ = const')
ax[0].plot(d['M'], d['pstar_v'], 'r--', label='c_p(T)')
ax[0].set_xlabel('M'); ax[0].set_ylabel('p*/p'); ax[0].grid(alpha=.3); ax[0].legend()
ax[0].set_title('p*/p при α=20°')
ax[1].plot(d['M'], d['sigma_c'], 'b-',  label='γ = const')
ax[1].plot(d['M'], d['sigma_v'], 'r--', label='c_p(T)')
ax[1].set_xlabel('M'); ax[1].set_ylabel('σ, град'); ax[1].grid(alpha=.3); ax[1].legend()
ax[1].set_title('Угол скачка σ при α=20°')
plt.tight_layout(); plt.savefig('10 Зависимость от числа Маха.png'); plt.close()

# 5. сравнение профилей cp(T) vs const (T и M)
dc = np.genfromtxt(RES + 'prof_cmp_const.csv', delimiter=',', names=True)
da = np.genfromtxt(RES + 'prof_cmp_air.csv',   delimiter=',', names=True)
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
ax[0].plot(dc['phi'], dc['T'], 'b-',  label='γ = const')
ax[0].plot(da['phi'], da['T'], 'r--', label='c_p(T)')
ax[0].set_xlabel('φ, град'); ax[0].set_ylabel('T, К'); ax[0].grid(alpha=.3); ax[0].legend()
ax[0].set_title('Температура (СС, M=4, α=18°)')
ax[1].plot(dc['phi'], dc['M'], 'b-',  label='γ = const')
ax[1].plot(da['phi'], da['M'], 'r--', label='c_p(T)')
ax[1].set_xlabel('φ, град'); ax[1].set_ylabel('M'); ax[1].grid(alpha=.3); ax[1].legend()
ax[1].set_title('Число Маха')
plt.tight_layout(); plt.savefig('11 Сравнение моделей газа.png'); plt.close()
print("plots done")
