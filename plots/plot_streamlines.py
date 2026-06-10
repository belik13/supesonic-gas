# plot_streamlines.py — линии тока: как газ втекает слева, проходит волны и вытекает справа.
# Дополнительный скрипт визуализации; исходный код ВКР не использует и не меняет.
# Поле строится из автомодельного решения: параметры зависят только от угла φ=atan2(y,x).
import os, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10, 'font.family':'DejaVu Sans', 'figure.dpi':130})
RES = '../results/'
D = np.pi/180.0

def draw(ax, csv, title):
    d = np.genfromtxt(RES+csv, delimiter=',', names=True)
    phi, p, T, V, th = d['phi'], d['p'], d['T'], d['V'], d['theta']
    o = np.argsort(phi); phi, p, T, V, th = phi[o], p[o], T[o], V[o], th[o]
    n = 420
    xs = np.linspace(-1.0, 1.45, n); ys = np.linspace(-1.05, 1.05, n)
    X, Y = np.meshgrid(xs, ys)
    PHI = np.degrees(np.arctan2(Y, X))            # автомодельный угол каждой точки
    # вне диапазона φ np.interp «зажимает» к краям = однородные набегающие потоки L и R
    Pg = np.interp(PHI, phi, p); THg = np.interp(PHI, phi, th); Vg = np.interp(PHI, phi, V)
    U = Vg*np.cos(THg*D); W = Vg*np.sin(THg*D)    # вектор скорости
    im = ax.pcolormesh(X, Y, Pg/1e5, cmap='plasma', shading='auto')
    # линии тока из «затравок» вдоль левой границы — это и есть набегающий газ
    sy = np.linspace(-0.95, 0.95, 21)
    seeds = np.column_stack([np.full_like(sy, -0.98), sy])
    ax.streamplot(xs, ys, U, W, color='white', density=2.2, linewidth=1.0,
                  arrowsize=1.0, start_points=seeds, zorder=2)
    # наложить волны: видно, как линии тока преломляются на скачке (сплошной) и контакте (пунктир)
    pmax, Tmax, Rm = p.max(), T.max(), 1.75
    for i in range(len(phi)-1):
        if phi[i+1]-phi[i] < 1e-3:
            a = 0.5*(phi[i]+phi[i+1])*D; x, y = Rm*np.cos(a), Rm*np.sin(a)
            if abs(p[i+1]-p[i])/pmax > 0.02:
                ax.plot([0,x],[0,y], '-', color='white', lw=4.5, zorder=4)
                ax.plot([0,x],[0,y], '-', color='black', lw=2.6, zorder=5)
            elif abs(T[i+1]-T[i])/Tmax > 0.02:
                ax.plot([0,x],[0,y], '-',  color='white', lw=4.0, zorder=4)
                ax.plot([0,x],[0,y], '--', color='black', lw=2.2, zorder=5)
    ax.plot(0, 0, 'o', color='cyan', mec='black', ms=6, zorder=6)
    ax.set_aspect('equal'); ax.set_xlim(-1.0, 1.45); ax.set_ylim(-1.05, 1.05)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title, fontsize=10)
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.02); cb.set_label('p, бар', fontsize=8)

fig, axes = plt.subplots(2, 2, figsize=(12.5, 10))
cases = [
    ('prof_SS.csv',      'Скачок–скачок: потоки сходятся → тормозятся, сжимаются'),
    ('prof_RR.csv',      'Разрежение–разрежение: потоки расходятся → ускоряются'),
    ('prof_SR.csv',      'Скачок–разрежение: перепад давления 2:1'),
    ('prof_custom2.csv', 'Несимметричный: смешанная конфигурация'),
]
for ax, (c, t) in zip(axes.flat, cases):
    if os.path.exists(RES+c): draw(ax, c, t)
    else: ax.axis('off')
fig.suptitle('Линии тока: газ втекает слева, проходит волны (скачки/веера) и вытекает справа',
             fontsize=12)
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig('13 Линии тока (как текут потоки).png'); plt.close()
print('ok')
