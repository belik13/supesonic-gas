# plot_geometry.py — коническая (автомодельная) картина течения: лучи из точки взаимодействия.
# Дополнительный скрипт визуализации; исходный код ВКР не использует и не меняет.
import os, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
from matplotlib import cm, colors
from matplotlib.lines import Line2D
plt.rcParams.update({'font.size':10, 'font.family':'DejaVu Sans', 'figure.dpi':130})
RES = '../results/'
D = np.pi/180.0
R = 1.0

def ray(ax, ang_deg, style):
    a = ang_deg*D; rr = 1.10*R; x, y = rr*np.cos(a), rr*np.sin(a)         # чуть за край веера
    ax.plot([0,x],[0,y], '-', color='white', lw=5.0, zorder=4)           # белый ореол для видимости
    ax.plot([0,x],[0,y], '-' if style=='shock' else '--',
            color='black', lw=2.6, zorder=5)

def flow_arrow(ax, phi_rep, theta_deg, label=None):
    a = phi_rep*D; r = 0.60
    x, y = r*np.cos(a), r*np.sin(a)
    dx, dy = 0.17*np.cos(theta_deg*D), 0.17*np.sin(theta_deg*D)
    ax.arrow(x-dx/2, y-dy/2, dx, dy, width=0.008, head_width=0.05,
             fc='white', ec='black', lw=0.8, length_includes_head=True, zorder=6)
    if label:
        ax.text(x-0.10*np.cos(a), y-0.10*np.sin(a), label, color='black',
                fontsize=9, ha='center', va='center', zorder=7,
                bbox=dict(boxstyle='circle', fc='white', ec='black', lw=0.7))

def draw_fan(ax, csv, title):
    d = np.genfromtxt(RES+csv, delimiter=',', names=True)
    phi, p, T, th = d['phi'], d['p'], d['T'], d['theta']
    o = np.argsort(phi); phi, p, T, th = phi[o], p[o], T[o], th[o]
    pmin, pmax, Tmax = p.min(), p.max(), T.max()
    norm = colors.Normalize(vmin=pmin/1e5, vmax=pmax/1e5); cmap = cm.plasma
    # 1) заливка секторов цветом по давлению (на каждом луче параметры постоянны)
    for i in range(len(phi)-1):
        if phi[i+1]-phi[i] < 1e-3:  # разрыв — сектора нет
            continue
        c = cmap(norm(0.5*(p[i]+p[i+1])/1e5))
        ax.add_patch(Wedge((0,0), R, phi[i], phi[i+1], width=R, facecolor=c, edgecolor='none'))
    # 2) разрывы: скачок (сплошной) либо контакт (пунктир)
    for i in range(len(phi)-1):
        if phi[i+1]-phi[i] < 1e-3:
            mid = 0.5*(phi[i]+phi[i+1])
            if abs(p[i+1]-p[i])/pmax > 0.02:      ray(ax, mid, 'shock')
            elif abs(T[i+1]-T[i])/Tmax > 0.02:    ray(ax, mid, 'contact')
    # 3) направления потока: набегающие L, R и зона взаимодействия
    flow_arrow(ax, phi[1],  th[0],  'L')
    flow_arrow(ax, phi[-2], th[-1], 'R')
    m = len(phi)//2; flow_arrow(ax, phi[m], th[m], '*')
    ax.plot(0,0,'o', color='white', mec='black', ms=5, zorder=8)
    ax.set_xlim(-0.1,1.12); ax.set_ylim(-1.05,1.05); ax.set_aspect('equal'); ax.axis('off')
    ax.set_title(title, fontsize=10)
    sm = cm.ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
    cb = plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.02); cb.set_label('p, бар', fontsize=8)

fig, axes = plt.subplots(2, 2, figsize=(12, 10))
cases = [
    ('prof_SS.csv',      'Скачок–скачок (M=3, схождение ±12°)'),
    ('prof_RR.csv',      'Разрежение–разрежение (M=3, расхождение ±10°)'),
    ('prof_SR.csv',      'Скачок–разрежение (перепад давления 2:1)'),
    ('prof_custom2.csv', 'Несимметричный: смешанная конфигурация'),
]
for ax, (c, t) in zip(axes.flat, cases):
    if os.path.exists(RES+c): draw_fan(ax, c, t)
    else: ax.axis('off'); ax.set_title(t+' — нет данных', fontsize=9)

leg = [Line2D([0],[0], color='black', lw=2.0, label='скачок уплотнения (один луч)'),
       Line2D([0],[0], color='black', lw=2.0, ls='--', label='контактный разрыв'),
       Line2D([0],[0], color='white', marker='>', mec='black', lw=0, label='направление потока θ (L, R, *)')]
fig.legend(handles=leg, loc='lower center', ncol=3, fontsize=9, frameon=False)
fig.suptitle('Коническое (автомодельное) течение: лучи из точки взаимодействия φ', fontsize=13)
plt.tight_layout(rect=[0,0.05,1,0.96])
plt.savefig('12 Геометрия течения (веер лучей).png'); plt.close()
print('ok')
