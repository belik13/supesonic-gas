# plot_prof.py — автомодельные профили трёх конфигураций (читает results/, пишет PNG рядом).
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':11,'font.family':'DejaVu Sans','figure.dpi':120})

RES = '../results/'

def plot(name, title, out):
    d = np.genfromtxt(RES + name, delimiter=',', names=True)
    fig, ax = plt.subplots(2, 2, figsize=(11, 7))
    ax[0, 0].plot(d['phi'], d['p']/1e5, 'b-'); ax[0, 0].set_ylabel('p, бар'); ax[0, 0].grid(alpha=.3)
    ax[0, 1].plot(d['phi'], d['M'], 'r-');     ax[0, 1].set_ylabel('M');     ax[0, 1].grid(alpha=.3)
    ax[1, 0].plot(d['phi'], d['T'], 'g-');     ax[1, 0].set_ylabel('T, К');  ax[1, 0].set_xlabel('φ, град'); ax[1, 0].grid(alpha=.3)
    ax[1, 1].plot(d['phi'], d['theta'], 'm-'); ax[1, 1].set_ylabel('θ, град'); ax[1, 1].set_xlabel('φ, град'); ax[1, 1].grid(alpha=.3)
    fig.suptitle(title); plt.tight_layout(); plt.savefig(out); plt.close()

plot('prof_SS.csv', 'Конфигурация скачок–скачок (схождение ±12°, M=3)', '01 Профиль скачок-скачок.png')
plot('prof_RR.csv', 'Конфигурация разрежение–разрежение (расхождение ±10°, M=3)', '02 Профиль разрежение-разрежение.png')
plot('prof_SR.csv', 'Конфигурация скачок — разрежение (перепад давления 2:1, M=3)', '03 Профиль скачок-разрежение.png')
import os
custom = [
    ('prof_custom1.csv', 'Случай 1: разрежение–разрежение (расхождение ±18°, M=3)', '04 Профиль пользовательский 1.png'),
    ('prof_custom2.csv', 'Случай 2: несимметричные потоки (смешанная конфигурация)', '05 Профиль пользовательский 2.png'),
    ('prof_custom3.csv', 'Случай 3: режим Маха (решения нет)', '06 Профиль пользовательский 3.png'),
]
for n, t, o in custom:
    if os.path.exists(RES + n):
        plot(n, t, o)
print("ok")
