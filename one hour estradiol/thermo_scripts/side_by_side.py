"""
Model beside the original figure panels, on the SAME time axis (days 13-18).

Top row  : untouched crop of each CBT-MA panel from mousy temp.png.
Bottom   : the fitted two-sided model, its 24 h waveform tiled across the same
           5 days, with the digitised mouse trace behind it for comparison.
Axes limits are matched panel-by-panel to the original (they differ!).
"""
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import minimize
exec(open('freephase.py').read().split('ALL = np.array')[0])

ALL = np.array([0, 1, 2])
x = np.load('best_dual_freephase.npy')
spec = MODELS['dual       ']
M = model(spec, x, ALL)                       # (3, 48) over 24 h

SRC = r"C:\Users\holly\source\garfinkel-xiong-menstrual-cycle\one hour estradiol\mousy temp.png"
img = np.asarray(Image.open(SRC).convert('RGB'))
Rc, Gc, Bc = [img[..., k].astype(float) for k in range(3)]
red = (Rc > 110) & (Rc - Gc > 45) & (Rc - Bc > 45)
red[:250, :] = False

XSPAN = [(116, 662), (739, 1285), (1362, 1907)]
YANCH = [[(326.5, 38.0), (365.5, 37.0), (403.5, 36.0)],
         [(312.5, 38.0), (367.5, 37.0), (423.0, 36.0)],
         [(320.5, 38.0), (372.0, 37.0), (423.0, 36.0)]]
ROW_TOP, ROW_BOT = 313, 423
TITLES = ['Non-pregnant (Cir8)', 'Pregnant-CTL (Di36)', 'Pregnant-Esr1i (Di7)']

days5 = np.linspace(13, 18, 5 * 48, endpoint=False)

fig, axes = plt.subplots(2, 3, figsize=(16, 6.6))
for i in range(3):
    x0, x1 = XSPAN[i]
    rows, temps = zip(*YANCH[i])
    m, c = np.polyfit(rows, temps, 1)
    ytop, ybot = m * ROW_TOP + c, m * ROW_BOT + c

    # ---- top: the original panel, untouched
    crop = img[ROW_TOP:ROW_BOT + 1, x0:x1 + 1]
    axes[0, i].imshow(crop, extent=[13, 18, ybot, ytop], aspect='auto')
    axes[0, i].set_title(f"{TITLES[i]}  -- original figure", fontsize=10)
    axes[0, i].set_ylim(ybot, ytop)

    # ---- bottom: model over the same 5 days, mouse behind it
    xs, ys = [], []
    for xc in range(x0, x1 + 1):
        col = np.nonzero(red[:, xc])[0]
        if col.size:
            xs.append(xc); ys.append(col.mean())
    tmouse = 13.0 + 5.0 * (np.array(xs, float) - x0) / (x1 - x0)
    Tmouse = m * np.array(ys) + c
    axes[1, i].plot(tmouse, Tmouse, lw=1.0, color='0.55', label='mouse')
    axes[1, i].plot(days5, np.tile(M[i], 5), lw=1.8, color='C3', label='model')
    axes[1, i].set_ylim(ybot, ytop); axes[1, i].set_xlim(13, 18)
    axes[1, i].grid(alpha=0.3)
    axes[1, i].set_title(f"model (amp {M[i].max()-M[i].min():.2f} C) vs mouse "
                         f"(amp {Tmouse.max()-Tmouse.min():.2f} C raw)", fontsize=9)
    axes[1, i].set_xlabel('day')
    for ax in (axes[0, i], axes[1, i]):
        ax.set_xticks(range(13, 19))
axes[0, 0].set_ylabel('CBT-MA (C)')
axes[1, 0].set_ylabel('CBT-MA (C)')
axes[1, 0].legend(fontsize=8, loc='lower left')
plt.tight_layout()
plt.savefig('fig_side_by_side.png', dpi=145)
print("wrote fig_side_by_side.png")
for i in range(3):
    print(f"  {LABELS[i]:9s} model daily amp {M[i].max()-M[i].min():.2f} C")
