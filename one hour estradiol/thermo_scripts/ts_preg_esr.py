"""Time series: two-sided model vs the pregnant and Esr1i mice."""
import numpy as np
from scipy.optimize import minimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
exec(open('freephase.py').read().split('ALL = np.array')[0])

ALL = np.array([0, 1, 2])
rng = np.random.default_rng(31)
spec = MODELS['dual       ']
bnds = spec['bounds'] + TAIL
lo = np.array([b[0] for b in bnds]); hi = np.array([b[1] for b in bnds])
best = None
for _ in range(200):
    x0 = lo + rng.random(len(bnds)) * (hi - lo)
    r = minimize(lambda z: loss(spec, z, ALL), x0, bounds=bnds,
                 method='L-BFGS-B', options=dict(maxiter=300, ftol=1e-14))
    if best is None or r.fun < best[0]:
        best = (float(r.fun), r.x)
e, x = best
M = model(spec, x, ALL)
np.save('best_dual_freephase.npy', x)

SD = np.vstack([fold[l][2] for l in LABELS])
show = [1, 2]                                   # pregnant, Esr1i
titles = ['Pregnant-CTL (Di36)', 'Pregnant-Esr1i (Di7)']

fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))
for ax, i, ti in zip(axes, show, titles):
    ax.fill_between(hrs, OBS[i] - SD[i], OBS[i] + SD[i], color='0.75',
                    alpha=0.45, lw=0, label='data +/- SD across days')
    ax.plot(hrs, OBS[i], 'o-', ms=4, lw=1.2, color='0.25', label='mouse (digitised)')
    ax.plot(hrs, M[i], '-', lw=2.6, color='C3', label='two-sided model')
    a_o = OBS[i].max() - OBS[i].min()
    a_m = M[i].max() - M[i].min()
    ax.set_title(f"{ti}\namplitude: mouse {a_o:.2f} C   model {a_m:.2f} C", fontsize=10)
    ax.set_xlabel('hour of day'); ax.grid(alpha=0.3); ax.set_xlim(0, 24)
    ax.set_xticks(range(0, 25, 4))
axes[0].set_ylabel('core body temperature (C)')
axes[0].legend(fontsize=8, loc='upper left')
lo_y = min(OBS[show].min(), M[show].min()) - 0.15
hi_y = max(OBS[show].max(), M[show].max()) + 0.15
for ax in axes:
    ax.set_ylim(lo_y, hi_y)
plt.tight_layout()
plt.savefig('fig_ts_preg_esr.png', dpi=145)

# same but all three, shared axis, for context
fig, axes = plt.subplots(1, 3, figsize=(15, 4.4), sharey=True)
for ax, i, ti in zip(axes, [0, 1, 2],
                     ['Non-pregnant (Cir8)'] + titles):
    ax.fill_between(hrs, OBS[i]-SD[i], OBS[i]+SD[i], color='0.75', alpha=0.45, lw=0)
    ax.plot(hrs, OBS[i], 'o-', ms=3.5, lw=1.1, color='0.25', label='mouse')
    ax.plot(hrs, M[i], '-', lw=2.4, color='C3', label='model')
    ax.set_title(f"{ti}   (mouse {OBS[i].max()-OBS[i].min():.2f} / "
                 f"model {M[i].max()-M[i].min():.2f} C)", fontsize=9.5)
    ax.set_xlabel('hour of day'); ax.grid(alpha=0.3); ax.set_xticks(range(0, 25, 6))
axes[0].set_ylabel('core body temperature (C)'); axes[0].legend(fontsize=8)
plt.tight_layout(); plt.savefig('fig_ts_all3.png', dpi=145)
print(f"loss {e:.4f}")
for i, l in enumerate(LABELS):
    print(f"  {l:9s} mouse amp {OBS[i].max()-OBS[i].min():.2f}  "
          f"model amp {M[i].max()-M[i].min():.2f}  "
          f"mouse mean {OBS[i].mean():.2f}  model mean {M[i].mean():.2f}")
print("wrote fig_ts_preg_esr.png, fig_ts_all3.png")
