"""REVERT to Holly's original form and establish an honest baseline against the RAW data.

    c(t) = cos(2*pi*(t-12)/24)          ONE clock, shared, PEAK LOCKED AT NOON
    loss = sig(bL + gL*c + wL*E2 + pL*P4)
    keep = sig(bK + gK*c + wK*E2 + pK*P4)
    CBT  = T0 + A_T*(keep - loss)
ALL COEFFICIENTS POSITIVE. E2 = 0/1/0, P4 = 0/line/line. No lambda, no two-phase, no dphi.

Judged against the RAW recording, not the 5-day fold that flattered every previous fit:
correlation per animal, and RMS relative to the animal's OWN day-to-day SD (the standard that
matters -- if RMS > that SD, the model predicts the animal worse than another day of the same
animal would).
"""
import numpy as np
from scipy.signal import find_peaks
from scipy.optimize import minimize

s = np.load('series.npy', allow_pickle=True).item()
fold = np.load('fold.npy', allow_pickle=True).item()
L = ['non-preg', 'pregnant', 'Esr1i']
hrs = fold[L[0]][0]
OBS = np.vstack([fold[l][1] for l in L])
AMP = OBS.max(1) - OBS.min(1)
E2 = np.array([0., 1., 0.])
midd = np.array([13, 14, 15, 16, 17], float) + 0.5
p4b = np.maximum(0., 1. - (midd - 13.) / 6.5).mean()
P4v = np.array([0., p4b, p4b])
c48 = np.cos(2 * np.pi * (hrs - 12.) / 24.)

# per-animal day-to-day SD = the bar the model must beat
SD = []
for l in L:
    t, y = s[l]
    D = []
    for d in range(13, 18):
        m = (t >= d) & (t < d + 1)
        if m.sum() >= 20:
            D.append(np.interp(hrs, (t[m] - d) * 24., y[m]))
    SD.append(float(np.mean(np.vstack(D).std(0))))
SD = np.array(SD)


def sig(z):
    return 1 / (1 + np.exp(-np.clip(z, -60, 60)))


def W_of(p):
    bK, gK, wK, pK, bL, gL, wL, pL, A_T = p
    keep = sig(bK + gK * c48[None, :] + wK * E2[:, None] + pK * P4v[:, None])
    loss = sig(bL + gL * c48[None, :] + wL * E2[:, None] + pL * P4v[:, None])
    return A_T * (keep - loss)


def score(p):
    if min(p[1], p[2], p[3], p[5], p[6], p[7]) < 0 or not (0 < p[8] <= 10):
        return 1e6
    W = W_of(p)
    W = W + (OBS.mean() - W.mean())
    a = W.max(1) - W.min(1)
    if a[0] < 1e-6:
        return 1e6
    rms = float(np.mean(np.sqrt(np.mean((W - OBS) ** 2, 1)) / AMP))
    aerr = float(np.mean(np.abs(a - AMP) / AMP))
    return rms + 0.5 * aerr


rng = np.random.default_rng(101)
N = 250_000
P = np.column_stack([rng.uniform(-15, 6, N), rng.uniform(0, 30, N), rng.uniform(0, 25, N),
                     rng.uniform(0, 15, N), rng.uniform(-15, 6, N), rng.uniform(0, 30, N),
                     rng.uniform(0, 25, N), rng.uniform(0, 15, N), rng.uniform(0.5, 8, N)])
vals = np.array([score(q) for q in P])
al = np.flatnonzero(vals < 1e5)
best = (float(vals[al].min()), P[al[np.argmin(vals[al])]].copy())
for i in al[np.argsort(vals[al])][:70]:
    z = P[i].copy()
    for _ in range(4):
        r = minimize(score, z, method='Nelder-Mead',
                     options=dict(maxiter=6000, maxfev=6000, xatol=1e-10, fatol=1e-13))
        z = r.x
    if float(r.fun) < best[0]:
        best = (float(r.fun), z.copy())

v, p = best
W = W_of(p); W = W + (OBS.mean() - W.mean())
print(f"=== HOLLY'S ORIGINAL FORM (noon-locked, all-positive, no lambda)  obj {v:.4f} ===",
      flush=True)
print(f'  KEEP b={p[0]:+.2f} g={p[1]:.2f} w(E2)={p[2]:.2f} p(P4)={p[3]:.2f}', flush=True)
print(f'  LOSS b={p[4]:+.2f} g={p[5]:.2f} w(E2)={p[6]:.2f} p(P4)={p[7]:.2f}   A_T={p[8]:.2f} C',
      flush=True)
print(f"\n  {'animal':10s}{'amp':>7}{'obs':>7}{'corr':>8}{'RMS':>8}{'daySD':>8}"
      f"{'vs spread':>12}   peaks / observed", flush=True)
obsp = ['8.2/15.2', '15.8/21.8', '13.2 only']
for i, l in enumerate(L):
    a = W[i].max() - W[i].min()
    ext = np.concatenate([W[i]] * 3)
    pk, _ = find_peaks(ext, prominence=0.20 * a)
    m = (pk >= 48) & (pk < 96)
    cc = np.corrcoef(W[i], OBS[i])[0, 1]
    rms = float(np.sqrt(np.mean((W[i] - OBS[i]) ** 2)))
    print(f'  {l:10s}{a:7.2f}{AMP[i]:7.2f}{cc:+8.3f}{rms:8.3f}{SD[i]:8.3f}'
          f'{rms/SD[i]:11.1f}x   {np.round(hrs[pk[m]-48],1).tolist()} / {obsp[i]}', flush=True)
np.save('best_revert_baseline.npy', p)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.0), sharex=True)
for i, l in enumerate(L):
    t, y = s[l]
    for d in range(13, 18):
        m = (t >= d) & (t < d + 1)
        if m.sum() >= 20:
            ax[i].plot((t[m] - d) * 24., y[m], '-', color='0.72', lw=.9)
    cc = np.corrcoef(W[i], OBS[i])[0, 1]
    ax[i].plot(hrs, OBS[i], 'k-', lw=2.2, label='mouse 5-day mean')
    ax[i].plot(hrs, W[i], 'r--', lw=2, label="Holly's original form")
    ax[i].set_title(f'{l}   amp {a if False else (W[i].max()-W[i].min()):.2f}/{AMP[i]:.2f} C'
                    f'   r={cc:+.2f}')
    ax[i].set_xlabel('hour of day'); ax[i].set_xticks([0, 6, 12, 18, 24])
ax[0].set_ylabel('CBT (C)'); ax[0].legend(fontsize=8)
fig.suptitle('BASELINE: original noon-locked two-arm model (grey = individual days)', y=1.03)
fig.tight_layout(); fig.savefig('fig_revert_baseline.png', dpi=130, bbox_inches='tight')
print('\nwrote fig_revert_baseline.png', flush=True)
