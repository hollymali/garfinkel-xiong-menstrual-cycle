"""Two-phase clock, fitted on PEAK POSITIONS with an exact prominence counter.

Four earlier topology tests were gamed by the optimizer, every time in the model's favour:
  _p4keep_strict.py    peaks 0.034 h apart scored as a "pair"
  _fit_p4_keep3.py     a monotone rise-then-fall scored as 2 peaks (trough tie)
  _twophase_fit2.py    pregnant and Esr1i dips parked at exactly the 10.0% cutoff
  _twophase_fit3.py    dip metric returned 0.0 whenever nmax != 2, so a THREE-peak Esr1i
                       passed the "unimodal" test
Common cause: a hand-rolled classifier with a single cutoff, never validated on the real data.

This file fixes the method, not just the symptom:
  * ONE peak routine, `count_peaks`, using scipy prominence on a tripled (circular) array.
  * VALIDATED ON THE OBSERVED DATA FIRST.  It reproduces the target 2/2/1 at thr=0.20 and only
    near there (0.15 -> 3/3/1, 0.25 -> 1/2/1).  So 2/2/1 is a narrow-window description of the
    fold, and demanding threshold-robustness the DATA does not have is what made every previous
    test gameable.
  * The objective therefore targets PEAK POSITIONS (observed 8.2/15.2, 15.8/21.8, 13.2), which
    a boundary-parked solution cannot satisfy, plus amplitudes.  Count mismatch is penalised but
    position is what carries the constraint.

Model: keep and loss arms driven by the SCN at DIFFERENT phases (separate M/E cell groups
projecting differentially); one extra parameter vs the shared noon-locked clock.
    keep = sig(bK + gK*cos(2pi(t-phiK)/24) + wK*E2 + pK*P4)
    loss = sig(bL + gL*cos(2pi(t-phiL)/24) + wL*E2 + pL*P4)
    CBT  = T0 + A_T*(keep - loss)
"""
import sys
import numpy as np
from scipy.signal import find_peaks
from scipy.optimize import minimize

fold = np.load('fold.npy', allow_pickle=True).item()
L = ['non-preg', 'pregnant', 'Esr1i']
hrs = fold[L[0]][0]
OBS = np.vstack([fold[l][1] for l in L])
AMP = OBS.max(1) - OBS.min(1)
R_OBS = AMP / AMP[0]

tf = np.arange(0, 24, 0.2)
E2 = np.array([0., 1., 0.])
mid = np.array([13, 14, 15, 16, 17], float) + 0.5
p4bar = np.maximum(0., 1. - (mid - 13.) / 6.5).mean()
P4v = np.array([0., p4bar, p4bar])

THR = 0.20                      # the threshold at which the DATA shows 2/2/1
OBS_PEAKS = [np.array([8.2, 15.2]), np.array([15.8, 21.8]), np.array([13.2])]


def pr(*a):
    print(*a, flush=True)


def sig(z):
    return 1 / (1 + np.exp(-np.clip(z, -60, 60)))


def count_peaks(W, t, thr=THR):
    W = np.asarray(W, float); n = len(W)
    amp = W.max() - W.min()
    if amp < 1e-12:
        return np.array([]), np.array([]), amp
    ext = np.concatenate([W, W, W])
    pk, props = find_peaks(ext, prominence=thr * amp)
    m = (pk >= n) & (pk < 2 * n)
    return t[pk[m] - n], props['prominences'][m], amp


def arms(p, e2, p4, t):
    bL, gL, phiL, wL, pL, bK, gK, phiK, wK, pK = [p[..., i] for i in range(10)]
    cK = np.cos(2 * np.pi * (t - phiK[..., None]) / 24.)
    cL = np.cos(2 * np.pi * (t - phiL[..., None]) / 24.)
    return (sig(bK[..., None] + gK[..., None] * cK + wK[..., None] * e2 + pK[..., None] * p4)
            - sig(bL[..., None] + gL[..., None] * cL + wL[..., None] * e2 + pL[..., None] * p4))


def circ_d(a, b):
    d = np.abs(a - b) % 24.
    return np.minimum(d, 24. - d)


def peak_penalty(p):
    """distance between model and observed peak SETS, per animal"""
    tot = 0.0
    info = []
    for i in range(3):
        W = arms(p[None, :], E2[i], P4v[i], tf)[0]
        pk, prom, amp = count_peaks(W, tf)
        obs = OBS_PEAKS[i]
        if len(pk) == 0:
            tot += 10.0; info.append((pk, prom)); continue
        # every observed peak must have a model peak near it
        d1 = np.mean([circ_d(o, pk).min() for o in obs])
        # and every model peak must be near an observed one (punishes extra peaks)
        d2 = np.mean([circ_d(m_, obs).min() for m_ in pk])
        tot += (d1 + d2) / 2.0 / 3.0          # hours -> ~O(1)
        tot += 1.5 * abs(len(pk) - len(obs))  # count mismatch
        info.append((pk, prom))
    return tot, info


def model_fold(p):
    return np.stack([arms(p[None, :], E2[i], P4v[i], hrs)[0] for i in range(3)])


def solve_lin(H):
    X = np.stack([np.ones(H.size), H.ravel()], 1)
    return np.linalg.lstsq(X, OBS.ravel(), rcond=None)[0]


def score(p):
    if np.any(p[[1, 3, 4, 6, 8, 9]] < 0):
        return 1e6
    H = model_fold(p)
    T0, A_T = solve_lin(H)
    if A_T <= 0 or A_T > 10.0:
        return 1e6
    W = T0 + A_T * H
    a = W.max(1) - W.min(1)
    if a[0] < 1e-6:
        return 1e6
    r = a / a[0]
    rms = float(np.mean(np.sqrt(np.mean((W - OBS) ** 2, 1)) / AMP))
    ratio = float(abs(r[1] - R_OBS[1]) + abs(r[2] - R_OBS[2]))
    absamp = float(abs(a[0] - AMP[0]) / AMP[0])
    pen, _ = peak_penalty(p)
    return rms + 4.0 * ratio + 3.0 * absamp + 3.0 * pen


# --------------------------------------------------- cheap vectorised seed screen
def screen(P, t):
    """non-preg must have 2 local maxima with a real dip and roughly the right placement"""
    W = arms(P, E2[0], P4v[0], t)
    amp = W.max(1) - W.min(1)
    lo = np.roll(W, 1, 1); hi = np.roll(W, -1, 1)
    imx = (W > lo) & (W > hi); imn = (W < lo) & (W < hi)
    ok = (imx.sum(1) == 2) & (imn.sum(1) == 2) & (amp > 1e-9)
    out = np.zeros(len(W), bool)
    idx = np.flatnonzero(ok)
    if not idx.size:
        return out
    _, cm = np.nonzero(imx[idx]); pk = cm.reshape(-1, 2)
    _, cn = np.nonzero(imn[idx]); vl = cn.reshape(-1, 2)
    Wt = W[idx]
    pv = np.take_along_axis(Wt, pk, 1); vv = np.take_along_axis(Wt, vl, 1)
    dip = (pv.min(1) - vv.max(1)) / amp[idx]
    t1, t2 = t[pk[:, 0]], t[pk[:, 1]]
    out[idx] = (dip >= 0.20) & (np.abs(t1 - 8.2) < 3.0) & (np.abs(t2 - 15.2) < 3.0)
    return out


rng = np.random.default_rng(271828)
N, CH = 4_000_000, 50_000
seeds = []
n1 = 0
for s in range(0, N, CH):
    P = np.stack([
        rng.uniform(-15, 8, CH), rng.uniform(0, 30, CH), rng.uniform(0, 24, CH),
        rng.uniform(0, 30, CH), rng.uniform(0, 15, CH),
        rng.uniform(-15, 8, CH), rng.uniform(0, 30, CH), rng.uniform(0, 24, CH),
        rng.uniform(0, 30, CH), rng.uniform(0, 15, CH)], 1)
    m = screen(P, tf)
    if m.any():
        n1 += int(m.sum()); seeds.append(P[m])
    if (s // CH) % 30 == 0:
        pr(f'  ... {s + CH:,}/{N:,}  seeds {n1:,}')
seeds = np.concatenate(seeds) if seeds else np.zeros((0, 10))
pr(f'draws {N:,}  seeds (non-preg 2 peaks, dip>=0.20, near 8.2/15.2): {n1:,}')
if not len(seeds):
    pr('none'); sys.exit()

vals = np.array([score(s) for s in seeds[:6000]])
alive = np.flatnonzero(vals < 1e5)
pr(f'  scored {len(vals):,} seeds, finite {len(alive):,}, best {vals[alive].min():.4f}')
best = (float(vals[alive].min()), seeds[alive[np.argmin(vals[alive])]].copy())
for k, i in enumerate(alive[np.argsort(vals[alive])][:100]):
    z = seeds[i].copy()
    for _ in range(3):
        r = minimize(score, z, method='Nelder-Mead',
                     options=dict(maxiter=4000, maxfev=4000, xatol=1e-9, fatol=1e-12))
        z = r.x
    if float(r.fun) < best[0]:
        best = (float(r.fun), z.copy())
        np.save('_twophase_final_ckpt.npy', best[1])
        pr(f'    [{k+1}/100] new best {best[0]:.4f}')

v, p = best
H = model_fold(p); T0, A_T = solve_lin(H)
W = T0 + A_T * H
a = W.max(1) - W.min(1)
pen, info = peak_penalty(p)
bL, gL, phiL, wL, pL, bK, gK, phiK, wK, pK = p
d = abs(phiK - phiL) % 24
pr(f'\n=== two-phase FINAL   objective {v:.4f}   peak-penalty {pen:.4f} ===')
pr(f'  LOSS  b={bL:+.2f} g={gL:.2f} phi={phiL:.2f} h  w(E2)={wL:.2f} p(P4)={pL:.2f}')
pr(f'  KEEP  b={bK:+.2f} g={gK:.2f} phi={phiK:.2f} h  w(E2)={wK:.2f} p(P4)={pK:.2f}')
pr(f'  ARM PHASE GAP {min(d, 24-d):.2f} h    A_T={A_T:.2f} C   T0={T0:.2f} C')
pr(f'\n  {"animal":10s}{"amp":>7}{"obs":>7}{"ratio":>8}{"obsR":>7}   model peaks (prominence%)   observed')
for i, l in enumerate(L):
    pk, prom = info[i]
    amp_i = a[i]
    s_pk = ', '.join(f'{x:.1f}({y/ (W[i].max()-W[i].min()) *100:.0f}%)'
                     for x, y in zip(pk, prom)) if len(pk) else '(none)'
    s_ob = ', '.join(f'{x:.1f}' for x in OBS_PEAKS[i])
    pr(f'  {l:10s}{amp_i:7.2f}{AMP[i]:7.2f}{amp_i/a[0]:8.2f}{R_OBS[i]:7.2f}   {s_pk:<28}{s_ob}')
np.save('best_twophase_final.npy', np.concatenate([p, [A_T, T0]]))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.0), sharex=True)
for i, l in enumerate(L):
    ax[i].plot(hrs, OBS[i], 'k-', lw=2.2, label='mouse')
    ax[i].plot(hrs, W[i], 'r--', lw=2, label='two-phase clock')
    for x in OBS_PEAKS[i]:
        ax[i].axvline(x, color='0.8', ls=':', lw=1, zorder=0)
    ax[i].set_title(f'{l}   amp {a[i]:.2f} / {AMP[i]:.2f} C')
    ax[i].set_xlabel('hour of day'); ax[i].set_xticks([0, 6, 12, 18, 24])
ax[0].set_ylabel('CBT (C)'); ax[0].legend(fontsize=8)
fig.suptitle('Two-phase clock fitted on PEAK POSITIONS (dotted = observed peaks); '
             f'arm phase gap {min(d,24-d):.2f} h', y=1.03)
fig.tight_layout(); fig.savefig('fig_twophase_final.png', dpi=130, bbox_inches='tight')
pr('wrote fig_twophase_final.png, best_twophase_final.npy')
