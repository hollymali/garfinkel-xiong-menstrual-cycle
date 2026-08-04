"""Two-phase clock, fitted against a MARGIN BAND instead of a single prominence cutoff.

Why this rewrite: _twophase_fit2.py's answer put pregnant's dip at exactly 10.0% of amplitude
and Esr1i's at exactly 10.0% -- both pinned to the PROM=0.10 classifier cutoff.  Esr1i's
"unimodality" would flip to bimodal at a 9% threshold.  That is the same boundary-parking that
already invalidated `_p4keep_strict.py` (peaks 0.034 h apart) and `_fit_p4_keep3.py` (monotone
curve scored as 2 peaks).  A single cutoff is always gameable: the optimizer gets the label for
free by sitting on it.

Fix = a BAND with a gap in the middle, so there is no boundary to park on:
    bimodal   requires dip >= BI_MIN (0.20) of amplitude
    unimodal  requires dip <= UNI_MAX (0.05) of amplitude   (or a single local max outright)
Nothing can satisfy both, and the forbidden zone between them is wide.

Second fix: seeds come from the LARGE non-preg-only pool (~16.8k per 6M) and the remaining
targets are reached by CONTINUOUS penalties, instead of hard-rejecting seeds that start
outside the target (which left the previous run optimising from only 22 usable seeds).
Strict band compliance is then re-checked on the final answer and reported explicitly.
"""
import sys
import numpy as np
from scipy.optimize import minimize

fold = np.load('fold.npy', allow_pickle=True).item()
L = ['non-preg', 'pregnant', 'Esr1i']
hrs = fold[L[0]][0]
OBS = np.vstack([fold[l][1] for l in L])
AMP = OBS.max(1) - OBS.min(1)
R_OBS = AMP / AMP[0]

tc = np.arange(0, 24, 0.5)
tf = np.arange(0, 24, 0.2)
E2 = np.array([0., 1., 0.])
mid = np.array([13, 14, 15, 16, 17], float) + 0.5
p4bar = np.maximum(0., 1. - (mid - 13.) / 6.5).mean()
P4v = np.array([0., p4bar, p4bar])

BI_MIN, UNI_MAX = 0.20, 0.05
T_NP_GAP, T_NP_MID, T_PR_MID = 7.0, 11.7, 18.8
TOL = 2.0


def pr(*a):
    print(*a, flush=True)


def sig(z):
    return 1 / (1 + np.exp(-np.clip(z, -60, 60)))


def arms(p, e2, p4, t):
    bL, gL, phiL, wL, pL, bK, gK, phiK, wK, pK = [p[..., i] for i in range(10)]
    cK = np.cos(2 * np.pi * (t - phiK[..., None]) / 24.)
    cL = np.cos(2 * np.pi * (t - phiL[..., None]) / 24.)
    return (sig(bK[..., None] + gK[..., None] * cK + wK[..., None] * e2 + pK[..., None] * p4)
            - sig(bL[..., None] + gL[..., None] * cL + wL[..., None] * e2 + pL[..., None] * p4))


def metrics(W, t):
    """-> nmax, dip_frac, gap, mid, amp   (dip_frac=0 when there is a single max)"""
    amp = W.max() - W.min()
    lo = np.roll(W, 1); hi = np.roll(W, -1)
    ismax = (W > lo) & (W > hi); ismin = (W < lo) & (W < hi)
    nmax = int(ismax.sum())
    if nmax != 2 or int(ismin.sum()) != 2 or amp < 1e-12:
        return nmax, 0.0, 0.0, 0.0, amp
    pv = W[ismax]; vv = W[ismin]
    dip = (pv.min() - vv.max()) / amp
    t1, t2 = t[ismax]
    g_lin = t2 - t1
    m = (t1 + t2) / 2. if g_lin <= 12. else ((t1 + t2) / 2. + 12.) % 24.
    return nmax, float(dip), float(min(g_lin, 24 - g_lin)), float(m), amp


def relu(x):
    return x if x > 0 else 0.0


def np_screen(P, t):
    """batch: non-preg 2 peaks, dip>=BI_MIN, gap and midpoint near target"""
    W = arms(P, E2[0], P4v[0], t)
    amp = W.max(1) - W.min(1)
    lo = np.roll(W, 1, 1); hi = np.roll(W, -1, 1)
    ismax = (W > lo) & (W > hi); ismin = (W < lo) & (W < hi)
    ok = (ismax.sum(1) == 2) & (ismin.sum(1) == 2) & (amp > 1e-9)
    out = np.zeros(len(W), bool)
    idx = np.flatnonzero(ok)
    if not idx.size:
        return out
    _, cm = np.nonzero(ismax[idx]); pk = cm.reshape(-1, 2)
    _, cn = np.nonzero(ismin[idx]); vl = cn.reshape(-1, 2)
    Wt = W[idx]
    pv = np.take_along_axis(Wt, pk, 1); vv = np.take_along_axis(Wt, vl, 1)
    dip = (pv.min(1) - vv.max(1)) / amp[idx]
    t1, t2 = t[pk[:, 0]], t[pk[:, 1]]
    g_lin = t2 - t1
    m_lin = (t1 + t2) / 2.
    m = np.where(g_lin <= 12., m_lin, (m_lin + 12.) % 24.)
    g = np.minimum(g_lin, 24 - g_lin)
    out[idx] = (dip >= BI_MIN) & (np.abs(g - T_NP_GAP) < TOL) & (np.abs(m - T_NP_MID) < TOL)
    return out


# ------------------------------------------------------------------ stage-1 seeds
rng = np.random.default_rng(1618)
N, CH = 6_000_000, 50_000
seeds = []
n1 = 0
for s in range(0, N, CH):
    P = np.stack([
        rng.uniform(-15, 8, CH), rng.uniform(0, 30, CH), rng.uniform(0, 24, CH),
        rng.uniform(0, 30, CH), rng.uniform(0, 15, CH),
        rng.uniform(-15, 8, CH), rng.uniform(0, 30, CH), rng.uniform(0, 24, CH),
        rng.uniform(0, 30, CH), rng.uniform(0, 15, CH)], 1)
    m = np_screen(P, tc)
    if m.any():
        Q = P[m]
        m2 = np_screen(Q, tf)
        if m2.any():
            n1 += int(m2.sum()); seeds.append(Q[m2])
    if (s // CH) % 30 == 0:
        pr(f'  ... {s + CH:,}/{N:,}  non-preg-band seeds {n1:,}')
seeds = np.concatenate(seeds) if seeds else np.zeros((0, 10))
pr(f'draws {N:,}  non-preg BAND seeds (dip>={BI_MIN}, gap/midpoint on target): {n1:,}')
if not len(seeds):
    pr('none'); sys.exit()
np.save('_twophase_seeds3.npy', seeds)


# ------------------------------------------------------------------ objective
def model_fold(p):
    return np.stack([arms(p[None, :], E2[i], P4v[i], hrs)[0] for i in range(3)])


def solve_lin(H):
    X = np.stack([np.ones(H.size), H.ravel()], 1)
    return np.linalg.lstsq(X, OBS.ravel(), rcond=None)[0]


def shape_pen(p):
    """continuous penalty pushing all three animals into the BAND"""
    m = [metrics(arms(p[None, :], E2[i], P4v[i], tf)[0], tf) for i in range(3)]
    (n0, d0, g0, m0, _), (n1_, d1, g1, m1, _), (n2, d2, _, _, _) = m
    pen = 0.0
    pen += 4.0 * relu(BI_MIN - d0) if n0 == 2 else 4.0 * BI_MIN + 1.0
    pen += 2.0 * relu(abs(g0 - T_NP_GAP) - TOL) + 2.0 * relu(abs(m0 - T_NP_MID) - TOL)
    pen += 4.0 * relu(BI_MIN - d1) if n1_ == 2 else 4.0 * BI_MIN + 1.0
    pen += 2.0 * relu(abs(m1 - T_PR_MID) - TOL) if n1_ == 2 else 2.0
    pen += 4.0 * relu(d2 - UNI_MAX)          # Esr1i: dip must be SMALL (n2==1 -> d2=0 -> free)
    return pen, m


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
    pen, _ = shape_pen(p)
    return rms + 6.0 * ratio + 3.0 * absamp + 10.0 * pen


vals = np.array([score(s) for s in seeds])
alive = np.flatnonzero(vals < 1e5)
pr(f'  finite seeds {len(alive):,}; best {vals[alive].min():.4f}')
best = (float(vals[alive].min()), seeds[alive[np.argmin(vals[alive])]].copy())
order = alive[np.argsort(vals[alive])][:250]
for k, i in enumerate(order):
    z = seeds[i].copy()
    for _ in range(3):
        r = minimize(score, z, method='Nelder-Mead',
                     options=dict(maxiter=5000, maxfev=5000, xatol=1e-9, fatol=1e-12))
        z = r.x
    if float(r.fun) < best[0]:
        best = (float(r.fun), z.copy())
        np.save('_twophase3_ckpt.npy', best[1])
        pr(f'    [{k+1}/{len(order)}] new best {best[0]:.4f}')

v, p = best
H = model_fold(p); T0, A_T = solve_lin(H)
W = T0 + A_T * H
a = W.max(1) - W.min(1)
pen, mm = shape_pen(p)
bL, gL, phiL, wL, pL, bK, gK, phiK, wK, pK = p
d = abs(phiK - phiL) % 24
pr(f'\n=== two-phase fit, MARGIN BAND (bimodal dip>={BI_MIN}, unimodal dip<={UNI_MAX})'
   f'   objective {v:.4f}  shape-penalty {pen:.4f} ===')
pr(f'  LOSS  b={bL:+.2f} g={gL:.2f} phi={phiL:.2f} h  w(E2)={wL:.2f} p(P4)={pL:.2f}')
pr(f'  KEEP  b={bK:+.2f} g={gK:.2f} phi={phiK:.2f} h  w(E2)={wK:.2f} p(P4)={pK:.2f}')
pr(f'  ARM PHASE GAP {min(d, 24-d):.2f} h     A_T={A_T:.2f} C   T0={T0:.2f} C')
obs_pk = ['8.2/15.2', '15.8/21.8', '13.2 only']
pr(f'\n  {"animal":10s}{"amp":>7}{"obs":>7}{"ratio":>8}{"obsR":>7}{"nmax":>6}{"dip%":>7}   verdict')
band_ok = True
for i, l in enumerate(L):
    n, dp, g, m, _ = mm[i]
    if i < 2:
        good = (n == 2 and dp >= BI_MIN)
        vd = 'BIMODAL ok' if good else f'FAILS band (needs dip>={BI_MIN})'
    else:
        good = (n == 1) or (dp <= UNI_MAX)
        vd = 'unimodal ok' if good else f'FAILS band (needs dip<={UNI_MAX})'
    band_ok &= good
    pr(f'  {l:10s}{a[i]:7.2f}{AMP[i]:7.2f}{a[i]/a[0]:8.2f}{R_OBS[i]:7.2f}{n:6d}{dp*100:7.1f}   {vd}')
pr(f'  observed peaks: ' + ' | '.join(f'{l} {o}' for l, o in zip(L, obs_pk)))
pr(f'\n  STRICT BAND SATISFIED: {band_ok}')
np.save('best_twophase3.npy', np.concatenate([p, [A_T, T0]]))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.0), sharex=True)
for i, l in enumerate(L):
    ax[i].plot(hrs, OBS[i], 'k-', lw=2.2, label='mouse')
    ax[i].plot(hrs, W[i], 'r--', lw=2, label='two-phase clock')
    ax[i].axvline(12, color='0.85', lw=1, zorder=0)
    ax[i].set_title(f'{l}   amp {a[i]:.2f} / {AMP[i]:.2f} C')
    ax[i].set_xlabel('hour of day'); ax[i].set_xticks([0, 6, 12, 18, 24])
ax[0].set_ylabel('CBT (C)'); ax[0].legend(fontsize=8)
fig.suptitle(f'Two-phase clock, margin band (bimodal dip>={BI_MIN}, unimodal dip<={UNI_MAX}); '
             f'arm phase gap {min(d,24-d):.2f} h', y=1.03)
fig.tight_layout(); fig.savefig('fig_twophase3.png', dpi=130, bbox_inches='tight')
pr('wrote fig_twophase3.png, best_twophase3.npy')
