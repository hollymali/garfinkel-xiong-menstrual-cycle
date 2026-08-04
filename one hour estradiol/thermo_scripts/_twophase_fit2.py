"""Two-phase clock fit -- staged/flushing rewrite of _twophase_fit.py.

Same model and same (corrected) target; three practical fixes after the first version was
killed mid-run and lost everything to output buffering:
  * STAGED SCREEN: non-preg's topology alone rejects ~99.7% of draws, so its waveform is
    tested FIRST on a coarse 0.5 h grid and pregnant/Esr1i are only computed for survivors.
  * seeds are SAVED TO DISK as soon as the scan ends, before any polishing, so a kill during
    Nelder-Mead cannot cost the scan.
  * every print flushes, and the best-so-far is checkpointed during polishing.

Model:
    keep = sig(bK + gK*cos(2pi(t-phiK)/24) + wK*E2 + pK*P4)
    loss = sig(bL + gL*cos(2pi(t-phiL)/24) + wL*E2 + pL*P4)
    CBT  = T0 + A_T*(keep - loss)
E2 = 0/1/0, P4 = 0/p4bar/p4bar for non-preg / pregnant / Esr1i.  All coeffs >= 0.

TARGET (all four, the last of which noon-lock forbids outright):
    non-preg 2 peaks, gap 7.0+-2 h, MIDPOINT 11.7+-2 h   <- midpoint was missing before
    Esr1i    1 peak
    pregnant 2 peaks, midpoint 18.8+-2 h
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

tc = np.arange(0, 24, 0.5)     # coarse screen
tf = np.arange(0, 24, 0.2)     # fine verify
E2 = np.array([0., 1., 0.])
mid = np.array([13, 14, 15, 16, 17], float) + 0.5
p4bar = np.maximum(0., 1. - (mid - 13.) / 6.5).mean()
P4v = np.array([0., p4bar, p4bar])

PROM = 0.10
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


def peak_stats(W, t):
    amp = W.max(1) - W.min(1)
    lo = np.roll(W, 1, 1); hi = np.roll(W, -1, 1)
    ismax = (W > lo) & (W > hi)
    ismin = (W < lo) & (W < hi)
    nmax = ismax.sum(1)
    gap = np.zeros(len(W)); midp = np.zeros(len(W)); nsig = nmax.copy()
    two = np.flatnonzero((nmax == 2) & (ismin.sum(1) == 2) & (amp > 1e-9))
    if two.size:
        _, cm = np.nonzero(ismax[two]); pk = cm.reshape(-1, 2)
        _, cn = np.nonzero(ismin[two]); vl = cn.reshape(-1, 2)
        Wt = W[two]
        pv = np.take_along_axis(Wt, pk, 1)
        vv = np.take_along_axis(Wt, vl, 1)
        real2 = (pv.min(1) - vv.max(1)) > PROM * amp[two]
        t1, t2 = t[pk[:, 0]], t[pk[:, 1]]
        g_lin = t2 - t1
        m_lin = (t1 + t2) / 2.
        midp[two] = np.where(g_lin <= 12., m_lin, (m_lin + 12.) % 24.)
        gap[two] = np.minimum(g_lin, 24 - g_lin)
        nsig[two] = np.where(real2, 2, 1)
    return nsig, gap, midp, amp


def np_ok(P, t):
    n, g, m, _ = peak_stats(arms(P, E2[0], P4v[0], t), t)
    return (n == 2) & (np.abs(g - T_NP_GAP) < TOL) & (np.abs(m - T_NP_MID) < TOL)


def full_ok(P, t):
    if not len(P):
        return np.zeros(0, bool)
    ok = np_ok(P, t)
    if not ok.any():
        return ok
    n_es, _, _, _ = peak_stats(arms(P, E2[2], P4v[2], t), t)
    n_pr, _, m_pr, _ = peak_stats(arms(P, E2[1], P4v[1], t), t)
    return ok & (n_es == 1) & (n_pr == 2) & (np.abs(m_pr - T_PR_MID) < TOL)


# ------------------------------------------------------------------ staged seed scan
rng = np.random.default_rng(4242)
N, CH = 6_000_000, 50_000
seeds = []
n_np, n_full = 0, 0
for s in range(0, N, CH):
    P = np.stack([
        rng.uniform(-15, 8, CH), rng.uniform(0, 30, CH), rng.uniform(0, 24, CH),
        rng.uniform(0, 30, CH), rng.uniform(0, 15, CH),
        rng.uniform(-15, 8, CH), rng.uniform(0, 30, CH), rng.uniform(0, 24, CH),
        rng.uniform(0, 30, CH), rng.uniform(0, 15, CH)], 1)
    m1 = np_ok(P, tc)                      # cheap: non-preg only, coarse grid
    n_np += int(m1.sum())
    if m1.any():
        Q = P[m1]
        m2 = full_ok(Q, tf)                # survivors only, fine grid
        if m2.any():
            n_full += int(m2.sum())
            seeds.append(Q[m2])
    if (s // CH) % 20 == 0:
        pr(f'  ... {s + CH:,}/{N:,} draws  non-preg-ok {n_np:,}  full-topology {n_full:,}')
seeds = np.concatenate(seeds) if seeds else np.zeros((0, 10))
pr(f'draws {N:,}')
pr(f'  non-preg 2 peaks + gap {T_NP_GAP}+-{TOL} + midpoint {T_NP_MID}+-{TOL} : {n_np:,}')
pr(f'  FULL topology (all 3 animals, incl. pregnant midpoint {T_PR_MID}) : {n_full:,}')
if len(seeds):
    np.save('_twophase_seeds2.npy', seeds)
    pr(f'  seeds saved -> _twophase_seeds2.npy')
else:
    pr('  NO full-topology seeds.')
    sys.exit()


# ------------------------------------------------------------------ objective
def model_fold(p):
    return np.stack([arms(p[None, :], E2[i], P4v[i], hrs)[0] for i in range(3)])


def solve_lin(H):
    X = np.stack([np.ones(H.size), H.ravel()], 1)
    return np.linalg.lstsq(X, OBS.ravel(), rcond=None)[0]


def score(p):
    if np.any(p[[1, 3, 4, 6, 8, 9]] < 0):
        return 1e6
    if not bool(full_ok(p[None, :], tf)[0]):
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
    return rms + 6.0 * ratio + 3.0 * absamp


vals = np.array([score(s) for s in seeds])
alive = np.flatnonzero(vals < 1e5)
pr(f'  finite-scoring seeds {len(alive):,}' + (f'; best {vals[alive].min():.4f}' if len(alive) else ''))
if not len(alive):
    sys.exit()
best = (float(vals[alive].min()), seeds[alive[np.argmin(vals[alive])]].copy())
order = alive[np.argsort(vals[alive])][:200]
for k, i in enumerate(order):
    z = seeds[i].copy()
    for _ in range(3):
        r = minimize(score, z, method='Nelder-Mead',
                     options=dict(maxiter=5000, maxfev=5000, xatol=1e-9, fatol=1e-12))
        z = r.x
    if float(r.fun) < best[0]:
        best = (float(r.fun), z.copy())
        np.save('_twophase_best_ckpt.npy', best[1])
        pr(f'    [{k+1}/{len(order)}] new best {best[0]:.4f}')

v, p = best
H = model_fold(p)
T0, A_T = solve_lin(H)
W = T0 + A_T * H
a = W.max(1) - W.min(1)
bL, gL, phiL, wL, pL, bK, gK, phiK, wK, pK = p
d = abs(phiK - phiL) % 24
pr(f'\n=== best two-phase fit   objective {v:.4f} ===')
pr(f'  LOSS  b={bL:+.2f} g={gL:.2f} phi={phiL:.2f} h  w(E2)={wL:.2f} p(P4)={pL:.2f}')
pr(f'  KEEP  b={bK:+.2f} g={gK:.2f} phi={phiK:.2f} h  w(E2)={wK:.2f} p(P4)={pK:.2f}')
pr(f'  ARM PHASE GAP {min(d, 24-d):.2f} h     A_T={A_T:.2f} C   T0={T0:.2f} C')
pr(f'\n  {"animal":10s}{"amp mod":>9}{"amp obs":>9}{"ratio mod":>11}{"ratio obs":>11}')
for i, l in enumerate(L):
    pr(f'  {l:10s}{a[i]:9.2f}{AMP[i]:9.2f}{a[i]/a[0]:11.2f}{R_OBS[i]:11.2f}')
pr('  observed peaks: non-preg 8.2/15.2 | pregnant 15.8/21.8 | Esr1i 13.2 only')
np.save('best_twophase.npy', np.concatenate([p, [A_T, T0]]))

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
fig.suptitle('Heat-KEEP and heat-LOSS arms driven by the SCN at DIFFERENT phases '
             f'(gap {min(d,24-d):.1f} h); E2=0 in Esr1i', y=1.03)
fig.tight_layout(); fig.savefig('fig_twophase.png', dpi=130, bbox_inches='tight')
pr('\nwrote fig_twophase.png, best_twophase.npy')
