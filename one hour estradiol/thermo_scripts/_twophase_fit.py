"""Two-phase clock, FITTED -- scan for topology-valid seeds then Nelder-Mead polish.

Fixes two problems with `_twophase_feasible.py`:
  1. it constrained non-preg's peak GAP but not its MIDPOINT, so the best-ratio draw parked
     non-preg's pair at 19.6 h instead of the observed 11.7 h.  Midpoint is now a target.
  2. it asked RANDOM DRAWS to land on the amplitude ratios (0 of 3M did).  Random draws are
     only good for seeds; the ratios have to be OPTIMISED for.

Model (Holly's arms, one change: the two arms see the clock at DIFFERENT phases, as separate
SCN morning/evening cell groups projecting differentially would give):
    keep = sig(bK + gK*cos(2pi(t-phiK)/24) + wK*E2 + pK*P4)
    loss = sig(bL + gL*cos(2pi(t-phiL)/24) + wL*E2 + pL*P4)
    CBT  = T0 + A_T*(keep - loss)
E2 = 0/1/0, P4 = 0/p4bar/p4bar for non-preg / pregnant / Esr1i.  All coeffs >= 0.
T0 and A_T enter linearly -> solved in closed form for every shape draw.

Objective follows the recorded method traps: amplitudes are weighted HEAVILY (RMS alone picks
the wrong model every time), and topology is a hard reject rather than a soft penalty.
Fold-only fit (no day-varying drift term) so there is no window-centring double-count.
"""
import numpy as np
from scipy.optimize import minimize

fold = np.load('fold.npy', allow_pickle=True).item()
L = ['non-preg', 'pregnant', 'Esr1i']
hrs = fold[L[0]][0]
OBS = np.vstack([fold[l][1] for l in L])
AMP = OBS.max(1) - OBS.min(1)
R_OBS = AMP / AMP[0]

tg = np.arange(0, 24, 0.2)
E2 = np.array([0., 1., 0.])
mid = np.array([13, 14, 15, 16, 17], float) + 0.5
p4bar = np.maximum(0., 1. - (mid - 13.) / 6.5).mean()
P4v = np.array([0., p4bar, p4bar])

PROM = 0.10
T_NP_GAP, T_NP_MID = 7.0, 11.7
T_PR_MID = 18.8
TOL_GAP, TOL_MID = 2.0, 2.0


def sig(z):
    return 1 / (1 + np.exp(-np.clip(z, -60, 60)))


def arms(p, e2, p4, t):
    bL, gL, phiL, wL, pL, bK, gK, phiK, wK, pK = [p[..., i] for i in range(10)]
    cK = np.cos(2 * np.pi * (t - phiK[..., None]) / 24.)
    cL = np.cos(2 * np.pi * (t - phiL[..., None]) / 24.)
    return (sig(bK[..., None] + gK[..., None] * cK + wK[..., None] * e2 + pK[..., None] * p4)
            - sig(bL[..., None] + gL[..., None] * cL + wL[..., None] * e2 + pL[..., None] * p4))


def peak_stats(W):
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
        depth = pv.min(1) - vv.max(1)
        real2 = depth > PROM * amp[two]
        t1, t2 = tg[pk[:, 0]], tg[pk[:, 1]]
        g_lin = t2 - t1
        m_lin = (t1 + t2) / 2.
        midp[two] = np.where(g_lin <= 12., m_lin, (m_lin + 12.) % 24.)
        gap[two] = np.minimum(g_lin, 24 - g_lin)
        nsig[two] = np.where(real2, 2, 1)
    return nsig, gap, midp, amp


def topo_ok_batch(P):
    st = [peak_stats(arms(P, E2[i], P4v[i], tg)) for i in range(3)]
    n_np, g_np, m_np, _ = st[0]
    n_pr, _, m_pr, _ = st[1]
    n_es, _, _, _ = st[2]
    return ((n_np == 2) & (np.abs(g_np - T_NP_GAP) < TOL_GAP)
            & (np.abs(m_np - T_NP_MID) < TOL_MID)
            & (n_es == 1) & (n_pr == 2)
            & (np.abs(m_pr - T_PR_MID) < TOL_MID))


# ------------------------------------------------------------------ seed scan
rng = np.random.default_rng(4242)
N, CH = 4_000_000, 20_000
seeds = []
n_hit = 0
for s in range(0, N, CH):
    P = np.stack([
        rng.uniform(-15, 8, CH), rng.uniform(0, 30, CH), rng.uniform(0, 24, CH),
        rng.uniform(0, 30, CH), rng.uniform(0, 15, CH),
        rng.uniform(-15, 8, CH), rng.uniform(0, 30, CH), rng.uniform(0, 24, CH),
        rng.uniform(0, 30, CH), rng.uniform(0, 15, CH)], 1)
    ok = topo_ok_batch(P)
    if ok.any():
        n_hit += int(ok.sum())
        seeds.append(P[ok])
seeds = np.concatenate(seeds) if seeds else np.zeros((0, 10))
print(f'draws {N:,}   full-topology seeds (incl. non-preg MIDPOINT): {n_hit:,}')
if not len(seeds):
    raise SystemExit('no seeds')


# ------------------------------------------------------------------ objective
def solve_lin(H):
    X = np.stack([np.ones(H.size), H.ravel()], 1)
    coef, *_ = np.linalg.lstsq(X, OBS.ravel(), rcond=None)
    return coef


def model_fold(p):
    return np.stack([arms(p[None, :], E2[i], P4v[i], hrs)[0] for i in range(3)])


def score(p):
    if np.any(p[[1, 3, 4, 6, 8, 9]] < 0):
        return 1e6
    if not bool(topo_ok_batch(p[None, :])[0]):
        return 1e6
    H = model_fold(p)
    T0, A_T = solve_lin(H)
    if A_T <= 0 or A_T > 10.0:            # physical bound on the Celsius scale factor
        return 1e6
    W = T0 + A_T * H
    a = W.max(1) - W.min(1)
    if a[0] < 1e-6:
        return 1e6
    r = a / a[0]
    rms = float(np.mean(np.sqrt(np.mean((W - OBS) ** 2, 1)) / AMP))
    ratio = float(np.abs(r[1] - R_OBS[1]) + np.abs(r[2] - R_OBS[2]))
    absamp = float(abs(a[0] - AMP[0]) / AMP[0])
    return rms + 6.0 * ratio + 3.0 * absamp


vals = np.array([score(s) for s in seeds])
alive = np.flatnonzero(vals < 1e5)
print(f'  seeds scoring finite: {len(alive):,}; best seed {vals[alive].min():.4f}'
      if len(alive) else '  no finite seed')
if not len(alive):
    raise SystemExit
best = (float(vals[alive].min()), seeds[alive[np.argmin(vals[alive])]].copy())
for i in alive[np.argsort(vals[alive])[:300]]:
    z = seeds[i].copy()
    for _ in range(4):
        r = minimize(score, z, method='Nelder-Mead',
                     options=dict(maxiter=6000, maxfev=6000, xatol=1e-9, fatol=1e-12))
        z = r.x
    if float(r.fun) < best[0]:
        best = (float(r.fun), z.copy())

v, p = best
H = model_fold(p)
T0, A_T = solve_lin(H)
W = T0 + A_T * H
a = W.max(1) - W.min(1)
bL, gL, phiL, wL, pL, bK, gK, phiK, wK, pK = p
d = abs(phiK - phiL) % 24
print(f'\n=== best two-phase fit   objective {v:.4f} ===')
print(f'  LOSS  b={bL:+.2f} g={gL:.2f} phi={phiL:.2f} h  w(E2)={wL:.2f} p(P4)={pL:.2f}')
print(f'  KEEP  b={bK:+.2f} g={gK:.2f} phi={phiK:.2f} h  w(E2)={wK:.2f} p(P4)={pK:.2f}')
print(f'  ARM PHASE GAP {min(d, 24-d):.2f} h        A_T={A_T:.2f} C   T0={T0:.2f} C')
print(f'\n  {"animal":10s}{"amp mod":>9}{"amp obs":>9}{"ratio mod":>11}{"ratio obs":>11}{"peaks (h)":>22}')
obs_pk = ['8.2, 15.2', '15.8, 21.8', '13.2']
for i, l in enumerate(L):
    n, g, m, _ = peak_stats(arms(p[None, :], E2[i], P4v[i], tg))
    Wf = arms(p[None, :], E2[i], P4v[i], tg)[0]
    lo = np.roll(Wf, 1); hi = np.roll(Wf, -1)
    pk = tg[(Wf > lo) & (Wf > hi)]
    pks = ', '.join(f'{x:.1f}' for x in pk)
    print(f'  {l:10s}{a[i]:9.2f}{AMP[i]:9.2f}{a[i]/a[0]:11.2f}{R_OBS[i]:11.2f}{pks:>22}')
    print(f'  {"":10s}{"":9s}{"":9s}{"":11s}{"obs:":>11}{obs_pk[i]:>22}')

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
print('\nwrote fig_twophase.png, best_twophase.npy')
