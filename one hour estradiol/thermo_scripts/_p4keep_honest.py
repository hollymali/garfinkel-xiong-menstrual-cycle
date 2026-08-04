"""Redo the P4-in-heat-KEEP-arm-only fit with a topology test that cannot be gamed.

Both earlier attempts "succeeded" via a loophole, not a real double hump:
  _p4keep_strict.py   (continuous sign-flip test)  -> peaks 0.034 h apart (single peak
                                                        sitting on the classifier boundary)
  _fit_p4_keep3.py     (discrete npeaks_grid test)   -> the "2nd peak" was a tie-breaking
                                                        artifact at the trough, non-preg
                                                        curve is monotone rise-then-fall

Fix: same continuous derivative test as _p4keep_strict.py, but require the resulting
pair's separation to clear an explicit floor (GAP_MIN hours) so the optimizer cannot
collapse the two peaks onto each other to cheat the amplitude loss.  Random draws show
genuine PAIR solutions are NOT rare (98% of all PAIR draws already clear 3h; median gap
~9.9h) -- so if a good amplitude fit still can't be found under this floor, that's a real
result, not a search-coverage problem.
"""
import numpy as np
import math
from scipy.optimize import minimize

fold = np.load('fold.npy', allow_pickle=True).item()
L = ['non-preg', 'pregnant', 'Esr1i']
hrs = fold[L[0]][0]
OBS = np.vstack([fold[l][1] for l in L])
AMP = OBS.max(1) - OBS.min(1)
omean = np.load('_omean.npy')
ocen = omean - omean.mean(1, keepdims=True)

c48 = np.cos(2 * np.pi * (hrs - 12.) / 24.)
cg = np.linspace(-1, 1, 801)
E2 = np.array([0., 1., 0.])
PREG = np.array([0., 1., 1.])
mid = np.array([13, 14, 15, 16, 17], float) + 0.5
P4d = np.maximum(0., 1. - (mid - 13.) / 6.5)
P4M = np.outer(PREG, P4d)
dd = mid - 15.5
GAP_MIN = 3.0   # hours; well below observed 7.0h, comfortably above the degenerate zone

PAIR, NOON, MIDN, OTHER = 0, 1, 2, 3


def sig(z):
    return 1 / (1 + np.exp(-np.clip(z, -60, 60)))


def dsig(z):
    s = sig(z)
    return s * (1 - s)


def shape_class_gap(bk, gK, bl, gL):
    """like shape_class but also returns the PAIR peak separation in hours (0 if not PAIR)."""
    uK = bk[..., None] + gK[..., None] * cg
    uL = bl[..., None] + gL[..., None] * cg
    d = gK[..., None] * dsig(uK) - gL[..., None] * dsig(uL)
    pos = d > 0
    flips = pos[..., :-1] & ~pos[..., 1:]
    drops = np.count_nonzero(flips, axis=-1)
    up_noon = pos[..., -1]
    dn_mid = ~pos[..., 0]
    cls = np.full(drops.shape, OTHER, np.int8)
    cls = np.where((drops == 1) & ~up_noon & ~dn_mid, PAIR, cls)
    cls = np.where((drops == 0) & up_noon & ~dn_mid, NOON, cls)
    cls = np.where((drops == 0) & dn_mid & ~up_noon, MIDN, cls)
    # gap: locate the (single) flip index where PAIR, convert cg[idx] -> hour separation
    idx = np.argmax(flips, axis=-1)
    cflip = cg[idx]
    ang = np.arccos(np.clip(cflip, -1, 1))
    gap_h = 2 * ang * 24 / (2 * np.pi)
    gap_h = np.where(cls == PAIR, gap_h, 0.0)
    f = sig(uK) - sig(uL)
    return cls, f.max(-1) - f.min(-1), gap_h


rng = np.random.default_rng(2718)
N, CH = 8_000_000, 100_000
p4bar = P4d.mean()
seeds = []
n_np_pair_wide = 0
n_es_noon = 0
n_both_wide = 0
for s in range(0, N, CH):
    bK = rng.uniform(-30, 10, CH); gK = rng.uniform(0, 60, CH)
    bL = rng.uniform(-30, 10, CH); gL = rng.uniform(0, 60, CH)
    pK = rng.uniform(0, 25, CH)
    wK = rng.uniform(0, 60, CH); wL = rng.uniform(0, 60, CH)
    c_np, a_np, gap_np = shape_class_gap(bK, gK, bL, gL)
    c_es, a_es, _ = shape_class_gap(bK + pK * p4bar, gK, bL, gL)
    m_np = (c_np == PAIR) & (gap_np >= GAP_MIN)
    m_es = c_es == NOON
    n_np_pair_wide += int(m_np.sum())
    n_es_noon += int(m_es.sum())
    both = m_np & m_es
    n_both_wide += int(both.sum())
    if both.any():
        seeds.append(np.stack([bL[both], gL[both], wL[both], bK[both], gK[both],
                                wK[both], pK[both]], 1))
seeds = np.concatenate(seeds) if seeds else np.zeros((0, 7))
print(f'draws {N:,}  GAP_MIN={GAP_MIN}h')
print(f'  non-preg PAIR with gap>={GAP_MIN}h    : {n_np_pair_wide:,}')
print(f'  Esr1i NOON                          : {n_es_noon:,}')
print(f'  BOTH simultaneously                 : {n_both_wide:,}')

if len(seeds) == 0:
    print('\nNO seeds found -- genuine wide-gap PAIR + Esr1i NOON may be structurally excluded.')
    raise SystemExit


def classes_gaps_of(sh):
    bL, gL, wL, bK, gK, wK, pK = sh
    bk = np.array([bK, bK + wK + pK * p4bar, bK + pK * p4bar])
    bl = np.array([bL, bL + wL, bL])
    g1 = np.full(3, gK); g2 = np.full(3, gL)
    cls, _, gap = shape_class_gap(bk, g1, bl, g2)
    return cls, gap


def H3_of(sh):
    bL, gL, wL, bK, gK, wK, pK = sh
    e2 = E2[:, None, None]; p4 = P4M[:, :, None]
    return (sig(bK + gK * c48 + wK * e2 + pK * p4) - sig(bL + gL * c48 + wL * e2))


def assemble(sh):
    H3 = H3_of(sh)
    Hf = H3.mean(0 if H3.ndim == 2 else 1)
    X = np.stack([np.ones(Hf.size), Hf.ravel()], 1)
    coef, *_ = np.linalg.lstsq(X, OBS.ravel(), rcond=None)
    T0, A_T = coef
    if A_T <= 0:
        return None
    hm = H3.mean(2); hcen = hm - hm.mean(1, keepdims=True)
    dr = float(((ocen - A_T * hcen) * dd[None, :]).sum() / (3.0 * float((dd ** 2).sum())))
    return T0 + A_T * H3 + dr * dd[None, :, None], (T0, A_T, dr)


def score(sh):
    if np.any(np.asarray(sh)[[1, 2, 4, 5, 6]] < 0):
        return 1e6, None, None
    cls, gap = classes_gaps_of(sh)
    if cls[0] != PAIR or gap[0] < GAP_MIN or cls[2] != NOON:
        return 1e6, None, None
    out = assemble(sh)
    if out is None:
        return 1e6, None, None
    W, lin = out
    Wm = W.mean(1)
    e = float(np.mean(np.sqrt(np.mean((Wm - OBS) ** 2, 1)) / AMP))
    fl = float(np.mean(((Wm.std(1) - OBS.std(1)) / OBS.std(1)) ** 2))
    dm = W.mean(2)
    de = float(np.sqrt(np.mean(((dm - dm.mean(1, keepdims=True)) - ocen) ** 2)))
    return e + 2.0 * fl + 1.5 * de, W, lin


vals = np.array([score(s)[0] for s in seeds])
alive = np.flatnonzero(vals < 1e5)
print(f'  seeds scored, {len(alive):,} valid; best seed {vals[alive].min():.4f}'
      if len(alive) else '  no valid seed')
best = (float(vals[alive].min()), seeds[alive[np.argmin(vals[alive])]].copy())
for i in alive[np.argsort(vals[alive])[:400]]:
    z = seeds[i].copy()
    for _ in range(3):
        r = minimize(lambda q: score(q)[0], z, method='Nelder-Mead',
                     options=dict(maxiter=4000, maxfev=4000, xatol=1e-8, fatol=1e-11))
        z = r.x
    if float(r.fun) < best[0]:
        best = (float(r.fun), z.copy())

v, sh = best
W, (T0, A_T, dr) = assemble(sh)
Wm = W.mean(1)
bL, gL, wL, bK, gK, wK, pK = sh
cls, gap = classes_gaps_of(sh)
print(f'\n=== best HONEST fit (gap>={GAP_MIN}h enforced)   objective {v:.4f} ===')
print(f'  LOSS  bL={bL:+.2f} gL={gL:.2f} wL={wL:.2f}   (no P4)')
print(f'  KEEP  bK={bK:+.2f} gK={gK:.2f} wK={wK:.2f} pK={pK:.2f}')
print(f'  A_T={A_T:.2f}  T0={T0:.2f}  drift={dr:+.3f} C/day')
print(f'  non-preg gap = {gap[0]:.2f} h   classes {[ ["PAIR","NOON","MIDN","OTHER"][k] for k in cls]}')
print(f'\n  {"animal":10s}{"amp mod":>9}{"amp obs":>9}{"SD mod":>8}{"SD obs":>8}{"peak h":>9}')
obs_pk = ['8.2,15.2', '15.8,21.8', '13.2']
for i, l in enumerate(L):
    pk = hrs[np.argmax(Wm[i])]
    print(f'  {l:10s}{np.ptp(Wm[i]):9.2f}{AMP[i]:9.2f}{Wm[i].std():8.3f}{OBS[i].std():8.3f}{pk:9.1f}')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.0), sharex=True)
for i, l in enumerate(L):
    ax[i].plot(hrs, OBS[i], 'k-', lw=2.2, label='mouse')
    ax[i].plot(hrs, Wm[i], 'r--', lw=2, label=f'P4 keep-arm, gap>={GAP_MIN}h enforced')
    ax[i].axvline(12, color='0.85', lw=1, zorder=0)
    ax[i].set_title(f'{l}   amp {np.ptp(Wm[i]):.2f} / {AMP[i]:.2f} C')
    ax[i].set_xlabel('hour of day'); ax[i].set_xticks([0, 6, 12, 18, 24])
ax[0].set_ylabel('CBT (C)'); ax[0].legend(fontsize=8)
fig.suptitle(f'HONEST topology test: non-preg PAIR forced to gap>={GAP_MIN}h (not a knife-edge); '
             'Esr1i forced NOON', y=1.03)
fig.tight_layout(); fig.savefig('fig_p4_keep_honest.png', dpi=130, bbox_inches='tight')
np.save('best_p4_keep_honest.npy', np.concatenate([sh, [A_T, T0, dr]]))
print('\nwrote fig_p4_keep_honest.png')
