"""P4 in the heat-KEEP arm only -- with the CORRECT bimodality criterion.

Earlier counters scored "peak at noon + peak at midnight" as bimodal.  That is the wrong
shape: the observed non-preg pair is 8.2/15.2 h, i.e. an INTERIOR maximum of f(c) visited
twice per clock sweep, with noon a local MINIMUM between the two humps.

Under noon-lock the waveform is f(c), c = cos(2pi(t-12)/24), so with
    f'(c) = gK*sig'(uK) - gL*sig'(uL)
    drops  = number of + -> - sign changes of f' inside (-1,1)   (interior maxima)
    up_noon = f'(1) > 0     dn_mid = f'(-1) < 0
the three shapes we care about are

    PAIR      drops == 1 and not up_noon and not dn_mid   -> 2 peaks straddling noon
    NOON      drops == 0 and up_noon and not dn_mid       -> 1 peak, at noon
    MIDNIGHT  drops == 0 and dn_mid and not up_noon       -> 1 peak, at midnight

Target: non-preg PAIR, Esr1i NOON.  (Pregnant's observed pair straddles 18.8 h, which
noon-lock forbids outright, so it is reported but never required.)
"""
import numpy as np
from scipy.optimize import minimize

fold = np.load('fold.npy', allow_pickle=True).item()
L = ['non-preg', 'pregnant', 'Esr1i']
hrs = fold[L[0]][0]
OBS = np.vstack([fold[l][1] for l in L])
AMP = OBS.max(1) - OBS.min(1)
omean = np.load('_omean.npy')
ocen = omean - omean.mean(1, keepdims=True)

c48 = np.cos(2 * np.pi * (hrs - 12.) / 24.)
cg = np.linspace(-1, 1, 401)
E2 = np.array([0., 1., 0.])
PREG = np.array([0., 1., 1.])
mid = np.array([13, 14, 15, 16, 17], float) + 0.5
P4d = np.maximum(0., 1. - (mid - 13.) / 6.5)
P4M = np.outer(PREG, P4d)
dd = mid - 15.5

PAIR, NOON, MIDN, OTHER = 0, 1, 2, 3


def sig(z):
    return 1 / (1 + np.exp(-np.clip(z, -60, 60)))


def dsig(z):
    s = sig(z)
    return s * (1 - s)


def shape_class(bk, gK, bl, gL):
    """classify f(c)=sig(bk+gK c)-sig(bl+gL c); arrays broadcast over leading axis."""
    uK = bk[..., None] + gK[..., None] * cg
    uL = bl[..., None] + gL[..., None] * cg
    d = gK[..., None] * dsig(uK) - gL[..., None] * dsig(uL)
    pos = d > 0
    drops = np.count_nonzero(pos[..., :-1] & ~pos[..., 1:], axis=-1)
    up_noon = pos[..., -1]
    dn_mid = ~pos[..., 0]
    cls = np.full(drops.shape, OTHER, np.int8)
    cls = np.where((drops == 1) & ~up_noon & ~dn_mid, PAIR, cls)
    cls = np.where((drops == 0) & up_noon & ~dn_mid, NOON, cls)
    cls = np.where((drops == 0) & dn_mid & ~up_noon, MIDN, cls)
    f = sig(uK) - sig(uL)
    return cls, f.max(-1) - f.min(-1)


# ---------------------------------------------------------------- structural scan
rng = np.random.default_rng(8080)
N, CH = 8_000_000, 200_000
p4bar = P4d.mean()
cnt = dict(np_pair=0, es_noon=0, both=0, both_amp=0)
seeds = []
for s in range(0, N, CH):
    bK = rng.uniform(-30, 10, CH); gK = rng.uniform(0, 60, CH)
    bL = rng.uniform(-30, 10, CH); gL = rng.uniform(0, 60, CH)
    pK = rng.uniform(0, 25, CH)
    wK = rng.uniform(0, 60, CH);  wL = rng.uniform(0, 60, CH)
    c_np, a_np = shape_class(bK, gK, bL, gL)
    c_es, a_es = shape_class(bK + pK * p4bar, gK, bL, gL)
    m1 = c_np == PAIR
    m2 = c_es == NOON
    both = m1 & m2
    ok_amp = both & (np.abs(a_es / np.maximum(a_np, 1e-12) - AMP[2] / AMP[0]) < 0.10) \
                  & (a_np > 0.05)
    cnt['np_pair'] += int(m1.sum()); cnt['es_noon'] += int(m2.sum())
    cnt['both'] += int(both.sum()); cnt['both_amp'] += int(ok_amp.sum())
    if ok_amp.any():
        seeds.append(np.stack([bL[ok_amp], gL[ok_amp], wL[ok_amp], bK[ok_amp],
                               gK[ok_amp], wK[ok_amp], pK[ok_amp]], 1))
seeds = np.concatenate(seeds) if seeds else np.zeros((0, 7))
print(f'draws {N:,}')
print(f'  non-preg = PAIR straddling noon              {cnt["np_pair"]:,}')
print(f'  Esr1i    = single NOON peak                  {cnt["es_noon"]:,}')
print(f'  BOTH                                         {cnt["both"]:,}')
print(f'  BOTH + amp ratio {AMP[2]/AMP[0]:.2f}+-0.10 & non-deg      {cnt["both_amp"]:,}')

# ---------------------------------------------------------------- fit the survivors
cb = c48


def H3_of(sh):
    bL, gL, wL, bK, gK, wK, pK = sh
    e2 = E2[:, None, None]; p4 = P4M[:, :, None]
    return (sig(bK + gK * cb + wK * e2 + pK * p4) - sig(bL + gL * cb + wL * e2))


def classes_of(sh):
    bL, gL, wL, bK, gK, wK, pK = sh
    bk = np.array([bK, bK + wK + pK * p4bar, bK + pK * p4bar])
    bl = np.array([bL, bL + wL, bL])
    g1 = np.full(3, gK); g2 = np.full(3, gL)
    return shape_class(bk, g1, bl, g2)[0]


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
    cls = classes_of(sh)
    if cls[0] != PAIR or cls[2] != NOON:
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
Wm = W.mean(1); dm = W.mean(2)
np.save('best_p4_keep_pairnoon.npy', np.concatenate([sh, [A_T, T0, dr]]))
bL, gL, wL, bK, gK, wK, pK = sh
print(f'\n=== best fit with non-preg=PAIR, Esr1i=NOON   objective {v:.4f} ===')
print(f'  LOSS  bL={bL:+.2f} gL={gL:.2f} wL={wL:.2f}   (no P4)')
print(f'  KEEP  bK={bK:+.2f} gK={gK:.2f} wK={wK:.2f} pK={pK:.2f}')
print(f'  A_T={A_T:.2f}  T0={T0:.2f}  drift={dr:+.3f} C/day')
print(f'  shape classes {[ ["PAIR","NOON","MIDN","OTHER"][k] for k in classes_of(sh)]}')
print(f'\n  {"animal":10s}{"amp mod":>9}{"amp obs":>9}{"SD mod":>8}{"SD obs":>8}'
      f'{"peak h":>9}{"obs peaks":>18}')
obs_pk = ['8.2,15.2', '15.8,21.8', '13.2']
for i, l in enumerate(L):
    pk = hrs[np.argmax(Wm[i])]
    print(f'  {l:10s}{np.ptp(Wm[i]):9.2f}{AMP[i]:9.2f}{Wm[i].std():8.3f}'
          f'{OBS[i].std():8.3f}{pk:9.1f}{obs_pk[i]:>18}')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.0), sharex=True)
for i, l in enumerate(L):
    ax[i].plot(hrs, OBS[i], 'k-', lw=2.2, label='mouse')
    ax[i].plot(hrs, Wm[i], 'r--', lw=2, label='P4 in keep arm only')
    ax[i].axvline(12, color='0.85', lw=1, zorder=0)
    ax[i].set_title(f'{l}   amp {np.ptp(Wm[i]):.2f} / {AMP[i]:.2f} C')
    ax[i].set_xlabel('hour of day'); ax[i].set_xticks([0, 6, 12, 18, 24])
ax[0].set_ylabel('CBT (C)'); ax[0].legend(fontsize=8)
fig.suptitle('P4 = vasoconstriction in the heat-KEEP arm only; E2 = 0 in Esr1i; '
             'non-preg forced to a noon-straddling pair, Esr1i to a single noon peak', y=1.03)
fig.tight_layout(); fig.savefig('fig_p4_keep_pairnoon.png', dpi=130, bbox_inches='tight')
print('\nwrote fig_p4_keep_pairnoon.png')
