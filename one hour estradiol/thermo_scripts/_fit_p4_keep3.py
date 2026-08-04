"""P4 in the heat-KEEP arm only, fitted with the observed peak topology as a HARD filter.

    c    = cos(2*pi*(t-12)/24)
    loss = sig(bL + gL*c + wL*E2)                <- no P4: P4 is not a heat-loss signal
    keep = sig(bK + gK*c + wK*E2 + pK*P4)        <- P4 = vasoconstriction / heat retention
    CBT  = T0 + A_T*(keep - loss) + drift*(day-15.5)
    E2 = 0/1/0,  P4 = 0/line/line   for non-preg / pregnant / Esr1i   (all coeffs >= 0)

Two things the earlier attempts got wrong and this fixes:
  * topology is counted DISCRETELY on the same 48-bin fold as the data, so it cannot be
    satisfied inside a flat tail that carries no amplitude (the smooth hinge let that through);
  * T0, A_T and drift enter the model LINEARLY, so for every shape draw they are solved in
    closed form instead of being searched -- the random scan then only spans the 7 shape
    parameters, which is what makes an exhaustive filtered search affordable.
Search = filtered random scan for seeds, then Nelder-Mead polish with topology rejection.
"""
import numpy as np
from scipy.optimize import minimize
from scipy.signal import find_peaks

fold = np.load('fold.npy', allow_pickle=True).item()
L = ['non-preg', 'pregnant', 'Esr1i']
hrs = fold[L[0]][0]
OBS = np.vstack([fold[l][1] for l in L])
AMP = OBS.max(1) - OBS.min(1)
NB = len(hrs)
omean = np.load('_omean.npy')
ocen = omean - omean.mean(1, keepdims=True)

c = np.cos(2 * np.pi * (hrs - 12.) / 24.)
E2 = np.array([0., 1., 0.])
PREG = np.array([0., 1., 1.])
mid = np.array([13, 14, 15, 16, 17], float) + 0.5
P4d = np.maximum(0., 1. - (mid - 13.) / 6.5)
P4M = np.outer(PREG, P4d)                       # (3,5)
dd = mid - 15.5
WANT = np.array([2, 2, 1])                      # peaks/day: non-preg, pregnant, Esr1i
PROM = 0.10
# Pregnant's peak PLACEMENT is unreachable under noon-lock whatever we do (its observed
# pair straddles 18.8 h, and mirror symmetry forces every pair to straddle noon), so the
# strict 2/2/1 filter buys an impossible feature at the price of all the amplitude.
# STRICT constrains all three; RELAXED constrains only the two animals whose peak counts
# this change is actually about -- non-preg bimodal, Esr1i unimodal.
FREE = np.array([False, True, False])           # pregnant's count left free in RELAXED


def sig(z):
    return 1 / (1 + np.exp(-np.clip(z, -60, 60)))


def H_of(sh, p4, e2):
    """heat balance keep-loss.  Two shapes are used:
       batch  sh (M,7), p4 (3,1),   e2 (3,1)   -> (M,3,48)
       days   sh (7,),  p4 (3,5,1), e2 (3,1,1) -> (3,5,48)
    """
    bL, gL, wL, bK, gK, wK, pK = [sh[..., i, None, None] for i in range(7)]
    cc = c
    return (sig(bK + gK * cc + wK * e2 + pK * p4)
            - sig(bL + gL * cc + wL * e2))


def npeaks_grid(W):
    """peaks/day on the circular 48-bin waveform, vectorised over leading axis."""
    a = np.ptp(W, axis=-1, keepdims=True)
    lo = np.roll(W, 1, -1)
    hi = np.roll(W, -1, -1)
    ismax = (W > lo) & (W >= hi)
    # prominence proxy: height above the lower of the two neighbouring troughs.
    # on 48 bins with <=4 extrema, depth-to-global-min scaled by amplitude is adequate
    depth = (W - W.min(-1, keepdims=True)) / np.maximum(a, 1e-12)
    return np.count_nonzero(ismax & (depth > PROM), axis=-1)


def solve_linear(H3):
    """given H (3,5,48) return best T0, A_T, drift and the fitted waveform set."""
    Hf = H3.mean(1)                                   # fold average, (3,48)
    X = np.stack([np.ones(Hf.size), Hf.ravel()], 1)
    coef, *_ = np.linalg.lstsq(X, OBS.ravel(), rcond=None)
    T0, A_T = coef
    if A_T <= 0:
        return None
    hm = H3.mean(2)                                   # (3,5) day means
    hcen = hm - hm.mean(1, keepdims=True)
    den = 3.0 * float((dd ** 2).sum())                # 3 animals share one drift
    dr = float(((ocen - A_T * hcen) * dd[None, :]).sum() / den)
    return T0, A_T, dr


def assemble(sh, lin=None):
    H3 = H_of(sh, P4M[:, :, None], E2[:, None, None])
    r = lin or solve_linear(H3)
    if r is None:
        return None
    T0, A_T, dr = r
    W = T0 + A_T * H3 + dr * dd[None, :, None]
    return W, (T0, A_T, dr)


def topo_ok(npk, relaxed):
    m = ~FREE if relaxed else np.ones(3, bool)
    return bool(np.all(npk[..., m] == WANT[m], axis=-1))


def score(sh, relaxed=False):
    out = assemble(sh)
    if out is None:
        return 1e6, None, None
    W, lin = out
    Wm = W.mean(1)
    if not topo_ok(npeaks_grid(Wm), relaxed):
        return 1e6, W, lin
    e = float(np.mean(np.sqrt(np.mean((Wm - OBS) ** 2, 1)) / AMP))
    fl = float(np.mean(((Wm.std(1) - OBS.std(1)) / OBS.std(1)) ** 2))
    dm = W.mean(2)
    de = float(np.sqrt(np.mean(((dm - dm.mean(1, keepdims=True)) - ocen) ** 2)))
    return e + 2.0 * fl + 1.5 * de, W, lin


# ---------------- filtered random scan over the 7 shape parameters ----------------
LOB = np.array([-30, 0, 0, -30, 0, 0, 0.])
HIB = np.array([10, 60, 60, 10, 60, 60, 25.])
rng = np.random.default_rng(31415)
N, CH = 8_000_000, 100_000
p4m = P4M.mean(1)[:, None]                     # day-averaged P4 for the filter
keepS, keepR = [], []
for s in range(0, N, CH):
    sh = LOB + rng.random((CH, 7)) * (HIB - LOB)
    npk = npeaks_grid(H_of(sh, p4m, E2[:, None]))   # topology is scale-free -> test on H
    hitR = np.all(npk[:, ~FREE] == WANT[~FREE][None, :], axis=1)
    hitS = hitR & (npk[:, 1] == WANT[1])
    if hitR.any():
        keepR.append(sh[hitR])
    if hitS.any():
        keepS.append(sh[hitS])
seedR = np.concatenate(keepR) if keepR else np.zeros((0, 7))
seedS = np.concatenate(keepS) if keepS else np.zeros((0, 7))
print(f'random draws {N:,}')
print(f'  non-preg bimodal + Esr1i unimodal          {len(seedR):,}')
print(f'  ... and pregnant bimodal too (strict 2/2/1) {len(seedS):,}')


def polish(seeds, relaxed, ntop=250):
    if not len(seeds):
        return None
    vals = np.array([score(sh, relaxed)[0] for sh in seeds])
    order = np.argsort(vals)[:ntop]

    def nm(sh):
        if np.any(sh[[1, 2, 4, 5, 6]] < 0):
            return 1e6
        return score(sh, relaxed)[0]

    best = (float(vals[order[0]]), seeds[order[0]].copy())
    for i in order:
        z = seeds[i].copy()
        for _ in range(3):
            r = minimize(nm, z, method='Nelder-Mead',
                         options=dict(maxiter=4000, maxfev=4000,
                                      xatol=1e-8, fatol=1e-11))
            z = r.x
        if float(r.fun) < best[0]:
            best = (float(r.fun), z.copy())
    return best


def peaks_of(w, p):
    a = np.ptp(w); ext = np.concatenate([w] * 3)
    pk, _ = find_peaks(ext, prominence=p * max(a, 1e-9))
    return [round(float(hrs[q - NB]), 1) for q in pk if NB <= q < 2 * NB]


def report(tag, best, fname):
    v, sh = best
    W, (T0, A_T, dr) = assemble(sh)
    Wm = W.mean(1); dm = W.mean(2)
    bL, gL, wL, bK, gK, wK, pK = sh
    print(f'\n================ {tag}   objective {v:.4f} ================')
    print(f'  LOSS  bL={bL:+.2f} gL={gL:.2f} wL={wL:.2f}    (P4 absent by construction)')
    print(f'  KEEP  bK={bK:+.2f} gK={gK:.2f} wK={wK:.2f} pK={pK:.2f}')
    print(f'  A_T={A_T:.2f}  T0={T0:.2f}  drift={dr:+.3f} C/day')
    print(f'\n  {"animal":10s}{"amp mod":>9}{"amp obs":>9}{"SD mod":>8}{"SD obs":>8}'
          f'{"drift mod":>11}{"drift obs":>10}')
    for i, l in enumerate(L):
        print(f'  {l:10s}{np.ptp(Wm[i]):9.2f}{AMP[i]:9.2f}{Wm[i].std():8.3f}'
              f'{OBS[i].std():8.3f}{dm[i][-1]-dm[i][0]:+11.2f}'
              f'{omean[i][-1]-omean[i][0]:+10.2f}')
    print('\n  peaks across prominence thresholds:')
    for i, l in enumerate(L):
        print(f'    {l:10s}' + '  '.join(f'{p:.2f}:{peaks_of(Wm[i], p)}'
                                         for p in (0.05, 0.10, 0.15, 0.25)))
    print('  OBSERVED    non-preg [8.2, 15.2] | pregnant [15.8, 21.8] | Esr1i [13.2] ONLY')
    np.save(fname.replace('fig_', 'best_').replace('.png', '.npy'),
            np.concatenate([sh, [A_T, T0, dr]]))
    return Wm


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

bR = polish(seedR, True)
bS = polish(seedS, False)
WR = report('RELAXED  (non-preg 2, Esr1i 1; pregnant count free)', bR, 'fig_p4_keep_only.png')
WS = report('STRICT   (2 / 2 / 1 all enforced)', bS, 'fig_p4_keep_strict.png') if bS else None

fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.0), sharex=True)
for i, l in enumerate(L):
    ax[i].plot(hrs, OBS[i], 'k-', lw=2.2, label='mouse')
    ax[i].plot(hrs, WR[i], 'r--', lw=2, label='P4 in keep arm (relaxed)')
    if WS is not None:
        ax[i].plot(hrs, WS[i], 'b:', lw=2, label='strict 2/2/1')
    ax[i].axvline(12, color='0.85', lw=1, zorder=0)
    ax[i].set_title(f'{l}   amp {np.ptp(WR[i]):.2f} / {AMP[i]:.2f} C')
    ax[i].set_xlabel('hour of day'); ax[i].set_xticks([0, 6, 12, 18, 24])
ax[0].set_ylabel('CBT (C)'); ax[0].legend(fontsize=7.5)
fig.suptitle('P4 in the heat-KEEP arm only (vasoconstriction, no heat-loss term); '
             'E2 = 0 in Esr1i', y=1.03)
fig.tight_layout(); fig.savefig('fig_p4_keep_only.png', dpi=130, bbox_inches='tight')
print('\nwrote fig_p4_keep_only.png')
