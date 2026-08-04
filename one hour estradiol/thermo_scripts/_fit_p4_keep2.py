"""Fit the P4-in-KEEP-arm-only model WITH the observed peak topology enforced.

Model (all coefficients positive, one noon-locked clock):
    c    = cos(2*pi*(t-12)/24)
    loss = sig(bL + gL*c + wL*E2)                 <- P4 absent, it is not a heat-loss signal
    keep = sig(bK + gK*c + wK*E2 + pK*P4)         <- P4 = cutaneous vasoconstriction / heat retention
    CBT  = T0 + A_T*(keep - loss) + drift*(day-15.5)
    E2 = 0/1/0 ,  P4 = 0/line/line   for non-preg / pregnant / Esr1i

Topology is a sign condition on f'(c) = gK*sig'(uK) - gL*sig'(uL):
  bimodal  <=> f' changes sign inside [-1,1]      (min f' < 0)
  unimodal <=> f' keeps one sign                  (min f' > 0, peak then at NOON)
so it goes into the objective as a smooth hinge instead of a peak count.
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

c = np.cos(2 * np.pi * (hrs - 12.) / 24.)
cg = np.linspace(-1, 1, 241)                      # grid for the sign condition
E2 = np.array([0., 1., 0.])
PREG = np.array([0., 1., 1.])
mid = np.array([13, 14, 15, 16, 17], float) + 0.5
P4d = np.maximum(0., 1. - (mid - 13.) / 6.5)
P4M = np.outer(PREG, P4d)[:, :, None]
dd = (mid - 15.5)[None, :, None]
cb = c[None, None, :]
WANT_UNIMODAL = np.array([False, False, True])    # non-preg, pregnant, Esr1i
MARGIN = 2e-3


def sig(z):
    return 1 / (1 + np.exp(-np.clip(z, -60, 60)))


def dsig(z):
    s = sig(z)
    return s * (1 - s)


def unpack(z):
    bL, gL, wL, bK, gK, wK, pK, A_T, T0, dr = z
    return bL, gL, wL, bK, gK, wK, pK, A_T, T0, dr


def waves(z):
    bL, gL, wL, bK, gK, wK, pK, A_T, T0, dr = unpack(z)
    e2 = E2[:, None, None]
    return (T0 + A_T * (sig(bK + gK * cb + wK * e2 + pK * P4M)
                        - sig(bL + gL * cb + wL * e2)) + dr * dd)


def min_slope(z):
    """min over c of f'(c), per animal, using the mid-window P4 level."""
    bL, gL, wL, bK, gK, wK, pK, A_T, T0, dr = unpack(z)
    p4 = PREG * P4d.mean()
    uK = (bK + wK * E2 + pK * p4)[:, None] + gK * cg[None, :]
    uL = (bL + wL * E2)[:, None] + gL * cg[None, :]
    return (gK * dsig(uK) - gL * dsig(uL)).min(1)


def obj(z, w_top):
    W = waves(z)
    Wm = W.mean(axis=1)
    r = Wm - OBS
    e = float(np.mean(np.sqrt(np.mean(r ** 2, 1)) / AMP))
    fl = float(np.mean(((Wm.std(1) - OBS.std(1)) / OBS.std(1)) ** 2))
    edge = sum(1.0 for i in range(3)
               for k in (int(np.argmax(Wm[i])), int(np.argmin(Wm[i])))
               if k < 2 or k > NB - 3)
    dm = W.mean(axis=2)
    de = float(np.sqrt(np.mean(((dm - dm.mean(1, keepdims=True))
                                - (omean - omean.mean(1, keepdims=True))) ** 2)))
    m = min_slope(z)
    top = np.where(WANT_UNIMODAL,
                   np.maximum(0., MARGIN - m),      # want min f' > 0
                   np.maximum(0., MARGIN + m))      # want min f' < 0
    return e + 2.0 * fl + 1.0 * edge + 1.5 * de + w_top * float(top.sum())


B = [(-30, 10), (0, 60), (0, 60), (-30, 10), (0, 60), (0, 60),
     (0, 25), (0.3, 8), (30, 42), (-0.4, 0.1)]
lo = np.array([q[0] for q in B]); hi = np.array([q[1] for q in B])
rng = np.random.default_rng(4242)
best = None
for _ in range(1200):
    z0 = lo + rng.random(len(B)) * (hi - lo)
    z = z0
    for w in (30., 300.):                          # ramp the topology weight
        r = minimize(obj, z, args=(w,), bounds=B, method='L-BFGS-B',
                     options=dict(maxiter=600, ftol=1e-13))
        z = r.x
    m = min_slope(z)
    if not (np.all(np.where(WANT_UNIMODAL, m > 0, m < 0))):
        continue                                    # topology not achieved, discard
    val = obj(z, 0.0)
    if best is None or val < best[0]:
        best = (float(val), z.copy())

if best is None:
    raise SystemExit('no start reached the requested topology')

e, z = best
np.save('best_p4_keep_topo.npy', z)
W = waves(z); Wm = W.mean(axis=1); dm = W.mean(axis=2)
np.save('_W_p4keep.npy', W)
bL, gL, wL, bK, gK, wK, pK, A_T, T0, dr = unpack(z)
print(f'fit loss (topology term excluded) {e:.4f}')
print(f'  LOSS  bL={bL:+.2f} gL={gL:.2f} wL={wL:.2f}   (no P4)')
print(f'  KEEP  bK={bK:+.2f} gK={gK:.2f} wK={wK:.2f} pK={pK:.2f}')
print(f'  A_T={A_T:.2f}  T0={T0:.2f}  drift={dr:+.3f} C/day')
print(f'\n  {"animal":10s}{"amp mod":>9}{"amp obs":>9}{"SD mod":>8}{"SD obs":>8}'
      f'{"drift mod":>11}{"drift obs":>10}')
for i, l in enumerate(L):
    print(f'  {l:10s}{np.ptp(Wm[i]):9.2f}{AMP[i]:9.2f}{Wm[i].std():8.3f}{OBS[i].std():8.3f}'
          f'{dm[i][-1]-dm[i][0]:+11.2f}{omean[i][-1]-omean[i][0]:+10.2f}')


def peaks_of(w, p):
    a = np.ptp(w); ext = np.concatenate([w] * 3)
    pk, _ = find_peaks(ext, prominence=p * max(a, 1e-9))
    return [round(float(hrs[q - NB]), 1) for q in pk if NB <= q < 2 * NB]


print('\n  peaks across prominence thresholds:')
for i, l in enumerate(L):
    print(f'    {l:10s}' + '  '.join(f'{p:.2f}:{peaks_of(Wm[i], p)}'
                                     for p in (0.05, 0.10, 0.15, 0.25)))
print('  OBSERVED    non-preg [8.2, 15.2] | pregnant [15.8, 21.8] | Esr1i [13.2] ONLY')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 3, figsize=(13, 3.8), sharex=True)
for i, l in enumerate(L):
    ax[i].plot(hrs, OBS[i], 'k-', lw=2, label='mouse')
    ax[i].plot(hrs, Wm[i], 'r--', lw=2, label='model')
    ax[i].set_title(f'{l}   amp {np.ptp(Wm[i]):.2f} vs {AMP[i]:.2f} C')
    ax[i].set_xlabel('hour of day'); ax[i].set_xticks([0, 6, 12, 18, 24])
    ax[i].axvline(12, color='0.8', lw=1, zorder=0)
ax[0].set_ylabel('CBT (C)'); ax[0].legend(fontsize=8)
fig.suptitle('P4 in the heat-KEEP arm only (vasoconstriction), E2=0 in Esr1i', y=1.02)
fig.tight_layout(); fig.savefig('fig_p4_keep_only.png', dpi=130, bbox_inches='tight')
print('\nwrote fig_p4_keep_only.png')
