"""Does P4-in-keep-only buy Esr1i's unimodality at the cost of the amplitudes?

Two fits of the SAME model disagree:
  unconstrained (_fit_p4_keep.py)  amps 2.33/0.51/1.24 vs obs 2.10/0.72/1.16  -- but Esr1i BIMODAL
  topology-forced (_p4keep_strict) Esr1i unimodal                             -- but amps collapse

So sweep the topology penalty from 0 to hard and watch both quantities at once.  If the
two are genuinely exclusive the curve is monotone: every unit of Esr1i unimodality costs
non-preg amplitude.  Also report, over ALL topology-valid random draws, the best
attainable amplitude-ratio match -- that separates "search failed" from "structurally out".
"""
import numpy as np

fold = np.load('fold.npy', allow_pickle=True).item()
L = ['non-preg', 'pregnant', 'Esr1i']
hrs = fold[L[0]][0]
OBS = np.vstack([fold[l][1] for l in L])
AMP = OBS.max(1) - OBS.min(1)
R_OBS = AMP / AMP[0]

cg = np.linspace(-1, 1, 401)
E2 = np.array([0., 1., 0.])
PREG = np.array([0., 1., 1.])
mid = np.array([13, 14, 15, 16, 17], float) + 0.5
p4bar = np.maximum(0., 1. - (mid - 13.) / 6.5).mean()
PAIR, NOON, MIDN, OTHER = 0, 1, 2, 3


def sig(z):
    return 1 / (1 + np.exp(-np.clip(z, -60, 60)))


def dsig(z):
    s = sig(z)
    return s * (1 - s)


def shape_class(bk, gK, bl, gL):
    uK = bk[..., None] + gK[..., None] * cg
    uL = bl[..., None] + gL[..., None] * cg
    d = gK[..., None] * dsig(uK) - gL[..., None] * dsig(uL)
    pos = d > 0
    drops = np.count_nonzero(pos[..., :-1] & ~pos[..., 1:], axis=-1)
    up, dn = pos[..., -1], ~pos[..., 0]
    cls = np.full(drops.shape, OTHER, np.int8)
    cls = np.where((drops == 1) & ~up & ~dn, PAIR, cls)
    cls = np.where((drops == 0) & up & ~dn, NOON, cls)
    cls = np.where((drops == 0) & dn & ~up, MIDN, cls)
    f = sig(uK) - sig(uL)
    return cls, f.max(-1) - f.min(-1)


rng = np.random.default_rng(5150)
N, CH = 12_000_000, 200_000
best_ratio_err = np.inf
best_rec = None
n_both = 0
amp_np_valid = []
for s in range(0, N, CH):
    bK = rng.uniform(-30, 10, CH); gK = rng.uniform(0, 60, CH)
    bL = rng.uniform(-30, 10, CH); gL = rng.uniform(0, 60, CH)
    pK = rng.uniform(0, 25, CH)
    wK = rng.uniform(0, 60, CH);  wL = rng.uniform(0, 60, CH)
    c_np, a_np = shape_class(bK, gK, bL, gL)
    c_es, a_es = shape_class(bK + pK * p4bar, gK, bL, gL)
    _, a_pr = shape_class(bK + wK + pK * p4bar, gK, bL + wL, gL)
    m = (c_np == PAIR) & (c_es == NOON)
    if not m.any():
        continue
    n_both += int(m.sum())
    an, ap, ae = a_np[m], a_pr[m], a_es[m]
    amp_np_valid.append(an)
    good = an > 1e-6
    err = np.where(good,
                   np.abs(ap / np.maximum(an, 1e-12) - R_OBS[1])
                   + np.abs(ae / np.maximum(an, 1e-12) - R_OBS[2]), np.inf)
    j = int(np.argmin(err))
    if err[j] < best_ratio_err:
        best_ratio_err = float(err[j])
        idx = np.flatnonzero(m)[j]
        best_rec = dict(bL=bL[idx], gL=gL[idx], wL=wL[idx], bK=bK[idx], gK=gK[idx],
                        wK=wK[idx], pK=pK[idx], a=(an[j], ap[j], ae[j]))

amp_np_valid = np.concatenate(amp_np_valid)
print(f'draws {N:,}')
print(f'  non-preg PAIR and Esr1i NOON simultaneously : {n_both:,}')
print(f'\n  arm-unit amplitude of the non-preg waveform among those:')
for q in (50, 90, 99, 99.9, 100):
    print(f'    p{q:<5} {np.percentile(amp_np_valid, q):.4f}')
print(f'  (arm units are the sigmoid difference, max possible 1.0;'
      f' Celsius = A_T x this)')

print(f'\n  best amplitude-RATIO match inside the topology-valid set:')
an, ap, ae = best_rec['a']
print(f'    model ratios 1.00/{ap/an:.2f}/{ae/an:.2f}   observed 1.00/'
      f'{R_OBS[1]:.2f}/{R_OBS[2]:.2f}   |err| {best_ratio_err:.3f}')
print(f'    at bL={best_rec["bL"]:+.2f} gL={best_rec["gL"]:.2f} wL={best_rec["wL"]:.2f}'
      f' | bK={best_rec["bK"]:+.2f} gK={best_rec["gK"]:.2f} wK={best_rec["wK"]:.2f}'
      f' pK={best_rec["pK"]:.2f}')
print(f'    non-preg arm amplitude there: {an:.4f}'
      f'  -> needs A_T = {AMP[0]/an:.1f} C to reach the observed 2.10 C')
