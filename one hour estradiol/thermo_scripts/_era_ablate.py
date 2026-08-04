"""ABLATION: is the ERalpha gate load-bearing, or would any 11th parameter have fitted?

Identical seeds, identical objective, identical polish budget. The ONLY difference is whether
rho is free (gate present) or pinned to 1.0 (gate absent -- Esr1i's loss arm keeps its full
circadian drive, i.e. exactly the model we had before the fix).
"""
import numpy as np
from scipy.signal import find_peaks
from scipy.optimize import minimize

fold = np.load('fold.npy', allow_pickle=True).item()
L = ['non-preg', 'pregnant', 'Esr1i']
hrs = fold[L[0]][0]
OBS = np.vstack([fold[l][1] for l in L])
AMP = OBS.max(1) - OBS.min(1); R_OBS = AMP / AMP[0]
tf = np.arange(0, 24, 0.2); NT = tf.size
E2 = np.array([0., 1., 0.])
mid = np.array([13, 14, 15, 16, 17], float) + 0.5
p4b = np.maximum(0., 1. - (mid - 13.) / 6.5).mean(); P4v = np.array([0., p4b, p4b])
DIP = np.array([0.359, 0.321, 0.088]); H_MIN = AMP[0] / 10.


def sig(z):
    return 1 / (1 + np.exp(-np.clip(z, -60, 60)))


def H(p, i, t):
    bL, gL, phiL, wL, pL, bK, gK, phiK, wK, pK, rho = p
    R = rho if i == 2 else 1.0
    cK = np.cos(2 * np.pi * (t - phiK) / 24.); cLl = np.cos(2 * np.pi * (t - phiL) / 24.)
    return (sig(bK + gK * cK + wK * E2[i] + pK * P4v[i])
            - sig(bL + gL * R * cLl + wL * E2[i] + pL * P4v[i]))


def ev(p):
    Hs = np.stack([H(p, i, tf) for i in range(3)])
    amp = Hs.max(1) - Hs.min(1); dips = []
    for i in range(3):
        a = amp[i]
        if a < 1e-12:
            dips.append(0.); continue
        ext = np.concatenate([Hs[i]] * 3)
        pk, pr_ = find_peaks(ext, prominence=1e-9)
        m = (pk >= NT) & (pk < 2 * NT)
        pv = np.sort(pr_['prominences'][m])[::-1]
        dips.append(pv[1] / a if len(pv) >= 2 else 0.)
    return Hs, np.array(dips), amp


def score(p):
    if np.any(np.asarray(p)[[1, 3, 4, 6, 8, 9]] < 0) or not (0 <= p[10] <= 1):
        return 1e6
    Hs, dip, amp = ev(p)
    if amp[0] < H_MIN:
        return 1e6
    A_T = AMP[0] / amp[0]
    Hd = np.stack([H(p, i, hrs) for i in range(3)])
    W = A_T * Hd; W = W + (OBS.mean() - W.mean())
    r = amp / amp[0]
    return (6 * np.sum(np.abs(dip - DIP))
            + 4 * (abs(r[1] - R_OBS[1]) + abs(r[2] - R_OBS[2]))
            + float(np.mean(np.sqrt(np.mean((W - OBS) ** 2, 1)) / AMP)))


def screen(P):
    bL, gL, phiL, wL, pL, bK, gK, phiK, wK, pK = [P[:, i:i + 1] for i in range(10)]
    cK = np.cos(2 * np.pi * (tf[None, :] - phiK) / 24.)
    cL = np.cos(2 * np.pi * (tf[None, :] - phiL) / 24.)
    W = sig(bK + gK * cK) - sig(bL + gL * cL)
    amp = W.max(1) - W.min(1)
    lo = np.roll(W, 1, 1); hi = np.roll(W, -1, 1)
    imx = (W > lo) & (W > hi); imn = (W < lo) & (W < hi)
    ok = (imx.sum(1) == 2) & (imn.sum(1) == 2) & (amp >= H_MIN)
    out = np.zeros(len(W), bool); idx = np.flatnonzero(ok)
    if not idx.size:
        return out
    _, cm = np.nonzero(imx[idx]); pk = cm.reshape(-1, 2)
    _, cn = np.nonzero(imn[idx]); vl = cn.reshape(-1, 2)
    Wt = W[idx]
    pv = np.take_along_axis(Wt, pk, 1); vv = np.take_along_axis(Wt, vl, 1)
    d = (pv.min(1) - vv.max(1)) / amp[idx]
    out[idx] = (d > 0.15) & (d < 0.60)
    return out


rng = np.random.default_rng(777)
N, CH = 1_000_000, 50_000
seeds = []
for s in range(0, N, CH):
    P = np.stack([rng.uniform(-15, 8, CH), rng.uniform(0, 30, CH), rng.uniform(0, 24, CH),
                  rng.uniform(0, 30, CH), rng.uniform(0, 15, CH), rng.uniform(-15, 8, CH),
                  rng.uniform(0, 30, CH), rng.uniform(0, 24, CH), rng.uniform(0, 30, CH),
                  rng.uniform(0, 15, CH)], 1)
    m = screen(P)
    if m.any():
        Q = P[m]
        seeds.append(np.column_stack([Q, rng.uniform(0, 0.6, len(Q))]))
seeds = np.concatenate(seeds)
sub = seeds[:2500]

for tag, pin in [('GATE FREE  (rho fitted)', False), ('GATE OFF   (rho pinned 1.0)', True)]:
    S = sub.copy()
    if pin:
        S[:, 10] = 1.0

    def f(q):
        if pin:
            q = np.concatenate([q[:10], [1.0]])
        return score(q)
    vals = np.array([f(s[:10] if pin else s) for s in S])
    alive = np.flatnonzero(vals < 1e5)
    best = (float(vals[alive].min()), None)
    for i in alive[np.argsort(vals[alive])][:45]:
        z = (S[i][:10] if pin else S[i]).copy()
        for _ in range(3):
            r = minimize(f, z, method='Nelder-Mead',
                         options=dict(maxiter=5000, maxfev=5000, xatol=1e-9, fatol=1e-12))
            z = r.x
        if float(r.fun) < best[0]:
            best = (float(r.fun), z.copy())
    q = best[1]
    if pin:
        q = np.concatenate([q[:10], [1.0]])
    Hs, dip, amp = ev(q)
    print(f'\n{tag}:  best objective {best[0]:.4f}', flush=True)
    print(f'   dips  {np.round(dip,3)}   target {DIP}', flush=True)
    print(f'   ratios 1.00/{amp[1]/amp[0]:.2f}/{amp[2]/amp[0]:.2f}  target 1.00/{R_OBS[1]:.2f}/{R_OBS[2]:.2f}', flush=True)
    for i, l in enumerate(L):
        a = Hs[i].max() - Hs[i].min()
        ext = np.concatenate([Hs[i]] * 3)
        pk, _ = find_peaks(ext, prominence=0.20 * a)
        m = (pk >= NT) & (pk < 2 * NT)
        print(f'   {l:10s} {int(m.sum())} peak(s)', flush=True)
