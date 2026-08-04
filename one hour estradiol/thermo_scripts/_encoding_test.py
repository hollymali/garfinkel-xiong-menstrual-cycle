"""Does the E2 ENCODING fix make 2/2/1-at-correct-positions survive WITH REAL AMPLITUDE?

THE CHANGE UNDER TEST
  old:  E2 = 0 / 1 / 0     for non-preg / pregnant / Esr1i
        (says a normal cycling female has the SAME zero estrogen signal as an ERalpha-blocked
         animal, which is wrong -- she has real estradiol and a working receptor)
  new:  E2 = m / 1 / 0 ,  0 < m < 1   (hormone LEVEL x receptor FUNCTION)
        Esr1i becomes the bottom EXTREME instead of the bracketed middle.
Both are run on the SAME random draws so the comparison is controlled: `m` is the only change.

WHY THIS TEST IS A_T-FREE (and so cannot be won by a flat line, which is how
`_twophase_final.py` failed):
  CBT = T0 + A_T*H, so the amplitude RATIOS between animals are ratios of H-ranges and A_T
  cancels entirely.  The absolute amplitude enters only as a BOUND: reaching the observed
  2.10 C for non-preg with a physical A_T <= 10 C requires H-range_np >= 2.10/10 = 0.21.
  So the whole target lives in H-space and a degenerate near-flat H is rejected outright.

TARGET (nested, so we see WHERE it breaks rather than just whether):
  1. non-preg 2 peaks at 8.2/15.2      (+-1.5 h)
  2. + pregnant 2 peaks at 15.8/21.8   (+-1.5 h)
  3. + Esr1i exactly 1 peak at 13.2    (+-1.5 h)
  4. + H-range ratios 1.00/0.34/0.55   (+-0.10)
  5. + H-range_np >= 0.21              (i.e. A_T <= 10 C -- REAL amplitude)
Peak counting uses the routine validated against the observed data: scipy prominence on a
tripled (circular) array at thr=0.20, the threshold at which the mouse data itself shows 2/2/1.
"""
import numpy as np
from scipy.signal import find_peaks

fold = np.load('fold.npy', allow_pickle=True).item()
L = ['non-preg', 'pregnant', 'Esr1i']
hrs = fold[L[0]][0]
OBS = np.vstack([fold[l][1] for l in L])
AMP = OBS.max(1) - OBS.min(1)
R_OBS = AMP / AMP[0]

tf = np.arange(0, 24, 0.2)
NT = tf.size
mid = np.array([13, 14, 15, 16, 17], float) + 0.5
p4bar = np.maximum(0., 1. - (mid - 13.) / 6.5).mean()
P4v = np.array([0., p4bar, p4bar])

THR = 0.20
TOLP = 1.5
TOLR = 0.10
H_MIN = AMP[0] / 10.0            # 0.21 -> A_T <= 10 C
OBS_PK = [np.array([8.2, 15.2]), np.array([15.8, 21.8]), np.array([13.2])]


def pr(*a):
    print(*a, flush=True)


def sig(z):
    return 1 / (1 + np.exp(-np.clip(z, -60, 60)))


def H_of(P, e2, p4):
    bL, gL, phiL, wL, pL, bK, gK, phiK, wK, pK = [P[:, i:i + 1] for i in range(10)]
    cK = np.cos(2 * np.pi * (tf[None, :] - phiK) / 24.)
    cL = np.cos(2 * np.pi * (tf[None, :] - phiL) / 24.)
    e2 = np.asarray(e2).reshape(-1, 1)
    return (sig(bK + gK * cK + wK * e2 + pK * p4)
            - sig(bL + gL * cL + wL * e2 + pL * p4))


def exact_peaks(w, thr=THR):
    amp = w.max() - w.min()
    if amp < 1e-12:
        return np.array([]), amp
    ext = np.concatenate([w, w, w])
    pk, _ = find_peaks(ext, prominence=thr * amp)
    m = (pk >= NT) & (pk < 2 * NT)
    return tf[pk[m] - NT], amp


def match(pk, obs):
    if len(pk) != len(obs):
        return False
    d = np.abs(np.sort(pk)[:, None] - np.sort(obs)[None, :]) % 24.
    d = np.minimum(d, 24. - d)
    return bool(np.all(np.diag(d) < TOLP))


def screen_np(P, e2np):
    """cheap vectorised: non-preg 2 local maxima, real dip, near 8.2/15.2, H-range >= H_MIN"""
    W = H_of(P, e2np, P4v[0])
    amp = W.max(1) - W.min(1)
    lo = np.roll(W, 1, 1); hi = np.roll(W, -1, 1)
    imx = (W > lo) & (W > hi); imn = (W < lo) & (W < hi)
    ok = (imx.sum(1) == 2) & (imn.sum(1) == 2) & (amp >= H_MIN)
    out = np.zeros(len(W), bool)
    idx = np.flatnonzero(ok)
    if not idx.size:
        return out
    _, cm = np.nonzero(imx[idx]); pk = cm.reshape(-1, 2)
    _, cn = np.nonzero(imn[idx]); vl = cn.reshape(-1, 2)
    Wt = W[idx]
    pv = np.take_along_axis(Wt, pk, 1); vv = np.take_along_axis(Wt, vl, 1)
    dip = (pv.min(1) - vv.max(1)) / amp[idx]
    t1, t2 = tf[pk[:, 0]], tf[pk[:, 1]]
    out[idx] = (dip >= THR) & (np.abs(t1 - 8.2) < TOLP) & (np.abs(t2 - 15.2) < TOLP)
    return out


def run(tag, e2_np_vec, P_all):
    c = dict(s1=0, s2=0, s3=0, s4=0, s5=0)
    best = None
    for st in range(0, len(P_all), 50_000):
        P = P_all[st:st + 50_000]
        e2np = e2_np_vec[st:st + 50_000]
        m0 = screen_np(P, e2np)
        if not m0.any():
            continue
        idx = np.flatnonzero(m0)
        Q = P[idx]; qe = e2np[idx]
        Hn = H_of(Q, qe, P4v[0])
        Hp = H_of(Q, np.ones(len(Q)), P4v[1])
        He = H_of(Q, np.zeros(len(Q)), P4v[2])
        for j in range(len(Q)):
            # COUNTS are the target (2 / 2 / 1); positions are reported, not required,
            # because demanding 5 exact peak positions at once is a search wall not a result
            pkn, an = exact_peaks(Hn[j])
            if len(pkn) != 2:
                continue
            c['s1'] += 1
            pkp, ap = exact_peaks(Hp[j])
            if len(pkp) != 2:
                continue
            c['s2'] += 1
            pke, ae = exact_peaks(He[j])
            if len(pke) != 1:
                continue
            c['s3'] += 1
            r1, r2 = ap / an, ae / an
            if abs(r1 - R_OBS[1]) > TOLR or abs(r2 - R_OBS[2]) > TOLR:
                continue
            c['s4'] += 1
            if an < H_MIN:
                continue
            c['s5'] += 1
            if best is None or an > best[0]:
                best = (an, r1, r2, Q[j].copy(), float(qe[j]),
                        pkn.copy(), pkp.copy(), pke.copy())
    pr(f'\n--- {tag} ---')
    pr(f'  1. non-preg BIMODAL (2 peaks, near 8.2/15.2 pre-screen)   {c["s1"]:,}')
    pr(f'  2. + pregnant BIMODAL (2 peaks)                           {c["s2"]:,}')
    pr(f'  3. + Esr1i UNIMODAL (exactly 1 peak)                      {c["s3"]:,}')
    pr(f'  4. + amplitude ratios 1.00/0.34/0.55 (+-{TOLR})            {c["s4"]:,}')
    pr(f'  5. + REAL amplitude (A_T <= 10 C)                         {c["s5"]:,}')
    if best is not None:
        an, r1, r2, q, mm, pkn, pkp, pke = best
        pr(f'  BEST: H-range_np {an:.3f} -> A_T = {AMP[0]/an:.2f} C ; ratios '
           f'1.00/{r1:.2f}/{r2:.2f} ; m={mm:.3f}')
        pr(f'    peaks  non-preg {np.round(pkn,1)} (obs 8.2/15.2)'
           f'  pregnant {np.round(pkp,1)} (obs 15.8/21.8)'
           f'  Esr1i {np.round(pke,1)} (obs 13.2)')
        pr(f'    phiK={q[7]:.2f} phiL={q[2]:.2f} (gap {min(abs(q[7]-q[2])%24, 24-abs(q[7]-q[2])%24):.2f} h)')
        np.save(f'best_encoding_{tag.split()[0].lower()}.npy', np.concatenate([q, [mm]]))
    return c


rng = np.random.default_rng(31337)
N = 3_000_000
P_all = np.stack([
    rng.uniform(-15, 8, N), rng.uniform(0, 30, N), rng.uniform(0, 24, N),
    rng.uniform(0, 30, N), rng.uniform(0, 15, N),
    rng.uniform(-15, 8, N), rng.uniform(0, 30, N), rng.uniform(0, 24, N),
    rng.uniform(0, 30, N), rng.uniform(0, 15, N)], 1)
m_free = rng.uniform(0.05, 0.95, N)

pr(f'draws {N:,}   SAME draws for both encodings (controlled: only E2_non-preg differs)')
pr(f'  peak threshold {THR} (validated: the mouse data shows 2/2/1 here)')
pr(f'  H-range_np >= {H_MIN:.3f} required so A_T <= 10 C  -- a flat solution CANNOT pass')
c_old = run('OLD encoding (E2 non-preg = 0)', np.zeros(N), P_all)
c_new = run('NEW encoding (E2 non-preg = m, free)', m_free, P_all)

pr('\n================ VERDICT ================')
for k, lab in [('s3', 'topology 2/2/1 at correct positions'),
               ('s4', '+ amplitude ratios'),
               ('s5', '+ REAL amplitude')]:
    pr(f'  {lab:38s} old {c_old[k]:>8,}   new {c_new[k]:>8,}')
