"""ERalpha-PRESENCE GATE on the heat-LOSS arm's circadian drive.

THE FIX (one new parameter):
    keep = sig(bK + gK*cos(2pi(t-phiK)/24) + wK*E2 + pK*P4)
    loss = sig(bL + gL*R*cos(2pi(t-phiL)/24) + wL*E2 + pL*P4)
    R    = 1 (non-preg)   1 (pregnant)   rho (Esr1i)
E2 = 0/1/0 and P4 = 0/p4bar/p4bar as before; the ONLY change is R.

WHY (grounded, not invented): the direct E2 -> MnPO -> MPA -> RPa HEAT-LOSS circuit is already
the established pathway in this project ([[ref-estrogen-thermoregulation-review]]).  Gating the
LOSS arm on receptor PRESENCE says: without ERalpha that limb loses its circadian drive, so the
loss arm stops carving a notch into the keep plateau and the day collapses to one hump.  This is
the "receptor-presence-gated term" already listed as one of the three escapes from impossibility
#3, and the one never tried.  It also explains why Esr1i differs from a merely LOW-E2 animal:
losing the receptor removes the pathway; low hormone only turns it down.

WHY IT SHOULD WORK (measured, from the per-day raw traces):
    amplitude   non-preg 2.45 > Esr1i 1.23 > pregnant 0.82  -> monotone in hormone LOAD, Esr1i MIDDLE
    dip/amp     non-preg 0.359, pregnant 0.321, Esr1i 0.088 -> tracks ERalpha FUNCTION, Esr1i EXTREME
Those are different quantities, so Esr1i can be middle in one and extreme in the other with no
contradiction.  The old "bracket" only bit because both were forced through a single bias axis.

TWO ANTI-TRAP MEASURES (both failure modes actually hit earlier in this work):
  * dip/amp is fitted as a CONTINUOUS target, so there is no binary label to park on.  Four
    successive hand-rolled peak classifiers were gamed by the optimizer; a continuous quantity
    removes that whole class of bug.
  * A_T is FIXED BY CONSTRUCTION as AMP_obs[0]/H-range(non-preg), so non-preg's amplitude is
    exactly 2.10 C and a near-flat solution is impossible (it would need A_T > 10 C and is
    rejected).  Previously a least-squares A_T collapsed to 0.12 C and "solved" the topology
    on a flat line.
"""
import numpy as np
from scipy.signal import find_peaks
from scipy.optimize import minimize

fold = np.load('fold.npy', allow_pickle=True).item()
L = ['non-preg', 'pregnant', 'Esr1i']
hrs = fold[L[0]][0]
OBS = np.vstack([fold[l][1] for l in L])
AMP = OBS.max(1) - OBS.min(1)
R_OBS = AMP / AMP[0]

tf = np.arange(0, 24, 0.2)
NT = tf.size
E2 = np.array([0., 1., 0.])
mid = np.array([13, 14, 15, 16, 17], float) + 0.5
p4bar = np.maximum(0., 1. - (mid - 13.) / 6.5).mean()
P4v = np.array([0., p4bar, p4bar])
RA = None                      # filled per-parameter: [1, 1, rho]

DIP_OBS = np.array([0.359, 0.321, 0.088])     # measured per-day, the continuous target
A_T_MAX = 10.0
H_MIN = AMP[0] / A_T_MAX                      # 0.21


def pr(*a):
    print(*a, flush=True)


def sig(z):
    return 1 / (1 + np.exp(-np.clip(z, -60, 60)))


def H_one(p, i, t=None):
    bL, gL, phiL, wL, pL, bK, gK, phiK, wK, pK, rho = p
    R = rho if i == 2 else 1.0
    t = tf if t is None else t
    cK = np.cos(2 * np.pi * (t - phiK) / 24.)
    cL = np.cos(2 * np.pi * (t - phiL) / 24.)
    return (sig(bK + gK * cK + wK * E2[i] + pK * P4v[i])
            - sig(bL + gL * R * cL + wL * E2[i] + pL * P4v[i]))


def dip_of(w):
    """continuous notch depth: prominence of the 2nd-most-prominent peak, / amplitude"""
    amp = w.max() - w.min()
    if amp < 1e-12:
        return 0.0, amp
    ext = np.concatenate([w, w, w])
    pk, props = find_peaks(ext, prominence=1e-9)
    m = (pk >= NT) & (pk < 2 * NT)
    pv = np.sort(props['prominences'][m])[::-1]
    return (float(pv[1] / amp) if len(pv) >= 2 else 0.0), amp


def evaluate(p):
    Hs, dips, amps = [], [], []
    for i in range(3):
        h = H_one(p, i)
        d, a = dip_of(h)
        Hs.append(h); dips.append(d); amps.append(a)
    return np.array(Hs), np.array(dips), np.array(amps)


def score(p):
    if np.any(np.asarray(p)[[1, 3, 4, 6, 8, 9]] < 0) or not (0.0 <= p[10] <= 1.0):
        return 1e6
    H, dip, amp = evaluate(p)
    if amp[0] < H_MIN:
        return 1e6                       # A_T would exceed 10 C -> flat/unphysical
    A_T = AMP[0] / amp[0]                # non-preg amplitude EXACT by construction
    Hd = np.stack([H_one(p, i, hrs) for i in range(3)])   # on the DATA grid for RMS
    W = A_T * Hd
    W = W + (np.mean(OBS) - np.mean(W))
    r = amp / amp[0]
    e_dip = float(np.sum(np.abs(dip - DIP_OBS)))
    e_rat = float(abs(r[1] - R_OBS[1]) + abs(r[2] - R_OBS[2]))
    e_rms = float(np.mean(np.sqrt(np.mean((W - OBS) ** 2, 1)) / AMP))
    return 6.0 * e_dip + 4.0 * e_rat + e_rms


# ---------------------------------------------------------------- seed scan
def screen(P):
    """vectorised: non-preg must have a real notch and enough H-range"""
    bL, gL, phiL, wL, pL, bK, gK, phiK, wK, pK = [P[:, i:i + 1] for i in range(10)]
    cK = np.cos(2 * np.pi * (tf[None, :] - phiK) / 24.)
    cL = np.cos(2 * np.pi * (tf[None, :] - phiL) / 24.)
    W = sig(bK + gK * cK) - sig(bL + gL * cL)
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
    d = (pv.min(1) - vv.max(1)) / amp[idx]
    out[idx] = (d > 0.15) & (d < 0.60)
    return out


rng = np.random.default_rng(777)
N, CH = 3_000_000, 50_000
seeds = []
for s in range(0, N, CH):
    P = np.stack([
        rng.uniform(-15, 8, CH), rng.uniform(0, 30, CH), rng.uniform(0, 24, CH),
        rng.uniform(0, 30, CH), rng.uniform(0, 15, CH),
        rng.uniform(-15, 8, CH), rng.uniform(0, 30, CH), rng.uniform(0, 24, CH),
        rng.uniform(0, 30, CH), rng.uniform(0, 15, CH)], 1)
    m = screen(P)
    if m.any():
        Q = P[m]
        seeds.append(np.column_stack([Q, rng.uniform(0, 0.6, len(Q))]))   # rho
seeds = np.concatenate(seeds) if seeds else np.zeros((0, 11))
pr(f'draws {N:,}  seeds with a real non-preg notch and H-range >= {H_MIN:.2f}: {len(seeds):,}')
if not len(seeds):
    raise SystemExit

sub = seeds[:8000]
vals = np.array([score(s) for s in sub])
alive = np.flatnonzero(vals < 1e5)
pr(f'  scored {len(sub):,}, finite {len(alive):,}, best seed {vals[alive].min():.4f}')
best = (float(vals[alive].min()), sub[alive[np.argmin(vals[alive])]].copy())
for k, i in enumerate(alive[np.argsort(vals[alive])][:150]):
    z = sub[i].copy()
    for _ in range(3):
        r = minimize(score, z, method='Nelder-Mead',
                     options=dict(maxiter=5000, maxfev=5000, xatol=1e-9, fatol=1e-12))
        z = r.x
    if float(r.fun) < best[0]:
        best = (float(r.fun), z.copy())
        np.save('_era_gate_ckpt.npy', best[1])
        pr(f'    [{k+1}/150] new best {best[0]:.4f}')

v, p = best
H, dip, amp = evaluate(p)
A_T = AMP[0] / amp[0]
W = A_T * H
W = W + (np.mean(OBS) - np.mean(W))
bL, gL, phiL, wL, pL, bK, gK, phiK, wK, pK, rho = p
d = abs(phiK - phiL) % 24
pr(f'\n================ ERalpha-GATE RESULT   objective {v:.4f} ================')
pr(f'  LOSS  b={bL:+.2f} g={gL:.2f} phi={phiL:.2f} h  w(E2)={wL:.2f} p(P4)={pL:.2f}')
pr(f'  KEEP  b={bK:+.2f} g={gK:.2f} phi={phiK:.2f} h  w(E2)={wK:.2f} p(P4)={pK:.2f}')
pr(f'  *** rho (ERalpha gate on loss-arm clock drive) = {rho:.4f} ***')
pr(f'  arm phase gap {min(d, 24-d):.2f} h    A_T={A_T:.2f} C (fixed by construction)')
pr(f'\n  {"animal":10s}{"amp mod":>9}{"amp obs":>9}{"ratio":>8}{"obsR":>7}{"dip mod":>9}{"dip obs":>9}')
for i, l in enumerate(L):
    pr(f'  {l:10s}{amp[i]*A_T:9.2f}{AMP[i]:9.2f}{amp[i]/amp[0]:8.2f}{R_OBS[i]:7.2f}'
       f'{dip[i]:9.3f}{DIP_OBS[i]:9.3f}')
# honest peak report at the data-validated threshold
pr('')
for i, l in enumerate(L):
    a = W[i].max() - W[i].min()
    ext = np.concatenate([W[i]] * 3)
    pk, _ = find_peaks(ext, prominence=0.20 * a)
    m = (pk >= NT) & (pk < 2 * NT)
    obs = ['8.2/15.2', '15.8/21.8', '13.2 only'][i]
    pr(f'  {l:10s} peaks@thr0.20: {np.round(tf[pk[m]-NT],1).tolist()}   observed {obs}')
np.save('best_era_gate.npy', np.concatenate([p, [A_T]]))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.0), sharex=True)
for i, l in enumerate(L):
    ax[i].plot(hrs, OBS[i], 'k-', lw=2.2, label='mouse')
    ax[i].plot(tf, W[i], 'r--', lw=2, label='ERalpha-gated loss arm')
    ax[i].set_title(f'{l}   amp {(W[i].max()-W[i].min()):.2f} / {AMP[i]:.2f} C'
                    f'   dip {dip[i]:.2f} / {DIP_OBS[i]:.2f}')
    ax[i].set_xlabel('hour of day'); ax[i].set_xticks([0, 6, 12, 18, 24])
ax[0].set_ylabel('CBT (C)'); ax[0].legend(fontsize=8)
fig.suptitle(f'ERalpha-presence gate on the heat-loss arm circadian drive (rho={rho:.3f})', y=1.03)
fig.tight_layout(); fig.savefig('fig_era_gate.png', dpi=130, bbox_inches='tight')
pr('\nwrote fig_era_gate.png, best_era_gate.npy')
