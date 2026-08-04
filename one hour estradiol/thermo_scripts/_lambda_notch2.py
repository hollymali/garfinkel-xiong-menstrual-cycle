"""ONE new number: let the heat-DUMP arm have a smaller maximum effect than the heat-KEEP arm.

    CBT = T0 + A_T*(keep - lam*R*loss)
    keep = sig(bK + gK*cos(2pi(t-phi)/24) + wK*E2 + pK*P4)
    loss = sig(bL + gL*cos(2pi(t-phi)/24) + wL*E2 + pL*P4)
    R    = 1 (non-preg, pregnant)   rho (Esr1i)     <- ERalpha-presence gate, optional

WHY (Holly spotted it): in `fig_era_gate.png` pregnant has two TROUGHS, not two peaks.  Cause:
both arms are full-range sigmoids (0->1), so when the dump arm fires mid-plateau it subtracts the
ENTIRE keep signal and carves a canyon straight to baseline; the "two peaks" are only the
shoulders of that canyon.

WHY lam IS EXACTLY THE RIGHT KNOB (analytic, not fitted intuition):
with both arms on one clock, f(c) = sig(bK+gK*c) - lam*sig(bL+gL*c).  Make keep BROAD (turns on
early) and loss NARROW (turns on only near the clock peak).  Then
    trough at the clock peak : f = 1 - lam
    flanking peaks           : f ~ 1
    baseline                 : f ~ 0
so amplitude stays ~1 (FULL) while the notch depth relative to the peak is exactly **lam**.
Observed dip/amp: 0.359 non-preg, 0.321 pregnant, 0.088 Esr1i -> lam ~ 0.36 with the receptor
present, and lam*rho ~ 0.09 without it (rho ~ 0.25).  That is the whole mechanism.
This also explains why every previous "fold" bimodality was amplitude-poor: they all had lam = 1
implicitly, which forces the notch to go all the way down and the shape to be a canyon.

Only ONE shared clock phase (no two-phase machinery needed): keep broad + loss narrow at the
same phase already gives a broad plateau with a mid-plateau notch, which is the observed
non-preg shape (peaks 8.2/15.2 with the dip at ~11.7 h).

METHOD (fixing what went wrong last time): the objective is dominated by WAVEFORM RMS.  The
previous attempt weighted dips 6x and ratios 4x against RMS 1x, matched every summary statistic,
and produced curves UNCORRELATED with the data (r = -0.12/+0.19/-0.17).  Here RMS carries the
weight and correlation is printed for every animal before anything is claimed.
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
tf = np.arange(0, 24, 0.2); NT = tf.size
E2 = np.array([0., 1., 0.])
mid = np.array([13, 14, 15, 16, 17], float) + 0.5
p4b = np.maximum(0., 1. - (mid - 13.) / 6.5).mean()
P4v = np.array([0., p4b, p4b])
DIP_OBS = np.array([0.359, 0.321, 0.088])
NP = 12          # bK,gK,wK,pK, bL,gL,wL,pL, phi, lam, rho, A_T


def pr(*a):
    print(*a, flush=True)


def sig(z):
    return 1 / (1 + np.exp(-np.clip(z, -60, 60)))


def W_of(p, i, t):
    bK, gK, wK, pK, bL, gL, wL, pL, phi, lam, rho, A_T = p
    R = rho if i == 2 else 1.0
    c = np.cos(2 * np.pi * (t - phi) / 24.)
    keep = sig(bK + gK * c + wK * E2[i] + pK * P4v[i])
    loss = sig(bL + gL * c + wL * E2[i] + pL * P4v[i])
    return A_T * (keep - lam * R * loss)


def score(p):
    bK, gK, wK, pK, bL, gL, wL, pL, phi, lam, rho, A_T = p
    if gK < 0 or gL < 0 or wK < 0 or wL < 0 or pK < 0 or pL < 0:
        return 1e6
    if not (0 < lam <= 1) or not (0 <= rho <= 1) or not (0 < A_T <= 10):
        return 1e6
    W = np.stack([W_of(p, i, hrs) for i in range(3)])
    W = W + (OBS.mean() - W.mean())
    a = W.max(1) - W.min(1)
    if a[0] < 1e-6:
        return 1e6
    # WAVEFORM is the objective; amplitude only anchors the scale
    rms = float(np.mean(np.sqrt(np.mean((W - OBS) ** 2, 1)) / AMP))
    aerr = float(np.mean(np.abs(a - AMP) / AMP))
    ext=[np.concatenate([W[i]]*3) for i in range(3)]
    dips=[]
    for i in range(3):
        aa=a[i]
        pk,pm=find_peaks(ext[i],prominence=1e-9)
        mm=(pk>=48)&(pk<96)
        pv=np.sort(pm["prominences"][mm])[::-1]
        dips.append(pv[1]/aa if len(pv)>=2 else 0.0)
    derr=float(np.sum(np.abs(np.array(dips)-DIP_OBS)))
    return rms + 0.5 * aerr + 1.2 * derr


rng = np.random.default_rng(20260803)
N = 400_000
P = np.column_stack([
    rng.uniform(-12, 4, N), rng.uniform(0, 20, N), rng.uniform(0, 20, N), rng.uniform(0, 12, N),
    rng.uniform(-12, 4, N), rng.uniform(0, 30, N), rng.uniform(0, 20, N), rng.uniform(0, 12, N),
    rng.uniform(0, 24, N), rng.uniform(0.05, 1.0, N), rng.uniform(0.0, 1.0, N),
    rng.uniform(0.5, 6.0, N)])
vals = np.array([score(q) for q in P])
alive = np.flatnonzero(vals < 1e5)
pr(f'random draws {N:,}, finite {len(alive):,}, best seed {vals[alive].min():.4f}')
best = (float(vals[alive].min()), P[alive[np.argmin(vals[alive])]].copy())
for k, i in enumerate(alive[np.argsort(vals[alive])][:220]):
    z = P[i].copy()
    for _ in range(4):
        r = minimize(score, z, method='Nelder-Mead',
                     options=dict(maxiter=8000, maxfev=8000, xatol=1e-10, fatol=1e-13))
        z = r.x
    if float(r.fun) < best[0]:
        best = (float(r.fun), z.copy())
        pr(f'   [{k+1}/220] new best {best[0]:.4f}')

v, p = best
bK, gK, wK, pK, bL, gL, wL, pL, phi, lam, rho, A_T = p
W = np.stack([W_of(p, i, hrs) for i in range(3)])
W = W + (OBS.mean() - W.mean())
Wf = np.stack([W_of(p, i, tf) for i in range(3)])
Wf = Wf + (OBS.mean() - Wf.mean())
pr(f'\n================ lambda-notch model   objective {v:.4f} ================')
pr(f'  KEEP  b={bK:+.2f} g={gK:.2f}  w(E2)={wK:.2f} p(P4)={pK:.2f}')
pr(f'  LOSS  b={bL:+.2f} g={gL:.2f}  w(E2)={wL:.2f} p(P4)={pL:.2f}')
pr(f'  *** lambda (dump-arm max effect) = {lam:.4f}   rho (ERa gate) = {rho:.4f} ***')
pr(f'  phase={phi:.2f} h   A_T={A_T:.2f} C')
pr(f'\n  {"animal":10s}{"amp":>7}{"obs":>7}{"dip":>8}{"obsdip":>8}{"corr":>8}{"RMS":>8}   peaks / observed')
obsp = ['8.2/15.2', '15.8/21.8', '13.2 only']
ok = True
for i, l in enumerate(L):
    a = Wf[i].max() - Wf[i].min()
    ext = np.concatenate([Wf[i]] * 3)
    pk, prm = find_peaks(ext, prominence=1e-9)
    m = (pk >= NT) & (pk < 2 * NT)
    pv = np.sort(prm['prominences'][m])[::-1]
    dip = pv[1] / a if len(pv) >= 2 else 0.0
    pk2, _ = find_peaks(ext, prominence=0.20 * a)
    m2 = (pk2 >= NT) & (pk2 < 2 * NT)
    cc = np.corrcoef(W[i], OBS[i])[0, 1]
    rms = np.sqrt(np.mean((W[i] - OBS[i]) ** 2))
    if cc < 0.5:
        ok = False
    pr(f'  {l:10s}{a:7.2f}{AMP[i]:7.2f}{dip:8.3f}{DIP_OBS[i]:8.3f}{cc:+8.3f}{rms:8.3f}'
       f'   {np.round(tf[pk2[m2]-NT],1).tolist()} / {obsp[i]}')
pr(f'\n  waveform correlation >= 0.5 in all three: {ok}')
np.save('best_lambda_notch2.npy', p)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.0), sharex=True)
for i, l in enumerate(L):
    cc = np.corrcoef(W[i], OBS[i])[0, 1]
    ax[i].plot(hrs, OBS[i], 'k-', lw=2.2, label='mouse')
    ax[i].plot(tf, Wf[i], 'r--', lw=2, label='keep - lambda*loss')
    ax[i].set_title(f'{l}   amp {(Wf[i].max()-Wf[i].min()):.2f}/{AMP[i]:.2f} C   r={cc:+.2f}')
    ax[i].set_xlabel('hour of day'); ax[i].set_xticks([0, 6, 12, 18, 24])
ax[0].set_ylabel('CBT (C)'); ax[0].legend(fontsize=8)
fig.suptitle(f'Dump arm scaled by lambda={lam:.2f} (= the notch depth); '
             f'ERalpha gate rho={rho:.2f} in Esr1i', y=1.03)
fig.tight_layout(); fig.savefig('fig_lambda_notch2.png', dpi=130, bbox_inches='tight')
pr('wrote fig_lambda_notch2.png, best_lambda_notch2.npy')
