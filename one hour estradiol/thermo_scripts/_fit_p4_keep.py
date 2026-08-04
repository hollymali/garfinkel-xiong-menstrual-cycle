"""P4 as a KEEP-arm-only term (cutaneous vasoconstriction / heat retention).

Holly 08-03: P4 is NOT part of the heat-loss arm.  It drives vasoconstriction and
raises heat, so it enters the heat-KEEP arm only:

    c(t) = cos(2*pi*(t-12)/24)                     one noon-locked clock
    loss = sig(bL + gL*c + wL*E2)                  <- no P4
    keep = sig(bK + gK*c + wK*E2 + pK*P4)          <- P4 here only
    CBT  = T0 + A_T*(keep - loss) + drift*(day-15.5)

    E2 = 0 / 1 / 0     non-preg / pregnant / Esr1i   (Esr1i transduces no estrogen)
    P4 = 0 / line / line                            (both pregnant animals)

All coefficients positive.  Question: does dropping P4 out of the loss arm let
Esr1i go UNIMODAL while the two ERalpha-intact animals stay bimodal?
Run with pL free as the control comparison.
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
E2 = np.array([0., 1., 0.])[:, None, None]
PREG = np.array([0., 1., 1.])
DAYS = [13, 14, 15, 16, 17]
mid = np.array(DAYS, float) + 0.5
P4M = np.outer(PREG, np.maximum(0., 1. - (mid - 13.) / 6.5))[:, :, None]
dd = (mid - 15.5)[None, :, None]
cb = c[None, None, :]


def sig(z):
    return 1 / (1 + np.exp(-np.clip(z, -60, 60)))


def waves(z, p4_in_loss):
    bL, gL, wL, bK, gK, wK, pK, A_T, T0, dr = z[:10]
    pL = z[10] if p4_in_loss else 0.0
    return (T0
            + A_T * (sig(bK + gK * cb + wK * E2 + pK * P4M)
                     - sig(bL + gL * cb + wL * E2 + pL * P4M))
            + dr * dd)


def obj(z, p4_in_loss):
    W = waves(z, p4_in_loss)
    Wm = W.mean(axis=1)
    r = Wm - OBS
    e = float(np.mean(np.sqrt(np.mean(r ** 2, 1)) / AMP))          # waveform
    fl = float(np.mean(((Wm.std(1) - OBS.std(1)) / OBS.std(1)) ** 2))  # anti-flatten
    edge = sum(1.0 for i in range(3)
               for k in (int(np.argmax(Wm[i])), int(np.argmin(Wm[i])))
               if k < 2 or k > NB - 3)                              # boundary cliffs
    dm = W.mean(axis=2)
    de = float(np.sqrt(np.mean(((dm - dm.mean(1, keepdims=True))
                                - (omean - omean.mean(1, keepdims=True))) ** 2)))
    return e + 2.0 * fl + 1.0 * edge + 1.5 * de


B0 = [(-30, 10), (0, 60), (0, 60), (-30, 10), (0, 60), (0, 60),
      (0, 25), (0.3, 8), (30, 42), (-0.4, 0.1)]


def peaks_of(w, prom_frac):
    a = np.ptp(w)
    ext = np.concatenate([w] * 3)
    pk, _ = find_peaks(ext, prominence=prom_frac * max(a, 1e-9))
    return [round(float(hrs[p - NB]), 1) for p in pk if NB <= p < 2 * NB]


def run(p4_in_loss, seed, ntry=400):
    B = B0 + ([(0, 25)] if p4_in_loss else [])
    lo = np.array([q[0] for q in B]); hi = np.array([q[1] for q in B])
    rng = np.random.default_rng(seed)
    best = None
    for _ in range(ntry):
        z0 = lo + rng.random(len(B)) * (hi - lo)
        r = minimize(obj, z0, args=(p4_in_loss,), bounds=B,
                     method='L-BFGS-B', options=dict(maxiter=500, ftol=1e-13))
        if best is None or r.fun < best[0]:
            best = (float(r.fun), r.x)
    return best


def report(tag, e, z, p4_in_loss):
    W = waves(z, p4_in_loss); Wm = W.mean(axis=1); dm = W.mean(axis=2)
    bL, gL, wL, bK, gK, wK, pK, A_T, T0, dr = z[:10]
    pL = z[10] if p4_in_loss else 0.0
    print(f'\n===== {tag}   loss={e:.4f} =====')
    print(f'  LOSS  bL={bL:+.2f} gL={gL:.2f} wL={wL:.2f} pL={pL:.2f}')
    print(f'  KEEP  bK={bK:+.2f} gK={gK:.2f} wK={wK:.2f} pK={pK:.2f}')
    print(f'  A_T={A_T:.2f}  T0={T0:.2f}  drift={dr:+.3f} C/day')
    print(f'  {"animal":10s}{"amp mod":>9}{"amp obs":>9}{"SD mod":>8}{"SD obs":>8}'
          f'{"drift mod":>11}{"drift obs":>10}')
    for i, l in enumerate(L):
        print(f'  {l:10s}{np.ptp(Wm[i]):9.2f}{AMP[i]:9.2f}{Wm[i].std():8.3f}'
              f'{OBS[i].std():8.3f}{dm[i][-1]-dm[i][0]:+11.2f}'
              f'{omean[i][-1]-omean[i][0]:+10.2f}')
    print('  PEAKS across prominence thresholds (fraction of amplitude):')
    for i, l in enumerate(L):
        row = '  '.join(f'{p:.2f}:{peaks_of(Wm[i], p)}' for p in (0.05, 0.10, 0.15, 0.25))
        print(f'    {l:10s}{row}')
    return Wm


print('OBSERVED  peaks: non-preg [8.2, 15.2] | pregnant [15.8, 21.8] | Esr1i [13.2] ONLY')
print(f'OBSERVED  amps: {AMP[0]:.2f} {AMP[1]:.2f} {AMP[2]:.2f}   '
      f'ratios 1.00/{AMP[1]/AMP[0]:.2f}/{AMP[2]/AMP[0]:.2f}')

eA, zA = run(False, 2026)
WA = report('P4 in KEEP arm only  (Holly 08-03)', eA, zA, False)
np.save('best_p4_keep_only.npy', zA)

eB, zB = run(True, 2026)
WB = report('control: P4 in BOTH arms (pL free)', eB, zB, True)

# ---- structural scan: is unimodal Esr1i reachable at all, at any pK? ----
print('\n===== pK scan (other params held at the keep-only fit) =====')
print(f'  {"pK":>6}{"npreg pk":>22}{"preg pk":>22}{"Esr1i pk":>22}{"E amp":>8}')
z = zA.copy()
for pK in [0, 0.25, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 20]:
    z[6] = pK
    Wm = waves(z, False).mean(axis=1)
    print(f'  {pK:6.2f}{str(peaks_of(Wm[0],0.10)):>22}{str(peaks_of(Wm[1],0.10)):>22}'
          f'{str(peaks_of(Wm[2],0.10)):>22}{np.ptp(Wm[2]):8.2f}')
