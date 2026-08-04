"""Is dphi ~ +7 h a BETTER optimum that the search simply missed?

The free fit chose dphi = +1.88 h and left pregnant broken (r=+0.18), which by the
pre-registered criterion rejects the E2-phase-shift hypothesis.  But the geometry says
otherwise: pregnant's fitted peaks (4.6/22.8 h) already have the RIGHT GAP (5.8 h vs observed
6.0 h) and are merely centred 6.9 h away from the observed midpoint of 18.8 h.  So dphi ~ +6.9
should align pregnant, and the free fit probably sat in a local optimum.

Test: PIN dphi on a grid, re-optimise everything else from the best known solution, and report
the objective and pregnant's correlation at each.  If the curve dips near +7 h, the mechanism is
real and the earlier rejection was a search artifact.  If it is flat or worst there, the
rejection stands.
"""
import numpy as np
from scipy.signal import find_peaks
from scipy.optimize import minimize

fold = np.load('fold.npy', allow_pickle=True).item()
L = ['non-preg', 'pregnant', 'Esr1i']
hrs = fold[L[0]][0]
OBS = np.vstack([fold[l][1] for l in L])
AMP = OBS.max(1) - OBS.min(1)
E2 = np.array([0., 1., 0.])
mid = np.array([13, 14, 15, 16, 17], float) + 0.5
p4b = np.maximum(0., 1. - (mid - 13.) / 6.5).mean()
P4v = np.array([0., p4b, p4b])
DIP_OBS = np.array([0.359, 0.321, 0.088])
p0 = np.load('best_lambda_phase.npy')


def sig(z):
    return 1 / (1 + np.exp(-np.clip(z, -60, 60)))


def W_of(p, i, t, dphi):
    bK, gK, wK, pK, bL, gL, wL, pL, phi, lam, rho, A_T = p[:12]
    R = rho if i == 2 else 1.0
    c = np.cos(2 * np.pi * (t - (phi + dphi * E2[i])) / 24.)
    return A_T * (sig(bK + gK * c + wK * E2[i] + pK * P4v[i])
                  - lam * R * sig(bL + gL * c + wL * E2[i] + pL * P4v[i]))


def make(dphi):
    def score(p):
        bK, gK, wK, pK, bL, gL, wL, pL, phi, lam, rho, A_T = p[:12]
        if min(gK, gL, wK, wL, pK, pL) < 0 or not (0 < lam <= 1) \
           or not (0 <= rho <= 1) or not (0 < A_T <= 10):
            return 1e6
        W = np.stack([W_of(p, i, hrs, dphi) for i in range(3)])
        W = W + (OBS.mean() - W.mean())
        a = W.max(1) - W.min(1)
        if a[0] < 1e-6:
            return 1e6
        rms = float(np.mean(np.sqrt(np.mean((W - OBS) ** 2, 1)) / AMP))
        aerr = float(np.mean(np.abs(a - AMP) / AMP))
        dips = []
        for i in range(3):
            ext = np.concatenate([W[i]] * 3)
            pk, pm = find_peaks(ext, prominence=1e-9)
            mm = (pk >= 48) & (pk < 96)
            pv = np.sort(pm['prominences'][mm])[::-1]
            dips.append(pv[1] / a[i] if len(pv) >= 2 else 0.0)
        derr = float(np.sum(np.abs(np.array(dips) - DIP_OBS)))
        return rms + 0.5 * aerr + 1.2 * derr
    return score


print(f'{"dphi":>7}{"objective":>12}{"r non-preg":>12}{"r pregnant":>12}{"r Esr1i":>10}'
      f'{"preg peaks":>22}', flush=True)
rows = []
for dphi in np.arange(-2, 12.1, 1.0):
    f = make(dphi)
    z = p0[:12].copy()
    best = (f(z), z.copy())
    for _ in range(5):
        r = minimize(f, z, method='Nelder-Mead',
                     options=dict(maxiter=6000, maxfev=6000, xatol=1e-10, fatol=1e-13))
        z = r.x
        if r.fun < best[0]:
            best = (float(r.fun), z.copy())
    v, z = best
    W = np.stack([W_of(z, i, hrs, dphi) for i in range(3)])
    W = W + (OBS.mean() - W.mean())
    cs = [np.corrcoef(W[i], OBS[i])[0, 1] for i in range(3)]
    a1 = W[1].max() - W[1].min()
    ext = np.concatenate([W[1]] * 3)
    pk, _ = find_peaks(ext, prominence=0.20 * a1)
    mm = (pk >= 48) & (pk < 96)
    pp = np.round(hrs[pk[mm] - 48], 1).tolist()
    rows.append((dphi, v, cs, z))
    print(f'{dphi:7.1f}{v:12.4f}{cs[0]:+12.3f}{cs[1]:+12.3f}{cs[2]:+10.3f}{str(pp):>22}',
          flush=True)
best = min(rows, key=lambda r: r[1])
print(f'\nBEST dphi = {best[0]:+.1f} h   objective {best[1]:.4f}   '
      f'r = {best[2][0]:+.2f}/{best[2][1]:+.2f}/{best[2][2]:+.2f}', flush=True)
print('observed pregnant peaks 15.8/21.8', flush=True)
np.save('best_dphi_scan.npy', np.concatenate([best[3][:12], [best[0]]]))
