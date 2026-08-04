"""
Final comparison, with the objective fixed.

Pooled RMSE in degrees C is dominated by the non-pregnant condition (amp 2.10,
and its waveform is strongly non-sinusoidal so it carries a large irreducible
residual). Under that loss the optimiser happily throws away the pregnant
amplitude, which is exactly the "over-compression" seen earlier -- an artefact
of the loss, not of the mechanism.

Fix: score each condition RELATIVE to its own observed amplitude, so a 0.1 C
error on the pregnant rhythm (amp 0.72) costs more than the same error on the
non-pregnant one (amp 2.10).
"""
import numpy as np
from scipy.optimize import minimize

fold = np.load('fold.npy', allow_pickle=True).item()
LABELS = ['non-preg', 'pregnant', 'Esr1i']
hrs = fold[LABELS[0]][0]
OBS = np.vstack([fold[l][1] for l in LABELS])
SD = np.vstack([fold[l][2] for l in LABELS])
NB = OBS.shape[1]
AMP = OBS.max(axis=1) - OBS.min(axis=1)
print(f"observed amps {np.round(AMP,2)}  ratios {np.round(AMP/AMP[0],2)}")

w_ang = 2 * np.pi * np.fft.fftfreq(NB, d=1440.0 / NB)
tgrid = hrs * 60.0

def sig(x):
    return 1.0 / (1.0 + np.exp(-np.clip(x, -60, 60)))

def plant(H, A_T, tau, T0):
    Tk = -A_T * np.fft.fft(H, axis=1) / (1.0 + 1j * w_ang * tau)[None, :]
    return np.real(np.fft.ifft(Tk, axis=1)) + T0

def clocks(phi):
    return np.cos(2 * np.pi * (tgrid - phi * 60.0) / 1440.0)[None, :]

def H_dual(p, E2, P4, c):
    b, g, w, pS, pI = p
    return sig(b + g*c + w*E2 + pS*P4) - sig(b - g*c + w*E2 + pI*P4)

def H_single(p, E2, P4, c):
    b, g, w, pS, pI = p
    return sig(b + g*c + w*E2 + pS*P4) - sig(b - g*c + pI*P4)

def H_amsh(p, E2, P4, c):
    R, Q, x0, wE, dP = p
    return sig(R*(1.0 - Q*E2)*c + x0 + wE*E2 + dP*P4)

def H_gain(p, E2, P4, c):
    """no sigmoid at all: pure gain scaling + additive offsets"""
    A, Q, _u1, wE, dP = p
    return A*(1.0 - Q*E2)*c + wE*E2 + dP*P4

WIDE = [(-8, 8), (0, 40), (-25, 25), (-8, 8), (-8, 8)]
MODELS = {
    'dual       ': dict(H=H_dual,   bounds=WIDE),
    'dual-forced': dict(H=H_dual,   bounds=[(-2, 2), (0, 40), (0, 25), (-8, 8), (-8, 8)]),
    'single     ': dict(H=H_single, bounds=WIDE),
    'amsh       ': dict(H=H_amsh,   bounds=[(0, 40), (-2, 2), (-8, 8), (-8, 8), (-8, 8)]),
    'linear-gain': dict(H=H_gain,   bounds=[(0, 10), (-2, 2), (0, 0), (-8, 8), (-8, 8)]),
}
TAIL = [(0.1, 15), (1, 1500), (30, 42), (0.0, 1.0), (0.0, 24.0)]

def model(spec, x, idx):
    core = x[:5]; A_T, tau, T0, rho, phi = x[5:10]
    E2 = np.array([0.0, 1.0, rho])[idx][:, None]
    P4 = np.array([0.0, 1.0, 1.0])[idx][:, None]
    return plant(spec['H'](core, E2, P4, clocks(phi)), A_T, tau, T0)

def loss(spec, x, idx):
    """RMS across conditions of each condition's RMSE / its observed amplitude"""
    r = model(spec, x, idx) - OBS[idx]
    per = np.sqrt(np.mean(r ** 2, axis=1)) / AMP[idx]
    return float(np.sqrt(np.mean(per ** 2)))

ALL = np.array([0, 1, 2]); TR = np.array([0, 1]); TE = np.array([2])
rng = np.random.default_rng(23)

def fit(spec, idx, n):
    bnds = spec['bounds'] + TAIL
    lo = np.array([b[0] for b in bnds]); hi = np.array([b[1] for b in bnds])
    out = []
    for _ in range(n):
        x0 = lo + rng.random(len(bnds)) * (hi - lo)
        r = minimize(lambda z: loss(spec, z, idx), x0, bounds=bnds,
                     method='L-BFGS-B', options=dict(maxiter=500, ftol=1e-15))
        out.append((float(r.fun), r.x))
    return sorted(out, key=lambda s: s[0])

print("\n" + "=" * 74)
print("IN-SAMPLE (normalised loss; lower = better)")
print("=" * 74)
best = {}
for name, spec in MODELS.items():
    e, x = fit(spec, ALL, 400)[0]
    best[name] = (e, x)
    M = model(spec, x, ALL)
    a = M.max(axis=1) - M.min(axis=1)
    pk = hrs[M.argmax(axis=1)]
    print(f"{name}  loss {e:.4f}   amps {np.round(a,2)}  ratios {np.round(a/a[0],2)}"
          f"   peak hrs {np.round(pk,1)}")
print(f"{'observed':11s}  {'':11s}  amps {np.round(AMP,2)}  "
      f"ratios {np.round(AMP/AMP[0],2)}   peak hrs "
      f"{np.round(hrs[OBS.argmax(axis=1)],1)}")

print("\n" + "=" * 74)
print("OUT-OF-SAMPLE: train non-preg+pregnant, predict Esr1i (rho only free)")
print("=" * 74)
for name, spec in MODELS.items():
    sols = fit(spec, TR, 300)
    thr = sols[0][0] * 1.10
    keep = [s for s in sols if s[0] <= thr]
    tests = np.array([min(loss(spec, np.r_[x[:8], r, x[9]], TE)
                          for r in np.linspace(0, 1, 101)) for _, x in keep])
    print(f"{name}  train {sols[0][0]:.4f} ({len(keep)} kept)   held-out Esr1i: "
          f"best {tests.min():.4f}  median {np.median(tests):.4f}")

np.save('best_final.npy', np.array(best, dtype=object), allow_pickle=True)
