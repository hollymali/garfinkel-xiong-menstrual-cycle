"""
Separate the two failures.

Every model is locked to a single clock phase, but the observed peak moves
8.2 -> 15.8 -> 13.2 h. A correctly-sized rhythm at the wrong phase scores worse
than a flat line, so the fitter kills the amplitude. That means the earlier
"all mechanisms over-compress" result confounds AMPLITUDE with PHASE.

Give each condition its own phase (phase is not an amplitude knob) and ask the
amplitude question cleanly: can the mechanism reproduce ratios 1.00/0.34/0.55?
"""
import numpy as np
from scipy.optimize import minimize

fold = np.load('fold.npy', allow_pickle=True).item()
LABELS = ['non-preg', 'pregnant', 'Esr1i']
hrs = fold[LABELS[0]][0]
OBS = np.vstack([fold[l][1] for l in LABELS])
NB = OBS.shape[1]
AMP = OBS.max(axis=1) - OBS.min(axis=1)
w_ang = 2 * np.pi * np.fft.fftfreq(NB, d=1440.0 / NB)
tgrid = hrs * 60.0

def sig(x): return 1.0 / (1.0 + np.exp(-np.clip(x, -60, 60)))
def plant(H, A_T, tau, T0):
    Tk = -A_T * np.fft.fft(H, axis=1) / (1.0 + 1j * w_ang * tau)[None, :]
    return np.real(np.fft.ifft(Tk, axis=1)) + T0
def clocks(ph):
    return np.cos(2*np.pi*(tgrid[None, :] - np.asarray(ph)[:, None]*60.0)/1440.0)

def H_dual(p, E2, P4, c):
    b, g, w, pS, pI = p
    return sig(b+g*c+w*E2+pS*P4) - sig(b-g*c+w*E2+pI*P4)
def H_single(p, E2, P4, c):
    b, g, w, pS, pI = p
    return sig(b+g*c+w*E2+pS*P4) - sig(b-g*c+pI*P4)
def H_amsh(p, E2, P4, c):
    R, Q, x0, wE, dP = p
    return sig(R*(1.0-Q*E2)*c + x0 + wE*E2 + dP*P4)
def H_gain(p, E2, P4, c):
    A, Q, _u, wE, dP = p
    return A*(1.0-Q*E2)*c + wE*E2 + dP*P4

WIDE = [(-8, 8), (0, 40), (-25, 25), (-8, 8), (-8, 8)]
MODELS = {
    'dual       ': dict(H=H_dual,   bounds=WIDE),
    'dual-forced': dict(H=H_dual,   bounds=[(-2, 2), (0, 40), (0, 25), (-8, 8), (-8, 8)]),
    'single     ': dict(H=H_single, bounds=WIDE),
    'amsh       ': dict(H=H_amsh,   bounds=[(0, 40), (-2, 2), (-8, 8), (-8, 8), (-8, 8)]),
    'linear-gain': dict(H=H_gain,   bounds=[(0, 10), (-2, 2), (0, 0), (-8, 8), (-8, 8)]),
}
TAIL = [(0.1, 15), (1, 1500), (30, 42), (0.0, 1.0)] + [(0.0, 24.0)]*3

def model(spec, x, idx):
    core = x[:5]; A_T, tau, T0, rho = x[5:9]; ph = x[9:12]
    E2 = np.array([0.0, 1.0, rho])[idx][:, None]
    P4 = np.array([0.0, 1.0, 1.0])[idx][:, None]
    return plant(spec['H'](core, E2, P4, clocks(ph)[idx]), A_T, tau, T0)

def loss(spec, x, idx):
    r = model(spec, x, idx) - OBS[idx]
    return float(np.sqrt(np.mean((np.sqrt(np.mean(r**2, axis=1))/AMP[idx])**2)))

ALL = np.array([0, 1, 2])
rng = np.random.default_rng(31)
print(f"observed ratios {np.round(AMP/AMP[0],2)}   peak hrs "
      f"{np.round(hrs[OBS.argmax(axis=1)],1)}\n")
print("PER-CONDITION FREE PHASE (amplitude question isolated)")
print("=" * 70)
for name, spec in MODELS.items():
    bnds = spec['bounds'] + TAIL
    lo = np.array([b[0] for b in bnds]); hi = np.array([b[1] for b in bnds])
    best = None
    for _ in range(200):
        x0 = lo + rng.random(len(bnds))*(hi-lo)
        r = minimize(lambda z: loss(spec, z, ALL), x0, bounds=bnds,
                     method='L-BFGS-B', options=dict(maxiter=300, ftol=1e-14))
        if best is None or r.fun < best[0]:
            best = (float(r.fun), r.x)
    e, x = best
    M = model(spec, x, ALL)
    a = M.max(axis=1) - M.min(axis=1)
    print(f"{name}  loss {e:.4f}   amps {np.round(a,2)}  "
          f"ratios {np.round(a/a[0],2)}   rho={x[8]:.2f}  "
          f"phases {np.round(x[9:12],1)}")
