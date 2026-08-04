"""
Definitive test on the digitised 24 h waveforms (144 points vs ~10 params).

Variants:
  dual        E2 pushes BOTH arms toward saturation, biases unconstrained
  dual-forced same, but |b| <= 2 and w >= 0, i.e. both arms are FORCED to sit
              in their responsive range so the two-sided mechanism is really
              used. (Unconstrained, the optimiser rails b to +-8 and effectively
              throws one arm away -- that is not a test of the hypothesis.)
  single      E2 reaches the stimulatory arm only
  amsh        E2 shrinks clock amplitude upstream of one sigmoid
"""
import numpy as np
from scipy.optimize import minimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fold = np.load('fold.npy', allow_pickle=True).item()
LABELS = ['non-preg', 'pregnant', 'Esr1i']
hrs = fold[LABELS[0]][0]
OBS = np.vstack([fold[l][1] for l in LABELS])
SD  = np.vstack([fold[l][2] for l in LABELS])
NB = OBS.shape[1]
NDAYS = 5
FLOOR = np.nanmean(SD) / np.sqrt(NDAYS)      # rough SEM noise floor

print("OBSERVED (digitised, folded to 24 h)")
for i, l in enumerate(LABELS):
    a = OBS[i].max() - OBS[i].min()
    print(f"  {l:9s} amp {a:.2f}  ratio {a/(OBS[0].max()-OBS[0].min()):.2f}  "
          f"mean {OBS[i].mean():.2f}  peak hour {hrs[OBS[i].argmax()]:.1f}")
print(f"  approx noise floor (mean bin SEM over {NDAYS} days): {FLOOR:.3f} C\n")

w_ang = 2 * np.pi * np.fft.fftfreq(NB, d=1440.0 / NB)
tgrid = hrs * 60.0

def sig(x):
    return 1.0 / (1.0 + np.exp(-np.clip(x, -60, 60)))

def plant(H, A_T, tau, T0):
    Tk = -A_T * np.fft.fft(H, axis=1) / (1.0 + 1j * w_ang * tau)[None, :]
    return np.real(np.fft.ifft(Tk, axis=1)) + T0

def clocks(phis):
    return np.cos(2 * np.pi * (tgrid[None, :] - phis[:, None] * 60.0) / 1440.0)

def H_dual(p, E2, P4, c):
    b, g, w, pS, pI = p
    return sig(b + g*c + w*E2 + pS*P4) - sig(b - g*c + w*E2 + pI*P4)

def H_single(p, E2, P4, c):
    b, g, w, pS, pI = p
    return sig(b + g*c + w*E2 + pS*P4) - sig(b - g*c + pI*P4)

def H_amsh(p, E2, P4, c):
    R, Q, x0, wE, dP = p
    return sig(R*(1.0 - Q*E2)*c + x0 + wE*E2 + dP*P4)

WIDE = [(-8, 8), (0, 40), (-25, 25), (-8, 8), (-8, 8)]
MODELS = {
    'dual       ': dict(H=H_dual,   bounds=WIDE,
                        names=['b', 'g', 'w', 'pS', 'pI']),
    'dual-forced': dict(H=H_dual,
                        bounds=[(-2, 2), (0, 40), (0, 25), (-8, 8), (-8, 8)],
                        names=['b', 'g', 'w', 'pS', 'pI']),
    'single     ': dict(H=H_single, bounds=WIDE,
                        names=['b', 'g', 'w', 'pS', 'pI']),
    'amsh       ': dict(H=H_amsh,
                        bounds=[(0, 40), (-2, 2), (-8, 8), (-8, 8), (-8, 8)],
                        names=['R', 'Q', 'x0', 'wE', 'dP']),
}
TAIL = [(0.1, 15), (1, 1500), (30, 42), (0.0, 1.0)]

def unpack(x, nphase):
    phis = x[9:9 + nphase]
    return x[:5], x[5], x[6], x[7], x[8], (np.repeat(phis, 3) if nphase == 1 else phis)

def model(spec, x, nphase, idx):
    core, A_T, tau, T0, rho, phis = unpack(x, nphase)
    E2 = np.array([0.0, 1.0, rho])[idx][:, None]
    P4 = np.array([0.0, 1.0, 1.0])[idx][:, None]
    return plant(spec['H'](core, E2, P4, clocks(phis)[idx]), A_T, tau, T0)

def rmse(spec, x, nphase, idx):
    return float(np.sqrt(np.mean((model(spec, x, nphase, idx) - OBS[idx]) ** 2)))

ALL = np.array([0, 1, 2]); TRAIN = np.array([0, 1]); TEST = np.array([2])
rng = np.random.default_rng(7)

def fit(spec, nphase, idx, n_starts):
    bnds = spec['bounds'] + TAIL + [(0.0, 24.0)] * nphase
    lo = np.array([b[0] for b in bnds]); hi = np.array([b[1] for b in bnds])
    out = []
    for _ in range(n_starts):
        x0 = lo + rng.random(len(bnds)) * (hi - lo)
        r = minimize(lambda z: rmse(spec, z, nphase, idx), x0, bounds=bnds,
                     method='L-BFGS-B', options=dict(maxiter=500, ftol=1e-14))
        out.append((float(r.fun), r.x))
    return sorted(out, key=lambda s: s[0])

print("=" * 78)
print("IN-SAMPLE FIT TO FULL WAVEFORMS")
print("=" * 78)
best = {}
for name, spec in MODELS.items():
    s1 = fit(spec, 1, ALL, 350)[0]
    s3 = fit(spec, 3, ALL, 350)[0]
    best[name] = s1
    core, A_T, tau, T0, rho, phis = unpack(s1[1], 1)
    M = model(spec, s1[1], 1, ALL)
    amps = M.max(axis=1) - M.min(axis=1)
    print(f"\n{name}  shared-phase {s1[0]:.4f} C   free-phase {s3[0]:.4f} C")
    print("   " + "  ".join(f"{n}={v:.2f}" for n, v in zip(spec['names'], core)) +
          f"  A_T={A_T:.2f} tau={tau:.0f} T0={T0:.2f} rho={rho:.2f} phi={phis[0]:.1f}h")
    print(f"   model amps {np.round(amps,2)}  ratios "
          f"{np.round(amps/amps[0],2)}   (obs ratios "
          f"{np.round((OBS.max(axis=1)-OBS.min(axis=1))/(OBS[0].max()-OBS[0].min()),2)})")

print("\n" + "=" * 78)
print("OUT-OF-SAMPLE: train on non-preg + pregnant, predict Esr1i (rho only free)")
print("=" * 78)
for name, spec in MODELS.items():
    sols = fit(spec, 1, TRAIN, 250)
    thr = sols[0][0] * 1.10
    train = [s for s in sols if s[0] <= thr]
    tests = np.array([min(rmse(spec, np.r_[x[:8], r, x[9:]], 1, TEST)
                          for r in np.linspace(0, 1, 101)) for _, x in train])
    print(f"\n{name}  train best {sols[0][0]:.4f}  ({len(train)} within 10%)")
    print(f"   held-out Esr1i RMSE: best {tests.min():.4f}  "
          f"median {np.median(tests):.4f}  worst {tests.max():.4f}")

fig, axes = plt.subplots(1, 4, figsize=(18, 4.2), sharey=True)
for j, (name, spec) in enumerate(MODELS.items()):
    M = model(spec, best[name][1], 1, ALL)
    for i, l in enumerate(LABELS):
        axes[j].plot(hrs, OBS[i], 'o', ms=3, color=f"C{i}", alpha=0.5)
        axes[j].fill_between(hrs, OBS[i]-SD[i], OBS[i]+SD[i], color=f"C{i}", alpha=0.10)
        axes[j].plot(hrs, M[i], '-', lw=2, color=f"C{i}", label=l)
    axes[j].set_title(f"{name.strip()}   RMSE={best[name][0]:.3f} C")
    axes[j].set_xlabel("hour of day"); axes[j].grid(alpha=0.3)
axes[0].set_ylabel("core temperature (C)"); axes[0].legend(fontsize=8)
plt.tight_layout(); plt.savefig("fig_waveform_fit.png", dpi=130)
print("\nwrote fig_waveform_fit.png")
