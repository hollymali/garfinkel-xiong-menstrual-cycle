"""Ragged E2 from MODEL.py -> two-arm heat-loss drive -> thermal plant -> CBT.

The thermo fit so far fed the two arms a CONSTANT per-condition E2 level, so the
only time-varying input was the SCN cosine.  Here E2 becomes a real time series
from the cascade (MODEL.py, band (55,65), tau=60, n=4), and the SCN gain is free
to drop because E2 now supplies part of the variance.

Arms (all three are the same animal model, differing in E2 signal):
  non-preg  : cycling, ERalpha intact -> FULL ragged cascade
  pregnant  : GnRH pulsatility suppressed -> smooth, elevated E2
  Esr1i     : pregnant + ERalpha inactivated -> smooth E2 AND reduced transduction (rho)

KEY DESIGN CHOICE (keeps the SD an honest prediction, not a fit):
the E2 fluctuation carries NO free scale parameter.  The cascade's fractional
swing is mapped into the thermo model's E2 units using one physiological number,
preg_ratio = (pregnant E2)/(cycling mean E2), because the thermo model defines
E2=0 at non-pregnant and E2=1 at pregnant.  So:

    E2_in(t) = L_arm + (E2_casc(t) - mean) / (mean * (preg_ratio - 1))

Only the thermo core is fitted, and it is fitted to the MEAN 24 h waveforms only.
The across-day SD is then compared as an out-of-sample prediction.
"""
import sys
import numpy as np
from scipy.optimize import minimize
from scipy.signal import lfilter

sys.path.insert(0, r"C:\Users\holly\source\garfinkel-xiong-menstrual-cycle\one hour estradiol")
import MODEL

DAY = 1440.0
NDAYS = 6                 # the CBT panels span ~5 days
DT = 0.25
NB = 48                   # 48 half-hour bins, matching the digitised fold

# ---------------- observed ----------------
fold = np.load("fold.npy", allow_pickle=True).item()
LABELS = ["non-preg", "pregnant", "Esr1i"]
hrs = fold[LABELS[0]][0]
OBS = np.vstack([fold[l][1] for l in LABELS])
OBS_SD = np.vstack([fold[l][2] for l in LABELS])
AMP = OBS.max(axis=1) - OBS.min(axis=1)


# ---------------- E2 time series per arm ----------------
def cascade_e2(n, tau, ndays=NDAYS, burn=4.0):
    """Ragged E2(t) at the late-follicular band. Returns (t_min, E2) after burn-in."""
    t_end = int((ndays + burn) * DAY)
    t, G, L, E = MODEL.run(n, n, tau=tau, t_end=t_end, dt=DT,
                           period_range=(55, 65), e0=15)
    m = t >= burn * DAY
    return t[m] - t[m][0], E[m]


def build_inputs(preg_ratio):
    """E2_in(t) for the three arms, on a common time grid."""
    t, e_ragged = cascade_e2(n=4, tau=60)             # cycling, ERalpha intact
    _, e_flat = cascade_e2(n=0, tau=0)                # feedback off -> regular pulses
    mean_r = e_ragged.mean()
    # fractional swing -> thermo E2 units (no free scale)
    frac_ragged = (e_ragged - mean_r) / (mean_r * (preg_ratio - 1.0))
    mean_f = e_flat.mean()
    frac_flat = (e_flat - mean_f) / (mean_f * (preg_ratio - 1.0))
    return t, frac_ragged, frac_flat


def sig(x):
    return 1.0 / (1.0 + np.exp(-np.clip(x, -60, 60)))


def H_dual(p, E2, P4, c):
    b, g, w, pS, pI = p
    return sig(b + g * c + w * E2 + pS * P4) - sig(b - g * c + w * E2 + pI * P4)


def run_arms(x, t, frac_ragged, frac_flat, rho_flat_scale=0.0):
    """Integrate the thermal plant over NDAYS for all three arms; return CBT(t)."""
    core = x[:5]
    A_T, tau_T, T0, rho = x[5:9]
    ph = x[9:12]
    L = np.array([0.0, 1.0, rho])
    # Esr1i: ERalpha inactivated -> E2 fluctuation is not transduced
    fr = [frac_ragged, frac_flat, frac_flat * rho_flat_scale]
    P4 = [0.0, 1.0, 1.0]

    out = []
    for i in range(3):
        c = np.cos(2 * np.pi * (t - ph[i] * 60.0) / DAY)
        H = H_dual(core, L[i] + fr[i], P4[i], c)
        drive = -A_T * H + T0
        # causal first-order plant y[k] = (1-a)y[k-1] + a x[k]  (dT/dt=(drive-T)/tau_T)
        a = min(DT / max(tau_T, 1e-6), 1.0)
        T = lfilter([a], [1.0, -(1.0 - a)], drive, zi=np.array([(1.0 - a) * drive[0]]))[0]
        out.append(T)
    return out


def fold_index(t):
    """Precompute the (day, clock-bin) flat index and its counts."""
    hour = (t / 60.0) % 24.0
    day = (t // DAY).astype(int)
    b = np.minimum((hour / (24.0 / NB)).astype(int), NB - 1)
    ndays = int(day.max()) + 1
    flat = day * NB + b
    counts = np.bincount(flat, minlength=ndays * NB).astype(float)
    return flat, counts, ndays


def fold_stats(T, flat, counts, ndays):
    """Fold to 48 half-hour clock bins across days -> (mean waveform, SD across days)."""
    tot = np.bincount(flat, weights=T, minlength=ndays * NB)
    with np.errstate(invalid="ignore", divide="ignore"):
        grid = (tot / counts).reshape(ndays, NB)
    return np.nanmean(grid, axis=0), np.nanstd(grid, axis=0, ddof=1)


def predict(x, t, fr_r, fr_f, idx, rho_flat_scale=0.0):
    Ts = run_arms(x, t, fr_r, fr_f, rho_flat_scale)
    flat, counts, ndays = idx
    means, sds = [], []
    for T in Ts:
        m, s = fold_stats(T, flat, counts, ndays)
        means.append(m)
        sds.append(s)
    return np.vstack(means), np.vstack(sds)


def loss_mean_only(x, t, fr_r, fr_f, idx):
    """Same normalised-RMS loss as freephase.py -- MEAN waveforms only."""
    M, _ = predict(x, t, fr_r, fr_f, idx)
    r = M - OBS
    return float(np.sqrt(np.mean((np.sqrt(np.mean(r ** 2, axis=1)) / AMP) ** 2)))


if __name__ == "__main__":
    PREG_RATIO = float(sys.argv[1]) if len(sys.argv) > 1 else 8.0
    t, fr_r, fr_f = build_inputs(PREG_RATIO)
    print(f"preg_ratio = {PREG_RATIO}  (pregnant E2 / cycling-mean E2)")
    print(f"cascade E2 fractional swing -> thermo units: "
          f"ragged {fr_r.min():+.3f}..{fr_r.max():+.3f}, flat {fr_f.min():+.3f}..{fr_f.max():+.3f}")
    print(f"grid: {len(t)} samples, {NDAYS} days\n")

    idx = fold_index(t)
    WIDE = [(-8, 8), (0, 40), (-25, 25), (-8, 8), (-8, 8)]
    TAIL = [(0.1, 15), (1, 1500), (30, 42), (0.0, 1.0)] + [(0.0, 24.0)] * 3
    bnds = WIDE + TAIL
    lo = np.array([b[0] for b in bnds]); hi = np.array([b[1] for b in bnds])
    rng = np.random.default_rng(31)
    best = None
    NSTART = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    for i in range(NSTART):
        x0 = lo + rng.random(len(bnds)) * (hi - lo)
        r = minimize(loss_mean_only, x0, args=(t, fr_r, fr_f, idx), bounds=bnds,
                     method="L-BFGS-B", options=dict(maxiter=120, ftol=1e-12))
        if best is None or r.fun < best[0]:
            best = (float(r.fun), r.x)
            print(f"  start {i:2d}: loss {r.fun:.4f}")
    e, x = best
    np.save("best_coupled.npy", x)

    M, S = predict(x, t, fr_r, fr_f, idx)
    print(f"\nfitted loss (mean waveforms only) = {e:.4f}")
    print(f"  SCN gain g = {x[1]:.2f}   E2 weight w = {x[2]:.2f}   "
          f"A_T = {x[5]:.2f}  tau = {x[6]:.1f} min  rho = {x[8]:.3f}")
    print(f"\n{'arm':10s} {'amp obs':>8} {'amp mod':>8} | {'SD obs':>7} {'SD pred':>8}   <- SD is OUT-OF-SAMPLE")
    for i, l in enumerate(LABELS):
        print(f"{l:10s} {AMP[i]:8.2f} {M[i].max()-M[i].min():8.2f} | "
              f"{OBS_SD[i].mean():7.2f} {S[i].mean():8.2f}")
