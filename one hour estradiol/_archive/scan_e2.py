"""
Scan the combined (n_amp, n_freq) plane, measuring irregularity of E2 ITSELF:
  E2 peak-height CV  = amplitude variation of estradiol
  E2 peak-period CV  = spacing variation of estradiol peaks
Goal: find where BOTH are high (E2 ragged in amplitude AND period).
"""
import numpy as np
from scipy.signal import find_peaks
from pi_cascade_combined import run

K = 14.8

def rcv(x):
    """Median-based coefficient of variation = 1.4826 * MAD / median."""
    x = np.asarray(x, float); med = np.median(x)
    if med == 0: return 0.0
    mad = np.median(np.abs(x - med))
    return 1.4826 * mad / med

def e2_metrics(n_amp, n_freq):
    T, G, F, E, pt, pa = run(n_amp, n_freq, K, K, t_end=14000.0)
    m = T > T[-1] - 5000
    Tt, Ee = T[m], E[m]
    rng = Ee.max() - Ee.min()
    if rng < 1: return 0.0, 0.0
    pk, _ = find_peaks(Ee, prominence=0.15*rng)   # real E2 peaks
    if len(pk) < 3: return 0.0, 0.0
    h = Ee[pk]; ivt = np.diff(Tt[pk])
    return rcv(h), rcv(ivt)                        # median-based CVs

print(" n_amp n_freq |  E2-height-CV  E2-period-CV   min(both)")
best = (0, None)
for na in [6, 8, 10, 12, 16, 20]:
    for nf in [0, 6, 10, 14, 18, 24]:
        hcv, icv = e2_metrics(na, nf)
        both = min(hcv, icv)
        if both > best[0]: best = (both, (na, nf, hcv, icv))
        print(f"   {na:2d}   {nf:2d}   |    {hcv:.3f}        {icv:.3f}       {both:.3f}")
print(f"\nBEST (both high): n_amp={best[1][0]}, n_freq={best[1][1]}  "
      f"-> height-CV {best[1][2]:.2f}, period-CV {best[1][3]:.2f}")
