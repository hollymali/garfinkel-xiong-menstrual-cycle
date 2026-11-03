"""
Shared E2-based variation metrics (median/MAD, no GnRH quantities).
  E2_amplitude_var = MAD/median of E2 peak HEIGHTS
  E2_period_var    = MAD/median of E2 peak SPACINGS
"""
import numpy as np
from scipy.signal import find_peaks

def rcv(x):
    """MAD/median: median absolute deviation as a fraction of the median."""
    x = np.asarray(x, float)
    if len(x) < 2: return 0.0
    med = np.median(x)
    return np.median(np.abs(x - med))/med if med else 0.0

def e2_peaks(T, E, tail=4000.0, prom=0.15):
    """Return (peak heights, peak-to-peak intervals) of E2 in the settled tail.
    prom = min peak prominence as a fraction of the E2 range: a bump must rise
    this much above its surroundings to count as its OWN peak (smaller shoulders
    merge into the neighbouring peak)."""
    m = T > T[-1] - tail
    Tt, Ee = T[m], E[m]
    rng = Ee.max() - Ee.min()
    if rng < 1: return np.array([]), np.array([])
    pk, _ = find_peaks(Ee, prominence=prom*rng)
    return Ee[pk], np.diff(Tt[pk])

def e2_variation(T, E, tail=4000.0):
    """Return (E2_amplitude_var, E2_period_var)."""
    h, iv = e2_peaks(T, E, tail)
    return rcv(h), rcv(iv)

def orbit_period(h, iv, tol=0.03, pmax=30):
    """Smallest p such that BOTH the E2 peak-height and peak-interval sequences
    repeat with period p (the orbit period / period-#). Returns:
       p (1,2,4,...)  = period-p orbit
      -1              = no repeat up to pmax (high-period / effectively irregular)
       0              = too few peaks to decide
    """
    m = min(len(h), len(iv))
    if m < 6: return 0
    h = np.asarray(h[-m:], float); iv = np.asarray(iv[-m:], float)
    hmed, imed = np.median(h), np.median(iv)
    for p in range(1, min(pmax, m//2) + 1):
        dh = np.mean(np.abs(h[p:] - h[:-p]))/hmed if hmed else 1.0
        di = np.mean(np.abs(iv[p:] - iv[:-p]))/imed if imed else 1.0
        if dh < tol and di < tol:
            return p
    return -1

def e2_orbit_period(T, E, tail=4000.0, tol=0.03, prom=0.15):
    h, iv = e2_peaks(T, E, tail, prom=prom)
    return orbit_period(h, iv, tol=tol)

def period_label(p):
    return "n/a" if p == 0 else (">30 (high-period)" if p < 0 else f"period-{p}")
