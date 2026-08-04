"""
We define "orbit-period" as the smallest period in which the peak heights and peak intervals match.
"""

# Import libraries
import numpy as np
from scipy.signal import find_peaks

"""
Transient is 3 days so 4320 minutes, mid-late follicular is around 6 days so 8640 minutes
"""
def e2_peaks_intervals(T, E, tail=8640, min_peak=0.05):
    # Keep only time and estrogen after 3 day transient
    after_transient = T > T[-1] - tail
    at_time, at_e2 = T[after_transient], E[after_transient]
    # Make sure E2 is not flat
    if at_e2.min() == at_e2.max():
        return np.array([]), np.array([])
    # Use scipy to find peaks that rise at least min_peak (pg/mL) above their surroundings
    peaks, _ = find_peaks(at_e2, prominence=min_peak)
    # Get height of peaks and intervals
    heights = at_e2[peaks]
    intervals = np.diff(at_time[peaks])
    # Return
    return heights, intervals

def orbit_period(h, iv, tol=0.03, pmax=30):
    # If there is only 2 peaks
    if len(h) <= 2:
        return 0
    # Make height and intervals floats
    h = np.asarray(h, float)
    iv = np.asarray(iv, float)
    # Set up loop to test period, with p as suspect # of peaks before repeat
    for p in range(1, min(pmax, len(h)//2) + 1):
        # Test maximal peak divergence
        max_h_var = np.max(np.abs(h[p:] - h[:-p]) / h[:-p])
        max_iv_var = np.max(np.abs(iv[p:] - iv[:-p]) / iv[:-p])
        # Make sure peaks and intervals do not diverge by more than 3%
        if max_h_var < tol and max_iv_var < tol:
            return p
    return -1

# Now put both functions together
def e2_orbit_period(T, E, tail=8640, tol=0.03, pmax=30):
    h, iv = e2_peaks_intervals(T, E, tail)
    return orbit_period(h, iv, tol, pmax)

# Return label
def period_label(p):
    if p == 0:
        return "Not enough peaks"
    if p < 0:
        return ">30 (high period)"
    else:
        return f"period-{p}"