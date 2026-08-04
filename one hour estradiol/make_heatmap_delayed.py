"""
Orbit-period heatmap over the PHYSIOLOGICAL range (n <= 5) WITH a physiological
E2->GnRH delay (tau = 120 min).  Robust settings.  Compare to the no-delay
version (which was all period-1/2).
"""
import os
import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from MODEL import run
from orbit_period import e2_orbit_period

HIGH = 40; TAU = 120.0
CACHE = "heatmap_delayed.npy"

n_amps  = np.arange(0.0, 5.01, 0.25)
n_freqs = np.arange(0.0, 5.01, 0.25)

def orbit(n_amp, n_freq):
    T, _, _, E = run(n_amp, n_freq, TAU)
    p = e2_orbit_period(T, E, tail=8000)   # your e2_peaks_intervals (min_peak) + orbit_period
    if p == 0: return np.nan
    return HIGH if p < 0 else float(p)

if os.path.exists(CACHE):
    P = np.load(CACHE); print("loaded cache")
else:
    P = np.full((len(n_amps), len(n_freqs)), np.nan)
    for i, na in enumerate(n_amps):
        for j, nf in enumerate(n_freqs):
            P[i, j] = orbit(na, nf)
        print(f"n_amp={na} row done")
    np.save(CACHE, P)

vals, counts = np.unique(P[~np.isnan(P)], return_counts=True)
print("period tally:", {int(v): int(c) for v, c in zip(vals, counts)})

bounds  = [0.5, 1.5, 2.5, 4.5, 8.5, 16.5, HIGH + 0.5]
centers = [1, 2, 3.5, 6.5, 12.5, (16.5 + HIGH + 0.5)/2]
labels  = ["period-1", "period-2", "period 3-4", "period 5-8", "period 9-16", ">16 (high)"]
cmap = plt.get_cmap("viridis", len(bounds)-1); cmap.set_bad("0.85")
norm = BoundaryNorm(bounds, cmap.N)
st = 0.25
ext = [n_freqs.min()-st/2, n_freqs.max()+st/2, n_amps.min()-st/2, n_amps.max()+st/2]
fig, ax = plt.subplots(figsize=(8.5, 7))
im = ax.imshow(P, origin="lower", aspect="auto", extent=ext, cmap=cmap, norm=norm)
ax.set_xlabel("n_freq  (positive frequency-feedback steepness)")
ax.set_ylabel("n_amp  (negative amplitude-feedback steepness)")
ax.set_title(f"E2 orbit period, PHYSIOLOGICAL range (n<=5) WITH delay tau={TAU:.0f} min\n"
             f"(compare: no delay was ALL period-1/2)")
cbar = fig.colorbar(im, ax=ax, ticks=centers, spacing="uniform")
cbar.ax.set_yticklabels(labels); cbar.set_label("orbit period")
fig.tight_layout(); fig.savefig("fig_heatmap_delayed.png", dpi=130)
print("wrote fig_heatmap_delayed.png")
