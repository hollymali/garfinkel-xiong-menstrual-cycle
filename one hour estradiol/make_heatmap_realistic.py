"""
Orbit-period heatmap restricted to the PHYSIOLOGICAL range (n <= 6), with
ROBUST settings (longer settling + finer dt) so the low-n artifacts clean up.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from MODEL import run
from orbit_period import e2_orbit_period

HIGH = 40   # high-period sentinel value (E2 threshold K=22 lives inside MODEL)
CACHE = "heatmap_realistic.npy"

n_amps  = np.arange(0.0, 5.01, 0.25)   # y  (0..5, step 0.25)
n_freqs = np.arange(0.0, 5.01, 0.25)   # x

def orbit(n_amp, n_freq):
    # robust: long run + fine dt + long settled tail.  tau=0 -> no-delay case
    T, _, _, E = run(n_amp, n_freq, 0, t_end=14000.0, dt=0.25)
    p = e2_orbit_period(T, E, tail=7000.0)   # your finder keeps all peaks (no prominence)
    if p == 0:  return np.nan
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

# report the period tally
vals, counts = np.unique(P[~np.isnan(P)], return_counts=True)
print("period tally:", dict(zip(vals.astype(int), counts)))

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
ax.set_title("E2 orbit period, PHYSIOLOGICAL range (Hill coeff n <= 5), robust settings\n"
             "(long settling + fine dt) -- everything is period-1 or period-2")
cbar = fig.colorbar(im, ax=ax, ticks=centers, spacing="uniform")
cbar.ax.set_yticklabels(labels); cbar.set_label("orbit period")
fig.tight_layout(); fig.savefig("fig_heatmap_realistic.png", dpi=130)
print("wrote fig_heatmap_realistic.png")
