"""
Heatmap over the (n_freq, n_amp) plane, colored by the E2 ORBIT PERIOD
(# of E2 peaks until the pattern repeats): period-1, -2, -4, ..., or high-period.
DISCRETE (chunky) colorbar -- one flat color per period band.
The 2-D complement to the single ragged time series (n_freq=14, n_amp=16).
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from pi_cascade_combined import run
from e2metrics import e2_orbit_period

K = 14.8
HIGH = 40   # display value for "high-period" (orbit_period returned -1)
CACHE = "heatmap_P.npy"

n_amps  = np.arange(2, 25)    # y
n_freqs = np.arange(0, 25)    # x

def orbit(n_amp, n_freq):
    T, G, F, E, pt, pa = run(n_amp, n_freq, K, K, t_end=10000.0, dt=0.5)
    p = e2_orbit_period(T, E, tail=5000.0)
    if p == 0:  return np.nan
    return HIGH if p < 0 else float(p)

# --- compute (or load cached) grid ----------------------------------------
if os.path.exists(CACHE):
    P = np.load(CACHE); print("loaded cached grid")
else:
    P = np.full((len(n_amps), len(n_freqs)), np.nan)
    for i, na in enumerate(n_amps):
        for j, nf in enumerate(n_freqs):
            P[i, j] = orbit(na, nf)
        print(f"n_amp={na} row done")
    np.save(CACHE, P)

# --- discrete colorbar -----------------------------------------------------
bounds  = [0.5, 1.5, 2.5, 4.5, 8.5, 16.5, HIGH + 0.5]        # 6 period bands
centers = [1, 2, 3.5, 6.5, 12.5, (16.5 + HIGH + 0.5)/2]
labels  = ["period-1", "period-2", "period 3-4", "period 5-8", "period 9-16", ">16 (high-period)"]
cmap = plt.get_cmap("viridis", len(bounds) - 1)
cmap.set_bad("0.85")                                         # NaN (undecidable) -> light gray
norm = BoundaryNorm(bounds, cmap.N)

ext = [n_freqs.min()-0.5, n_freqs.max()+0.5, n_amps.min()-0.5, n_amps.max()+0.5]
fig, ax = plt.subplots(figsize=(9.5, 7))
im = ax.imshow(P, origin="lower", aspect="auto", extent=ext, cmap=cmap, norm=norm)
ax.set_xlabel("n_freq  (positive frequency-feedback steepness)")
ax.set_ylabel("n_amp  (negative amplitude-feedback steepness)")
ax.set_title("E2 orbit period over the (n_freq x n_amp) plane\n(# of E2 peaks until the pattern repeats)")

cbar = fig.colorbar(im, ax=ax, ticks=centers, spacing="uniform")
cbar.ax.set_yticklabels(labels)
cbar.set_label("orbit period")

ax.plot(14, 16, marker="*", ms=16, color="red", mec="white", mew=1.0)
ax.annotate("time series", (14, 16), textcoords="offset points", xytext=(8, 6),
            color="white", fontsize=9, weight="bold")

fig.tight_layout()
fig.savefig("fig_heatmaps.png", dpi=130)
print("wrote fig_heatmaps.png")
