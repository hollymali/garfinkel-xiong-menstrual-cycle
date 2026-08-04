"""Compare the orbit-period heatmap at tolerance 3% vs 1% (side by side)."""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from pi_cascade_combined import run
from e2metrics import e2_orbit_period

K = 14.8; HIGH = 40
n_amps  = np.arange(2, 25)
n_freqs = np.arange(0, 25)

def grid(tol):
    cache = "heatmap_P.npy" if abs(tol-0.03) < 1e-9 else f"heatmap_P_tol{int(round(tol*1000))}.npy"
    if os.path.exists(cache):
        print(f"loaded {cache}"); return np.load(cache)
    P = np.full((len(n_amps), len(n_freqs)), np.nan)
    for i, na in enumerate(n_amps):
        for j, nf in enumerate(n_freqs):
            T, G, F, E, pt, pa = run(na, nf, K, K, t_end=10000.0, dt=0.5)
            p = e2_orbit_period(T, E, tail=5000.0, tol=tol)
            P[i, j] = np.nan if p == 0 else (HIGH if p < 0 else float(p))
        print(f"tol={tol} n_amp={na} done")
    np.save(cache, P); return P

P3 = grid(0.03)
P1 = grid(0.01)
print(f"high-period cells: tol3%={np.sum(P3==HIGH)}  tol1%={np.sum(P1==HIGH)}")

bounds  = [0.5, 1.5, 2.5, 4.5, 8.5, 16.5, HIGH + 0.5]
centers = [1, 2, 3.5, 6.5, 12.5, (16.5 + HIGH + 0.5)/2]
labels  = ["period-1", "period-2", "period 3-4", "period 5-8", "period 9-16", ">16 (high-period)"]
cmap = plt.get_cmap("viridis", len(bounds)-1); cmap.set_bad("0.85")
norm = BoundaryNorm(bounds, cmap.N)
ext = [n_freqs.min()-0.5, n_freqs.max()+0.5, n_amps.min()-0.5, n_amps.max()+0.5]

fig, axes = plt.subplots(1, 2, figsize=(15, 6.2))
for ax, P, t in zip(axes, [P3, P1], ["tolerance = 3%", "tolerance = 1%"]):
    im = ax.imshow(P, origin="lower", aspect="auto", extent=ext, cmap=cmap, norm=norm)
    ax.plot(14, 16, marker="*", ms=14, color="red", mec="white")
    ax.set_xlabel("n_freq"); ax.set_ylabel("n_amp"); ax.set_title(f"orbit period, {t}")
cbar = fig.colorbar(im, ax=axes, ticks=centers, spacing="uniform")
cbar.ax.set_yticklabels(labels); cbar.set_label("orbit period")
fig.savefig("fig_tol_compare.png", dpi=130, bbox_inches="tight")
print("wrote fig_tol_compare.png")
