"""
Orbit-period heatmap with the HONEST two-lag budget, over the physiological Hill range.

    tau_prod = 30 min   LH -> E2 steroidogenesis (ovarian response to the LH pulse)
    tau_fb   = 30 min   E2 -> KNDy genomic feedback read
    TOTAL loop lag = 60 min

Both lags come from MODEL.py's single `tau`, which sits on BOTH arms of the loop: the
E2 read that sets pulse amplitude and frequency, and the LH -> E2 production step.  So
tau=30 IS the 30 + 30 budget, and the loop lag this figure tests is 2*tau = 60 min.

Only the TOTAL loop lag matters (placement was proven irrelevant), and the requirement is
~1 PULSE INTERVAL, not a fixed number of minutes. So a 60 min budget only closes if the pulse
interval is ~60 min -- i.e. the LATE-FOLLICULAR band (55, 65) min (Sollenberger 1990 measures
~53 min late-follicular). At the early-follicular band (70, 80) the same 60 min lag is BELOW
threshold and everything collapses to period-1.

Two panels, identical except for the lag, so the lag is the only thing that can explain a
difference:
    LEFT   no lag        (tau_prod = tau_fb = 0)
    RIGHT  realistic lag (tau_prod = tau_fb = 30)
both at period_range = (55, 65).
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from MODEL import run
from orbit_period import e2_orbit_period

HIGH = 40
BAND = (55.0, 65.0)
TAU = 30.0          # per arm; loop lag = 2*TAU = 60 min
TAU_PROD = TAU_FB = TAU
CACHE = "heatmap_realdelays.npy"

# 0.5 step -> 11 x 11 = 121 cells per panel, matching the grid the earlier lag work used
# (its headline was "55/121 ragged with lag vs 0/121 without"), so the numbers are comparable.
STEP = 0.5
n_amps = np.arange(0.0, 5.01, STEP)
n_freqs = np.arange(0.0, 5.01, STEP)


def orbit(n_amp, n_freq, tau):
    T, _, _, E = run(n_amp, n_freq, tau=tau,
                     t_end=12960, dt=0.25, period_range=BAND)
    p = e2_orbit_period(T, E, tail=8640)
    if p == 0:
        return np.nan
    return HIGH if p < 0 else float(p)


def sweep(tau, tag, P, other):
    """row-level checkpointing: a kill costs one row, not the whole sweep"""
    for i, na in enumerate(n_amps):
        if not np.all(np.isnan(P[i])):      # row already done on a previous run
            continue
        for j, nf in enumerate(n_freqs):
            P[i, j] = orbit(na, nf, tau)
        np.save(CACHE, np.stack([P, other] if tag == "no lag" else [other, P]))
        print(f"  [{tag}] n_amp={na:.2f} row done (checkpointed)", flush=True)
    return P


if os.path.exists(CACHE):
    d = np.load(CACHE)
    P0, P1 = d[0], d[1]
    print(f"resuming from cache: no-lag {int(np.sum(~np.isnan(P0)))}/121 cells, "
          f"lag {int(np.sum(~np.isnan(P1)))}/121 cells", flush=True)
else:
    P0 = np.full((len(n_amps), len(n_freqs)), np.nan)
    P1 = np.full((len(n_amps), len(n_freqs)), np.nan)

if np.all(np.isnan(P0)) or np.any(np.all(np.isnan(P0), axis=1)):
    print("no-lag sweep:", flush=True)
    P0 = sweep(0.0, "no lag", P0, P1)
if np.all(np.isnan(P1)) or np.any(np.all(np.isnan(P1), axis=1)):
    print("realistic-lag sweep:", flush=True)
    P1 = sweep(TAU, "lag 60", P1, P0)
np.save(CACHE, np.stack([P0, P1]))


def tally(P, name):
    v, c = np.unique(P[~np.isnan(P)], return_counts=True)
    d = {int(a): int(b) for a, b in zip(v, c)}
    n = int(np.sum(~np.isnan(P)))
    ragged = int(np.sum((P > 2) & ~np.isnan(P)))
    print(f"{name}: {d}")
    print(f"   ragged (period>2 incl. high): {ragged}/{n}"
          f"   period-1: {d.get(1,0)}   period-2: {d.get(2,0)}")
    return ragged, n


r0, n0 = tally(P0, "no lag        ")
r1, n1 = tally(P1, "realistic lag ")

bounds = [0.5, 1.5, 2.5, 4.5, 8.5, 16.5, HIGH + 0.5]
centers = [1, 2, 3.5, 6.5, 12.5, (16.5 + HIGH + 0.5) / 2]
labels = ["period-1", "period-2", "period 3-4", "period 5-8", "period 9-16", ">16 (high)"]
cmap = plt.get_cmap("viridis", len(bounds) - 1)
cmap.set_bad("0.85")
norm = BoundaryNorm(bounds, cmap.N)
st = STEP
ext = [n_freqs.min() - st / 2, n_freqs.max() + st / 2,
       n_amps.min() - st / 2, n_amps.max() + st / 2]

fig, ax = plt.subplots(1, 2, figsize=(14.5, 6.4), sharey=True)
for a, P, ttl, rr, nn in [
        (ax[0], P0, "NO lag", r0, n0),
        (ax[1], P1, f"steroidogenesis {TAU_PROD:.0f} + genomic {TAU_FB:.0f} = 60 min", r1, n1)]:
    im = a.imshow(P, origin="lower", aspect="auto", extent=ext, cmap=cmap, norm=norm)
    a.set_xlabel("n_freq  (frequency-feedback steepness)")
    a.set_title(f"{ttl}\nragged (period>2): {rr}/{nn}")
ax[0].set_ylabel("n_amp  (amplitude-feedback steepness)")
cbar = fig.colorbar(im, ax=ax, ticks=centers, spacing="uniform", fraction=0.035)
cbar.ax.set_yticklabels(labels)
cbar.set_label("E2 orbit period")
fig.suptitle(f"E2 orbit period, physiological Hill range, GnRH pulse band "
             f"{BAND[0]:.0f}-{BAND[1]:.0f} min (late follicular)", y=0.99)
fig.savefig("fig_heatmap_realdelays.png", dpi=130, bbox_inches="tight")
print("\nwrote fig_heatmap_realdelays.png, heatmap_realdelays.npy")
