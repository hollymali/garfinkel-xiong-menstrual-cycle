"""
Spot-check a single parameter set.

Edit the values in the EDIT block below, run `python spotcheck.py`, and it:
  - runs MODEL.run with those parameters
  - prints a short summary (E2 range + orbit period-#)
  - saves the LAST two days of the E2 time series to fig_spotcheck.png

Everything MODEL.run accepts is exposed here, so you can sweep by hand.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")          # save to file, don't open a window
import matplotlib.pyplot as plt

import MODEL
from orbit_period import e2_orbit_period, period_label

# ======================= EDIT THESE =======================
n_amp        = 4            # amplitude-feedback steepness (Hill exponent)
n_freq       = 4            # frequency-feedback steepness (Hill exponent)
tau          = 120          # E2 -> GnRH delay in minutes (0 = no delay)
e0           = 15           # starting E2 (pg/mL)
t_end        = 16000        # total run length (min)  [~11 days]
dt           = 0.25         # integration step (min)
period_range = (70, 80)     # GnRH pulse spacing bounds (min): (fastest, slowest) — mid-late follicular
circ         = False       # daily circadian envelope on pulse amplitude (peaks 2pm)
# ==========================================================

# --- run the model ---
t, G, L, E = MODEL.run(n_amp, n_freq, tau=tau, e0=e0,
                       t_end=t_end, dt=dt, period_range=period_range, circ=circ)

# --- summary over the settled tail (last 6 days) ---
tail = 8640                       # 6 days, matches the period-finder default
settled = E[t > t[-1] - tail]
p = e2_orbit_period(t, E)
print(f"n_amp={n_amp}  n_freq={n_freq}  tau={tau} min  e0={e0} pg/mL")
print(f"settled E2 range: {settled.min():.1f} .. {settled.max():.1f} pg/mL")
print(f"orbit period: {period_label(p)}")

# --- plot the LAST two days of E2 ---
two_days = 2 * 1440
m = t > t[-1] - two_days
hours = (t[m] - t[m][0]) / 60.0

fig, ax = plt.subplots(figsize=(12, 4))
ax.plot(hours, E[m], lw=1.0, color="C4")
for x in range(0, 49, 12):        # day/half-day gridlines
    ax.axvline(x, color="0.85", lw=0.6)
ax.text(0.02, 0.95, period_label(p), transform=ax.transAxes, fontsize=13,
        fontweight="bold", va="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.85))
ax.set_xlabel("time (h)  —  last 2 days of the run")
ax.set_ylabel("E2 (pg/mL)")
ax.set_title(f"spotcheck:  n_amp={n_amp}, n_freq={n_freq}, tau={tau}, e0={e0}  ->  {period_label(p)}")
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig("fig_spotcheck.png", dpi=130)
print("saved fig_spotcheck.png")
