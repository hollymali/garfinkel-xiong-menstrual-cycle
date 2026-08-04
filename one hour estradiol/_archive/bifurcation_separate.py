"""
Two 1-D bifurcations that vary n_freq and n_amp SEPARATELY.
  Left  column: sweep n_freq, n_amp fixed = 8   (matches figure [4])
  Right column: sweep n_amp,  n_freq fixed = 8   (matches figure [5])
  Top row: E2 peak HEIGHTS (amplitude attractor)
  Bot row: E2 peak-to-peak PERIODS in min (timing attractor)
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
from pi_cascade_combined import run
from e2metrics import e2_peaks

K = 14.8

def attractor(n_amp, n_freq):
    T, G, F, E, pt, pa = run(n_amp, n_freq, K, K, t_end=9000.0, dt=0.5)
    return e2_peaks(T, E, tail=4000.0)   # (heights, intervals)

nf_vals = np.linspace(0.0, 26.0, 105)    # sweep n_freq (n_amp=8)
na_vals = np.linspace(2.0, 26.0, 105)    # sweep n_amp  (n_freq=8)

fig, ax = plt.subplots(2, 2, figsize=(15, 9), sharex="col")

# --- left column: sweep n_freq, n_amp = 8 ---
for nf in nf_vals:
    h, iv = attractor(8, nf)
    ax[0, 0].plot([nf]*len(h),  h,  ".", ms=1.6, color="C4", alpha=0.5)
    ax[1, 0].plot([nf]*len(iv), iv, ".", ms=1.6, color="C3", alpha=0.5)
# --- right column: sweep n_amp, n_freq = 8 ---
for na in na_vals:
    h, iv = attractor(na, 8)
    ax[0, 1].plot([na]*len(h),  h,  ".", ms=1.6, color="C4", alpha=0.5)
    ax[1, 1].plot([na]*len(iv), iv, ".", ms=1.6, color="C3", alpha=0.5)

ax[0, 0].set_title("sweep n_freq   (n_amp = 8 fixed)")
ax[0, 1].set_title("sweep n_amp   (n_freq = 8 fixed)")
ax[0, 0].set_ylabel("E2 peak heights\n(attractor)")
ax[1, 0].set_ylabel("E2 peak-to-peak\nperiods (min)")
ax[1, 0].set_xlabel("n_freq"); ax[1, 1].set_xlabel("n_amp")

for a in ax.ravel():
    a.xaxis.set_major_locator(MultipleLocator(2)); a.xaxis.set_minor_locator(AutoMinorLocator(2))
    a.yaxis.set_minor_locator(AutoMinorLocator(2))
    a.grid(alpha=0.3); a.grid(which="minor", alpha=0.12)
for a in ax[0, :]: a.yaxis.set_major_locator(MultipleLocator(2))
for a in ax[1, :]: a.yaxis.set_major_locator(MultipleLocator(20))

fig.suptitle("Separate bifurcations: vary n_freq and n_amp independently "
             "(1 dot=period-1, 2=period-2, band=high-period)", fontsize=12)
fig.tight_layout()
fig.savefig("fig_bifurcation_separate.png", dpi=130)
print("wrote fig_bifurcation_separate.png")
