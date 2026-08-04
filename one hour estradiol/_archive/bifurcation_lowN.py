"""
Bifurcation over the PHYSIOLOGICAL range n <= 5 (realistic Hill coefficients).
Tied diagonal n_amp = n_freq = n, fine resolution.  Top: E2 peak heights;
bottom: E2 peak-to-peak periods.
"""
import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
from bifurcation import run, e2_peaks

ns = np.linspace(0.0, 6.0, 160)
fig, ax = plt.subplots(2, 1, figsize=(11, 8), sharex=True)
for n in ns:
    T, E = run(n, t_end=12000.0)
    h, iv = e2_peaks(T, E, tail=6000.0)
    ax[0].plot([n]*len(h),  h,  ".", ms=2.0, color="C4", alpha=0.55)
    ax[1].plot([n]*len(iv), iv, ".", ms=2.0, color="C3", alpha=0.55)
ax[0].set_ylabel("E2 peak heights\n(attractor)")
ax[0].set_title("Bifurcation over the PHYSIOLOGICAL range  n = n_amp = n_freq <= 6\n"
                "(realistic Hill coefficients ~ 2-5)")
ax[1].set_ylabel("E2 peak-to-peak\nperiods (min)")
ax[1].set_xlabel("feedback steepness  n (= n_amp = n_freq)")
ax[0].axvspan(0, 5, color="green", alpha=0.06)   # shade the physiological band
ax[1].axvspan(0, 5, color="green", alpha=0.06)
for a in ax:
    a.xaxis.set_major_locator(MultipleLocator(0.5)); a.xaxis.set_minor_locator(AutoMinorLocator(2))
    a.yaxis.set_minor_locator(AutoMinorLocator(2))
    a.grid(alpha=0.3); a.grid(which="minor", alpha=0.12)
fig.tight_layout(); fig.savefig("fig_bifurcation_lowN.png", dpi=130)
print("wrote fig_bifurcation_lowN.png")
