"""What does the delay-induced ragged E2 look like on the DAY and HOUR scale?
n=3 (physiological), tau=120 min."""
import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from MODEL import run

T, _, _, E = run(3, 3, 120)
m = T > T[-1] - 6000
Tt, Ee = T[m], E[m]
pk, _ = find_peaks(Ee, prominence=0.15*(Ee.max()-Ee.min()))
iv = np.diff(Tt[pk])
print(f"n=3, tau=120:  E2 range {Ee.min():.1f}..{Ee.max():.1f} pg/mL")
print(f"  peak-to-peak intervals: median {np.median(iv):.0f} min, range {iv.min():.0f}-{iv.max():.0f} min")
print(f"  ({np.median(iv)/60:.1f} h typical spacing -> {24*60/np.median(iv):.0f} peaks/day)")

fig, ax = plt.subplots(2, 1, figsize=(12, 6))
# day scale: 2 days
d = (Tt - Tt[0]) < 2*1440
ax[0].plot((Tt[d]-Tt[0])/60, Ee[d], lw=1.0, color="C4")
ax[0].set_title(f"n=3 (physiological) + tau=120 min:  E2 over 2 DAYS")
ax[0].set_xlabel("time (h)"); ax[0].set_ylabel("E2 (pg/mL)"); ax[0].grid(alpha=0.3)
for x in range(0, 49, 12): ax[0].axvline(x, color="0.8", lw=0.6)
# hour scale: 6 h zoom
z = (Tt - Tt[0]) < 6*60
ax[1].plot((Tt[z]-Tt[0])/60, Ee[z], lw=1.4, color="C4", marker=".", ms=3)
ax[1].set_title("zoom: 6 HOURS (are pulses still ~hourly?)")
ax[1].set_xlabel("time (h)"); ax[1].set_ylabel("E2 (pg/mL)"); ax[1].grid(alpha=0.3)
fig.tight_layout(); fig.savefig("fig_delay_dynamics.png", dpi=130)
print("wrote fig_delay_dynamics.png")
