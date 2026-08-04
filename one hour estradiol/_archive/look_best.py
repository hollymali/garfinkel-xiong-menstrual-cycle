"""Time series of E2 at the regime with both amplitude AND period variation."""
import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from pi_cascade_combined import run

def rcv(x):
    """MAD/median: median absolute deviation as a fraction of the median."""
    x = np.asarray(x, float); med = np.median(x)
    return np.median(np.abs(x-med))/med if med else 0.0

K = 14.8; n_amp, n_freq = 16, 14
T, G, F, E, pt, pa = run(n_amp, n_freq, K, K, t_end=16000.0)
m = T > T[-1] - 60*60
Tt, Ee = T[m]/60, E[m]
pk, _ = find_peaks(Ee, prominence=0.15*(Ee.max()-Ee.min()))
h = Ee[pk]; ivt = np.diff(Tt[pk])*60
print(f"n_amp={n_amp}, n_freq={n_freq}: E2 peaks n={len(pk)}")
print(f"  peak heights (MAD/median) = {rcv(h):.2f}   -> AMPLITUDE variation")
print(f"  peak periods (MAD/median) = {rcv(ivt):.2f}   -> PERIOD variation")
print(f"  heights: {np.round(h[-10:],0)}")
print(f"  periods(min): {np.round(ivt[-10:],0)}")

fig, ax = plt.subplots(figsize=(12, 4.2))
ax.plot(Tt, Ee, lw=1.1, color="C4")
ax.plot(Tt[pk], h, "v", ms=6, color="C3")
ax.set_xlabel("time (h)"); ax.set_ylabel("E2 (pg/mL)"); ax.grid(alpha=0.3)
ax.set_title(f"E2 with BOTH amplitude & period variation  (n_amp={n_amp}, n_freq={n_freq})\n"
             f"E2 peak-height MAD/median = {rcv(h):.2f}   |   E2 peak-period MAD/median = {rcv(ivt):.2f}",
             fontsize=10)
fig.tight_layout(); fig.savefig("fig_look_best.png", dpi=130)
print("wrote fig_look_best.png")
