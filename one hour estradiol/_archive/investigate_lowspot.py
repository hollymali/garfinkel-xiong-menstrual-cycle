"""Investigate the yellow (high-period) cell at low n (~ n_amp=3, n_freq=5):
real high-period orbit, or a peak-detection artifact?"""
import numpy as np
import matplotlib.pyplot as plt
from pi_cascade_combined import run
from e2metrics import e2_peaks, e2_orbit_period, period_label

K = 14.8
print("orbit period in the low-n neighborhood:")
for na in [2, 3, 4]:
    for nf in [4, 5, 6]:
        T, G, F, E, pt, pa = run(na, nf, K, K, t_end=10000.0, dt=0.5)
        print(f"  n_amp={na} n_freq={nf}: {period_label(e2_orbit_period(T, E, tail=5000.0))}")

na, nf = 3, 5
T, G, F, E, pt, pa = run(na, nf, K, K, t_end=10000.0, dt=0.5)
h, iv = e2_peaks(T, E, tail=5000.0)
print(f"\nn_amp={na}, n_freq={nf}:  E2 range {E[-6000:].min():.2f}..{E[-6000:].max():.2f} pg/mL")
print(f"  {len(h)} peaks;  heights(last 10) = {np.round(h[-10:],2)}")
print(f"  intervals(last 10) = {np.round(iv[-10:],0)}")

m = T > T[-1] - 3000
fig, ax = plt.subplots(figsize=(11, 3.6))
ax.plot(T[m]/60, E[m], lw=1.0, color="C4")
mm = (T > T[-1]-3000)
from scipy.signal import find_peaks
Ee = E[mm]; pk, _ = find_peaks(Ee, prominence=0.15*(Ee.max()-Ee.min()))
ax.plot((T[mm]/60)[pk], Ee[pk], "v", ms=6, color="C3")
ax.set_xlabel("time (h)"); ax.set_ylabel("E2 (pg/mL)")
ax.set_title(f"Low-n 'yellow spot' n_amp={na}, n_freq={nf}:  what is E2 actually doing?")
ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig("fig_lowspot.png", dpi=130)
print("wrote fig_lowspot.png")
