"""Close look at n_amp=8, n_freq=4: is it irregular in timing AND still alternating in amplitude?"""
import numpy as np
import matplotlib.pyplot as plt
from pi_cascade_combined import run

K = 14.8
T, G, F, E, pt, pa = run(8, 4, K, K)

# E2 peaks (local maxima) in the settled tail
m = T > T[-1] - 3000
Tt, Ee = T[m], E[m]
pk = np.where((Ee[1:-1] > Ee[:-2]) & (Ee[1:-1] > Ee[2:]))[0] + 1
pk_t, pk_h = Tt[pk], Ee[pk]
iv = np.diff(pk_t)
print(f"peak heights (last 12): {np.round(pk_h[-12:],1)}")
print(f"intervals   (last 12):  {np.round(iv[-12:],0)}")
print(f"height CV {pk_h.std()/pk_h.mean():.2f}   interval CV {iv.std()/iv.mean():.2f}")

fig, ax = plt.subplots(2, 1, figsize=(11, 6.5))
w = T > T[-1] - 40*60
ax[0].plot(T[w]/60, E[w], lw=1.1, color="C4")
pw = (pk_t > T[-1]-40*60)
ax[0].plot(pk_t[pw]/60, pk_h[pw], "v", ms=6, color="C3")
ax[0].set_ylabel("E2 (pg/mL)"); ax[0].grid(alpha=0.3)
ax[0].set_title("n_amp=8, n_freq=4:  E2 time series (triangles = peaks)")

idx = np.arange(len(pk_h[-16:]))
ax[1].plot(idx, pk_h[-16:], "-o", color="C3")
ax[1].set_ylabel("peak height"); ax[1].set_xlabel("pulse number")
ax[1].grid(alpha=0.3)
ax[1].set_title("peak-height sequence  (does it still alternate tall/short?)")
fig.tight_layout(); fig.savefig("fig_look_n4n8.png", dpi=130)
print("wrote fig_look_n4n8.png")
