"""Quick check: figure-4 setup (n_amp=8 fixed) but push n_freq higher, to 24."""
import numpy as np
from pi_cascade_combined import run

K = 14.8
print("n_amp = 8 fixed, sweeping n_freq:")
for nf in [0, 2, 4, 8, 12, 16, 24]:
    T, G, F, E, pt, pa = run(8, nf, K, K)
    iv = np.diff(pt[pt > T[-1]-4000])
    cv = iv.std()/iv.mean() if len(iv) > 1 else 0.0
    m = T > T[-1]-3000
    print(f"  n_freq={nf:2d}:  interval CV {cv:.3f}   "
          f"interval {iv.min():.0f}-{iv.max():.0f} min   E2 {E[m].min():.0f}-{E[m].max():.0f}")
