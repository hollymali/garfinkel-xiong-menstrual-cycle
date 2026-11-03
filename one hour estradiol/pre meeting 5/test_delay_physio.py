"""Does a time delay in E2->GnRH create complexity at PHYSIOLOGICAL steepness (n<=5)?
Delayed feedback: GnRH pulse period & amplitude respond to E2(t - tau)."""
import numpy as np
from scipy.signal import find_peaks
from cascade_model import run_delay
from orbit_period import orbit_period, period_label

def op(na, nf, tau):
    T, E = run_delay(na, nf, tau, t_end=14000.0)
    m = T > T[-1] - 7000
    Tt, Ee = T[m], E[m]; rng = Ee.max() - Ee.min()
    if rng < 1: return 0
    pk, _ = find_peaks(Ee, prominence=0.15*rng)
    return orbit_period(Ee[pk], np.diff(Tt[pk]))

print("PHYSIOLOGICAL steepness (n_amp = n_freq = n), vs delay tau (min):")
print(f"  {'':6s}" + "".join(f"tau={t:<4d}" for t in [0, 15, 30, 60, 120, 240]))
for n in [2, 3, 4, 5]:
    row = f"  n={n}:  "
    for tau in [0, 15, 30, 60, 120, 240]:
        row += f"{period_label(op(n, n, tau)):<8s}"
    print(row)
