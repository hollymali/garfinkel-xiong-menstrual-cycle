"""Does the delay-induced high-period survive LOOSE tolerances (real ragged),
or collapse to low-period (harsh-tolerance artifact)?
Control: the tau=0 low-n 'yellow spot' (known fragile)."""
import numpy as np
from scipy.signal import find_peaks
from cascade_model import run_delay
from orbit_period import orbit_period, period_label

def op(na, nf, tau, tol):
    T, E = run_delay(na, nf, tau)
    m = T > T[-1] - 8000
    Tt, Ee = T[m], E[m]; rng = Ee.max()-Ee.min()
    pk, _ = find_peaks(Ee, prominence=0.15*rng)
    h = Ee[pk]
    hspread = (h.max()-h.min())/np.median(h)
    return orbit_period(h, np.diff(Tt[pk]), tol=tol), hspread

tols = [0.01, 0.03, 0.05, 0.10, 0.15, 0.20]
print("DELAY case  n=3, tau=120:")
for tol in tols:
    p, spr = op(3, 3, 120, tol)
    print(f"  tol {tol*100:>2.0f}% : {period_label(p):<18s}  (peak-height spread = {spr*100:.0f}% of median)")

print("\nCONTROL (fragile tau=0 spot)  n_amp=3, n_freq=5:")
for tol in tols:
    p, spr = op(3, 5, 0, tol)
    print(f"  tol {tol*100:>2.0f}% : {period_label(p):<18s}  (peak-height spread = {spr*100:.0f}% of median)")
