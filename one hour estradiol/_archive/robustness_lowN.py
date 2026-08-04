"""Are the low-n 'interesting spots' robust, or artifacts?
Test the yellow spot (n_amp=3, n_freq=5) against: initial conditions, longer
settling + finer dt, and the peak-detection prominence threshold."""
import numpy as np
from scipy.signal import find_peaks
from pi_cascade_combined import run
from e2metrics import orbit_period, period_label

K = 14.8

def op(na, nf, t_end=10000.0, dt=0.5, e0=15.0, tol=0.03, prom=0.15):
    T, G, F, E, pt, pa = run(na, nf, K, K, t_end=t_end, dt=dt, y0=(0.0, 5.0, e0))
    m = T > T[-1] - t_end*0.5
    Tt, Ee = T[m], E[m]
    rng = Ee.max() - Ee.min()
    if rng < 1: return 0
    pk, _ = find_peaks(Ee, prominence=prom*rng)
    return orbit_period(Ee[pk], np.diff(Tt[pk]), tol=tol)

print("=== YELLOW SPOT (n_amp=3, n_freq=5) ===")
tests = [
    ("default (t=10k, dt=.5, e0=15)", dict()),
    ("longer + finer (t=24k, dt=.25)", dict(t_end=24000.0, dt=0.25)),
    ("initial cond. e0=8",            dict(e0=8.0)),
    ("initial cond. e0=25",           dict(e0=25.0)),
    ("initial cond. e0=40",           dict(e0=40.0)),
    ("prominence 0.10 (looser)",      dict(prom=0.10)),
    ("prominence 0.20 (stricter)",    dict(prom=0.20)),
    ("tol 1%",                        dict(tol=0.01)),
]
for label, kw in tests:
    print(f"  {label:34s}: {period_label(op(3, 5, **kw))}")

print("\n=== is it isolated? fine grid around (3,5) ===")
for na in [2.6, 2.8, 3.0, 3.2, 3.4]:
    row = f"  n_amp={na}: "
    for nf in [4.6, 4.8, 5.0, 5.2, 5.4]:
        row += f"{period_label(op(na, nf)):>10s}"
    print(row)
