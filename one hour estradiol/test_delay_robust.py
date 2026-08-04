"""Is the delay-induced high-period at physiological steepness ROBUST?
Stress-test n=3, tau=120 min against ICs, settling, dt, prominence."""
import numpy as np
from scipy.signal import find_peaks
from cascade_model import run_delay
from orbit_period import orbit_period, period_label

def op(tau=120, t_end=14000.0, dt=0.25, e0=15.0, prom=0.15, tol=0.03):
    T, E = run_delay(3, 3, tau, t_end=t_end, dt=dt, e0=e0)
    m = T > T[-1] - t_end*0.5
    Tt, Ee = T[m], E[m]; rng = Ee.max() - Ee.min()
    if rng < 1: return 0
    pk, _ = find_peaks(Ee, prominence=prom*rng)
    return orbit_period(Ee[pk], np.diff(Tt[pk]), tol=tol)

print("=== n=3, tau=120 robustness ===")
for label, kw in [
    ("default",                  dict()),
    ("longer+finer (t=28k,dt=.125)", dict(t_end=28000.0, dt=0.125)),
    ("IC e0=8",                  dict(e0=8.0)),
    ("IC e0=30",                 dict(e0=30.0)),
    ("IC e0=45",                 dict(e0=45.0)),
    ("prominence 0.10",          dict(prom=0.10)),
    ("prominence 0.20",          dict(prom=0.20)),
]:
    print(f"  {label:28s}: {period_label(op(**kw))}")

# also: does a small IC change diverge? (sensitive dependence?)
T1, E1 = run_delay(3, 3, 120, t_end=14000.0, e0=15.0)
T2, E2 = run_delay(3, 3, 120, t_end=14000.0, e0=15.000001)
sep = np.abs(E1 - E2)
print(f"\nIC-sensitivity (1e-6 perturb): final |dE2| = {sep[-4000:].mean():.4f}  "
      f"(>0.5 = CHAOS; ~0 = periodic/quasi)")
